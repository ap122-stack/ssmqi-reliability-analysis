function run_ssmqi_registered_analysis(mode,varargin)
% RUN_SSMQI_REGISTERED_ANALYSIS  Controlled entry point for the preregistered analysis.
%
%   run_ssmqi_registered_analysis('freeze')   run the tests, then record SHA-256 hashes of the
%                                             analysis code. Do this once, after all data are
%                                             collected and before the final run.
%   run_ssmqi_registered_analysis('run')      verify the frozen code, check that every session is
%                                             finished, run the registered analysis once, and
%                                             write a manifest. A second final run is refused.
%
% Controls that keep the analysis from being run repeatedly until it looks good:
%   1. The code must match the frozen hashes or the run stops.
%   2. Every execution, including incomplete-data previews, is appended to
%      Analysis_Run_Log.csv and can never be edited by this code.
%   3. After one complete final run a lock file is written; a second run needs a written
%      reason, is saved in its own folder, and is logged as a rerun.
%   4. The data checked at run time are fingerprinted in Input_Files_SHA256.csv.
%
% Options (name-value): 'nBoot' (default 2000), 'seed' (default 20261024),
% 'allowIncomplete' (preview only), 'rerunReason' (text), 'testMode' with 'tableOverride'
% and 'rootOverride' (used only by the automated tests).

p = inputParser;
addParameter(p,'nBoot',2000);
addParameter(p,'seed',20261024);
addParameter(p,'allowIncomplete',false);
addParameter(p,'rerunReason','');
addParameter(p,'testMode',false);
addParameter(p,'tableOverride',table());
addParameter(p,'rootOverride','');
parse(p,varargin{:});
o = p.Results;

scriptsDir = fileparts(mfilename('fullpath'));
projectRoot = fileparts(scriptsDir);
if ~isempty(o.rootOverride), projectRoot = o.rootOverride; end
studyDir = fullfile(projectRoot,'Reliability Study','Inter_Intra_Rater_Study');
analysisRoot = fullfile(studyDir,'Analysis','Registered');
if ~isfolder(analysisRoot), mkdir(analysisRoot); end
freezePath = fullfile(analysisRoot,'Analysis_Code_Freeze.csv');
logPath = fullfile(analysisRoot,'Analysis_Run_Log.csv');
lockPath = fullfile(analysisRoot,'ANALYSIS_FINAL_RUN_LOCK.json');
codeFiles = {'ssmqi_stats_core.m','ssmqi_registered_analysis.m', ...
    'ssmqi_load_registered_measurements.m','run_ssmqi_registered_analysis.m', ...
    'ssmqi_formal_artifact_mask.m','ssmqi_read_csv_strings.m','ssmqi_sha256_file.m'};

switch lower(string(mode))
    case "freeze"
        if isfile(freezePath)
            error('The code is already frozen (%s). Delete nothing; document any change as a deviation.',freezePath);
        end
        fprintf('Running tests before freezing...\n');
        test_ssmqi_stats_core();
        test_ssmqi_registered_pipeline();
        rows = cell(numel(codeFiles),5);
        stamp = string(datetime('now','Format','yyyy-MM-dd HH:mm:ss'));
        for k = 1:numel(codeFiles)
            rows(k,:) = {codeFiles{k},ssmqi_sha256_file(fullfile(scriptsDir,codeFiles{k})),char(stamp), ...
                version,'https://osf.io/gfd6t'};
        end
        writetable(cell2table(rows,'VariableNames',{'File','SHA256','Frozen_At','MATLAB_Version','Registration'}),freezePath);
        append_log(logPath,'freeze','','tests passed; code hashes recorded');
        fprintf('Code frozen. Hashes written to %s\n',freezePath);

    case "run"
        verify_freeze(freezePath,scriptsDir);
        locked = isfile(lockPath);
        if locked && isempty(o.rerunReason)
            error(['A final run was already completed (%s). To run again, supply a written reason with ' ...
                '''rerunReason''; the rerun will be logged and must be reported in the manuscript.'],lockPath);
        end

        if o.testMode && ~isempty(o.tableOverride)
            T = o.tableOverride;
        else
            T = ssmqi_load_registered_measurements(studyDir);
        end
        complete = true; note = 'all four sessions finished';
        for s = [1 1; 1 2; 2 1; 2 2]'
            rows = T(T.Rater == s(1) & T.Session == s(2),:);
            if isempty(rows)
                complete = false; note = sprintf('Rater %d Session %d is missing',s(1),s(2)); break;
            end
            if o.testMode, continue; end
            listPath = fullfile(studyDir,sprintf('Rater_%02d',s(1)),sprintf('Session_%d',s(2)), ...
                sprintf('Session_%d_Case_List.csv',s(2)));
            st = lower(strtrim(ssmqi_read_csv_strings(listPath).Status));
            if any(st == "pending")
                complete = false; note = sprintf('Rater %d Session %d still has pending images',s(1),s(2)); break;
            end
        end
        if ~complete && ~o.allowIncomplete
            append_log(logPath,'run-refused','',note);
            error('Not all data are complete: %s. Finish collection before the final run.',note);
        end

        stamp = char(datetime('now','Format','yyyyMMdd_HHmmss'));
        if ~complete, tag = 'INCOMPLETE_PREVIEW'; elseif locked, tag = 'RERUN'; else, tag = 'FINAL'; end
        outDir = fullfile(analysisRoot,[tag '_' stamp]);
        mkdir(outDir);

        flagged = table();
        flagPath = fullfile(studyDir,'Admin_Private','Adjudication_Flagged_Outlines.csv');
        if isfile(flagPath)
            fl = ssmqi_read_csv_strings(flagPath);
            flagged = table(double(fl.Rater),double(fl.Session),fl.Source_Image, ...
                'VariableNames',{'Rater','Session','Source_Image'});
        end

        res = ssmqi_registered_analysis(T,fullfile(outDir,'Primary_AllTracings'),'nBoot',o.nBoot,'seed',o.seed); %#ok<NASGU>
        if height(flagged) > 0
            ssmqi_registered_analysis(T,fullfile(outDir,'Sensitivity_ExcludingAdjudicatedOutlines'),'nBoot',o.nBoot, ...
                'seed',o.seed,'flagged',flagged);
        end

        inputs = input_fingerprint(studyDir,o.testMode);
        if ~isempty(inputs), writetable(inputs,fullfile(outDir,'Input_Files_SHA256.csv')); end
        manifest = struct('run_type',tag,'timestamp',char(datetime('now','Format','yyyy-MM-dd HH:mm:ss')), ...
            'matlab',version,'seed',o.seed,'nBoot',o.nBoot,'registration','https://osf.io/gfd6t', ...
            'data_status',note,'tracings_total',height(T),'tracings_accepted',sum(T.Complete), ...
            'participants',numel(unique(T.Participant)),'sensitivity_run',height(flagged) > 0, ...
            'rerun_reason',o.rerunReason,'code_files',{codeFiles});
        fid = fopen(fullfile(outDir,'Run_Manifest.json'),'w'); fwrite(fid,jsonencode(manifest,'PrettyPrint',true)); fclose(fid);
        append_log(logPath,tag,outDir,note);
        if strcmp(tag,'FINAL')
            fid = fopen(lockPath,'w'); fwrite(fid,jsonencode(struct('final_run_folder',outDir,'timestamp',manifest.timestamp),'PrettyPrint',true)); fclose(fid);
        end
        fprintf('%s analysis complete. Output: %s\n',tag,outDir);

    otherwise
        error('Mode must be ''freeze'' or ''run''.');
end
end

function verify_freeze(freezePath,scriptsDir)
if ~isfile(freezePath)
    error('The analysis code has not been frozen. Run run_ssmqi_registered_analysis(''freeze'') first.');
end
fz = ssmqi_read_csv_strings(freezePath);
for k = 1:height(fz)
    current = ssmqi_sha256_file(fullfile(scriptsDir,char(fz.File(k))));
    if ~strcmpi(current,char(fz.SHA256(k)))
        error('Code-freeze violation: %s changed after the freeze. Document the change as a deviation and decide how to proceed.',char(fz.File(k)));
    end
end
end

function append_log(logPath,kind,outDir,note)
row = table(string(datetime('now','Format','yyyy-MM-dd HH:mm:ss')),string(kind),string(outDir),string(note), ...
    'VariableNames',{'Timestamp','Event','Output_Folder','Note'});
if isfile(logPath)
    writetable([ssmqi_read_csv_strings(logPath); row],logPath);
else
    writetable(row,logPath);
end
end

function tbl = input_fingerprint(studyDir,testMode)
tbl = table();
if testMode, return; end
files = {};
for r = {'Rater_01','Rater_02'}
    for s = 1:2
        d = fullfile(studyDir,r{1},sprintf('Session_%d',s));
        if isfolder(d)
            files = [files; {fullfile(d,sprintf('Session_%d_Case_List.csv',s))}]; %#ok<AGROW>
            m = dir(fullfile(d,'Trace_Results','*_manual_trace.mat'));
            files = [files; fullfile({m.folder},{m.name})']; %#ok<AGROW>
        end
        k = fullfile(studyDir,'Admin_Private',sprintf('%s_Session_%d_Randomization_Key.csv',r{1},s));
        if isfile(k), files = [files; {k}]; end %#ok<AGROW>
    end
end
files = files(cellfun(@isfile,files));
h = cell(numel(files),1);
for i = 1:numel(files), h{i} = ssmqi_sha256_file(files{i}); end
tbl = table(string(files),string(h),'VariableNames',{'File','SHA256'});
end
