function test_ssmqi_runner_guards()
% TEST_SSMQI_RUNNER_GUARDS  Confirms the controls in run_ssmqi_registered_analysis work.
% Uses a throwaway project folder and synthetic data; no study files are read or changed.

root = fullfile(tempdir,'ssmqi_guard_test');
if isfolder(root), rmdir(root,'s'); end
mkdir(root);
T = ssmqi_synthetic_study_table(3);
Tshort = ssmqi_synthetic_study_table(3,[2 2]);
problems = {};

    function expectError(name,fcn,pattern)
        try
            fcn();
            problems{end+1} = [name ': expected an error but none occurred']; %#ok<AGROW>
        catch ME
            if ~contains(ME.message,pattern)
                problems{end+1} = sprintf('%s: wrong error "%s"',name,ME.message); %#ok<AGROW>
            end
        end
    end

expectError('run before freeze',@() run_ssmqi_registered_analysis('run','testMode',true,'tableOverride',T,'rootOverride',root,'nBoot',10),'not been frozen');

run_ssmqi_registered_analysis('freeze','rootOverride',root);
freezePath = fullfile(root,'Reliability Study','Inter_Intra_Rater_Study','Analysis','Registered','Analysis_Code_Freeze.csv');
if ~isfile(freezePath), problems{end+1} = 'freeze file not written'; end
expectError('second freeze',@() run_ssmqi_registered_analysis('freeze','rootOverride',root),'already frozen');

expectError('incomplete data refused',@() run_ssmqi_registered_analysis('run','testMode',true,'tableOverride',Tshort,'rootOverride',root,'nBoot',10),'Not all data are complete');
run_ssmqi_registered_analysis('run','testMode',true,'tableOverride',Tshort,'rootOverride',root,'nBoot',10,'allowIncomplete',true);

run_ssmqi_registered_analysis('run','testMode',true,'tableOverride',T,'rootOverride',root,'nBoot',10);
reg = fullfile(root,'Reliability Study','Inter_Intra_Rater_Study','Analysis','Registered');
if ~isfile(fullfile(reg,'ANALYSIS_FINAL_RUN_LOCK.json')), problems{end+1} = 'lock not written after the final run'; end
if isempty(dir(fullfile(reg,'FINAL_*','Run_Manifest.json'))), problems{end+1} = 'final manifest missing'; end
if isempty(dir(fullfile(reg,'INCOMPLETE_PREVIEW_*'))), problems{end+1} = 'incomplete preview not saved separately'; end

expectError('second final run',@() run_ssmqi_registered_analysis('run','testMode',true,'tableOverride',T,'rootOverride',root,'nBoot',10),'already completed');
run_ssmqi_registered_analysis('run','testMode',true,'tableOverride',T,'rootOverride',root,'nBoot',10,'rerunReason','test of rerun path');
if isempty(dir(fullfile(reg,'RERUN_*'))), problems{end+1} = 'rerun folder missing'; end

lg = readtable(fullfile(reg,'Analysis_Run_Log.csv'),'TextType','string');
need = ["freeze","run-refused","INCOMPLETE_PREVIEW","FINAL","RERUN"];
for k = 1:numel(need)
    if ~any(lg.Event == need(k)), problems{end+1} = ['log missing event ' char(need(k))]; end %#ok<AGROW>
end

% Tamper: change a recorded hash and confirm the run stops
fz = readtable(freezePath,'TextType','string');
fz.SHA256(1) = "0000000000000000000000000000000000000000000000000000000000000000";
writetable(fz,freezePath);
expectError('tampered freeze',@() run_ssmqi_registered_analysis('run','testMode',true,'tableOverride',T,'rootOverride',root,'nBoot',10,'rerunReason','x'),'Code-freeze violation');

if isempty(problems)
    fprintf('PASS  runner controls: freeze required, tamper detected, incomplete refused, one final run, rerun logged\n');
else
    for k = 1:numel(problems), fprintf('FAIL  %s\n',problems{k}); end
    error('Runner guard test found %d problem(s).',numel(problems));
end
rmdir(root,'s');
end
