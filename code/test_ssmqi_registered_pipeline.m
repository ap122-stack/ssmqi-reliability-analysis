function test_ssmqi_registered_pipeline()
% TEST_SSMQI_REGISTERED_PIPELINE  Runs the full registered analysis on synthetic data.
% The synthetic study has the same layout as the real one: 80 images, 17 participants,
% 57 long-axis and 23 short-axis images, four sessions, some images untraceable.
% Known truth is built in so the output can be checked, not only that it runs.

rng(7,'twister');
nImg = 80; nPart = 17;
part = string(1 + mod(0:nImg-1,nPart))';
part = "G" + pad(part,3,'left','0');
axisLab = [repmat("Long axis",57,1); repmat("Short axis",23,1)];
imgName = "IMG" + string((1:nImg)');
trueEI = 110 + 12*randn(nPart,1);      % participant-level echo intensity
baseEI = trueEI(double(extractAfter(part,1))) + 6*randn(nImg,1);

sessions = [1 1; 1 2; 2 1; 2 2];
pTrace = [1.00 0.99 0.85 0.85];
raterShift = [0 2.5];
rows = {};
for s = 1:4
    for i = 1:nImg
        complete = rand < pTrace(s);
        mask = []; ei = NaN; nPix = 0; mm = 0.081;
        if complete
            [X,Y] = meshgrid(1:260,1:200);
            cx = 130 + 4*randn; cy = 100 + 4*randn;
            ax = 70*(1+0.05*randn); by = 35*(1+0.05*randn);
            if axisLab(i) == "Short axis", ax = 55*(1+0.05*randn); by = 40*(1+0.05*randn); end
            mask = ((X-cx)/ax).^2 + ((Y-cy)/by).^2 <= 1;
            nPix = nnz(mask);
            ei = baseEI(i) + raterShift(sessions(s,1)) + 4*randn;
        end
        rows(end+1,:) = {sessions(s,1),sessions(s,2),"C"+string(s)+"_"+string(i),imgName(i),part(i), ...
            axisLab(i),complete,{mask},mm,ei,nPix}; %#ok<AGROW>
    end
end
T = cell2table(rows,'VariableNames',{'Rater','Session','Case_ID','Source_Image','Participant', ...
    'Axis','Complete','Mask','mmPerPx','EI','EI_Pixels'});

outDir = fullfile(tempdir,'ssmqi_pipeline_test');
if isfolder(outDir), rmdir(outDir,'s'); end
out = ssmqi_registered_analysis(T,outDir,'nBoot',40,'seed',1);

problems = {};
need = {'Table2_Outline_Agreement.csv','Table3_EchoIntensity_Agreement.csv','Table4_Variance_Components_ICC.csv', ...
    'Table5_Decision_Study.csv','Table6_Traceability_Agreement.csv','Pair_Level_Metrics.csv'};
for k = 1:numel(need)
    if ~isfile(fullfile(outDir,need{k})), problems{end+1} = ['missing ' need{k}]; end %#ok<AGROW>
end
if isempty(dir(fullfile(outDir,'BlandAltman_*.png'))), problems{end+1} = 'no Bland-Altman figures'; end
if height(out.outline) ~= 10, problems{end+1} = sprintf('expected 10 outline rows, got %d',height(out.outline)); end
if height(out.echoIntensity) ~= 10, problems{end+1} = sprintf('expected 10 echo rows, got %d',height(out.echoIntensity)); end
if height(out.traceability) ~= 10, problems{end+1} = sprintf('expected 10 traceability rows, got %d',height(out.traceability)); end
if height(out.variance) ~= 2, problems{end+1} = sprintf('expected 2 variance rows, got %d',height(out.variance)); end
if height(out.decisionStudy) ~= 8, problems{end+1} = sprintf('expected 8 decision rows, got %d',height(out.decisionStudy)); end
if any(out.outline.Dice_Mean < 0.80 | out.outline.Dice_Mean > 1), problems{end+1} = 'Dice outside the expected range for jittered ellipses'; end
if any(out.outline.Dice_CI_Low > out.outline.Dice_Mean | out.outline.Dice_CI_High < out.outline.Dice_Mean), problems{end+1} = 'Dice CI does not bracket the mean'; end
if any(out.outline.HD95_mm_Mean <= 0), problems{end+1} = 'HD95 not positive'; end
inter = out.echoIntensity(contains(out.echoIntensity.Comparison,'Inter-rater, primary'),:);
if any(inter.Bias_Percent < 0) , problems{end+1} = 'inter-rater bias should be positive (rater 2 built in +2.5)'; end
if any(out.echoIntensity.SEM <= 0), problems{end+1} = 'SEM not positive'; end
if any(out.variance.ICC_A1 < 0.2 | out.variance.ICC_A1 > 1), problems{end+1} = 'ICC outside the plausible range'; end
if any(out.traceability.PercentAgreement < 0 | out.traceability.PercentAgreement > 100), problems{end+1} = 'agreement out of range'; end

% Sensitivity path: flag a few tracings and confirm they are dropped from the metrics only
flagged = table(uint8([1;1;2]),uint8([1;1;1]),["IMG3";"IMG8";"IMG12"], ...
    'VariableNames',{'Rater','Session','Source_Image'});
flagged.Rater = double(flagged.Rater); flagged.Session = double(flagged.Session);
out2 = ssmqi_registered_analysis(T,fullfile(tempdir,'ssmqi_pipeline_test_sens'),'nBoot',20,'seed',1,'flagged',flagged);
if height(out2.traceability) ~= height(out.traceability), problems{end+1} = 'flagging changed the traceability table'; end
nA = sum(out.outline.Pairs); nB = sum(out2.outline.Pairs);
if ~(nB <= nA), problems{end+1} = 'flagging did not reduce pairs'; end

if isempty(problems)
    fprintf('PASS  full registered pipeline ran to completion on synthetic data\n');
    fprintf('      outline rows=%d, echo rows=%d, variance rows=%d, decision rows=%d, traceability rows=%d\n', ...
        height(out.outline),height(out.echoIntensity),height(out.variance),height(out.decisionStudy),height(out.traceability));
    disp(out.variance(:,{'Plane','Tracings','Participants','ICC_A1','ICC_CI_Low','ICC_CI_High'}));
else
    for k = 1:numel(problems), fprintf('FAIL  %s\n',problems{k}); end
    error('Pipeline test found %d problem(s).',numel(problems));
end
end
