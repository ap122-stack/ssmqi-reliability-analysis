function T = ssmqi_synthetic_study_table(seed,dropSession)
% SSMQI_SYNTHETIC_STUDY_TABLE  Synthetic study with the real layout, used only by the tests.
% 80 images, 17 participants, 57 long-axis and 23 short-axis, four sessions.
% dropSession is [rater session] to leave out of the table, or [] to keep all.
if nargin < 2, dropSession = []; end
rng(seed,'twister');
nImg = 80; nPart = 17;
part = "G" + pad(string(1 + mod(0:nImg-1,nPart))',3,'left','0');
axisLab = [repmat("Long axis",57,1); repmat("Short axis",23,1)];
imgName = "IMG" + string((1:nImg)');
trueEI = 110 + 12*randn(nPart,1);
baseEI = trueEI(double(extractAfter(part,1))) + 6*randn(nImg,1);
sessions = [1 1; 1 2; 2 1; 2 2];
pTrace = [1.00 0.99 0.85 0.85];
raterShift = [0 2.5];
[X,Y] = meshgrid(1:260,1:200);
rows = {};
for s = 1:4
    if ~isempty(dropSession) && isequal(sessions(s,:),dropSession), continue; end
    for i = 1:nImg
        complete = rand < pTrace(s);
        mask = []; ei = NaN; nPix = 0;
        if complete
            cx = 130 + 4*randn; cy = 100 + 4*randn;
            ax = 70*(1+0.05*randn); by = 35*(1+0.05*randn);
            if axisLab(i) == "Short axis", ax = 55*(1+0.05*randn); by = 40*(1+0.05*randn); end
            mask = ((X-cx)/ax).^2 + ((Y-cy)/by).^2 <= 1;
            nPix = nnz(mask);
            ei = baseEI(i) + raterShift(sessions(s,1)) + 4*randn;
        end
        rows(end+1,:) = {sessions(s,1),sessions(s,2),"C"+string(s)+"_"+string(i),imgName(i),part(i), ...
            axisLab(i),complete,{mask},0.081,ei,nPix}; %#ok<AGROW>
    end
end
T = cell2table(rows,'VariableNames',{'Rater','Session','Case_ID','Source_Image','Participant', ...
    'Axis','Complete','Mask','mmPerPx','EI','EI_Pixels'});
end
