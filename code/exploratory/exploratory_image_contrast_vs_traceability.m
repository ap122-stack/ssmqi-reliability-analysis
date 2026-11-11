function exploratory_image_contrast_vs_traceability()
% EXPLORATORY (not preregistered; written after the final registered run).
% Question: do objective image-contrast measures differ between images Rater 2 marked
% untraceable in Session 2 and images Rater 2 traced? This tests the rater's stated
% reason (contrast too high) against the pixels instead of relying on the statement.
% Reads blinded DICOM copies and the Session 2 case list; writes one CSV; does not touch
% any registered output. Uses the same displayed B-mode crop and artifact mask as the
% registered pipeline. Features are computed BEFORE the per-image 0-255 rescale used for echo intensity.
root = 'C:\Users\Owner\OneDrive\Rotator Cuff MVI\Reliability Study\Inter_Intra_Rater_Study';
outDir = fullfile(root,'Analysis','Registered','FINAL_20261106_190031','Exploratory');
if ~exist(outDir,'dir'), mkdir(outDir); end
addpath('C:\Users\Owner\OneDrive\Rotator Cuff MVI\Scripts');
cases = readtable(fullfile(root,'Rater_02','Session_2','Session_2_Case_List.csv'),'TextType','string','VariableNamingRule','preserve');
key = readtable(fullfile(root,'Admin_Private','Rater_02_Session_2_Randomization_Key.csv'),'TextType','string','VariableNamingRule','preserve');
n = height(cases);
rows = cell(n,10);
for i = 1:n
    cid = cases.Case_ID(i);
    info = dicominfo(fullfile(root,'Rater_02','Session_2','Blinded_Images',char(cid)+".dcm"));
    img = squeeze(dicomread(info));
    [bmode,originXY] = crop_left(img,info);
    art = ssmqi_formal_artifact_mask(bmode,info,originXY);
    if ndims(bmode) == 3, g = double(rgb2gray(bmode)); else, g = double(bmode); end
    v = g(~art);
    k = find(key.Case_ID == cid,1);
    traced = lower(strtrim(cases.Status(i))) == "complete";
    p = prctile(v,[5 50 95]);
    rows(i,:) = {char(key.Source_Image(k)), string(key.Patient_Group_Inferred(k)), char(cases.Axis(i)), traced, ...
        mean(v), std(v), p(3)-p(1), mean(v >= 250), mean(v <= 5), numel(v)};
end
T = cell2table(rows,'VariableNames',{'Source_Image','Participant','Axis','Traced','MeanGray','SDGray','P95_minus_P5','FracSaturated','FracBlack','Pixels'});
writetable(T,fullfile(outDir,'Image_Contrast_Features_R2S2.csv'));
% Per-feature comparison with participant-cluster bootstrap (seed fixed, 2000 resamples)
feats = {'MeanGray','SDGray','P95_minus_P5','FracSaturated','FracBlack'};
rng(20261108);
parts = unique(T.Participant);
res = cell(numel(feats),6);
for f = 1:numel(feats)
    x = T.(feats{f}); tr = T.Traced;
    auc = auc_fun(x(~tr),x(tr));          % P(untraceable value > traced value)
    bs = nan(2000,1);
    for b = 1:2000
        pick = parts(randi(numel(parts),numel(parts),1));
        idx = []; for q = 1:numel(pick), idx = [idx; find(T.Participant == pick(q))]; end %#ok<AGROW>
        xb = x(idx); tb = tr(idx);
        if any(tb) && any(~tb), bs(b) = auc_fun(xb(~tb),xb(tb)); end
    end
    ci = prctile(bs(~isnan(bs)),[2.5 97.5]);
    res(f,:) = {feats{f}, median(x(~tr)), median(x(tr)), auc, ci(1), ci(2)};
end
R = cell2table(res,'VariableNames',{'Feature','Median_Untraceable','Median_Traced','AUC_untraceable_higher','CI_low','CI_high'});
writetable(R,fullfile(outDir,'Image_Contrast_vs_Traceability_R2S2.csv'));
fprintf('n images %d, untraceable %d, traced %d\n',n,nnz(~T.Traced),nnz(T.Traced));
disp(R);
end

function a = auc_fun(pos,neg)
% Probability that a random "pos" (untraceable) value exceeds a random "neg" (traced) value; ties count half.
[P,N] = meshgrid(pos,neg);
a = mean((P(:) > N(:)) + 0.5*(P(:) == N(:)));
end

function [cropped,originXY] = crop_left(img,info)
[h,w,~] = size(img); originXY = [1 1]; bestArea = -Inf; bestBounds = [];
if isfield(info,'SequenceOfUltrasoundRegions')
    names = fieldnames(info.SequenceOfUltrasoundRegions);
    for i = 1:numel(names)
        r = info.SequenceOfUltrasoundRegions.(names{i});
        f = {'RegionLocationMinX0','RegionLocationMinY0','RegionLocationMaxX1','RegionLocationMaxY1'};
        if ~all(isfield(r,f)), continue; end
        x0 = double(r.RegionLocationMinX0); x1 = double(r.RegionLocationMaxX1);
        y0 = double(r.RegionLocationMinY0); y1 = double(r.RegionLocationMaxY1);
        area = max(0,x1-x0+1)*max(0,y1-y0+1);
        if (x0+x1)/2 < w/2 && area > bestArea, bestArea = area; bestBounds = [x0 y0 x1 y1]; end
    end
end
if ~isempty(bestBounds)
    x0 = max(0,round(bestBounds(1))); y0 = max(0,round(bestBounds(2)));
    x1 = min(w-1,round(bestBounds(3))); y1 = min(h-1,round(bestBounds(4)));
    cropped = img(y0+1:y1+1,x0+1:x1+1,:); originXY = [x0+1 y0+1]; return;
end
if w >= 1.15*h, cropped = img(:,1:floor(w/2),:); else, cropped = img; end
end
