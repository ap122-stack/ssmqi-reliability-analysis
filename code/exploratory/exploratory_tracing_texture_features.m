function exploratory_tracing_texture_features()
% EXPLORATORY (not preregistered; written after the final registered run).
% Computes a fixed, declared list of grayscale and texture features inside every accepted tracing,
% on the same displayed B-mode panel, same artifact exclusion, and same per-image 0-255 rescale
% as the registered echo-intensity measure. One row per tracing; no clinical data are read here.
% Output: FINAL_20261106_190031\Exploratory\Tracing_Texture_Features.csv
studyDir = 'C:\Users\Owner\OneDrive\Rotator Cuff MVI\Reliability Study\Inter_Intra_Rater_Study';
addpath('C:\Users\Owner\OneDrive\Rotator Cuff MVI\Scripts');
T = ssmqi_load_registered_measurements(studyDir);
specs = {1,1,'Rater_01'; 1,2,'Rater_01'; 2,1,'Rater_02'; 2,2,'Rater_02'};
names = {'EI_mean','Gray_SD','Gray_CV','Gray_skew','Gray_kurt','Gray_P10','Gray_P90','Gray_IQR','Hist_entropy', ...
    'Local_SD5','GLCM_contrast','GLCM_homog','GLCM_energy','GLCM_corr','EI_ratio_to_panel'};
n = height(T); F = nan(n,numel(names));
for i = 1:n
    if ~T.Complete(i), continue; end
    folder = specs{[specs{:,1}] == T.Rater(i) & [specs{:,2}] == T.Session(i), 3};
    imgPath = fullfile(studyDir,folder,sprintf('Session_%d',T.Session(i)),'Blinded_Images',[char(T.Case_ID(i)) '.dcm']);
    [gray,artifact] = displayed_image(imgPath);
    use = logical(T.Mask{i}) & ~artifact;
    v = double(gray(use));
    if numel(v) < 50, continue; end
    mu = mean(v); sd = std(v);
    cnt = histcounts(v,0:256); pr = cnt/sum(cnt); pr = pr(pr>0);
    % local standard deviation in 5x5 windows, averaged inside the use mask
    g = double(gray); m1 = imboxfilt(g,5); m2 = imboxfilt(g.^2,5); lsd = sqrt(max(m2 - m1.^2,0));
    % GLCM (horizontal neighbour, 16 levels) from pixel pairs that are both inside the use mask
    q = floor(g/16); q(q>15) = 15; A = use(:,1:end-1) & use(:,2:end);
    a = q(:,1:end-1); b = q(:,2:end);
    P = accumarray([a(A)+1, b(A)+1],1,[16 16]); P = P + P'; P = P/sum(P(:));
    [I,J] = meshgrid(0:15,0:15); mi = sum(sum(I.*P)); mj = sum(sum(J.*P));
    si = sqrt(sum(sum((I-mi).^2.*P))); sj = sqrt(sum(sum((J-mj).^2.*P)));
    glcmCorr = NaN; if si>0 && sj>0, glcmCorr = sum(sum((I-mi).*(J-mj).*P))/(si*sj); end
    panelMean = mean(g(~artifact));
    F(i,:) = [mu, sd, sd/mu, skew(v), kurt(v), prctile(v,10), prctile(v,90), iqr(v), -sum(pr.*log2(pr)), ...
        mean(lsd(use)), sum(sum((I-J).^2.*P)), sum(sum(P./(1+abs(I-J)))), sum(P(:).^2), glcmCorr, mu/panelMean];
end
out = table(T.Rater,T.Session,T.Source_Image,T.Participant,T.Axis,T.Complete,'VariableNames',{'Rater','Session','Source_Image','Participant','Axis','Complete'});
out = [out array2table(F,'VariableNames',names)];
outDir = fullfile(studyDir,'Analysis','Registered','FINAL_20261106_190031','Exploratory');
if ~exist(outDir,'dir'), mkdir(outDir); end
writetable(out,fullfile(outDir,'Tracing_Texture_Features.csv'));
fprintf('rows %d, with features %d\n',n,nnz(~isnan(F(:,1))));
end

function s = skew(v), m = mean(v); d = v-m; s = mean(d.^3)/max(std(v,1)^3,eps); end
function k = kurt(v), m = mean(v); d = v-m; k = mean(d.^4)/max(std(v,1)^4,eps); end

function [gray,artifact] = displayed_image(path)
info = dicominfo(path);
img = squeeze(dicomread(info));
[bmode,originXY] = crop_left(img,info);
artifact = ssmqi_formal_artifact_mask(bmode,info,originXY);
if ndims(bmode) == 3, g = rgb2gray(bmode); else, g = bmode; end
gray = uint8(rescale(double(g),0,255));
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
