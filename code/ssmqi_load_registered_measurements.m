function T = ssmqi_load_registered_measurements(studyDir)
% SSMQI_LOAD_REGISTERED_MEASUREMENTS  Collect every tracing into one table.
%
% Reads each rater session's case list, randomization key and saved masks and
% computes echo intensity (mean grayscale inside the outline, excluding scanner
% overlay pixels) on the image that session's rater saw.
%
% Output columns
%   Rater       1 or 2
%   Session     1 or 2 (current folder numbering)
%   Case_ID     blinded identifier
%   Source_Image, Participant, Axis
%   Complete    true when the rater accepted a tracing
%   Mask        logical mask in the cropped B-mode frame (empty if not traced)
%   mmPerPx     millimeters per pixel from the DICOM header
%   EI          echo intensity, 0 to 255 scale (NaN if not traced)
%   EI_Pixels   number of pixels averaged

specs = {1,1,'Rater_01'; 1,2,'Rater_01'; 2,1,'Rater_02'; 2,2,'Rater_02'};
adminDir = fullfile(studyDir,'Admin_Private');
rows = {};
for s = 1:size(specs,1)
    raterNum = specs{s,1}; sessNum = specs{s,2}; folder = specs{s,3};
    sessDir = fullfile(studyDir,folder,sprintf('Session_%d',sessNum));
    listPath = fullfile(sessDir,sprintf('Session_%d_Case_List.csv',sessNum));
    keyPath = fullfile(adminDir,sprintf('%s_Session_%d_Randomization_Key.csv',folder,sessNum));
    if ~isfile(listPath) || ~isfile(keyPath), continue; end
    cases = ssmqi_read_csv_strings(listPath);
    key = ssmqi_read_csv_strings(keyPath);
    for i = 1:height(cases)
        caseId = char(cases.Case_ID(i));
        k = find(strcmp(string(key.Case_ID),cases.Case_ID(i)),1);
        if isempty(k), error('Case %s has no entry in %s.',caseId,keyPath); end
        status = lower(strtrim(string(cases.Status(i))));
        isComplete = status == "complete";
        mask = []; mmPerPx = NaN; ei = NaN; nPix = 0;
        if isComplete
            matPath = fullfile(sessDir,'Trace_Results',[caseId '_manual_trace.mat']);
            if ~isfile(matPath), error('Missing saved trace: %s',matPath); end
            S = load(matPath,'mask','cmPerPixel');
            mask = logical(S.mask);
            mmPerPx = 10 * S.cmPerPixel;
            imgPath = fullfile(sessDir,'Blinded_Images',[caseId '.dcm']);
            [gray,artifact] = displayed_image(imgPath);
            if ~isequal(size(gray),size(mask))
                error('Mask and displayed image sizes differ for %s.',caseId);
            end
            use = mask & ~artifact;
            nPix = nnz(use);
            if nPix > 0, ei = mean(double(gray(use))); end
        end
        rows(end+1,:) = {raterNum,sessNum,string(caseId),key.Source_Image(k), ...
            key.Patient_Group_Inferred(k),string(cases.Axis(i)),isComplete,{mask},mmPerPx,ei,nPix}; %#ok<AGROW>
    end
end
T = cell2table(rows,'VariableNames',{'Rater','Session','Case_ID','Source_Image', ...
    'Participant','Axis','Complete','Mask','mmPerPx','EI','EI_Pixels'});
end

function [gray,artifact] = displayed_image(path)
info = dicominfo(path);
img = squeeze(dicomread(info));
[bmode,originXY] = crop_left_bmode_region(img,info);
artifact = ssmqi_formal_artifact_mask(bmode,info,originXY);
if ndims(bmode) == 3, g = rgb2gray(bmode); else, g = bmode; end
gray = uint8(rescale(double(g),0,255));
end

function [cropped,originXY] = crop_left_bmode_region(img,info)
% Same rule as the tracing interface (run_ssmqi_formal_reliability_session.m).
[h,w,~] = size(img);
originXY = [1 1];
bestArea = -Inf;
bestBounds = [];
if isfield(info,'SequenceOfUltrasoundRegions')
    names = fieldnames(info.SequenceOfUltrasoundRegions);
    for i = 1:numel(names)
        r = info.SequenceOfUltrasoundRegions.(names{i});
        fields = {'RegionLocationMinX0','RegionLocationMinY0', ...
            'RegionLocationMaxX1','RegionLocationMaxY1'};
        if ~all(isfield(r,fields)), continue; end
        x0 = double(r.RegionLocationMinX0); x1 = double(r.RegionLocationMaxX1);
        y0 = double(r.RegionLocationMinY0); y1 = double(r.RegionLocationMaxY1);
        area = max(0,x1-x0+1)*max(0,y1-y0+1);
        if (x0+x1)/2 < w/2 && area > bestArea
            bestArea = area;
            bestBounds = [x0 y0 x1 y1];
        end
    end
end
if ~isempty(bestBounds)
    x0 = max(0,round(bestBounds(1))); y0 = max(0,round(bestBounds(2)));
    x1 = min(w-1,round(bestBounds(3))); y1 = min(h-1,round(bestBounds(4)));
    cropped = img(y0+1:y1+1,x0+1:x1+1,:);
    originXY = [x0+1 y0+1];
    return;
end
if w >= 1.15*h
    cropped = img(:,1:floor(w/2),:);
else
    cropped = img;
end
end
