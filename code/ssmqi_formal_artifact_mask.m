function [artifactMask,summary] = ssmqi_formal_artifact_mask(image,info,cropOriginXY)
% Mask scanner overlays while preserving the native B-mode pixel geometry.
% The right SWE panel is removed before this function is called. The faint
% box in the left B-mode panel is mapped from its DICOM ultrasound-region
% metadata, avoiding unreliable visual rectangle detection.

if nargin < 2, info = struct(); end
if nargin < 3, cropOriginXY = [1 1]; end
if ndims(image) == 3, gray = rgb2gray(image); else, gray = image; end
gray = im2uint8(rescale(double(gray),0,1));
[h,w] = size(gray);
artifactMask = false(h,w);
labels = strings(3,1); labelCount = 0;

[sweBox,boxFound] = dicom_swe_box_in_left_panel(info,cropOriginXY,h,w);
if boxFound
    artifactMask = artifactMask | rectangle_border_mask(sweBox,h,w,4);
    labelCount = labelCount+1;
    labels(labelCount) = "DICOM-mapped SWE box border";
end

% Scanner text, the right depth scale, and the bottom body-position diagram
% are bright peripheral drawings. Mask only their pixels, not entire bands.
highThreshold = 180;
bright = double(gray) >= highThreshold;
peripheral = false(h,w);
peripheral(1:max(1,round(0.065*h)),:) = true;
peripheral(min(h,round(0.76*h)):end,1:max(1,round(0.27*w))) = true;
peripheralMask = imdilate(bright & peripheral,strel('disk',3));
peripheralMask = bwareaopen(peripheralMask,6);
if any(peripheralMask(:))
    artifactMask = artifactMask | peripheralMask;
    labelCount = labelCount+1;
    labels(labelCount) = "top labels/bottom body diagram";
end

% A depth scale is a thin, nearly vertical line near the right border. Use
% line geometry instead of brightness so bright lateral anatomy is retained.
rightZone = false(h,w);
rightZone(:,min(w,round(0.92*w)):end) = true;
scaleEdges = edge(gray,'canny',[0.04 0.18]) & rightZone;
scaleLines = imopen(scaleEdges,strel('line',max(25,round(0.05*h)),90));
scaleMask = imdilate(scaleLines,strel('disk',2));
if any(scaleMask(:))
    artifactMask = artifactMask | scaleMask;
    labelCount = labelCount+1;
    labels(labelCount) = "right depth-scale line";
end

artifactMask = imdilate(artifactMask,strel('disk',1));
labels = labels(1:labelCount);
if isempty(labels)
    summary = 'none detected';
else
    summary = char(strjoin(labels,'; '));
end
end

function [box,found] = dicom_swe_box_in_left_panel(info,cropOriginXY,h,w)
box = []; found = false;
if ~isstruct(info) || ~isfield(info,'SequenceOfUltrasoundRegions'), return; end
names = fieldnames(info.SequenceOfUltrasoundRegions);
regions = struct('x0',{},'y0',{},'x1',{},'y1',{},'area',{},'cx',{});
for i = 1:numel(names)
    region = info.SequenceOfUltrasoundRegions.(names{i});
    required = {'RegionLocationMinX0','RegionLocationMinY0', ...
        'RegionLocationMaxX1','RegionLocationMaxY1'};
    if ~all(isfield(region,required)), continue; end
    item.x0 = double(region.RegionLocationMinX0);
    item.y0 = double(region.RegionLocationMinY0);
    item.x1 = double(region.RegionLocationMaxX1);
    item.y1 = double(region.RegionLocationMaxY1);
    item.area = max(0,item.x1-item.x0+1)*max(0,item.y1-item.y0+1);
    item.cx = (item.x0+item.x1)/2;
    regions(end+1) = item; %#ok<AGROW>
end
if numel(regions) < 3, return; end

imageMidX = double(info.Columns)/2;
leftCandidates = find([regions.cx] < imageMidX);
rightCandidates = find([regions.cx] >= imageMidX);
if isempty(leftCandidates) || isempty(rightCandidates), return; end
[~,leftRelative] = max([regions(leftCandidates).area]);
leftMain = regions(leftCandidates(leftRelative));
[~,rightRelative] = max([regions(rightCandidates).area]);
rightMain = regions(rightCandidates(rightRelative));

smallCandidates = rightCandidates([regions(rightCandidates).area] < 0.50*rightMain.area);
inside = false(size(smallCandidates));
for i = 1:numel(smallCandidates)
    candidate = regions(smallCandidates(i));
    inside(i) = candidate.x0 >= rightMain.x0 && candidate.x1 <= rightMain.x1 && ...
        candidate.y0 >= rightMain.y0 && candidate.y1 <= rightMain.y1;
end
smallCandidates = smallCandidates(inside);
if isempty(smallCandidates), return; end
[~,largestSmall] = max([regions(smallCandidates).area]);
sweRegion = regions(smallCandidates(largestSmall));

% Transfer the right-panel ROI's relative coordinates to the left panel.
leftFullX0 = leftMain.x0 + (sweRegion.x0-rightMain.x0);
leftFullY0 = leftMain.y0 + (sweRegion.y0-rightMain.y0);
localX = leftFullX0-cropOriginXY(1)+1;
localY = leftFullY0-cropOriginXY(2)+1;
boxWidth = sweRegion.x1-sweRegion.x0+1;
boxHeight = sweRegion.y1-sweRegion.y0+1;
if localX < 1 || localY < 1 || localX+boxWidth-1 > w || localY+boxHeight-1 > h
    return;
end
box = [localX localY boxWidth boxHeight];
found = true;
end

function mask = rectangle_border_mask(box,h,w,pad)
x1 = max(1,floor(box(1))-pad); y1 = max(1,floor(box(2))-pad);
x2 = min(w,ceil(box(1)+box(3)-1)+pad); y2 = min(h,ceil(box(2)+box(4)-1)+pad);
mask = false(h,w); band = max(3,pad);
mask(y1:min(y1+band,y2),x1:x2) = true;
mask(max(y2-band,y1):y2,x1:x2) = true;
mask(y1:y2,x1:min(x1+band,x2)) = true;
mask(y1:y2,max(x2-band,x1):x2) = true;
end
