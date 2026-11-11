function core = ssmqi_stats_core()
% SSMQI_STATS_CORE  Statistical building blocks for the registered analysis.
%
% Returns a struct of function handles so each calculation can be tested on
% its own. No file input or output happens here.
%
%   core.diceJaccard(a,b)            overlap of two logical masks
%   core.boundaryDist(a,b,mmPerPx)   HD95 and mean surface distance (mm)
%   core.baPercent(x1,x2)            percent difference and pair mean
%   core.agreementStats(x1,x2)       SEM, MDC95 and within-subject CV
%   core.fitComponents(tbl)          mixed-model variance components and ICC(A,1)
%   core.dstudy(comp,k,nRaters)      projected ICC when tracings or raters are averaged
%   core.kappaPabak(s1,s2)           percent agreement, Cohen kappa and PABAK
%   core.clusterBoot(cluster,fun,B)  percentile bootstrap resampling whole clusters

core.diceJaccard = @dice_jaccard;
core.boundaryDist = @boundary_dist_mm;
core.baPercent = @ba_percent;
core.agreementStats = @agreement_stats;
core.fitComponents = @fit_components;
core.dstudy = @dstudy;
core.kappaPabak = @kappa_pabak;
core.clusterBoot = @cluster_boot;
end

function [dice,jaccard] = dice_jaccard(a,b)
a = logical(a); b = logical(b);
inter = nnz(a & b);
sa = nnz(a); sb = nnz(b);
if sa + sb == 0
    dice = NaN; jaccard = NaN; return;
end
dice = 2*inter/(sa+sb);
jaccard = inter/nnz(a | b);
end

function [hd95,msd] = boundary_dist_mm(a,b,mmPerPx)
a = logical(a); b = logical(b);
if ~any(a(:)) || ~any(b(:))
    hd95 = NaN; msd = NaN; return;
end
pa = bwperim(a); pb = bwperim(b);
dToB = bwdist(pb);
dToA = bwdist(pa);
dab = double(dToB(pa)) * mmPerPx;
dba = double(dToA(pb)) * mmPerPx;
hd95 = max(prctile(dab,95), prctile(dba,95));
msd = (mean(dab) + mean(dba)) / 2;
end

function [d,m] = ba_percent(x1,x2)
x1 = x1(:); x2 = x2(:);
m = (x1 + x2) / 2;
d = 100 * (x2 - x1) ./ m;
end

function s = agreement_stats(x1,x2)
x1 = x1(:); x2 = x2(:);
d = x1 - x2;
m = (x1 + x2) / 2;
s.wsd = sqrt(mean(d.^2) / 2);
s.sem = s.wsd;
s.mdc95 = 1.96 * sqrt(2) * s.sem;
s.cvPercent = 100 * sqrt(mean((d ./ m).^2 / 2));
end

function comp = fit_components(tbl)
% tbl variables: EI (double), Rater (categorical, 2 levels), Participant, Image
% Model: EI ~ Rater + (1|Participant) + (1|Image), REML.
% sigma2_rater is the variance of the two rater means, delta^2 / 2, where delta is
% the fitted rater contrast (the (k-1)-divisor estimate used by the ICC(A,1) ANOVA formula).
tbl.Rater = categorical(tbl.Rater);
tbl.Participant = categorical(tbl.Participant);
tbl.Image = categorical(tbl.Image);
lme = fitlme(tbl,'EI ~ Rater + (1|Participant) + (1|Image)','FitMethod','REML');
[psi,mse,stats] = covarianceParameters(lme);
comp.sigma2_participant = NaN; comp.sigma2_image = NaN;
for i = 1:numel(psi)
    g = strtrim(char(stats{i}.Group(1,:)));
    if strcmp(g,'Participant'), comp.sigma2_participant = psi{i}(1,1); end
    if strcmp(g,'Image'), comp.sigma2_image = psi{i}(1,1); end
end
comp.sigma2_residual = mse;
nm = lme.Coefficients.Name;
k = find(startsWith(nm,'Rater_'),1);
if isempty(k), delta = 0; else, delta = lme.Coefficients.Estimate(k); end
comp.delta = delta;
comp.sigma2_rater = delta^2 / 2;
num = comp.sigma2_participant + comp.sigma2_image;
comp.icc = num / (num + comp.sigma2_rater + comp.sigma2_residual);
end

function r = dstudy(comp,k,nRaters)
% Projected ICC(A) when k tracings by the same rater are averaged (rater effect unchanged)
% and when nRaters raters each contribute one tracing (rater and residual variance shrink).
num = comp.sigma2_participant + comp.sigma2_image;
r.sameRater = num / (num + comp.sigma2_rater + comp.sigma2_residual / k);
r.multiRater = num / (num + comp.sigma2_rater / nRaters + comp.sigma2_residual / nRaters);
end

function r = kappa_pabak(s1,s2)
s1 = logical(s1(:)); s2 = logical(s2(:));
po = mean(s1 == s2);
p1 = mean(s1); p2 = mean(s2);
pe = p1*p2 + (1-p1)*(1-p2);
if pe >= 1, r.kappa = NaN; else, r.kappa = (po - pe) / (1 - pe); end
r.percentAgreement = 100 * po;
r.pabak = 2*po - 1;
end

function [est,ci,boots] = cluster_boot(cluster,fun,B)
% Resample whole clusters (participants) with replacement B times.
% fun(rowIndices) returns a scalar statistic computed on those rows.
cluster = cluster(:);
[uc,~,g] = unique(cluster);
rows = accumarray(g,(1:numel(g))',[],@(x){x});
est = fun((1:numel(g))');
boots = NaN(B,1);
n = numel(uc);
for b = 1:B
    pick = randi(n,n,1);
    idx = vertcat(rows{pick});
    try, boots(b) = fun(idx); catch, boots(b) = NaN; end
end
ci = prctile(boots(~isnan(boots)),[2.5 97.5]);
end
