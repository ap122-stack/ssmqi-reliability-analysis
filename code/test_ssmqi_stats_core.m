function test_ssmqi_stats_core()
% TEST_SSMQI_STATS_CORE  Checks each statistic against known answers on synthetic data.
% Run before any real analysis. Uses no study data.

core = ssmqi_stats_core();
rng(20261024,'twister');
results = {};
function check(name,ok,detail)
    results(end+1,:) = {name,ok,detail}; %#ok<AGROW>
end

% Dice and Jaccard
a = false(100); a(11:20,11:30) = true;
[d,j] = core.diceJaccard(a,a);
check('Dice identical = 1',abs(d-1)<1e-12,sprintf('%.4f',d));
b = false(100); b(11:20,21:40) = true;
[d,j] = core.diceJaccard(a,b);
check('Dice half overlap = 0.5',abs(d-0.5)<1e-12,sprintf('%.4f',d));
check('Jaccard half overlap = 1/3',abs(j-1/3)<1e-12,sprintf('%.4f',j));
c = false(100); c(60:80,60:80) = true;
d = core.diceJaccard(a,c);
check('Dice disjoint = 0',d==0,sprintf('%.4f',d));

% Boundary distance on concentric squares: boundaries are 5 px apart
A = false(120); A(41:80,41:80) = true;
Bm = false(120); Bm(36:85,36:85) = true;
[hd,msd] = core.boundaryDist(A,A,0.1);
check('HD95 identical = 0',hd==0 && msd==0,sprintf('%.3f %.3f',hd,msd));
[hd,msd] = core.boundaryDist(A,Bm,0.1);
% Directed distances from the outer square to the inner one are 5 px along the edges and 5.1 to 7.07 px
% within 4 px of each corner. Corners make up about 18% of that perimeter, so the 95th percentile is 6.40 px.
check('HD95 concentric squares = 0.64 mm',abs(hd-0.64)<0.01,sprintf('%.3f mm',hd));
check('MSD concentric squares about 0.5 mm',abs(msd-0.5)<0.02,sprintf('%.3f mm',msd));

% Bland-Altman percent difference
[dp,m] = core.baPercent(100,110);
check('Percent difference 100 vs 110',abs(dp-100*10/105)<1e-9 && abs(m-105)<1e-9,sprintf('%.4f',dp));

% SEM, MDC, CV
s = core.agreementStats([10 20 30],[11 19 31]);
check('Within-subject SD = 0.7071',abs(s.wsd-sqrt(0.5))<1e-9,sprintf('%.4f',s.wsd));
check('MDC95 = 1.96*sqrt(2)*SEM',abs(s.mdc95-1.96*sqrt(2)*sqrt(0.5))<1e-9,sprintf('%.4f',s.mdc95));

% Kappa and PABAK
r = core.kappaPabak([1 1 0 0],[1 0 0 0]);
check('Kappa = 0.5',abs(r.kappa-0.5)<1e-12,sprintf('%.4f',r.kappa));
check('PABAK = 0.5',abs(r.pabak-0.5)<1e-12,sprintf('%.4f',r.pabak));
check('Percent agreement = 75',abs(r.percentAgreement-75)<1e-12,sprintf('%.1f',r.percentAgreement));

% Variance components recovered from simulated data with known truth
nP = 200; nImgPer = 3;
sp2 = 100; si2 = 50; delta = 4; se2 = 25;
rows = {};
for p = 1:nP
    pe = sqrt(sp2)*randn;
    for k = 1:nImgPer
        ie = sqrt(si2)*randn;
        img = sprintf('P%d_I%d',p,k);
        base = 120 + pe + ie;
        rows(end+1,:) = {base - delta/2 + sqrt(se2)*randn,'R1',sprintf('P%d',p),img}; %#ok<AGROW>
        rows(end+1,:) = {base - delta/2 + sqrt(se2)*randn,'R1',sprintf('P%d',p),img}; %#ok<AGROW>
        rows(end+1,:) = {base + delta/2 + sqrt(se2)*randn,'R2',sprintf('P%d',p),img}; %#ok<AGROW>
    end
end
tbl = cell2table(rows,'VariableNames',{'EI','Rater','Participant','Image'});
comp = core.fitComponents(tbl);
trueIcc = (sp2+si2)/(sp2+si2+delta^2/2+se2);
check('Participant variance near 100',abs(comp.sigma2_participant-sp2)/sp2<0.25,sprintf('%.1f',comp.sigma2_participant));
check('Image variance near 50',abs(comp.sigma2_image-si2)/si2<0.25,sprintf('%.1f',comp.sigma2_image));
check('Residual variance near 25',abs(comp.sigma2_residual-se2)/se2<0.15,sprintf('%.1f',comp.sigma2_residual));
check('Rater variance near delta^2/2 = 8',abs(comp.sigma2_rater-8)<2.5,sprintf('%.2f',comp.sigma2_rater));
check(sprintf('ICC(A,1) near truth %.3f',trueIcc),abs(comp.icc-trueIcc)<0.03,sprintf('%.3f',comp.icc));

% Decision study
ds = core.dstudy(comp,2,2);
manual = (comp.sigma2_participant+comp.sigma2_image) / ...
    (comp.sigma2_participant+comp.sigma2_image+comp.sigma2_rater+comp.sigma2_residual/2);
check('D-study averaging 2 tracings matches formula',abs(ds.sameRater-manual)<1e-12,sprintf('%.3f',ds.sameRater));
check('Averaging more tracings never lowers ICC',ds.sameRater >= comp.icc && ds.multiRater >= comp.icc,sprintf('%.3f %.3f',ds.sameRater,ds.multiRater));

% Cluster bootstrap: interval should bracket the estimate and be wider with more spread
cl = repelem((1:30)',3,1);
x = randn(90,1) + repelem(randn(30,1)*2,3,1);
[est,ci] = core.clusterBoot(cl,@(ix) mean(x(ix)),500);
check('Bootstrap CI brackets estimate',ci(1)<=est && est<=ci(2),sprintf('%.2f [%.2f %.2f]',est,ci(1),ci(2)));

% Report
ok = cell2mat(results(:,2));
for i = 1:size(results,1)
    fprintf('%s  %-52s %s\n',ternary(ok(i),'PASS','FAIL'),results{i,1},results{i,3});
end
fprintf('\n%d of %d checks passed.\n',sum(ok),numel(ok));
if ~all(ok), error('Some checks failed.'); end
end

function s = ternary(c,a,b)
if c, s = a; else, s = b; end
end
