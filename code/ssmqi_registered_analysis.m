function out = ssmqi_registered_analysis(T,outDir,varargin)
% SSMQI_REGISTERED_ANALYSIS  Runs the seven analyses preregistered at https://osf.io/gfd6t.
%
%   out = ssmqi_registered_analysis(T,outDir,'nBoot',2000,'seed',20261024)
%
% T comes from ssmqi_load_registered_measurements. Optional name-value pair
% 'flagged' is a table with Rater, Session and Source_Image columns listing tracings
% to leave out of analyses 1 to 6 (sensitivity analysis only); traceability always
% uses the original accepted or untraceable status.
%
% Comparisons, planes, statistics and confidence intervals follow the registration:
% long-axis and short-axis images are analyzed separately; image pairs where either
% tracing is untraceable are left out of analyses 1 to 6; confidence intervals come
% from a percentile bootstrap that resamples whole participants.

p = inputParser;
addParameter(p,'nBoot',2000);
addParameter(p,'seed',20261024);
addParameter(p,'flagged',table());
parse(p,varargin{:});
nBoot = p.Results.nBoot;
seed = p.Results.seed;
flagged = p.Results.flagged;

core = ssmqi_stats_core();
if ~isfolder(outDir), mkdir(outDir); end
rng(seed,'twister');

T.Traceable = T.Complete;
T.Usable = T.Complete;
if ~isempty(flagged) && height(flagged) > 0
    for i = 1:height(flagged)
        m = T.Rater == flagged.Rater(i) & T.Session == flagged.Session(i) & ...
            T.Source_Image == string(flagged.Source_Image(i));
        T.Usable(m) = false;
    end
end

comps = {
    'Intra-rater, Rater 1 (S1 vs S2)',          [1 1], [1 2];
    'Intra-rater, Rater 2 (S1 vs S2)',          [2 1], [2 2];
    'Inter-rater, primary (R1 S1 vs R2 S1)',    [1 1], [2 1];
    'Inter-rater, secondary (R1 S2 vs R2 S2)',  [1 2], [2 2];
    'Sensitivity (R1 S2 vs R2 S1)',             [1 2], [2 1]};
planes = {'Long axis','Short axis'};

outline = table(); ei = table(); trace = table(); pairData = table();
for c = 1:size(comps,1)
    for pl = 1:numel(planes)
        [A,B] = join_pairs(T,comps{c,2},comps{c,3},planes{pl});
        if isempty(A), continue; end
        label = comps{c,1}; plane = planes{pl};

        % Analysis 7: traceability (all images present in both sessions)
        if height(A) >= 2
            part = A.Participant;
            kp = @(ix) core.kappaPabak(A.Traceable(ix),B.Traceable(ix));
            r0 = kp((1:height(A))');
            [~,ciPa] = core.clusterBoot(part,@(ix) pick(kp(ix),'percentAgreement'),nBoot);
            [~,ciKa] = core.clusterBoot(part,@(ix) pick(kp(ix),'kappa'),nBoot);
            [~,ciPb] = core.clusterBoot(part,@(ix) pick(kp(ix),'pabak'),nBoot);
            trace = [trace; table(string(label),string(plane),height(A),r0.percentAgreement,ciPa(1),ciPa(2), ...
                r0.kappa,ciKa(1),ciKa(2),r0.pabak,ciPb(1),ciPb(2), ...
                'VariableNames',{'Comparison','Plane','Images','PercentAgreement','PA_CI_Low','PA_CI_High', ...
                'Kappa','Kappa_CI_Low','Kappa_CI_High','PABAK','PABAK_CI_Low','PABAK_CI_High'})]; %#ok<AGROW>
        end

        both = A.Usable & B.Usable;
        A = A(both,:); B = B(both,:);
        n = height(A);
        if n < 3, continue; end
        part = A.Participant;
        dice = NaN(n,1); jac = NaN(n,1); hd = NaN(n,1); msd = NaN(n,1);
        for i = 1:n
            [dice(i),jac(i)] = core.diceJaccard(A.Mask{i},B.Mask{i});
            [hd(i),msd(i)] = core.boundaryDist(A.Mask{i},B.Mask{i},A.mmPerPx(i));
        end
        x1 = A.EI; x2 = B.EI;
        pairData = [pairData; table(repmat(string(label),n,1),repmat(string(plane),n,1),A.Source_Image, ...
            part,dice,jac,hd,msd,x1,x2,'VariableNames',{'Comparison','Plane','Source_Image', ...
            'Participant','Dice','Jaccard','HD95_mm','MSD_mm','EI_First','EI_Second'})]; %#ok<AGROW>

        % Analyses 1 and 2: overlap and boundary distance
        row = {string(label),string(plane),n};
        for v = {dice,jac,hd,msd}
            vec = v{1};
            [est,ci] = core.clusterBoot(part,@(ix) mean(vec(ix),'omitnan'),nBoot);
            row = [row, {est,ci(1),ci(2),median(vec,'omitnan')}]; %#ok<AGROW>
        end
        outline = [outline; cell2table(row,'VariableNames',{'Comparison','Plane','Pairs', ...
            'Dice_Mean','Dice_CI_Low','Dice_CI_High','Dice_Median', ...
            'Jaccard_Mean','Jaccard_CI_Low','Jaccard_CI_High','Jaccard_Median', ...
            'HD95_mm_Mean','HD95_mm_CI_Low','HD95_mm_CI_High','HD95_mm_Median', ...
            'MSD_mm_Mean','MSD_mm_CI_Low','MSD_mm_CI_High','MSD_mm_Median'})]; %#ok<AGROW>

        % Analyses 3 and 4: Bland-Altman of echo intensity, SEM, MDC95, CV
        [d,m] = core.baPercent(x1,x2);
        bias = @(ix) mean(d(ix));
        lo = @(ix) mean(d(ix)) - 1.96*std(d(ix));
        hi = @(ix) mean(d(ix)) + 1.96*std(d(ix));
        [b0,cb] = core.clusterBoot(part,bias,nBoot);
        [l0,cl] = core.clusterBoot(part,lo,nBoot);
        [h0,ch] = core.clusterBoot(part,hi,nBoot);
        lm = fitlm(m,d);
        s0 = core.agreementStats(x1,x2);
        [~,cs] = core.clusterBoot(part,@(ix) getfield(core.agreementStats(x1(ix),x2(ix)),'sem'),nBoot);
        [~,cm] = core.clusterBoot(part,@(ix) getfield(core.agreementStats(x1(ix),x2(ix)),'mdc95'),nBoot);
        [~,cv] = core.clusterBoot(part,@(ix) getfield(core.agreementStats(x1(ix),x2(ix)),'cvPercent'),nBoot);
        ei = [ei; table(string(label),string(plane),n,b0,cb(1),cb(2),l0,cl(1),cl(2),h0,ch(1),ch(2), ...
            lm.Coefficients.Estimate(2),lm.Coefficients.pValue(2), ...
            s0.sem,cs(1),cs(2),s0.mdc95,cm(1),cm(2),s0.cvPercent,cv(1),cv(2), ...
            'VariableNames',{'Comparison','Plane','Pairs','Bias_Percent','Bias_CI_Low','Bias_CI_High', ...
            'LoA_Lower','LoA_Lower_CI_Low','LoA_Lower_CI_High','LoA_Upper','LoA_Upper_CI_Low','LoA_Upper_CI_High', ...
            'ProportionalBias_Slope','ProportionalBias_p', ...
            'SEM','SEM_CI_Low','SEM_CI_High','MDC95','MDC95_CI_Low','MDC95_CI_High', ...
            'CV_Percent','CV_CI_Low','CV_CI_High'})]; %#ok<AGROW>

        make_ba_plot(m,d,b0,l0,h0,sprintf('%s, %s',label,plane),outDir);
    end
end

% Analyses 5 and 6: mixed model, ICC(A,1) and decision study (all tracings, per plane)
varRows = table(); dsRows = table();
for pl = 1:numel(planes)
    R = T(T.Usable & T.Axis == planes{pl},:);
    if numel(unique(R.Rater)) < 2 || height(R) < 20, continue; end
    tbl = table(R.EI,"R"+string(R.Rater),R.Participant,R.Source_Image, ...
        'VariableNames',{'EI','Rater','Participant','Image'});
    comp = core.fitComponents(tbl);
    boots = bootstrap_model(core,tbl,nBoot);
    cols = {'sigma2_participant','sigma2_image','sigma2_rater','sigma2_residual','icc'};
    row = {string(planes{pl}),height(R),numel(unique(R.Participant))};
    for k = 1:numel(cols)
        ci = prctile(boots.(cols{k})(~isnan(boots.(cols{k}))),[2.5 97.5]);
        row = [row,{comp.(cols{k}),ci(1),ci(2)}]; %#ok<AGROW>
    end
    varRows = [varRows; cell2table(row,'VariableNames',{'Plane','Tracings','Participants', ...
        'Participant_Var','Participant_Var_CI_Low','Participant_Var_CI_High', ...
        'Image_Var','Image_Var_CI_Low','Image_Var_CI_High', ...
        'Rater_Var','Rater_Var_CI_Low','Rater_Var_CI_High', ...
        'Residual_Var','Residual_Var_CI_Low','Residual_Var_CI_High', ...
        'ICC_A1','ICC_CI_Low','ICC_CI_High'})]; %#ok<AGROW>
    specs = {'1 tracing',1,1; '2 tracings, same rater',2,1; '3 tracings, same rater',3,1; '2 raters, 1 tracing each',1,2};
    for s = 1:size(specs,1)
        ds = core.dstudy(comp,specs{s,2},specs{s,3});
        if specs{s,3} == 1, est = ds.sameRater; else, est = ds.multiRater; end
        bs = NaN(numel(boots.icc),1);
        for b = 1:numel(bs)
            cb = struct('sigma2_participant',boots.sigma2_participant(b),'sigma2_image',boots.sigma2_image(b), ...
                'sigma2_rater',boots.sigma2_rater(b),'sigma2_residual',boots.sigma2_residual(b));
            dsb = core.dstudy(cb,specs{s,2},specs{s,3});
            if specs{s,3} == 1, bs(b) = dsb.sameRater; else, bs(b) = dsb.multiRater; end
        end
        ci = prctile(bs(~isnan(bs)),[2.5 97.5]);
        dsRows = [dsRows; table(string(planes{pl}),string(specs{s,1}),est,ci(1),ci(2), ...
            'VariableNames',{'Plane','Design','Projected_ICC','CI_Low','CI_High'})]; %#ok<AGROW>
    end
end

writetable(outline,fullfile(outDir,'Table2_Outline_Agreement.csv'));
writetable(ei,fullfile(outDir,'Table3_EchoIntensity_Agreement.csv'));
writetable(varRows,fullfile(outDir,'Table4_Variance_Components_ICC.csv'));
writetable(dsRows,fullfile(outDir,'Table5_Decision_Study.csv'));
writetable(trace,fullfile(outDir,'Table6_Traceability_Agreement.csv'));
writetable(pairData,fullfile(outDir,'Pair_Level_Metrics.csv'));
out = struct('outline',outline,'echoIntensity',ei,'variance',varRows,'decisionStudy',dsRows, ...
    'traceability',trace,'pairs',pairData);
end

function [A,B] = join_pairs(T,a,b,plane)
A = T(T.Rater==a(1) & T.Session==a(2) & T.Axis==plane,:);
B = T(T.Rater==b(1) & T.Session==b(2) & T.Axis==plane,:);
if isempty(A) || isempty(B), A = A([],:); B = B([],:); return; end
[~,ia,ib] = intersect(A.Source_Image,B.Source_Image,'stable');
A = A(ia,:); B = B(ib,:);
end

function v = pick(s,name)
v = s.(name);
end

function boots = bootstrap_model(core,tbl,nBoot)
parts = unique(tbl.Participant);
f = {'sigma2_participant','sigma2_image','sigma2_rater','sigma2_residual','icc'};
for k = 1:numel(f), boots.(f{k}) = NaN(nBoot,1); end
rows = arrayfun(@(q) find(tbl.Participant == parts(q)),1:numel(parts),'UniformOutput',false);
n = numel(parts);
warnState = warning('off','all');
cleanup = onCleanup(@() warning(warnState));
for b = 1:nBoot
    pick = randi(n,n,1);
    pieces = cell(n,1);
    for q = 1:n
        piece = tbl(rows{pick(q)},:);
        piece.Participant = repmat("b"+string(q),height(piece),1);
        piece.Image = piece.Image + "_b" + string(q);
        pieces{q} = piece;
    end
    try
        comp = core.fitComponents(vertcat(pieces{:}));
        for k = 1:numel(f), boots.(f{k})(b) = comp.(f{k}); end
    catch
    end
end
end

function make_ba_plot(m,d,bias,lo,hi,titleText,outDir)
name = regexprep(titleText,'[^A-Za-z0-9]+','_');
target = fullfile(outDir,['BlandAltman_' name '.png']);
% Headless MATLAB can invalidate an invisible figure before export. Retry with a fresh
% figure; if the plot still cannot be written, warn and continue (the plot is a display
% aid only and no registered statistic depends on it).
for attempt = 1:3
    fig = [];
    try
        fig = figure('Visible','off','Color','w','Position',[100 100 700 500]);
        ax = axes(fig);
        scatter(ax,m,d,28,'filled','MarkerFaceAlpha',0.6); hold(ax,'on');
        yline(ax,bias,'-','Bias'); yline(ax,lo,'--'); yline(ax,hi,'--');
        xlabel(ax,'Mean echo intensity of the two tracings (0 to 255)');
        ylabel(ax,'Difference (% of pair mean)');
        title(ax,titleText,'Interpreter','none');
        drawnow;
        exportgraphics(fig,target,'Resolution',200);
        if isgraphics(fig), close(fig); end
        return
    catch err
        if ~isempty(fig) && isgraphics(fig), close(fig); end
        if attempt == 3
            warning('ssmqi:plotFailed','Bland-Altman plot not written (%s): %s',name,err.message);
        end
    end
end
end
