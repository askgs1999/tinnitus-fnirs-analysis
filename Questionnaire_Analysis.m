%  QUESTIONNAIRE ANALYSIS
%  THI and TFI scores, and their correlation with cortical measures
%
%  Satish G, Taylor LM, Arnold MP, Gallagher-Shale J, Krishnamurthy K,
%  Basura GJ. Subjective tinnitus distress correlates with cortical
%  hemodynamic changes in human auditory cortex as measured by functional
%  near-infrared spectroscopy. PLOS ONE.
%
%  WHAT THIS SCRIPT DOES
%    1. Reads the REDCap export, separates pre- and post-session records,
%       and scores the Tinnitus Handicap Inventory and the Tinnitus
%       Functional Index for each participant.
%    2. Tests the within-session change in each instrument.
%    3. Generates Figure 3, the participant-wise pre-to-post change plots.
%    4. Correlates both instruments, as absolute and change scores, with
%       the hemodynamic and connectivity measures produced by the other two
%       scripts, applying Bonferroni correction within each prespecified
%       family of tests.
%    5. Generates Figure 10, the association between within-session change
%       in TFI and late-window HbO during interstimulus rest.
%
%  TEST SELECTION
%    TFI scores were normally distributed and are analysed with Pearson
%    correlations. THI scores were not, and are analysed with Spearman
%    correlations. This rule is applied consistently throughout.
%
%  REQUIREMENTS
%    MATLAB R2024b
%    Statistics and Machine Learning Toolbox
%
%  INPUT   all three files must sit in the same folder
%    SubjectiveAndSomatic_DATA_2026-01-15_1203.csv   REDCap questionnaire export
%    RSFC_outputs.mat    written by rsfc_analysis.m
%    HR_outputs.mat      written by hemodynamic_analysis.m
%
%  OUTPUT
%    TINN_HbO_THI_TFI_merged.mat and .csv   merged participant-level table
%    Fig3.tif Fig10_dTFI_vs_lateISR.tif     written to FIG_DIR
%
%  HOW TO RUN
%    Run rsfc_analysis.m and hemodynamic_analysis.m first, then set
%    MATLAB's current folder to the one holding the three input files and
%    run this script.

clearvars; close all; clc;

%% ---- PATHS ----
% The REDCap export is located automatically using WHICH, which searches the
% current folder and the MATLAB path. To set the folder by hand, assign
% DATA_DIR directly and delete the block below.

INPUT_FILE = 'SubjectiveAndSomatic_DATA_2026-01-15_1203.csv';

found = which(INPUT_FILE);
if isempty(found) && isfile(fullfile(pwd, INPUT_FILE))
    found = fullfile(pwd, INPUT_FILE);
end
if isempty(found)
    error(['Cannot find %s.\n' ...
           'Current folder is %s.\n' ...
           'Either cd to the folder holding the data, add it to the MATLAB ' ...
           'path, or set DATA_DIR by hand near the top of this script.'], ...
           INPUT_FILE, pwd);
end

DATA_DIR    = fileparts(found);
FIG_DIR     = fullfile(DATA_DIR, 'figures');
RESULTS_DIR = fullfile(DATA_DIR, 'results');

if ~exist(FIG_DIR,'dir'),     mkdir(FIG_DIR);     end
if ~exist(RESULTS_DIR,'dir'), mkdir(RESULTS_DIR); end

for f = {'RSFC_outputs.mat','HR_outputs.mat'}
    if ~isfile(fullfile(DATA_DIR, f{1}))
        error(['%s not found in %s.\n' ...
               'Run rsfc_analysis.m and hemodynamic_analysis.m first.'], ...
               f{1}, DATA_DIR);
    end
end

fprintf('Reading data from : %s\n', DATA_DIR);
fprintf('Figures to        : %s\n', FIG_DIR);


%%
clc
close all;



% 1) LOAD NEW CSV

csvPath = fullfile(DATA_DIR, 'SubjectiveAndSomatic_DATA_2026-01-15_1203.csv');
data = readtable(csvPath);


% REMOVE EXCLUDED PARTICIPANT(S)

excludeIDs = ["Tinnitus11","Tinnitus71","Tinnitus72","Tinnitus37","Tinnitus55","Tinnitus33"]; % removing the subjects excluded from analysis
rmIdx = ismember(string(data.record_id), excludeIDs);
data(rmIdx,:) = [];


% 2) SPLIT PRE vs POST ROWS

isPre  = strcmpi(string(data.redcap_event_name), "pre_test_survey_arm_1");
isPost = strcmpi(string(data.redcap_event_name), "post_test_survey_arm_1");

preTbl  = data(isPre,  :);
postTbl = data(isPost, :);

[ids, iPre, iPost] = intersect(preTbl.record_id, postTbl.record_id, 'stable');
preTbl  = preTbl(iPre, :);
postTbl = postTbl(iPost,:);

% 3) COMPUTE THI totals

thi_pre_cols  = startsWith(preTbl.Properties.VariableNames, "thi_") & ...
                ~contains(preTbl.Properties.VariableNames, "_v2") & ...
                ~contains(preTbl.Properties.VariableNames, "complete") & ...
                ~contains(preTbl.Properties.VariableNames, "score");

thi_post_cols = startsWith(postTbl.Properties.VariableNames, "thi_") & ...
                contains(postTbl.Properties.VariableNames, "_v2") & ...
                ~contains(postTbl.Properties.VariableNames, "complete") & ...
                ~contains(postTbl.Properties.VariableNames, "score");

thi_pre_items  = preTbl{:, thi_pre_cols};
thi_post_items = postTbl{:, thi_post_cols};

THI_pre_all  = sum(thi_pre_items,  2, 'omitnan');
THI_post_all = sum(thi_post_items, 2, 'omitnan');
dTHI_all     = THI_post_all - THI_pre_all;


% 4) COMPUTE TFI totals

tfi_pre_cols  = ~contains(preTbl.Properties.VariableNames, "_v2") & ...
                startsWith(preTbl.Properties.VariableNames, "tfi_") & ...
                ~contains(preTbl.Properties.VariableNames, "complete");

tfi_post_cols = contains(postTbl.Properties.VariableNames, "_v2") & ...
                startsWith(postTbl.Properties.VariableNames, "tfi_") & ...
                ~contains(postTbl.Properties.VariableNames, "complete");

tfi_pre_items  = preTbl{:, tfi_pre_cols};
tfi_post_items = postTbl{:, tfi_post_cols};

max_per_item = 10;
num_items_pre  = sum(~isnan(tfi_pre_items),  2);
num_items_post = sum(~isnan(tfi_post_items), 2);
sum_pre  = sum(tfi_pre_items,  2, 'omitnan');
sum_post = sum(tfi_post_items, 2, 'omitnan');

TFI_pre_all  = (sum_pre  ./ (num_items_pre  * max_per_item)) * 100;
TFI_post_all = (sum_post ./ (num_items_post * max_per_item)) * 100;
dTFI_all     = TFI_post_all - TFI_pre_all;


% 5) KEEP VALID SUBJECTS, drop (dTHI=0 & dTFI=0)

validBoth = ~isnan(THI_pre_all) & ~isnan(THI_post_all) & ...
            ~isnan(TFI_pre_all) & ~isnan(TFI_post_all);
zeroBoth  = (dTHI_all == 0) & (abs(dTFI_all) < 1e-12);
keepIdx   = validBoth & ~zeroBoth;

fprintf('\nTotal paired IDs: %d\n', numel(ids));
fprintf('Valid THI+TFI pairs: %d\n', sum(validBoth));
fprintf('Removed (dTHI=0 AND dTFI=0): %d\n', sum(validBoth & zeroBoth));
fprintf('Remaining: %d\n', sum(keepIdx));

removedIDs = ids(validBoth & zeroBoth);
if ~isempty(removedIDs)
    disp('Removed IDs (dTHI=0 AND dTFI=0):'); disp(removedIDs(:));
end

ids_final = ids(keepIdx);
THI_pre  = THI_pre_all(keepIdx);  THI_post = THI_post_all(keepIdx);
dTHI     = dTHI_all(keepIdx);     nTHI = numel(THI_pre);
TFI_pre  = TFI_pre_all(keepIdx);  TFI_post = TFI_post_all(keepIdx);
dTFI     = dTFI_all(keepIdx);     nTFI = numel(TFI_pre);


% 6) PLOT THI

figure('Color','w'); hold on;
y_ranges = [0 16; 18 36; 38 56; 58 76; 78 100];
colors = [0.85 0.93 1.00; 0.90 1.00 0.90; 1.00 0.95 0.80; 1.00 0.88 0.88; 1.00 0.80 0.80];
for i = 1:size(y_ranges,1)
    rectangle('Position',[0.7 y_ranges(i,1) 1.6 diff(y_ranges(i,:))], ...
              'FaceColor', colors(i,:), 'EdgeColor','none');
end
for i = 1:nTHI
    if dTHI(i) < 0, col = [0.1 0.6 0.1]; else, col = [0.8 0.1 0.1]; end
    plot([1 2], [THI_pre(i) THI_post(i)], '-', 'Color', col, 'LineWidth', 2);
    plot(1, THI_pre(i),  'o', 'MarkerSize', 8, 'MarkerFaceColor', col, 'MarkerEdgeColor','k');
    plot(2, THI_post(i), 'o', 'MarkerSize', 8, 'MarkerFaceColor', col, 'MarkerEdgeColor','k');
    text(2.03, THI_post(i), sprintf('D%d', round(dTHI(i))), 'FontWeight','bold','FontSize',10,'Color',col);
end
xlim([0.8 2.3]); ylim([0 100]);
set(gca,'XTick',[1 2],'XTickLabel',{'Pre-test','Post-test'},'FontSize',14,'LineWidth',1.4);
ylabel('THI Score (0-100)','FontSize',16,'FontWeight','bold');
title(sprintf('Participant-wise THI Change (Pre vs Post)  |  N=%d', nTHI),'FontSize',18,'FontWeight','bold');
grid on; box on;


% 7) PLOT TFI

figure('Color','w','Position',[200 100 1400 600]); hold on;
tfi_ranges = [0 17; 18 31; 32 53; 54 72; 73 100];
tfi_colors = [0.85 0.95 1.00; 0.90 1.00 0.90; 1.00 0.98 0.80; 1.00 0.88 0.88; 1.00 0.75 0.75];
for i = 1:size(tfi_ranges,1)
    rectangle('Position',[0.7 tfi_ranges(i,1) 1.6 diff(tfi_ranges(i,:))], ...
              'FaceColor', tfi_colors(i,:), 'EdgeColor','none');
end
for i = 1:nTFI
    if dTFI(i) < 0, col = [0.1 0.6 0.1]; else, col = [0.8 0.1 0.1]; end
    plot([1 2], [TFI_pre(i) TFI_post(i)], '-', 'Color', col, 'LineWidth', 2);
    plot(1, TFI_pre(i),  'o', 'MarkerSize', 8, 'MarkerFaceColor', col, 'MarkerEdgeColor','k');
    plot(2, TFI_post(i), 'o', 'MarkerSize', 8, 'MarkerFaceColor', col, 'MarkerEdgeColor','k');
    text(2.05, TFI_post(i), sprintf('D%.1f', dTFI(i)), 'FontSize',10,'FontWeight','bold','Color',col);
end
xlim([0.8 2.3]); ylim([0 100]);
set(gca,'XTick',[1 2], 'XTickLabel',{'Pre-test','Post-test'},'FontSize',14,'LineWidth',1.4);
ylabel('TFI Score (0-100)', 'FontWeight','bold','FontSize',16);
title(sprintf('Participant-wise TFI Change (Pre vs Post)  |  N=%d', nTFI),'FontSize',18,'FontWeight','bold');
grid on; box on;


%%  THI stats  (paired t-test + Wilcoxon verification)

fprintf('\n================ THI PRE vs POST ================\n');
fprintf('N (paired) = %d\n', nTHI);
mean_pre_thi=mean(THI_pre); mean_post_thi=mean(THI_post);
sd_pre_thi=std(THI_pre);    sd_post_thi=std(THI_post);
[~, p_thi, ci_thi, stats_thi] = ttest(THI_post, THI_pre);
p_thi_w = signrank(THI_post, THI_pre);
d_thi = mean(dTHI)/std(dTHI); mean_change_thi = mean(dTHI);
fprintf('Pre  mean +/- SD = %.2f +/- %.2f\n', mean_pre_thi, sd_pre_thi);
fprintf('Post mean +/- SD = %.2f +/- %.2f\n', mean_post_thi, sd_post_thi);
fprintf('Mean change (Post-Pre) = %.2f\n', mean_change_thi);
fprintf('Paired t-test: t(%d)=%.3f, p=%.6f\n', stats_thi.df, stats_thi.tstat, p_thi);
fprintf('95%% CI of mean change = [%.2f, %.2f]\n', ci_thi(1), ci_thi(2));
fprintf('Cohen''s d_z = %.3f\n', d_thi);
fprintf('Wilcoxon signed-rank p = %.6f\n', p_thi_w);

%% --- TFI stats (paired t-test) ---
fprintf('\n================ TFI PRE vs POST ================\n');
fprintf('N (paired) = %d\n', nTFI);
mean_pre=mean(TFI_pre); mean_post=mean(TFI_post);
sd_pre=std(TFI_pre);    sd_post=std(TFI_post);
[~, p_tfi, ci_tfi, stats_tfi] = ttest(TFI_post, TFI_pre);
p_tfi_w = signrank(TFI_post, TFI_pre);
d_tfi = mean(dTFI)/std(dTFI); mean_change_tfi = mean(dTFI);
fprintf('Pre  mean +/- SD = %.2f +/- %.2f\n', mean_pre, sd_pre);
fprintf('Post mean +/- SD = %.2f +/- %.2f\n', mean_post, sd_post);
fprintf('Mean change (Post-Pre) = %.2f\n', mean_change_tfi);
fprintf('Paired t-test: t(%d)=%.3f, p=%.6f\n', stats_tfi.df, stats_tfi.tstat, p_tfi);
fprintf('95%% CI of mean change = [%.2f, %.2f]\n', ci_tfi(1), ci_tfi(2));
fprintf('Cohen''s d_z = %.3f\n', d_tfi);
fprintf('Wilcoxon signed-rank p = %.6f\n', p_tfi_w);


%% DELTA SCORE TABLE

DeltaTbl = table(string(ids_final(:)), THI_pre(:), THI_post(:), dTHI(:), ...
    TFI_pre(:), TFI_post(:), dTFI(:), ...
    'VariableNames', {'SubjectID','THI_Pre','THI_Post','Delta_THI','TFI_Pre','TFI_Post','Delta_TFI'});
disp(DeltaTbl)


%% MERGE HbO + THI/TFI

% Written by hemodynamic_analysis.m
load(fullfile(DATA_DIR, 'HR_outputs.mat'));   % loads struct HR
T_CorrData_TINN = HR.Tinnitus.ROI;
MergedTbl = innerjoin(T_CorrData_TINN, DeltaTbl, 'Keys', 'SubjectID');
fprintf('\nMerged table: %d subjects with both HbO and THI/TFI data\n', height(MergedTbl));
disp(MergedTbl);
% save(fullfile(DATA_DIR, 'TINN_HbO_THI_TFI_merged.mat'), 'MergedTbl');
% writetable(MergedTbl, fullfile(RESULTS_DIR, 'TINN_HbO_THI_TFI_merged.csv'));
% fprintf('Saved merged file.\n');


%% TEST-TYPE RULE (used everywhere below)
%   THI-based score -> Spearman ; TFI-based score -> Pearson

testForScore = @(scoreName) ternaryTest(contains(lower(scoreName),'thi'));


%% HbO BBN correlations  (THI->Spearman, TFI->Pearson)

thi_avg = (MergedTbl.THI_Pre + MergedTbl.THI_Post)/2;
tfi_avg = (MergedTbl.TFI_Pre + MergedTbl.TFI_Post)/2;

bbn_plot = {
    MergedTbl.early_BBN, MergedTbl.Delta_THI, 'dTHI',     'early', [0.2 0.4 0.8];
    MergedTbl.late_BBN,  MergedTbl.Delta_THI, 'dTHI',     'late',  [0.2 0.4 0.8];
    MergedTbl.late_BBN,  thi_avg,             'mean THI', 'late',  [0.2 0.4 0.8];
    MergedTbl.early_BBN, MergedTbl.Delta_TFI, 'dTFI',     'early', [0.8 0.2 0.2];
    MergedTbl.late_BBN,  MergedTbl.Delta_TFI, 'dTFI',     'late',  [0.8 0.2 0.2];
    MergedTbl.late_BBN,  tfi_avg,             'mean TFI', 'late',  [0.8 0.2 0.2];
};
fprintf('\n================ HbO BBN CORRELATIONS ================\n');
for k = 1:size(bbn_plot,1)
    xall=bbn_plot{k,1}; yall=bbn_plot{k,2}; yLab=bbn_plot{k,3};
    win=bbn_plot{k,4};  col=bbn_plot{k,5};
    tType = testForScore(yLab);
    ok=isfinite(xall)&isfinite(yall); x=xall(ok); y=yall(ok); n=sum(ok);
    [rval,pval]=corr(x,y,'Type',tType);
    if strcmp(tType,'Spearman'), stat='rho'; else, stat='r'; end
    fprintf('\n=== %s vs HbO BBN (ROI, %s) [%s] ===\nN=%d | %s=%.4f | p=%.4g\n', yLab,win,tType,n,stat,rval,pval);
    figure; hold on;
    scatter(x,y,60,'filled','MarkerFaceColor',col,'MarkerEdgeColor','k');
    c=polyfit(x,y,1); xf=linspace(min(x),max(x),100);
    plot(xf,polyval(c,xf),'r-','LineWidth',2);
    xlabel(sprintf('Mean HbO BBN - ROI (%s)',win),'FontSize',13);
    ylabel(yLab,'FontSize',13);
    title(sprintf('%s vs HbO BBN (%s)  |  %s=%.3f, p=%.4g, n=%d',yLab,win,stat,rval,pval,n),'FontSize',13);
    yline(0,'--','Color',[0.5 0.5 0.5]); xline(0,'--','Color',[0.5 0.5 0.5]);
    grid on; box on;
end


%% HbO ISR correlations  (THI->Spearman, TFI->Pearson)

isr_early = MergedTbl.early_ISR;  isr_late = MergedTbl.late_ISR;
vars = {MergedTbl.Delta_THI,'dTHI'; thi_avg,'mean THI'; MergedTbl.Delta_TFI,'dTFI'; tfi_avg,'mean TFI'};
windows = {isr_early,'early'; isr_late,'late'};
isrColor=[0.2 0.6 0.2];
fprintf('\n================ HbO ISR CORRELATIONS ================\n');
for w = 1:2
    x_all=windows{w,1}; winLabel=windows{w,2};
    for v = 1:4
        y_all=vars{v,1}; yLabel=vars{v,2};
        tType=testForScore(yLabel);
        ok=isfinite(x_all)&isfinite(y_all); x=x_all(ok); y=y_all(ok); n=sum(ok);
        [rval,pval]=corr(x,y,'Type',tType);
        if strcmp(tType,'Spearman'), stat='rho'; else, stat='r'; end
        fprintf('\n=== %s vs HbO ISR (ROI, %s) [%s] ===\nN=%d | %s=%.4f | p=%.4g\n',yLabel,winLabel,tType,n,stat,rval,pval);
        figure; hold on;
        scatter(x,y,60,'filled','MarkerFaceColor',isrColor,'MarkerEdgeColor','k');
        c=polyfit(x,y,1); xf=linspace(min(x),max(x),100);
        plot(xf,polyval(c,xf),'r-','LineWidth',2);
        xlabel(sprintf('Mean HbO ISR - ROI (%s)',winLabel),'FontSize',13);
        ylabel(yLabel,'FontSize',13);
        title(sprintf('%s vs HbO ISR (%s)  |  %s=%.3f, p=%.4g, n=%d',yLabel,winLabel,stat,rval,pval,n),'FontSize',13);
        yline(0,'--','Color',[0.5 0.5 0.5]); xline(0,'--','Color',[0.5 0.5 0.5]);
        grid on; box on;
    end
end

%% LOAD RSFC VARIABLES

% Written by rsfc_analysis.m
load(fullfile(DATA_DIR, 'RSFC_outputs.mat'));   % loads struct RSFC

Tinn_SubjectID  = RSFC.SubjectID;
RSFC_ROI_Z_pre  = RSFC.ROI.Z_pre;
RSFC_ROI_Z_post = RSFC.ROI.Z_post;
RSFC_ROI_Z_diff = RSFC.ROI.Z_diff;
RSFC_Lobe_Z_pre  = RSFC.Lobe.Z_pre;
RSFC_Lobe_Z_post = RSFC.Lobe.Z_post;
RSFC_Lobe_Z_diff = RSFC.Lobe.Z_diff;
RSFC_WholeBrain  = RSFC.WholeBrain;
Lobes            = RSFC.Lobes;

fprintf('RSFC outputs loaded\n');
fprintf('   ROI subjects        : %d\n', numel(Tinn_SubjectID));
fprintf('   Whole-brain subjects: %d\n', height(RSFC_WholeBrain));

%% --- Align ROI-seed RSFC with THI/TFI ---
ids_rsfc_roi = string(Tinn_SubjectID(:));
ids_thi      = string(DeltaTbl.SubjectID(:));
[commonIDs_roi, iROI, iTHI_roi] = intersect(ids_rsfc_roi, ids_thi, 'stable');
nROI = numel(commonIDs_roi);
fprintf('\nROI-seed: subjects with BOTH RSFC and THI/TFI = %d\n', nROI);
if nROI == 0, error('No overlapping subjects for ROI-seed.'); end
disp(commonIDs_roi);

rsfc_roi_pre=RSFC_ROI_Z_pre(iROI); rsfc_roi_post=RSFC_ROI_Z_post(iROI); rsfc_roi_diff=RSFC_ROI_Z_diff(iROI);
Zlob_pre_aligned=RSFC_Lobe_Z_pre(iROI,:); Zlob_post_aligned=RSFC_Lobe_Z_post(iROI,:); Zlob_diff_aligned=RSFC_Lobe_Z_diff(iROI,:);
thi_pre_roi=DeltaTbl.THI_Pre(iTHI_roi); thi_post_roi=DeltaTbl.THI_Post(iTHI_roi); d_thi_roi=DeltaTbl.Delta_THI(iTHI_roi);
tfi_pre_roi=DeltaTbl.TFI_Pre(iTHI_roi); tfi_post_roi=DeltaTbl.TFI_Post(iTHI_roi); d_tfi_roi=DeltaTbl.Delta_TFI(iTHI_roi);
mean_thi_roi=(thi_pre_roi+thi_post_roi)/2; mean_tfi_roi=(tfi_pre_roi+tfi_post_roi)/2;

%% --- Align Whole-Brain RSFC with THI/TFI ---
ids_rsfc_wb = string(RSFC_WholeBrain.SubjectID(:));
[commonIDs_wb, iWB, iTHI_wb] = intersect(ids_rsfc_wb, ids_thi, 'stable');
nWB = numel(commonIDs_wb);
fprintf('\nWhole-brain: subjects with BOTH RSFC and THI/TFI = %d\n', nWB);
disp(commonIDs_wb);
rsfc_wb_pre=RSFC_WholeBrain.RSFC_pre(iWB); rsfc_wb_post=RSFC_WholeBrain.RSFC_post(iWB); rsfc_wb_delta=RSFC_WholeBrain.RSFC_delta(iWB);
thi_pre_wb=DeltaTbl.THI_Pre(iTHI_wb); thi_post_wb=DeltaTbl.THI_Post(iTHI_wb); d_thi_wb=DeltaTbl.Delta_THI(iTHI_wb);
tfi_pre_wb=DeltaTbl.TFI_Pre(iTHI_wb); tfi_post_wb=DeltaTbl.TFI_Post(iTHI_wb); d_tfi_wb=DeltaTbl.Delta_TFI(iTHI_wb);
mean_thi_wb=(thi_pre_wb+thi_post_wb)/2; mean_tfi_wb=(tfi_pre_wb+tfi_post_wb)/2;

%% --- ROI-seed RSFC x THI/TFI (THI->Spearman, TFI->Pearson) ---
fprintf('\n============================================================\n');
fprintf('ROI-SEED OVERALL RSFC x THI/TFI  (N=%d)\n', nROI);
fprintf('============================================================\n');
fprintf('%-22s  %-12s  %-9s  %6s  %8s\n','RSFC metric','Score','Test','stat','p');
fprintf('%s\n', repmat('-',1,70));

corrPairs_roi = {
    rsfc_roi_diff,  'RSFC dZ (ROI)',    d_thi_roi,   'dTHI';
    rsfc_roi_diff,  'RSFC dZ (ROI)',    d_tfi_roi,   'dTFI';
    rsfc_roi_diff,  'RSFC dZ (ROI)',    mean_thi_roi,'Mean THI';
    rsfc_roi_diff,  'RSFC dZ (ROI)',    mean_tfi_roi,'Mean TFI';
    rsfc_roi_pre,   'RSFC Z PRE (ROI)', thi_pre_roi, 'THI Pre';
    rsfc_roi_pre,   'RSFC Z PRE (ROI)', tfi_pre_roi, 'TFI Pre';
    rsfc_roi_post,  'RSFC Z POST (ROI)',thi_post_roi,'THI Post';
    rsfc_roi_post,  'RSFC Z POST (ROI)',tfi_post_roi,'TFI Post';
};
ROI_CorrResults = table();
for k = 1:size(corrPairs_roi,1)
    x=corrPairs_roi{k,1}; xLab=corrPairs_roi{k,2}; y=corrPairs_roi{k,3}; yLab=corrPairs_roi{k,4};
    tType = testForScore(yLab);
    ok=isfinite(x)&isfinite(y);
    if sum(ok)<5, fprintf('%-22s  %-12s  SKIPPED (n=%d)\n',xLab,yLab,sum(ok)); continue; end
    [rval,pval]=corr(x(ok),y(ok),'Type',tType);
    sig=pval<0.05;
    fprintf('%-22s  %-12s  %-9s  %6.3f  %8.4g  %s\n',xLab,yLab,tType,rval,pval,repmat('*',1,sig));
    ROI_CorrResults=[ROI_CorrResults; table({xLab},{yLab},{tType},rval,pval,sig,sum(ok), ...
        'VariableNames',{'RSFC_metric','Score','Test','r','p','sig','n'})];
end

%% --- Whole-brain RSFC x THI/TFI (THI->Spearman, TFI->Pearson) ---
fprintf('\n============================================================\n');
fprintf('WHOLE-BRAIN RSFC x THI/TFI  (N=%d)\n', nWB);
fprintf('============================================================\n');
fprintf('%-26s  %-12s  %-9s  %6s  %8s\n','RSFC metric','Score','Test','stat','p');
fprintf('%s\n', repmat('-',1,74));

corrPairs_wb = {
    rsfc_wb_delta, 'RSFC dr (WB)',   d_thi_wb,   'dTHI';
    rsfc_wb_delta, 'RSFC dr (WB)',   d_tfi_wb,   'dTFI';
    rsfc_wb_delta, 'RSFC dr (WB)',   mean_thi_wb,'Mean THI';
    rsfc_wb_delta, 'RSFC dr (WB)',   mean_tfi_wb,'Mean TFI';
    rsfc_wb_pre,   'RSFC r PRE (WB)',thi_pre_wb,'THI Pre';
    rsfc_wb_pre,   'RSFC r PRE (WB)',tfi_pre_wb,'TFI Pre';
    rsfc_wb_post,  'RSFC r POST (WB)',thi_post_wb,'THI Post';
    rsfc_wb_post,  'RSFC r POST (WB)',tfi_post_wb,'TFI Post';
};
WB_CorrResults = table();
for k = 1:size(corrPairs_wb,1)
    x=corrPairs_wb{k,1}; xLab=corrPairs_wb{k,2}; y=corrPairs_wb{k,3}; yLab=corrPairs_wb{k,4};
    tType = testForScore(yLab);
    ok=isfinite(x)&isfinite(y);
    if sum(ok)<5, fprintf('%-26s  %-12s  SKIPPED (n=%d)\n',xLab,yLab,sum(ok)); continue; end
    [rval,pval]=corr(x(ok),y(ok),'Type',tType);
    sig=pval<0.05;
    fprintf('%-26s  %-12s  %-9s  %6.3f  %8.4g  %s\n',xLab,yLab,tType,rval,pval,repmat('*',1,sig));
    WB_CorrResults=[WB_CorrResults; table({xLab},{yLab},{tType},rval,pval,sig,sum(ok), ...
        'VariableNames',{'RSFC_metric','Score','Test','r','p','sig','n'})];
end

%% --- Lobe-wise RSFC dZ x THI/TFI (THI->Spearman, TFI->Pearson) ---
fprintf('\n============================================================\n');
fprintf('LOBE-WISE RSFC dZ x THI/TFI  (N=%d)\n', nROI);
fprintf('============================================================\n');
scoreVecs  = {d_thi_roi, d_tfi_roi, mean_thi_roi, mean_tfi_roi};
scoreNames = {'dTHI','dTFI','Mean THI','Mean TFI'};
scoreTests = {testForScore('dTHI'), testForScore('dTFI'), testForScore('Mean THI'), testForScore('Mean TFI')};

lobeCorr_r = NaN(numel(Lobes), numel(scoreVecs));
lobeCorr_p = NaN(numel(Lobes), numel(scoreVecs));
for sc = 1:numel(scoreVecs)
    score=scoreVecs{sc}; sLab=scoreNames{sc}; tType=scoreTests{sc};
    fprintf('\n--- %s [%s] ---\n', sLab, tType);
    fprintf('%-6s  %6s  %8s\n','Lobe','stat','p');
    for L = 1:numel(Lobes)
        x=Zlob_diff_aligned(:,L); ok=isfinite(x)&isfinite(score);
        if sum(ok)<5, continue; end
        [rval,pval]=corr(x(ok),score(ok),'Type',tType);
        lobeCorr_r(L,sc)=rval; lobeCorr_p(L,sc)=pval;
        sig=pval<0.05; trend=pval<0.10 && ~sig; marker='';
        if sig, marker=' *'; end
        if trend, marker=' (trend)'; end
        fprintf('%-6s  %6.3f  %8.4g%s\n',Lobes{L},rval,pval,marker);
    end
end

%  SAVE QUESTIONNAIRE OUTPUTS

QUEST = struct();

% ---- merged participant-level table ----
QUEST.Merged   = MergedTbl;
QUEST.DeltaTbl = DeltaTbl;

% ---- item-level records, needed for subscale scoring ----
QUEST.Items.preTbl  = preTbl;
QUEST.Items.postTbl = postTbl;
QUEST.Items.ids     = ids;
QUEST.Items.keepIdx = keepIdx;

% ---- total scores ----
QUEST.Scores.THI_pre  = THI_pre;
QUEST.Scores.THI_post = THI_post;
QUEST.Scores.dTHI     = dTHI;
QUEST.Scores.TFI_pre  = TFI_pre;
QUEST.Scores.TFI_post = TFI_post;
QUEST.Scores.dTFI     = dTFI;

% ---- connectivity aligned to questionnaire row order ----
QUEST.Aligned.Lobes             = Lobes;
QUEST.Aligned.rsfc_roi_pre      = rsfc_roi_pre;
QUEST.Aligned.rsfc_roi_post     = rsfc_roi_post;
QUEST.Aligned.rsfc_roi_diff     = rsfc_roi_diff;
QUEST.Aligned.Zlob_pre_aligned  = Zlob_pre_aligned;
QUEST.Aligned.Zlob_post_aligned = Zlob_post_aligned;
QUEST.Aligned.Zlob_diff_aligned = Zlob_diff_aligned;
QUEST.Aligned.rsfc_wb_pre       = rsfc_wb_pre;
QUEST.Aligned.rsfc_wb_post      = rsfc_wb_post;
QUEST.Aligned.rsfc_wb_delta     = rsfc_wb_delta;

% ---- questionnaire scores matched to each connectivity set ----
QUEST.Aligned.d_thi_roi    = d_thi_roi;
QUEST.Aligned.d_tfi_roi    = d_tfi_roi;
QUEST.Aligned.mean_thi_roi = mean_thi_roi;
QUEST.Aligned.mean_tfi_roi = mean_tfi_roi;
QUEST.Aligned.thi_pre_roi  = thi_pre_roi;
QUEST.Aligned.tfi_pre_roi  = tfi_pre_roi;
QUEST.Aligned.thi_post_roi = thi_post_roi;
QUEST.Aligned.tfi_post_roi = tfi_post_roi;
QUEST.Aligned.d_thi_wb     = d_thi_wb;
QUEST.Aligned.d_tfi_wb     = d_tfi_wb;
QUEST.Aligned.mean_thi_wb  = mean_thi_wb;
QUEST.Aligned.mean_tfi_wb  = mean_tfi_wb;
QUEST.Aligned.thi_pre_wb   = thi_pre_wb;
QUEST.Aligned.tfi_pre_wb   = tfi_pre_wb;
QUEST.Aligned.thi_post_wb  = thi_post_wb;
QUEST.Aligned.tfi_post_wb  = tfi_post_wb;

QUEST.Info.Created = datetime('now');

save(fullfile(DATA_DIR, 'QUEST_outputs.mat'), 'QUEST');
fprintf('\nQuestionnaire outputs saved to %s\n', fullfile(DATA_DIR, 'QUEST_outputs.mat'));
%% BONFERRONI CORRECTION FOR ALL FAMILIES
% (each pair corrected with its own test type)

% --- 1) THI/TFI pre vs post (2 tests) ---
p_thi_bonf = min(p_thi*2,1); p_tfi_bonf = min(p_tfi*2,1);
fprintf('\n=== PRE vs POST (Bonferroni, 2 tests) ===\n');
fprintf('THI: uncorrected p=%.4g | Bonferroni p=%.4g\n', p_thi, p_thi_bonf);
fprintf('TFI: uncorrected p=%.4g | Bonferroni p=%.4g\n', p_tfi, p_tfi_bonf);

% --- 2) HbO BBN family (8 tests) ---
bbn_pairs = {
    MergedTbl.early_BBN, MergedTbl.Delta_THI, 'dTHI ~ BBN early';
    MergedTbl.late_BBN,  MergedTbl.Delta_THI, 'dTHI ~ BBN late';
    MergedTbl.early_BBN, thi_avg,             'mean THI ~ BBN early';
    MergedTbl.late_BBN,  thi_avg,             'mean THI ~ BBN late';
    MergedTbl.early_BBN, MergedTbl.Delta_TFI, 'dTFI ~ BBN early';
    MergedTbl.late_BBN,  MergedTbl.Delta_TFI, 'dTFI ~ BBN late';
    MergedTbl.early_BBN, tfi_avg,             'mean TFI ~ BBN early';
    MergedTbl.late_BBN,  tfi_avg,             'mean TFI ~ BBN late';
};
runFamily('HbO BBN', bbn_pairs, testForScore);

% --- 3) HbO ISR family (8 tests) ---
isr_pairs = {
    MergedTbl.early_ISR, MergedTbl.Delta_THI, 'dTHI ~ ISR early';
    MergedTbl.late_ISR,  MergedTbl.Delta_THI, 'dTHI ~ ISR late';
    MergedTbl.early_ISR, thi_avg,             'mean THI ~ ISR early';
    MergedTbl.late_ISR,  thi_avg,             'mean THI ~ ISR late';
    MergedTbl.early_ISR, MergedTbl.Delta_TFI, 'dTFI ~ ISR early';
    MergedTbl.late_ISR,  MergedTbl.Delta_TFI, 'dTFI ~ ISR late';
    MergedTbl.early_ISR, tfi_avg,             'mean TFI ~ ISR early';
    MergedTbl.late_ISR,  tfi_avg,             'mean TFI ~ ISR late';
};
runFamily('HbO ISR', isr_pairs, testForScore);

% --- 4) ROI-seed RSFC family (8 tests) ---
nTests_rsfc_roi = height(ROI_CorrResults);
rsfc_roi_p_bonf = min(ROI_CorrResults.p * nTests_rsfc_roi, 1);
fprintf('\n=== ROI-seed RSFC Correlations (Bonferroni, %d tests) ===\n', nTests_rsfc_roi);
fprintf('%-22s  %-12s  %-9s  %6s  %8s  %12s  %s\n','RSFC metric','Score','Test','stat','p_uncorr','p_bonf','sig');
for k = 1:height(ROI_CorrResults)
    sig = rsfc_roi_p_bonf(k) < 0.05;
    fprintf('%-22s  %-12s  %-9s  %6.3f  %8.4g  %12.4g  %s\n', ...
        ROI_CorrResults.RSFC_metric{k}, ROI_CorrResults.Score{k}, ROI_CorrResults.Test{k}, ...
        ROI_CorrResults.r(k), ROI_CorrResults.p(k), rsfc_roi_p_bonf(k), repmat('*',1,sig));
end

% --- 5) Whole-brain RSFC family (8 tests) ---
nTests_rsfc_wb = height(WB_CorrResults);
rsfc_wb_p_bonf = min(WB_CorrResults.p * nTests_rsfc_wb, 1);
fprintf('\n=== Whole-brain RSFC Correlations (Bonferroni, %d tests) ===\n', nTests_rsfc_wb);
fprintf('%-26s  %-12s  %-9s  %6s  %8s  %12s  %s\n','RSFC metric','Score','Test','stat','p_uncorr','p_bonf','sig');
for k = 1:height(WB_CorrResults)
    sig = rsfc_wb_p_bonf(k) < 0.05;
    fprintf('%-26s  %-12s  %-9s  %6.3f  %8.4g  %12.4g  %s\n', ...
        WB_CorrResults.RSFC_metric{k}, WB_CorrResults.Score{k}, WB_CorrResults.Test{k}, ...
        WB_CorrResults.r(k), WB_CorrResults.p(k), rsfc_wb_p_bonf(k), repmat('*',1,sig));
end

% --- 6) Lobe-wise RSFC family (9 lobes x 4 scores = 36 tests) ---
nTests_lobe = numel(Lobes) * numel(scoreVecs);
lobe_p_bonf = min(lobeCorr_p(:) * nTests_lobe, 1);
lobeCorr_p_bonf = reshape(lobe_p_bonf, size(lobeCorr_p));
fprintf('\n=== Lobe-wise RSFC Correlations (Bonferroni, %d tests) ===\n', nTests_lobe);
fprintf('%-6s  %-9s  %-9s  %6s  %8s  %12s  %s\n','Lobe','Score','Test','stat','p_uncorr','p_bonf','sig');
for sc = 1:numel(scoreVecs)
    for L = 1:numel(Lobes)
        if ~isfinite(lobeCorr_p(L,sc)), continue; end
        sig = lobeCorr_p_bonf(L,sc) < 0.05;
        fprintf('%-6s  %-9s  %-9s  %6.3f  %8.4g  %12.4g  %s\n', ...
            Lobes{L}, scoreNames{sc}, scoreTests{sc}, lobeCorr_r(L,sc), ...
            lobeCorr_p(L,sc), lobeCorr_p_bonf(L,sc), repmat('*',1,sig));
    end
end

fprintf('\nScript complete. THI correlations = Spearman, TFI correlations = Pearson.\n');
close all;

%% FIG 10: dTFI vs late ISR-evoked HbO (ROI)  -- surviving correlation
% Saved at 300 DPI, no title

saveDir = FIG_DIR;
if ~exist(saveDir, 'dir'), mkdir(saveDir); end

% --- data ---
x_all = MergedTbl.late_ISR;      % late ISR HbO (ROI)
y_all = MergedTbl.Delta_TFI;     % dTFI
ok = isfinite(x_all) & isfinite(y_all);
x = x_all(ok); y = y_all(ok); n = sum(ok);

% --- Pearson (TFI is normal) ---
[r, p] = corr(x, y, 'Type', 'Pearson');
fprintf('Fig10: dTFI vs late ISR  | r=%.4f, p=%.4g, n=%d\n', r, p, n);

% --- figure ---
fig = figure('Color','w','Units','inches','Position',[1 1 5 4]);
hold on;

scatter(x, y, 70, 'filled', ...
    'MarkerFaceColor', [0.2 0.6 0.2], 'MarkerEdgeColor', 'k');

% regression line
coeffs = polyfit(x, y, 1);
xFit = linspace(min(x), max(x), 100);
plot(xFit, polyval(coeffs, xFit), 'r-', 'LineWidth', 2);

% reference lines at 0
yline(0,'--','Color',[0.5 0.5 0.5]);
xline(0,'--','Color',[0.5 0.5 0.5]);

xlabel('Mean HbO ISR – ROI (late window, \muM)', 'FontSize', 13);
ylabel('\DeltaTFI (Post – Pre)', 'FontSize', 13);

% NO title (per request)
set(gca, 'FontSize', 12, 'LineWidth', 1.2);
box on; grid on;
hold off;

% --- save at 300 DPI ---
% outFile = fullfile(saveDir, 'Fig10_dTFI_vs_lateISR');
% print(fig, outFile, '-dtiff', '-r300');   % TIFF, 300 dpi
% % also save PNG version if wanted:
% print(fig, outFile, '-dpng', '-r300');     % PNG, 300 dpi
% 
% fprintf('Saved Fig 10 to:\n  %s.tif\n  %s.png\n', outFile, outFile);


%%  LOCAL FUNCTIONS
%  MATLAB requires all local functions to appear at the end of a script.


function t = ternaryTest(isTHI)
    if isTHI, t = 'Spearman'; else, t = 'Pearson'; end
end

function runFamily(famName, pairs, testForScore)
    nT = size(pairs,1);
    rv = zeros(nT,1); pv = zeros(nT,1); tt = cell(nT,1);
    for k = 1:nT
        lbl = pairs{k,3};
        tType = testForScore(lbl);
        tt{k} = tType;
        ok = isfinite(pairs{k,1}) & isfinite(pairs{k,2});
        [rv(k), pv(k)] = corr(pairs{k,1}(ok), pairs{k,2}(ok), 'Type', tType);
    end
    pb = min(pv * nT, 1);
    fprintf('\n=== %s Correlations (Bonferroni, %d tests) ===\n', famName, nT);
    fprintf('%-22s  %-9s  %6s  %8s  %12s  %s\n','Pair','Test','stat','p_uncorr','p_bonf','sig');
    for k = 1:nT
        sig = pb(k) < 0.05;
        fprintf('%-22s  %-9s  %6.3f  %8.4g  %12.4g  %s\n', pairs{k,3}, tt{k}, rv(k), pv(k), pb(k), repmat('*',1,sig));
    end
end