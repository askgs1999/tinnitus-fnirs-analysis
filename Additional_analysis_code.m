%%  Supplementary analyses performed  (v1.1.0, PTA-4 corrected)

%  Satish G, Taylor LM, Arnold MP, Gallagher-Shale J, Krishnamurthy K,
%  Basura GJ. Subjective tinnitus distress correlates with cortical
%  hemodynamic changes in human auditory cortex as measured by functional
%  near-infrared spectroscopy. PLOS ONE.
%
%  WHAT THIS SCRIPT CONTAINS
%
%  Analysis 1   Audiometric thresholds, and whether sensation level
%               predicts response amplitude. Also reports the between-group
%               HbO comparisons adjusted for PTA-4 and age, and a
%               sensitivity analysis restricted to PTA-4 at or below
%               25 dB HL.
%  Analysis 2   Baseline distress entered as a covariate on the
%               within-group BBN-versus-ISR contrast.
%  Analysis 3   Influence diagnostics for the change-in-TFI correlation:
%               Spearman and Kendall coefficients, Cook's distances,
%               bisquare robust regression, leave-one-out refitting, and a
%               bootstrap confidence interval.
%  Analysis 4   THI and TFI subscale scoring, and subscale associations
%               with HbO in both analysis windows.
%  Analysis 5   Between-group comparisons adjusted for age and PTA-4, for
%               connectivity and for the hemodynamic response, plus the
%               sex comparison and the covariate-adjusted change-in-TFI
%               correlation.
%  Analysis 6   Distributions of the pre-session THI and TFI scores,
%               including severity-band frequencies.
%  Analysis 7   Baseline distress against noise-evoked HbO, and whether
%               individual connectivity change tracks individual
%               questionnaire change.
%  Analysis 8   Sensitivity analysis giving the smallest effects
%               detectable at 80% power with this sample.
%
%  REQUIREMENTS
%    MATLAB R2024b
%    Statistics and Machine Learning Toolbox
%
%  INPUT   all three files are written by the main analysis scripts
%    RSFC_outputs.mat    written by rsfc_analysis.m
%    HR_outputs.mat      written by hemodynamic_analysis.m
%    QUEST_outputs.mat   written by questionnaire_analysis.m
%
%  HOW TO RUN
%    Run the three main scripts first, in any order except that
%    questionnaire_analysis.m must come last. Then set MATLAB's current
%    folder to the one holding the three .mat files and run this script.
%
%    Sections within this script must run in the order given, because later
%    sections use tables built by earlier ones. In particular Analysis 1
%    builds the audiometric table that Analysis 5 needs.
%
%  AUDIOMETRIC DATA
%    Four-frequency pure-tone averages (500, 1000, 2000, 4000 Hz) and ages
%    are entered directly in Analysis 1, because they are not held in the
%    fNIRS data files. PTA-4 values were recalculated from the original
%    audiograms and match S1_Table.xlsx. Tinnitus08 has no frequency-specific
%    audiogram, so its PTA-4 is NaN and it is left out of every PTA analysis.
%    Control ages appear a second time in the Point 13 section, ordered to
%    match the hemodynamic curve matrices.
%
%  OUTPUT
%    Result tables are written to RESULTS_DIR as CSV files and assigned to
%    the base workspace. Analysis 3 also draws one diagnostic figure.


clearvars; close all; clc;

%% ---- PATHS ----


INPUT_FILE = 'QUEST_outputs.mat';

found = which(INPUT_FILE);
if isempty(found) && isfile(fullfile(pwd, INPUT_FILE))
    found = fullfile(pwd, INPUT_FILE);
end
if isempty(found)
    error(['Cannot find %s.\n' ...
           'Current folder is %s.\n' ...
           'Run rsfc_analysis.m, hemodynamic_analysis.m and ' ...
           'questionnaire_analysis.m first, then cd to the folder holding ' ...
           'their output files.'], INPUT_FILE, pwd);
end

DATA_DIR    = fileparts(found);
RESULTS_DIR = fullfile(DATA_DIR, 'results');
if ~exist(RESULTS_DIR,'dir'), mkdir(RESULTS_DIR); end

for f = {'RSFC_outputs.mat','HR_outputs.mat','QUEST_outputs.mat'}
    if ~isfile(fullfile(DATA_DIR, f{1}))
        error('%s not found in %s.', f{1}, DATA_DIR);
    end
end

load(fullfile(DATA_DIR, 'RSFC_outputs.mat'));    % struct RSFC
load(fullfile(DATA_DIR, 'HR_outputs.mat'));      % struct HR
load(fullfile(DATA_DIR, 'QUEST_outputs.mat'));   % struct QUEST

fprintf('Reading data from : %s\n', DATA_DIR);
fprintf('Result tables to  : %s\n', RESULTS_DIR);

% ---- unpack into the variable names the sections below use ----

% questionnaire
MergedTbl = QUEST.Merged;
DeltaTbl  = QUEST.DeltaTbl;
preTbl    = QUEST.Items.preTbl;
postTbl   = QUEST.Items.postTbl;
ids       = QUEST.Items.ids;
keepIdx   = QUEST.Items.keepIdx;
THI_pre   = QUEST.Scores.THI_pre;
THI_post  = QUEST.Scores.THI_post;
dTHI      = QUEST.Scores.dTHI;
TFI_pre   = QUEST.Scores.TFI_pre;
TFI_post  = QUEST.Scores.TFI_post;
dTFI      = QUEST.Scores.dTFI;

% hemodynamic window means and curves
T_CorrData_TINN = HR.Tinnitus.ROI;
tEpoch          = HR.Curves.tEpoch;
CTRL_ROI_BBN    = HR.Curves.CTRL_ROI_BBN;
CTRL_ROI_ISR    = HR.Curves.CTRL_ROI_ISR;
CTRL_NROI_BBN   = HR.Curves.CTRL_NROI_BBN;
CTRL_NROI_ISR   = HR.Curves.CTRL_NROI_ISR;
TINN_ROI_BBN    = HR.Curves.TINN_ROI_BBN;
TINN_ROI_ISR    = HR.Curves.TINN_ROI_ISR;
TINN_NROI_BBN   = HR.Curves.TINN_NROI_BBN;
TINN_NROI_ISR   = HR.Curves.TINN_NROI_ISR;

% connectivity, pooled ROI seed
zPre_exp   = RSFC.ROI.Z_pre;
zPost_exp  = RSFC.ROI.Z_post;
zDiff_exp  = RSFC.ROI.Z_diff;
zPre_ctrl  = RSFC.Control.ROI_Z_pre;
zPost_ctrl = RSFC.Control.ROI_Z_post;
zDiff_ctrl = RSFC.Control.ROI_Z_diff;

% connectivity, lobe-wise
Zlob_pre_exp   = RSFC.Lobe.Z_pre;
Zlob_post_exp  = RSFC.Lobe.Z_post;
Zlob_diff_exp  = RSFC.Lobe.Z_diff;
Zlob_pre_ctrl  = RSFC.Control.Lobe_Z_pre;
Zlob_post_ctrl = RSFC.Control.Lobe_Z_post;
Zlob_diff_ctrl = RSFC.Control.Lobe_Z_diff;
lobes          = RSFC.Lobes;

% connectivity aligned to questionnaire row order
Lobes             = QUEST.Aligned.Lobes;
rsfc_roi_diff     = QUEST.Aligned.rsfc_roi_diff;
Zlob_diff_aligned = QUEST.Aligned.Zlob_diff_aligned;
rsfc_wb_delta     = QUEST.Aligned.rsfc_wb_delta;
d_thi_roi    = QUEST.Aligned.d_thi_roi;
d_tfi_roi    = QUEST.Aligned.d_tfi_roi;
mean_thi_roi = QUEST.Aligned.mean_thi_roi;
mean_tfi_roi = QUEST.Aligned.mean_tfi_roi;
d_thi_wb     = QUEST.Aligned.d_thi_wb;
d_tfi_wb     = QUEST.Aligned.d_tfi_wb;
mean_thi_wb  = QUEST.Aligned.mean_thi_wb;
mean_tfi_wb  = QUEST.Aligned.mean_tfi_wb;

% test-selection rule: THI is analysed with Spearman, TFI with Pearson
testForScore = @(scoreName) ternaryTest(contains(lower(scoreName),'thi'));



%% Analysis 1
%  POINT 4 — COMPLETE
%  (1) report thresholds
%  (2) does sensation level predict response amplitude?
%      - within tinnitus group
%      - within control group
%      - between-group HbO adjusted for PTA + age
%      - sensitivity: exclude PTA > 25 dB HL
%
%  SL = 65 dB SPL - PTA, so SL is a constant minus PTA.
%  A NEGATIVE PTA slope = lower SL -> smaller response, which
%  is the peripheral account the reviewer raises.
% =========================================================

STOP_ON_ID_MISMATCH = true;

% PTA AND AGE LOOKUPS
% PTA-4 recalculated from the original audiograms (matches S1_Table.xlsx).
% Tinnitus08 has no frequency-specific audiogram -> NaN.

tinnIDs = {'Tinnitus01','Tinnitus12','Tinnitus13','Tinnitus18','Tinnitus19', ...
           'Tinnitus21','Tinnitus23','Tinnitus28','Tinnitus32','Tinnitus34', ...
           'Tinnitus02','Tinnitus38','Tinnitus39','Tinnitus42','Tinnitus46', ...
           'Tinnitus47','Tinnitus59','Tinnitus61','Tinnitus62','Tinnitus69', ...
           'Tinnitus03','Tinnitus70','Tinnitus04','Tinnitus05','Tinnitus06', ...
           'Tinnitus07','Tinnitus08'};
tinnR   = [15 8.75 17.5 25 26.25 15 23.75 21.25 21.25 30 6.25 13.75 11.25 27.5 ...
           2.5 20 10 12.5 16.25 18.75 16.25 17.5 12.5 13.75 30 30 NaN];
tinnL   = [20 7.5 15 28.75 31.25 18.75 17.5 22.5 20 30 5 12.5 20 27.5 ...
           5 16.25 15 15 17.5 21.25 18.75 15 13.75 18.75 27.5 28.75 NaN];
tinnAge = [45 63 63 70 52 63 61 65 69 72 49 59 70 68 65 68 66 63 65 60 62 64 46 64 54 65 53];

ctrlIDs = {'Tinnitus09','Tinnitus26','Tinnitus25','Tinnitus27','Tinnitus31', ...
           'Tinnitus40','Tinnitus41','Tinnitus43','Tinnitus44','Tinnitus45', ...
           'Tinnitus48','Tinnitus50','Tinnitus51','Tinnitus52','Tinnitus54', ...
           'Tinnitus60','Tinnitus63','Tinnitus64','Tinnitus65','Tinnitus66', ...
           'Tinnitus67','Tinnitus68','Tinnitus73','Tinnitus16','Tinnitus17', ...
           'Tinnitus20','Tinnitus24'};
ctrlR   = [8.75 10 8.75 6.25 12.5 5 18.75 8.75 20 6.25 17.5 8.75 2.5 11.25 ...
           11.25 21.25 16.25 18.75 11.25 7.5 10 8.75 22.5 8.75 7.5 25 7.5];
ctrlL   = [10 6.25 10 11.25 15 11.25 22.5 10 25 7.5 17.5 8.75 2.5 16.25 ...
           10 18.75 13.75 21.25 12.5 16.25 6.25 13.75 13.75 12.5 6.25 26.25 2.5];
ctrlAge = [45 41 33 41 39 35 30 69 62 62 60 62 69 53 69 62 52 69 62 60 65 51 59 59 62 50 59];

assert(isequal(numel(tinnIDs),numel(tinnR),numel(tinnL),numel(tinnAge)), 'tinnitus vector length mismatch');
assert(isequal(numel(ctrlIDs),numel(ctrlR),numel(ctrlL),numel(ctrlAge)), 'control vector length mismatch');

pT = (tinnR + tinnL)/2;
pC = (ctrlR + ctrlL)/2;

% ---- ID normalizer: strip spaces, case-fold ----
nrm = @(s) upper(string(regexprep(cellstr(string(s)), '\s+', '')));

AudTbl = table( nrm([tinnIDs(:); ctrlIDs(:)]), ...
                [pT(:); pC(:)], ...
                [tinnR(:); ctrlR(:)], [tinnL(:); ctrlL(:)], ...
                [tinnAge(:); ctrlAge(:)], ...
                [ones(numel(pT),1); zeros(numel(pC),1)], ...
    'VariableNames', {'SubjectID','PTA','PTA_R','PTA_L','Age','Group'});

fprintf('\nAge check (should match Table 1): tinn %.2f (%.2f) | ctrl %.2f (%.2f)\n', ...
    mean(tinnAge), std(tinnAge), mean(ctrlAge), std(ctrlAge));


% 1) THRESHOLDS BETWEEN GROUPS

fprintf('\n=================================================================\n');
fprintf(' 1) AUDIOMETRIC THRESHOLDS (PTA-4: 0.5, 1, 2, 4 kHz)\n');
fprintf('=================================================================\n');
fprintf('%-14s %-28s %-28s %s\n','Measure','Tinnitus','Control','Welch test');
fprintf('%s\n', repmat('-',1,96));
for v = {'PTA','PTA_R','PTA_L','Age'}
    a = AudTbl.(v{1})(AudTbl.Group==1);
    b = AudTbl.(v{1})(AudTbl.Group==0);
    a = a(isfinite(a)); b = b(isfinite(b));      % drop participants without PTA-4
    [~,pv,~,sv] = ttest2(a, b, 'Vartype','unequal');
    spool = sqrt(((numel(a)-1)*var(a) + (numel(b)-1)*var(b)) / (numel(a)+numel(b)-2));
    fprintf('%-14s %6.2f (%5.2f) %5.2f-%5.2f  %6.2f (%5.2f) %5.2f-%5.2f  t(%.1f)=%+.3f p=%.4f d=%.3f (n=%d vs %d)\n', ...
        v{1}, mean(a), std(a), min(a), max(a), mean(b), std(b), min(b), max(b), ...
        sv.df, sv.tstat, pv, (mean(a)-mean(b))/spool, numel(a), numel(b));
end
fprintf('\nPTA > 25 dB HL : tinnitus %d/%d, control %d/%d\n', ...
    sum(pT>25), sum(isfinite(pT)), sum(pC>25), sum(isfinite(pC)));
fprintf('All participants <= 30 dB HL : tinnitus %d/%d, control %d/%d\n', ...
    sum(pT<=30), sum(isfinite(pT)), sum(pC<=30), sum(isfinite(pC)));

% 2) LOAD CONTROL HbO

hboCols = {'early_BBN','early_ISR','late_BBN','late_ISR'};
assert(all(ismember(hboCols, MergedTbl.Properties.VariableNames)), ...
    'MergedTbl is missing HbO columns: %s', ...
    strjoin(setdiff(hboCols, MergedTbl.Properties.VariableNames), ', '));

% Control HbO window means come from HR_outputs.mat (hemodynamic_analysis.m)
Tc = HR.Control.ROI;
assert(all(ismember(hboCols, Tc.Properties.VariableNames)), ...
    'HR.Control.ROI is missing HbO columns: %s', ...
    strjoin(setdiff(hboCols, Tc.Properties.VariableNames), ', '));



% 2b) ID CHECK

Dt = MergedTbl;  Dt.SubjectID = nrm(Dt.SubjectID);
Dc = Tc;         Dc.SubjectID = nrm(Dc.SubjectID);

allHbO   = [Dt.SubjectID; Dc.SubjectID];
onlyHbO  = setdiff(allHbO, AudTbl.SubjectID);
onlyAud  = setdiff(AudTbl.SubjectID, allHbO);

fprintf('\n=================================================================\n');
fprintf(' 2b) ID CHECK (after stripping spaces / case)\n');
fprintf('=================================================================\n');
fprintf(' HbO rows: %d tinnitus + %d control = %d\n', height(Dt), height(Dc), numel(allHbO));
fprintf(' AudTbl rows: %d\n', height(AudTbl));
if isempty(onlyHbO) && isempty(onlyAud)
    fprintf(' All IDs matched.\n');
else
    fprintf(' In HbO but not AudTbl: %s\n', strjoin(cellstr(onlyHbO), ', '));
    fprintf(' In AudTbl but not HbO: %s\n', strjoin(cellstr(onlyAud), ', '));
    if STOP_ON_ID_MISMATCH
        error(['ID mismatch. These are different people, not a formatting issue. ' ...
               'Resolve which participants are actually in the analysis before ' ...
               'trusting Table 1 or the covariate models.']);
    else
        warning('Proceeding with unmatched IDs — those rows will be dropped.');
    end
end


% 3) PTA vs HbO, TINNITUS GROUP

[tf, loc] = ismember(Dt.SubjectID, AudTbl.SubjectID);
Dt.PTA = nan(height(Dt),1);  Dt.Age = nan(height(Dt),1);
Dt.PTA(tf) = AudTbl.PTA(loc(tf));
Dt.Age(tf) = AudTbl.Age(loc(tf));
Dt.delta_early = Dt.early_BBN - Dt.early_ISR;
Dt.delta_late  = Dt.late_BBN  - Dt.late_ISR;

quant = [hboCols, {'delta_early','delta_late'}];

fprintf('\n=================================================================\n');
fprintf(' 3) PTA vs HbO AMPLITUDE — TINNITUS GROUP\n');
fprintf('=================================================================\n');
fprintf(' negative slope = lower sensation level -> smaller response\n\n');
fprintf('%-12s %4s %9s %9s %9s %9s\n','Quantity','n','r','p','rho','p_rho');
fprintf('%s\n', repmat('-',1,56));
P4t = table();
for q = 1:numel(quant)
    y = Dt.(quant{q}); x = Dt.PTA;
    ok = isfinite(x) & isfinite(y);
    if sum(ok) < 5, continue; end
    [r,p]      = corr(x(ok), y(ok), 'Type','Pearson');
    [rho,prho] = corr(x(ok), y(ok), 'Type','Spearman');
    fprintf('%-12s %4d %+9.3f %9.4f %+9.3f %9.4f%s\n', quant{q}, sum(ok), r, p, rho, prho, p4bstar(p));
    P4t = [P4t; table({'tinnitus'}, quant(q), sum(ok), r, p, rho, prho, ...
        'VariableNames',{'Group','Quantity','n','r','p','rho','p_rho'})]; %#ok<AGROW>
end


% 4) PTA vs HbO, CONTROL GROUP

[tf2, loc2] = ismember(Dc.SubjectID, AudTbl.SubjectID);
Dc.PTA = nan(height(Dc),1);  Dc.Age = nan(height(Dc),1);
Dc.PTA(tf2) = AudTbl.PTA(loc2(tf2));
Dc.Age(tf2) = AudTbl.Age(loc2(tf2));
Dc.delta_early = Dc.early_BBN - Dc.early_ISR;
Dc.delta_late  = Dc.late_BBN  - Dc.late_ISR;

fprintf('\n=================================================================\n');
fprintf(' 4) PTA vs HbO AMPLITUDE — CONTROL GROUP\n');
fprintf('=================================================================\n');
fprintf('%-12s %4s %9s %9s %9s %9s\n','Quantity','n','r','p','rho','p_rho');
fprintf('%s\n', repmat('-',1,56));
P4c = table();
for q = 1:numel(quant)
    y = Dc.(quant{q}); x = Dc.PTA;
    ok = isfinite(x) & isfinite(y);
    if sum(ok) < 5, continue; end
    [r,p]      = corr(x(ok), y(ok), 'Type','Pearson');
    [rho,prho] = corr(x(ok), y(ok), 'Type','Spearman');
    fprintf('%-12s %4d %+9.3f %9.4f %+9.3f %9.4f%s\n', quant{q}, sum(ok), r, p, rho, prho, p4bstar(p));
    P4c = [P4c; table({'control'}, quant(q), sum(ok), r, p, rho, prho, ...
        'VariableNames',{'Group','Quantity','n','r','p','rho','p_rho'})]; %#ok<AGROW>
end


% 5) BETWEEN-GROUP HbO, UNADJUSTED vs ADJUSTED FOR PTA + AGE

keep = ['SubjectID', hboCols];
A = [Dt(:,keep); Dc(:,keep)];
A.Group = [ones(height(Dt),1); zeros(height(Dc),1)];
A.delta_early = A.early_BBN - A.early_ISR;
A.delta_late  = A.late_BBN  - A.late_ISR;
[tf3, loc3] = ismember(A.SubjectID, AudTbl.SubjectID);
A.PTA = nan(height(A),1);  A.Age = nan(height(A),1);
A.PTA(tf3) = AudTbl.PTA(loc3(tf3));
A.Age(tf3) = AudTbl.Age(loc3(tf3));
cc = isfinite(A.PTA) & isfinite(A.Age);

fprintf('\n=================================================================\n');
fprintf(' 5) BETWEEN-GROUP HbO — unadjusted vs adjusted (group|PTA|age)\n');
fprintf('=================================================================\n');
fprintf(' n = %d tinnitus + %d control; complete covariates for %d\n\n', ...
    sum(A.Group==1), sum(A.Group==0), sum(cc));
fprintf('%-12s | %-24s | %s\n','Quantity','Unadjusted (Welch)','Adjusted');
fprintf('%s\n', repmat('-',1,94));
P4b = table();
for q = 1:numel(quant)
    nm = quant{q};  y = A.(nm);
    ok0 = isfinite(y);
    [~,pu,~,su] = ttest2(y(ok0 & A.Group==1), y(ok0 & A.Group==0), 'Vartype','unequal');
    ok = isfinite(y) & cc;
    mdl = fitlm([A.Group(ok), A.PTA(ok)-mean(A.PTA(ok)), A.Age(ok)-mean(A.Age(ok))], y(ok), ...
                'VarNames',{'Group','PTA','Age',nm});
    bg = mdl.Coefficients.Estimate(2);
    pg = mdl.Coefficients.pValue(2);
    pp = mdl.Coefficients.pValue(3);
    pa = mdl.Coefficients.pValue(4);
    fprintf('%-12s | t=%+.3f p=%.4f%-4s | grp b=%+.4f p=%.4f%-4s PTA p=%.3f Age p=%.3f (n=%d)\n', ...
        nm, su.tstat, pu, p4bstar(pu), bg, pg, p4bstar(pg), pp, pa, sum(ok));
    P4b = [P4b; table(quant(q), sum(ok0), pu, su.tstat, sum(ok), bg, pg, pp, pa, ...
        'VariableNames',{'Quantity','n_unadj','p_unadj','t_unadj','n_adj', ...
                         'b_group','p_group_adj','p_PTA','p_Age'})]; %#ok<AGROW>
end


% 6) SENSITIVITY: drop everyone above 25 dB HL

sub = A(A.PTA <= 25, :);
[~,pS] = ttest2(sub.PTA(sub.Group==1), sub.PTA(sub.Group==0), 'Vartype','unequal');
fprintf('\n=================================================================\n');
fprintf(' 6) SENSITIVITY — participants with PTA <= 25 dB HL only\n');
fprintf('=================================================================\n');
fprintf(' retained: %d tinnitus, %d control | PTA now %.2f vs %.2f, p = %.4f\n\n', ...
    sum(sub.Group==1), sum(sub.Group==0), ...
    mean(sub.PTA(sub.Group==1)), mean(sub.PTA(sub.Group==0)), pS);
fprintf('%-12s %4s %9s %9s %9s\n','Quantity','n','t','p','Cohen d');
fprintf('%s\n', repmat('-',1,46));
for q = 1:numel(quant)
    y = sub.(quant{q});
    a = y(sub.Group==1 & isfinite(y));  b = y(sub.Group==0 & isfinite(y));
    if numel(a)<5 || numel(b)<5, continue; end
    [~,pv,~,sv] = ttest2(a, b, 'Vartype','unequal');
    spool = sqrt(((numel(a)-1)*var(a) + (numel(b)-1)*var(b)) / (numel(a)+numel(b)-2));
    fprintf('%-12s %4d %+9.3f %9.4f %+9.3f%s\n', quant{q}, numel(a)+numel(b), ...
        sv.tstat, pv, (mean(a)-mean(b))/spool, p4bstar(pv));
end


% 7) VERDICT

fprintf('\n=================================================================\n');
fprintf(' 7) VERDICT FOR POINT 4\n');
fprintf('=================================================================\n');
neg = P4t(P4t.p < 0.05 & P4t.r < 0, :);
if isempty(neg)
    fprintf('\nTinnitus group: no negative PTA-amplitude association in any\n');
    fprintf('window or condition. The peripheral account predicts one and\n');
    fprintf('it is absent.\n');
else
    fprintf('\nTinnitus group: NEGATIVE PTA association present — report it.\n');
    disp(neg);
end
flip = P4b(xor(P4b.p_unadj<0.05, P4b.p_group_adj<0.05), ...
           {'Quantity','p_unadj','p_group_adj'});
if isempty(flip)
    fprintf('\nBetween-group: no effect changed significance after adjustment.\n');
else
    fprintf('\nBetween-group effects that changed significance:\n');
    disp(flip);
end
fprintf('Was PTA ever a significant predictor? %s\n', p4bYN(any(P4b.p_PTA < 0.05)));
fprintf('Was age ever a significant predictor? %s\n', p4bYN(any(P4b.p_Age < 0.05)));

P4all = [P4t; P4c];
writetable(AudTbl, fullfile(RESULTS_DIR, 'Point4_audiometry.csv'));
writetable(P4all,  fullfile(RESULTS_DIR, 'Point4_PTA_vs_HbO_bothgroups.csv'));
writetable(P4b,    fullfile(RESULTS_DIR, 'Point4_between_group_adjusted.csv'));
assignin('base','AudTbl',AudTbl);
assignin('base','P4all',P4all);
assignin('base','P4b',P4b);
assignin('base','A',A);
fprintf('\nSaved: Point4_audiometry.csv, Point4_PTA_vs_HbO_bothgroups.csv, Point4_between_group_adjusted.csv\n');
fprintf('#################################################################\n');

%%  ANALYSIS 2
%  Baseline distress as a covariate

%  POINT 6 (CORRECTED): BASELINE DISTRESS AS A COVARIATE
%  ON THE CONTRAST THE PAPER ACTUALLY TESTS
%


%  Model:  y = b0 + b1 * (distress - mean(distress))
%    b0 = adjusted mean of y at average distress
%         -> does the effect REMAIN after adjustment?
%    b1 = does distress predict y?
%  The covariate is centred so b0 is the adjusted mean rather
%  than the value at distress = 0.


fprintf('\n\n#################################################################\n');
fprintf('# POINT 6 (CORRECTED) — DISTRESS COVARIATE ON BBN vs ISR\n');
fprintf('#################################################################\n');

D = MergedTbl;
req = {'late_ISR','late_BBN','early_ISR','early_BBN','THI_Pre','TFI_Pre'};
miss = req(~ismember(req, D.Properties.VariableNames));
if ~isempty(miss)
    error('MergedTbl is missing: %s', strjoin(miss, ', '));
end


% Build the contrasts the paper reports
%   delta = BBN - ISR  (the quantity in the between-group test)

D.delta_late  = D.late_BBN  - D.late_ISR;
D.delta_early = D.early_BBN - D.early_ISR;

contrasts = {
    'delta_late',  'BBN - ISR, late window  (paper''s contrast)';
    'delta_early', 'BBN - ISR, early window (paper''s contrast)';
    'late_ISR',    'ISR alone, late window  (the elevation)';
    'late_BBN',    'BBN alone, late window  (the suppression)';
};

% 1) UNADJUSTED — what the paper reports

fprintf('\n=================================================================\n');
fprintf(' 1) UNADJUSTED (paired t-test vs zero, tinnitus group)\n');
fprintf('=================================================================\n');
fprintf('%-13s %4s %10s %10s %9s %9s   %s\n', ...
    'Quantity','n','mean','SD','t','p','description');
fprintf('%s\n', repmat('-',1,96));
for c = 1:size(contrasts,1)
    y = D.(contrasts{c,1});  y = y(isfinite(y));
    [~,p,~,st] = ttest(y);
    fprintf('%-13s %4d %+10.4f %10.4f %+9.3f %9.4f%-4s %s\n', ...
        contrasts{c,1}, numel(y), mean(y), std(y), st.tstat, p, ...
        p6bstar(p), contrasts{c,2});
end
fprintf('\nNote: delta_late / delta_early vs zero is equivalent to the\n');
fprintf('paired within-group BBN vs ISR comparison reported in the paper.\n');


% 2) ADJUSTED for baseline distress

distVars = {'THI_Pre','TFI_Pre'};
R = table();

for d = 1:numel(distVars)
    dv = distVars{d};
    fprintf('\n=================================================================\n');
    fprintf(' 2.%d) ADJUSTED for %s\n', d, dv);
    fprintf('=================================================================\n');
    fprintf('%-13s %4s | %-26s | %-26s\n', ...
        'Quantity','n','ADJUSTED MEAN (b0)','DISTRESS SLOPE (b1)');
    fprintf('%s\n', repmat('-',1,76));

    for c = 1:size(contrasts,1)
        y = D.(contrasts{c,1});
        x = D.(dv);
        ok = isfinite(y) & isfinite(x);
        if sum(ok) < 5, continue; end
        y = y(ok);  x = x(ok);

        mdl = fitlm(x - mean(x), y);
        b0 = mdl.Coefficients.Estimate(1);  p0 = mdl.Coefficients.pValue(1);
        b1 = mdl.Coefficients.Estimate(2);  p1 = mdl.Coefficients.pValue(2);
        r2 = mdl.Rsquared.Ordinary;

        % Spearman as a robustness check on the slope
        [rho, pRho] = corr(x, y, 'Type','Spearman');

        fprintf('%-13s %4d | %+.4f  p=%7.4f%-4s | %+.5f p=%7.4f%-4s R2=%.3f (rho=%+.3f p=%.3f)\n', ...
            contrasts{c,1}, sum(ok), b0, p0, p6bstar(p0), ...
            b1, p1, p6bstar(p1), r2, rho, pRho);

        R = [R; table({contrasts{c,1}}, {dv}, sum(ok), b0, p0, b1, p1, r2, rho, pRho, ...
            'VariableNames', {'Quantity','Covariate','n','adjMean','p_adjMean', ...
                              'slope','p_slope','R2','rho','p_rho'})]; %#ok<AGROW>
    end
end


% 3) VERDICT

fprintf('\n=================================================================\n');
fprintf(' 3) VERDICT\n');
fprintf('=================================================================\n');

for c = 1:2   % the two contrasts the paper reports
    q = contrasts{c,1};
    y = D.(q); y = y(isfinite(y));
    [~,pu] = ttest(y);
    fprintf('\n%s (unadjusted p = %.4f)\n', q, pu);
    rows = R(strcmp(R.Quantity,q), :);
    for k = 1:height(rows)
        if pu < 0.05 && rows.p_adjMean(k) < 0.05
            v = 'significant before AND after adjustment';
        elseif pu < 0.05 && rows.p_adjMean(k) >= 0.05
            v = 'significant before, NOT after adjustment';
        elseif pu >= 0.05
            v = 'was not significant unadjusted either';
        else
            v = '';
        end
        fprintf('   adj. for %-8s : p = %.4f  -> %s\n', ...
            rows.Covariate{k}, rows.p_adjMean(k), v);
    end
end

fprintf('\nDistress measures that significantly predicted a quantity:\n');
sig = R(R.p_slope < 0.05, {'Quantity','Covariate','slope','p_slope','R2','p_rho'});
if isempty(sig)
    fprintf('  none\n');
else
    disp(sig);
    fprintf('(check p_rho: if Spearman disagrees, the slope is outlier-sensitive)\n');
end

fprintf('\n-----------------------------------------------------------------\n');
fprintf('DESIGN LIMITATION, to state in the response:\n');
fprintf('Controls were not administered the THI or TFI, so baseline\n');
fprintf('distress cannot be entered as a covariate in the BETWEEN-group\n');
fprintf('model. These analyses are within the tinnitus group only.\n');
fprintf('-----------------------------------------------------------------\n');

writetable(R,fullfile(RESULTS_DIR, 'Point6_corrected.csv'));
assignin('base','P6c',R);
fprintf('\nSaved: Point6_corrected.csv  (workspace: P6c)\n');
fprintf('#################################################################\n');




%%  ANALYSIS 3
%  Influence diagnostics for the change-in-TFI correlation

%  POINT 10 ROBUSTNESS CHECKS: dTFI vs late ISR HbO (ROI)

% ---- SET THE VARIABLES
x = MergedTbl.late_ISR;      % late-window ISR HbO, auditory ROI
y = MergedTbl.Delta_TFI;     % change in TFI (Post - Pre)
ok = isfinite(x) & isfinite(y);
x  = x(ok);
y  = y(ok);
n  = numel(x);

% ---- SANITY CHECK:
fprintf('\n--- SANITY CHECK ---\n');
fprintf('n            = %d          (expected 27)\n', n);
fprintf('x range      = %+.3f to %+.3f   (expected about -0.32 to 0.30)\n', min(x), max(x));
fprintf('y range      = %+.2f to %+.2f   (expected about -12 to +8)\n', min(y), max(y));
fprintf('y mean       = %+.2f            (expected about -3.6)\n', mean(y));
if mean(y) > 5
    warning('y looks like ABSOLUTE TFI, not change scores. Check Delta_TFI.');
end

fprintf('\n=========== FIG 10 ROBUSTNESS CHECKS (n=%d) ===========\n', n);

% ---- 1) Pearson (what the paper reports) ----
[r_p, p_p] = corr(x, y, 'Type', 'Pearson');
fprintf('\n1) Pearson       : r    = %.3f, p = %.4g\n', r_p, p_p);

% ---- 2) Spearman (rank-based; robust to outliers) ----
[r_s, p_s] = corr(x, y, 'Type', 'Spearman');
fprintf('2) Spearman      : rho  = %.3f, p = %.4g\n', r_s, p_s);

% ---- 3) Kendall (second rank-based check) ----
[r_k, p_k] = corr(x, y, 'Type', 'Kendall');
fprintf('3) Kendall       : tau  = %.3f, p = %.4g\n', r_k, p_k);

% ---- 4) Cook's distance ----
mdl  = fitlm(x, y);
cd   = mdl.Diagnostics.CooksDistance;
thr  = 4/n;
infl = find(cd > thr);
fprintf('\n4) Cook''s distance (threshold 4/n = %.3f)\n', thr);
fprintf('   max Cook''s D = %.3f (observation %d)\n', max(cd), find(cd==max(cd),1));
if isempty(infl)
    fprintf('   No observations exceed threshold.\n');
else
    fprintf('   %d observation(s) exceed threshold: %s\n', numel(infl), mat2str(infl'));
    for k = 1:numel(infl)
        i = infl(k);
        fprintf('     obs %2d: x = %+.4f, y = %+.2f, Cook''s D = %.3f\n', ...
                i, x(i), y(i), cd(i));
    end
end

% ---- 5) Robust (bisquare) regression ----
[b_rob, st_rob] = robustfit(x, y);
b_ols = polyfit(x, y, 1);
fprintf('\n5) Robust fit (bisquare)\n');
fprintf('   slope = %.2f, SE = %.2f, t = %.3f, p = %.4g\n', ...
        b_rob(2), st_rob.se(2), st_rob.t(2), st_rob.p(2));
fprintf('   (OLS slope for comparison = %.2f)\n', b_ols(1));

% ---- 6) Leave-one-out range ----
r_loo = nan(n,1); p_loo = nan(n,1);
for i = 1:n
    keep = true(n,1); keep(i) = false;
    [r_loo(i), p_loo(i)] = corr(x(keep), y(keep), 'Type', 'Pearson');
end
fprintf('\n6) Leave-one-out (Pearson)\n');
fprintf('   r ranges %.3f to %.3f\n', min(r_loo), max(r_loo));
fprintf('   p ranges %.4g to %.4g\n', min(p_loo), max(p_loo));
fprintf('   p < 0.05 in %d of %d leave-one-out fits\n', sum(p_loo<0.05), n);
[~, iWorst] = min(abs(r_loo));
fprintf('   Weakest fit obtained by dropping obs %d (r = %.3f, p = %.4g)\n', ...
        iWorst, r_loo(iWorst), p_loo(iWorst));

% ---- 7) Refit with influential observations removed ----
if ~isempty(infl)
    keep = true(n,1); keep(infl) = false;
    [r_x,  p_x ] = corr(x(keep), y(keep), 'Type', 'Pearson');
    [rs_x, ps_x] = corr(x(keep), y(keep), 'Type', 'Spearman');
    fprintf('\n7) With %d influential obs removed (n=%d)\n', numel(infl), sum(keep));
    fprintf('   Pearson  r   = %.3f, p = %.4g\n', r_x,  p_x);
    fprintf('   Spearman rho = %.3f, p = %.4g\n', rs_x, ps_x);
end

% ---- 8) Bootstrap CI on Pearson r ----
rng(1);
nboot = 10000;
bootr = bootstrp(nboot, @(idx) corr(x(idx), y(idx)), (1:n)');
ci    = prctile(bootr, [2.5 97.5]);
fprintf('\n8) Bootstrap (%d resamples)\n', nboot);
fprintf('   Pearson r = %.3f, 95%% CI [%.3f, %.3f]\n', r_p, ci(1), ci(2));
fprintf('   Proportion of resamples with r <= 0: %.3f\n', mean(bootr<=0));

% ---- 9) Multiple-comparison context (Reviewer point 9) ----
fprintf('\n9) Bonferroni context for this single correlation\n');
fprintf('   uncorrected p        = %.4g\n', p_p);
fprintf('   x 8  (ISR family)    = %.4g\n', min(p_p*8, 1));
fprintf('   x 16 (ISR + BBN)     = %.4g\n', min(p_p*16,1));

% ---- 10) Diagnostic figure ----
figure('Color','w','Units','inches','Position',[1 1 10 3.2]);

subplot(1,3,1);
stem(1:n, cd, 'filled', 'MarkerSize', 4); hold on;
yline(thr, 'r--', 'LineWidth', 1.2);
xlabel('Observation'); ylabel('Cook''s distance');
title('Influence'); box on; grid on;

subplot(1,3,2);
scatter(x, y, 55, 'filled', 'MarkerFaceColor',[0.2 0.6 0.2], ...
        'MarkerEdgeColor','k'); hold on;
xf = linspace(min(x), max(x), 100);
plot(xf, polyval(b_ols, xf), 'r-',  'LineWidth', 2);
plot(xf, b_rob(1) + b_rob(2)*xf,  'b--', 'LineWidth', 2);
legend({'data','OLS','Robust'}, 'Location','best', 'FontSize', 8);
xlabel('Late ISR HbO (\muM)'); ylabel('\DeltaTFI');
title('OLS vs robust fit'); box on; grid on;

subplot(1,3,3);
histogram(bootr, 40); hold on;
xline(0, 'r--', 'LineWidth', 1.2);
xlabel('Bootstrap r'); ylabel('Count');
title('Bootstrap distribution'); box on; grid on;

fprintf('\n============ END ROBUSTNESS CHECKS ============\n');


%%  ANALYSIS 4
%  THI and TFI subscales against HbO


%  POINT 16 — BASELINE DISTRESS vs CORTICAL RESPONSE
%  Subscale and total-score analysis, early and late windows
%

%
%  SCORING (verified against source papers):
%    TFI  Meikle et al. 2012, Ear Hear 33:153-176
%         25 items, each 0-10. Subscales are consecutive blocks.
%         Score = mean of answered items x 10 -> 0-100.
%    THI  Newman et al. 1996, Arch Otolaryngol 122:143-148
%         25 items, yes=4 / sometimes=2 / no=0. Score = item sum.
%         Functional 11 (max 44), Emotional 9 (max 36),
%         Catastrophic 5 (max 20). Item 14 is EMOTIONAL.
% clc;   % disabled so earlier output stays visible


fprintf('\n\n#################################################################\n');
fprintf('# POINT 16 — BASELINE DISTRESS vs HbO\n');
fprintf('#################################################################\n');


% CONFIGURATION

hbo_set = {'early_ISR','early_BBN','late_ISR','late_BBN'};

% Which correction populates the single 'p_adj' column.
% Both are computed and printed regardless.
%   'bonf' = consistent with the rest of the manuscript
%   'fdr'  = less conservative, appropriate for correlated tests
CORRECTION = 'bonf';


% 1) DEFINE SUBSCALES

TFI_names = {'Intrusive','Control','Cognitive','Sleep', ...
             'Auditory','Relaxation','QoL','Emotional'};
TFI_items = {[1 2 3],[4 5 6],[7 8 9],[10 11 12], ...
             [13 14 15],[16 17 18],[19 20 21 22],[23 24 25]};

THI_names = {'Functional','Emotional','Catastrophic'};
THI_items = {[1 2 4 7 9 12 13 15 18 20 24], ...
             [3 6 10 14 16 17 21 22 25], ...
             [5 8 11 19 23]};
THI_max   = [44 36 20];

assert(isequal(sort([TFI_items{:}]), 1:25), 'TFI item assignment error');
assert(isequal(sort([THI_items{:}]), 1:25), 'THI item assignment error');
assert(isequal(cellfun(@numel, THI_items), [11 9 5]), 'THI subscale sizes wrong');
assert(isequal(cellfun(@numel, TFI_items), [3 3 3 3 3 3 4 3]), 'TFI subscale sizes wrong');
assert(all(cellfun(@numel, THI_items) * 4 == THI_max), 'THI maxima inconsistent');
fprintf('\nItem assignments verified: TFI 3/3/3/3/3/3/4/3, THI 11/9/5.\n');


% 2) DATA INTEGRITY CHECKS

tfiPre  = cellstr(compose('tfi_%d',    1:25));
tfiPost = cellstr(compose('tfi_%d_v2', 1:25));
thiPre  = cellstr(compose('thi_%d',    1:25));
thiPost = cellstr(compose('thi_%d_v2', 1:25));

for c = [tfiPre thiPre]
    assert(ismember(c{1}, preTbl.Properties.VariableNames), ...
        'preTbl is missing column %s', c{1});
end
for c = [tfiPost thiPost]
    assert(ismember(c{1}, postTbl.Properties.VariableNames), ...
        'postTbl is missing column %s', c{1});
end
assert(numel(ids) == height(preTbl), ...
    'ids has %d entries but preTbl has %d rows', numel(ids), height(preTbl));

nFixed = 0;

V = preTbl{:, tfiPre};
bad = (V > 10 | V < 0) & ~isnan(V);
bad(all(V == 11, 2), :) = false;
nFixed = nFixed + reportBad(bad, V, tfiPre, 'pre TFI');
V(bad) = NaN;  preTbl{:, tfiPre} = V;

V = preTbl{:, thiPre};
bad = ~ismember(V, [0 2 4]) & ~isnan(V);
bad(all(V == 1, 2), :) = false;
nFixed = nFixed + reportBad(bad, V, thiPre, 'pre THI');
V(bad) = NaN;  preTbl{:, thiPre} = V;

V = postTbl{:, tfiPost};
bad = (V > 10 | V < 0) & ~isnan(V);
bad(all(V == 11, 2), :) = false;
nFixed = nFixed + reportBad(bad, V, tfiPost, 'post TFI');
V(bad) = NaN;  postTbl{:, tfiPost} = V;

V = postTbl{:, thiPost};
bad = ~ismember(V, [0 2 4]) & ~isnan(V);
bad(all(V == 1, 2), :) = false;
nFixed = nFixed + reportBad(bad, V, thiPost, 'post THI');
V(bad) = NaN;  postTbl{:, thiPost} = V;

if nFixed == 0
    fprintf('No out-of-range responses found.\n');
else
    fprintf('Blanked %d out-of-range response(s) in total.\n', nFixed);
end

Xk = preTbl{keepIdx, tfiPre};
Tk = preTbl{keepIdx, thiPre};
nSent = sum(all(Xk == 11, 2)) + sum(all(Tk == 1, 2));
assert(nSent == 0, ...
    ['%d analysed participant(s) carry sentinel questionnaire values ' ...
     '(all TFI = 11 / all THI = 1). Fix keepIdx before scoring.'], nSent);

mxTFI = max(Xk(:), [], 'omitnan');
assert(mxTFI <= 10, ...
    ['TFI items reach %g in the analysed rows. Values > 10 mean the ' ...
     'percentage items (1, 3) were not recoded to 0-10.'], mxTFI);

uT = unique(Tk(~isnan(Tk)));
assert(all(ismember(uT, [0 2 4])), ...
    'THI responses include values other than 0/2/4: %s', mat2str(uT'));

fprintf('Integrity checks passed on %d analysed participants.\n', sum(keepIdx));

% 3) COMPUTE SUBSCALE AND TOTAL SCORES

Sub = table(string(ids(:)), 'VariableNames', {'SubjectID'});

for s = 1:numel(TFI_names)
    it = TFI_items{s};
    pre_v  = scoreScale(preTbl{:,  cellstr(compose('tfi_%d',    it))}, 'mean');
    post_v = scoreScale(postTbl{:, cellstr(compose('tfi_%d_v2', it))}, 'mean');
    Sub.(['TFI_' TFI_names{s} '_pre'])  = pre_v;
    Sub.(['TFI_' TFI_names{s} '_post']) = post_v;
    Sub.(['TFI_' TFI_names{s} '_d'])    = post_v - pre_v;
end

pre_v  = scoreScale(preTbl{:,  tfiPre},  'mean');
post_v = scoreScale(postTbl{:, tfiPost}, 'mean');
Sub.TFI_Total_pre  = pre_v;
Sub.TFI_Total_post = post_v;
Sub.TFI_Total_d    = post_v - pre_v;

for s = 1:numel(THI_names)
    it = THI_items{s};
    pre_v  = scoreScale(preTbl{:,  cellstr(compose('thi_%d',    it))}, 'sum');
    post_v = scoreScale(postTbl{:, cellstr(compose('thi_%d_v2', it))}, 'sum');
    Sub.(['THI_' THI_names{s} '_pre'])  = pre_v;
    Sub.(['THI_' THI_names{s} '_post']) = post_v;
    Sub.(['THI_' THI_names{s} '_d'])    = post_v - pre_v;
end

pre_v  = scoreScale(preTbl{:,  thiPre},  'sum');
post_v = scoreScale(postTbl{:, thiPost}, 'sum');
Sub.THI_Total_pre  = pre_v;
Sub.THI_Total_post = post_v;
Sub.THI_Total_d    = post_v - pre_v;

Sub = Sub(keepIdx, :);
fprintf('Score table built for %d participants.\n', height(Sub));

nPror = sum(sum(isnan(preTbl{keepIdx, [tfiPre thiPre]})));
fprintf('Missing pre-session items across analysed rows: %d\n', nPror);


% 4) BASELINE DESCRIPTIVES

fprintf('\n--- BASELINE (PRE-SESSION) SCORES ---\n');
fprintf('%-24s %5s %4s %8s %8s %16s\n', 'Measure','max','n','mean','SD','range');
fprintf('%s\n', repmat('-',1,70));

descRows = [ ...
    cellfun(@(x) ['TFI_' x], [TFI_names {'Total'}], 'uni', 0), ...
    cellfun(@(x) ['THI_' x], [THI_names {'Total'}], 'uni', 0)];
descMax  = [repmat(100, 1, numel(TFI_names)+1), THI_max, 100];

for k = 1:numel(descRows)
    v = Sub.([descRows{k} '_pre']);
    fprintf('%-24s %5d %4d %8.2f %8.2f   %6.1f - %6.1f\n', ...
        strrep(descRows{k},'_',' '), descMax(k), sum(isfinite(v)), ...
        mean(v,'omitnan'), std(v,'omitnan'), min(v), max(v));
end


% 5) MERGE WITH HbO DATA

SubMerged = innerjoin(MergedTbl, Sub, 'Keys', 'SubjectID');
fprintf('\nMerged with HbO data: %d participants\n', height(SubMerged));

missHbO = setdiff(hbo_set, SubMerged.Properties.VariableNames);
if ~isempty(missHbO)
    warning('HbO column(s) not found and will be skipped: %s', strjoin(missHbO, ', '));
    hbo_set = intersect(hbo_set, SubMerged.Properties.VariableNames, 'stable');
end


% 6) CORRELATIONS

R = table();

for m = {'TFI_Total','THI_Total'}
    for h = 1:numel(hbo_set)
        R = [R; corrRow(SubMerged, [m{1} '_pre'], hbo_set{h}, 'Total', 'Totals')]; %#ok<AGROW>
    end
end
for s = 1:numel(TFI_names)
    for h = 1:numel(hbo_set)
        R = [R; corrRow(SubMerged, ['TFI_' TFI_names{s} '_pre'], hbo_set{h}, ...
                        'TFI subscale', ['TFIsub_' hbo_set{h}])]; %#ok<AGROW>
    end
end
for s = 1:numel(THI_names)
    for h = 1:numel(hbo_set)
        R = [R; corrRow(SubMerged, ['THI_' THI_names{s} '_pre'], hbo_set{h}, ...
                        'THI subscale', ['THIsub_' hbo_set{h}])]; %#ok<AGROW>
    end
end

% -- compute BOTH corrections within each family --
R.nFamily    = nan(height(R),1);
R.p_fdr      = nan(height(R),1);
R.p_bonf     = nan(height(R),1);
R.p_bonf_rho = nan(height(R),1);

fams = unique(R.Family);
for f = 1:numel(fams)
    m  = strcmp(R.Family, fams{f});
    nf = sum(m);
    R.nFamily(m)    = nf;
    R.p_fdr(m)      = bhFDR(R.p_pearson(m));
    R.p_bonf(m)     = min(R.p_pearson(m)  * nf, 1);
    R.p_bonf_rho(m) = min(R.p_spearman(m) * nf, 1);
end

switch lower(CORRECTION)
    case 'fdr',  R.p_adj = R.p_fdr;
    case 'bonf', R.p_adj = R.p_bonf;
    otherwise,   error('Unknown CORRECTION: %s', CORRECTION);
end

fprintf('\nFamily structure (n tests per family):\n');
for f = 1:numel(fams)
    fprintf('  %-18s %d\n', fams{f}, sum(strcmp(R.Family, fams{f})));
end


% 7) PRINT

fprintf('\n=================================================================\n');
fprintf(' BASELINE DISTRESS vs HbO  (both corrections shown)\n');
fprintf(' p_adj column follows CORRECTION = %s\n', upper(CORRECTION));
fprintf('=================================================================\n');

blocks = {'Total','TFI subscale','THI subscale'};
for b = 1:numel(blocks)
    rows = strcmp(R.Kind, blocks{b});
    if ~any(rows), continue; end
    fprintf('\n%s\n%s\n', upper(blocks{b}), repmat('=',1,104));
    fprintf('%-20s %-11s %4s %5s %8s %8s %8s %8s %8s %11s\n', ...
        'Measure','HbO','n','nFam','r','p','p_fdr','p_bonf','rho','p_bonf_rho');
    fprintf('%s\n', repmat('-',1,104));
    idx = find(rows);
    for k = idx'
        fprintf('%-20s %-11s %4d %5d %+8.3f %8.4f %8.4f %8.4f %+8.3f %11.4f%s\n', ...
            strrep(strrep(R.Measure{k},'_pre',''),'_',' '), R.HbO{k}, R.n(k), ...
            R.nFamily(k), R.r_pearson(k), R.p_pearson(k), R.p_fdr(k), ...
            R.p_bonf(k), R.rho_spearman(k), R.p_bonf_rho(k), sigStar(R.p_pearson(k)));
    end
end
fprintf('\n(* uncorrected p < .05, ** < .01, *** < .001)\n');

% -- sign consistency --
fprintf('\n--- SIGN CONSISTENCY ACROSS TFI SUBSCALES ---\n');
for h = 1:numel(hbo_set)
    m = strcmp(R.Kind,'TFI subscale') & strcmp(R.HbO, hbo_set{h});
    rr = R.r_pearson(m);
    fprintf('%-11s  %d/%d positive, %d/%d negative   (median r = %+.3f)\n', ...
        hbo_set{h}, sum(rr>0), numel(rr), sum(rr<0), numel(rr), median(rr));
end

% -- Pearson / Spearman divergence --
fprintf('\n--- ROBUSTNESS: where Pearson and Spearman disagree ---\n');
flag = xor(R.p_pearson < 0.05, R.p_spearman < 0.05);
if ~any(flag)
    fprintf('  none (no result depends on the choice of coefficient)\n');
else
    disp(R(flag, {'Measure','HbO','r_pearson','p_pearson','rho_spearman','p_spearman'}));
end

% -- what survives, under each policy --
fprintf('\n--- SURVIVES FDR (p_fdr < .05) ---\n');
s1 = R(R.p_fdr < 0.05, {'Measure','HbO','n','nFamily','r_pearson','p_pearson','p_fdr'});
if isempty(s1), fprintf('  none\n'); else, disp(s1); end

fprintf('\n--- SURVIVES BONFERRONI, PEARSON (p_bonf < .05) ---\n');
s2 = R(R.p_bonf < 0.05, {'Measure','HbO','n','nFamily','r_pearson','p_pearson','p_bonf'});
if isempty(s2), fprintf('  none\n'); else, disp(s2); end

fprintf('\n--- SURVIVES BONFERRONI, SPEARMAN (p_bonf_rho < .05) ---\n');
s3 = R(R.p_bonf_rho < 0.05, {'Measure','HbO','n','nFamily','rho_spearman','p_spearman','p_bonf_rho'});
if isempty(s3), fprintf('  none\n'); else, disp(s3); end

fprintf('\n--- NOMINALLY SIGNIFICANT (uncorrected p < .05) ---\n');
nom = R(R.p_pearson < 0.05, {'Kind','Measure','HbO','r_pearson','p_pearson','p_fdr','p_bonf'});
if isempty(nom), fprintf('  none\n'); else, disp(nom); end


% 8) SAVE
%
% writetable(Sub,       fullfile(RESULTS_DIR, 'Point16_scores.csv'));
% writetable(SubMerged, fullfile(RESULTS_DIR, 'Point16_scores_with_HbO.csv'));
% writetable(R,         fullfile(RESULTS_DIR, 'Point16_correlations.csv'));
% assignin('base','ScoreTbl',  Sub);
% assignin('base','SubMerged', SubMerged);
% assignin('base','P16',       R);
% fprintf('\nSaved: Point16_scores.csv, Point16_scores_with_HbO.csv,\n');
% fprintf('       Point16_correlations.csv\n');
% fprintf('#################################################################\n');


%%  ANALYSIS 5a
%  Connectivity adjusted for age and PTA-4


%  POINT 13 — RSFC ADJUSTED FOR AGE + PTA
%  Requires: AudTbl (from the Point 4 block)
%            T_CorrData_TINN / HR.Control.ROI  (row order)
%            zPre_exp/zPost_exp/zDiff_exp       (ROI seed, tinnitus)
%            zPre_ctrl/zPost_ctrl/zDiff_ctrl    (ROI seed, control)
%            Zlob_pre_exp/post/diff  [27 x 9]   (lobe-wise, tinnitus)
%            Zlob_pre_ctrl/post/diff [27 x 9]   (lobe-wise, control)
%            lobes                              (1 x 9 names)
% =========================================================

fprintf('\n\n#################################################################\n');
fprintf('# POINT 13 — RSFC BETWEEN GROUPS, ADJUSTED FOR AGE + PTA\n');
fprintf('#################################################################\n');

assert(exist('AudTbl','var')==1, 'Run the Point 4 block first (AudTbl missing).');

nrm13 = @(s) upper(string(regexprep(cellstr(string(s)), '\s+', '')));

idT = nrm13(T_CorrData_TINN.SubjectID);
idC = nrm13(HR.Control.ROI.SubjectID);

% ---- attach covariates in array row order ----
[tfT, locT] = ismember(idT, AudTbl.SubjectID);
[tfC, locC] = ismember(idC, AudTbl.SubjectID);
assert(all(tfT), 'Unmatched tinnitus IDs: %s', strjoin(cellstr(idT(~tfT)),', '));
assert(all(tfC), 'Unmatched control IDs: %s',  strjoin(cellstr(idC(~tfC)),', '));

ptaT = AudTbl.PTA(locT);  ageT = AudTbl.Age(locT);
ptaC = AudTbl.PTA(locC);  ageC = AudTbl.Age(locC);

fprintf('\nCovariates attached: %d tinnitus, %d control (all matched)\n', ...
    numel(ptaT), numel(ptaC));
fprintf('Participants without PTA-4 (dropped from adjusted models): %d tinnitus, %d control\n', ...
    sum(~isfinite(ptaT)), sum(~isfinite(ptaC)));

PTA = [ptaT; ptaC];
AGE = [ageT; ageC];
GRP = [ones(numel(ptaT),1); zeros(numel(ptaC),1)];

% ---- 1) ROI-SEED POOLED ----

fprintf(' 1) ROI-SEED POOLED CONNECTIVITY\n');

R13 = table();
R13 = [R13; bg13(zPre_exp,  zPre_ctrl,  GRP, PTA, AGE, 'ROI seed PRE')];
R13 = [R13; bg13(zPost_exp, zPost_ctrl, GRP, PTA, AGE, 'ROI seed POST')];
R13 = [R13; bg13(zDiff_exp, zDiff_ctrl, GRP, PTA, AGE, 'ROI seed dZ')];
disp(R13);

% ---- 2) LOBE-WISE, THREE WINDOWS ----

fprintf(' 2) LOBE-WISE CONNECTIVITY (Bonferroni x9 within each window)\n');

winNames = {'PRE','POST','dZ'};
Eexp = {Zlob_pre_exp,  Zlob_post_exp,  Zlob_diff_exp};
Ectl = {Zlob_pre_ctrl, Zlob_post_ctrl, Zlob_diff_ctrl};

L13new = table();
for w = 1:3
    fprintf('\n--- %s ---\n', winNames{w});
    fprintf('%-5s %10s %10s | %11s %11s %8s %8s\n', ...
        'Lobe','t_unadj','p_unadj','b_group','p_adj','p_PTA','p_Age');
    fprintf('%s\n', repmat('-',1,72));
    for L = 1:numel(lobes)
        row = bg13(Eexp{w}(:,L), Ectl{w}(:,L), GRP, PTA, AGE, ...
                   sprintf('%s %s', winNames{w}, lobes{L}));
        row.Window = string(winNames{w});
        row.Lobe   = string(lobes{L});
        row.p_unadj_bonf    = min(row.p_unadj*9, 1);
        row.p_group_adj_bonf = min(row.p_group_adj*9, 1);
        fprintf('%-5s %+10.3f %10.4f | %+11.4f %11.4f %8.3f %8.3f%s\n', ...
            lobes{L}, row.t_unadj, row.p_unadj, row.b_group, row.p_group_adj, ...
            row.p_PTA, row.p_Age, p13s(row.p_group_adj));
        L13new = [L13new; row]; %#ok<AGROW>
    end
end

% ---- 3) VERDICT ----

fprintf(' 3) VERDICT\n');


All13 = [R13; L13new(:, R13.Properties.VariableNames)];

flip13 = All13(xor(All13.p_unadj<0.05, All13.p_group_adj<0.05), ...
               {'Measure','p_unadj','p_group_adj'});
if isempty(flip13)
    fprintf('\nNo RSFC effect changed significance after adjusting for age and PTA.\n');
else
    fprintf('\nEffects that changed significance after adjustment:\n');
    disp(flip13);
end

fprintf('\nPTA significant in %d of %d models.\n', sum(All13.p_PTA<0.05), height(All13));
fprintf('Age significant in %d of %d models.\n', sum(All13.p_Age<0.05), height(All13));

fprintf('\nSurviving Bonferroni x9 after adjustment (lobe-wise only):\n');
srv13 = L13new(L13new.p_group_adj_bonf < 0.05, ...
               {'Window','Lobe','p_unadj','p_group_adj','p_group_adj_bonf'});
if isempty(srv13)
    fprintf('  none\n');
else
    disp(srv13);
end

writetable(R13,    fullfile(RESULTS_DIR, 'Point13_RSFC_ROIseed_adjusted.csv'));
writetable(L13new, fullfile(RESULTS_DIR, 'Point13_RSFC_lobewise_adjusted.csv'));
assignin('base','R13',R13);
assignin('base','L13new',L13new);
fprintf('\nSaved: Point13_RSFC_ROIseed_adjusted.csv, Point13_RSFC_lobewise_adjusted.csv\n');
fprintf('#################################################################\n');


%%  ANALYSIS 5b
%  Hemodynamic response adjusted for age

%  POINT 13: AGE AS COVARIATE — BETWEEN-GROUP HEMODYNAMIC RESPONSE

%  Requires in workspace:
%    tEpoch
%    CTRL_ROI_BBN,  CTRL_ROI_ISR,  CTRL_NROI_BBN,  CTRL_NROI_ISR
%    TINN_ROI_BBN,  TINN_ROI_ISR,  TINN_NROI_BBN,  TINN_NROI_ISR


fprintf(' BETWEEN-GROUP HbO WITH AGE AS COVARIATE (ANCOVA)\n');


% ---------------------------------------------------------
% 1) AGE LOOKUP TABLES (ID -> age)
% ---------------------------------------------------------
% Tinnitus ages (same master list used in the RSFC script)
tinnMapIDs = ["Tinnitus01";"Tinnitus12";"Tinnitus13";"Tinnitus18";"Tinnitus19"; ...
              "Tinnitus21";"Tinnitus23";"Tinnitus28";"Tinnitus32";"Tinnitus34"; ...
              "Tinnitus02";"Tinnitus38";"Tinnitus39";"Tinnitus42";"Tinnitus46"; ...
              "Tinnitus47";"Tinnitus59";"Tinnitus61";"Tinnitus62";"Tinnitus69"; ...
              "Tinnitus03";"Tinnitus70";"Tinnitus04";"Tinnitus05";"Tinnitus06"; ...
              "Tinnitus07";"Tinnitus08"];
tinnMapAge = [45;63;63;70;52;63;61;65;69;72;49;59;70;68;65;68;66;63;65;60; ...
              62;64;46;64;54;65;53];
tinnAgeMap = containers.Map(cellstr(tinnMapIDs), num2cell(tinnMapAge));



% 2) SUBJECT ID LISTS (from HR_outputs.mat, in curve column order)

ids_ctrl_roi  = HR.Control.SubjectID;
ids_ctrl_nroi = HR.Control.SubjectID;
ids_tinn_roi  = HR.Tinnitus.SubjectID;
ids_tinn_nroi = HR.Tinnitus.SubjectID;

fprintf('\nColumn-count checks:\n');
fprintf('  CTRL ROI  : %2d IDs vs %2d columns\n', numel(ids_ctrl_roi),  size(CTRL_ROI_BBN,2));
fprintf('  CTRL nROI : %2d IDs vs %2d columns\n', numel(ids_ctrl_nroi), size(CTRL_NROI_BBN,2));
fprintf('  TINN ROI  : %2d IDs vs %2d columns\n', numel(ids_tinn_roi),  size(TINN_ROI_BBN,2));
fprintf('  TINN nROI : %2d IDs vs %2d columns\n', numel(ids_tinn_nroi), size(TINN_NROI_BBN,2));

assert(numel(ids_ctrl_roi)  == size(CTRL_ROI_BBN,2)  && ...
       numel(ids_ctrl_nroi) == size(CTRL_NROI_BBN,2) && ...
       numel(ids_tinn_roi)  == size(TINN_ROI_BBN,2)  && ...
       numel(ids_tinn_nroi) == size(TINN_NROI_BBN,2), ...
    ['Subject ID count does not match curve column count. The ROI and ' ...
     'non-ROI channel gating must have excluded different participants, ' ...
     'so the ID lists cannot be shared.']);

% Control ages, in the order of HR.Control.SubjectID
age_ctrl_ordered = [45,41,33,41,39,35,30,69,62,62,60,62,69,53,69,62,52,69, ...
                    62,60,65,51,59,59,62,50,59]';
assert(numel(age_ctrl_ordered) == numel(ids_ctrl_roi), ...
    'age_ctrl_ordered (%d) does not match control n (%d)', ...
    numel(age_ctrl_ordered), numel(ids_ctrl_roi));
ctrlAgeMap = containers.Map(cellstr(ids_ctrl_roi), num2cell(age_ctrl_ordered));

% 3) RUN THE FOUR BETWEEN-GROUP TESTS, UNADJUSTED AND AGE-ADJUSTED

windows = {'EARLY', [5 12]; 'LATE', [13 20]};
regions = {
    'ROI',    CTRL_ROI_BBN,  CTRL_ROI_ISR,  ids_ctrl_roi,  TINN_ROI_BBN,  TINN_ROI_ISR,  ids_tinn_roi;
    'nonROI', CTRL_NROI_BBN, CTRL_NROI_ISR, ids_ctrl_nroi, TINN_NROI_BBN, TINN_NROI_ISR, ids_tinn_nroi;
};

hboResults = table();

for w = 1:size(windows,1)
    wname = windows{w,1};
    win   = windows{w,2};

    fprintf('\n-----------------------------------------------------------------\n');
    fprintf('%s WINDOW [%g-%g s] — between-group delta (BBN - ISR)\n', wname, win(1), win(2));
    fprintf('-----------------------------------------------------------------\n');

    for r = 1:size(regions,1)
        rname = regions{r,1};

        [dC, aC] = localDeltaAndAge(regions{r,2}, regions{r,3}, regions{r,4}, ...
                                    tEpoch, win, ctrlAgeMap);
        [dT, aT] = localDeltaAndAge(regions{r,5}, regions{r,6}, regions{r,7}, ...
                                    tEpoch, win, tinnAgeMap);

        if numel(dC) < 3 || numel(dT) < 3
            fprintf('%-8s | insufficient data (nC=%d, nT=%d)\n', rname, numel(dC), numel(dT));
            continue;
        end

        % --- Unadjusted Welch t-test (as currently reported) ---
        [~, p_raw, ~, s_raw] = ttest2(dC, dT, 'Vartype','unequal');

        % --- Age-adjusted linear model: delta ~ group + age ---
        y   = [dC; dT];
        grp = [zeros(numel(dC),1); ones(numel(dT),1)];   % 0 = control, 1 = tinnitus
        age = [aC; aT];
        mdl = fitlm([grp age], y, 'VarNames', {'Group','Age','DeltaHbO'});

        t_adj = mdl.Coefficients.tStat( strcmp(mdl.CoefficientNames,'Group'));
        p_adj = mdl.Coefficients.pValue(strcmp(mdl.CoefficientNames,'Group'));
        p_age = mdl.Coefficients.pValue(strcmp(mdl.CoefficientNames,'Age'));

        fprintf('%-8s | nC=%2d nT=%2d | CTRL mean=%+.4f  TINN mean=%+.4f\n', ...
            rname, numel(dC), numel(dT), mean(dC), mean(dT));
        fprintf('%-8s |   unadjusted: t=%+6.3f  p=%7.4f%s\n', ...
            '', s_raw.tstat, p_raw, localStar(p_raw));
        fprintf('%-8s |   age-adjusted: t=%+6.3f  p=%7.4f%s   (age p=%.4f)\n', ...
            '', t_adj, p_adj, localStar(p_adj), p_age);

        hboResults = [hboResults; table({wname}, {rname}, ...
            numel(dC), numel(dT), mean(dC), mean(dT), ...
            s_raw.tstat, p_raw, t_adj, p_adj, p_age, ...
            'VariableNames', {'Window','Region','n_ctrl','n_tinn', ...
                              'mean_ctrl','mean_tinn','t_raw','p_raw', ...
                              't_adj','p_adj','p_age'})]; %#ok<AGROW>
    end
end


% 4) SUMMARY

fprintf('\n=================================================================\n');
fprintf(' SUMMARY\n');
fprintf('=================================================================\n');
disp(hboResults);

changed = hboResults( (hboResults.p_raw < 0.05 & hboResults.p_adj >= 0.05) | ...
                      (hboResults.p_raw >= 0.05 & hboResults.p_adj < 0.05), :);
if isempty(changed)
    fprintf('\nNo between-group HbO effect changed significance status after\n');
    fprintf('adjusting for age.\n');
else
    fprintf('\nChanged significance status after age adjustment:\n');
    disp(changed(:, {'Window','Region','t_raw','p_raw','t_adj','p_adj'}));
end

if all(hboResults.p_age >= 0.05)
    fprintf('\nAge was not a significant predictor in any model (all p >= 0.05).\n');
end

assignin('base','hbo_age_adjusted', hboResults);
fprintf('\nTable saved to workspace as ''hbo_age_adjusted''.\n');
fprintf('=================================================================\n');



%%  ANALYSIS 5c
%  Sex comparison, and change-in-TFI adjusted for age and PTA-4

[~, p_sex] = fishertest([13 14; 14 13]);
fprintf('Sex, Fisher exact p = %.4f\n', p_sex);

%%  AGE-ADJUSTED dTFI CORRELATION
%  Requires: Dt (from the Point 4 block) with late_ISR, Delta_TFI, Age


fprintf('\n\n#################################################################\n');
fprintf('# dTFI vs LATE ISR HbO, ADJUSTED FOR AGE\n');
fprintf('#################################################################\n');

y = Dt.late_ISR;
x = Dt.Delta_TFI;
a = Dt.Age;
ok = isfinite(y) & isfinite(x) & isfinite(a);
fprintf('\nn = %d\n', sum(ok));

% ---- unadjusted ----
[r,p]      = corr(x(ok), y(ok), 'Type','Pearson');
[rho,prho] = corr(x(ok), y(ok), 'Type','Spearman');
fprintf('\nUnadjusted:  r = %+.3f (p = %.4f)   rho = %+.3f (p = %.4f)\n', r, p, rho, prho);

% ---- partial correlation, controlling age ----
[rpart, ppart] = partialcorr(x(ok), y(ok), a(ok));
fprintf('Partial (age controlled):  r = %+.3f (p = %.4f)\n', rpart, ppart);

% ---- multiple regression ----
mdl = fitlm([x(ok), a(ok)-mean(a(ok))], y(ok), 'VarNames',{'dTFI','Age','late_ISR'});
fprintf('\nModel: late_ISR ~ dTFI + Age\n');
disp(mdl.Coefficients);
fprintf('Adjusted R-squared = %.3f\n', mdl.Rsquared.Adjusted);

% ---- is age itself related to either variable? ----
[ra,pa] = corr(a(ok), y(ok));   fprintf('\nAge vs late_ISR:  r = %+.3f (p = %.4f)\n', ra, pa);
[rb,pb] = corr(a(ok), x(ok));   fprintf('Age vs dTFI:      r = %+.3f (p = %.4f)\n', rb, pb);

fprintf('\n#################################################################\n');

%% ---- dTFI vs LATE ISR HbO, ADJUSTED FOR AGE + PTA ----
y = Dt.late_ISR;  x = Dt.Delta_TFI;  a = Dt.Age;  pt = Dt.PTA;
ok = isfinite(y) & isfinite(x) & isfinite(a) & isfinite(pt);
fprintf('\nn = %d\n', sum(ok));

[r,p] = corr(x(ok), y(ok));
fprintf('Unadjusted:  r = %+.3f (p = %.4f)\n', r, p);

[rp2, pp2] = partialcorr(x(ok), y(ok), [a(ok), pt(ok)]);
fprintf('Partial (age + PTA controlled):  r = %+.3f (p = %.4f)\n', rp2, pp2);

mdl2 = fitlm([x(ok), a(ok)-mean(a(ok)), pt(ok)-mean(pt(ok))], y(ok), ...
             'VarNames',{'dTFI','Age','PTA','late_ISR'});
fprintf('\nModel: late_ISR ~ dTFI + Age + PTA\n');
disp(mdl2.Coefficients);
fprintf('Adjusted R-squared = %.3f\n', mdl2.Rsquared.Adjusted);

[rc,pc] = corr(pt(ok), y(ok));  fprintf('\nPTA vs late_ISR:  r = %+.3f (p = %.4f)\n', rc, pc);
[rd,pd] = corr(pt(ok), x(ok));  fprintf('PTA vs dTFI:      r = %+.3f (p = %.4f)\n', rd, pd);


%%  ANALYSIS 6
%  Distributions of the pre-session THI and TFI scores
% ---- DISTRIBUTION STATS FOR POINT 15 ----
fprintf('\n================ DISTRIBUTIONS (pre-session) ================\n');

for k = 1:2
    if k==1
        v = THI_pre; nm = 'THI';
        edges = [0 16; 17 36; 37 56; 57 76; 77 100];
        labs  = {'slight','mild','moderate','severe','catastrophic'};
    else
        v = TFI_pre; nm = 'TFI';
        edges = [0 17; 17.01 31; 31.01 53; 53.01 72; 72.01 100];
        labs  = {'not a problem','small','moderate','big','very big'};
    end

    fprintf('\n--- %s (n=%d) ---\n', nm, numel(v));
    fprintf('  mean   = %.2f  (SD %.2f)\n', mean(v), std(v));
    fprintf('  median = %.2f  (IQR %.1f - %.1f)\n', median(v), prctile(v,25), prctile(v,75));
    fprintf('  range  = %.1f - %.1f\n', min(v), max(v));
    fprintf('  skewness = %.2f\n', skewness(v));
    fprintf('  severity bands:\n');
    for b = 1:5
        n = sum(v >= edges(b,1) & v <= edges(b,2));
        fprintf('    %-16s n=%2d (%4.1f%%)\n', labs{b}, n, 100*n/numel(v));
    end
end


%%  ANALYSIS 7
%  Baseline distress against HbO, and change against change
%  POINT 16: TWO ADDITIONAL ANALYSES REQUESTED BY REVIEWER 2
%
%  Analysis 1: baseline THI/TFI vs noise-evoked (BBN) HbO   [NEW]
%  Analysis 2: individual RSFC change vs individual score change
%              [already computed; reported explicitly here]


fprintf('\n\n');

fprintf('# POINT 16 — ADDITIONAL ANALYSES\n');



% ANALYSIS 1: BASELINE DISTRESS vs BBN-EVOKED HbO
% Does pre-session distress predict the cortical response
% to a neutral broadband stimulus?

fprintf('\n=================================================================\n');
fprintf(' ANALYSIS 1: baseline THI / TFI  vs  BBN-evoked HbO (ROI)\n');
fprintf('=================================================================\n');
fprintf('%-14s %-8s %-9s %4s %8s %9s\n', 'Score','Window','Test','n','stat','p');
fprintf('%s\n', repmat('-',1,58));

base_pairs = {
    MergedTbl.early_BBN, MergedTbl.THI_Pre, 'THI baseline', 'early';
    MergedTbl.late_BBN,  MergedTbl.THI_Pre, 'THI baseline', 'late';
    MergedTbl.early_BBN, MergedTbl.TFI_Pre, 'TFI baseline', 'early';
    MergedTbl.late_BBN,  MergedTbl.TFI_Pre, 'TFI baseline', 'late';
};

A1 = table();
for k = 1:size(base_pairs,1)
    xa = base_pairs{k,1};  ya = base_pairs{k,2};
    lab = base_pairs{k,3}; win = base_pairs{k,4};
    tType = testForScore(lab);                    % THI -> Spearman, TFI -> Pearson
    ok = isfinite(xa) & isfinite(ya);
    [rv, pv] = corr(xa(ok), ya(ok), 'Type', tType);
    if strcmp(tType,'Spearman'), st='rho'; else, st='r  '; end
    fprintf('%-14s %-8s %-9s %4d  %s=%+.3f %8.4f%s\n', ...
        lab, win, tType, sum(ok), st, rv, pv, p16star(pv));
    A1 = [A1; table({lab},{win},{tType},sum(ok),rv,pv, ...
        'VariableNames',{'Score','Window','Test','n','stat','p'})]; %#ok<AGROW>
end

% Bonferroni within this new exploratory family of 4
A1.p_bonf = min(A1.p * height(A1), 1);
fprintf('\nBonferroni within this family of %d tests:\n', height(A1));
for k = 1:height(A1)
    fprintf('  %-14s %-6s  p_uncorr=%.4f  p_bonf=%.4f%s\n', ...
        A1.Score{k}, A1.Window{k}, A1.p(k), A1.p_bonf(k), p16star(A1.p_bonf(k)));
end

% Optional: also test baseline scores against ISR HbO for completeness
fprintf('\n--- (supplementary) baseline scores vs ISR HbO ---\n');
isr_base = {
    MergedTbl.early_ISR, MergedTbl.THI_Pre, 'THI baseline', 'early';
    MergedTbl.late_ISR,  MergedTbl.THI_Pre, 'THI baseline', 'late';
    MergedTbl.early_ISR, MergedTbl.TFI_Pre, 'TFI baseline', 'early';
    MergedTbl.late_ISR,  MergedTbl.TFI_Pre, 'TFI baseline', 'late';
};
for k = 1:size(isr_base,1)
    xa = isr_base{k,1}; ya = isr_base{k,2};
    lab = isr_base{k,3}; win = isr_base{k,4};
    tType = testForScore(lab);
    ok = isfinite(xa) & isfinite(ya);
    [rv, pv] = corr(xa(ok), ya(ok), 'Type', tType);
    fprintf('  %-14s %-6s %-9s n=%2d  stat=%+.3f  p=%.4f%s\n', ...
        lab, win, tType, sum(ok), rv, pv, p16star(pv));
end


% ANALYSIS 2: INDIVIDUAL RSFC CHANGE vs INDIVIDUAL SCORE CHANGE

fprintf('\n=================================================================\n');
fprintf(' ANALYSIS 2: RSFC change  vs  questionnaire change (per participant)\n');
fprintf('=================================================================\n');

fprintf('\n--- ROI-seed RSFC (dZ) ---\n');
d_pairs = {
    rsfc_roi_diff, d_thi_roi, 'dTHI';
    rsfc_roi_diff, d_tfi_roi, 'dTFI';
};
A2 = table();
for k = 1:size(d_pairs,1)
    xa = d_pairs{k,1}; ya = d_pairs{k,2}; lab = d_pairs{k,3};
    tType = testForScore(lab);
    ok = isfinite(xa) & isfinite(ya);
    [rv, pv] = corr(xa(ok), ya(ok), 'Type', tType);
    fprintf('  RSFC dZ (ROI) x %-6s [%-8s] n=%2d  stat=%+.3f  p=%.4f%s\n', ...
        lab, tType, sum(ok), rv, pv, p16star(pv));
    A2 = [A2; table({'ROI-seed dZ'},{lab},{tType},sum(ok),rv,pv, ...
        'VariableNames',{'Metric','Score','Test','n','stat','p'})]; %#ok<AGROW>
end

fprintf('\n--- Whole-brain RSFC (dr) ---\n');
dwb_pairs = {
    rsfc_wb_delta, d_thi_wb, 'dTHI';
    rsfc_wb_delta, d_tfi_wb, 'dTFI';
};
for k = 1:size(dwb_pairs,1)
    xa = dwb_pairs{k,1}; ya = dwb_pairs{k,2}; lab = dwb_pairs{k,3};
    tType = testForScore(lab);
    ok = isfinite(xa) & isfinite(ya);
    [rv, pv] = corr(xa(ok), ya(ok), 'Type', tType);
    fprintf('  RSFC dr (WB) x %-6s [%-8s] n=%2d  stat=%+.3f  p=%.4f%s\n', ...
        lab, tType, sum(ok), rv, pv, p16star(pv));
    A2 = [A2; table({'Whole-brain dr'},{lab},{tType},sum(ok),rv,pv, ...
        'VariableNames',{'Metric','Score','Test','n','stat','p'})]; %#ok<AGROW>
end

fprintf('\n--- Lobe-wise RSFC (dZ) x score change ---\n');
fprintf('%-6s %-6s %8s %9s %11s\n','Lobe','Score','stat','p','p_bonf(9)');
fprintf('%s\n', repmat('-',1,46));
for sc = 1:2
    if sc==1, score = d_thi_roi; lab='dTHI'; else, score = d_tfi_roi; lab='dTFI'; end
    tType = testForScore(lab);
    for L = 1:numel(Lobes)
        xa = Zlob_diff_aligned(:,L);
        ok = isfinite(xa) & isfinite(score);
        if sum(ok) < 5, continue; end
        [rv, pv] = corr(xa(ok), score(ok), 'Type', tType);
        pb = min(pv*numel(Lobes),1);
        fprintf('%-6s %-6s %+8.3f %9.4f %11.4f%s\n', ...
            Lobes{L}, lab, rv, pv, pb, p16star(pb));
    end
end

assignin('base','point16_baseline', A1);
assignin('base','point16_deltadelta', A2);
fprintf('\nTables saved to workspace: point16_baseline, point16_deltadelta\n');


%%  ANALYSIS 8
%  Minimum detectable effects at 80% power
%  SENSITIVITY ANALYSIS — minimum detectable effect sizes
%  80% power, alpha = 0.05 two-tailed, n = 27 per group
%  Requires Statistics and Machine Learning Toolbox

fprintf('# SENSITIVITY ANALYSIS — MINIMUM DETECTABLE EFFECTS\n');
alpha = 0.05;
target_power = 0.80;
n1 = 27; n2 = 27;

% ---- 1) BETWEEN-GROUP (two independent samples) ----
df  = n1 + n2 - 2;
crit = tinv(1 - alpha/2, df);
pw_bg = @(d) 1 - nctcdf(crit, df, d*sqrt(n1*n2/(n1+n2))) ...
               + nctcdf(-crit, df, d*sqrt(n1*n2/(n1+n2)));
d_bg = fzero(@(d) pw_bg(d) - target_power, [0.05 1.5]);

% ---- 2) WITHIN-GROUP (paired) ----
dfp  = n1 - 1;
critp = tinv(1 - alpha/2, dfp);
pw_wg = @(d) 1 - nctcdf(critp, dfp, d*sqrt(n1)) ...
               + nctcdf(-critp, dfp, d*sqrt(n1));
dz_wg = fzero(@(d) pw_wg(d) - target_power, [0.05 1.5]);

% ---- 3) CORRELATION (Fisher z) ----
se = 1/sqrt(n1 - 3);
zc = norminv(1 - alpha/2);
pw_r = @(r) normcdf(abs(0.5*log((1+r)/(1-r)))/se - zc) ...
          + normcdf(-abs(0.5*log((1+r)/(1-r)))/se - zc);
r_min = fzero(@(r) pw_r(r) - target_power, [0.05 0.95]);

fprintf('\nAt %.0f%% power, alpha = %.2f (two-tailed):\n\n', target_power*100, alpha);
fprintf('  Between-group (n = %d vs %d)   d  = %.3f\n', n1, n2, d_bg);
fprintf('  Within-group paired (n = %d)    dz = %.3f\n', n1, dz_wg);
fprintf('  Correlation (n = %d)            r  = %.3f\n', n1, r_min);

% ---- achieved power for the effects actually observed ----
fprintf('\nAchieved power for observed effects:\n\n');
obs = { 'early ΔHbO between-group', 'bg', 0.601; ...
        'late ΔHbO between-group',  'bg', 0.384; ...
        'late BBN vs ISR, tinnitus','wg', 0.465; ...
        'early BBN vs ISR, control','wg', 1.163; ...
        'ΔTFI vs late ISR HbO',     'r',  0.515 };
for k = 1:size(obs,1)
    switch obs{k,2}
        case 'bg', p = pw_bg(obs{k,3});
        case 'wg', p = pw_wg(obs{k,3});
        case 'r',  p = pw_r(obs{k,3});
    end
    fprintf('  %-28s %5.3f -> power = %.2f\n', obs{k,1}, obs{k,3}, p);
end
fprintf('\n#################################################################\n');


%%  LOCAL FUNCTIONS


function s = p4bstar(p)
    if ~isfinite(p), s=''; elseif p<0.001, s=' ***';
    elseif p<0.01, s=' **'; elseif p<0.05, s=' *'; else, s=''; end
end

function s = p4bYN(c)
    if c, s = 'YES — see the p_PTA / p_Age column'; else, s = 'No (all p >= 0.05)'; end
end

function s = p6bstar(p)
    if     ~isfinite(p), s = '';
    elseif p < 0.001,    s = ' ***';
    elseif p < 0.01,     s = ' **';
    elseif p < 0.05,     s = ' *';
    else,                s = '';
    end
end

function n = reportBad(bad, V, cols, label)
    n = sum(bad(:));
    if n == 0, return; end
    [ri, ci] = find(bad);
    for j = 1:numel(ri)
        fprintf('  blanked %s row %d: %s = %g\n', ...
            label, ri(j), cols{ci(j)}, V(ri(j), ci(j)));
    end
end

function out = scoreScale(X, mode)
    nItems = size(X, 2);
    nAns   = sum(isfinite(X), 2);
    if nItems <= 3
        minAns = 2;
    else
        minAns = ceil(0.75 * nItems);
    end
    itemMean = mean(X, 2, 'omitnan');
    switch mode
        case 'mean', out = itemMean * 10;
        case 'sum',  out = itemMean * nItems;
        otherwise,   error('scoreScale: unknown mode "%s"', mode);
    end
    out(nAns < minAns) = NaN;
end

function T = corrRow(D, xname, yname, kind, family)
    xv = D.(xname);
    yv = D.(yname);
    ok = isfinite(xv) & isfinite(yv);
    if sum(ok) < 4
        [rp, pp, rs, ps] = deal(NaN, NaN, NaN, NaN);
    else
        [rp, pp] = corr(xv(ok), yv(ok), 'Type', 'Pearson');
        [rs, ps] = corr(xv(ok), yv(ok), 'Type', 'Spearman');
    end
    T = table({xname}, {yname}, {kind}, {family}, sum(ok), rp, pp, rs, ps, ...
        'VariableNames', {'Measure','HbO','Kind','Family','n', ...
                          'r_pearson','p_pearson','rho_spearman','p_spearman'});
end

function q = bhFDR(p)
    q  = nan(size(p));
    ok = isfinite(p);
    pv = p(ok);
    m  = numel(pv);
    if m == 0, return; end
    [ps, ord] = sort(pv);
    adj = ps .* m ./ (1:m)';
    adj = min(1, flipud(cummin(flipud(adj))));
    out = nan(m,1);
    out(ord) = adj;
    q(ok) = out;
end

function s = sigStar(p)
    if     ~isfinite(p), s = '';
    elseif p < 0.001,    s = ' ***';
    elseif p < 0.01,     s = ' **';
    elseif p < 0.05,     s = ' *';
    else,                s = '';
    end
end

function s = p16star(p)
    if p < 0.001,     s = ' ***';
    elseif p < 0.01,  s = ' **';
    elseif p < 0.05,  s = ' *';
    else,             s = '';
    end
end

function out = bg13(yT, yC, GRP, PTA, AGE, label)
    y  = [yT(:); yC(:)];
    ok = isfinite(y) & isfinite(PTA) & isfinite(AGE);
    [~,pu,~,su] = ttest2(y(GRP==1 & isfinite(y)), y(GRP==0 & isfinite(y)), ...
                         'Vartype','unequal');
    mdl = fitlm([GRP(ok), PTA(ok)-mean(PTA(ok)), AGE(ok)-mean(AGE(ok))], y(ok), ...
                'VarNames',{'Group','PTA','Age','Y'});
    out = table({label}, su.tstat, pu, ...
                mdl.Coefficients.Estimate(2), mdl.Coefficients.pValue(2), ...
                mdl.Coefficients.pValue(3), mdl.Coefficients.pValue(4), sum(ok), ...
        'VariableNames',{'Measure','t_unadj','p_unadj','b_group', ...
                         'p_group_adj','p_PTA','p_Age','n'});
end

function s = p13s(p)
    if ~isfinite(p), s=''; elseif p<0.001, s=' ***';
    elseif p<0.01, s=' **'; elseif p<0.05, s=' *'; else, s=''; end
end

function [delta, ages] = localDeltaAndAge(BBN, ISR, ids, tEpoch, win, ageMap)
% Subject-wise window means, delta = BBN - ISR, with matched ages.

    delta = []; ages = [];
    if isempty(BBN) || isempty(ISR) || isempty(ids), return; end

    nT = min([size(BBN,1), size(ISR,1), numel(tEpoch)]);
    t  = tEpoch(1:nT);
    iWin = find(t >= win(1) & t <= win(2));
    if isempty(iWin), return; end

    nSub = min([size(BBN,2), size(ISR,2), numel(ids)]);
    BBN  = BBN(1:nT, 1:nSub);
    ISR  = ISR(1:nT, 1:nSub);
    ids  = ids(1:nSub);

    mB = mean(BBN(iWin,:), 1, 'omitnan')';
    mI = mean(ISR(iWin,:), 1, 'omitnan')';
    d  = mB - mI;

    % Ages by ID lookup
    a = nan(nSub,1);
    for k = 1:nSub
        key = char(ids(k));
        if isKey(ageMap, key)
            a(k) = ageMap(key);
        else
            warning('No age found for %s — excluded from ANCOVA', key);
        end
    end

    ok    = isfinite(d) & isfinite(a);
    delta = d(ok);
    ages  = a(ok);
end

function s = localStar(p)
    if p < 0.001,     s = ' ***';
    elseif p < 0.01,  s = ' **';
    elseif p < 0.05,  s = ' *';
    else,             s = '';
    end
end


function t = ternaryTest(isTHI)
    if isTHI, t = 'Spearman'; else, t = 'Pearson'; end
end