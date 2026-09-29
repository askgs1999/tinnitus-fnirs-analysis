%%  RSFC ANALYSIS
%  Resting-state functional connectivity from fNIRS in tinnitus and controls
%
%  Satish G, Taylor LM, Arnold MP, Gallagher-Shale J, Krishnamurthy K,
%  Basura GJ. Subjective tinnitus distress correlates with cortical
%  hemodynamic changes in human auditory cortex as measured by functional
%  near-infrared spectroscopy. PLOS ONE.
%
%  WHAT THIS SCRIPT DOES
%    1. Preprocesses fNIRS recordings: scalp coupling index screening,
%       optical density conversion, TDDR motion correction, short-separation
%       channel regression, 0.01-0.08 Hz band-pass filtering, and the
%       modified Beer-Lambert law.
%    2. Registers the probe to the Colin27 atlas and assigns each channel an
%       AAL cortical label.
%    3. Computes resting-state functional connectivity during the pre- and
%       post-stimulation silence periods at three levels: ROI seed-based,
%       ROI-seed-to-lobe, and whole-brain channel-pair.
%    4. Runs the within-group, between-group and lobe-wise statistics
%       reported in the manuscript, with Bonferroni correction across the
%       nine cortical lobes.
%    5. Generates Figures 1, 7, 8 and 9.
%    6. Writes every connectivity variable needed by the downstream
%       hemodynamic and questionnaire analyses to RSFC_outputs.mat.
%
%  REQUIREMENTS
%    MATLAB R2024b
%    Statistics and Machine Learning Toolbox
%    NIRS Brain AnalyzIR Toolbox   https://github.com/huppertt/nirs-toolbox
%
%  INPUT
%    Aim1_raw_loaded_new_03_17_new.mat   preprocessed fNIRS data and
%                                        participant demographics table
%
%  OUTPUT
%    RSFC_outputs.mat   struct RSFC, written to DATA_DIR
%    Fig1.tif Fig7.tif Fig8.tif Fig9.tif   written to FIG_DIR
%
%  HOW TO RUN
%    Keep Aim1_raw_loaded_new_03_17_new.mat in the same folder as this
%    script, then run the script.
%
%  LICENSE   MIT (see LICENSE in the repository root)

clearvars; close all; clc;

%% PATHS

INPUT_FILE = 'Aim1_raw_loaded_new_03_17_new.mat';

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

fprintf('Reading data from : %s\n', DATA_DIR);
fprintf('Figures to        : %s\n', FIG_DIR);
fprintf('Result tables to  : %s\n', RESULTS_DIR);

%% loading data script 
% rootFolder = fullfile(DATA_DIR, 'groupSubject');
% % Map subfolder names 
% folderHierarchy = {'group','subject'};   % matches your tree: Control/Experimental → Sub#
% % If every Sub* contains NIRx folders (not files), default loader is fine:
% raw = nirs.io.loadDirectory(rootFolder, folderHierarchy);
% 
% % Sanity check: you should see group & subject populated
% tbl_demo = nirs.createDemographicsTable(raw);
% 
% disp(tbl_demo(:,{'group','subject'}));
% % % %
% % populate the missing SubjectID
% tbl_demo.SubjectID(35)  = {'Tinnitus01'};
% tbl_demo.SubjectID(46) = {'Tinnitus02'};
% tbl_demo.SubjectID(56) = {'Tinnitus03'};

% 
% % raw(1).draw
% % title("raw data")
% % 
% % raw(1).probe.link % to see if the short channels are labeled (mostly not)
% 
%% save the loaded data for future use 
%save('Aim1_raw_loaded_new_03_17_new.mat', 'raw', 'tbl_demo', '-v7.3'); %Aim1_raw_loaded_new_03_17_new.mat - FINAL

%% load raw data 
load(fullfile(DATA_DIR, 'Aim1_raw_loaded_new_03_17_new.mat'));   % loads raw + tbl

%% fall back 
for i = 1:height(tbl_demo)
    if ~isempty(tbl_demo.SubjectID{i})
        raw(i).demographics('SubjectID') = tbl_demo.SubjectID{i};
    else
        raw(i).demographics('SubjectID') = tbl_demo.subject{i}; % fallback
    end
end
%% apply short channels to all the subjects 
for i = 1:length(raw)
    labelShort = nirs.modules.LabelShortSeperation();
    labelShort.max_distance = 10; % mm
    raw(i) = labelShort.run(raw(i));
end

%% SCI 
% STEP 1: Compute SCI per subject

band = [0.2 2.5];
sciThr = 0.7;

SCI_tables = cell(numel(raw),1);

for i = 1:numel(raw)
    SCI_tables{i} = scalp_coupling_index(raw(i), band, false);

    fprintf('Subject %d (%s): mean SCI = %.3f\n', ...
        i, raw(i).demographics('SubjectID'), mean(SCI_tables{i}.sci));
end


%% REMOVE SUBJECTS WITH MEAN SCI < 0.65

sciCutoff = 0.65;

% compute mean SCI per subject
meanSCI = zeros(numel(raw),1);
for i = 1:numel(raw)
    meanSCI(i) = mean(SCI_tables{i}.sci);
end

% identify good vs bad subjects
keepIdx   = meanSCI >= sciCutoff;
removeIdx = meanSCI <  sciCutoff;

% print removed subjects
fprintf('\n=== REMOVED SUBJECTS (mean SCI < %.2f) ===\n', sciCutoff);
for i = find(removeIdx)'
    fprintf('%s | mean SCI = %.3f\n', ...
        raw(i).demographics('SubjectID'), meanSCI(i));
end
fprintf('Removed %d / %d subjects\n', sum(removeIdx), numel(raw));

% remove from analysis
raw        = raw(keepIdx);
SCI_tables = SCI_tables(keepIdx);
tbl_demo  = tbl_demo(keepIdx, :);

fprintf('Remaining subjects for analysis: %d\n', numel(raw));

%% STEP 2: Optical Density conversion (after SCI)

odJob = nirs.modules.OpticalDensity();
ODdata = odJob.run(raw);   % OD is nirs.core.Data (same size as raw)

% quick sanity check
% ODdata(1).draw
% title('Optical Density (OD)')

%% Step 3 : Motion Correction using TDDR
job = nirs.modules.TDDR();
TDDR_data = job.run(ODdata);

% STEP 3B: Remove channels with NaNs after TDDR

TDDR_clean = TDDR_data;  % copy

for subj = 1:numel(TDDR_clean)

    bad = any(~isfinite(TDDR_clean(subj).data), 1);  % columns with NaN/Inf

    if any(bad)
        link = TDDR_clean(subj).probe.link;

        % Identify bad SD pairs
        bad_pairs = unique([link.source(bad), link.detector(bad)], 'rows');

        remove_idx = false(height(link),1);

        for k = 1:size(bad_pairs,1)
            s = bad_pairs(k,1);
            d = bad_pairs(k,2);
            pair_idx = (link.source == s & link.detector == d);
            remove_idx = remove_idx | pair_idx;
        end

        % Apply removal
        TDDR_clean(subj).data(:,remove_idx) = [];
        TDDR_clean(subj).probe.link(remove_idx,:) = [];

        fprintf('Subject %s: removed %d channels due to TDDR NaNs\n', ...
            TDDR_clean(subj).demographics('SubjectID'), sum(remove_idx));
    end
end


%% Step 4: Short channel regression

ODdata_ssr = TDDR_clean;   % copy

for subj = 1:numel(TDDR_clean)

    link = TDDR_clean(subj).probe.link;

    % identify short and long channels
    isShort = link.ShortSeperation == 1;
    isLong  = link.ShortSeperation == 0;

    shortIdx = find(isShort);
    longIdx  = find(isLong);

    % sanity check
    if isempty(shortIdx)
        error('Subject %s has NO short channels — SSR impossible', ...
            TDDR_clean(subj).demographics('SubjectID'));
    end

    % ---- GLOBAL SHORT REGRESSOR ----
    X = mean(TDDR_clean(subj).data(:, shortIdx), 2);
    X = X - mean(X);

    % ---- REGRESS FROM EACH LONG CHANNEL ----
    for ii = 1:numel(longIdx)
        L = longIdx(ii);

        Y = TDDR_clean(subj).data(:, L);
        Y = Y - mean(Y);

        beta = X \ Y;
        ODdata_ssr(subj).data(:, L) = Y - beta * X;
    end

    fprintf('Subject %s: GLOBAL SSR applied (%d shorts → %d longs)\n', ...
        TDDR_clean(subj).demographics('SubjectID'), ...
        numel(shortIdx), numel(longIdx));
end

%% Step 5: Remove short channels 
ODdata_long = ODdata_ssr;

for i = 1:numel(ODdata_long)
    keep = ~ODdata_long(i).probe.link.ShortSeperation;
    ODdata_long(i).data = ODdata_long(i).data(:,keep);
    ODdata_long(i).probe.link = ODdata_long(i).probe.link(keep,:);
end

%% Step 6 : Band Pass filter
bp = nirs.modules.BandPassFilter();
bp.highpass = 0.01;
bp.lowpass  = 0.08;
bp.do_downsample = false;

ODdata_filt = bp.run(ODdata_long);

%% Step 7 : MBLL


mbll = nirs.modules.BeerLambertLaw();
mbll.PPF = 6;   % standard adult (OK for 760/850 too)

HaemoData = mbll.run(ODdata_filt);
% HaemoData - both HbO & HbR data 

%% plotting montage/ probe
figure('Color','w');
raw(1).probe.draw1020;
sc = findobj(gca,'Type','Scatter');
for k = 1:length(sc)
    if numel(sc(k).XData) > 100   % only the gray grid has hundreds of points
        delete(sc(k));
    end
end
axis equal
axis off

% Save as TIFF at 300 dpi - Fig 1 - PLOS ONE
% output_path = fullfile(FIG_DIR, 'Fig1.tif');
% print(gcf, output_path, '-dtiff', '-r300');

%% Step 8 : Map channels to anatomical ROIs (using Colin27) 
% Register probe of first subject (all subjects share probe geometry)
mesh = nirs.registration.Colin27.mesh();
probe1 = HaemoData(1).probe.register_mesh2probe(mesh);

% Use all channels
link_all = probe1.link;
nCh = height(link_all);

% Compute channel midpoints
xyz = zeros(nCh,3);
for ch = 1:nCh
    s = link_all.source(ch);
    d = link_all.detector(ch);
    xyz(ch,:) = (probe1.srcPos3D(s,:) + probe1.detPos3D(d,:)) / 2;
end

% Map to atlas ROIs
labels   = mesh(end).labels;
keysList = labels.keys;

roiNamesFull  = cell(nCh,1);
roiNamesShort = cell(nCh,1);

for ch = 1:nCh
    minDist = inf;
    closestROI_full  = 'Unknown';
    closestROI_short = 'Unknown';

    for k = 1:length(keysList)
        lbl = labels(keysList{k});
        for j = 1:numel(lbl.VertexIndex)
            verts = lbl.VertexIndex{j};
            nodes = mesh(end).nodes(verts,:);
            d = min(sum((nodes - xyz(ch,:)).^2,2));
            if d < minDist
                minDist = d;
                closestROI_full  = lbl.Region{j};
                closestROI_short = lbl.Label{j};
            end
        end
    end

    roiNamesFull{ch}  = closestROI_full;
    roiNamesShort{ch} = closestROI_short;
end

% Save mapping into a table (include type to make keys unique)
regionTbl = table(link_all.source, link_all.detector, link_all.type, ...
                  roiNamesFull, roiNamesShort, ...
                  xyz(:,1), xyz(:,2), xyz(:,3), ...
    'VariableNames', {'source','detector','type','ROI_code','ROI_full','MNI_X','MNI_Y','MNI_Z'});

fprintf('✅ Region mapping table created for %d channels\n', nCh);

% Attach ROI info back to every subject (unique join)
for i = 1:length(HaemoData)
    HaemoData(i).probe.link = innerjoin(HaemoData(i).probe.link, ...
        regionTbl(:,{'source','detector','type','ROI_code','ROI_full'}), ...
        'Keys',{'source','detector','type'}, ...
        'RightVariables',{'ROI_code','ROI_full'});
end


%% Visualization : Channel-wise plots (HbO + HbR, all channels: long + short)
subj = 1;   % choose subject index
tbl_HaemoData = HaemoData(subj).probe.link;

% logical indices
is_hbo = strcmpi(tbl_HaemoData.type,'hbo');   % all HbO (long + short)
is_hbr = strcmpi(tbl_HaemoData.type,'hbr');   % all HbR (long + short)

hbo_data_all = HaemoData(subj).data(:, is_hbo);
hbr_data_all = HaemoData(subj).data(:, is_hbr);
t = HaemoData(subj).time(:);

% loop through channels
nCh = size(hbo_data_all,2);
for ch = 1
    figure;
    plot(t, hbo_data_all(:,ch), 'r'); hold on;
    plot(t, hbr_data_all(:,ch), 'b');
    legend({'HbO','HbR'});
    
    % label: get source-detector & SS flag from probe.link
    src = tbl_HaemoData.source(is_hbo | is_hbr);
    det = tbl_HaemoData.detector(is_hbo | is_hbr);
    ssFlag = tbl_HaemoData.ShortSeperation(is_hbo | is_hbr);

    if ssFlag(ch)
        ssLabel = ' (Short)';
    else
        ssLabel = ' (Long)';
    end
    
    title(sprintf('Subject %s | Channel %d (S%d-D%d)%s', ...
        HaemoData(subj).demographics('SubjectID'), ch, src(ch), det(ch), ssLabel));
    xlabel('Time (s)'); ylabel('\Delta\muM');
    
    pause; close;
end


%% Split HaemoData into HbO and HbR datasets

HaemoData_hbo = HaemoData;   % copy structure
HaemoData_hbr = HaemoData;   % copy structure

for i = 1:length(HaemoData)
    link = HaemoData(i).probe.link;

    % normalize type column (lowercase, strip spaces)
    types = lower(strtrim(link.type));

    % --- HbO channels ---
    keep_hbo = strcmp(types, 'hbo');
    HaemoData_hbo(i).data = HaemoData(i).data(:, keep_hbo);
    HaemoData_hbo(i).probe.link = link(keep_hbo,:);

    % --- HbR channels ---
    keep_hbr = strcmp(types, 'hbr');   % will now catch properly
    HaemoData_hbr(i).data = HaemoData(i).data(:, keep_hbr);
    HaemoData_hbr(i).probe.link = link(keep_hbr,:);
end

fprintf('✅ Split into HbO (%d channels subj1) and HbR (%d channels subj1)\n', ...
    size(HaemoData_hbo(1).data,2), size(HaemoData_hbr(1).data,2));

%% Step 9 - Run RSFC - Seed-based
% Control Group
% Pre/Post = silence trigger '1'
% Window length = 180 s
% Corr → Fisher Z → pooled seed (n-1 logic)
% GROUP DIFF: per-subject ΔZ then average, then tanh

clc

groupName = 'Control';
durSec    = 180;
silKey    = '1';

roiIdx_hbo_ctrl    = [2 8 21 41];
nonRoiIdx_hbo_ctrl = [60 9 58 18];

fprintf('\n=============================\n');
fprintf('Running CONTROL RSFC (paper-identical)\n');
fprintf('=============================\n');

% ---------- Select CONTROL subjects ----------
group_idx_ctrl  = strcmpi(tbl_demo.group, groupName);
group_data_ctrl = HaemoData_hbo(group_idx_ctrl);

nSub_ctrl = numel(group_data_ctrl);
fprintf('N Control subjects = %d\n', nSub_ctrl);

if nSub_ctrl == 0
    error('No Control subjects found.');
end

% ---------- HbO montage (locked to first subject) ----------
link0_ctrl  = group_data_ctrl(1).probe.link;
isHbO_ctrl  = strcmpi(link0_ctrl.type,'hbo');
nCh_ctrl    = sum(isHbO_ctrl);

assert(all(roiIdx_hbo_ctrl>=1 & roiIdx_hbo_ctrl<=nCh_ctrl));
assert(all(nonRoiIdx_hbo_ctrl>=1 & nonRoiIdx_hbo_ctrl<=nCh_ctrl));

% ---------- Allocate subject-level seed maps ----------
Z_ROI_pre_all_ctrl  = NaN(nSub_ctrl, nCh_ctrl);
Z_ROI_post_all_ctrl = NaN(nSub_ctrl, nCh_ctrl);
Z_NON_pre_all_ctrl  = NaN(nSub_ctrl, nCh_ctrl);
Z_NON_post_all_ctrl = NaN(nSub_ctrl, nCh_ctrl);


% SUBJECT LOOP
for i = 1:nSub_ctrl

    subj = group_data_ctrl(i);
    sid  = subj.demographics('SubjectID');

    if ~ismember(silKey, subj.stimulus.keys)
        warning('%s missing silence trigger', sid);
        continue;
    end

    on = subj.stimulus(silKey).onset;
    if numel(on) < 2
        warning('%s <2 silence blocks', sid);
        continue;
    end

    t = subj.time(:);
    X = subj.data(:, isHbO_ctrl);

    idxPre  = t >= on(1) & t < on(1) + durSec;
    idxPost = t >= on(2) & t < on(2) + durSec;

    if nnz(idxPre) < 10 || nnz(idxPost) < 10
        warning('%s insufficient samples', sid);
        continue;
    end

    % ---------- Channel x Channel connectivity ----------
    Rpre  = corr(X(idxPre,:),  'Rows','pairwise');
    Rpost = corr(X(idxPost,:), 'Rows','pairwise');

    % zero diagonal
    Rpre(1:nCh_ctrl+1:end)  = 0;
    Rpost(1:nCh_ctrl+1:end) = 0;

    % clamp for Fisher
    Rpre  = max(min(Rpre,  0.999999), -0.999999);
    Rpost = max(min(Rpost, 0.999999), -0.999999);

    Zpre  = atanh(Rpre);
    Zpost = atanh(Rpost);

    % ---------- ROI pooled seed (n-1 logic) ----------
    for ch = 1:nCh_ctrl
        roiSet = roiIdx_hbo_ctrl;
        if ismember(ch, roiIdx_hbo_ctrl)
            roiSet = setdiff(roiIdx_hbo_ctrl, ch);
        end
        Z_ROI_pre_all_ctrl(i,ch)  = mean(Zpre(roiSet,ch),  'omitnan');
        Z_ROI_post_all_ctrl(i,ch) = mean(Zpost(roiSet,ch), 'omitnan');
    end

    % ---------- Non-ROI pooled seed (n-1 logic) ----------
    for ch = 1:nCh_ctrl
        nonSet = nonRoiIdx_hbo_ctrl;
        if ismember(ch, nonRoiIdx_hbo_ctrl)
            nonSet = setdiff(nonRoiIdx_hbo_ctrl, ch);
        end
        Z_NON_pre_all_ctrl(i,ch)  = mean(Zpre(nonSet,ch),  'omitnan');
        Z_NON_post_all_ctrl(i,ch) = mean(Zpost(nonSet,ch), 'omitnan');
    end

    fprintf('✓ %s\n', sid);
end


% GROUP MAPS
% pre/post: tanh(meanZ)
% diff: per-subject ΔZ then average then tanh

% --- Mean Z maps ---
meanZ_ROI_pre_ctrl  = nanmean(Z_ROI_pre_all_ctrl,  1);
meanZ_ROI_post_ctrl = nanmean(Z_ROI_post_all_ctrl, 1);

meanZ_NON_pre_ctrl  = nanmean(Z_NON_pre_all_ctrl,  1);
meanZ_NON_post_ctrl = nanmean(Z_NON_post_all_ctrl, 1);

% --- Convert meanZ -> r for pre/post maps ---
R_ROI_pre_ctrl  = tanh(meanZ_ROI_pre_ctrl);
R_ROI_post_ctrl = tanh(meanZ_ROI_post_ctrl);

R_NON_pre_ctrl  = tanh(meanZ_NON_pre_ctrl);
R_NON_post_ctrl = tanh(meanZ_NON_post_ctrl);

% --- DIFF: per-subject ΔZ then average then tanh ---
dZ_ROI_ctrl = Z_ROI_post_all_ctrl - Z_ROI_pre_all_ctrl;
dZ_NON_ctrl = Z_NON_post_all_ctrl - Z_NON_pre_all_ctrl;

R_ROI_diff_ctrl = tanh(nanmean(dZ_ROI_ctrl, 1));
R_NON_diff_ctrl = tanh(nanmean(dZ_NON_ctrl, 1));

fprintf('✓ CONTROL group RSFC maps computed\n');

%% Experimental Group
clc

groupName = 'Experimental';
durSec    = 180;
silKey    = '1';


roiIdx_hbo_exp    = [2 8 21 41];
nonRoiIdx_hbo_exp = [60 9 58 18];

fprintf('\n=============================\n');
fprintf('Running EXPERIMENTAL RSFC (paper-identical)\n');
fprintf('=============================\n');

% ---------- Select EXPERIMENTAL subjects ----------
group_idx_exp  = strcmpi(tbl_demo.group, groupName);
group_data_exp = HaemoData_hbo(group_idx_exp);

nSub_exp = numel(group_data_exp);
fprintf('N Experimental subjects = %d\n', nSub_exp);

if nSub_exp == 0
    error('No Experimental subjects found.');
end

% ---------- HbO montage (locked to first subject) ----------
link0_exp  = group_data_exp(1).probe.link;
isHbO_exp  = strcmpi(link0_exp.type,'hbo');
nCh_exp    = sum(isHbO_exp);

assert(all(roiIdx_hbo_exp>=1 & roiIdx_hbo_exp<=nCh_exp));
assert(all(nonRoiIdx_hbo_exp>=1 & nonRoiIdx_hbo_exp<=nCh_exp));

% ---------- Allocate subject-level seed maps ----------
Z_ROI_pre_all_exp  = NaN(nSub_exp, nCh_exp);
Z_ROI_post_all_exp = NaN(nSub_exp, nCh_exp);
Z_NON_pre_all_exp  = NaN(nSub_exp, nCh_exp);
Z_NON_post_all_exp = NaN(nSub_exp, nCh_exp);


% SUBJECT LOOP
for i = 1:nSub_exp

    subj = group_data_exp(i);
    sid  = subj.demographics('SubjectID');

    if ~ismember(silKey, subj.stimulus.keys)
        warning('%s missing silence trigger', sid);
        continue;
    end

    on = subj.stimulus(silKey).onset;
    if numel(on) < 2
        warning('%s <2 silence blocks', sid);
        continue;
    end

    t = subj.time(:);
    X = subj.data(:, isHbO_exp);

    idxPre  = t >= on(1) & t < on(1) + durSec;
    idxPost = t >= on(2) & t < on(2) + durSec;

    if nnz(idxPre) < 10 || nnz(idxPost) < 10
        warning('%s insufficient samples', sid);
        continue;
    end

    % ---------- Channel x Channel connectivity ----------
    Rpre  = corr(X(idxPre,:),  'Rows','pairwise');
    Rpost = corr(X(idxPost,:), 'Rows','pairwise');

    % zero diagonal
    Rpre(1:nCh_exp+1:end)  = 0;
    Rpost(1:nCh_exp+1:end) = 0;

    % clamp for Fisher
    Rpre  = max(min(Rpre,  0.999999), -0.999999);
    Rpost = max(min(Rpost, 0.999999), -0.999999);

    Zpre  = atanh(Rpre);
    Zpost = atanh(Rpost);

    % ---------- ROI pooled seed (n-1 logic) ----------
    for ch = 1:nCh_exp
        roiSet = roiIdx_hbo_exp;
        if ismember(ch, roiIdx_hbo_exp)
            roiSet = setdiff(roiIdx_hbo_exp, ch);
        end
        Z_ROI_pre_all_exp(i,ch)  = mean(Zpre(roiSet,ch),  'omitnan');
        Z_ROI_post_all_exp(i,ch) = mean(Zpost(roiSet,ch), 'omitnan');
    end

    % ---------- Non-ROI pooled seed (n-1 logic) ----------
    for ch = 1:nCh_exp
        nonSet = nonRoiIdx_hbo_exp;
        if ismember(ch, nonRoiIdx_hbo_exp)
            nonSet = setdiff(nonRoiIdx_hbo_exp, ch);
        end
        Z_NON_pre_all_exp(i,ch)  = mean(Zpre(nonSet,ch),  'omitnan');
        Z_NON_post_all_exp(i,ch) = mean(Zpost(nonSet,ch), 'omitnan');
    end

    fprintf('✓ %s\n', sid);
end


% GROUP MAPS
% pre/post: tanh(meanZ)
% diff: per-subject ΔZ then average then tanh


% --- Mean Z maps ---
meanZ_ROI_pre_exp  = nanmean(Z_ROI_pre_all_exp,  1);
meanZ_ROI_post_exp = nanmean(Z_ROI_post_all_exp, 1);

meanZ_NON_pre_exp  = nanmean(Z_NON_pre_all_exp,  1);
meanZ_NON_post_exp = nanmean(Z_NON_post_all_exp, 1);

% --- Convert meanZ -> r for pre/post maps ---
R_ROI_pre_exp  = tanh(meanZ_ROI_pre_exp);
R_ROI_post_exp = tanh(meanZ_ROI_post_exp);

R_NON_pre_exp  = tanh(meanZ_NON_pre_exp);
R_NON_post_exp = tanh(meanZ_NON_post_exp);

% --- DIFF: per-subject ΔZ then average then tanh ---
dZ_ROI_exp = Z_ROI_post_all_exp - Z_ROI_pre_all_exp;
dZ_NON_exp = Z_NON_post_all_exp - Z_NON_pre_all_exp;

R_ROI_diff_exp = tanh(nanmean(dZ_ROI_exp, 1));
R_NON_diff_exp = tanh(nanmean(dZ_NON_exp, 1));

fprintf('✓ EXPERIMENTAL group RSFC maps computed\n');


%% Fig 7
close all;
% generate_Fig7_PLOS.m


% Panels (3 rows x 2 cols):
%   (a) Control Pre         (b) Control Post
%   (c) Tinnitus Pre        (d) Tinnitus Post
%   (e) Control Difference  (f) Tinnitus Difference

% PREREQUISITE: Run your full RSFC pipeline first so the workspace contains:
%   R_ROI_pre_ctrl,  R_ROI_post_ctrl,  R_ROI_diff_ctrl
%   R_ROI_pre_exp,   R_ROI_post_exp,   R_ROI_diff_exp
%   roiIdx_hbo_ctrl, roiIdx_hbo_exp
%   group_data_ctrl, group_data_exp   (for probe registration)


% ----- Output -----
outDir = FIG_DIR;
if ~exist(outDir, 'dir'); mkdir(outDir); end
outFile = fullfile(outDir, 'Fig7.tif');

% ----- Figure size: full page width, 3 rows x 2 cols -----
figWidthIn  = 7.5;
figHeightIn = 8.5;       % near max (8.75) — needed for brain renderings
dpi         = 300;

% ----- Font settings (Arial 8-12 pt per PLOS) -----
axFont    = 'Arial';
panelSize = 11;          % (a), (b), (c)... bold panel letters
subSize   = 10;          % "Pre", "Post", "Difference" sub-labels
cbSize    = 9;           % colorbar tick + label
chLabSize = 7;           % channel number labels (small to fit)

% ----- Visual parameters -----
doLabels     = true;     % show channel number labels on brain
rSphere      = 4;
clim_prepost = [-1 1];
clim_diff    = [-0.5 0.5];
cmap         = jet(256);

% ----- Panel definitions (data, group/condition, panel letter, sub-label) -----
panels = {
    R_ROI_pre_ctrl,   roiIdx_hbo_ctrl, group_data_ctrl, '(a)', 'Control Pre',         clim_prepost;
    R_ROI_post_ctrl,  roiIdx_hbo_ctrl, group_data_ctrl, '(b)', 'Control Post',        clim_prepost;
    R_ROI_pre_exp,    roiIdx_hbo_exp,  group_data_exp,  '(c)', 'Tinnitus Pre',        clim_prepost;
    R_ROI_post_exp,   roiIdx_hbo_exp,  group_data_exp,  '(d)', 'Tinnitus Post',       clim_prepost;
    R_ROI_diff_ctrl,  roiIdx_hbo_ctrl, group_data_ctrl, '(e)', 'Control Difference',  clim_diff;
    R_ROI_diff_exp,   roiIdx_hbo_exp,  group_data_exp,  '(f)', 'Tinnitus Difference', clim_diff;
};


% Build figure
fig = figure('Units','inches', ...
             'Position',[1 1 figWidthIn figHeightIn], ...
             'Color','w', ...
             'PaperUnits','inches', ...
             'PaperPosition',[0 0 figWidthIn figHeightIn], ...
             'PaperSize',[figWidthIn figHeightIn]);

% ----- Manual axis positions (3 rows x 2 cols) -----
leftMargin   = 0.04;
rightMargin  = 0.12;     % wider to fit colorbars
bottomMargin = 0.03;
topMargin    = 0.03;
hGap         = 0.02;
vGap         = 0.04;

panelW = (1 - leftMargin - rightMargin - hGap) / 2;
panelH = (1 - bottomMargin - topMargin - 2*vGap) / 3;

positions = cell(6,1);
for r = 0:2
    for c = 0:1
        idx = 2*r + c + 1;
        left   = leftMargin + c*(panelW + hGap);
        bottom = 1 - topMargin - panelH - r*(panelH + vGap);
        positions{idx} = [left bottom panelW panelH];
    end
end

% ----- Pre-load mesh once (shared across panels) -----
mesh = nirs.registration.Colin27.mesh();
for ii = 1:length(mesh)
    if istable(mesh(ii).fiducials)
        mesh(ii).fiducials.Draw(:) = false;
    end
end

% Loop through panels
for d = 1:size(panels,1)
    vals       = panels{d,1};
    roiIdx     = panels{d,2};
    group_data = panels{d,3};
    panelLab   = panels{d,4};
    subLab     = panels{d,5};
    clim       = panels{d,6};

    % --- Register probe to mesh (subject-specific geometry) ---
    probe = group_data(1).probe.register_mesh2probe(mesh);
    probe = probe.SetFiducials_Visibility(false);

    linkP      = probe.link;
    isHbO_plot = strcmpi(linkP.type,'hbo');
    link_hbo   = linkP(isHbO_plot,:);

    srcPos3D = probe.srcPos3D;
    detPos3D = probe.detPos3D;

    nCh_plot = height(link_hbo);
    xyz = zeros(nCh_plot, 3);
    for ch = 1:nCh_plot
        xyz(ch,:) = (srcPos3D(link_hbo.source(ch),:) + ...
                     detPos3D(link_hbo.detector(ch),:)) / 2;
    end

    % --- Zero out the ROI seed channels (they're the seed, no self-corr) ---
    valsPlot = vals;
    valsPlot(roiIdx) = 0;

    % --- Create axes at manual position ---
    ax = axes('Parent', fig, 'Position', positions{d});

    % --- Draw mesh into this axes ---
    hP = mesh.draw([],[],[],[],ax);
    delete(hP(1:4));
    set(hP(5),'FaceAlpha',0.6,'EdgeColor','none');
    hold(ax,'on');

    % --- Draw channel spheres ---
    [sx, sy, sz] = sphere(20);
    for ch = 1:nCh_plot
        val = valsPlot(ch);
        if ~isfinite(val), continue; end

        idx = round(((val - clim(1)) / (clim(2) - clim(1))) ...
                    * (size(cmap,1)-1)) + 1;
        idx = min(max(idx,1), size(cmap,1));

        surf(ax, rSphere*sx + xyz(ch,1), ...
                 rSphere*sy + xyz(ch,2), ...
                 rSphere*sz + xyz(ch,3), ...
                 'FaceColor', cmap(idx,:), ...
                 'EdgeColor', 'none');

        if doLabels
            v      = xyz(ch,:);
            cam    = campos(ax);
            dirVec = cam - v;
            dirVec = dirVec / norm(dirVec);
            labelPos = v + (rSphere * 0.6) * dirVec;
            text(ax, labelPos(1), labelPos(2), labelPos(3), ...
                sprintf('%d', ch), ...
                'FontName',axFont, 'FontSize', chLabSize, ...
                'FontWeight', 'bold', 'Color', 'k', ...
                'HorizontalAlignment', 'center', ...
                'VerticalAlignment', 'middle', ...
                'Clipping', 'off');
        end
    end

    axis(ax, 'equal', 'off');
    set(ax, 'SortMethod','childorder');
    camlight(ax, 'headlight');
    lighting(ax, 'gouraud');

    colormap(ax, cmap);
    caxis(ax, clim);

    % % --- Sub-title above panel ---
    % title(ax, subLab, ...
    %     'FontName',axFont, 'FontSize',subSize, 'FontWeight','normal');

    % --- Panel label (a), (b)... outside top-left ---
    text(ax, -0.02, 1.05, panelLab, ...
        'Units','normalized', ...
        'FontName',axFont, 'FontSize',panelSize, 'FontWeight','bold');
end

% Two shared colorbars on the right side
%   - Top colorbar: Pre/Post range [-1, 1]   (rows 1 & 2)
%   - Bottom colorbar: Diff range  [-0.5, 0.5] (row 3)

% --- Pre/Post colorbar (covers rows 1 and 2) ---
cbAx_pp = axes('Parent', fig, ...
    'Position', [1-rightMargin+0.015, ...
                 positions{4}(2), ...               % bottom of row 2
                 0.025, ...
                 positions{2}(2)+positions{2}(4) - positions{4}(2)]);  % spans rows 1-2
colormap(cbAx_pp, cmap);
caxis(cbAx_pp, clim_prepost);
cb_pp = colorbar(cbAx_pp, 'Position', cbAx_pp.Position);
ylabel(cb_pp, 'Correlation (r)', 'FontName',axFont, 'FontSize',cbSize);
set(cb_pp, 'FontName',axFont, 'FontSize',cbSize-1);
axis(cbAx_pp, 'off');

% --- Diff colorbar (row 3 only) ---
cbAx_d = axes('Parent', fig, ...
    'Position', [1-rightMargin+0.015, ...
                 positions{6}(2), ...
                 0.025, ...
                 positions{6}(4)]);
colormap(cbAx_d, cmap);
caxis(cbAx_d, clim_diff);
cb_d = colorbar(cbAx_d, 'Position', cbAx_d.Position);
ylabel(cb_d, '\Delta Correlation', 'FontName',axFont, 'FontSize',cbSize);
set(cb_d, 'FontName',axFont, 'FontSize',cbSize-1);
axis(cbAx_d, 'off');


% Export as PLOS-compliant TIFF
% delete(findall(fig,'Type','hggroup'));
% delete(findall(fig,'Tag','DataTipMarker'));
% 
% try
%     exportgraphics(fig, outFile, 'Resolution', dpi, 'BackgroundColor','white');
%     fprintf('Saved (exportgraphics): %s\n', outFile);
% catch
%     print(fig, outFile, '-dtiff', sprintf('-r%d', dpi));
%     fprintf('Saved (print): %s\n', outFile);
% end

% Compliance check
% info = imfinfo(outFile);
% fprintf('\n--- Fig7.tif specs ---\n');
% fprintf('  Width:       %d px\n', info.Width);
% fprintf('  Height:      %d px\n', info.Height);
% fprintf('  XRes:        %g dpi\n', info.XResolution);
% fprintf('  YRes:        %g dpi\n', info.YResolution);
% fprintf('  ColorType:   %s, BitDepth: %d\n', info.ColorType, info.BitDepth);
% fprintf('  Compression: %s\n', info.Compression);
% fprintf('  File size:   %.2f MB\n', info.FileSize / 1024 / 1024);
% 
% fprintf('\n--- PLOS ONE compliance check ---\n');
% ok = true;
% if info.Width < 789 || info.Width > 2250
%     fprintf('  [WARN] Width %d outside 789-2250 px range\n', info.Width); ok = false;
% end
% if info.Height > 2625
%     fprintf('  [WARN] Height %d exceeds 2625 px max\n', info.Height); ok = false;
% end
% if info.XResolution < 300 || info.XResolution > 600
%     fprintf('  [WARN] DPI %g outside 300-600 range\n', info.XResolution); ok = false;
% end
% if (info.FileSize / 1024 / 1024) > 10
%     fprintf('  [WARN] File size exceeds 10 MB\n'); ok = false;
% end
% if ok
%     fprintf('  PASS: dimensions, resolution, and file size all within PLOS specs.\n');
% end
% fprintf('  NOTE: If compression is not LZW, re-save in GIMP/Photoshop with LZW.\n');

%% Step 10 : ROI-SEED RSFC STATS 
%   Fisher Z for each subject, then ΔZ = Zpost - Zpre,
%   stats in Z-space, and (if desired) report Δr = tanh(mean(ΔZ)).

excludeSeedTargets = true;

% ---- ROI indices (use group-specific, safer) ----
roiIdxC = roiIdx_hbo_ctrl;
roiIdxE = roiIdx_hbo_exp;  

% ---- Columns to pool across (optionally exclude seed targets) ----
if excludeSeedTargets
    keepColsC = setdiff(1:size(Z_ROI_pre_all_ctrl,2), roiIdxC);
    keepColsE = setdiff(1:size(Z_ROI_pre_all_exp,2),  roiIdxE);
else
    keepColsC = 1:size(Z_ROI_pre_all_ctrl,2);
    keepColsE = 1:size(Z_ROI_pre_all_exp,2);
end

% ---- Subject-level pooled ROI connectivity (in Fisher Z) ----
zPre_ctrl_all  = nanmean(Z_ROI_pre_all_ctrl(:,  keepColsC), 2);
zPost_ctrl_all = nanmean(Z_ROI_post_all_ctrl(:, keepColsC), 2);

zPre_exp_all   = nanmean(Z_ROI_pre_all_exp(:,  keepColsE), 2);
zPost_exp_all  = nanmean(Z_ROI_post_all_exp(:, keepColsE), 2);

% ---- Keep only subjects with BOTH pre and post (paired) ----
maskC = isfinite(zPre_ctrl_all) & isfinite(zPost_ctrl_all);
zPre_ctrl  = zPre_ctrl_all(maskC);
zPost_ctrl = zPost_ctrl_all(maskC);
zDiff_ctrl = zPost_ctrl - zPre_ctrl;      % ΔZ 

maskE = isfinite(zPre_exp_all) & isfinite(zPost_exp_all);
zPre_exp  = zPre_exp_all(maskE);
zPost_exp = zPost_exp_all(maskE);
zDiff_exp = zPost_exp - zPre_exp;         % ΔZ 

fprintf('\n================ ROI-SEED RSFC STATS ================\n');

% (1) One-sample t-tests vs zero (per group, per window) IN Z-SPACE
fprintf('\n--- ONE-SAMPLE t-tests in Fisher-Z (H0: mean(Z)=0) ---\n');

[~, pC_pre,  ~, sC_pre]  = ttest(zPre_ctrl,  0);
[~, pC_post, ~, sC_post] = ttest(zPost_ctrl, 0);

[~, pE_pre,  ~, sE_pre]  = ttest(zPre_exp,   0);
[~, pE_post, ~, sE_post] = ttest(zPost_exp,  0);

fprintf('\nCONTROL:\n');
fprintf('  PRE   n=%d | meanZ=%.4f | t=%.3f | p=%.4g | reported r=%.4f\n', ...
    numel(zPre_ctrl),  mean(zPre_ctrl,'omitnan'),  sC_pre.tstat,  pC_pre,  tanh(mean(zPre_ctrl,'omitnan')));
fprintf('  POST  n=%d | meanZ=%.4f | t=%.3f | p=%.4g | reported r=%.4f\n', ...
    numel(zPost_ctrl), mean(zPost_ctrl,'omitnan'), sC_post.tstat, pC_post, tanh(mean(zPost_ctrl,'omitnan')));

fprintf('\nEXPERIMENTAL (Tinnitus):\n');
fprintf('  PRE   n=%d | meanZ=%.4f | t=%.3f | p=%.4g | reported r=%.4f\n', ...
    numel(zPre_exp),   mean(zPre_exp,'omitnan'),   sE_pre.tstat,  pE_pre,  tanh(mean(zPre_exp,'omitnan')));
fprintf('  POST  n=%d | meanZ=%.4f | t=%.3f | p=%.4g | reported r=%.4f\n', ...
    numel(zPost_exp),  mean(zPost_exp,'omitnan'),  sE_post.tstat, pE_post, tanh(mean(zPost_exp,'omitnan')));

% (2) Two-sample group comparisons (Control vs Tinnitus) at PRE and POST (Z-SPACE)
fprintf('\n--- TWO-SAMPLE t-tests (Control vs Tinnitus) in Fisher-Z ---\n');

[~, p_pre,  ~, s_pre]  = ttest2(zPre_ctrl,  zPre_exp);
[~, p_post, ~, s_post] = ttest2(zPost_ctrl, zPost_exp);

fprintf('\nPRE (Control vs Tinnitus):\n');
fprintf('  Control meanZ=%.4f | Tinnitus meanZ=%.4f | t=%.3f | p=%.4g\n', ...
    mean(zPre_ctrl,'omitnan'), mean(zPre_exp,'omitnan'), s_pre.tstat, p_pre);

fprintf('\nPOST (Control vs Tinnitus):\n');
fprintf('  Control meanZ=%.4f | Tinnitus meanZ=%.4f | t=%.3f | p=%.4g\n', ...
    mean(zPost_ctrl,'omitnan'), mean(zPost_exp,'omitnan'), s_post.tstat, p_post);

% (3) Within-group Post–Pre change: one-sample t-test on ΔZ vs 0
fprintf('\n--- WITHIN-GROUP POST−PRE (one-sample t-test of ΔZ vs 0) ---\n');

[~, p_dC, ~, s_dC] = ttest(zDiff_ctrl, 0);
[~, p_dE, ~, s_dE] = ttest(zDiff_exp,  0);

fprintf('\nCONTROL ΔZ:\n');
fprintf('  n=%d | meanΔZ=%.4f | t=%.3f | p=%.4g | paper Δr=tanh(meanΔZ)=%.4f\n', ...
    numel(zDiff_ctrl), mean(zDiff_ctrl,'omitnan'), s_dC.tstat, p_dC, tanh(mean(zDiff_ctrl,'omitnan')));

fprintf('\nTINNITUS ΔZ:\n');
fprintf('  n=%d | meanΔZ=%.4f | t=%.3f | p=%.4g | paper Δr=tanh(meanΔZ)=%.4f\n', ...
    numel(zDiff_exp), mean(zDiff_exp,'omitnan'), s_dE.tstat, p_dE, tanh(mean(zDiff_exp,'omitnan')));

% (4) Difference-of-differences: compare ΔZ between groups 
fprintf('\n--- DIFFERENCE-OF-DIFFERENCES: ΔZ(Control) vs ΔZ(Tinnitus) ---\n');

[~, p_dd, ~, s_dd] = ttest2(zDiff_ctrl, zDiff_exp);
fprintf('  t=%.3f | p=%.4g\n', s_dd.tstat, p_dd);

%% Step 11 : ROI → LOBE-WISE RSFC

% Uses existing:
%   Z_ROI_pre_all_ctrl, Z_ROI_post_all_ctrl
%   Z_ROI_pre_all_exp,  Z_ROI_post_all_exp
%   roiIdx_hbo_ctrl (and optionally roiIdx_hbo_exp)

%  Lobe dictionary

lobeDict = containers.Map;

% Frontal (F)
F = {'Precentral_L','Precentral_R', ...
     'Frontal_Sup_L','Frontal_Sup_R', ...
     'Frontal_Mid_L','Frontal_Mid_R'};
for k = 1:numel(F), lobeDict(F{k}) = 'F'; end

% Fronto-Temporal (FT)
FT = {'Frontal_Inf_Oper_L', ...
      'Frontal_Inf_Tri_R', ...
      'Rolandic_Oper_L'};
for k = 1:numel(FT), lobeDict(FT{k}) = 'FT'; end

% FP
FP = {'Postcentral_L','Postcentral_R','Paracentral_Lobule_L','Paracentral_Lobule_R'};
for k = 1:numel(FP), lobeDict(FP{k}) = 'FP'; end

% TP
TP = {'Angular_L','Angular_R','SupraMarginal_L','SupraMarginal_R'};
for k = 1:numel(TP), lobeDict(TP{k}) = 'TP'; end

% T
T = {'Temporal_Sup_L','Temporal_Sup_R','Temporal_Pole_Sup_L','Temporal_Pole_Sup_R', ...
     'Temporal_Mid_L','Temporal_Mid_R','Temporal_Pole_Mid_L','Temporal_Pole_Mid_R'};
for k = 1:numel(T), lobeDict(T{k}) = 'T'; end

% P
P = {'Precuneus_L','Precuneus_R','Parietal_Sup_L','Parietal_Sup_R','Parietal_Inf_L','Parietal_Inf_R'};
for k = 1:numel(P), lobeDict(P{k}) = 'P'; end

% OP
OP = {'Cuneus_R','Calcarine_L'};
for k = 1:numel(OP), lobeDict(OP{k}) = 'OP'; end

% OT
OT = {'Fusiform_L','Fusiform_R','Temporal_Inf_L','Temporal_Inf_R','Lingual_R'};
for k = 1:numel(OT), lobeDict(OT{k}) = 'OT'; end

% O
O = {'Occipital_Sup_L','Occipital_Sup_R','Occipital_Mid_L','Occipital_Mid_R', ...
     'Occipital_Inf_L','Occipital_Inf_R'};
for k = 1:numel(O), lobeDict(O{k}) = 'O'; end


% 1) Per-HbO-channel ROI labels

roiLabel_hbo = HaemoData_hbo(1).probe.link.ROI_full;

if isstring(roiLabel_hbo), roiLabel_hbo = cellstr(roiLabel_hbo); end
if iscell(roiLabel_hbo) && size(roiLabel_hbo,2) == 1
    roiLabel_hbo = roiLabel_hbo(:)';   % force 1 x nCh
end

nCh = size(Z_ROI_pre_all_ctrl,2);
assert(numel(roiLabel_hbo)==nCh, ...
    'Mismatch: roiLabel_hbo length (%d) vs HbO channels in Z (%d).', numel(roiLabel_hbo), nCh);


% 2) Map ROI label → lobe code

chanLobe_hbo = repmat({'UNK'}, 1, nCh);

for ch = 1:nCh
    lab = roiLabel_hbo{ch};
    if isstring(lab), lab = char(lab); end
    if isempty(lab) || ~ischar(lab)
        chanLobe_hbo{ch} = 'UNK';
        continue;
    end

    if isKey(lobeDict, lab)
        chanLobe_hbo{ch} = lobeDict(lab);
    else
        chanLobe_hbo{ch} = 'UNK';
    end
end


% 3) Coverage report

fprintf('\n=== LOBE MAPPING COVERAGE (HbO channels) ===\n');
lobesAll = {'F','FT','FP','TP','T','P','OP','OT','O','UNK'};
for k = 1:numel(lobesAll)
    fprintf('%3s : %d\n', lobesAll{k}, sum(strcmpi(chanLobe_hbo, lobesAll{k})));
end


% 4) Exclude ROI target channels from pooling (recommended)

excludeSeedTargets = true;

% If your ROI indices are identical for both groups, keep this:
roiIdx = roiIdx_hbo_ctrl;


validTarget = true(1,nCh);
if excludeSeedTargets
    validTarget(roiIdx) = false;
end


% 5) Compute ROI→lobe subject means (Z-space)

lobes = {'F','FT','FP','TP','T','P','OP','OT','O'};
nL = numel(lobes);

Zlob_pre_ctrl  = computeLobeMeans(Z_ROI_pre_all_ctrl,  chanLobe_hbo, lobes, validTarget);
Zlob_post_ctrl = computeLobeMeans(Z_ROI_post_all_ctrl, chanLobe_hbo, lobes, validTarget);

Zlob_pre_exp   = computeLobeMeans(Z_ROI_pre_all_exp,   chanLobe_hbo, lobes, validTarget);
Zlob_post_exp  = computeLobeMeans(Z_ROI_post_all_exp,  chanLobe_hbo, lobes, validTarget);

fprintf('\n✓ ROI → LOBE matrices computed (Z-space)\n');


% 6) Stats per lobe (Fisher-Z, SUBJECT is unit)


alpha = 0.05;
fprintf('\n================ ROI → LOBE-WISE PAPER-STYLE STATS ================\n');
fprintf('Lobes order: %s\n', strjoin(lobes, ', '));

for L = 1:nL
    name = lobes{L};

    % ---------- CONTROL paired subjects ----------
    zPreC_all  = Zlob_pre_ctrl(:,L);
    zPostC_all = Zlob_post_ctrl(:,L);
    mC = isfinite(zPreC_all) & isfinite(zPostC_all);

    zPreC  = zPreC_all(mC);
    zPostC = zPostC_all(mC);
    zDiffC = zPostC - zPreC;        % ΔZ per subject 

    % ---------- EXPERIMENTAL paired subjects ----------
    zPreE_all  = Zlob_pre_exp(:,L);
    zPostE_all = Zlob_post_exp(:,L);
    mE = isfinite(zPreE_all) & isfinite(zPostE_all);

    zPreE  = zPreE_all(mE);
    zPostE = zPostE_all(mE);
    zDiffE = zPostE - zPreE;        % ΔZ per subject 

    % one-sample vs 0 (pre/post)
    [~,pCpre,~,sCpre]   = ttest(zPreC,  0);
    [~,pCpost,~,sCpost] = ttest(zPostC, 0);
    [~,pEpre,~,sEpre]   = ttest(zPreE,  0);
    [~,pEpost,~,sEpost] = ttest(zPostE, 0);

    % between-group (pre/post)
    [~,pPre,~,sPre]   = ttest2(zPreC,  zPreE);
    [~,pPost,~,sPost] = ttest2(zPostC, zPostE);

    % within-group ΔZ vs 0
    [~,pDC,~,sDC] = ttest(zDiffC, 0);
    [~,pDE,~,sDE] = ttest(zDiffE, 0);

    % diff-of-diff on ΔZ
    [~,pDD,~,sDD] = ttest2(zDiffC, zDiffE);

    fprintf('\n[%s]\n', name);

    fprintf('  One-sample vs 0:\n');
    fprintf('    CTRL PRE   n=%2d meanZ=% .3f t=% .2f p=% .4g r=% .3f\n', ...
        numel(zPreC),  mean(zPreC,'omitnan'),  sCpre.tstat,  pCpre,  tanh(mean(zPreC,'omitnan')));
    fprintf('    CTRL POST  n=%2d meanZ=% .3f t=% .2f p=% .4g r=% .3f\n', ...
        numel(zPostC), mean(zPostC,'omitnan'), sCpost.tstat, pCpost, tanh(mean(zPostC,'omitnan')));
    fprintf('    TINN PRE   n=%2d meanZ=% .3f t=% .2f p=% .4g r=% .3f\n', ...
        numel(zPreE),  mean(zPreE,'omitnan'),  sEpre.tstat,  pEpre,  tanh(mean(zPreE,'omitnan')));
    fprintf('    TINN POST  n=%2d meanZ=% .3f t=% .2f p=% .4g r=% .3f\n', ...
        numel(zPostE), mean(zPostE,'omitnan'), sEpost.tstat, pEpost, tanh(mean(zPostE,'omitnan')));

    fprintf('  Between-group (CTRL vs TINN):\n');
    fprintf('    PRE  t=% .2f p=% .4g\n',  sPre.tstat,  pPre);
    fprintf('    POST t=% .2f p=% .4g\n', sPost.tstat, pPost);

    fprintf('  Within-group ΔZ (POST−PRE):\n');
    fprintf('    CTRL  n=%2d meanΔZ=% .3f t=% .2f p=% .4g | paper Δr=tanh(meanΔZ)=% .3f\n', ...
        numel(zDiffC), mean(zDiffC,'omitnan'), sDC.tstat, pDC, tanh(mean(zDiffC,'omitnan')));
    fprintf('    TINN  n=%2d meanΔZ=% .3f t=% .2f p=% .4g | paper Δr=tanh(meanΔZ)=% .3f\n', ...
        numel(zDiffE), mean(zDiffE,'omitnan'), sDE.tstat, pDE, tanh(mean(zDiffE,'omitnan')));

    fprintf('  Diff-of-diff on ΔZ:\n');
    fprintf('    t=% .2f p=% .4g\n', sDD.tstat, pDD);

    if pDD < alpha
        fprintf('    ** Interaction significant at alpha=%.2f **\n', alpha);
    end
end

fprintf('\n===============================================================\n');


% 7)  lobe diffs as r (for plotting)
%   Δr_lobe = tanh(mean(ΔZ_lobe))  where ΔZ_lobe is subject-level

Zlob_diff_ctrl = Zlob_post_ctrl - Zlob_pre_ctrl;
Zlob_diff_exp  = Zlob_post_exp  - Zlob_pre_exp;

Rlob_pre_ctrl   = tanh(nanmean(Zlob_pre_ctrl,  1));
Rlob_post_ctrl  = tanh(nanmean(Zlob_post_ctrl, 1));
Rlob_diff_ctrl  = tanh(nanmean(Zlob_diff_ctrl, 1));   % paper-style diff

Rlob_pre_exp    = tanh(nanmean(Zlob_pre_exp,   1));
Rlob_post_exp   = tanh(nanmean(Zlob_post_exp,  1));
Rlob_diff_exp   = tanh(nanmean(Zlob_diff_exp,  1));   % paper-style diff

% Now use Rlob_* for any r-space plots.



%% LOBE WISE SUMMARY 
fprintf('\n=== BASELINE (PRE) meanZ BY LOBE ===\n');
fprintf('Lobe | CTRL meanZ (r)        | TINN meanZ (r)\n');
fprintf('-----|------------------------|------------------------\n');

for k = 1:numel(lobes)
    mC = mean(Zlob_pre_ctrl(:,k), 'omitnan');
    mE = mean(Zlob_pre_exp(:,k),  'omitnan');

    fprintf('%-3s  | %8.4f (%6.3f)     | %8.4f (%6.3f)\n', ...
        lobes{k}, mC, tanh(mC), mE, tanh(mE));
end

fprintf('\n=== POST (POST-stimulation) meanZ BY LOBE ===\n');
fprintf('Lobe | CTRL meanZ (r)        | TINN meanZ (r)\n');
fprintf('-----|------------------------|------------------------\n');

for k = 1:numel(lobes)
    mC = mean(Zlob_post_ctrl(:,k), 'omitnan');
    mE = mean(Zlob_post_exp(:,k),  'omitnan');

    fprintf('%-3s  | %8.4f (%6.3f)     | %8.4f (%6.3f)\n', ...
        lobes{k}, mC, tanh(mC), mE, tanh(mE));
end


%% BONFERRONI CORRECTION — LOBE-WISE RSFC
% Applied to within-group ΔZ (primary findings)
% and between-group comparisons


nLobes     = numel(lobes);       % 9
alpha_nom  = 0.05;
alpha_bonf = alpha_nom / nLobes; % 0.05/9 = 0.0056

fprintf('\n================ BONFERRONI CORRECTION (alpha/9 = %.4f) ================\n', alpha_bonf);
fprintf('%-4s | %-20s | %-10s | %-10s | %-8s || %-20s | %-10s | %-10s | %-8s\n', ...
    'Lobe', 'TINN deltaZ', 'p_unc', 'p_bonf', 'sig*', ...
    'CTRL deltaZ', 'p_unc', 'p_bonf', 'sig*');
fprintf('%s\n', repmat('-',1,105));

% Storage for summary table
bonf_results = table();

for L = 1:nLobes
    name = lobes{L};

    % --- Tinnitus ΔZ ---
    zPreE_all  = Zlob_pre_exp(:,L);
    zPostE_all = Zlob_post_exp(:,L);
    mE = isfinite(zPreE_all) & isfinite(zPostE_all);
    zDiffE = zPostE_all(mE) - zPreE_all(mE);

    % --- Control ΔZ ---
    zPreC_all  = Zlob_pre_ctrl(:,L);
    zPostC_all = Zlob_post_ctrl(:,L);
    mC = isfinite(zPreC_all) & isfinite(zPostC_all);
    zDiffC = zPostC_all(mC) - zPreC_all(mC);

    % --- Between-group ΔZ ---
    [~, pBG, ~, sBG] = ttest2(zDiffC, zDiffE);

    % --- Within-group tests ---
    [~, pE, ~, sE] = ttest(zDiffE, 0);
    [~, pC, ~, sC] = ttest(zDiffC, 0);

    pE_bonf = min(pE * nLobes, 1);
    pC_bonf = min(pC * nLobes, 1);
    pBG_bonf = min(pBG * nLobes, 1);

    sigE  = pE_bonf  < alpha_nom;
    sigC  = pC_bonf  < alpha_nom;
    sigBG = pBG_bonf < alpha_nom;

    fprintf('%-4s | meanΔZ=%6.3f t=%5.2f | p=%7.4f | p_b=%7.4f | %-5s || meanΔZ=%6.3f t=%5.2f | p=%7.4f | p_b=%7.4f | %-5s\n', ...
        name, ...
        mean(zDiffE,'omitnan'), sE.tstat, pE, pE_bonf, string(sigE), ...
        mean(zDiffC,'omitnan'), sC.tstat, pC, pC_bonf, string(sigC));

    newRow = table({name}, ...
        mean(zDiffE,'omitnan'), sE.tstat, pE, pE_bonf, sigE, ...
        mean(zDiffC,'omitnan'), sC.tstat, pC, pC_bonf, sigC, ...
        sBG.tstat, pBG, pBG_bonf, sigBG, ...
        'VariableNames', {'Lobe', ...
        'TINN_meanDZ','TINN_t','TINN_p_unc','TINN_p_bonf','TINN_sig', ...
        'CTRL_meanDZ','CTRL_t','CTRL_p_unc','CTRL_p_bonf','CTRL_sig', ...
        'BG_t','BG_p_unc','BG_p_bonf','BG_sig'});

    bonf_results = [bonf_results; newRow];
end

fprintf('\n--- SUMMARY: Tinnitus ΔZ surviving Bonferroni ---\n');
disp(bonf_results(bonf_results.TINN_sig == true, {'Lobe','TINN_meanDZ','TINN_t','TINN_p_unc','TINN_p_bonf'}));

fprintf('\n--- SUMMARY: Between-group ΔZ surviving Bonferroni ---\n');
disp(bonf_results(bonf_results.BG_sig == true, {'Lobe','BG_t','BG_p_unc','BG_p_bonf'}));

assignin('base', 'bonf_results_rsfc', bonf_results);
fprintf('\n✓ Bonferroni correction done. Results saved to bonf_results_rsfc.\n');


%% Step 12 : WHOLE BRAIN CONNECTIVITY
durSec = 180;
silKey = '1';

groupsToRun = {'Control','Experimental'};

% ---- 9-lobe dictionary ----
lobes = {'F','FT','FP','TP','T','P','OP','OT','O'};
nL = numel(lobes);

lobeDict = containers.Map;

F = {'Precentral_L','Precentral_R', ...
     'Frontal_Sup_L','Frontal_Sup_R', ...
     'Frontal_Mid_L','Frontal_Mid_R'};
for k = 1:numel(F), lobeDict(F{k}) = 'F'; end

FT = {'Frontal_Inf_Oper_L','Frontal_Inf_Tri_R','Rolandic_Oper_L'};
for k = 1:numel(FT), lobeDict(FT{k}) = 'FT'; end

FP = {'Postcentral_L','Postcentral_R','Paracentral_Lobule_L','Paracentral_Lobule_R'};
for k = 1:numel(FP), lobeDict(FP{k}) = 'FP'; end

TP = {'Angular_L','Angular_R','SupraMarginal_L','SupraMarginal_R'};
for k = 1:numel(TP), lobeDict(TP{k}) = 'TP'; end

T = {'Temporal_Sup_L','Temporal_Sup_R','Temporal_Pole_Sup_L','Temporal_Pole_Sup_R', ...
     'Temporal_Mid_L','Temporal_Mid_R','Temporal_Pole_Mid_L','Temporal_Pole_Mid_R'};
for k = 1:numel(T), lobeDict(T{k}) = 'T'; end

P = {'Precuneus_L','Precuneus_R','Parietal_Sup_L','Parietal_Sup_R', ...
     'Parietal_Inf_L','Parietal_Inf_R'};
for k = 1:numel(P), lobeDict(P{k}) = 'P'; end

OP = {'Cuneus_L','Cuneus_R','Calcarine_L'};
for k = 1:numel(OP), lobeDict(OP{k}) = 'OP'; end

OT = {'Fusiform_L','Fusiform_R','Temporal_Inf_L','Temporal_Inf_R','Lingual_R'};
for k = 1:numel(OT), lobeDict(OT{k}) = 'OT'; end

O = {'Occipital_Sup_L','Occipital_Sup_R','Occipital_Mid_L','Occipital_Mid_R', ...
     'Occipital_Inf_L','Occipital_Inf_R'};
for k = 1:numel(O), lobeDict(O{k}) = 'O'; end


% =========================================================
% RUN PER GROUP
% =========================================================
dR_stack_control = [];
dR_stack_exp     = [];
chanLobe_global  = [];

for gg = 1:numel(groupsToRun)

    groupName = groupsToRun{gg};
    fprintf('\n==================== %s ====================\n', upper(groupName));

    group_idx  = strcmpi(tbl_demo.group, groupName);
    group_data = HaemoData_hbo(group_idx);

    nSub = numel(group_data);
    fprintf('N subjects = %d\n', nSub);
    if nSub == 0
        warning('No subjects found for %s', groupName);
        continue;
    end

    % ---- HbO channel selection ----
    link0 = group_data(1).probe.link;
    isHbO = strcmpi(link0.type,'hbo');
    nCh   = sum(isHbO);

    % ---- ROI_full labels ----
    if ~ismember('ROI_full', link0.Properties.VariableNames)
        error('probe.link does not contain ROI_full.');
    end
    roi_full_hbo = link0.ROI_full(isHbO);

    % ---- Map channels to lobes ----
    chanLobe = repmat("UNK", nCh, 1);
    for ch = 1:nCh
        lab = roi_full_hbo{ch};
        if isstring(lab), lab = char(lab); end
        if ischar(lab) && isKey(lobeDict, lab)
            chanLobe(ch) = string(lobeDict(lab));
        end
    end

    if isempty(chanLobe_global)
        chanLobe_global = chanLobe;
    end

    % ---- Coverage report ----
    fprintf('\n=== LOBE COVERAGE ===\n');
    allLobesUNK = [lobes, {'UNK'}];
    for k = 1:numel(allLobesUNK)
        fprintf('%3s : %d\n', allLobesUNK{k}, sum(chanLobe == allLobesUNK{k}));
    end

    % ---- Subject-level CHANGE matrices ----
    dR_stack = NaN(nCh, nCh, nSub);
    kept     = false(nSub, 1);

    for s = 1:nSub
        subj = group_data(s);
        sid  = subj.demographics('SubjectID');

        if ~ismember(silKey, subj.stimulus.keys)
            warning('%s: missing silence trigger', sid);
            continue;
        end

        on = subj.stimulus(silKey).onset;
        if numel(on) < 2
            warning('%s: <2 silence blocks', sid);
            continue;
        end

        t = subj.time(:);
        X = subj.data(:, isHbO);

        idxPre  = t >= on(1) & t < on(1) + durSec;
        idxPost = t >= on(2) & t < on(2) + durSec;

        if nnz(idxPre) < 10 || nnz(idxPost) < 10
            warning('%s: insufficient samples', sid);
            continue;
        end

        Rpre  = corr(X(idxPre,:),  'Rows','pairwise');
        Rpost = corr(X(idxPost,:), 'Rows','pairwise');

        Rpre  = max(min(Rpre,  0.999999), -0.999999);
        Rpost = max(min(Rpost, 0.999999), -0.999999);

        Zpre  = atanh(Rpre);
        Zpost = atanh(Rpost);
        dZ    = Zpost - Zpre;
        dR    = tanh(dZ);

        dR(1:nCh+1:end) = NaN;

        dR_stack(:,:,s) = dR;
        kept(s) = true;

        fprintf('✓ %s\n', sid);
    end

    dR_stack = dR_stack(:,:,kept);
    nKept    = size(dR_stack, 3);
    fprintf('Kept subjects = %d\n', nKept);

    % ---- Group mean CHANGE ----
    dR_avg = mean(dR_stack, 3, 'omitnan');

    % ---- Large-change threshold on group-average ----
    maskUT  = triu(true(nCh), 1);
    absVals = abs(dR_avg(maskUT));
    mu      = mean(absVals, 'omitnan');
    sd      = std(absVals,  'omitnan');
    thr     = mu + sd;

    strongMask = maskUT & (abs(dR_avg) > thr);

    fprintf('Large-change threshold (|Dr| > mean+SD): %.4f\n', thr);
    fprintf('Kept edges = %d (out of %d)\n', nnz(strongMask), nnz(maskUT));

    % ---- Edge list ----
    [ii, jj] = find(strongMask);
    dR_edges = dR_avg(strongMask);

    Tedges = table(ii, jj, dR_edges, ...
        'VariableNames', {'chA','chB','dR'});

    Tedges.lobeA = chanLobe(Tedges.chA);
    Tedges.lobeB = chanLobe(Tedges.chB);

    Tok = Tedges(~(Tedges.lobeA=="UNK" | Tedges.lobeB=="UNK"), :);
    fprintf('Selected edges (all) = %d | after dropping UNK = %d\n', ...
        height(Tedges), height(Tok));

    % ---- Lobe-pair summary ----
    if height(Tok) == 0
        warning('No edges survived threshold for %s.', groupName);
        Tsum = table();
    else
        Tok.lobeA = categorical(Tok.lobeA, lobes);
        Tok.lobeB = categorical(Tok.lobeB, lobes);

        a = string(Tok.lobeA); b = string(Tok.lobeB);
        swap = a > b;
        a(swap) = string(Tok.lobeB(swap));
        b(swap) = string(Tok.lobeA(swap));
        Tok.lobeA2 = categorical(a, lobes);
        Tok.lobeB2 = categorical(b, lobes);

        G      = findgroups(Tok.lobeA2, Tok.lobeB2);
        nEdge  = splitapply(@numel, Tok.dR, G);
        meanDR = splitapply(@(x) mean(x,'omitnan'), Tok.dR, G);

        lobeA_u = splitapply(@(x) x(1), Tok.lobeA2, G);
        lobeB_u = splitapply(@(x) x(1), Tok.lobeB2, G);

        Tsum = table(string(lobeA_u), string(lobeB_u), nEdge, meanDR, ...
            'VariableNames', {'lobeA','lobeB','n','mean_dR'});

        Tsum.pct = 100 * Tsum.n / sum(Tsum.n);
        Tsum = sortrows(Tsum, 'n', 'descend');

        disp(Tsum(1:min(15,height(Tsum)), :));
    end

    % ---- Lobe x Lobe mean dR matrix ----
    LobeMat = NaN(nL, nL);
    for aL = 1:nL
        for bL = 1:nL
            idxA = chanLobe == lobes{aL};
            idxB = chanLobe == lobes{bL};
            subM = dR_avg(idxA, idxB);
            if aL == bL
                subM = subM(triu(true(sum(idxA)), 1));
            end
            LobeMat(aL,bL) = mean(subM(:), 'omitnan');
        end
    end

    % ---- Visualizations ----
    figure('Color','w');
    imagesc(dR_avg);
    axis image; colorbar;
    title([upper(groupName) ': Whole-brain \Deltar (post-pre)'], ...
          'FontSize', 14, 'FontWeight', 'bold');
    xlabel('Channel'); ylabel('Channel');

    figure('Color','w');
    imagesc(LobeMat);
    axis image; colorbar;
    set(gca,'XTick',1:nL,'XTickLabel',lobes,'XTickLabelRotation',45, ...
            'YTick',1:nL,'YTickLabel',lobes);
    title([upper(groupName) ': Lobe x Lobe mean \Deltar'], ...
          'FontSize', 14, 'FontWeight', 'bold');

    % ---- Save per group ----
    tag = lower(groupName);
    assignin('base', ['dR_avg_'      tag], dR_avg);
    assignin('base', ['thr_'         tag], thr);
    assignin('base', ['strongMask_'  tag], strongMask);
    assignin('base', ['Tedges_'      tag], Tedges);
    assignin('base', ['Tok_'         tag], Tok);
    assignin('base', ['Tsum_'        tag], Tsum);
    assignin('base', ['LobeMat_'     tag], LobeMat);

    if strcmpi(groupName, 'Control')
        dR_stack_control = dR_stack;
    else
        dR_stack_exp = dR_stack;
    end

    fprintf('✓ %s done.\n', upper(groupName));
end



% BETWEEN-GROUP STATS ON LOBE-PAIR LEVEL


fprintf('\n================ BETWEEN-GROUP LOBE-PAIR STATS ================\n');

chanLobe = chanLobe_global;
alpha    = 0.05;
results  = table();

for aL = 1:nL
    for bL = aL:nL

        lobeA = lobes{aL};
        lobeB = lobes{bL};

        idxA = find(chanLobe == lobeA);
        idxB = find(chanLobe == lobeB);

        if isempty(idxA) || isempty(idxB)
            continue;
        end

        % Subject-level mean dR for Control
        dR_lobepair_ctrl = NaN(size(dR_stack_control, 3), 1);
        for s = 1:size(dR_stack_control, 3)
            subM = dR_stack_control(idxA, idxB, s);
            if aL == bL
                subM = subM(triu(true(numel(idxA)), 1));
            end
            dR_lobepair_ctrl(s) = mean(subM(:), 'omitnan');
        end

        % Subject-level mean dR for Experimental
        dR_lobepair_exp = NaN(size(dR_stack_exp, 3), 1);
        for s = 1:size(dR_stack_exp, 3)
            subM = dR_stack_exp(idxA, idxB, s);
            if aL == bL
                subM = subM(triu(true(numel(idxA)), 1));
            end
            dR_lobepair_exp(s) = mean(subM(:), 'omitnan');
        end

        dR_ctrl = dR_lobepair_ctrl(isfinite(dR_lobepair_ctrl));
        dR_exp  = dR_lobepair_exp(isfinite(dR_lobepair_exp));

        if numel(dR_ctrl) < 2 || numel(dR_exp) < 2
            continue;
        end

        [~, pC,  ~, sC]  = ttest(dR_ctrl,  0);
        [~, pE,  ~, sE]  = ttest(dR_exp,   0);
        [~, pBG, ~, sBG] = ttest2(dR_ctrl, dR_exp);

        fprintf('\n[%s - %s]\n', lobeA, lobeB);
        fprintf('  CTRL  n=%d | mean Dr=%.4f | t=%.3f | p=%.4g\n', ...
            numel(dR_ctrl), mean(dR_ctrl,'omitnan'), sC.tstat, pC);
        fprintf('  TINN  n=%d | mean Dr=%.4f | t=%.3f | p=%.4g\n', ...
            numel(dR_exp),  mean(dR_exp,'omitnan'),  sE.tstat, pE);
        fprintf('  CTRL vs TINN: t=%.3f | p=%.4g', sBG.tstat, pBG);
        if pBG < alpha, fprintf('  **'); end
        fprintf('\n');

        newRow = table({lobeA}, {lobeB}, ...
            numel(dR_ctrl), mean(dR_ctrl,'omitnan'), sC.tstat, pC, ...
            numel(dR_exp),  mean(dR_exp,'omitnan'),  sE.tstat, pE, ...
            sBG.tstat, pBG, pBG < alpha, ...
            'VariableNames', {'lobeA','lobeB', ...
            'n_ctrl','mean_dR_ctrl','t_ctrl','p_ctrl', ...
            'n_exp', 'mean_dR_exp', 't_exp', 'p_exp', ...
            't_between','p_between','sig'});

        results = [results; newRow];
    end
end

results = sortrows(results, 'p_between', 'ascend');

fprintf('\n================ SIGNIFICANT LOBE PAIRS (p < %.2f) ================\n', alpha);
disp(results(results.sig == true, :));

fprintf('\n================ FULL RESULTS TABLE ================\n');
disp(results);

assignin('base', 'lobe_pair_stats', results);
fprintf('\n✓ Between-group lobe-pair stats done.\n');


%% BUILD LeftMat and RightMat (required by your Fig 6-style plotting)
% Uses your already-computed lobe-wise matrices:
%   Zlob_pre_ctrl,  Zlob_post_ctrl
%   Zlob_pre_exp,   Zlob_post_exp


% --- sanity checks ---
assert(exist('Zlob_pre_ctrl','var')==1,  'Missing Zlob_pre_ctrl');
assert(exist('Zlob_post_ctrl','var')==1, 'Missing Zlob_post_ctrl');
assert(exist('Zlob_pre_exp','var')==1,   'Missing Zlob_pre_exp');
assert(exist('Zlob_post_exp','var')==1,  'Missing Zlob_post_exp');

% --- Mean across subjects in Z-space (paper) ---
meanZ_ctrl_pre  = nanmean(Zlob_pre_ctrl,  1);   % 1x9
meanZ_ctrl_post = nanmean(Zlob_post_ctrl, 1);   % 1x9
meanZ_exp_pre   = nanmean(Zlob_pre_exp,   1);   % 1x9
meanZ_exp_post  = nanmean(Zlob_post_exp,  1);   % 1x9

% --- Convert to Pearson r for PRE/POST display ---
Ctrl_pre  = tanh(meanZ_ctrl_pre);      % 1x9
Ctrl_post = tanh(meanZ_ctrl_post);     % 1x9
Tinn_pre  = tanh(meanZ_exp_pre);       % 1x9
Tinn_post = tanh(meanZ_exp_post);      % 1x9

% --- PAPER change (Δr) = tanh( meanZ_post - meanZ_pre ) ---
Ctrl_change = tanh(meanZ_ctrl_post - meanZ_ctrl_pre);   % 1x9
Tinn_change = tanh(meanZ_exp_post  - meanZ_exp_pre);    % 1x9

% --- Matrices expected by your plotting block ---
LeftMat  = [Tinn_post; Tinn_pre; Ctrl_post; Ctrl_pre];      % 4x9
RightMat = [Tinn_change; Ctrl_change];                      % 2x9

% --- optional: print quick check ---
disp('LeftMat size:'), disp(size(LeftMat))
disp('RightMat size:'), disp(size(RightMat))

%% Fig 8a and 8b 
% generate_Fig8_PLOS.m
% PLOS ONE compliant 2-panel figure: ROI connectivity heatmaps across lobes
%
% Panels (2 rows, stacked):
%   (a) Mean ROI-seed connectivity (Pearson r) — 4 rows x 9 lobes
%   (b) Change in connectivity (Δr = tanh(ΔZ)) — 2 rows x 9 lobes
%
% PREREQUISITE: workspace contains:
%   LeftMat   = [Tinn_post; Tinn_pre; Ctrl_post; Ctrl_pre]   (4 x 9)
%   RightMat  = [Tinn_change; Ctrl_change]                   (2 x 9)
%   lobes     = {'F','FT','FP','TP','T','P','OP','OT','O'}


% ----- Safety checks -----
if ~exist('LeftMat','var') || isempty(LeftMat)
    error('LeftMat not found in workspace.');
end
if ~exist('RightMat','var') || isempty(RightMat)
    error('RightMat not found in workspace.');
end
if ~exist('lobes','var') || isempty(lobes)
    lobes = {'F','FT','FP','TP','T','P','OP','OT','O'};
end
if istable(LeftMat),  LeftMat  = table2array(LeftMat);  end
if istable(RightMat), RightMat = table2array(RightMat); end

assert(isequal(size(LeftMat),  [4 9]), 'LeftMat must be 4x9.');
assert(isequal(size(RightMat), [2 9]), 'RightMat must be 2x9.');

% ----- Output -----
outDir = FIG_DIR;
if ~exist(outDir, 'dir'); mkdir(outDir); end
outFile = fullfile(outDir, 'Fig8.tif');

% ----- Figure size -----
% (a) 4x9 + (b) 2x9 stacked. With square cells and shared x-axis, this needs
% room for both. Full page width 7.5"; height ~5.5" gives good cell aspect.
figWidthIn  = 7.5;
figHeightIn = 5.5;
dpi         = 300;

% ----- Font settings (Arial 8-12 pt per PLOS) -----
axFont    = 'Arial';
axSize    = 10;     % tick labels
labSize   = 11;     % axis labels (we have none, but in colorbar)
panelSize = 12;     % (a), (b) panel letters
numSize   = 8;      % numbers inside cells (must be small to fit)
cbSize    = 9;      % colorbar tick labels

% ----- Row labels -----
rowLabelsLeft  = {'Tinnitus post','Tinnitus pre','Control post','Control pre'};
rowLabelsRight = {'Tinnitus change','Control change'};

% ----- Color limits and maps -----
climLeft  = [0 1];          % Pearson r: 0 to max (sequential)
climRight = [-0.2 0.2];     % Δr: diverging around 0
cmapLeft  = jet(256);
cmapRight = jet(256);       % keeping jet for consistency with your code


% Build figure

fig = figure('Units','inches', ...
             'Position',[1 1 figWidthIn figHeightIn], ...
             'Color','w', ...
             'PaperUnits','inches', ...
             'PaperPosition',[0 0 figWidthIn figHeightIn], ...
             'PaperSize',[figWidthIn figHeightIn]);

% ----- Layout -----
% Panel (a): 4 rows, taller
% Panel (b): 2 rows, shorter
% Both share x-axis position so columns align visually.
%
% Heights apportioned roughly by row count: (a) gets ~2x the height of (b).

leftMargin   = 0.16;     % wider for long row labels ("Tinnitus post")
rightMargin  = 0.12;     % colorbar space
bottomMargin = 0.10;
topMargin    = 0.06;
vGap         = 0.10;     % gap between panel a and b

% Total height available for both panels
totalH = 1 - bottomMargin - topMargin - vGap;

% Apportion height by row count: 4 + 2 = 6 row-units
heightA = totalH * (4/6);
heightB = totalH * (2/6);

% Width is shared
panelW  = 1 - leftMargin - rightMargin;

% Panel positions
posA = [leftMargin, 1 - topMargin - heightA, panelW, heightA];
posB = [leftMargin, posA(2) - vGap - heightB, panelW, heightB];


% PANEL (a): Pearson r heatmap (4x9)

ax1 = axes('Parent', fig, 'Position', posA);
imagesc(ax1, LeftMat);
set(ax1, 'YDir','reverse');
set(ax1, 'TickLength',[0 0]);

% Tick labels
set(ax1, 'XTick', 1:numel(lobes), 'XTickLabel', lobes, ...
         'YTick', 1:numel(rowLabelsLeft), 'YTickLabel', rowLabelsLeft, ...
         'FontName', axFont, 'FontSize', axSize);

xlim(ax1, [0.5 numel(lobes)+0.5]);
ylim(ax1, [0.5 numel(rowLabelsLeft)+0.5]);

colormap(ax1, cmapLeft);
caxis(ax1, climLeft);

% Colorbar for (a)
cbA = colorbar(ax1, 'eastoutside');
ylabel(cbA, 'r', 'FontName',axFont, 'FontSize',labSize, 'FontWeight','bold');
set(cbA, 'FontName',axFont, 'FontSize',cbSize);

% Grid lines
hold(ax1, 'on');
[nR, nC] = size(LeftMat);
for k = 0:nC
    plot(ax1, [k+0.5 k+0.5], [0.5 nR+0.5], 'k-', 'LineWidth', 0.8);
end
for k = 0:nR
    plot(ax1, [0.5 nC+0.5], [k+0.5 k+0.5], 'k-', 'LineWidth', 0.8);
end

% Numbers in cells
for r = 1:nR
    for c = 1:nC
        v = LeftMat(r,c);
        if isfinite(v)
            % Choose text color for readability against background
            normVal = (v - climLeft(1)) / (climLeft(2) - climLeft(1));
            if normVal > 0.6 && normVal < 0.8
                txtColor = 'k';
            else
                txtColor = 'k';
            end
            text(ax1, c, r, sprintf('%.2f', v), ...
                'HorizontalAlignment','center', ...
                'VerticalAlignment','middle', ...
                'FontName',axFont, 'FontSize',numSize, ...
                'FontWeight','bold', 'Color',txtColor);
        end
    end
end

% Panel label (a)
text(ax1, -0.13, 1.08, '(a)', ...
    'Units','normalized', ...
    'FontName',axFont, 'FontSize',panelSize, 'FontWeight','bold');


% PANEL (b): Δr heatmap (2x9)

ax2 = axes('Parent', fig, 'Position', posB);
imagesc(ax2, RightMat);
set(ax2, 'YDir','reverse');
set(ax2, 'TickLength',[0 0]);

set(ax2, 'XTick', 1:numel(lobes), 'XTickLabel', lobes, ...
         'YTick', 1:numel(rowLabelsRight), 'YTickLabel', rowLabelsRight, ...
         'FontName', axFont, 'FontSize', axSize);

xlim(ax2, [0.5 numel(lobes)+0.5]);
ylim(ax2, [0.5 numel(rowLabelsRight)+0.5]);

colormap(ax2, cmapRight);
caxis(ax2, climRight);

% Colorbar for (b)
cbB = colorbar(ax2, 'eastoutside');
ylabel(cbB, '\Deltar', 'FontName',axFont, 'FontSize',labSize, 'FontWeight','bold');
set(cbB, 'FontName',axFont, 'FontSize',cbSize);

% Grid lines
hold(ax2, 'on');
[nR, nC] = size(RightMat);
for k = 0:nC
    plot(ax2, [k+0.5 k+0.5], [0.5 nR+0.5], 'k-', 'LineWidth', 0.8);
end
for k = 0:nR
    plot(ax2, [0.5 nC+0.5], [k+0.5 k+0.5], 'k-', 'LineWidth', 0.8);
end

% Numbers in cells
for r = 1:nR
    for c = 1:nC
        v = RightMat(r,c);
        if isfinite(v)
            text(ax2, c, r, sprintf('%.2f', v), ...
                'HorizontalAlignment','center', ...
                'VerticalAlignment','middle', ...
                'FontName',axFont, 'FontSize',numSize, ...
                'FontWeight','bold', 'Color','k');
        end
    end
end

% Panel label (b)
text(ax2, -0.13, 1.15, '(b)', ...
    'Units','normalized', ...
    'FontName',axFont, 'FontSize',panelSize, 'FontWeight','bold');


% Export as PLOS-compliant TIFF

% delete(findall(fig,'Type','hggroup'));
% delete(findall(fig,'Tag','DataTipMarker'));
% 
% try
%     exportgraphics(fig, outFile, 'Resolution', dpi, 'BackgroundColor','white');
%     fprintf('Saved (exportgraphics): %s\n', outFile);
% catch
%     print(fig, outFile, '-dtiff', sprintf('-r%d', dpi));
%     fprintf('Saved (print): %s\n', outFile);
% end
% 
% % =========================================================================
% % Compliance check
% % =========================================================================
% info = imfinfo(outFile);
% fprintf('\n--- Fig8.tif specs ---\n');
% fprintf('  Width:       %d px\n', info.Width);
% fprintf('  Height:      %d px\n', info.Height);
% fprintf('  XRes:        %g dpi\n', info.XResolution);
% fprintf('  YRes:        %g dpi\n', info.YResolution);
% fprintf('  ColorType:   %s, BitDepth: %d\n', info.ColorType, info.BitDepth);
% fprintf('  Compression: %s\n', info.Compression);
% fprintf('  File size:   %.2f MB\n', info.FileSize / 1024 / 1024);
% 
% fprintf('\n--- PLOS ONE compliance check ---\n');
% ok = true;
% if info.Width < 789 || info.Width > 2250
%     fprintf('  [WARN] Width %d outside 789-2250 px range\n', info.Width); ok = false;
% end
% if info.Height > 2625
%     fprintf('  [WARN] Height %d exceeds 2625 px max\n', info.Height); ok = false;
% end
% if info.XResolution < 300 || info.XResolution > 600
%     fprintf('  [WARN] DPI %g outside 300-600 range\n', info.XResolution); ok = false;
% end
% if (info.FileSize / 1024 / 1024) > 10
%     fprintf('  [WARN] File size exceeds 10 MB\n'); ok = false;
% end
% if ok
%     fprintf('  PASS: dimensions, resolution, and file size all within PLOS specs.\n');
% end
% fprintf('  NOTE: If compression is not LZW, re-save in GIMP/Photoshop with LZW.\n');


%% generate_Fig9_PLOS.m
% PLOS ONE compliant version of the triangular lobe-pair connectivity figure.
% Upper-left triangle = Control, lower-right triangle = Tinnitus.
% Circle COLOR = mean Δr per lobe-pair, Circle SIZE = % of selected edges.
%
% PREREQUISITE: Run your full Fig9 analysis pipeline first so the workspace contains:
%   HaemoData_hbo
%   Tok_ctrl (or Tok_control / Tok_Control)
%   Tok_exp  (or Tok_experimental / Tok_tinn)
% with Tok tables containing: chA, chB, lobeA, lobeB, and one of {dZ, meanZ, dR}


% ----- Output -----
outDir = FIG_DIR;
if ~exist(outDir, 'dir'); mkdir(outDir); end
outFile = fullfile(outDir, 'Fig9.tif');
% ----- Figure size: square-ish, fits well within full page

figWidthIn  = 6.5;
figHeightIn = 6.5;
dpi         = 300;

% ----- Font settings 
axFont    = 'Arial';
axSize    = 10;     % tick labels (was 14)
labSize   = 11;     % axis labels (was 14)
sideSize  = 9;      % top/right reversed labels (was 12)
cbSize    = 9;      % colorbar tick labels


%  Lobe order (paper)

lobes = {'F','FT','FP','TP','T','P','OP','OT','O'};
nL    = numel(lobes);
nGrid = nL + 1;
idxOf = containers.Map(lobes, num2cell(1:nL));

% 1) Grab Tok tables from workspace 

if exist('Tok_ctrl','var'); TokC = Tok_ctrl;
elseif exist('Tok_control','var'); TokC = Tok_control;
elseif exist('Tok_Control','var'); TokC = Tok_Control;
else
    error('Could not find Tok_ctrl or Tok_control in workspace.');
end

if exist('Tok_exp','var'); TokT = Tok_exp;
elseif exist('Tok_experimental','var'); TokT = Tok_experimental;
elseif exist('Tok_tinn','var'); TokT = Tok_tinn;
else
    error('Could not find Tok_exp or Tok_experimental in workspace.');
end

needCols = {'chA','chB'};
for k = 1:numel(needCols)
    if ~ismember(needCols{k}, TokC.Properties.VariableNames)
        error('Control Tok missing column: %s', needCols{k});
    end
    if ~ismember(needCols{k}, TokT.Properties.VariableNames)
        error('Tinnitus/Experimental Tok missing column: %s', needCols{k});
    end
end

% Edge Z values
if ismember('dZ', TokC.Properties.VariableNames)
    edgeZ_C = TokC.dZ;
elseif ismember('meanZ', TokC.Properties.VariableNames)
    edgeZ_C = TokC.meanZ;
elseif ismember('dR', TokC.Properties.VariableNames)
    edgeZ_C = atanh(max(min(TokC.dR,0.999999),-0.999999));
else
    error('Control Tok must have dZ OR dR OR meanZ.');
end


if ismember('dZ', TokT.Properties.VariableNames)
    edgeZ_T = TokT.dZ;
elseif ismember('meanZ', TokT.Properties.VariableNames)
    edgeZ_T = TokT.meanZ;
elseif ismember('dR', TokT.Properties.VariableNames)
    edgeZ_T = atanh(max(min(TokT.dR,0.999999),-0.999999));
else
    error('Tinnitus Tok must have dZ OR dR OR meanZ.');
end

if ~ismember('lobeA', TokC.Properties.VariableNames) || ~ismember('lobeB', TokC.Properties.VariableNames)
    error('Control Tok must contain lobeA and lobeB.');
end
if ~ismember('lobeA', TokT.Properties.VariableNames) || ~ismember('lobeB', TokT.Properties.VariableNames)
    error('Tinnitus Tok must contain lobeA and lobeB.');
end

TokC.lobeA = string(TokC.lobeA); TokC.lobeB = string(TokC.lobeB);
TokT.lobeA = string(TokT.lobeA); TokT.lobeB = string(TokT.lobeB);


% 2) Hemisphere labels for isBilateral (kept for table only)

link0 = HaemoData_hbo(1).probe.link;
isHbO = strcmpi(link0.type,'hbo');

if ~ismember('ROI_full', link0.Properties.VariableNames)
    error('HaemoData_hbo(1).probe.link.ROI_full not found.');
end

roi_full_hbo = string(link0.ROI_full(isHbO));
nHbO = numel(roi_full_hbo);

chanHemi = strings(nHbO,1);
chanHemi(:) = "UNK";
chanHemi(endsWith(roi_full_hbo,"_L")) = "L";
chanHemi(endsWith(roi_full_hbo,"_R")) = "R";

if any(TokC.chA<1 | TokC.chA>nHbO | TokC.chB<1 | TokC.chB>nHbO)
    error('TokC channel indices out of range for HbO channels (%d).', nHbO);
end
if any(TokT.chA<1 | TokT.chA>nHbO | TokT.chB<1 | TokT.chB>nHbO)
    error('TokT channel indices out of range for HbO channels (%d).', nHbO);
end

TokC.isBilateral = (chanHemi(TokC.chA) ~= chanHemi(TokC.chB)) ...
                 & (chanHemi(TokC.chA) ~= "UNK") ...
                 & (chanHemi(TokC.chB) ~= "UNK");

TokT.isBilateral = (chanHemi(TokT.chA) ~= chanHemi(TokT.chB)) ...
                 & (chanHemi(TokT.chA) ~= "UNK") ...
                 & (chanHemi(TokT.chB) ~= "UNK");


% 3) Normalize lobe-pair ordering

aC = TokC.lobeA; bC = TokC.lobeB;
aT = TokT.lobeA; bT = TokT.lobeB;

swapC = aC > bC; tmp = aC(swapC); aC(swapC)=bC(swapC); bC(swapC)=tmp;
swapT = aT > bT; tmp = aT(swapT); aT(swapT)=bT(swapT); bT(swapT)=tmp;

TokC.lobeA2 = categorical(aC, lobes);
TokC.lobeB2 = categorical(bC, lobes);
TokT.lobeA2 = categorical(aT, lobes);
TokT.lobeB2 = categorical(bT, lobes);


% 4) Build summary tables

S_ctrl = buildSummary(TokC, edgeZ_C, lobes);
S_tinn = buildSummary(TokT, edgeZ_T, lobes);


% 5) Colormap + limits

cmap = jet(256);
clim = 0.1;


% 6) Circle size scaling

minRad = 0.06;
maxRad = 0.45;
pctMax = max([S_ctrl.pct(:); S_tinn.pct(:)]);
if ~isfinite(pctMax) || pctMax <= 0, pctMax = 1; end
pct_to_rad = @(p) (minRad + (maxRad-minRad) * sqrt(max(p,0)/pctMax));


% 7) DRAW FIGURE (PLOS-compliant)

fig = figure('Units','inches', ...
             'Position',[1 1 figWidthIn figHeightIn], ...
             'Color','w', ...
             'PaperUnits','inches', ...
             'PaperPosition',[0 0 figWidthIn figHeightIn], ...
             'PaperSize',[figWidthIn figHeightIn]);

% Axes positioned to leave room for colorbar on right
ax = axes('Parent', fig, 'Position',[0.11 0.10 0.64 0.80]);
hold(ax, 'on'); axis(ax,'equal','ij');
xlim(ax, [0 nGrid]); ylim(ax, [0 nGrid]);

set(ax, ...
    'XTick', (0.5:1:nL-0.5), ...
    'YTick', (1.5:1:nL+0.5), ...
    'XTickLabel', lobes, ...
    'YTickLabel', lobes, ...
    'XTickLabelRotation', 45, ...
    'FontName', axFont, ...
    'FontSize', axSize, ...
    'LineWidth', 1.0, ...
    'TickLength', [0 0], ...
    'XColor', 'k', ...
    'YColor', 'k');

xlabel(ax, 'Tinnitus', 'FontName',axFont, 'FontWeight','bold', 'FontSize',labSize);
ylabel(ax, 'Tinnitus', 'FontName',axFont, 'FontWeight','bold', 'FontSize',labSize);

colormap(ax, cmap);
caxis(ax, [-clim clim]);

% Colorbar
cb = colorbar(ax, 'Position', [0.90 0.10 0.025 0.80]);
ylabel(cb, '\Delta r (Post–Pre)', ...
    'FontName',axFont, 'FontWeight','bold', 'FontSize',cbSize);
set(cb, 'FontName',axFont, 'FontSize',cbSize-1);

% Diagonal gray blocks
for i = 1:nL
    rectangle(ax,'Position',[i i 1 1], ...
        'FaceColor',[0.7 0.7 0.7], 'EdgeColor','none');
end
% Top-left padding square
rectangle(ax,'Position',[0 0 1 1], ...
    'FaceColor',[0.7 0.7 0.7], 'EdgeColor','none');

% Grid lines
for i = 0:nGrid
    plot(ax,[i i],[0 nGrid],'k-','LineWidth',0.5);
    plot(ax,[0 nGrid],[i i],'k-','LineWidth',0.5);
end

% Top-row reversed labels (for Control quadrant)
for i = 1:nL
    text(ax, i+0.5, -0.35, lobes{nL - i + 1}, ...
        'HorizontalAlignment','center', ...
        'VerticalAlignment','middle', ...
        'FontName',axFont, 'FontSize',sideSize, ...
        'FontWeight','bold', 'Rotation',90);
end
% Right-column reversed labels (for Control quadrant)
for i = 1:nL
    text(ax, nGrid+0.25, i-0.5, lobes{nL - i + 1}, ...
        'HorizontalAlignment','left', ...
        'VerticalAlignment','middle', ...
        'FontName',axFont, 'FontSize',sideSize, ...
        'FontWeight','bold');
end


% 8) Draw circles (color=meanR, size=pct)

drawCircleSimple = @(cx,cy,rad,val) drawCircleFilled( ...
    ax, cx, cy, rad, valToRGB(val, cmap, [-clim clim]) );

% TINNITUS triangle (bottom-right)
for r = 1:height(S_tinn)
    a = S_tinn.lobeA(r); b = S_tinn.lobeB(r);
    if ~isKey(idxOf, char(a)) || ~isKey(idxOf, char(b)), continue; end
    ia = idxOf(char(a)); ib = idxOf(char(b));

    i = max(ia,ib);
    j = min(ia,ib);
    if i==j, continue; end

    cx = (j-1) + 0.5;
    cy = i + 0.5;

    rad = pct_to_rad(S_tinn.pct(r));
    drawCircleSimple(cx, cy, rad, S_tinn.meanR(r));
end

% CONTROL triangle (top-left, reversed mapping)
for r = 1:height(S_ctrl)
    a = S_ctrl.lobeA(r); b = S_ctrl.lobeB(r);
    if ~isKey(idxOf, char(a)) || ~isKey(idxOf, char(b)), continue; end
    ia = idxOf(char(a)); ib = idxOf(char(b));
    if ia==ib, continue; end

    row = nL - ib + 1;
    col = nL - ia + 1;

    if row >= col
        row = nL - ia + 1;
        col = nL - ib + 1;
    end
    if row >= col, continue; end

    cx = col + 0.5;
    cy = (row-1) + 0.5;

    rad = pct_to_rad(S_ctrl.pct(r));
    drawCircleSimple(cx, cy, rad, S_ctrl.meanR(r));
end


% 9) Save outputs to base workspace

assignin('base','Tsum_control',S_ctrl);
assignin('base','Tsum_experimental',S_tinn);

fprintf('\n✓ Fig 9 figure built.\n');
fprintf('  S_ctrl (Tsum_control) rows=%d\n', height(S_ctrl));
fprintf('  S_tinn (Tsum_experimental) rows=%d\n', height(S_tinn));


% 10) Export as PLOS-compliant TIFF

% delete(findall(fig,'Type','hggroup'));
% delete(findall(fig,'Tag','DataTipMarker'));
% 
% try
%     exportgraphics(fig, outFile, 'Resolution', dpi, 'BackgroundColor','white');
%     fprintf('Saved (exportgraphics): %s\n', outFile);
% catch
%     print(fig, outFile, '-dtiff', sprintf('-r%d', dpi));
%     fprintf('Saved (print): %s\n', outFile);
% end
% 
% 
% % 11) Compliance check
% 
% info = imfinfo(outFile);
% fprintf('\n--- Fig9.tif specs ---\n');
% fprintf('  Width:       %d px\n', info.Width);
% fprintf('  Height:      %d px\n', info.Height);
% fprintf('  XRes:        %g dpi\n', info.XResolution);
% fprintf('  YRes:        %g dpi\n', info.YResolution);
% fprintf('  ColorType:   %s, BitDepth: %d\n', info.ColorType, info.BitDepth);
% fprintf('  Compression: %s\n', info.Compression);
% fprintf('  File size:   %.2f MB\n', info.FileSize / 1024 / 1024);
% 
% fprintf('\n--- PLOS ONE compliance check ---\n');
% ok = true;
% if info.Width < 789 || info.Width > 2250
%     fprintf('  [WARN] Width %d outside 789-2250 px range\n', info.Width); ok = false;
% end
% if info.Height > 2625
%     fprintf('  [WARN] Height %d exceeds 2625 px max\n', info.Height); ok = false;
% end
% if info.XResolution < 300 || info.XResolution > 600
%     fprintf('  [WARN] DPI %g outside 300-600 range\n', info.XResolution); ok = false;
% end
% if (info.FileSize / 1024 / 1024) > 10
%     fprintf('  [WARN] File size exceeds 10 MB\n'); ok = false;
% end
% if ok
%     fprintf('  PASS: dimensions, resolution, and file size all within PLOS specs.\n');
% end
% fprintf('  NOTE: If compression is not LZW, re-save in GIMP/Photoshop with LZW.\n');



%% WHOLE-BRAIN RSFC (PRE / POST / DELTA) — TINNITUS GROUP

silKey    = '1';
durSec    = 180;
groupName = 'Experimental';

group_idx  = strcmpi(tbl_demo.group, groupName);
group_data = HaemoData_hbo(group_idx);
nSub       = numel(group_data);

SubjectID  = strings(nSub, 1);
RSFC_pre   = NaN(nSub, 1);
RSFC_post  = NaN(nSub, 1);
RSFC_delta = NaN(nSub, 1);

% HbO channel mask (locked to first subject)
link0 = group_data(1).probe.link;
isHbO = strcmpi(link0.type, 'hbo');
nCh   = sum(isHbO);
UT    = triu(true(nCh), 1);   % upper triangle mask

for s = 1:nSub

    subj = group_data(s);
    sid  = subj.demographics('SubjectID');
    SubjectID(s) = string(sid);

    if ~ismember(silKey, subj.stimulus.keys)
        warning('%s: missing silence trigger', sid);
        continue;
    end

    on = subj.stimulus(silKey).onset;
    if numel(on) < 2
        warning('%s: fewer than 2 silence blocks', sid);
        continue;
    end

    t = subj.time(:);
    X = subj.data(:, isHbO);

    idxPre  = t >= on(1) & t < on(1) + durSec;
    idxPost = t >= on(2) & t < on(2) + durSec;

    if nnz(idxPre) < 10 || nnz(idxPost) < 10
        warning('%s: insufficient samples', sid);
        continue;
    end

    Rpre  = corr(X(idxPre,:),  'Rows','pairwise');
    Rpost = corr(X(idxPost,:), 'Rows','pairwise');

    Rpre  = max(min(Rpre,  0.999999), -0.999999);
    Rpost = max(min(Rpost, 0.999999), -0.999999);

    Zpre  = atanh(Rpre);
    Zpost = atanh(Rpost);
    dZ    = Zpost - Zpre;

    RSFC_pre(s)   = tanh(mean(Zpre(UT),  'omitnan'));
    RSFC_post(s)  = tanh(mean(Zpost(UT), 'omitnan'));
    RSFC_delta(s) = tanh(mean(dZ(UT),    'omitnan'));

    fprintf('✓ %s\n', sid);
end

% --- Build table ---
RSFC_WholeBrain = table(SubjectID, RSFC_pre, RSFC_post, RSFC_delta);

% --- Drop subjects that failed (all NaN rows) ---
failed = ~all(isfinite(RSFC_WholeBrain{:, 2:end}), 2);
if any(failed)
    fprintf('\nDropped from whole-brain (failed): %d subject(s)\n', sum(failed));
    disp(RSFC_WholeBrain.SubjectID(failed));
end
RSFC_WholeBrain = RSFC_WholeBrain(~failed, :);

disp(RSFC_WholeBrain);


%%  ASSEMBLE AND SAVE RSFC OUTPUTS
%
%  Everything downstream analyses need is written to a single struct, RSFC,
%  in RSFC_outputs.mat. Rows of every subject-level field are in the order
%  given by RSFC.SubjectID.
% =========================================================================

idxExp = strcmpi(tbl_demo.group, 'Experimental');

RSFC = struct();
RSFC.SubjectID = string(tbl_demo.SubjectID(idxExp));
RSFC.Lobes     = lobes;

% ---- ROI-seed connectivity, Fisher Z (tinnitus group) ----
RSFC.ROI.Z_pre  = zPre_exp_all;
RSFC.ROI.Z_post = zPost_exp_all;
RSFC.ROI.Z_diff = zPost_exp_all - zPre_exp_all;
RSFC.ROI.r_diff = tanh(RSFC.ROI.Z_diff);

% ---- ROI-seed-to-lobe connectivity, Fisher Z (nSub x 9) ----
RSFC.Lobe.Z_pre  = Zlob_pre_exp;
RSFC.Lobe.Z_post = Zlob_post_exp;
RSFC.Lobe.Z_diff = Zlob_diff_exp;

% ---- Whole-brain connectivity, Pearson r (tinnitus group) ----
RSFC.WholeBrain = RSFC_WholeBrain;

% ---- Control-group equivalents, for between-group work downstream ----
RSFC.Control.SubjectID = string(tbl_demo.SubjectID(strcmpi(tbl_demo.group,'Control')));
RSFC.Control.ROI_Z_pre  = zPre_ctrl_all;
RSFC.Control.ROI_Z_post = zPost_ctrl_all;
RSFC.Control.ROI_Z_diff = zPost_ctrl_all - zPre_ctrl_all;
RSFC.Control.Lobe_Z_pre  = Zlob_pre_ctrl;
RSFC.Control.Lobe_Z_post = Zlob_post_ctrl;
RSFC.Control.Lobe_Z_diff = Zlob_diff_ctrl;

% ---- provenance ----
RSFC.Info.ROI_seed_channels = roiIdx_hbo_exp;
RSFC.Info.RestWindow_s      = durSec;
RSFC.Info.Units             = 'Fisher Z unless a field name says r';
RSFC.Info.Created           = datetime('now');

save(fullfile(DATA_DIR, 'RSFC_outputs.mat'), 'RSFC');

fprintf('\nRSFC outputs saved to %s\n', fullfile(DATA_DIR, 'RSFC_outputs.mat'));
fprintf('  tinnitus n = %d | control n = %d | lobes = %d\n', ...
    numel(RSFC.SubjectID), numel(RSFC.Control.SubjectID), numel(RSFC.Lobes));


%% ALL LOCAL FUNCTIONS TO RUN RSFC CODE 


function Zlob = computeLobeMeans(Zsubj_by_ch, chanLobe_hbo, lobes, validTarget)
    nSub = size(Zsubj_by_ch,1);
    nL   = numel(lobes);
    Zlob = NaN(nSub, nL);

    for L = 1:nL
        lb = lobes{L};
        idx = strcmpi(chanLobe_hbo, lb) & validTarget;

        if ~any(idx)
            continue;
        end

        Zlob(:,L) = nanmean(Zsubj_by_ch(:, idx), 2);
    end
end

function S = buildSummary(Tok, edgeZ, ~)
    if numel(edgeZ) ~= height(Tok)
        error('edgeZ length (%d) does not match Tok height (%d).', numel(edgeZ), height(Tok));
    end

    G = findgroups(Tok.lobeA2, Tok.lobeB2);

    nEdge   = splitapply(@numel, edgeZ, G);
    meanZ   = splitapply(@(x) mean(x,'omitnan'), edgeZ, G);
    meanR   = tanh(meanZ);

    pctBil  = 100 * splitapply(@(x) mean(x,'omitnan'), double(Tok.isBilateral), G);

    lobeA_u = splitapply(@(x) x(1), Tok.lobeA2, G);
    lobeB_u = splitapply(@(x) x(1), Tok.lobeB2, G);

    S = table(string(lobeA_u), string(lobeB_u), nEdge, meanZ, meanR, pctBil, ...
        'VariableNames', {'lobeA','lobeB','n','meanZ','meanR','pctBilateral'});

    S.pct = 100 * S.n / sum(S.n);
    S = sortrows(S, 'n', 'descend');
end

function rgb = valToRGB(v, cmap, clim)
    if ~isfinite(v), rgb = [1 1 1]; return; end
    v = max(min(v, clim(2)), clim(1));
    t = (v - clim(1)) / (clim(2) - clim(1));
    idx = 1 + round(t * (size(cmap,1)-1));
    idx = max(1, min(size(cmap,1), idx));
    rgb = cmap(idx, :);
end

function drawCircleFilled(ax, cx, cy, rad, faceRGB)
    rectangle(ax,'Position',[cx-rad, cy-rad, 2*rad, 2*rad], ...
        'Curvature',[1 1], ...
        'FaceColor',faceRGB, ...
        'EdgeColor','k', ...
        'LineWidth',0.6);
end