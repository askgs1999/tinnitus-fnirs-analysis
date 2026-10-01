%%  EVENT-LOCKED HEMODYNAMIC RESPONSE ANALYSIS
%  Stimulus-evoked HbO responses to broadband noise and interstimulus rest
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
%    3. Defines BBN and ISR events from the recorded triggers and extracts
%       block-averaged, baseline-corrected HbO epochs.
%    4. Validates the auditory ROI functionally in the control group, then
%       derives the early and late analysis windows from the control
%       group-mean response.
%    5. Computes group-mean HbO curves for the auditory ROI, the non-ROI
%       spatial control channels, and the frontal and occipital subsets.
%    6. Runs the within-group and between-group statistics reported in the
%       manuscript for the early (5-12 s) and late (13-20 s) windows.
%    7. Generates Figures 4, 5 and 6.
%    8. Writes subject-wise window means to HR_outputs.mat for the
%       questionnaire analysis.
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
%    HR_outputs.mat   struct HR, written to DATA_DIR
%    Fig4.tif Fig5.tif Fig6.tif   written to FIG_DIR
%
%  HOW TO RUN
%    Set MATLAB's current folder to the one holding
%    Aim1_raw_loaded_new_03_17_new.mat, then run this script.
%



clearvars; close all; clc;

%% ---- PATHS ----
% The input file is located automatically using WHICH, which searches the
% current folder and the MATLAB path, so this works whether or not the
% script has been saved. To set the folder by hand, assign DATA_DIR
% directly and delete the block below.

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


%% loadAim1_raw_loaded.mat directly   
load(fullfile(DATA_DIR, 'Aim1_raw_loaded_new_03_17_new.mat'));   % loads raw + tbl_demo
%Aim1_raw_loaded_new_03_17_new.mat - FINAL

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

%% Step 1: SCI 


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

    % identify long and short channels
    isShort = link.ShortSeperation == 1;
    isLong  = link.ShortSeperation == 0;

    % skip subjects without short channels
    if ~any(isShort)
        warning('Subject %s has no short channels — skipping SSR', ...
            TDDR_clean(subj).demographics('SubjectID'));
        continue;
    end

    % Channel midpoints
    sPos = TDDR_clean(subj).probe.srcPos(link.source,:);
    dPos = TDDR_clean(subj).probe.detPos(link.detector,:);
    mid  = (sPos + dPos)/2;

    shortIdx = find(isShort);
    longIdx  = find(isLong);

    for ii = 1:numel(longIdx)
        L = longIdx(ii);

        % find closest short channel
        d = vecnorm(mid(shortIdx,:) - mid(L,:), 2, 2);
        [~, k] = min(d);
        S = shortIdx(k);

        y = TDDR_clean(subj).data(:,L);   % long OD
        x = TDDR_clean(subj).data(:,S);   % short OD

        % mean-center
        y = y - mean(y);
        x = x - mean(x);

        % regression coefficient
        beta = x \ y;

        % subtract scaled short channel
        ODdata_ssr(subj).data(:,L) = y - beta*x;
    end

    fprintf('Subject %s: SSR applied to %d long channels\n', ...
        TDDR_clean(subj).demographics('SubjectID'), numel(longIdx));
end

%% Remove short channels 
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
%% Step 8: Map channels to anatomical ROIs (using Colin27)
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
disp(HaemoData(1).probe.link.Properties.VariableNames)


%% === FIX: unify ROI_full columns ===
for i = 1:length(HaemoData)
    link = HaemoData(i).probe.link;

    % Check if columns like ROI_full_left / ROI_full_right exist
    if ismember('ROI_full_left', link.Properties.VariableNames) && ...
       ismember('ROI_full_right', link.Properties.VariableNames)

        % Create a clean ROI_full column from left/right
        link.ROI_full = link.ROI_full_left;  % start with left hemisphere

        % If left side is empty, replace with right
        emptyIdx = cellfun(@isempty, link.ROI_full);
        link.ROI_full(emptyIdx) = link.ROI_full_right(emptyIdx);

        % --- Remove duplicate ROI columns safely (backward-compatible) ---
        colsToRemove = {'ROI_full_left','ROI_full_right',...
                        'ROI_code_left','ROI_code_right',...
                        'ROI_full_left_1','ROI_full_right_1',...
                        'ROI_code_left_1','ROI_code_right_1'};

        colsToRemove = intersect(colsToRemove, link.Properties.VariableNames);
        if ~isempty(colsToRemove)
            link = removevars(link, colsToRemove);
        end

        % Store back
        HaemoData(i).probe.link = link;
    end
end

disp("✅ ROI_full unified successfully. Sample columns:")
disp(HaemoData(1).probe.link.Properties.VariableNames)
disp(unique(HaemoData(1).probe.link.ROI_full))


%% Step 9 : Split HaemoData into HbO and HbR datasets safely

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

%% Step 10: Chopping fragments of data to BBN and ISR 
for i = 1:length(HaemoData_hbo)
    dat = HaemoData_hbo(i);

    % --- Extract correct trigger ("2") as BBN
    stim_bbn = dat.stimulus('2');   % not '1'
    stim_bbn.name = 'BBN';
    stim_bbn.dur  = 12 * ones(size(stim_bbn.onset));  % enforce 12s

    % --- Define ISR 
    stim_isr = nirs.design.StimulusEvents();
    stim_isr.name  = 'ISR';
    stim_isr.onset = stim_bbn.onset + 12;
    stim_isr.dur   = 17 * ones(size(stim_bbn.onset));
    stim_isr.amp   = ones(size(stim_bbn.onset));

    % --- Save back
    dat.stimulus('BBN') = stim_bbn;
    dat.stimulus('ISR') = stim_isr;
    HaemoData_hbo(i) = dat;
end

dat = HaemoData_hbo(1);


control_idx   = strcmpi(tbl_demo.group,'Control');
tinnitus_idx  = strcmpi(tbl_demo.group,'Experimental');

ControlData   = HaemoData_hbo(control_idx);
TinnitusData  = HaemoData_hbo(tinnitus_idx);

%% --- CHECK BBN & ISR ONSETS AND DURATIONS ---
for i = 1:length(HaemoData_hbo)
    subjID = HaemoData_hbo(i).demographics('SubjectID');
    dat = HaemoData_hbo(i);

    % --- Get list of available stimuli
    stim_names = dat.stimulus.keys;

    % --- Proceed only if both exist
    if ismember('BBN', stim_names) && ismember('ISR', stim_names)
        stim_bbn = dat.stimulus('BBN');
        stim_isr = dat.stimulus('ISR');
        
        % --- Print summary
        fprintf('\n===== Subject %s =====\n', subjID);
        fprintf('Total BBN blocks: %d\n', length(stim_bbn.onset));
        fprintf('Total ISR blocks: %d\n', length(stim_isr.onset));
        fprintf('First BBN onset: %.2f s | Duration: %.2f s\n', stim_bbn.onset(1), stim_bbn.dur(1));
        fprintf('First ISR onset: %.2f s | Duration: %.2f s\n', stim_isr.onset(1), stim_isr.dur(1));
        fprintf('Last BBN onset:  %.2f s\n', stim_bbn.onset(end));
        fprintf('Last ISR onset:  %.2f s\n', stim_isr.onset(end));

        % % --- Visualize timeline
        % figure('Name', ['Stimulus timeline - ' subjID]);
        % hold on;
        % % Plot BBN (red bars)
        % for j = 1:length(stim_bbn.onset)
        %     rectangle('Position', [stim_bbn.onset(j), 0, stim_bbn.dur(j), 1], ...
        %               'FaceColor', [1 0 0 0.4], 'EdgeColor', 'none');
        % end
        % % Plot ISR (blue bars)
        % for j = 1:length(stim_isr.onset)
        %     rectangle('Position', [stim_isr.onset(j), 0, stim_isr.dur(j), 1], ...
        %               'FaceColor', [0 0 1 0.4], 'EdgeColor', 'none');
        % end
        % 
        % % Format
    
        % ylim([0, 1]);
        % xlabel('Time (s)');
        % ylabel('Stimulus block');
        % title(sprintf('BBN (red) vs ISR (blue) – %s', subjID));
        % legend({'BBN','ISR'});
        % grid on;

        % --- Optional: check gap consistency
        gap = stim_isr.onset - (stim_bbn.onset + stim_bbn.dur);
        fprintf('Mean gap between BBN→ISR: %.2f s (range %.2f–%.2f)\n', ...
                mean(gap), min(gap), max(gap));

    else
        fprintf('\n⚠️ Subject %s missing BBN/ISR events\n', subjID);
        fprintf('Available stimuli: %s\n', strjoin(stim_names, ', '));
    end
end


%% Step 11: DEFINE ROI AND NON ROI LABELS
roi_labels = {
    'Temporal_Sup_L', 'Temporal_Mid_L', 'Temporal_Inf_L', ...
    'Temporal_Sup_R', 'Temporal_Mid_R', 'Temporal_Inf_R', ...
    'Supramarginal_L', 'Supramarginal_R', ...
    'Postcentral_L', 'Postcentral_R', ...
    'Angular_L', 'Angular_R'
};

nonroi_labels = {
    'Occipital_Sup_L', 'Occipital_Mid_L', 'Occipital_Sup_R', 'Occipital_Mid_R', ...
    'Cuneus_L', 'Cuneus_R', 'Lingual_R', 'Calcarine_L', ...
    'Frontal_Mid_L', 'Frontal_Mid_R', ...
    'Frontal_Sup_L', 'Frontal_Sup_R', ...
    'Frontal_Inf_Tri_R', 'Frontal_Inf_Oper_L'
};

%% Step 12: Control Group:  FUNCTIONAL VALIDATION OF ROI 
fs               = ControlData(1).Fs;
prestim          = 2;
poststim_val     = 12;
prestim_samp     = round(prestim * fs);
poststim_samp    = round(poststim_val * fs);
sampOffsets_val  = (-prestim_samp : poststim_samp)';
timevec_len      = prestim_samp + poststim_samp + 1;
roi_win          = (prestim_samp + 1) : timevec_len;

roi_final = struct();

for s = 1:length(ControlData)
    subj = ControlData(s);
    link = subj.probe.link;

    is_hbo  = strcmpi(link.type, 'hbo');
    roi_idx = find(is_hbo & ismember(link.ROI_full, roi_labels));

    if isempty(roi_idx)
        fprintf('⚠️ No anatomical ROI channels for subject %s\n', subj.demographics('SubjectID'));
        roi_final(s).SubjectID = subj.demographics('SubjectID');
        roi_final(s).channels  = [];
        continue;
    end

    stim_bbn = subj.stimulus('BBN');
    stim_isr = subj.stimulus('ISR');

    if isempty(stim_bbn.onset) || isempty(stim_isr.onset)
        fprintf('⚠️ Missing stimulus for subject %s\n', subj.demographics('SubjectID'));
        roi_final(s).SubjectID = subj.demographics('SubjectID');
        roi_final(s).channels  = [];
        continue;
    end

    bbn_onsets = stim_bbn.onset(:);
    isr_onsets = stim_isr.onset(:);

    meanBBN = zeros(length(roi_idx), 1);
    meanISR = zeros(length(roi_idx), 1);

    for ch_i = 1:length(roi_idx)
        ch         = roi_idx(ch_i);
        bbn_trials = [];
        isr_trials = [];

        for b = 1:min(length(bbn_onsets), length(isr_onsets))
            idx_bbn = round(bbn_onsets(b) * fs) + sampOffsets_val;
            idx_isr = round(isr_onsets(b) * fs) + sampOffsets_val;

            if idx_bbn(1) < 1 || idx_bbn(end) > size(subj.data,1), continue; end
            if idx_isr(1) < 1 || idx_isr(end) > size(subj.data,1), continue; end

            trial_bbn = subj.data(idx_bbn, ch);
            trial_bbn = trial_bbn - mean(trial_bbn(1:prestim_samp));
            bbn_trials(:, end+1) = trial_bbn; %#ok<SAGROW>

            trial_isr = subj.data(idx_isr, ch);
            trial_isr = trial_isr - mean(trial_isr(1:prestim_samp));
            isr_trials(:, end+1) = trial_isr; %#ok<SAGROW>
        end

        if isempty(bbn_trials) || isempty(isr_trials)
            continue;
        end

        meanBBN(ch_i) = mean(mean(bbn_trials(roi_win, :), 1), 'omitnan');
        meanISR(ch_i) = mean(mean(isr_trials(roi_win, :), 1), 'omitnan');
    end

    is_valid      = meanBBN > meanISR;
    roi_final_idx = roi_idx(is_valid);

    roi_final(s).SubjectID = subj.demographics('SubjectID');
    roi_final(s).channels  = roi_final_idx;

    fprintf('✅ %s → %d validated ROI channels\n', ...
        roi_final(s).SubjectID, length(roi_final_idx));
end

% === NON-ROI CHANNELS (CONTROLS) ===
nonroi_final = struct();

for s = 1:length(ControlData)
    subj = ControlData(s);
    link = subj.probe.link;

    is_hbo     = strcmpi(link.type, 'hbo');
    nonroi_idx = find(is_hbo & ismember(link.ROI_full, nonroi_labels));

    if isempty(nonroi_idx)
        fprintf('⚠️ No anatomical non-ROI channels for subject %s\n', subj.demographics('SubjectID'));
        nonroi_final(s).SubjectID = subj.demographics('SubjectID');
        nonroi_final(s).channels  = [];
        continue;
    end

    nonroi_final(s).SubjectID = subj.demographics('SubjectID');
    nonroi_final(s).channels  = nonroi_idx;

    fprintf('✅ %s → %d anatomical non-ROI channels\n', ...
        nonroi_final(s).SubjectID, length(nonroi_idx));
end




%% FUNCTIONAL VALIDATION STEPS 

% 1. For every BBN block, take from -2 s before BBN onset (baseline
% correction)
%2. Define duration of BBN and ISR 
%3. Define search windows 
%4. For each participant u build a hbo response cure, extract epochs around every BBN onset, for each ROI channel
% 5. Average trials → channel response; Average channels → subject ROI response
%6. average the per sub hbo response curve to form one mean hbo response
%curve 
% 7. take the peak +ve hbo amplitude and what time its happening ; take the
% peak -ve amplitude of hbo curve and what time its happening..thats the
% mena time of max and min 
% 8 ; now define window as [ mean - 1 sd, mean, mean + 1 sd]

% Final perfectly working version 
% CONTROLS ROI HbO: PEAK (MAX) + ISR LOWEST (MIN) WINDOWS
% Centers: from GROUP MEAN curve
% Widths:  ± 1 SD of SUBJECT-wise peak/trough times
% Plot:
%  - Shade BBN (0–12s) in red, ISR (12–27s) in blue
%  - Red dotted lines = peak window boundaries
%  - Blue dotted lines = trough window boundaries


prestim  = 2;
poststim = 30;

BBN_dur = 12;
ISR_dur = 17;              % enforced in your stim creation
ISR_end = BBN_dur + ISR_dur;

% Search windows (broad, only to avoid nonsense)
peakSearch   = [3 16];     % peak should be here
troughSearch = [12 25];    % trough should be here (ISR/recovery)

fs = ControlData(1).Fs;
tEpoch = (-prestim : 1/fs : poststim)';     % time axis
sampOffsets = round(tEpoch * fs);

iBase   = find(tEpoch >= -prestim & tEpoch < 0);
iPeakS  = find(tEpoch >= peakSearch(1) & tEpoch <= peakSearch(2));
iTrS    = find(tEpoch >= troughSearch(1) & tEpoch <= troughSearch(2));

% per-subject outputs (for SD in time)
subj_peak_t   = nan(numel(ControlData),1);
subj_trough_t = nan(numel(ControlData),1);

% store subject ROI timecourses for group plot (time x subjects)
subj_tc = [];
subj_used = strings(0);

for s = 1:numel(ControlData)
    subj = ControlData(s);

    % check BBN exists
    stim_names = subj.stimulus.keys;
    if ~ismember('BBN', stim_names)
        fprintf('⚠️ %s: no BBN key\n', subj.demographics('SubjectID'));
        continue;
    end

    % check roi_final exists for subject
    if s > numel(roi_final) || ~isfield(roi_final(s),'channels') || isempty(roi_final(s).channels)
        fprintf('⚠️ %s: roi_final empty (no validated ROI channels)\n', subj.demographics('SubjectID'));
        continue;
    end

    onsets = subj.stimulus('BBN').onset(:);
    if isempty(onsets)
        fprintf('⚠️ %s: BBN onsets empty\n', subj.demographics('SubjectID'));
        continue;
    end

    roi_ch = roi_final(s).channels(:)';  % validated ROI channels

    % ---- subject ROI block-average curve ----
    roi_ch_tc = nan(numel(tEpoch), numel(roi_ch)); % time x ROIch

    for ii = 1:numel(roi_ch)
        ch = roi_ch(ii);
        trials = nan(numel(tEpoch), numel(onsets)); % time x blocks

        for b = 1:numel(onsets)
            idx0 = round(onsets(b) * fs);
            idx  = idx0 + sampOffsets;

            if idx(1) < 1 || idx(end) > size(subj.data,1)
                continue;
            end

            x = subj.data(idx, ch);
            x = x - mean(x(iBase), 'omitnan');   % baseline (-2..0)
            trials(:,b) = x;
        end

        roi_ch_tc(:,ii) = mean(trials, 2, 'omitnan'); % avg across blocks
    end

    tc = mean(roi_ch_tc, 2, 'omitnan'); % avg across ROI channels
    if all(isnan(tc))
        fprintf('⚠️ %s: tc all NaN\n', subj.demographics('SubjectID'));
        continue;
    end

    subj_tc(:, end+1) = tc; %#ok<SAGROW>
    subj_used(end+1)  = string(subj.demographics('SubjectID')); %#ok<SAGROW>

    % ---- SUBJECT peak time (MAX HbO) ----
    [~, irelP] = max(tc(iPeakS));
    subj_peak_t(s) = tEpoch(iPeakS(irelP));

    % ---- SUBJECT trough time (MIN HbO during ISR) ----
    [~, irelT] = min(tc(iTrS));
    subj_trough_t(s) = tEpoch(iTrS(irelT));

    fprintf('✅ %s | subj peak @ %.2fs | subj trough @ %.2fs\n', ...
        subj.demographics('SubjectID'), subj_peak_t(s), subj_trough_t(s));
end

% --- VALID subjects for SD calculation
validP = ~isnan(subj_peak_t);
validT = ~isnan(subj_trough_t);

if isempty(subj_tc)
    error('No valid subjects made it into subj_tc. Check roi_final and stimulus keys.');
end

% --- GROUP MEAN curve (center comes from this)
muCurve = mean(subj_tc, 2, 'omitnan');

% ======== PEAK center from GROUP MEAN curve ========
[~, irelP_mu] = max(muCurve(iPeakS));
t_peak = tEpoch(iPeakS(irelP_mu));       % center time of peak from group mean

sdPeak = std(subj_peak_t(validP));       % width from subject variability
muPeak = t_peak;
peakWin = [muPeak - sdPeak, muPeak + sdPeak];

% ======== TROUGH center from GROUP MEAN curve ========
[~, irelT_mu] = min(muCurve(iTrS));
t_lowest = tEpoch(iTrS(irelT_mu));       % center time of trough from group mean

sdTr = std(subj_trough_t(validT));       % width from subject variability
muTr = t_lowest;
troughWin = [muTr - sdTr, muTr + sdTr];

fprintf('\n=== WINDOWS (center from GROUP MEAN, width = SD of subject times) ===\n');
fprintf('Peak time (group-mean):   %.2f ± %.2f s  => [%.2f, %.2f]\n', ...
    muPeak, sdPeak, peakWin(1), peakWin(2));
fprintf('Trough time (group-mean): %.2f ± %.2f s  => [%.2f, %.2f]\n', ...
    muTr, sdTr, troughWin(1), troughWin(2));

% --- Plot group mean ± SEM with shading + dotted boundaries
sem = std(subj_tc, 0, 2, 'omitnan') ./ sqrt(size(subj_tc,2));

figure; hold on;

% Plot SEM band + mean
fill([tEpoch; flipud(tEpoch)], [muCurve-sem; flipud(muCurve+sem)], 'k', ...
    'FaceAlpha', 0.12, 'EdgeColor', 'none');
plot(tEpoch, muCurve, 'k', 'LineWidth', 2);
yline(0,'--');

yl = ylim;

% Shading: BBN red, ISR blue
patch([0 BBN_dur BBN_dur 0], [yl(1) yl(1) yl(2) yl(2)], [1 0 0], ...
    'FaceAlpha', 0.08, 'EdgeColor', 'none');

patch([BBN_dur min(ISR_end,poststim) min(ISR_end,poststim) BBN_dur], ...
      [yl(1) yl(1) yl(2) yl(2)], [0 0 1], ...
      'FaceAlpha', 0.06, 'EdgeColor', 'none');

% Re-draw on top of shading
fill([tEpoch; flipud(tEpoch)], [muCurve-sem; flipud(muCurve+sem)], 'k', ...
    'FaceAlpha', 0.12, 'EdgeColor', 'none');
plot(tEpoch, muCurve, 'k', 'LineWidth', 2);
yline(0,'--');

% Stimulus boundaries
xline(0,':');
xline(BBN_dur,':');

% Peak window dotted RED lines (centered on group mean peak)
xline(peakWin(1), ':', 'Color', [1 0 0], 'LineWidth', 2);
xline(peakWin(2), ':', 'Color', [1 0 0], 'LineWidth', 2);

% Trough window dotted BLUE lines (centered on group mean trough)
xline(troughWin(1), ':', 'Color', [0 0 1], 'LineWidth', 2);
xline(troughWin(2), ':', 'Color', [0 0 1], 'LineWidth', 2);



xlabel('Time (s) from BBN onset');
ylabel('\Delta HbO (\muM, baseline-corrected)');
title(sprintf('Controls ROI HbO (mean \\pm SEM), n=%d', size(subj_tc,2)));
grid on;

legend({'SEM','Mean','Baseline', ...
        'BBN (shaded)','ISR (shaded)', ...
        'BBN onset','BBN offset', ...
        'Peak win start','Peak win end', ...
        'Trough win start','Trough win end'}, ...
        'Location','best');

%% generate_Fig4_PLOS.m

% % PREREQUISITE: Run your main analysis script first so the workspace contains:
% %   tEpoch       (time vector)
% %   subj_tc      (time x subjects matrix of subject ROI timecourses)
% %   muCurve      (group-mean curve)
% %   peakWin      (1x2 vector: [peak_lo peak_hi])
% %   troughWin    (1x2 vector: [trough_lo trough_hi])
% %   BBN_dur      (BBN duration, e.g. 12)
% %   ISR_end      (ISR end time, e.g. 27)
% %   poststim     (post-stimulus axis end, e.g. 30)


outDir = FIG_DIR;
if ~exist(outDir, 'dir'); mkdir(outDir); end
outFile = fullfile(outDir, 'Fig4.tif');

% ----- Figure size: full page width, single panel -----
% PLOS max width = 19.05 cm (7.5 in). For a single-panel time series,
% a wide-but-short aspect works well (e.g., 7.5 x 4 inches).
figWidthIn  = 7.5;
figHeightIn = 4.0;
dpi         = 300;

% ----- Font settings -----
axFont    = 'Arial';
axSize    = 10;
labSize   = 11;
legSize   = 9;

% ----- Recompute SEM in case it isn't in workspace -----
sem = std(subj_tc, 0, 2, 'omitnan') ./ sqrt(size(subj_tc,2));
nSubj = size(subj_tc, 2);


% Build figure

fig = figure('Units','inches', ...
             'Position',[1 1 figWidthIn figHeightIn], ...
             'Color','w', ...
             'PaperUnits','inches', ...
             'PaperPosition',[0 0 figWidthIn figHeightIn], ...
             'PaperSize',[figWidthIn figHeightIn]);

ax = axes('Parent', fig); hold(ax,'on');

% --- Determine y-limits from data so shading patches fill correctly ---
yMin = min(muCurve - sem) * 1.15;
yMax = max(muCurve + sem) * 1.30;   % give headroom for legend
ylim(ax, [yMin yMax]);

% --- Shaded stimulus periods (drawn FIRST, so they sit behind the data) ---
hBBN = patch(ax, [0 BBN_dur BBN_dur 0], [yMin yMin yMax yMax], ...
    [1 0.85 0.85], 'FaceAlpha', 0.5, 'EdgeColor','none', ...
    'DisplayName','BBN stimulus');

isrEnd = min(ISR_end, poststim);
hISR = patch(ax, [BBN_dur isrEnd isrEnd BBN_dur], [yMin yMin yMax yMax], ...
    [0.85 0.85 1], 'FaceAlpha', 0.5, 'EdgeColor','none', ...
    'DisplayName','ISR period');

% --- SEM band ---
hSEM = fill(ax, [tEpoch; flipud(tEpoch)], [muCurve-sem; flipud(muCurve+sem)], ...
    [0.5 0.5 0.5], 'FaceAlpha', 0.30, 'EdgeColor','none', ...
    'DisplayName','SEM');

% --- Group mean curve ---
hMean = plot(ax, tEpoch, muCurve, 'k-', 'LineWidth', 2, ...
    'DisplayName','Group mean');

% --- Zero baseline ---
plot(ax, [0 poststim], [0 0], 'k--', 'LineWidth', 0.5, ...
    'HandleVisibility','off');

% --- Peak window (red dotted verticals) ---
hPeak = plot(ax, [peakWin(1) peakWin(1)], [yMin yMax], ':', ...
    'Color',[0.85 0 0], 'LineWidth', 1.8, ...
    'DisplayName','Peak window');
plot(ax, [peakWin(2) peakWin(2)], [yMin yMax], ':', ...
    'Color',[0.85 0 0], 'LineWidth', 1.8, 'HandleVisibility','off');

% --- Trough window (blue dotted verticals) ---
hTr = plot(ax, [troughWin(1) troughWin(1)], [yMin yMax], ':', ...
    'Color',[0 0 0.85], 'LineWidth', 1.8, ...
    'DisplayName','Trough window');
plot(ax, [troughWin(2) troughWin(2)], [yMin yMax], ':', ...
    'Color',[0 0 0.85], 'LineWidth', 1.8, 'HandleVisibility','off');

% --- Axes formatting ---

xlabel(ax, 'Time (s) from BBN onset', 'FontName',axFont, 'FontSize',labSize);
ylabel(ax, '\Delta HbO (\muM, baseline-corrected)', 'FontName',axFont, 'FontSize',labSize);

set(ax, 'FontName', axFont, 'FontSize', axSize, ...
        'LineWidth', 1, 'TickDir','out', 'Box','on');

% --- Legend (only meaningful entries, in logical order) ---
lgd = legend(ax, [hMean hSEM hBBN hISR hPeak hTr], ...
    'Location','northeast', ...
    'FontName',axFont, 'FontSize',legSize, ...
    'Box','on');
% Make legend background semi-opaque so it doesn't obscure data
set(lgd, 'Color', [1 1 1 0.9]);

% Annotate N inside the plot (top-left), so it's reproducible without title
% text(ax, 0.02, 0.97, sprintf('n = %d', nSubj), ...
%     'Units','normalized', ...
%     'FontName',axFont, 'FontSize',legSize, ...
%     'VerticalAlignment','top', ...
%     'BackgroundColor',[1 1 1 0.7]);

% =========================================================================
% Export as PLOS-compliant TIFF
% =========================================================================
% try
%     exportgraphics(fig, outFile, 'Resolution', dpi, 'BackgroundColor','white');
%     fprintf('Saved (exportgraphics): %s\n', outFile);
% catch
%     print(fig, outFile, '-dtiff', sprintf('-r%d', dpi));
%     fprintf('Saved (print): %s\n', outFile);
% end

% =========================================================================
% Compliance check
% =========================================================================
% info = imfinfo(outFile);
% fprintf('\n--- Fig4.tif specs ---\n');
% fprintf('  Width:     %d px\n', info.Width);
% fprintf('  Height:    %d px\n', info.Height);
% fprintf('  XRes:      %g dpi\n', info.XResolution);
% fprintf('  YRes:      %g dpi\n', info.YResolution);
% fprintf('  ColorType: %s, BitDepth: %d\n', info.ColorType, info.BitDepth);
% fprintf('  Compression: %s\n', info.Compression);
% fprintf('  File size: %.2f MB\n', info.FileSize / 1024 / 1024);
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
% 

%% Visualizations

% CONTROL ROI
fs = ControlData(1).Fs;
prestim  = 2;
poststim = 20;

tEpoch = (-prestim:1/fs:poststim)';
sampOffsets = round(tEpoch * fs);

iBase = find(tEpoch >= -1 & tEpoch <= 0);
iPlot = find(tEpoch >= 0 & tEpoch <= poststim);

% Colors (match figure)
cBBN = [1 0 0];      % red
cISR = [0 1 0];      % green

ROI_BBN_ctrl = [];
ROI_ISR_ctrl = [];

for s = 1:numel(ControlData)
    subj = ControlData(s);

    if s > numel(roi_final) || isempty(roi_final(s).channels)
        continue;
    end
    roi_ch = roi_final(s).channels(:)';

    onBBN = subj.stimulus('BBN').onset(:);
    onISR = subj.stimulus('ISR').onset(:);

    % ----- BBN -----
    tc = nan(numel(tEpoch), numel(onBBN));
    for b = 1:numel(onBBN)
        idx = round(onBBN(b)*fs) + sampOffsets;
        if idx(1)<1 || idx(end)>size(subj.data,1), continue; end
        x = mean(subj.data(idx,roi_ch),2);
        x = x - mean(x(iBase));
        tc(:,b) = x;
    end
    ROI_BBN_ctrl(:,end+1) = mean(tc,2,'omitnan'); %#ok<SAGROW>

    % ----- ISR -----
    tc = nan(numel(tEpoch), numel(onISR));
    for b = 1:numel(onISR)
        idx = round(onISR(b)*fs) + sampOffsets;
        if idx(1)<1 || idx(end)>size(subj.data,1), continue; end
        x = mean(subj.data(idx,roi_ch),2);
        x = x - mean(x(iBase));
        tc(:,b) = x;
    end
    ROI_ISR_ctrl(:,end+1) = mean(tc,2,'omitnan'); %#ok<SAGROW>
end

muBBN = mean(ROI_BBN_ctrl,2,'omitnan');
muISR = mean(ROI_ISR_ctrl,2,'omitnan');

seBBN = std(ROI_BBN_ctrl,0,2,'omitnan') / sqrt(size(ROI_BBN_ctrl,2));
seISR = std(ROI_ISR_ctrl,0,2,'omitnan') / sqrt(size(ROI_ISR_ctrl,2));

figure; hold on;

fill([tEpoch(iPlot); flipud(tEpoch(iPlot))], ...
     [muBBN(iPlot)-seBBN(iPlot); flipud(muBBN(iPlot)+seBBN(iPlot))], ...
     cBBN, 'FaceAlpha',0.15,'EdgeColor','none');
plot(tEpoch(iPlot), muBBN(iPlot), 'Color', cBBN, 'LineWidth',2);

fill([tEpoch(iPlot); flipud(tEpoch(iPlot))], ...
     [muISR(iPlot)-seISR(iPlot); flipud(muISR(iPlot)+seISR(iPlot))], ...
     cISR, 'FaceAlpha',0.15,'EdgeColor','none');
plot(tEpoch(iPlot), muISR(iPlot), 'Color', cISR, 'LineWidth',2);

yline(0,'--','Color',[0.5 0.5 0.5]);
xlabel('Time (s)');
ylabel('\Delta HbO (\muM)');
title('Control Group — Mean HbO (ROI Channels)');
legend({'BBN','BBN mean','ISR','ISR mean'},'Location','southwest');
grid on;

% CONTROL NON ROI (CTRL vars)

NROI_BBN_ctrl = [];
NROI_ISR_ctrl = [];

for s = 1:numel(ControlData)
    subj = ControlData(s);

    if s > numel(nonroi_final) || isempty(nonroi_final(s).channels)
        continue;
    end
    nonroi_ch = nonroi_final(s).channels(:)';

    onBBN = subj.stimulus('BBN').onset(:);
    onISR = subj.stimulus('ISR').onset(:);

    % ----- BBN -----
    tc = nan(numel(tEpoch), numel(onBBN));
    for b = 1:numel(onBBN)
        idx = round(onBBN(b)*fs) + sampOffsets;
        if idx(1)<1 || idx(end)>size(subj.data,1), continue; end
        x = mean(subj.data(idx,nonroi_ch),2);
        x = x - mean(x(iBase));
        tc(:,b) = x;
    end
    NROI_BBN_ctrl(:,end+1) = mean(tc,2,'omitnan'); %#ok<SAGROW>

    % ----- ISR -----
    tc = nan(numel(tEpoch), numel(onISR));
    for b = 1:numel(onISR)
        idx = round(onISR(b)*fs) + sampOffsets;
        if idx(1)<1 || idx(end)>size(subj.data,1), continue; end
        x = mean(subj.data(idx,nonroi_ch),2);
        x = x - mean(x(iBase));
        tc(:,b) = x;
    end
    NROI_ISR_ctrl(:,end+1) = mean(tc,2,'omitnan'); %#ok<SAGROW>
end

muBBN = mean(NROI_BBN_ctrl,2,'omitnan');
muISR = mean(NROI_ISR_ctrl,2,'omitnan');

seBBN = std(NROI_BBN_ctrl,0,2,'omitnan') / sqrt(size(NROI_BBN_ctrl,2));
seISR = std(NROI_ISR_ctrl,0,2,'omitnan') / sqrt(size(NROI_ISR_ctrl,2));

figure; hold on;

fill([tEpoch(iPlot); flipud(tEpoch(iPlot))], ...
     [muBBN(iPlot)-seBBN(iPlot); flipud(muBBN(iPlot)+seBBN(iPlot))], ...
     cBBN, 'FaceAlpha',0.15,'EdgeColor','none');
plot(tEpoch(iPlot), muBBN(iPlot), 'Color', cBBN, 'LineWidth',2);

fill([tEpoch(iPlot); flipud(tEpoch(iPlot))], ...
     [muISR(iPlot)-seISR(iPlot); flipud(muISR(iPlot)+seISR(iPlot))], ...
     cISR, 'FaceAlpha',0.15,'EdgeColor','none');
plot(tEpoch(iPlot), muISR(iPlot), 'Color', cISR, 'LineWidth',2);

yline(0,'--','Color',[0.5 0.5 0.5]);
xlabel('Time (s)');
ylabel('\Delta HbO (\muM)');
title('Control Group — Mean HbO (Non-ROI Channels)');
legend({'BBN','BBN mean','ISR','ISR mean'},'Location','southwest');
grid on;
%% OPTIONAL: freeze into paper-stats variable names I gave you earlier
CTRL_ROI_BBN = ROI_BBN_ctrl;
CTRL_ROI_ISR = ROI_ISR_ctrl;
CTRL_NROI_BBN = NROI_BBN_ctrl;
CTRL_NROI_ISR = NROI_ISR_ctrl;


%% TINNITUS — Group HbO curve (BBN-locked) using fixed channel set


clc;
roi_fixed = [38 41 21 8 24 2];

% Diagnostic — what these channels map to

link     = TinnitusData(1).probe.link;
is_hbo   = strcmpi(link.type, 'hbo');
hbo_link = link(is_hbo, :);

fprintf('\n=== roi_fixed channel mapping ===\n');
disp(hbo_link(roi_fixed, {'source','detector','type','ROI_full'}));


% Epoch params

prestim  = 2;
poststim = 30;
BBN_dur  = 12;
ISR_dur  = 17;
ISR_end  = BBN_dur + ISR_dur;

troughSearch = [12 25];

fs          = TinnitusData(1).Fs;
tEpoch      = (-prestim : 1/fs : poststim)';
sampOffsets = round(tEpoch * fs);

iBase = find(tEpoch >= -prestim & tEpoch < 0);
iTrS  = find(tEpoch >= troughSearch(1) & tEpoch <= troughSearch(2));


% Build subject-level ROI block-averaged curves

subj_trough_t = nan(numel(TinnitusData), 1);
subj_tc       = [];
subj_used     = strings(0);

for s = 1:numel(TinnitusData)
    subj = TinnitusData(s);

    if ~ismember('BBN', subj.stimulus.keys)
        fprintf('⚠️ %s: no BBN key\n', subj.demographics('SubjectID'));
        continue;
    end

    onsets = subj.stimulus('BBN').onset(:);
    if isempty(onsets), continue; end

    nCh    = size(subj.data, 2);
    roi_ch = roi_fixed(roi_fixed >= 1 & roi_fixed <= nCh);
    if isempty(roi_ch), continue; end

    % per-channel block average
    roi_ch_tc = nan(numel(tEpoch), numel(roi_ch));

    for ii = 1:numel(roi_ch)
        ch     = roi_ch(ii);
        trials = nan(numel(tEpoch), numel(onsets));

        for b = 1:numel(onsets)
            idx0 = round(onsets(b) * fs);
            idx  = idx0 + sampOffsets;

            if idx(1) < 1 || idx(end) > size(subj.data, 1), continue; end

            x = subj.data(idx, ch);
            x = x - mean(x(iBase), 'omitnan');
            trials(:, b) = x;
        end

        roi_ch_tc(:, ii) = mean(trials, 2, 'omitnan');
    end

    tc = mean(roi_ch_tc, 2, 'omitnan');
    if all(isnan(tc)), continue; end

    subj_tc(:, end+1) = tc; %#ok<SAGROW>
    subj_used(end+1)  = string(subj.demographics('SubjectID')); %#ok<SAGROW>

    [~, irelT] = min(tc(iTrS));
    subj_trough_t(s) = tEpoch(iTrS(irelT));

    fprintf('✅ %s | trough @ %.2fs\n', ...
        subj.demographics('SubjectID'), subj_trough_t(s));
end


% Group-level trough window

validT = ~isnan(subj_trough_t);

if isempty(subj_tc)
    error('No valid tinnitus subjects with these channels.');
end

muCurve = mean(subj_tc, 2, 'omitnan');

[~, irelT_mu] = min(muCurve(iTrS));
muTr      = tEpoch(iTrS(irelT_mu));
sdTr      = std(subj_trough_t(validT));
troughWin = [muTr - sdTr, muTr + sdTr];

fprintf('\n=== TINNITUS TROUGH WINDOW (fixed channels) ===\n');
fprintf('Trough time: %.2f ± %.2f s  => [%.2f, %.2f]\n', ...
    muTr, sdTr, troughWin(1), troughWin(2));

% Plot

sem = std(subj_tc, 0, 2, 'omitnan') ./ sqrt(size(subj_tc, 2));

figure('Color', 'w'); hold on;

fill([tEpoch; flipud(tEpoch)], [muCurve - sem; flipud(muCurve + sem)], 'k', ...
     'FaceAlpha', 0.12, 'EdgeColor', 'none');
plot(tEpoch, muCurve, 'k', 'LineWidth', 2);
yline(0, '--');

yl = ylim;

% BBN red shading
patch([0 BBN_dur BBN_dur 0], [yl(1) yl(1) yl(2) yl(2)], [1 0 0], ...
      'FaceAlpha', 0.08, 'EdgeColor', 'none');

% ISR blue shading
patch([BBN_dur min(ISR_end, poststim) min(ISR_end, poststim) BBN_dur], ...
      [yl(1) yl(1) yl(2) yl(2)], [0 0 1], ...
      'FaceAlpha', 0.06, 'EdgeColor', 'none');

% redraw mean on top of shading
fill([tEpoch; flipud(tEpoch)], [muCurve - sem; flipud(muCurve + sem)], 'k', ...
     'FaceAlpha', 0.12, 'EdgeColor', 'none');
plot(tEpoch, muCurve, 'k', 'LineWidth', 2);
yline(0, '--');

xline(0,       ':');
xline(BBN_dur, ':');

% trough window (blue)
xline(troughWin(1), ':', 'Color', [0 0 1], 'LineWidth', 2);
xline(troughWin(2), ':', 'Color', [0 0 1], 'LineWidth', 2);

xlabel('Time (s) from BBN onset');

ylabel('\Delta HbO (\muM, baseline-corrected)');
title(sprintf('Tinnitus ROI HbO  (mean \\pm SEM), n=%d', ...
      size(subj_tc, 2)));
grid on;

legend({'SEM','Mean','Baseline', ...
        'BBN (shaded)','ISR (shaded)', ...
        'BBN onset','BBN offset', ...
        'Trough win start','Trough win end'}, ...
        'Location','best');


%% generate_Fig5_PLOS.m
% PLOS ONE compliant version of the Tinnitus ROI HbO group-mean figure
% (trough-only analysis).
%
% PREREQUISITE: Run your main analysis script first so the workspace contains:
%   tEpoch       (time vector)
%   subj_tc      (time x subjects matrix of subject ROI timecourses)
%   muCurve      (group-mean curve)
%   troughWin    (1x2 vector: [trough_lo trough_hi])
%   BBN_dur      (BBN duration, e.g. 12)
%   ISR_end      (ISR end time, e.g. 29)
%   poststim     (post-stimulus axis end)


% ----- Output -----
outDir = FIG_DIR;
if ~exist(outDir, 'dir'); mkdir(outDir); end
outFile = fullfile(outDir, 'Fig5.tif');

% ----- Figure size: full page width, single panel -----
figWidthIn  = 7.5;
figHeightIn = 4.0;
dpi         = 300;

% ----- Font settings -----
axFont    = 'Arial';
axSize    = 10;
labSize   = 11;
legSize   = 9;

% ----- Recompute SEM -----
sem = std(subj_tc, 0, 2, 'omitnan') ./ sqrt(size(subj_tc,2));
nSubj = size(subj_tc, 2);

% Build figure

fig = figure('Units','inches', ...
             'Position',[1 1 figWidthIn figHeightIn], ...
             'Color','w', ...
             'PaperUnits','inches', ...
             'PaperPosition',[0 0 figWidthIn figHeightIn], ...
             'PaperSize',[figWidthIn figHeightIn]);

ax = axes('Parent', fig); hold(ax,'on');

% --- Determine y-limits from data so shading patches fill correctly ---
yMin = min(muCurve - sem) * 1.15;
yMax = max(muCurve + sem) * 1.30;   % headroom for legend
ylim(ax, [yMin yMax]);

% --- Shaded stimulus periods (drawn FIRST, behind data) ---
hBBN = patch(ax, [0 BBN_dur BBN_dur 0], [yMin yMin yMax yMax], ...
    [1 0.85 0.85], 'FaceAlpha', 0.5, 'EdgeColor','none', ...
    'DisplayName','BBN stimulus');

isrEnd = min(ISR_end, poststim);
hISR = patch(ax, [BBN_dur isrEnd isrEnd BBN_dur], [yMin yMin yMax yMax], ...
    [0.85 0.85 1], 'FaceAlpha', 0.5, 'EdgeColor','none', ...
    'DisplayName','ISR period');

% --- SEM band ---
hSEM = fill(ax, [tEpoch; flipud(tEpoch)], [muCurve-sem; flipud(muCurve+sem)], ...
    [0.5 0.5 0.5], 'FaceAlpha', 0.30, 'EdgeColor','none', ...
    'DisplayName','SEM');

% --- Group mean curve ---
hMean = plot(ax, tEpoch, muCurve, 'k-', 'LineWidth', 2, ...
    'DisplayName','Group mean');

% --- Zero baseline ---
plot(ax, [0 poststim], [0 0], 'k--', 'LineWidth', 0.5, ...
    'HandleVisibility','off');

% --- Trough window (blue dotted verticals) ---
hTr = plot(ax, [troughWin(1) troughWin(1)], [yMin yMax], ':', ...
    'Color',[0 0 0.85], 'LineWidth', 1.8, ...
    'DisplayName','Trough window');
plot(ax, [troughWin(2) troughWin(2)], [yMin yMax], ':', ...
    'Color',[0 0 0.85], 'LineWidth', 1.8, 'HandleVisibility','off');

% --- Axes formatting ---

xlabel(ax, 'Time (s) from BBN onset', 'FontName',axFont, 'FontSize',labSize);
xlim([0 30]);
ylabel(ax, '\Delta HbO (\muM, baseline-corrected)', 'FontName',axFont, 'FontSize',labSize);

set(ax, 'FontName', axFont, 'FontSize', axSize, ...
        'LineWidth', 1, 'TickDir','out', 'Box','on');

% --- Legend (only meaningful entries) ---
lgd = legend(ax, [hMean hSEM hBBN hISR hTr], ...
    'Location','northeast', ...
    'FontName',axFont, 'FontSize',legSize, ...
    'Box','on');
set(lgd, 'Color', [1 1 1 0.9]);


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
% fprintf('\n--- Fig5.tif specs ---\n');
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
% fprintf('  NOTE: n = %d -- include this in the figure caption, not the figure.\n', nSubj);




%% TINNITUS GROUP — PAPER-STYLE HbO PLOTS (ROI + nonROI)

% Uses the SAME logic as your Control plots:

clc;

%  Define ROI + nonROI labels
roi_fixed = [38 41 21 8 24 2];
nonroi_labels = {
    'Occipital_Sup_L', 'Occipital_Mid_L', 'Occipital_Sup_R', 'Occipital_Mid_R', ...
    'Cuneus_L', 'Cuneus_R', 'Lingual_R', 'Calcarine_L', ...
    'Frontal_Mid_L', 'Frontal_Mid_R', ...
    'Frontal_Sup_L', 'Frontal_Sup_R', ...
    'Frontal_Inf_Tri_R', 'Frontal_Inf_Oper_L'
};


% 1) Build nonROI struct (TINNITUS)
nonroi_final_tinn = struct();

for s = 1:numel(TinnitusData)
    subj = TinnitusData(s);
    link = subj.probe.link;

    is_hbo = strcmpi(link.type,'hbo');
    nonroi_idx = find(is_hbo & ismember(link.ROI_full, nonroi_labels));

    if isempty(nonroi_idx)
        fprintf("⚠️ No anatomical non-ROI channels for tinnitus subject %s\n", ...
            subj.demographics('SubjectID'));
        continue;
    end

    nonroi_final_tinn(s).SubjectID = subj.demographics('SubjectID');
    nonroi_final_tinn(s).channels  = nonroi_idx;

    fprintf("✅ %s → %d anatomical non-ROI channels (tinnitus)\n", ...
        nonroi_final_tinn(s).SubjectID, numel(nonroi_idx));
end


% 2) Plot params

fs = TinnitusData(1).Fs;
prestim  = 2;
poststim = 20;

tEpoch = (-prestim:1/fs:poststim)';
sampOffsets = round(tEpoch * fs);

iBase = find(tEpoch >= -1 & tEpoch <= 0);
iPlot = find(tEpoch >= 0 & tEpoch <= poststim);

cBBN = [1 0 0];
cISR = [0 1 0];


% 3) TINNITUS ROI (fixed channels)

ROI_BBN_tinn = [];
ROI_ISR_tinn = [];

for s = 1:numel(TinnitusData)
    subj = TinnitusData(s);

    stim_names = subj.stimulus.keys;
    if ~ismember('BBN',stim_names) || ~ismember('ISR',stim_names)
        continue;
    end

    onBBN = subj.stimulus('BBN').onset(:);
    onISR = subj.stimulus('ISR').onset(:);
    if isempty(onBBN) || isempty(onISR), continue; end

    nCh = size(subj.data,2);
    roi_ch = roi_fixed(roi_fixed >= 1 & roi_fixed <= nCh);
    if isempty(roi_ch), continue; end

    % ----- BBN -----
    tc = nan(numel(tEpoch), numel(onBBN));
    for b = 1:numel(onBBN)
        idx = round(onBBN(b)*fs) + sampOffsets;
        if idx(1)<1 || idx(end)>size(subj.data,1), continue; end
        x = mean(subj.data(idx,roi_ch),2);
        x = x - mean(x(iBase),'omitnan');
        tc(:,b) = x;
    end
    ROI_BBN_tinn(:,end+1) = mean(tc,2,'omitnan'); %#ok<SAGROW>

    % ----- ISR -----
    tc = nan(numel(tEpoch), numel(onISR));
    for b = 1:numel(onISR)
        idx = round(onISR(b)*fs) + sampOffsets;
        if idx(1)<1 || idx(end)>size(subj.data,1), continue; end
        x = mean(subj.data(idx,roi_ch),2);
        x = x - mean(x(iBase),'omitnan');
        tc(:,b) = x;
    end
    ROI_ISR_tinn(:,end+1) = mean(tc,2,'omitnan'); %#ok<SAGROW>
end

muBBN = mean(ROI_BBN_tinn,2,'omitnan');
muISR = mean(ROI_ISR_tinn,2,'omitnan');

seBBN = std(ROI_BBN_tinn,0,2,'omitnan') / sqrt(size(ROI_BBN_tinn,2));
seISR = std(ROI_ISR_tinn,0,2,'omitnan') / sqrt(size(ROI_ISR_tinn,2));

figure; hold on;
fill([tEpoch(iPlot); flipud(tEpoch(iPlot))], ...
     [muBBN(iPlot)-seBBN(iPlot); flipud(muBBN(iPlot)+seBBN(iPlot))], ...
     cBBN,'FaceAlpha',0.15,'EdgeColor','none');
plot(tEpoch(iPlot), muBBN(iPlot),'Color',cBBN,'LineWidth',2);

fill([tEpoch(iPlot); flipud(tEpoch(iPlot))], ...
     [muISR(iPlot)-seISR(iPlot); flipud(muISR(iPlot)+seISR(iPlot))], ...
     cISR,'FaceAlpha',0.15,'EdgeColor','none');
plot(tEpoch(iPlot), muISR(iPlot),'Color',cISR,'LineWidth',2);

yline(0,'--');
xlabel('Time (s)');
ylabel('\Delta HbO (\muM)');

title(sprintf('Tinnitus — Mean HbO (ROI fixed), n=%d', size(ROI_BBN_tinn,2)));
legend({'Time-locked to BBN onset ','BBN mean','Time-locked to ISR onset','ISR mean'},'Location','southwest');
grid on;


% 4) TINNITUS nonROI (anatomical)

NROI_BBN_tinn = [];
NROI_ISR_tinn = [];

for s = 1:numel(TinnitusData)
    subj = TinnitusData(s);

    if s > numel(nonroi_final_tinn) || isempty(nonroi_final_tinn(s).channels)
        continue;
    end

    onBBN = subj.stimulus('BBN').onset(:);
    onISR = subj.stimulus('ISR').onset(:);
    if isempty(onBBN) || isempty(onISR), continue; end

    nonroi_ch = nonroi_final_tinn(s).channels(:)';

    % ----- BBN -----
    tc = nan(numel(tEpoch), numel(onBBN));
    for b = 1:numel(onBBN)
        idx = round(onBBN(b)*fs) + sampOffsets;
        if idx(1)<1 || idx(end)>size(subj.data,1), continue; end
        x = mean(subj.data(idx,nonroi_ch),2);
        x = x - mean(x(iBase),'omitnan');
        tc(:,b) = x;
    end
    NROI_BBN_tinn(:,end+1) = mean(tc,2,'omitnan'); %#ok<SAGROW>

    % ----- ISR -----
    tc = nan(numel(tEpoch), numel(onISR));
    for b = 1:numel(onISR)
        idx = round(onISR(b)*fs) + sampOffsets;
        if idx(1)<1 || idx(end)>size(subj.data,1), continue; end
        x = mean(subj.data(idx,nonroi_ch),2);
        x = x - mean(x(iBase),'omitnan');
        tc(:,b) = x;
    end
    NROI_ISR_tinn(:,end+1) = mean(tc,2,'omitnan'); %#ok<SAGROW>
end

muBBN = mean(NROI_BBN_tinn,2,'omitnan');
muISR = mean(NROI_ISR_tinn,2,'omitnan');

seBBN = std(NROI_BBN_tinn,0,2,'omitnan') / sqrt(size(NROI_BBN_tinn,2));
seISR = std(NROI_ISR_tinn,0,2,'omitnan') / sqrt(size(NROI_ISR_tinn,2));

figure; hold on;
fill([tEpoch(iPlot); flipud(tEpoch(iPlot))], ...
     [muBBN(iPlot)-seBBN(iPlot); flipud(muBBN(iPlot)+seBBN(iPlot))], ...
     cBBN,'FaceAlpha',0.15,'EdgeColor','none');
plot(tEpoch(iPlot), muBBN(iPlot),'Color',cBBN,'LineWidth',2);

fill([tEpoch(iPlot); flipud(tEpoch(iPlot))], ...
     [muISR(iPlot)-seISR(iPlot); flipud(muISR(iPlot)+seISR(iPlot))], ...
     cISR,'FaceAlpha',0.15,'EdgeColor','none');
plot(tEpoch(iPlot), muISR(iPlot),'Color',cISR,'LineWidth',2);

yline(0,'--');
xlabel('Time (s)');
ylabel('\Delta HbO (\muM)');

title(sprintf('Tinnitus — Mean HbO (Non-ROI anatomical), n=%d', size(NROI_BBN_tinn,2)));
legend({'BBN','BBN mean','ISR','ISR mean'},'Location','southwest');
grid on;


% Freeze for stats

TINN_ROI_BBN  = ROI_BBN_tinn;
TINN_ROI_ISR  = ROI_ISR_tinn;
TINN_NROI_BBN = NROI_BBN_tinn;
TINN_NROI_ISR = NROI_ISR_tinn;


%% NON-ROI LOBE-SPLIT (FRONTAL vs OCCIPITAL) — TINNITUS ONLY


groupName = 'Tinnitus';
Data      = TinnitusData;

frontal_labels = {
    'Frontal_Mid_L','Frontal_Mid_R', ...
    'Frontal_Sup_L','Frontal_Sup_R', ...
    'Frontal_Inf_Tri_R','Frontal_Inf_Oper_L'
};

occipital_labels = {
    'Occipital_Sup_L','Occipital_Mid_L','Occipital_Sup_R','Occipital_Mid_R', ...
     'Cuneus_L','Cuneus_R','Lingual_R','Calcarine_L'
};

fs       = Data(1).Fs;
prestim  = 2;
poststim = 20;

tEpoch      = (-prestim:1/fs:poststim)';
sampOffsets = round(tEpoch * fs);

iBase = find(tEpoch >= -1 & tEpoch <= 1);
iPlot = find(tEpoch >= 0  & tEpoch <= poststim);

cBBN = [1 0 0];
cISR = [0 1 0];

% ---- TINNITUS FRONTAL ----
[TINN_FR_BBN, TINN_FR_ISR, usedFR_tinn] = local_compute_lobe_curves( ...
    Data, frontal_labels, fs, tEpoch, sampOffsets, iBase);

muBBN = mean(TINN_FR_BBN, 2, 'omitnan');
muISR = mean(TINN_FR_ISR, 2, 'omitnan');
seBBN = std(TINN_FR_BBN, 0, 2, 'omitnan') / sqrt(size(TINN_FR_BBN, 2));
seISR = std(TINN_FR_ISR, 0, 2, 'omitnan') / sqrt(size(TINN_FR_ISR, 2));

figure; hold on;
fill([tEpoch(iPlot); flipud(tEpoch(iPlot))], ...
     [muBBN(iPlot)-seBBN(iPlot); flipud(muBBN(iPlot)+seBBN(iPlot))], ...
     cBBN, 'FaceAlpha', 0.15, 'EdgeColor', 'none');
plot(tEpoch(iPlot), muBBN(iPlot), 'Color', cBBN, 'LineWidth', 2);
fill([tEpoch(iPlot); flipud(tEpoch(iPlot))], ...
     [muISR(iPlot)-seISR(iPlot); flipud(muISR(iPlot)+seISR(iPlot))], ...
     cISR, 'FaceAlpha', 0.15, 'EdgeColor', 'none');
plot(tEpoch(iPlot), muISR(iPlot), 'Color', cISR, 'LineWidth', 2);
yline(0, '--', 'Color', [0.5 0.5 0.5]);
xlabel('Time (s)'); ylabel('\Delta HbO (\muM)');
title(sprintf('%s — Frontal nonROI HbO (n=%d)', groupName, size(TINN_FR_BBN, 2)));
legend({'BBN SEM','BBN mean','ISR SEM','ISR mean'}, 'Location', 'southwest');
grid on;

% ---- TINNITUS OCCIPITAL ----
[TINN_OC_BBN, TINN_OC_ISR, usedOC_tinn] = local_compute_lobe_curves( ...
    Data, occipital_labels, fs, tEpoch, sampOffsets, iBase);

muBBN = mean(TINN_OC_BBN, 2, 'omitnan');
muISR = mean(TINN_OC_ISR, 2, 'omitnan');
seBBN = std(TINN_OC_BBN, 0, 2, 'omitnan') / sqrt(size(TINN_OC_BBN, 2));
seISR = std(TINN_OC_ISR, 0, 2, 'omitnan') / sqrt(size(TINN_OC_ISR, 2));

figure; hold on;
fill([tEpoch(iPlot); flipud(tEpoch(iPlot))], ...
     [muBBN(iPlot)-seBBN(iPlot); flipud(muBBN(iPlot)+seBBN(iPlot))], ...
     cBBN, 'FaceAlpha', 0.15, 'EdgeColor', 'none');
plot(tEpoch(iPlot), muBBN(iPlot), 'Color', cBBN, 'LineWidth', 2);
fill([tEpoch(iPlot); flipud(tEpoch(iPlot))], ...
     [muISR(iPlot)-seISR(iPlot); flipud(muISR(iPlot)+seISR(iPlot))], ...
     cISR, 'FaceAlpha', 0.15, 'EdgeColor', 'none');
plot(tEpoch(iPlot), muISR(iPlot), 'Color', cISR, 'LineWidth', 2);
yline(0, '--', 'Color', [0.5 0.5 0.5]);
xlabel('Time (s)'); ylabel('\Delta HbO (\muM)');
title(sprintf('%s — Occipital nonROI HbO (n=%d)', groupName, size(TINN_OC_BBN, 2)));
legend({'BBN SEM','BBN mean','ISR SEM','ISR mean'}, 'Location', 'southwest');
grid on;

fprintf('\nTinnitus Frontal subjects used:\n');  disp(usedFR_tinn(:));
fprintf('Tinnitus Occipital subjects used:\n'); disp(usedOC_tinn(:));


%% ALL STATS for Control and Tinnitus group

% Uses saved matrices:
%   CTRL_ROI_BBN,  CTRL_ROI_ISR,  CTRL_NROI_BBN,  CTRL_NROI_ISR
%   TINN_ROI_BBN,  TINN_ROI_ISR,  TINN_NROI_BBN,  TINN_NROI_ISR

% OUTPUTS:
%   1) Within-group paired t-test (BBN vs ISR) for ROI + nonROI
%   2) Between-group t-test on DELTA = (BBN-ISR) for ROI + nonROI
%   3) Optional: Cohen's dz (paired) + Cohen's d (between on delta)
%   4) Summary table StatsTbl

%  Windows

win = [5 12];   % seconds
iWin = find(tEpoch >= win(1) & tEpoch <= win(2));

fprintf('\n================ ALL STATS (HbO), %g–%gs window ================\n', win(1), win(2));


% 1) Collect datasets

D = struct();

D(1).name = 'CONTROL ROI';
D(1).BBN  = CTRL_ROI_BBN;
D(1).ISR  = CTRL_ROI_ISR;

D(2).name = 'CONTROL nonROI';
D(2).BBN  = CTRL_NROI_BBN;
D(2).ISR  = CTRL_NROI_ISR;

D(3).name = 'TINNITUS ROI';
D(3).BBN  = TINN_ROI_BBN;
D(3).ISR  = TINN_ROI_ISR;

D(4).name = 'TINNITUS nonROI';
D(4).BBN  = TINN_NROI_BBN;
D(4).ISR  = TINN_NROI_ISR;


% 2) Within-group paired stats

StatsTbl = table();

row = 0;

for k = 1:numel(D)

    [subjBBN, subjISR, nUsed] = local_subject_window_means(D(k).BBN, D(k).ISR, iWin);

    row = row + 1;
    StatsTbl.GroupRow{row,1} = D(k).name;
    StatsTbl.n(row,1)        = nUsed;

    if nUsed < 2
        StatsTbl.BBN_mean(row,1) = NaN;  StatsTbl.BBN_se(row,1) = NaN;
        StatsTbl.ISR_mean(row,1) = NaN;  StatsTbl.ISR_se(row,1) = NaN;
        StatsTbl.t_paired(row,1) = NaN;  StatsTbl.p_paired(row,1) = NaN; StatsTbl.df_paired(row,1)=NaN;
        StatsTbl.dz_paired(row,1)= NaN;
        fprintf('\n--- %s ---\nNot enough subjects (n=%d)\n', D(k).name, nUsed);
        continue;
    end

    mBBN = mean(subjBBN);
    mISR = mean(subjISR);

    seBBN = std(subjBBN,0) / sqrt(nUsed);
    seISR = std(subjISR,0) / sqrt(nUsed);

    [~,p,~,st] = ttest(subjBBN, subjISR);   % paired

    dz = (mean(subjBBN - subjISR)) / std(subjBBN - subjISR, 0);  % Cohen's dz

    StatsTbl.BBN_mean(row,1) = mBBN;
    StatsTbl.BBN_se(row,1)   = seBBN;
    StatsTbl.ISR_mean(row,1) = mISR;
    StatsTbl.ISR_se(row,1)   = seISR;
    StatsTbl.t_paired(row,1) = st.tstat;
    StatsTbl.p_paired(row,1) = p;
    StatsTbl.df_paired(row,1)= st.df;
    StatsTbl.dz_paired(row,1)= dz;

    fprintf('\n--- %s ---\n', D(k).name);
    fprintf('BBN mean=%.4f | SE=%.4f\n', mBBN, seBBN);
    fprintf('ISR mean=%.4f | SE=%.4f\n', mISR, seISR);
    fprintf('Paired t-test (BBN vs ISR): t=%.3f | p=%.4g | df=%d | n=%d | dz=%.3f\n', ...
        st.tstat, p, st.df, nUsed, dz);

end


% 3) Between-group stats on delta (BBN-ISR)
%    (Control vs Tinnitus), separately for ROI and nonROI

fprintf('\n================ BETWEEN-GROUP (Control vs Tinnitus) on Δ=(BBN−ISR) ================\n');

% ROI delta
[ctrlBBN, ctrlISR, nC_roi] = local_subject_window_means(CTRL_ROI_BBN, CTRL_ROI_ISR, iWin);
[tinnBBN, tinnISR, nT_roi] = local_subject_window_means(TINN_ROI_BBN, TINN_ROI_ISR, iWin);

dC_roi = ctrlBBN - ctrlISR;
dT_roi = tinnBBN - tinnISR;

[~,p_roi,~,st_roi] = ttest2(dC_roi, dT_roi, 'Vartype','unequal');  % Welch
d_between_roi = (mean(dC_roi)-mean(dT_roi)) / sqrt( ((var(dC_roi,0)) + (var(dT_roi,0))) / 2 ); % Cohen d (pooled-ish)

fprintf('\nROI Δ (BBN−ISR): CTRL n=%d | TINN n=%d\n', nC_roi, nT_roi);
fprintf('CTRL meanΔ=%.4f | TINN meanΔ=%.4f\n', mean(dC_roi), mean(dT_roi));
fprintf('Welch t-test: t=%.3f | p=%.4g | df=%.2f | d=%.3f\n', st_roi.tstat, p_roi, st_roi.df, d_between_roi);

% nonROI delta
[ctrlBBN, ctrlISR, nC_nroi] = local_subject_window_means(CTRL_NROI_BBN, CTRL_NROI_ISR, iWin);
[tinnBBN, tinnISR, nT_nroi] = local_subject_window_means(TINN_NROI_BBN, TINN_NROI_ISR, iWin);

dC_nroi = ctrlBBN - ctrlISR;
dT_nroi = tinnBBN - tinnISR;

[~,p_nroi,~,st_nroi] = ttest2(dC_nroi, dT_nroi, 'Vartype','unequal'); % Welch
d_between_nroi = (mean(dC_nroi)-mean(dT_nroi)) / sqrt( ((var(dC_nroi,0)) + (var(dT_nroi,0))) / 2 );

fprintf('\nnonROI Δ (BBN−ISR): CTRL n=%d | TINN n=%d\n', nC_nroi, nT_nroi);
fprintf('CTRL meanΔ=%.4f | TINN meanΔ=%.4f\n', mean(dC_nroi), mean(dT_nroi));
fprintf('Welch t-test: t=%.3f | p=%.4g | df=%.2f | d=%.3f\n', st_nroi.tstat, p_nroi, st_nroi.df, d_between_nroi);


% 4) Add between-group rows to table

% ROI row
row = row + 1;
StatsTbl.GroupRow{row,1} = 'BETWEEN (Δ ROI) CTRL vs TINN';
StatsTbl.n(row,1)        = min(numel(dC_roi), numel(dT_roi));  % not exact; just a placeholder
StatsTbl.BBN_mean(row,1) = mean(dC_roi);   % store CTRL meanΔ in BBN_mean column
StatsTbl.ISR_mean(row,1) = mean(dT_roi);   % store TINN meanΔ in ISR_mean column
StatsTbl.t_paired(row,1) = st_roi.tstat;
StatsTbl.p_paired(row,1) = p_roi;
StatsTbl.df_paired(row,1)= st_roi.df;
StatsTbl.dz_paired(row,1)= d_between_roi;

% nonROI row
row = row + 1;
StatsTbl.GroupRow{row,1} = 'BETWEEN (Δ nonROI) CTRL vs TINN';
StatsTbl.n(row,1)        = min(numel(dC_nroi), numel(dT_nroi));
StatsTbl.BBN_mean(row,1) = mean(dC_nroi);
StatsTbl.ISR_mean(row,1) = mean(dT_nroi);
StatsTbl.t_paired(row,1) = st_nroi.tstat;
StatsTbl.p_paired(row,1) = p_nroi;
StatsTbl.df_paired(row,1)= st_nroi.df;
StatsTbl.dz_paired(row,1)= d_between_nroi;

disp('================ SUMMARY TABLE ================');
disp(StatsTbl);



% Local helper: subject-wise window means
%   Inputs:  BBN, ISR are (time x subjects)
%   Output:  subjBBN, subjISR are (subjects x 1)



%% SECONDARY (EXPLORATORY) STATS — LATE window - doing window averaging 


lateWin = [13 20];   % seconds (adjustable)

fprintf('\n================ SECONDARY (EXPLORATORY) STATS: %g–%gs ================\n', ...
    lateWin(1), lateWin(2));


% 1) WITHIN-GROUP (paired BBN vs ISR)

fprintf('\n--- WITHIN-GROUP (paired) ---\n');

groups = { ...
    'CONTROL ROI',     CTRL_ROI_BBN,     CTRL_ROI_ISR; ...
    'CONTROL nonROI',  CTRL_NROI_BBN,    CTRL_NROI_ISR; ...
    'TINNITUS ROI',    TINN_ROI_BBN,     TINN_ROI_ISR; ...
    'TINNITUS nonROI', TINN_NROI_BBN,    TINN_NROI_ISR };

for g = 1:size(groups,1)

    label = groups{g,1};
    BBN   = groups{g,2};
    ISR   = groups{g,3};

    fprintf('\n%s (late window)\n', label);

    [xBBN, xISR, nUsed, usedWin] = local_subject_window_means_timeaware(BBN, ISR, tEpoch, lateWin);

    if isempty(usedWin)
        fprintf('Late window not present in this epoch (available t=%.2f..%.2fs)\n', tEpoch(1), tEpoch(min(numel(tEpoch), size(BBN,1))));
        continue;
    end

    if nUsed < 2
        fprintf('Not enough subjects (n=%d)\n', nUsed);
        continue;
    end

    [~,p,~,st] = ttest(xBBN, xISR);
    dz = mean(xBBN - xISR) / std(xBBN - xISR,0);

    fprintf('Used window: %g–%gs\n', usedWin(1), usedWin(2));
    fprintf('BBN mean=%.4f | ISR mean=%.4f\n', mean(xBBN), mean(xISR));
    fprintf('Paired t-test: t=%.3f | p=%.4g | df=%d | n=%d | dz=%.3f\n', ...
        st.tstat, p, st.df, nUsed, dz);
end


% 2) BETWEEN-GROUP (Δ = BBN − ISR), Welch t-test

fprintf('\n--- BETWEEN-GROUP (Δ late window) ---\n');

% ROI
[cBBN, cISR, nC] = local_subject_window_means_timeaware(CTRL_ROI_BBN, CTRL_ROI_ISR, tEpoch, lateWin);
[tBBN, tISR, nT] = local_subject_window_means_timeaware(TINN_ROI_BBN, TINN_ROI_ISR, tEpoch, lateWin);

dC = cBBN - cISR;
dT = tBBN - tISR;

if numel(dC) >= 2 && numel(dT) >= 2
    [~,p_roi,~,st_roi] = ttest2(dC, dT, 'Vartype','unequal');
    d_roi = (mean(dC)-mean(dT)) / sqrt((var(dC)+var(dT))/2);
    fprintf('\nROI late Δ (BBN−ISR): CTRL n=%d | TINN n=%d\n', numel(dC), numel(dT));
    fprintf('CTRL meanΔ=%.4f | TINN meanΔ=%.4f\n', mean(dC), mean(dT));
    fprintf('Welch t-test: t=%.3f | p=%.4g | df=%.2f | d=%.3f\n', ...
        st_roi.tstat, p_roi, st_roi.df, d_roi);
else
    fprintf('\nROI late Δ: Not enough subjects for between-group test.\n');
end

% nonROI
[cBBN, cISR, nC] = local_subject_window_means_timeaware(CTRL_NROI_BBN, CTRL_NROI_ISR, tEpoch, lateWin);
[tBBN, tISR, nT] = local_subject_window_means_timeaware(TINN_NROI_BBN, TINN_NROI_ISR, tEpoch, lateWin);

dC = cBBN - cISR;
dT = tBBN - tISR;

if numel(dC) >= 2 && numel(dT) >= 2
    [~,p_nroi,~,st_nroi] = ttest2(dC, dT, 'Vartype','unequal');
    d_nroi = (mean(dC)-mean(dT)) / sqrt((var(dC)+var(dT))/2);
    fprintf('\nnonROI late Δ (BBN−ISR): CTRL n=%d | TINN n=%d\n', numel(dC), numel(dT));
    fprintf('CTRL meanΔ=%.4f | TINN meanΔ=%.4f\n', mean(dC), mean(dT));
    fprintf('Welch t-test: t=%.3f | p=%.4g | df=%.2f | d=%.3f\n', ...
        st_nroi.tstat, p_nroi, st_nroi.df, d_nroi);
else
    fprintf('\nnonROI late Δ: Not enough subjects for between-group test.\n');
end

fprintf('\nNOTE: These results are EXPLORATORY and reflect late post-stimulus dynamics.\n');

% Local helper: time-aware window means (safe for different epoch lengths)
% Put at END of your script (below all code), or as a separate file.



%% generate_Fig6_PLOS.m
% Panels:
%   (a) Control ROI
%   (b) Control Non-ROI
%   (c) Tinnitus ROI
%   (d) Tinnitus Non-ROI
%   (e) Tinnitus Frontal Non-ROI
%   (f) Tinnitus Occipital Non-ROI

% PREREQUISITE: Run your full analysis pipeline first so the workspace contains:
%   tEpoch
%   CTRL_ROI_BBN,  CTRL_ROI_ISR
%   CTRL_NROI_BBN, CTRL_NROI_ISR
%   TINN_ROI_BBN,  TINN_ROI_ISR
%   TINN_NROI_BBN, TINN_NROI_ISR
%   TINN_FR_BBN,   TINN_FR_ISR
%   TINN_OC_BBN,   TINN_OC_ISR
% -------------------------------------------------------------------------

% ----- Output -----
outDir = FIG_DIR;
if ~exist(outDir, 'dir'); mkdir(outDir); end
outFile = fullfile(outDir, 'Fig6.tif');

% ----- Figure size: full page width, 3 rows x 2 cols -----
figWidthIn  = 7.5;
figHeightIn = 7.2;
dpi         = 300;

% ----- Font settings (Arial 8-12 pt per PLOS) -----
axFont    = 'Arial';
axSize    = 9;
labSize   = 10;
panelSize = 11;
legSize   = 8;

% ----- Colors -----
cBBN = [0.85 0.10 0.10];
cISR = [0.10 0.65 0.10];

% ----- Datasets in panel order (a,b,c,d,e,f) -----
datasets = {
    CTRL_ROI_BBN,  CTRL_ROI_ISR,  'Control ROI',          '(a)';
    CTRL_NROI_BBN, CTRL_NROI_ISR, 'Control Non-ROI',      '(b)';
    TINN_ROI_BBN,  TINN_ROI_ISR,  'Tinnitus ROI',         '(c)';
    TINN_NROI_BBN, TINN_NROI_ISR, 'Tinnitus Non-ROI',     '(d)';
    TINN_FR_BBN,   TINN_FR_ISR,   'Tinnitus Frontal',     '(e)';
    TINN_OC_BBN,   TINN_OC_ISR,   'Tinnitus Occipital',   '(f)';
};


% Build figure

fig = figure('Units','inches', ...
             'Position',[1 1 figWidthIn figHeightIn], ...
             'Color','w', ...
             'PaperUnits','inches', ...
             'PaperPosition',[0 0 figWidthIn figHeightIn], ...
             'PaperSize',[figWidthIn figHeightIn]);

% ----- Manual axis positions for tight, even layout -----
leftMargin   = 0.085;
rightMargin  = 0.025;
bottomMargin = 0.07;
topMargin    = 0.05;
hGap         = 0.09;
vGap         = 0.075;

panelW = (1 - leftMargin - rightMargin - hGap) / 2;
panelH = (1 - bottomMargin - topMargin - 2*vGap) / 3;

% Compute (left, bottom) for each panel
positions = cell(6,1);
for r = 0:2
    for c = 0:1
        idx = 2*r + c + 1;
        left   = leftMargin + c*(panelW + hGap);
        bottom = 1 - topMargin - panelH - r*(panelH + vGap);
        positions{idx} = [left bottom panelW panelH];
    end
end


% Loop through panels

for d = 1:size(datasets,1)
    BBN_data = datasets{d,1};
    ISR_data = datasets{d,2};
    titleStr = datasets{d,3};
    panelLab = datasets{d,4};

    % --- Local time vector matching this dataset's length ---
    nT_local = size(BBN_data, 1);
    t_local  = tEpoch(1:nT_local);

    % --- Group mean and SEM ---
    muBBN = mean(BBN_data, 2, 'omitnan');
    muISR = mean(ISR_data, 2, 'omitnan');
    seBBN = std(BBN_data, 0, 2, 'omitnan') ./ sqrt(sum(~isnan(BBN_data), 2));
    seISR = std(ISR_data, 0, 2, 'omitnan') ./ sqrt(sum(~isnan(ISR_data), 2));

    % --- Pointwise paired t-test at each second ---
    sig_times = [];
    sig_y     = [];

    for k = 0:20
        if k < t_local(1) || k > t_local(end)
            continue;
        end
        [~, ik] = min(abs(t_local - k));

        bbn_vals = BBN_data(ik, :);
        isr_vals = ISR_data(ik, :);
        valid    = ~isnan(bbn_vals) & ~isnan(isr_vals);

        if sum(valid) >= 2
            [~, p] = ttest(bbn_vals(valid)', isr_vals(valid)');
            if p < 0.05
                sig_times(end+1) = t_local(ik); %#ok<SAGROW>
                y_star = max(muBBN(ik), muISR(ik)) + ...
                         max(seBBN(ik), seISR(ik)) + ...
                         0.005;
                sig_y(end+1) = y_star; %#ok<SAGROW>
            end
        end
    end

    % --- Plot range: 0-20 s ---
    iPlot = find(t_local >= 0 & t_local <= 20);
    tAll  = t_local(iPlot);

    % --- Create axes at manual position ---
    ax = axes('Parent', fig, 'Position', positions{d});
    hold(ax,'on');

    % BBN band + line
    fill(ax, [tAll; flipud(tAll)], ...
         [muBBN(iPlot)-seBBN(iPlot); flipud(muBBN(iPlot)+seBBN(iPlot))], ...
         cBBN, 'FaceAlpha', 0.18, 'EdgeColor', 'none');
    hBBN = plot(ax, tAll, muBBN(iPlot), 'Color', cBBN, 'LineWidth', 1.5);

    % ISR band + line
    fill(ax, [tAll; flipud(tAll)], ...
         [muISR(iPlot)-seISR(iPlot); flipud(muISR(iPlot)+seISR(iPlot))], ...
         cISR, 'FaceAlpha', 0.18, 'EdgeColor', 'none');
    hISR = plot(ax, tAll, muISR(iPlot), 'Color', cISR, 'LineWidth', 1.5);

    % Significance stars
    if ~isempty(sig_times)
        plot(ax, sig_times, sig_y, 'k*', 'MarkerSize', 5);
    end

    % Zero baseline
    plot(ax, [0 20], [0 0], '--', 'Color',[0.5 0.5 0.5], 'LineWidth', 0.5);

    % --- Axes formatting ---
    xlim(ax, [0 20]);

    set(ax, 'FontName', axFont, 'FontSize', axSize, ...
            'LineWidth', 0.8, 'TickDir','out', 'Box','on');

    % X-label only on bottom row
    if d >= 5
        xlabel(ax, 'Time (s)', 'FontName',axFont, 'FontSize',labSize);
    end
    % Y-label only on left column
    if mod(d,2) == 1
        ylabel(ax, '\Delta HbO (\muM)', 'FontName',axFont, 'FontSize',labSize);
    end

    % Sub-title above each panel
    title(ax, titleStr, ...
    'FontName',axFont, 'FontSize',axSize, 'FontWeight','normal');

    % Panel label (a), (b), etc.
    text(ax, -0.18, 1.12, panelLab, ...
        'Units','normalized', ...
        'FontName',axFont, 'FontSize',panelSize, 'FontWeight','bold');

    % Legend on panel (b) where there's empty space
    if d == 2
        lgd = legend(ax, [hBBN hISR], {'Time-locked to BBN onset','Time-locked to ISR onset'}, ...
            'Location','northwest', ...
            'FontName',axFont, 'FontSize',legSize, ...
            'Box','on');
        set(lgd, 'Color', [1 1 1 0.9]);
    end
end


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
% fprintf('\n--- Fig6.tif specs ---\n');
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


%% =========================================================================
%  ASSEMBLE AND SAVE HEMODYNAMIC OUTPUTS
%
%  Subject-wise window means for every region and condition are written to
%  a single struct, HR, in HR_outputs.mat. These are the values the
%  questionnaire analysis correlates against THI and TFI scores.
%
%  Early window 5-12 s, late window 13-20 s, both relative to BBN onset.
% =========================================================================

earlyWin = [5 12];
lateWin  = [13 20];

iEarly = find(tEpoch >= earlyWin(1) & tEpoch <= earlyWin(2));
iLate  = find(tEpoch >= lateWin(1)  & tEpoch <= lateWin(2));

% ---- subject IDs in the same column order as the curve matrices ----
ids_ctrl = strings(0,1);
for s = 1:numel(ControlData)
    if s > numel(roi_final) || isempty(roi_final(s).channels), continue; end
    ids_ctrl(end+1,1) = string(ControlData(s).demographics('SubjectID'));   %#ok<SAGROW>
end

ids_tinn = strings(0,1);
for s = 1:numel(TinnitusData)
    subj = TinnitusData(s);
    k = subj.stimulus.keys;
    if ~ismember('BBN',k) || ~ismember('ISR',k), continue; end
    if isempty(subj.stimulus('BBN').onset) || isempty(subj.stimulus('ISR').onset), continue; end
    if isempty(roi_fixed(roi_fixed >= 1 & roi_fixed <= size(subj.data,2))), continue; end
    ids_tinn(end+1,1) = string(subj.demographics('SubjectID'));             %#ok<SAGROW>
end

% ---- helper: window means for one BBN/ISR matrix pair ----
winMeans = @(M, idx) mean(M(idx, :), 1, 'omitnan')';

HR = struct();
HR.Info.EarlyWindow_s      = earlyWin;
HR.Info.LateWindow_s       = lateWin;
HR.Info.TinnitusROIChannels = roi_fixed;
HR.Info.Units              = 'micromolar change in HbO, baseline corrected';
HR.Info.Created            = datetime('now');

% ---- TINNITUS GROUP ----
HR.Tinnitus.SubjectID = ids_tinn;

HR.Tinnitus.ROI = table( ids_tinn, ...
    winMeans(TINN_ROI_BBN, iEarly), winMeans(TINN_ROI_ISR, iEarly), ...
    winMeans(TINN_ROI_BBN, iLate),  winMeans(TINN_ROI_ISR, iLate), ...
    'VariableNames', {'SubjectID','early_BBN','early_ISR','late_BBN','late_ISR'});
HR.Tinnitus.ROI.early_delta = HR.Tinnitus.ROI.early_BBN - HR.Tinnitus.ROI.early_ISR;
HR.Tinnitus.ROI.late_delta  = HR.Tinnitus.ROI.late_BBN  - HR.Tinnitus.ROI.late_ISR;

HR.Tinnitus.NonROI = table( ...
    winMeans(TINN_NROI_BBN, iEarly), winMeans(TINN_NROI_ISR, iEarly), ...
    winMeans(TINN_NROI_BBN, iLate),  winMeans(TINN_NROI_ISR, iLate), ...
    'VariableNames', {'early_BBN','early_ISR','late_BBN','late_ISR'});

HR.Tinnitus.Frontal = table( ...
    winMeans(TINN_FR_BBN, iEarly), winMeans(TINN_FR_ISR, iEarly), ...
    winMeans(TINN_FR_BBN, iLate),  winMeans(TINN_FR_ISR, iLate), ...
    'VariableNames', {'early_BBN','early_ISR','late_BBN','late_ISR'});

HR.Tinnitus.Occipital = table( ...
    winMeans(TINN_OC_BBN, iEarly), winMeans(TINN_OC_ISR, iEarly), ...
    winMeans(TINN_OC_BBN, iLate),  winMeans(TINN_OC_ISR, iLate), ...
    'VariableNames', {'early_BBN','early_ISR','late_BBN','late_ISR'});

% ---- CONTROL GROUP ----
HR.Control.SubjectID = ids_ctrl;

HR.Control.ROI = table( ids_ctrl, ...
    winMeans(CTRL_ROI_BBN, iEarly), winMeans(CTRL_ROI_ISR, iEarly), ...
    winMeans(CTRL_ROI_BBN, iLate),  winMeans(CTRL_ROI_ISR, iLate), ...
    'VariableNames', {'SubjectID','early_BBN','early_ISR','late_BBN','late_ISR'});
HR.Control.ROI.early_delta = HR.Control.ROI.early_BBN - HR.Control.ROI.early_ISR;
HR.Control.ROI.late_delta  = HR.Control.ROI.late_BBN  - HR.Control.ROI.late_ISR;

HR.Control.NonROI = table( ...
    winMeans(CTRL_NROI_BBN, iEarly), winMeans(CTRL_NROI_ISR, iEarly), ...
    winMeans(CTRL_NROI_BBN, iLate),  winMeans(CTRL_NROI_ISR, iLate), ...
    'VariableNames', {'early_BBN','early_ISR','late_BBN','late_ISR'});

% ---- group-mean curves, kept so the figures can be redrawn ----
HR.Curves.tEpoch        = tEpoch;
HR.Curves.CTRL_ROI_BBN  = CTRL_ROI_BBN;
HR.Curves.CTRL_ROI_ISR  = CTRL_ROI_ISR;
HR.Curves.CTRL_NROI_BBN = CTRL_NROI_BBN;
HR.Curves.CTRL_NROI_ISR = CTRL_NROI_ISR;
HR.Curves.TINN_ROI_BBN  = TINN_ROI_BBN;
HR.Curves.TINN_ROI_ISR  = TINN_ROI_ISR;
HR.Curves.TINN_NROI_BBN = TINN_NROI_BBN;
HR.Curves.TINN_NROI_ISR = TINN_NROI_ISR;
HR.Curves.TINN_FR_BBN   = TINN_FR_BBN;
HR.Curves.TINN_FR_ISR   = TINN_FR_ISR;
HR.Curves.TINN_OC_BBN   = TINN_OC_BBN;
HR.Curves.TINN_OC_ISR   = TINN_OC_ISR;

save(fullfile(DATA_DIR, 'HR_outputs.mat'), 'HR');

fprintf('\nHemodynamic outputs saved to %s\n', fullfile(DATA_DIR, 'HR_outputs.mat'));
fprintf('  tinnitus n = %d | control n = %d\n', ...
    numel(HR.Tinnitus.SubjectID), numel(HR.Control.SubjectID));


%% =========================================================================
%  LOCAL FUNCTIONS
%  MATLAB requires all local functions to appear at the end of a script.
% =========================================================================

function [LOB_BBN, LOB_ISR, usedIDs] = local_compute_lobe_curves(Data, labelList, fs, tEpoch, sampOffsets, iBase)

LOB_BBN = [];
LOB_ISR = [];
usedIDs = {};   % use cell array instead of string array to avoid size issues

for s = 1:numel(Data)
    subj = Data(s);

    stim_names = subj.stimulus.keys;
    if ~ismember('BBN', stim_names) || ~ismember('ISR', stim_names)
        continue;
    end

    onBBN = subj.stimulus('BBN').onset(:);
    onISR = subj.stimulus('ISR').onset(:);
    if isempty(onBBN) || isempty(onISR)
        continue;
    end

    link   = subj.probe.link;
    is_hbo = strcmpi(link.type, 'hbo');

    if ~ismember('ROI_full', link.Properties.VariableNames)
        error('probe.link does not contain ROI_full.');
    end

    lob_idx = find(is_hbo & ismember(link.ROI_full, labelList));
    if isempty(lob_idx)
        continue;
    end

    % ----- BBN -----
    tc = nan(numel(tEpoch), numel(onBBN));
    for b = 1:numel(onBBN)
        idx = round(onBBN(b) * fs) + sampOffsets;
        if idx(1) < 1 || idx(end) > size(subj.data, 1), continue; end
        x = mean(subj.data(idx, lob_idx), 2, 'omitnan');
        x = x - mean(x(iBase), 'omitnan');
        tc(:, b) = x;
    end
    subjBBN = mean(tc, 2, 'omitnan');

    % ----- ISR -----
    tc = nan(numel(tEpoch), numel(onISR));
    for b = 1:numel(onISR)
        idx = round(onISR(b) * fs) + sampOffsets;
        if idx(1) < 1 || idx(end) > size(subj.data, 1), continue; end
        x = mean(subj.data(idx, lob_idx), 2, 'omitnan');
        x = x - mean(x(iBase), 'omitnan');
        tc(:, b) = x;
    end
    subjISR = mean(tc, 2, 'omitnan');

    if all(isnan(subjBBN)) || all(isnan(subjISR))
        continue;
    end

    LOB_BBN(:, end+1) = subjBBN; %#ok<SAGROW>
    LOB_ISR(:, end+1) = subjISR; %#ok<SAGROW>

    % ---- Safe SubjectID retrieval ----
    sid = sprintf('Subject_%d', s);   % default fallback
    demKeys = subj.demographics.keys;
    for kk = 1:numel(demKeys)
        if strcmpi(demKeys{kk}, 'SubjectID')
            val = subj.demographics(demKeys{kk});
            if ischar(val) || isstring(val)
                sid = char(val);
            end
            break;
        end
    end
    usedIDs{end+1} = sid; %#ok<AGROW>

end

% convert cell to string array at the end
usedIDs = string(usedIDs);

end

function [subjBBN, subjISR, nUsed] = local_subject_window_means(BBN, ISR, iWin)

% empty guard
if isempty(BBN) || isempty(ISR)
    subjBBN = []; subjISR = []; nUsed = 0;
    return;
end

% match #subjects
nB = size(BBN,2);
nI = size(ISR,2);
n  = min(nB,nI);
BBN = BBN(:,1:n);
ISR = ISR(:,1:n);

% subject-wise mean in window
subjBBN = mean(BBN(iWin,:), 1, 'omitnan')';
subjISR = mean(ISR(iWin,:), 1, 'omitnan')';

% remove NaN subjects
ok = isfinite(subjBBN) & isfinite(subjISR);
subjBBN = subjBBN(ok);
subjISR = subjISR(ok);

nUsed = numel(subjBBN);

end

function [subjBBN, subjISR, nUsed, usedWin] = local_subject_window_means_timeaware(BBN, ISR, tEpoch, win)

subjBBN = [];
subjISR = [];
nUsed   = 0;
usedWin = [];

% Guard
if isempty(BBN) || isempty(ISR)
    return;
end

% Match time length safely
nT = min([size(BBN,1), size(ISR,1), numel(tEpoch)]);
t  = tEpoch(1:nT);

% Window indices within available t
iWin = find(t >= win(1) & t <= win(2));
if isempty(iWin)
    return;
end
usedWin = [t(iWin(1)) t(iWin(end))];

% Match subjects
nSub = min(size(BBN,2), size(ISR,2));
BBN  = BBN(1:nT, 1:nSub);
ISR  = ISR(1:nT, 1:nSub);

% Subject means in window
subjBBN = mean(BBN(iWin,:), 1, 'omitnan')';
subjISR = mean(ISR(iWin,:), 1, 'omitnan')';

% Remove NaN subjects
ok = isfinite(subjBBN) & isfinite(subjISR);
subjBBN = subjBBN(ok);
subjISR = subjISR(ok);
nUsed   = numel(subjBBN);

end