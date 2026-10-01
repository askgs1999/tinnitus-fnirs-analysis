
%   1. computes PTA-4 (500, 1000, 2000, 4000 Hz) for every participant
%      and compares it with the values in the sheet
%   2. computes the Mean, SD and n rows at the bottom of each sheet
%   3. computes group statistics for thresholds, PTA-4, SRT and WRS,
%      with Mann-Whitney (ranksum) p values, Holm and FDR (Benjamini-Hochberg)
%
% Requires the Statistics and Machine Learning Toolbox (for ranksum, ttest2).
% Put the Excel file in the same folder as this script, or give the full path.

clear; clc;

file   = 'Audiograms_by_group.xlsx';
groups = {'Controls', 'Tinnitus'};
freqs  = [250 500 1000 2000 3000 4000 6000 8000];
ptaIdx = [2 3 4 6];          % positions of 500, 1000, 2000, 4000 Hz in freqs
tol    = 1e-6;               % tolerance for numeric comparisons
nFail  = 0;

% Column layout (relative to column B):
%   1-16  thresholds (250 R, 250 L, 500 R, 500 L, ... 8000 R, 8000 L)
%   17-19 PTA-4 R, PTA-4 L, PTA-4 mean
%   20-21 SRT R, SRT L
%   22-25 WRS R (%), WRS R level, WRS L (%), WRS L level
colNames = [reshape(cellstr([compose('%d R',freqs); compose('%d L',freqs)]),1,[]), ...
            {'PTA-4 R','PTA-4 L','PTA-4 mean','SRT R','SRT L', ...
             'WRS R (%)','WRS R level','WRS L (%)','WRS L level'}];

D = struct();

for g = 1:numel(groups)
    sh  = groups{g};
    ids = readcell(file, 'Sheet', sh, 'Range', 'A4:A30');
    M   = cellsToNum(readcell(file, 'Sheet', sh, 'Range', 'B4:Z30'));   % 27 x 25
    S   = cellsToNum(readcell(file, 'Sheet', sh, 'Range', 'B32:Z34'));  % Mean, SD, n

    R = M(:, 1:2:16);   % right-ear thresholds, 27 x 8
    L = M(:, 2:2:16);   % left-ear thresholds,  27 x 8

    fprintf('================ %s (%d participants) ================\n', sh, size(M,1));

    % ---- 1. PTA-4 per participant --------------------------------------
    ptaR = pta4(R, ptaIdx);
    ptaL = pta4(L, ptaIdx);
    ptaM = mean([ptaR ptaL], 2);            % NaN unless both ears available
    mine = [ptaR ptaL ptaM];
    shown = M(:, 17:19);

    bad = ~(abs(mine - shown) < tol | (isnan(mine) & isnan(shown)));
    if any(bad(:))
        [r, c] = find(bad);
        for k = 1:numel(r)
            fprintf('  FAIL  %s  %s: sheet = %g, recomputed = %g\n', ...
                ids{r(k)}, colNames{16 + c(k)}, shown(r(k), c(k)), mine(r(k), c(k)));
        end
        nFail = nFail + nnz(bad);
    else
        fprintf('  PASS  PTA-4 matches for every participant\n');
    end

    % ---- 2. Summary rows (Mean, SD, n) ---------------------------------
    myMean = mean(M, 1, 'omitnan');
    mySD   = std(M, 0, 1, 'omitnan');       % sample SD, same as Excel STDEV
    myN    = sum(~isnan(M), 1);
    summaryBad = false;
    for c = 1:25
        checks = [myMean(c) S(1,c); mySD(c) S(2,c); myN(c) S(3,c)];
        labels = {'Mean', 'SD', 'n'};
        for q = 1:3
            a = checks(q,1); b = checks(q,2);
            if ~(abs(a - b) < tol || (isnan(a) && isnan(b)))
                fprintf('  FAIL  %s of %s: sheet = %g, recomputed = %g\n', ...
                    labels{q}, colNames{c}, b, a);
                summaryBad = true; nFail = nFail + 1;
            end
        end
    end
    if ~summaryBad
        fprintf('  PASS  Mean, SD and n rows match for all 25 columns\n');
    end

    % store for group comparisons
    D.(sh).R    = R;
    D.(sh).L    = L;
    D.(sh).avg  = squeeze(mean(cat(3, R, L), 3, 'omitnan'));  % ear-averaged, one ear if other missing
    D.(sh).pta  = ptaM;
    D.(sh).srt  = mean(M(:, 20:21), 2);     % NaN unless both ears available
    D.(sh).wrsR = M(:, 22);
    D.(sh).wrsL = M(:, 24);
    fprintf('\n');
end

%% ---- 3. Group comparisons --------------------------------------------
C = D.Controls; T = D.Tinnitus;

fprintf('================ Thresholds by frequency (ear-averaged, dB HL) ================\n');
p = nan(1, numel(freqs));
fprintf('%-8s %-20s %-20s %-10s\n', 'Freq', 'Controls', 'Tinnitus', 'p (MWU)');
for k = 1:numel(freqs)
    c = rmNaN(C.avg(:,k)); t = rmNaN(T.avg(:,k));
    p(k) = ranksum(c, t);
    fprintf('%-8d %-20s %-20s %.4f\n', freqs(k), msd(c), msd(t), p(k));
end
pHolm = holm(p);
pBH   = bh(p);
fprintf('\n%-8s %-12s %-12s %-12s\n', 'Freq', 'Uncorrected', 'Holm', 'FDR (BH)');
for k = 1:numel(freqs)
    fprintf('%-8d %-12.4f %-12.4f %-12.4f\n', freqs(k), p(k), pHolm(k), pBH(k));
end

fprintf('\n================ PTA-4, SRT and WRS ================\n');
report('PTA-4 (mean of ears)', C.pta,  T.pta);
report('SRT (mean of ears)',   C.srt,  T.srt);
report('WRS right (%)',        C.wrsR, T.wrsR);
report('WRS left (%)',         C.wrsL, T.wrsL);

fprintf('\n================ Summary ================\n');
if nFail == 0
    fprintf('All checks passed.\n');
else
    fprintf('%d check(s) failed. See FAIL lines above.\n', nFail);
end

%% ---- Local functions --------------------------------------------------
function X = cellsToNum(Cc)
    % Convert a cell array from readcell to numbers. Blank or text cells become NaN.
    X = nan(size(Cc));
    for i = 1:numel(Cc)
        v = Cc{i};
        if isnumeric(v) && isscalar(v)
            X(i) = v;
        end
    end
end

function out = pta4(A, idx)
    % Mean of the four PTA frequencies, NaN if any of them is missing.
    sub = A(:, idx);
    out = mean(sub, 2);
    out(any(isnan(sub), 2)) = NaN;
end

function x = rmNaN(x)
    x = x(~isnan(x));
end

function s = msd(x)
    s = sprintf('%.1f +/- %.1f (n=%d)', mean(x), std(x), numel(x));
end

function report(name, c, t)
    c = rmNaN(c); t = rmNaN(t);
    pM = ranksum(c, t);
    [~, pT] = ttest2(c, t);
    fprintf('%-22s Controls %-22s Tinnitus %-22s MWU p = %.4f, t-test p = %.4f\n', ...
        name, msd(c), msd(t), pM, pT);
end

function padj = holm(p)
    % Holm-Bonferroni adjusted p values
    m = numel(p); [ps, o] = sort(p); padj = zeros(size(p)); run = 0;
    for i = 1:m
        run = max(run, (m - i + 1) * ps(i));
        padj(o(i)) = min(run, 1);
    end
end

function padj = bh(p)
    % Benjamini-Hochberg (FDR) adjusted p values
    m = numel(p); [ps, o] = sort(p); padj = zeros(size(p)); run = 1;
    for i = m:-1:1
        run = min(run, ps(i) * m / i);
        padj(o(i)) = run;
    end
end