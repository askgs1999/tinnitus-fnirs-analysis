# Tinnitus fNIRS analysis code

MATLAB analysis code for:

Satish G, Taylor LM, Arnold MP, Gallagher-Shale J, Krishnamurthy K, Basura GJ. Subjective tinnitus distress correlates with cortical hemodynamic changes in human auditory cortex as measured by functional near-infrared spectroscopy. *PLOS ONE*.

The study compared 27 adults with chronic bilateral non-pulsatile tinnitus against 27 non-tinnitus controls. Participants heard 20 blocks of broadband noise separated by randomised silent intervals, with three minutes of silence before and after the block sequence. Oxygenated haemoglobin was recorded with an 83-channel NIRSport2 system.

## Requirements

- MATLAB R2024b
- Statistics and Machine Learning Toolbox
- [NIRS Brain AnalyzIR Toolbox](https://github.com/huppertt/nirs-toolbox) at commit `995f95e8`, on the MATLAB path. Later versions of the toolbox register the probe to the Colin27 atlas differently and assign different atlas labels to some channels, which changes the control-group ROI, the non-ROI channel sets and the lobe groupings. The published values are reproduced only with this commit.

## Input data

Neither input file is held in this repository.

**fNIRS recordings.** The per-participant NIRx recordings are deposited on OpenNeuro under accession [ds007990](https://openneuro.org/datasets/ds007990). The scripts do not read those directories directly. They read `Aim1_raw_loaded_new_03_17_new.mat`, which holds the recordings as an array of `nirs.core.Data` objects together with the participant demographics table. To build that file, download the OpenNeuro dataset, then uncomment and run the loading block at the top of `RSFC_Code.m`, setting `rootFolder` to the downloaded folder. It calls `nirs.io.loadDirectory` with a group and subject folder hierarchy, populates three missing subject identifiers, and saves the result. All preprocessing is performed afterwards by the scripts.

**Questionnaire responses.** `SubjectiveAndSomatic_DATA_2026-01-15_1203.csv` is the REDCap export of the THI and TFI responses, published with the article as Supporting Information.

Place both files in the same folder as the scripts, or in any folder on the MATLAB path. Each script locates its input automatically and reports the folders it is reading from and writing to before doing anything else.

Four-frequency pure-tone averages and participant ages are not held in the data files. They are entered directly in `Additional_analysis_code`, transcribed from the audiologist's records.

## Scripts

Run these in order. Each writes a `.mat` file that the later scripts read, so the order matters.

### 1. RSFC_Code.m

Preprocesses the recordings, using scalp coupling index screening, optical density conversion, TDDR motion correction, short-separation channel regression, 0.01 to 0.08 Hz band-pass filtering and the modified Beer-Lambert law. Registers the probe to the Colin27 atlas and assigns each channel an AAL cortical label. Computes resting-state functional connectivity during the pre- and post-stimulation silence periods at three levels: ROI seed-based, ROI-seed-to-lobe, and whole-brain channel-pair.

Writes `RSFC_outputs.mat`. Generates Figures 1, 7, 8 and 9.

### 2. Hemodynamic_response_code.m

Repeats the preprocessing, then defines broadband-noise and interstimulus-rest events from the recorded triggers and extracts block-averaged, baseline-corrected epochs. Validates the auditory ROI functionally in the control group and derives the early window, 5 to 12 s, and the late window, 13 to 20 s, from the control group-mean response. Computes group-mean curves for the auditory ROI, the non-ROI spatial control channels, and the frontal and occipital subsets.

Writes `HR_outputs.mat`. Generates Figures 4, 5 and 6.

### 3. Questionnaire_Analysis.m

Scores the Tinnitus Handicap Inventory and the Tinnitus Functional Index, tests the within-session change in each, and correlates both instruments, as absolute scores and as change scores, against the haemodynamic and connectivity measures from the two scripts above. Bonferroni correction is applied within each prespecified family of tests.

Writes `QUEST_outputs.mat`. Generates Figures 3 and 10.

### 4. Additional_analysis_code.m

Supplementary analyses added during peer review: audiometric thresholds and sensation level, baseline distress as a covariate, influence diagnostics for the headline correlation, questionnaire subscale associations, covariate adjustment for age and pure-tone average, score distributions, and a sensitivity analysis of minimum detectable effects.

Writes result tables as CSV files. Produces no figures.

## Output

Each script creates two subfolders beside the data. Figures go to `figures/` as 300 dpi TIFFs and result tables go to `results/` as CSV files.

The TIFF export blocks in the figure sections are commented out, so figures are displayed rather than written. Uncomment those blocks to save the files.

## Test selection

TFI scores were normally distributed and are analysed with Pearson correlations. THI scores were not, and are analysed with Spearman correlations. This rule is applied consistently throughout, and the test used is printed alongside every coefficient.

## Notes

Preprocessing is repeated in scripts 1 and 2 rather than shared, because the two analyses were developed separately and the duplication keeps each script runnable on its own. The preprocessing parameters are the same in both.

To install the required toolbox version, clone the toolbox and check out the commit before adding it to the MATLAB path:

    git clone https://github.com/huppertt/nirs-toolbox.git
    cd nirs-toolbox
    git checkout 995f95e8

To confirm the correct version is in use, run `Hemodynamic_response_code.m` and check the `roi_fixed channel mapping` printout. Channel 8 (source 3, detector 7) should be labelled `Angular_L`. If it is labelled `Occipital_Mid_L`, a different toolbox version is on the path.


## License

MIT. See `LICENSE`.
