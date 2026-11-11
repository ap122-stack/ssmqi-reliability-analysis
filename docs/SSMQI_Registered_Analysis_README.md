# Registered analysis: how it works and how to run it

Preregistration: https://osf.io/gfd6t (registered October 24, 2026).
Study folder: `Reliability Study\Inter_Intra_Rater_Study`.
Code folder: `Scripts`.

## What it does

Reads every saved outline from the four rater sessions (Rater 1 Sessions 1 and 2, Rater 2 Sessions 1 and 2), computes echo intensity inside each outline, and runs the seven registered analyses separately for long-axis and short-axis images:

1. Dice and Jaccard overlap (Dice is the primary outcome).
2. 95th-percentile Hausdorff distance and mean surface distance, in millimeters.
3. Bland-Altman agreement of echo intensity as percent difference, with a test for proportional bias.
4. Standard error of measurement, minimal detectable change (1.96 x sqrt(2) x SEM) and within-subject coefficient of variation.
5. Mixed-model variance components and ICC(A,1).
6. Generalizability decision study (projected reliability when tracings or raters are averaged).
7. Agreement on whether an image could be traced (percent agreement, Cohen kappa, PABAK).

Confidence intervals are 95% percentile bootstrap intervals from 2,000 resamples that draw whole participants (images from one participant are not independent). Pairs where either tracing was untraceable are left out of analyses 1 to 6 and reported in analysis 7.

## Files

| File | Purpose |
|---|---|
| `ssmqi_stats_core.m` | The statistics, one function each, so each can be tested alone. |
| `ssmqi_load_registered_measurements.m` | Collects all tracings into one table and computes echo intensity on the image each rater saw. Scanner-overlay pixels (text, body-position diagram, depth scale) are excluded. |
| `ssmqi_registered_analysis.m` | Runs analyses 1 to 7 and writes Tables 2 to 6, pair-level metrics and Bland-Altman figures. |
| `run_ssmqi_registered_analysis.m` | The controlled entry point (freeze, run, locks, logs). |
| `test_ssmqi_stats_core.m` | 21 checks of each statistic against answers worked out by hand and simulated data with known variance components. |
| `test_ssmqi_registered_pipeline.m` | Runs the whole analysis on synthetic data shaped like the study. |
| `test_ssmqi_runner_guards.m` | Confirms the controls below work. |
| `ssmqi_synthetic_study_table.m` | Builds the synthetic test data. |

## How to run it (in this order, once all four sessions are finished)

```matlab
addpath('C:\Users\Owner\OneDrive\Rotator Cuff MVI\Scripts')
run_ssmqi_registered_analysis('freeze')   % runs the tests, then records code hashes
run_ssmqi_registered_analysis('run')      % the one final run
```

Output goes to `Reliability Study\Inter_Intra_Rater_Study\Analysis\Registered\FINAL_<date>`:
Tables 2 to 6 as CSV files, `Pair_Level_Metrics.csv`, Bland-Altman figures, `Run_Manifest.json` and `Input_Files_SHA256.csv`.

Optional sensitivity analysis: if a physician, blind to the agreement results, flags outlines that are not supraspinatus, list them in `Admin_Private\Adjudication_Flagged_Outlines.csv` (columns `Rater`, `Session`, `Source_Image`) before the run. The runner then also writes a `Sensitivity_ExcludingAdjudicatedOutlines` folder. The registered analysis (all tracings) stays the primary result.

## Controls

- **Freeze:** the code is hashed (SHA-256) after the tests pass. If any analysis file changes afterward, the run stops.
- **Completeness:** the final run refuses to start if any session has pending images.
- **One final run:** after a complete run a lock file is written. A second run needs a written reason, is saved in its own `RERUN_` folder, and is logged.
- **Log:** every execution, including refused runs and incomplete previews, is appended to `Analysis_Run_Log.csv`. Anything viewed before the final run must be disclosed in the manuscript.
- **Fingerprint:** hashes of the case lists, keys and saved tracings used are written with each run.

## Definitions used

- Echo intensity: mean grayscale (0 to 255) inside the outline on the image the rater saw, excluding scanner-overlay pixels.
- HD95: the larger of the two directed 95th-percentile boundary distances between the outlines. Mean surface distance: the average of the two directed mean distances. Distances use the DICOM pixel spacing.
- ICC(A,1): (participant variance + image variance) / (participant + image + rater + residual variance), where the rater variance is the squared rater contrast divided by 2.
- Decision study: reliability with k tracings by the same rater averaged, or with two raters averaged.

## Re-testing

```matlab
test_ssmqi_stats_core
test_ssmqi_registered_pipeline
test_ssmqi_runner_guards
```

The project history of deviations from the registered plan is in `Protocol_Deviation_Log.md`.
