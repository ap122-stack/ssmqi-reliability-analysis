# SSMQI registered analysis: the algorithm, step by step

Describes exactly what the frozen code does (hashes in `Analysis\Registered\Analysis_Code_Freeze.csv`). Source files: `ssmqi_load_registered_measurements.m`, `ssmqi_stats_core.m`, `ssmqi_registered_analysis.m`, `run_ssmqi_registered_analysis.m`. Registered plan: https://osf.io/gfd6t. MATLAB R2026a (26.1.0.3251617) with the Image Processing and Statistics and Machine Learning Toolboxes. Random seed 20261024, Mersenne Twister, set once at the start of the run.

## Step 0. Controls (runner)

1. `freeze`: runs the three test suites; if all pass, writes SHA-256 hashes of the seven analysis files.
2. `run`: refuses if any hash changed, if any case in any session is still pending, or if a final run already exists (lock file). Every action is appended to `Analysis_Run_Log.csv`. Hashes of case lists, randomization keys and saved tracings are written with the output.

## Step 1. Load every tracing into one table (loader)

For each of the four rater sessions (Rater 1 Sessions 1 and 2, Rater 2 Sessions 1 and 2):

1. Read the session case list (blinded ID, plane, status) and the coordinator's randomization key; every case must have a key row or the run stops.
2. Status "Complete" means the rater accepted a tracing; anything else (Untraceable) has no mask.
3. For a complete case, load the saved binary mask and the pixel spacing (cm per pixel, converted to mm).
4. Rebuild the image the rater saw: read the DICOM, choose the left B-mode region from the ultrasound-region metadata (same rule as the tracing tool; fallback: left half if the image is wider than 1.15 times its height), build the scanner-overlay mask (text, body-position diagram, depth scale), convert to grayscale and **rescale that image linearly to 0 to 255** (per-image min-max).
5. Echo intensity (EI) = mean grayscale of pixels inside the mask and outside the overlay mask. If the mask and image sizes differ, the run stops.
6. Output row: rater, session, case, source image, participant (from the key), plane, complete flag, mask, mm per pixel, EI, pixel count.

## Step 2. Pair tracings and choose comparisons

Five registered comparisons, each run separately for long-axis and short-axis images:

| Label | A | B |
|---|---|---|
| Intra-rater, Rater 1 | R1 S1 | R1 S2 |
| Intra-rater, Rater 2 | R2 S1 | R2 S2 |
| Inter-rater, primary | R1 S1 | R2 S1 |
| Inter-rater, secondary | R1 S2 | R2 S2 |
| Sensitivity | R1 S2 | R2 S1 |

Pairs are formed by matching the same source image in A and B. Analysis 7 uses every matched image. Analyses 1 to 4 use only pairs where both tracings were accepted; a comparison with fewer than 3 usable pairs is skipped.

## Step 3. The seven analyses

**1. Overlap (per usable pair, then summarized).**
Dice = 2|A∩B| / (|A|+|B|). Jaccard = |A∩B| / |A∪B|. Reported as mean with 95% CI and the median.

**2. Boundary distance (mm).**
Take the perimeter pixels of each mask. For each perimeter pixel of A, the Euclidean distance to the nearest perimeter pixel of B (distance transform), and the reverse. HD95 = the larger of the two 95th percentiles of those directed distances. Mean surface distance (MSD) = the average of the two directed means. Distances multiplied by the DICOM pixel spacing.

**3. Bland-Altman of echo intensity, as percent difference.**
d = 100 × (EI_B − EI_A) / mean(EI_A, EI_B); m = pair mean. Bias = mean(d). Limits of agreement = mean(d) ± 1.96 × SD(d). Each has a bootstrap 95% CI. Proportional bias: ordinary least squares of d on m; the slope and its p-value are reported (the p-value does not account for clustering by participant, so it is descriptive).

**4. Reliability of echo intensity.**
Within-subject SD (= SEM) = sqrt(mean(difference²) / 2) from the pair differences (A − B in raw EI units). MDC95 = 1.96 × sqrt(2) × SEM. Within-subject CV (%) = 100 × sqrt(mean((difference / pair mean)² / 2)). Each with a bootstrap 95% CI.

**5. Variance components and ICC(A,1).**
For each plane, all usable tracings from both raters and both sessions are pooled. Linear mixed-effects model fitted by REML: EI ~ Rater + (1 | Participant) + (1 | Image). Variance components: participant, image, residual (which also absorbs session-to-session variation), and rater = δ² / 2, where δ is the fitted rater contrast. ICC(A,1) = (σ²participant + σ²image) / (σ²participant + σ²image + σ²rater + σ²residual), the absolute-agreement, single-measurement form. Thresholds read with Koo and Li: below 0.50 poor, 0.50 to 0.75 moderate, 0.75 to 0.90 good, above 0.90 excellent.

**6. Decision study (G-theory).**
From the same components. Averaging k tracings by one rater: num / (num + σ²rater + σ²residual / k). Two raters, one tracing each: num / (num + σ²rater / 2 + σ²residual / 2), where num = σ²participant + σ²image. Designs reported: 1 tracing; 2 and 3 tracings by the same rater; 2 raters with 1 tracing each. CIs from the bootstrap draws of the components.

**7. Traceability agreement.**
For every matched image, traceable (accepted) or not in each session. Percent agreement = 100 × proportion equal. Cohen kappa = (po − pe) / (1 − pe), where pe is the chance agreement from the two marginal rates (NaN when pe = 1, which happens when both tracings are always traceable). PABAK = 2 × po − 1.

## Step 4. Confidence intervals: cluster bootstrap

Images from one participant are not independent, so resampling draws **whole participants** with replacement. 2,000 resamples. For each: recompute the statistic on all rows of the drawn participants; the 95% CI is the 2.5th and 97.5th percentiles of the 2,000 values (percentile method, no bias correction). For the mixed model, each resample renames the drawn participants and their images so repeated draws count as separate clusters, refits the model, and any resample whose fit fails is skipped (the CI then uses the successful fits). No multiplicity correction (an estimation study; the registration says so).

## Step 5. Outputs

`Table2_Outline_Agreement.csv`, `Table3_EchoIntensity_Agreement.csv`, `Table4_Variance_Components_ICC.csv`, `Table5_Decision_Study.csv`, `Table6_Traceability_Agreement.csv`, `Pair_Level_Metrics.csv`, one Bland-Altman figure per comparison and plane, `Run_Manifest.json`, `Input_Files_SHA256.csv`. Tables 2 to 6 map to the manuscript tables of the same name; `Scripts\build_results_draft_from_final.py` formats them without recomputing anything.

## Step 6. Things a reader should know (not changes; plain properties of the method)

- EI is relative, not absolute: each image is min-max scaled to 0 to 255 before averaging, so EI compares a muscle to the rest of that image's displayed range.
- The ICC residual includes session-to-session variation, because sessions are not a separate random factor in the registered model.
- With 17 participants the percentile intervals will be wide and the bootstrap can be unstable for variance components; failed fits are skipped.
- The proportional-bias p-value ignores clustering.
- Untraceable images drop out of analyses 1 to 6 and appear only in analysis 7. Rater 2 Session 2 has 31 of 80 untraceable, versus 12 in his Session 1, so his comparisons use fewer pairs.
- Traced area is recorded but is not an outcome. The SSMQI composite score is not part of this study.

## Step 7. Tests that guard the code

- `test_ssmqi_stats_core.m`: 21 checks, each statistic against hand-worked answers and simulated data with known variance components (for example ICC near the true 0.82, participant variance near 100, image 50, residual 25, rater 8).
- `test_ssmqi_registered_pipeline.m`: runs the whole pipeline on synthetic data shaped like the study.
- `test_ssmqi_runner_guards.m`: freeze required, tampering detected, incomplete data refused, a single final run, reruns logged.
