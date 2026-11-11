# SSMQI reliability study: analysis code

Code for the preregistered inter- and intra-rater reliability analysis of supraspinatus outlining and echo intensity on 2D shoulder ultrasound.

- Preregistration: https://osf.io/gfd6t (registered October 24, 2026)
- Software: MATLAB R2026a with the Image Processing and Statistics and Machine Learning Toolboxes; Python 3 (standard library only) for the helper scripts.
- No images, tracings, randomization keys, or participant-level data are included. The de-identified image set and tracing masks are available from the corresponding author on reasonable request, subject to the governing institutional review board protocol. The operator-guided tracing software is also available from the corresponding author on request; its method is described step by step in the manuscript and in `docs/SSMQI_Analysis_Algorithm.md`.

## Layout
| Path | What it holds |
|---|---|
| `code/` | The frozen registered analysis (`ssmqi_stats_core.m` holds the Dice, Hausdorff, Bland-Altman, mixed-model ICC, decision-study, and kappa code; `ssmqi_registered_analysis.m`, `ssmqi_load_registered_measurements.m`, `run_ssmqi_registered_analysis.m`, and helpers), three synthetic-data test suites, and `build_results_draft_from_final.py`. |
| `code/exploratory/` | Analyses written after the final registered run and labeled exploratory in the manuscript: selection sensitivity, Rater 2 untraceable-image comparison, and image contrast versus traceability. |
| `docs/` | `SSMQI_Analysis_Algorithm.md` (algorithm step by step), `SSMQI_Registered_Analysis_README.md` (how to run), `Protocol_Deviation_Log.md`, `Analysis_Code_Freeze.csv` (SHA-256 of the frozen files, 2026-11-06 18:48:35), `provenance/` (append-only run log and final-run lock). |
| `results/registered_final/` | Aggregate tables from the one-shot final run (Tables 2 to 6 as CSV) and the exploratory summaries. No pair-level or image-level rows. |

## Check the frozen code
`sha256sum` each of the seven files in `docs/Analysis_Code_Freeze.csv` (all in `code/`); they should match. This repository was assembled after the freeze from unchanged copies, and `.gitattributes` keeps line endings byte-for-byte.

## Tests (synthetic data only)
```matlab
test_ssmqi_stats_core
test_ssmqi_registered_pipeline
test_ssmqi_runner_guards
```
Paths inside the scripts point to the original study folder and must be edited to run elsewhere.

## What was and was not preregistered
Seven analyses were preregistered (see the algorithm document). Everything in `code/exploratory/` was added after the final run, is reported as exploratory, and does not change any registered result.

## License
MIT (see `LICENSE`). The preregistration is at https://osf.io/gfd6t.
