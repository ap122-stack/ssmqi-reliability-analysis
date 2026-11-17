# Reproducibility check, 2026-11-16

The frozen registered analysis (`ssmqi_registered_analysis.m` with `ssmqi_load_registered_measurements.m`) was run five more times on the final data with the registered seed (20261024) and 2000 bootstrap resamples. Before the runs, the SHA-256 hashes of all seven frozen files were checked against `Analysis_Code_Freeze.csv`; all matched.

Each run's six output files (Tables 2 to 6 and the pair-level metrics) were compared byte for byte with the files from the one-shot final run (`FINAL_20261106_190031`). All six were identical in all five runs.

The runs wrote to a scratch folder and called the analysis function directly, so they are not entries in `Analysis_Run_Log.csv` and did not change the final-run folder or the lock file. `Exploratory_Selection_Sensitivity.csv` comes from the separate exploratory script and was not part of this check.
