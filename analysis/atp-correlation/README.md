# atp-correlation

`atp_transport_correlation.py` correlates the ATP/ADP ratio (488/405 sensor) with nuclear import and
export rates from `../nct-analysis/nct_analysis.py`. Settings are in the PARAMETERS block at the top of
the script; the docstring describes the inputs, matching modes and outputs.

```
python atp_transport_correlation.py
```

Requires numpy, pandas and matplotlib.

## Batch runs (NCT and ATP in different folders)

Each experiment needs its NCT results and its ATP file paired. There are two ways to do it:

- **Pairs table.** Set `PAIRS_FILE` to a CSV in the same format as the repo's `folder_pairs.csv`
  (`batch,imaging_folders,kinetics_folders,include`). `imaging_folders` hold
  `all_images_measurements.csv`, and `kinetics_folders` hold the `nct_results/` written by
  `nct_analysis.py`. Folder names are relative to `PAIRS_ROOT`, or to the pairs file's own folder
  when `PAIRS_ROOT` is empty. Separate several folders with `;` to pool them, and set `include` to
  `no` to skip a row.
- **`BATCH_RUNS` list.** Write the pairs straight into the script as
  `{"name": ..., "nct": <path>, "atp": <path>}`. The paths can be absolute and can point anywhere.

Each pair writes its results to `atp_results/<batch>/`.
