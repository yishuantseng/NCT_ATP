# atp-correlation

`atp_transport_correlation.py` correlates the ATP/ADP ratio (488/405 sensor) with nuclear import and
export rates from `../nct-analysis/nct_analysis.py`. Settings are in the PARAMETERS block at the top of
the script; the docstring describes the inputs, matching modes and outputs.

```
python atp_transport_correlation.py
```

Requires numpy, pandas and matplotlib.

## Other treatments

By default (`GROUP_BY = "dose"`), ATP and NCT data are matched by the mM number in the treatment text,
which suits a 2DG dose series. For any other mix of treatments (controls, other drugs, µM doses), set
`GROUP_BY = "treatment"`. Conditions are then matched by their full treatment label, ignoring case and
spaces. Add `"treat_time"` or `"pyruvate"` to `TREATMENT_COLUMNS` to keep those apart, and use
`CONDITION_ORDER` to set the plot order. `TREATMENT_ALIASES` renames whole words before matching, so with the default
`{"glu": "gluc"}` "25mM glu" and "25mM gluc" count as the same condition ("glucose" is left alone).
`IGNORE_REPORTERS` (default `[185]`) leaves reporters out of the analysis.

## Bulk run (no pairing)

Set `BULK_ATP` and `BULK_NCT` to one or more folders. The script searches them recursively for every
`all_images_measurements.csv` and every `nct_results/rates_per_cell.csv` (folders named
`nct_results_pooled/` are skipped so their cells aren't counted twice). It then pools all the cells and
compares ATP/ADP with the transport rates by dose. The results go to `atp_results/bulk/`.

## Excluding bad data

`EXCLUDE_RATES` and `EXCLUDE_ATP` drop rows from specific files. For example,
`{"path_contains": "20250709_MN_glu_110,111", "reporter": 111, "treatment": "control"}` drops the 111
control cells from that experiment only. Leave out a column to match any value.

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
