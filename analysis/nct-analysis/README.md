# What's in MN_no-glu_110 and _111

## What the files are
Each file is an export of Fiji/ImageJ **Measure** results with experiment labels added. There is **one row per cell (ROI) per time frame**. The labels in the file names and the `Label` column are `com_<reporter>_MN_no-glu_50%_1s_30s_60x_<field>_XY<n>_561.tif`. *Inferred* reading: 561 nm (red) channel, 50% laser, 1 s exposure, a frame every 30 s, 60x objective, glucose-free medium.

| | file 110 | file 111 |
|---|---|---|
| reporter | 110 | 111 |
| cells (ROIs) | 426 | 450 |
| frames | 33 (0–960 s) | 28 (0–810 s) |
| rows | 14,058 | 12,600 |

Both files use the same design: a 2-deoxyglucose (2DG) dose of 0, 2, 6, 10 or 50 mM, a 1 h pretreatment, pyruvate present, replicate A, and 8 XY fields per dose with about 75–96 cells per dose.

## Columns
- **ImageJ measurements:** `Area` (pixels), `Mean`, `StdDev`, `Min`, `Max`, `Median`, `Kurt`, `IntDen`, `RawIntDen`. Intensities are raw camera grey levels. Area changes slightly from frame to frame (about 2% CV), so ROIs were likely re-segmented or tracked each frame.
- **Labels:** `Image`, `roi`, `Slice` (frame number), `file_id` (0, 2, 4, 6, 8 map to 50, 10, 6, 2, 0 mM, so the 50 mM dish was imaged first), `category`=MN, `position`=XY field, `reporter`, `treatment`, `treat_time`, `pyruvate`, `replicate`, `group`, `img_id`, `order`, `roi_uni` (unique cell id within a file).
- **Derived columns (checked against the data):**
  - `bg_mean` is the background intensity, one value per image frame (about 102).
  - `time` = (Slice−1)×30 s.
  - `Mean_SN` = Mean − bg_mean, and `Median_SN` = Median − bg_mean.
  - `pre_mean` and `pre_median` are the averages of the first 4 frames (0–90 s), i.e. the baseline.
  - **`norm_mean`** = Mean_SN / (pre_mean − bg_mean) is the background-subtracted signal relative to baseline, with 1.0 meaning no change. `norm_median` is built the same way.
  - `count_1` = the number of frames where `norm_median` is exactly 1, i.e. where Median equals pre_median. This is confirmed by the user's code and matches 100% of cells in both files.
    - It is an artifact of integer rounding. Median is a whole number, while pre_median is the average of 4 medians. Cells whose pre_median isn't a whole number (about 62%) can never score, which is why ~63% of cells have 0.
    - High counts mostly flag dim cells whose median barely moves (a range of 1 to 1.5 grey levels), and 32 to 39% of the hits fall in the baseline frames.
    - It is not a reliable "unchanged frames" metric. A tolerance-based count after 90 s would be a better alternative.
- The unnamed first column is a row id from a larger master table. File 110 covers ids 1–24714 and 111 covers 24721–47486.

## Main patterns
- Both reporters are flat for frames 1–4 and then change abruptly from about 120 s. That pattern suggests something was added right after the 90 s frame (*inferred*: probably a stimulus or a glucose add-back).
- **Reporter 110 goes up.** Per-cell plateau over the last 120 s: 0 mM gives 1.19×, 2 mM 1.14×, 6 mM 1.12×, 10 mM 1.08× and 50 mM 1.04×. That is a clear dose trend (Kruskal–Wallis p≈4e-11, Spearman ρ=−0.35). About half the cells rise by more than 10%.
- **Reporter 111 goes down.** Plateau values are 0.60–0.68× for 0–10 mM and 0.77× for 50 mM (KW p≈1e-5). 77% of cells fall by more than 10%. Only 50 mM clearly dampens the drop. Below that, the doses are not in order.
- **So 2DG blunts both responses, and most strongly at 50 mM.** The pattern is consistent across the 8 fields (see per-field table in `per_cell_summary.csv`).
- **Caveat for 111:** baseline brightness strongly predicts the drop (Spearman ρ=−0.66). Dim cells sit close to background, so their ratios are noisy. Consider a brightness cutoff.
- **Caveat for both:** dishes were imaged in sequence (50 mM first, 0 mM last), so dose is confounded with imaging order and time in no-glucose medium.

## Files here
- `combined_long.csv`: both files stacked, plus `dose_mM` and `cell_id` (= reporter_roi_uni).
- `per_cell_summary.csv`: one row per cell with baseline, area, plateau_norm (mean of the last 5 frames), peak, min, value at 180 s, and count_1.
- `timecourses.png`: mean ± SEM of norm_mean per dose for each reporter.

## Transport rates (added after the user explained the experiment)
The user explained the setup. These are nuclear intensities of photoactivatable reporters, with activation at frame 4 (t = 90 s). Reporter 110 reports import and 111 reports export. Both rates use `norm_mean` with t0 = 90 s.

- **Initial rate (model-free):** a linear slope over frames 4–7 (90–180 s), in % of baseline per minute. This is the headline for import, because the import traces rise slowly and many don't plateau within the movie.
- **Exponential fit:** y = y_inf + (1 − y_inf)·exp(−k(t − 90)). The outputs are k, t½ = ln2/k and amplitude = y_inf − 1.
  - k is the relaxation rate to a new steady state (k_in + k_out), not a pure one-way rate.
  - This works well for export (median R² = 0.99 in bright cells) but poorly for import.
- **QC:**
  - `bright` = baseline − bg ≥ 20 grey levels. Below that, R² collapses.
  - `responder` = |amplitude| ≥ 0.05 in the expected direction.
  - `good_fit` = bright & responder & R² ≥ 0.8 & k not at a fit bound.

| 2DG mM | import slope %/min | % import responders | export slope %/min | export t½ s | export amplitude |
|---|---|---|---|---|---|
| 0 | 6.3 | 96 | −27.0 | 44 | −0.58 |
| 2 | 4.7 | 89 | −27.0 | 44 | −0.56 |
| 6 | 2.2 | 77 | −23.0 | 53 | −0.53 |
| 10 | 2.7 | 62 | −18.3 | 63 | −0.50 |
| 50 | 0.7 | 39 | −10.9 | 107 | −0.43 |

All values are medians of bright cells (t½ and amplitude use good fits). The dose effect is significant for both reporters, on cells and also on field medians (KW p ≤ 1e-4 for slopes).

**Caveats:** there is one replicate (A) per condition, so cells are not independent biological replicates. Dose is also confounded with imaging order.

Files: `rates_per_cell.csv` (one row per cell with QC flags), `rates_by_dose.png`.

### Log-linear fits (`log_fits.png`, `log_fits_plateau_subtracted.png`)
- A plain ln(norm) vs time fit is linear only over the first ~90 s after activation. After that it bends, because both reporters plateau instead of going to 0.
- The plateau-subtracted log ln((norm − y_inf)/(1 − y_inf)) is linear. Its slope is −k, stored per cell in `rates_per_cell.csv` as `k_lnz_per_min`. The fit window is 90–600 s for import and 90–300 s for export, using good-fit cells only.

| 2DG mM | import k (/min) | n | export k (/min) | n |
|---|---|---|---|---|
| 0 | 0.26 | 38 | 0.87 | 47 |
| 2 | 0.14 | 42 | 0.82 | 50 |
| 6 | 0.12 | 34 | 0.72 | 46 |
| 10 | 0.21 | 22 | 0.66 | 52 |
| 50 | 0.12 | 15 | 0.37 | 48 |

- Import k uses only the cells that import well enough to fit. At high dose most cells don't import at all, so this k is biased toward the best importers. For import, initial rate and % responders capture the dose effect better.

### Better import readouts (all bright cells, no fit filter)
k from well-fit cells hides the inhibition, because the non-importing cells at 50 mM are excluded. The measures below include every bright cell (medians, n = 54/63/66/63/69):

| 2DG mM | gain at 600 s (norm − 1) | AUC above baseline (fraction·min) | population-fit initial flux k·amp (%/min) |
|---|---|---|---|
| 0 | 0.21 | 2.57 | 7.4 |
| 2 | 0.17 | 2.02 | 4.5 |
| 6 | 0.10 | 1.28 | 2.5 |
| 10 | 0.07 | 0.92 | 3.8 |
| 50 | 0.02 | 0.27 | 1.2 |

KW p ≈ 1e-18 to 1e-20; every dose differs from 0 mM (MWU p < 0.015).

### Initial flux = k × amplitude (recommended rate for both reporters)
- **Definition.** This is the slope of the fitted exponential at activation, in % of baseline per min. It is stored per cell in `rates_per_cell.csv` as `flux_pct_per_min`.
- **Fit.** It comes from a refit with k capped at 1/30 s, the frame-rate limit. Without the cap, flat cells could fit a spurious "instant jump" and get a large flux.
- **Coverage.** It is computed for all bright cells, with no responder or R² filter. Flat cells naturally get a flux near 0.

| 2DG mM | import flux median (% of 0 mM) | export flux median (% of 0 mM) |
|---|---|---|
| 0 | 7.3 (100) | −56 (100) |
| 2 | 5.2 (72) | −53 (95) |
| 6 | 3.1 (42) | −38 (69) |
| 10 | 4.0 (55) | −31 (56) |
| 50 | 0.5 (7) | −13 (23) |

- KW p ≈ 1e-19 for both reporters, and the trend holds on field medians (p ≤ 2e-5).
- Every import dose differs from 0 mM. Export differs from 6 mM upward.
- Figure: `flux_by_dose.png`.

### Photobleaching correction (option 3: bleaching term in the fit)
- **Model.** Each cell is fitted over all frames with y(t) = A·e^(−b·t)·S(t), where S(t) = 1 before 90 s and y∞ + (1 − y∞)·e^(−k(t − 90)) after.
  - b is one bleaching rate shared by the cells it covers.
  - A, y∞ and k are per cell (A ≈ 1 absorbs the baseline normalisation).
- **Estimating b.** b is found by minimising the summed squared error over a grid of b values. This was done two ways: per dish (`bleach_rates_per_dish.csv`), and pooled over all dishes with b ≥ 0.
- **Per-dish b is not identifiable.** It ranges from −0.28 to +0.90 %/min across dishes that were imaged with the same settings, and two values are negative, which bleaching cannot be. The term is absorbing slow transport as well.
- **Pooled b = 0.46 %/min** (0.23 %/frame). It is the same for both reporters (110: 0.51, 111: 0.42), so the pooled value is the one used.
- **Effect on initial flux** (median, % of 0 mM, uncorrected → pooled-corrected):
  - Export: 100/95/70/56/22 → 100/96/70/57/22, essentially unchanged.
  - Import: 100/71/42/55/7 → 100/70/36/53/16. Only 50 mM moves: 0.48 → 1.08 %/min.
- **Import plateau.** After correction the import traces at low dose keep rising to 960 s without levelling off, so the import plateau and k get less certain.
- **Files.** Per-cell results are in `rates_bleach_corrected.csv`: `*_nob` = same model with b = 0, `*_corr` = per-dish b, `*_pool` = pooled b. Figure: `bleach_correction.png`.

### Per-reporter bleach correction (user's preferred approach; reporters are fluorescent before activation)
- **Estimating b.** For each reporter, b = −(mean per-cell slope of ln(norm) over the pre-activation frames), across all bright cells. Results are in `bleach_rate_per_reporter.csv`.
  - Frames 2–4 (used): import −0.04 %/min (95% CI −0.45 to 0.37), export 0.07 (−0.30 to 0.44). Both are effectively 0, so the import value is clipped to 0.
  - Frames 1–4: import −0.35, export 1.02. The export value is driven by a one-off drop from frame 1 to frame 2, not steady bleaching.
  - Fit-based per reporter (option 3, for comparison): import 0.51, export 0.42 %/min.
- **Correction.** corrected = norm(t)·e^(b·(t − 45 s)), then the standard exponential refit with k capped at 1/30 s, on all bright cells. Results are in `flux_bleach_per_reporter.csv`, and the figure is `bleach_per_reporter.png`.
- **Result.** The pre-activation b leaves every flux unchanged. The fit-based b changes only 50 mM import, from 7% to 18% of control. Export is unchanged under every option.

## Reproducible script: `nct_analysis.py`
Runs the whole pipeline from the raw inputs (`combined_Results.csv`, `mCherry_BGInt_Results.csv`, `fileid_treatment.csv`). It does not compute count_1.

Edit the PARAMETERS block at the top of the file, then run `python3 nct_analysis.py`. Settings include the data folder, an optional `BATCH_FOLDERS` list (each folder is analysed separately, with an optional pooled run via `COMBINE_BATCH`), the bleach method, the QC thresholds and `SAVE_CELL_PLOTS`.

- **ROI filter:** an ROI is kept if it is tracked in every frame and has area ≥ 1000 px in every frame. This reproduces the original 876 cells exactly.
- **Background:** each frame uses its own image's background. The original `MN_no-glu_110_30s_all.txt` used reporter 111's background for 110 frames 1–28, which looks like a merge-key collision. The difference is about 0.2 grey levels and the results are unchanged.
- **Outputs** (in `<folder>/nct_results/`; the current run is in `inputs/nct_results/`):
  - tables: `timecourses_long.csv`, `timecourse_mean_by_dose.csv`, `timecourse_<reporter>_wide.csv`, `rates_per_cell.csv`, `summary_by_dose.csv`, `statistics.csv`, `bleach_rates.csv`
  - figures per reporter: `timecourse_`, `flux_`, `slope_`, `half_time_`, `log_fits_<reporter>.png`, plus `bleach_correction.png`. Stars are Mann-Whitney vs 0 mM, Holm-adjusted.
  - dot plots are superplots: grey dots = single cells, open box = all cells, large markers = median of each experiment (folder) + replicate. A missing `replicate` column or value counts as "A". These medians are also saved in `replicate_medians.csv`.
  - bleaching rate is estimated per experiment × reporter, so the pooled run keeps each folder's own rate.
  - `export_vs_import_k.png`, `_flux.png`, `_slope.png`: export rate (x) vs import rate (y) per dose, as median ± SD. k uses well-fit cells; flux and model-free slope use bright cells, with export's sign flipped so both axes are positive. The metrics are set by `RATE_SCATTER_METRICS` and the error bars by `RATE_SCATTER_ERROR`. Faint dots are experiment + replicate medians. Once every dose has ≥ 2 of them, the large dot becomes their mean ± SD.
  - `cell_plots/<reporter>_<dose>mM/` when `SAVE_CELL_PLOTS = True`: each cell's trace with its fitted exponential.
- `script_output/` holds the output of the earlier command-line version.
