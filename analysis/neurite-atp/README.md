# neurite-atp

`neurite_atp.py` measures the ATP/ADP sensor ratio (vLH92, 488/405 excitation) along neurites.
First used on `20261002_MN_92` (`neurite_*` images: 7 dishes, 6–7 XY fields each).

## Run

Edit the PARAMETERS block at the top of the script, then `python neurite_atp.py`.
Needs numpy, pandas, scipy, scikit-image, tifffile and matplotlib.

- Batch: list experiment folders in `EXPERIMENT_FOLDERS`, and/or parent folders in `BATCH_ROOTS`.
  Every folder below a root (up to `SEARCH_DEPTH` levels) with `neurite_*` images in `tif/` is
  added; `SKIP_FOLDERS` leaves folders out by name. Each folder is one experiment. Results of all
  folders are also pooled (`POOLED_OUTPUT`, default `<first root>/neurite_atp_pooled/`), with one
  big dot per experiment + replicate.
- Images are read from `<folder>/tif/` (`IMAGE_SUBFOLDER`), only files starting with `neurite_`
  (`FILE_PREFIX`).
- `TREATMENT_FILE`: `analysis/fileid_treatment.csv` in each experiment folder. Its `image` column
  is the image name without `_XY..` (for example `neurite_vLH92_MN_60x-1xzoom_0000`). An optional
  `replicate` column sets the replicate; otherwise it is "A".
- Outputs go to `<folder>/analysis/neurite_atp/` (without `BATCH_ROOTS`, the pooled results go to
  `<first folder>/analysis/neurite_atp/pooled/`).

## Method

1. Background: per channel, a smooth map (grey opening with a 15 px disk, then Gaussian blur)
   is subtracted.
2. Neurite mask: Sato ridge filter on the summed signal, Li threshold, signal above 3× noise.
   Blobs wider than ~12 px (somata, clumps, debris) and pieces shorter than 40 px are removed.
3. Segments: the mask skeleton is cut at branch points and crossings. Each piece ≥ 20 px
   (≈ 3.7 µm) is one neurite segment, and mask pixels go to the nearest segment. Single cells
   cannot be followed in these dense networks, so the segment is the per-"cell" unit.
4. Ratio: only pixels above 3× noise in both channels are used.
   `ratio_median` = median of the per-pixel 488/405 ratios, which is what the plots use
   (`RATIO_COLUMN`). `ratio_sum` = Σ(488 − bg) / Σ(405 − bg) (intensity-weighted) is also saved.
5. Normalisation: each value is divided by the median of its batch's control. Controls follow
   `ATP_CONTROLS` (2DG → "0mM 2DG", MMS → "0mM MMS", gluc → "2.5mM gluc"). If "0mM 2DG" is
   missing, "2.5mM gluc" is used instead (`CONTROL_FALLBACK`). The control with the same
   `treat_time` is used when there is one (`MATCH_TREAT_TIME`).

## Outputs

| File | Contents |
|---|---|
| `neurite_segments.csv` | one row per segment: image, treatment, length (px, µm), intensities, `ratio_sum`, `ratio_median`, `ratio_norm` |
| `neurite_images.csv` | one row per image: whole-mask ratio, background, noise, mask coverage, saturated pixels, n segments |
| `summary_by_treatment.csv` | per condition: n images and segments, raw and normalised medians |
| `stats_vs_control.csv` | Mann-Whitney U against the control, at image level and at segment level |
| `superplot_segments.png` | segments (faint dots), boxplots, image medians (open circles), experiment + replicate median (big dot) |
| `superplot_images.png` | the same with one dot per image |
| `ratio_vs_intensity.png` | QC: segment ratio vs 405 intensity |
| `overlays/` | QC per image: mask outline, and segments coloured by ratio |

With one experiment, use the image-level p-values: segments from the same field are not
independent, so segment-level p-values are far too small.
