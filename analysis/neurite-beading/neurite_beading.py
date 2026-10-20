#!/usr/bin/env python3
"""
Neurite beading (focal swellings / varicosities along neurites), vLH92 motor-neuron neurite fields.

HOW TO USE
  1. Edit the PARAMETERS block below (experiment folders, file names, settings).
  2. Run the file:  python neurite_beading.py   (or press Run in your editor).
  Requires: numpy, pandas, scipy, scikit-image, tifffile, matplotlib.

WHAT IT DOES (per image = one XY field)
  1. Reads the channels in CHANNEL_SUFFIXES (488 and 405 excitation of vLH92) and sums them after
     subtracting a smooth background map from each. vLH92 fills the cytosol, so the summed signal
     follows neurite volume; summing both excitations makes it insensitive to the ATP/ADP ratio.
  2. Neurite mask (same as neurite_atp.py): ridge filter, threshold, remove blobs wider than
     THICK_RADIUS_PX (somata, clumps, debris) and short pieces.
  3. Neurite segments: the skeleton is cut at branch points and crossings; each piece at least
     MIN_SEGMENT_LENGTH_PX long is one segment (a stretch of one neurite between crossings).
  4. Beads: local maxima of the (smoothed) signal on the mask that are at least BEAD_FOLD x the
     neurite shaft around them. The shaft level is the median along the skeleton in a ring
     (BEAD_ANNULUS_PX) around the peak, so it is local and does not depend on how the network was
     cut into segments. The bead outline is the connected area above half-way between shaft and
     peak. Kept if at least BEAD_MIN_AREA_PX, no longer than BEAD_MAX_LENGTH_UM (longer bright
     stretches are thick neurites or bundles), at least BEAD_MIN_SNR x noise above the shaft, and
     not a narrow spot at a crossing (two overlapping neurites are bright too).
     Each bead is assigned to the segment it sits on (for the segment-level table).

BEADING METRICS
  beads_per_100um      number of beads per 100 um of neurite (skeleton) length   <- main readout
  bead_area_fraction   fraction of the neurite mask area that is bead
  intensity_cv         SD / mean of the smoothed intensity along the skeleton (threshold-free:
                       a beaded neurite has an uneven profile, a smooth one a flat profile)
  Image-level values pool all segments of the field (total beads / total length), so they are not
  dominated by short segments; segment-level values are also saved and plotted.

OUTPUTS (in <experiment>/OUTPUT_SUBFOLDER; with several folders, POOLED_OUTPUT gets the combined
         tables and plots, one big dot per experiment + replicate)
  beading_segments.csv     one row per segment: image, treatment, length, n beads, metrics
  beading_beads.csv        one row per bead: position, area, length, width, peak / shaft
  beading_images.csv       one row per image: pooled metrics, total length, blobs removed
  summary_by_treatment.csv per treatment: n images / segments, medians of each metric
  stats_vs_control.csv     Mann-Whitney U of each treatment vs its control (image and segment level)
  superplot_images.png     image-level metrics (faint dots = images), big dots = experiment+replicate
  superplot_segments.png   segment-level metrics (faint dots = segments), blue circles = image medians
  overlays/<image>.png     QC: neurite mask and detected beads
"""

# =============================================================================
# PARAMETERS — edit these
# =============================================================================

# --- Which data to analyse ---------------------------------------------------
# Batch: one folder per experiment. Images are read from <folder>/IMAGE_SUBFOLDER.
# List experiment folders here, and/or put parent folders in BATCH_ROOTS: every folder below them
# (up to SEARCH_DEPTH levels) whose IMAGE_SUBFOLDER holds FILE_PREFIX images is added.
# Folders whose name contains a word in SKIP_FOLDERS are left out.
EXPERIMENT_FOLDERS = [
    "/Volumes/holtl02lab/holtl02labspace/Holt_Lab_Members/Yi-Shuan/microscope/2026/10-12/20261015_MN_92",
]
BATCH_ROOTS = [
    # "/Volumes/holtl02lab/holtl02labspace/Holt_Lab_Members/Yi-Shuan/microscope/2026/10-12",
]
SEARCH_DEPTH = 2
SKIP_FOLDERS = []
IMAGE_SUBFOLDER = "tif"
FILE_PREFIX = "neurite_"                      # only images whose name starts with this
CHANNEL_SUFFIXES = ["_488.tif", "_Dapi-ex_Green-em.tif"]   # summed; the first one names the image
# Treatment table (CSV, relative to the experiment folder). Column IMAGE_COLUMN holds the image
# name without "_XY..." (e.g. neurite_vLH92_MN_60x-1xzoom_0000). Optional column "replicate".
TREATMENT_FILE = "analysis/fileid_treatment.csv"
IMAGE_COLUMN = "image"
TREATMENT_COLUMN = "treatment"
TIME_COLUMN = "treat_time"                    # "" if there is none
RECOVERY_COLUMN = "recovery time"             # "" if there is none
DEFAULT_REPLICATE = "A"

OUTPUT_SUBFOLDER = "analysis/neurite_beading" # per experiment, created inside the experiment folder
POOLED_OUTPUT = ""                            # folder for the combined results of all experiments (batch)
                                              # ("" = <first BATCH_ROOT>/neurite_beading_pooled, or
                                              #  <first experiment>/OUTPUT_SUBFOLDER/pooled)

# --- Controls (for statistics and plot order) -------------------------------
# Treatment family -> its control label. A treatment's family is the first key found in its label.
# "glu" is read as "gluc".
CONTROLS = {
    "2DG": "0mM 2DG",
    "MMS": "0mM MMS",
    "gluc": "2.5mM gluc",
}
CONTROL_FALLBACK = {"0mM 2DG": "2.5mM gluc"}  # used when a family's control is not in the batch
# Control matching: same treat_time and recovery time if that control exists, else same treat_time,
# else all control images of the experiment.
MATCH_TREAT_TIME = True

# --- Segmentation (as in neurite_atp.py) -------------------------------------
BG_RADIUS_PX = 15                             # grey-opening disk radius for the background map
SMOOTH_SIGMA = 1.0                            # Gaussian blur before the ridge filter
RIDGE_SIGMAS = (1, 2, 3)                      # Sato filter scales (px); neurites are ~2-8 px wide
RIDGE_THRESHOLD = "li"                        # "li" (keeps dim neurites), "otsu" (bright ones only) or a number
MIN_SNR = 3.0                                 # mask pixels need signal > MIN_SNR x background noise SD
THICK_RADIUS_PX = 15                          # structures surviving an opening with this disk are blobs
                                              # (diameter > ~5.5 um: cell bodies, clumps; beads are smaller)
MIN_OBJECT_LENGTH_PX = 40                     # drop mask pieces whose skeleton is shorter
MIN_SEGMENT_LENGTH_PX = 20                    # shortest neurite segment kept (0.18 um/px at 60x)
SEGMENT_MAX_DIST_PX = 12                      # mask pixels further than this from a segment are unused
                                              # (large enough to include the body of a bead)

# --- Bead detection ----------------------------------------------------------
BEAD_SIGMA = 1.0                              # Gaussian blur (px) of the summed signal for bead detection
BEAD_FOLD = 2.0                               # bead peak >= BEAD_FOLD x the shaft around it
BEAD_ANNULUS_PX = (8, 22)                     # shaft = median along the skeleton in this ring (px) around
                                              # the peak (inner radius > bead radius; 22 px ~ 4 um)
BEAD_MIN_SHAFT_PX = 8                         # need at least this many skeleton pixels in the ring
BEAD_PEAK_SEP_PX = 3                          # peaks closer than this are one candidate
BEAD_MIN_AREA_PX = 6                          # smallest bead (px), outline = half-way between shaft and peak
BEAD_MAX_LENGTH_UM = 5.0                      # longest bead (major axis); longer = thick stretch, not a bead
BEAD_MIN_SNR = 5.0                            # bead peak must be this many noise SDs above the shaft
JUNCTION_EXCLUDE_PX = 6                       # a bright spot this close to a crossing / branch point is
                                              # taken as overlapping neurites unless it is wide (below);
                                              # a crossing is where >= 3 kept segments meet (the short
                                              # spurs a bead adds to the skeleton do not count)
JUNCTION_MERGE_PX = 5                         # skeleton nodes this close are merged into one crossing
BEAD_MIN_WIDTH_AT_CROSSING_UM = 1.5           # wider spots at a crossing still count as beads
PIXEL_SIZE_UM = 0.183333                      # used if the TIFF has no OME pixel size

# --- Plots -------------------------------------------------------------------
PLOT_METRICS = ["beads_per_100um", "bead_area_fraction", "intensity_cv"]
MIN_PLOT_SEGMENT_UM = 10                      # segment-level plot: only segments at least this long
TREATMENT_ORDER = []                          # plot order of condition labels; [] = automatic
SAVE_OVERLAYS = True
FIG_DPI = 200

# =============================================================================
# END OF PARAMETERS
# =============================================================================

import os
import re
import warnings

import numpy as np
import pandas as pd
import tifffile
from scipy import ndimage as ndi
from scipy.stats import mannwhitneyu
from skimage import filters, morphology, measure, segmentation

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

warnings.filterwarnings("ignore", category=RuntimeWarning)
warnings.filterwarnings("ignore", category=FutureWarning)

METRIC_LABELS = {
    "beads_per_100um": "beads per 100 µm neurite",
    "bead_area_fraction": "fraction of neurite area in beads",
    "intensity_cv": "intensity CV along neurite",
}


def here(path):
    return path if os.path.isabs(path) else os.path.join(os.path.dirname(os.path.abspath(__file__)), path)


def apply_aliases(text):
    return re.sub(r"\bglu\b", "gluc", str(text), flags=re.I).strip()


def norm_key(text):
    return re.sub(r"\s+", "", apply_aliases(text)).lower()


def pixel_size_um(path):
    try:
        with tifffile.TiffFile(path) as t:
            m = re.search(r'PhysicalSizeX="([0-9.]+)"', t.ome_metadata or "")
            return float(m.group(1)) if m else PIXEL_SIZE_UM
    except Exception:
        return PIXEL_SIZE_UM


# --- image processing --------------------------------------------------------
def background_map(img):
    opened = ndi.grey_opening(img, footprint=morphology.disk(BG_RADIUS_PX))
    return ndi.gaussian_filter(opened, BG_RADIUS_PX)


def noise_sd(sig, mask):
    v = sig[~mask]
    return 1.4826 * np.median(np.abs(v - np.median(v)))


def neurite_mask(total_sn):
    total = ndi.gaussian_filter(total_sn, SMOOTH_SIGMA)
    ridge = filters.sato(total, sigmas=RIDGE_SIGMAS, black_ridges=False)
    auto = {"li": filters.threshold_li, "otsu": filters.threshold_otsu}
    thr = auto[RIDGE_THRESHOLD](ridge) if RIDGE_THRESHOLD in auto else float(RIDGE_THRESHOLD)
    mask = ridge > thr
    sd = noise_sd(total, mask)
    mask &= total > MIN_SNR * sd
    # beads are brighter and wider than the ridge response; add bright pixels touching the mask
    mask |= ndi.binary_dilation(mask, morphology.disk(2)) & (total > 2 * MIN_SNR * sd)
    mask = morphology.remove_small_holes(ndi.binary_closing(mask, morphology.disk(1)), 20)
    thick = ndi.binary_opening(mask, morphology.disk(THICK_RADIUS_PX))
    thick = ndi.binary_dilation(thick, morphology.disk(2))
    mask &= ~thick
    lab = measure.label(mask)
    skel = morphology.skeletonize(mask)
    lengths = np.bincount(lab[skel], minlength=lab.max() + 1)
    keep = lengths >= MIN_OBJECT_LENGTH_PX
    keep[0] = False
    return keep[lab], thick


def segments(mask):
    skel = morphology.skeletonize(mask)
    nb = ndi.convolve(skel.astype(int), np.ones((3, 3), int), mode="constant") - 1
    junction = skel & (nb >= 3)
    junction_zone = ndi.binary_dilation(junction, morphology.disk(2))
    branches = skel & ~junction_zone
    lab = measure.label(branches, connectivity=2)
    lengths = np.bincount(lab.ravel(), minlength=lab.max() + 1)
    lengths[0] = 0
    lab[lengths[lab] < MIN_SEGMENT_LENGTH_PX] = 0
    dist, (iy, ix) = ndi.distance_transform_edt(lab == 0, return_indices=True)
    seg = lab[iy, ix]
    seg[~mask | (dist > SEGMENT_MAX_DIST_PX)] = 0
    # real crossings / branch points: junction zones touched by >= 3 kept segments
    # (nodes closer than 2 x JUNCTION_MERGE_PX are one crossing: an X skeletonises into two nodes)
    zlab = measure.label(ndi.binary_dilation(junction, morphology.disk(JUNCTION_MERGE_PX)))
    grown = segmentation.expand_labels(zlab, 2)
    both = (grown > 0) & (lab > 0)
    pairs = np.unique(np.stack([grown[both], lab[both]]), axis=1)
    n_touch = np.bincount(pairs[0], minlength=zlab.max() + 1)
    crossing = (n_touch >= 3)[zlab] & (zlab > 0)
    return seg, lab, lengths, crossing, skel


def find_beads(smooth, mask, skel, seg, crossing, sd, um):
    """Beads = local maxima of the signal that are BEAD_FOLD x brighter than the neurite shaft
    around them (median along the skeleton in an annulus), grown to their half-maximum outline."""
    peaks = (smooth == ndi.maximum_filter(smooth, size=2 * BEAD_PEAK_SEP_PX + 1)) & mask
    near_crossing = ndi.distance_transform_edt(~crossing) <= JUNCTION_EXCLUDE_PX
    r_in, r_out = BEAD_ANNULUS_PX
    yy, xx = np.mgrid[-r_out:r_out + 1, -r_out:r_out + 1]
    rr = np.hypot(yy, xx)
    ring = (rr >= r_in) & (rr <= r_out)
    disk = rr <= r_out
    h, w = smooth.shape
    pad = r_out
    sm_p = np.pad(smooth, pad, constant_values=0)
    sk_p = np.pad(skel, pad, constant_values=False)
    mk_p = np.pad(mask, pad, constant_values=False)
    bead_mask = np.zeros(smooth.shape, bool)
    beads = []
    for py, px in sorted(zip(*np.nonzero(peaks)), key=lambda t: -smooth[t]):
        if bead_mask[py, px]:
            continue
        win = (slice(py, py + 2 * pad + 1), slice(px, px + 2 * pad + 1))
        shaft = sm_p[win][ring & sk_p[win]]
        if shaft.size < BEAD_MIN_SHAFT_PX:
            continue
        base = np.median(shaft)
        peak = smooth[py, px]
        if base <= 0 or peak < BEAD_FOLD * base or peak - base < BEAD_MIN_SNR * sd:
            continue
        # outline: connected pixels above half-way between shaft and peak, within r_out
        region = (sm_p[win] >= base + 0.5 * (peak - base)) & mk_p[win] & disk
        lab = measure.label(region)
        region = lab == lab[pad, pad]
        if region[0, :].any() or region[-1, :].any() or region[:, 0].any() or region[:, -1].any():
            continue                                    # runs out of the window: a thick stretch
        props = measure.regionprops(region.astype(int))[0]
        length, width = props.axis_major_length * um, props.axis_minor_length * um
        if props.area < BEAD_MIN_AREA_PX or length > BEAD_MAX_LENGTH_UM:
            continue
        if near_crossing[py, px] and width < BEAD_MIN_WIDTH_AT_CROSSING_UM:
            continue                                    # two neurites overlapping, not a bead
        ys, xs = np.nonzero(region)
        ys, xs = ys + py - pad, xs + px - pad
        ok = (ys >= 0) & (ys < h) & (xs >= 0) & (xs < w)
        bead_mask[ys[ok], xs[ok]] = True
        beads.append({"segment": int(seg[py, px]), "x": int(px), "y": int(py), "area_px": int(props.area),
                      "area_um2": props.area * um * um, "length_um": length, "width_um": width,
                      "peak": peak, "shaft": base, "peak_over_shaft": peak / base})
    return beads, bead_mask


def process_image(paths, overlay_path=None):
    imgs = [tifffile.imread(p).astype(float) for p in paths]
    total = sum(im - background_map(im) for im in imgs)
    raw = sum(imgs)
    um = pixel_size_um(paths[0])
    mask, thick = neurite_mask(total)
    seg, branch_lab, lengths, crossing, skel = segments(mask)
    smooth = ndi.gaussian_filter(total, BEAD_SIGMA)
    sd = noise_sd(smooth, mask | thick)
    beads, bead_mask = find_beads(smooth, mask, skel, seg, crossing, sd, um)
    ids = np.unique(branch_lab[branch_lab > 0])
    baseline = np.zeros(branch_lab.max() + 1)
    if len(ids):
        baseline[ids] = ndi.median(smooth, branch_lab, ids)
    edt = ndi.distance_transform_edt(mask)

    n_beads = np.bincount([b["segment"] for b in beads], minlength=branch_lab.max() + 1)
    rows = []
    for s in np.unique(branch_lab[branch_lab > 0]):
        on_skel = branch_lab == s
        in_seg = seg == s
        prof = smooth[on_skel]
        length_um = lengths[s] * um
        rows.append({
            "segment": int(s), "length_px": int(lengths[s]), "length_um": length_um,
            "area_px": int(in_seg.sum()), "n_beads": int(n_beads[s]),
            "beads_per_100um": 100 * n_beads[s] / length_um,
            "bead_area_fraction": (bead_mask & in_seg).sum() / max(in_seg.sum(), 1),
            "intensity_cv": prof.std() / prof.mean() if prof.mean() > 0 else np.nan,
            "shaft_intensity": baseline[s], "shaft_width_um": 2 * np.median(edt[on_skel]) * um,
            "x": float(np.nonzero(on_skel)[1].mean()), "y": float(np.nonzero(on_skel)[0].mean()),
        })
    total_len = sum(r["length_um"] for r in rows)
    blobs = measure.label(thick & (total > MIN_SNR * sd))
    info = {
        "total_length_um": total_len, "n_segments": len(rows), "n_beads": len(beads),
        "beads_per_100um": 100 * len(beads) / total_len if total_len else np.nan,
        "bead_area_fraction": bead_mask.sum() / mask.sum() if mask.any() else np.nan,
        "intensity_cv": np.average([r["intensity_cv"] for r in rows], weights=[r["length_um"] for r in rows])
                        if rows else np.nan,
        "bead_area_um2_median": np.median([b["area_um2"] for b in beads]) if beads else np.nan,
        "mask_fraction": mask.mean(), "noise_sd": sd,
        "n_large_blobs": int(blobs.max()), "large_blob_fraction": thick.mean(),
    }
    if overlay_path:
        save_overlay(raw, mask, bead_mask, beads, overlay_path, info)
    return info, rows, beads


def save_overlay(raw, mask, bead_mask, beads, path, info):
    fig, ax = plt.subplots(1, 2, figsize=(10, 7.5))
    lo, hi = np.percentile(raw, [1, 99.7])
    for a in ax:
        a.imshow(raw, cmap="gray", vmin=lo, vmax=hi)
        a.axis("off")
    ax[0].contour(mask, levels=[0.5], colors="yellow", linewidths=0.3)
    ax[0].set_title("neurite mask", fontsize=9)
    ax[1].contour(bead_mask, levels=[0.5], colors="red", linewidths=0.5)
    ax[1].set_title(f"{len(beads)} beads, {info['beads_per_100um']:.1f} per 100 µm, "
                    f"{100 * info['bead_area_fraction']:.1f}% of area", fontsize=9)
    fig.suptitle(os.path.basename(path)[:-4], fontsize=9)
    fig.tight_layout()
    fig.savefig(path, dpi=130)
    plt.close(fig)


# --- experiment --------------------------------------------------------------
def load_treatments(folder):
    path = os.path.join(folder, TREATMENT_FILE)
    if not os.path.exists(path):
        print(f"  ! no treatment file {path}; all images get treatment 'unknown'")
        return None
    t = pd.read_csv(path, encoding="utf-8-sig")
    t.columns = [c.strip() for c in t.columns]
    t[IMAGE_COLUMN] = t[IMAGE_COLUMN].astype(str).str.strip()
    return t.drop_duplicates(IMAGE_COLUMN).set_index(IMAGE_COLUMN)


def condition_label(tr, tt, rec):
    parts = [tr] + ([tt] if tt else []) + ([f"rec {rec}"] if rec else [])
    return " | ".join(parts)


def run_experiment(folder):
    name = os.path.basename(os.path.normpath(folder))
    img_dir = os.path.join(folder, IMAGE_SUBFOLDER)
    out = os.path.join(folder, OUTPUT_SUBFOLDER)
    os.makedirs(out, exist_ok=True)
    if SAVE_OVERLAYS:
        os.makedirs(os.path.join(out, "overlays"), exist_ok=True)
    treat = load_treatments(folder)
    first = CHANNEL_SUFFIXES[0]
    files = sorted(f for f in os.listdir(img_dir) if f.startswith(FILE_PREFIX) and f.endswith(first))
    print(f"{name}: {len(files)} images")
    seg_rows, img_rows, bead_rows = [], [], []
    for f in files:
        base = f[: -len(first)]
        paths = [os.path.join(img_dir, base + s) for s in CHANNEL_SUFFIXES]
        missing = [p for p in paths if not os.path.exists(p)]
        if missing:
            print(f"  ! missing {os.path.basename(missing[0])}")
            continue
        stem = re.sub(r"_XY\d+$", "", base)
        meta = {"experiment": name, "image": base, "dish": stem}
        if treat is not None and stem in treat.index:
            row = treat.loc[stem]
            meta["treatment"] = apply_aliases(row[TREATMENT_COLUMN])
            meta["treat_time"] = str(row[TIME_COLUMN]).strip() if TIME_COLUMN and TIME_COLUMN in row else ""
            meta["recovery"] = (str(row[RECOVERY_COLUMN]).strip()
                                if RECOVERY_COLUMN and RECOVERY_COLUMN in row else "")
            meta["replicate"] = str(row["replicate"]) if "replicate" in row else DEFAULT_REPLICATE
        else:
            meta.update({"treatment": "unknown", "treat_time": "", "recovery": "",
                         "replicate": DEFAULT_REPLICATE})
        meta["condition"] = condition_label(meta["treatment"], meta["treat_time"], meta["recovery"])
        ov = os.path.join(out, "overlays", base + ".png") if SAVE_OVERLAYS else None
        info, rows, beads = process_image(paths, ov)
        img_rows.append({**meta, **info})
        seg_rows += [{**meta, **r} for r in rows]
        bead_rows += [{**meta, **b} for b in beads]
        print(f"  {base}: {info['n_segments']} segments, {info['n_beads']} beads, "
              f"{info['beads_per_100um']:.2f} per 100 um")
    segs, imgs, beads = pd.DataFrame(seg_rows), pd.DataFrame(img_rows), pd.DataFrame(bead_rows)
    segs.to_csv(os.path.join(out, "beading_segments.csv"), index=False)
    imgs.to_csv(os.path.join(out, "beading_images.csv"), index=False)
    beads.to_csv(os.path.join(out, "beading_beads.csv"), index=False)
    write_summary_and_plots(segs, imgs, out)
    return segs, imgs, beads


# --- controls, summary and plots ---------------------------------------------
def family(text):
    for fam in CONTROLS:
        if fam.lower() in text.lower():
            return fam
    return None


def hours(text):
    m = re.match(r"\s*([0-9.]+)", str(text))
    return float(m.group(1)) if m else 0.0


def control_pool(e, g):
    """Rows of experiment table e that are the control for condition group g."""
    tr, tt, rec = g["treatment"].iloc[0], g["treat_time"].iloc[0], g["recovery"].iloc[0]
    fam = family(tr)
    if fam is None:
        return None, ""
    ctrl = CONTROLS[fam]
    pool = e[e["treatment"].map(norm_key) == norm_key(ctrl)]
    if pool.empty and ctrl in CONTROL_FALLBACK:
        ctrl = CONTROL_FALLBACK[ctrl]
        pool = e[e["treatment"].map(norm_key) == norm_key(ctrl)]
    if pool.empty:
        return None, ""
    if MATCH_TREAT_TIME:
        both = pool[(pool["treat_time"] == tt) & (pool["recovery"] == rec)]
        same_t = pool[pool["treat_time"] == tt]
        pool = both if len(both) else same_t if len(same_t) else pool
    labels = sorted(pool["condition"].unique())
    return pool, " + ".join(labels)


def condition_order(conds, df):
    if TREATMENT_ORDER:
        return [c for c in TREATMENT_ORDER if c in conds] + [c for c in conds if c not in TREATMENT_ORDER]
    controls = {norm_key(v) for v in list(CONTROLS.values()) + list(CONTROL_FALLBACK.values())}
    first = df.drop_duplicates("condition").set_index("condition")

    def key(c):
        r = first.loc[c]
        fam = family(r["treatment"])
        dose = hours(r["treatment"])
        rank = 0 if norm_key(r["treatment"]) in controls else 1
        return (hours(r["treat_time"]), hours(r["recovery"]), rank,
                -dose if fam == "gluc" else dose, c)
    return sorted(conds, key=key)


def stats_vs_control(df, level):
    rows = []
    for exp, e in df.groupby("experiment"):
        for cond, g in e.groupby("condition"):
            pool, ctrl = control_pool(e, g)
            if pool is None or set(pool["condition"]) == {cond}:
                continue
            for m in PLOT_METRICS:
                x, y = g[m].dropna(), pool[m].dropna()
                p = mannwhitneyu(x, y).pvalue if len(x) > 1 and len(y) > 1 else np.nan
                rows.append({"experiment": exp, "level": level, "metric": m, "condition": cond,
                             "control": ctrl, "n": len(x), "n_control": len(y), "median": x.median(),
                             "median_control": y.median(),
                             "fold_vs_control": x.median() / y.median() if y.median() else np.nan,
                             "mannwhitney_p": p})
    return rows


def write_summary_and_plots(segs, imgs, out):
    if imgs.empty:
        print("  ! no images; nothing to plot")
        return
    rows = []
    for (exp, cond), gi in imgs.groupby(["experiment", "condition"]):
        g = segs[(segs["experiment"] == exp) & (segs["condition"] == cond)]
        r = {"experiment": exp, "condition": cond, "n_images": len(gi), "n_segments": len(g),
             "total_length_um": gi["total_length_um"].sum(), "n_beads": gi["n_beads"].sum()}
        for m in PLOT_METRICS:
            r[f"image_{m}_median"] = gi[m].median()
            r[f"segment_{m}_median"] = g[m].median() if len(g) else np.nan
        r["bead_area_um2_median"] = gi["bead_area_um2_median"].median()
        rows.append(r)
    pd.DataFrame(rows).to_csv(os.path.join(out, "summary_by_treatment.csv"), index=False)
    st = stats_vs_control(imgs, "image") + stats_vs_control(segs, "segment")
    pd.DataFrame(st).to_csv(os.path.join(out, "stats_vs_control.csv"), index=False)
    superplot(imgs, os.path.join(out, "superplot_images.png"), image_medians=False,
              title="Images (all segments of a field pooled)")
    long_segs = segs[segs["length_um"] >= MIN_PLOT_SEGMENT_UM]
    superplot(long_segs, os.path.join(out, "superplot_segments.png"), image_medians=True,
              title=f"Neurite segments ≥ {MIN_PLOT_SEGMENT_UM} µm")


def superplot(df, path, image_medians, title):
    if df.empty:
        return
    conds = condition_order(list(df["condition"].unique()), df)
    groups = sorted(df.groupby(["experiment", "replicate"]).groups)
    markers = "osD^v<>p*h"
    fig, axes = plt.subplots(len(PLOT_METRICS), 1, figsize=(1.4 + 1.0 * len(conds), 3.3 * len(PLOT_METRICS)),
                             squeeze=False)
    for ax, m in zip(axes[:, 0], PLOT_METRICS):
        rng = np.random.default_rng(0)
        for i, c in enumerate(conds):
            g = df[(df["condition"] == c) & df[m].notna()]
            ax.scatter(i + rng.uniform(-0.25, 0.25, len(g)), g[m], s=6, color="0.4", alpha=0.25,
                       linewidths=0, zorder=1)
            ax.boxplot(g[m], positions=[i], widths=0.55, showfliers=False,
                       medianprops={"color": "k"}, boxprops={"color": "k"},
                       whiskerprops={"color": "k"}, capprops={"color": "k"}, zorder=2)
            if image_medians:
                im = g.groupby(["experiment", "image"])[m].median()
                ax.scatter(np.full(len(im), i + 0.33), im, s=14, facecolors="none", edgecolors="tab:blue",
                           linewidths=0.8, zorder=3, label="image median" if i == 0 else None)
            for j, key in enumerate(groups):
                gg = g[(g["experiment"] == key[0]) & (g["replicate"] == key[1])]
                if len(gg):
                    ax.scatter(i, gg[m].median(), s=90, marker=markers[j % len(markers)],
                               color=plt.cm.tab10(j % 10), edgecolors="k", zorder=4,
                               label=f"{key[0]} {key[1]}" if i == 0 else None)
        ax.set_xticks(range(len(conds)))
        ax.set_xticklabels([c.replace(" | ", "\n") for c in conds], fontsize=7)
        ax.set_ylabel(METRIC_LABELS.get(m, m), fontsize=8)
        ax.spines[["top", "right"]].set_visible(False)
        ax.legend(fontsize=6, frameon=False, loc="best")
    axes[0, 0].set_title(f"{title} (n per box: " + ", ".join(str((df['condition'] == c).sum()) for c in conds)
                         + ")", fontsize=9)
    fig.tight_layout()
    fig.savefig(path, dpi=FIG_DPI)
    plt.close(fig)


def has_images(folder):
    d = os.path.join(folder, IMAGE_SUBFOLDER)
    return os.path.isdir(d) and any(f.startswith(FILE_PREFIX) and f.endswith(CHANNEL_SUFFIXES[0])
                                    for f in os.listdir(d))


def find_experiments():
    found = [here(f) for f in EXPERIMENT_FOLDERS]
    for root in BATCH_ROOTS:
        root = here(root)
        base = root.rstrip(os.sep).count(os.sep)
        for d, subs, _ in os.walk(root):
            if d.count(os.sep) - base >= SEARCH_DEPTH:
                subs[:] = []
            subs[:] = [x for x in subs if not x.startswith(".") and x != IMAGE_SUBFOLDER]
            if has_images(d):
                found.append(d)
    out = []
    for f in found:
        name = os.path.basename(os.path.normpath(f))
        if f in out or any(w in name for w in SKIP_FOLDERS):
            continue
        if not has_images(f):
            print(f"! skipped {f}: no {FILE_PREFIX}*{CHANNEL_SUFFIXES[0]} in {IMAGE_SUBFOLDER}/")
            continue
        out.append(f)
    return out


def main():
    folders = find_experiments()
    print(f"{len(folders)} experiment folder(s):")
    for f in folders:
        print("  " + f)
    all_s, all_i, all_b = [], [], []
    for folder in folders:
        s, i, b = run_experiment(folder)
        if len(i):
            all_s.append(s)
            all_i.append(i)
            all_b.append(b)
    if len(all_i) > 1:
        if POOLED_OUTPUT:
            pooled = here(POOLED_OUTPUT)
        elif BATCH_ROOTS:
            pooled = os.path.join(here(BATCH_ROOTS[0]), "neurite_beading_pooled")
        else:
            pooled = os.path.join(folders[0], OUTPUT_SUBFOLDER, "pooled")
        os.makedirs(pooled, exist_ok=True)
        segs, imgs = pd.concat(all_s, ignore_index=True), pd.concat(all_i, ignore_index=True)
        segs.to_csv(os.path.join(pooled, "beading_segments.csv"), index=False)
        imgs.to_csv(os.path.join(pooled, "beading_images.csv"), index=False)
        pd.concat(all_b, ignore_index=True).to_csv(os.path.join(pooled, "beading_beads.csv"), index=False)
        write_summary_and_plots(segs, imgs, pooled)
        print(f"pooled results: {pooled}")


if __name__ == "__main__":
    main()
