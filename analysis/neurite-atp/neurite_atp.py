#!/usr/bin/env python3
"""
ATP/ADP (vLH92 / PercevalHR-type sensor, 488/405 excitation ratio) in neurites.

HOW TO USE
  1. Edit the PARAMETERS block below (experiment folders, file names, settings).
  2. Run the file:  python neurite_atp.py   (or press Run in your editor).
  Requires: numpy, pandas, scipy, scikit-image, tifffile, matplotlib.

WHAT IT DOES (per image = one XY field)
  1. Reads the 488-excitation image (ATP) and the 405-excitation image ("Dapi-ex_Green-em", ADP).
  2. Background: a smooth background map per channel (grey opening with a disk wider than any
     neurite, then Gaussian smoothing), subtracted from each channel.
  3. Neurite mask: ridge filter (Sato) on the summed, background-subtracted signal, thresholded,
     then thick blobs (somata, clumps, debris wider than THICK_RADIUS_PX) and small pieces
     (skeleton shorter than MIN_OBJECT_LENGTH_PX) are removed.
  4. Neurite segments: the mask is skeletonised and cut at branch points and crossings; each piece
     of skeleton at least MIN_SEGMENT_LENGTH_PX long is one "segment" (a stretch of one neurite
     between crossings). Mask pixels go to their nearest segment.
  5. Ratio: pixels with enough signal in both channels (MIN_SNR x background noise) are used.
     ratio = sum(488 - bg) / sum(405 - bg)  (intensity-weighted, robust for dim pixels);
     the median of the per-pixel ratio is also saved.
  6. Normalisation: each value is divided by the median of its batch's control (per experiment
     folder, see ATP_CONTROLS / MATCH_TREAT_TIME).

  Neurites in these fields form a dense network with no somata, so single cells cannot be
  followed; the per-"cell" unit here is the neurite segment.

OUTPUTS (in <experiment>/OUTPUT_SUBFOLDER; POOLED_OUTPUT gets the combined tables and plots)
  neurite_segments.csv     one row per segment: image, treatment, length, intensities, ratios
  neurite_images.csv       one row per image: whole-mask ratio, n segments, background, coverage
  summary_by_treatment.csv per treatment: n images / segments, medians of raw and normalised ratio
  superplot_segments.png   normalised ratio per segment (faint dots), boxplots, per-image medians,
                           big dots = experiment + replicate medians
  superplot_images.png     same, one dot per image (whole-mask ratio)
  ratio_vs_intensity.png   QC: segment ratio vs 405 intensity (dim segments should not drift)
  overlays/<image>.png     QC: mask outline and segments coloured by ratio
"""

# =============================================================================
# PARAMETERS — edit these
# =============================================================================

# --- Which data to analyse ---------------------------------------------------
# One folder per experiment. Images are read from <folder>/IMAGE_SUBFOLDER.
EXPERIMENT_FOLDERS = [
    "/Volumes/holtl02lab/holtl02labspace/Holt_Lab_Members/Yi-Shuan/microscope/2026/10-12/20261002_MN_92",
]
IMAGE_SUBFOLDER = "tif"
FILE_PREFIX = "neurite_"                      # only images whose name starts with this
ATP_SUFFIX = "_488.tif"                       # 488 nm excitation (ATP-bound)
ADP_SUFFIX = "_Dapi-ex_Green-em.tif"          # 405 nm excitation (ADP-bound)
# Treatment table (CSV, relative to the experiment folder). Column IMAGE_COLUMN holds the image
# name without "_XY..." (e.g. neurite_vLH92_MN_60x-1xzoom_0000). Optional column "replicate".
TREATMENT_FILE = "analysis/fileid_treatment.csv"
IMAGE_COLUMN = "image"
TREATMENT_COLUMN = "treatment"
TIME_COLUMN = "treat_time"                    # "" if there is none
DEFAULT_REPLICATE = "A"

OUTPUT_SUBFOLDER = "analysis/neurite_atp"     # per experiment, created inside the experiment folder
POOLED_OUTPUT = ""                            # folder for the combined results of all experiments
                                              # ("" = the first experiment's output folder)

# --- Normalisation -----------------------------------------------------------
NORMALIZE = True
# Treatment family -> its control label. A treatment's family is the first key found as a word in
# its label. "glu" is read as "gluc".
ATP_CONTROLS = {
    "2DG": "0mM 2DG",
    "MMS": "0mM MMS",
    "gluc": "2.5mM gluc",
}
# If a family's control is not in the batch, use this one instead (2DG is added to 2.5 mM glucose
# medium, so "2.5mM gluc" is the 0 mM 2DG control).
CONTROL_FALLBACK = {"0mM 2DG": "2.5mM gluc"}
MATCH_TREAT_TIME = True                       # use the control with the same treat_time when there is one
NORMALIZE_LEVEL = "segment"                   # control median taken over "segment" or "image" values

# --- Segmentation ------------------------------------------------------------
BG_RADIUS_PX = 15                             # grey-opening disk radius for the background map
SMOOTH_SIGMA = 1.0                            # Gaussian blur before the ridge filter
RIDGE_SIGMAS = (1, 2, 3)                      # Sato filter scales (px); neurites are ~2-8 px wide
RIDGE_THRESHOLD = "li"                        # "li" (keeps dim neurites), "otsu" (bright ones only) or a number
MIN_SNR = 3.0                                 # pixel needs signal > MIN_SNR x noise SD in both channels
THICK_RADIUS_PX = 6                           # structures surviving an opening with this disk are blobs
MIN_OBJECT_LENGTH_PX = 40                     # drop mask pieces whose skeleton is shorter
MIN_SEGMENT_LENGTH_PX = 20                    # shortest neurite segment kept (0.18 um/px at 60x)
SEGMENT_MAX_DIST_PX = 6                       # mask pixels further than this from a segment are unused
MIN_SEGMENT_PIXELS = 30                       # segments with fewer usable pixels are dropped

# --- Plots -------------------------------------------------------------------
RATIO_COLUMN = "ratio_sum"                    # "ratio_sum" (sum 488 / sum 405) or "ratio_median"
TREATMENT_ORDER = []                          # plot order of "treatment | treat_time" labels; [] = automatic
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
from skimage import filters, morphology, measure

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

warnings.filterwarnings("ignore", category=RuntimeWarning)


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
            return float(m.group(1)) if m else np.nan
    except Exception:
        return np.nan


# --- image processing --------------------------------------------------------
def background_map(img):
    opened = ndi.grey_opening(img, footprint=morphology.disk(BG_RADIUS_PX))
    return ndi.gaussian_filter(opened, BG_RADIUS_PX)


def noise_sd(sig, mask):
    v = sig[~mask]
    return 1.4826 * np.median(np.abs(v - np.median(v)))


def neurite_mask(atp_sn, adp_sn):
    total = ndi.gaussian_filter(atp_sn + adp_sn, SMOOTH_SIGMA)
    ridge = filters.sato(total, sigmas=RIDGE_SIGMAS, black_ridges=False)
    auto = {"li": filters.threshold_li, "otsu": filters.threshold_otsu}
    thr = auto[RIDGE_THRESHOLD](ridge) if RIDGE_THRESHOLD in auto else float(RIDGE_THRESHOLD)
    mask = ridge > thr
    # signal must also be above background noise (ridges also fire on noise)
    sd = noise_sd(total, mask)
    mask &= total > MIN_SNR * sd
    mask = ndi.binary_closing(mask, morphology.disk(1))
    # thick blobs: somata, clumps, debris
    thick = ndi.binary_opening(mask, morphology.disk(THICK_RADIUS_PX))
    thick = ndi.binary_dilation(thick, morphology.disk(2))
    mask &= ~thick
    # drop short pieces
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
    junction = ndi.binary_dilation(junction, morphology.disk(2))
    branches = skel & ~junction
    lab = measure.label(branches, connectivity=2)
    lengths = np.bincount(lab.ravel(), minlength=lab.max() + 1)
    lengths[0] = 0
    lab[lengths[lab] < MIN_SEGMENT_LENGTH_PX] = 0
    # each mask pixel -> nearest kept branch pixel
    dist, (iy, ix) = ndi.distance_transform_edt(lab == 0, return_indices=True)
    seg = lab[iy, ix]
    seg[~mask | (dist > SEGMENT_MAX_DIST_PX)] = 0
    return seg, lab, lengths


def ratio_stats(a, b):
    return {
        "atp_sum": a.sum(), "adp_sum": b.sum(),
        "atp_mean": a.mean(), "adp_mean": b.mean(),
        "ratio_sum": a.sum() / b.sum(),
        "ratio_median": np.median(a / b),
        "n_pixels": a.size,
    }


def process_image(atp_path, adp_path, overlay_path=None):
    atp = tifffile.imread(atp_path).astype(float)
    adp = tifffile.imread(adp_path).astype(float)
    bg_a, bg_b = background_map(atp), background_map(adp)
    atp_sn, adp_sn = atp - bg_a, adp - bg_b
    mask, thick = neurite_mask(atp_sn, adp_sn)
    sd_a, sd_b = noise_sd(atp_sn, mask | thick), noise_sd(adp_sn, mask | thick)
    usable = mask & (atp_sn > MIN_SNR * sd_a) & (adp_sn > MIN_SNR * sd_b)
    seg, branch_lab, lengths = segments(mask)

    img = ratio_stats(atp_sn[usable], adp_sn[usable]) if usable.sum() >= MIN_SEGMENT_PIXELS else None
    info = {
        "bg_488": np.median(bg_a), "bg_405": np.median(bg_b),
        "noise_488": sd_a, "noise_405": sd_b,
        "mask_fraction": mask.mean(), "usable_pixels": int(usable.sum()),
        "saturated_pixels": int(((atp >= 4095) | (adp >= 4095))[mask].sum()),
    }
    rows = []
    for s in np.unique(seg[usable]):
        if s == 0:
            continue
        sel = usable & (seg == s)
        if sel.sum() < MIN_SEGMENT_PIXELS:
            continue
        ys, xs = np.nonzero(branch_lab == s)
        r = ratio_stats(atp_sn[sel], adp_sn[sel])
        r.update({"segment": int(s), "length_px": int(lengths[s]),
                  "x": float(xs.mean()), "y": float(ys.mean())})
        rows.append(r)

    if overlay_path:
        save_overlay(atp, mask, seg, rows, overlay_path)
    return img, info, rows


def save_overlay(atp, mask, seg, rows, path):
    fig, ax = plt.subplots(1, 2, figsize=(10, 7.5))
    lo, hi = np.percentile(atp, [1, 99.5])
    for a in ax:
        a.imshow(atp, cmap="gray", vmin=lo, vmax=hi)
        a.axis("off")
    ax[0].contour(mask, levels=[0.5], colors="yellow", linewidths=0.4)
    ax[0].set_title("neurite mask", fontsize=9)
    if rows:
        vals = {r["segment"]: r[RATIO_COLUMN] for r in rows}
        rmap = np.full(seg.shape, np.nan)
        for s, v in vals.items():
            rmap[seg == s] = v
        v = np.array(list(vals.values()))
        im = ax[1].imshow(rmap, cmap="viridis", vmin=np.percentile(v, 5), vmax=np.percentile(v, 95))
        fig.colorbar(im, ax=ax[1], fraction=0.04, label="488/405 (raw)")
    ax[1].set_title(f"segments ({len(rows)}) coloured by ratio", fontsize=9)
    fig.suptitle(os.path.basename(path)[:-4], fontsize=9)
    fig.tight_layout()
    fig.savefig(path, dpi=110)
    plt.close(fig)


# --- experiment --------------------------------------------------------------
def load_treatments(folder):
    path = os.path.join(folder, TREATMENT_FILE)
    if not os.path.exists(path):
        print(f"  ! no treatment file {path}; all images get treatment 'unknown'")
        return None
    t = pd.read_csv(path, encoding="utf-8-sig")
    t[IMAGE_COLUMN] = t[IMAGE_COLUMN].astype(str).str.strip()
    return t.set_index(IMAGE_COLUMN)


def run_experiment(folder):
    name = os.path.basename(os.path.normpath(folder))
    img_dir = os.path.join(folder, IMAGE_SUBFOLDER)
    out = os.path.join(folder, OUTPUT_SUBFOLDER)
    os.makedirs(out, exist_ok=True)
    if SAVE_OVERLAYS:
        os.makedirs(os.path.join(out, "overlays"), exist_ok=True)
    treat = load_treatments(folder)
    files = sorted(f for f in os.listdir(img_dir) if f.startswith(FILE_PREFIX) and f.endswith(ATP_SUFFIX))
    print(f"{name}: {len(files)} images")
    seg_rows, img_rows = [], []
    for f in files:
        base = f[: -len(ATP_SUFFIX)]
        adp = os.path.join(img_dir, base + ADP_SUFFIX)
        if not os.path.exists(adp):
            print(f"  ! no 405 image for {f}")
            continue
        stem = re.sub(r"_XY\d+$", "", base)
        meta = {"experiment": name, "image": base, "dish": stem}
        if treat is not None and stem in treat.index:
            row = treat.loc[stem]
            meta["treatment"] = apply_aliases(row[TREATMENT_COLUMN])
            meta["treat_time"] = str(row[TIME_COLUMN]) if TIME_COLUMN and TIME_COLUMN in row else ""
            meta["replicate"] = str(row["replicate"]) if "replicate" in row else DEFAULT_REPLICATE
        else:
            meta.update({"treatment": "unknown", "treat_time": "", "replicate": DEFAULT_REPLICATE})
        ov = os.path.join(out, "overlays", base + ".png") if SAVE_OVERLAYS else None
        img, info, rows = process_image(os.path.join(img_dir, f), adp, ov)
        um = pixel_size_um(os.path.join(img_dir, f))
        img_rows.append({**meta, **info, **(img or {}), "n_segments": len(rows)})
        for r in rows:
            r["length_um"] = r["length_px"] * um
            seg_rows.append({**meta, **r})
        print(f"  {base}: {len(rows)} segments, ratio {img['ratio_sum'] if img else np.nan:.3f}")
    segs, imgs = pd.DataFrame(seg_rows), pd.DataFrame(img_rows)
    for df in (segs, imgs):
        df["condition"] = df["treatment"] + (" | " + df["treat_time"]).where(df["treat_time"] != "", "")
    segs, imgs = normalize(segs, imgs)
    segs.to_csv(os.path.join(out, "neurite_segments.csv"), index=False)
    imgs.to_csv(os.path.join(out, "neurite_images.csv"), index=False)
    write_summary_and_plots(segs, imgs, out)
    return segs, imgs


# --- normalisation -----------------------------------------------------------
def family(text):
    for fam in ATP_CONTROLS:
        if re.search(rf"\b{re.escape(fam)}\b", text, flags=re.I) or fam.lower() in text.lower():
            return fam
    return None


def normalize(segs, imgs):
    ref = segs if NORMALIZE_LEVEL == "segment" else imgs
    for df in (segs, imgs):
        df["ratio_norm"] = np.nan
        df["control_used"] = ""
    if not NORMALIZE:
        return segs, imgs
    for exp in imgs["experiment"].unique():
        r = ref[ref["experiment"] == exp]
        for cond, g in imgs[imgs["experiment"] == exp].groupby("condition"):
            tr, tt = g["treatment"].iloc[0], g["treat_time"].iloc[0]
            fam = family(tr)
            if fam is None:
                print(f"  ! {exp}: '{tr}' has no family in ATP_CONTROLS; not normalised")
                continue
            ctrl = ATP_CONTROLS[fam]
            pool = r[r["treatment"].map(norm_key) == norm_key(ctrl)]
            if pool.empty and ctrl in CONTROL_FALLBACK:
                ctrl = CONTROL_FALLBACK[ctrl]
                pool = r[r["treatment"].map(norm_key) == norm_key(ctrl)]
            if pool.empty:
                print(f"  ! {exp}: no control '{ctrl}' for '{tr}'; not normalised")
                continue
            label = ctrl
            if MATCH_TREAT_TIME and (pool["treat_time"] == tt).any():
                pool = pool[pool["treat_time"] == tt]
                label = f"{ctrl} | {tt}"
            med = pool[RATIO_COLUMN].median()
            for df in (segs, imgs):
                sel = (df["experiment"] == exp) & (df["condition"] == cond)
                df.loc[sel, "ratio_norm"] = df.loc[sel, RATIO_COLUMN] / med
                df.loc[sel, "control_used"] = label
    return segs, imgs


# --- summary and plots ------------------------------------------------------
def condition_order(conds):
    if TREATMENT_ORDER:
        return [c for c in TREATMENT_ORDER if c in conds] + [c for c in conds if c not in TREATMENT_ORDER]

    # by treat_time; within a time: the control, then gluc (high to low), then 2DG / MMS (low to high)
    controls = {norm_key(v) for v in list(ATP_CONTROLS.values()) + list(CONTROL_FALLBACK.values())}

    def key(c):
        tr, _, tt = c.partition(" | ")
        fam = family(tr)
        dose = re.match(r"\s*([0-9.]+)", tr)
        dose = float(dose.group(1)) if dose else 0
        hours = re.match(r"\s*([0-9.]+)", tt)
        rank = 0 if norm_key(tr) in controls else 1 if fam == "gluc" else 2
        return (float(hours.group(1)) if hours else 0, rank, -dose if fam == "gluc" else dose, c)
    return sorted(conds, key=key)


def stats_vs_control(df, level):
    """Mann-Whitney U of each condition against its control, within each experiment."""
    rows = []
    for exp, e in df.groupby("experiment"):
        for cond, g in e.groupby("condition"):
            ctrl = g["control_used"].iloc[0]
            if not ctrl or ctrl == cond:
                continue
            tr, _, tt = ctrl.partition(" | ")
            c = e[(e["treatment"].map(norm_key) == norm_key(tr)) & ((e["treat_time"] == tt) if tt else True)]
            if c.empty:
                continue
            x, y = g[RATIO_COLUMN].dropna(), c[RATIO_COLUMN].dropna()
            p = mannwhitneyu(x, y).pvalue if len(x) > 1 and len(y) > 1 else np.nan
            rows.append({"experiment": exp, "level": level, "condition": cond, "control": ctrl,
                         "n": len(x), "n_control": len(y), "median": x.median(),
                         "median_control": y.median(), "fold_vs_control": x.median() / y.median(),
                         "mannwhitney_p": p})
    return rows


def write_summary_and_plots(segs, imgs, out):
    if segs.empty:
        print("  ! no segments; nothing to plot")
        return
    val = "ratio_norm" if NORMALIZE and segs["ratio_norm"].notna().any() else RATIO_COLUMN
    rows = []
    for (exp, cond), g in segs.groupby(["experiment", "condition"]):
        gi = imgs[(imgs["experiment"] == exp) & (imgs["condition"] == cond)]
        rows.append({
            "experiment": exp, "condition": cond, "control_used": g["control_used"].iloc[0],
            "n_images": gi["image"].nunique(), "n_segments": len(g),
            "segment_ratio_median": g[RATIO_COLUMN].median(),
            "segment_ratio_norm_median": g["ratio_norm"].median(),
            "image_ratio_median": gi[RATIO_COLUMN].median(),
            "image_ratio_norm_median": gi["ratio_norm"].median(),
            "image_ratio_norm_sd": gi["ratio_norm"].std(),
            "segment_length_um_median": g["length_um"].median(),
        })
    pd.DataFrame(rows).to_csv(os.path.join(out, "summary_by_treatment.csv"), index=False)
    st = stats_vs_control(imgs, "image") + stats_vs_control(segs, "segment")
    pd.DataFrame(st).to_csv(os.path.join(out, "stats_vs_control.csv"), index=False)
    ylab = "ATP/ADP (488/405), norm. to control" if val == "ratio_norm" else "ATP/ADP (488/405)"
    superplot(segs, val, ylab, os.path.join(out, "superplot_segments.png"), image_medians=True,
              title="Neurite segments")
    superplot(imgs, val, ylab, os.path.join(out, "superplot_images.png"), image_medians=False,
              title="Images (whole neurite mask)")
    qc_ratio_vs_intensity(segs, os.path.join(out, "ratio_vs_intensity.png"))


def superplot(df, val, ylab, path, image_medians, title):
    df = df[df[val].notna()]
    conds = condition_order(list(df["condition"].unique()))
    groups = sorted(df.groupby(["experiment", "replicate"]).groups)
    markers = "osD^v<>p*h"
    rng = np.random.default_rng(0)
    fig, ax = plt.subplots(figsize=(1.2 + 1.1 * len(conds), 4.2))
    for i, c in enumerate(conds):
        g = df[df["condition"] == c]
        ax.scatter(i + rng.uniform(-0.25, 0.25, len(g)), g[val], s=6, color="0.4", alpha=0.25,
                   linewidths=0, zorder=1)
        ax.boxplot(g[val], positions=[i], widths=0.55, showfliers=False,
                   medianprops={"color": "k"}, boxprops={"color": "k"},
                   whiskerprops={"color": "k"}, capprops={"color": "k"}, zorder=2)
        if image_medians:
            im = g.groupby("image")[val].median()
            ax.scatter(np.full(len(im), i + 0.33), im, s=14, facecolors="none", edgecolors="tab:blue",
                       linewidths=0.8, zorder=3, label="image median" if i == 0 else None)
        for j, key in enumerate(groups):
            gg = g[(g["experiment"] == key[0]) & (g["replicate"] == key[1])]
            if len(gg):
                ax.scatter(i, gg[val].median(), s=90, marker=markers[j % len(markers)],
                           color=plt.cm.tab10(j % 10), edgecolors="k", zorder=4,
                           label=f"{key[0]} {key[1]}" if i == 0 else None)
    if val == "ratio_norm":
        ax.axhline(1, color="0.6", lw=0.8, ls="--", zorder=0)
    ax.set_xticks(range(len(conds)))
    ax.set_xticklabels([c.replace(" | ", "\n") for c in conds], fontsize=8)
    ax.set_ylabel(ylab)
    ax.set_title(f"{title} (n per box: " + ", ".join(str((df['condition'] == c).sum()) for c in conds) + ")",
                 fontsize=9)
    ax.spines[["top", "right"]].set_visible(False)
    ax.legend(fontsize=7, frameon=False, loc="best")
    fig.tight_layout()
    fig.savefig(path, dpi=FIG_DPI)
    plt.close(fig)


def qc_ratio_vs_intensity(segs, path):
    conds = condition_order(list(segs["condition"].unique()))
    fig, ax = plt.subplots(figsize=(6, 4.2))
    for i, c in enumerate(conds):
        g = segs[segs["condition"] == c]
        ax.scatter(g["adp_mean"], g[RATIO_COLUMN], s=5, alpha=0.4, color=plt.cm.tab10(i % 10), label=c)
    ax.set_xscale("log")
    ax.set_xlabel("405 intensity, background-subtracted (segment mean)")
    ax.set_ylabel("488/405 (raw)")
    ax.legend(fontsize=6, frameon=False, markerscale=2)
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    fig.savefig(path, dpi=FIG_DPI)
    plt.close(fig)


def main():
    all_s, all_i = [], []
    for folder in EXPERIMENT_FOLDERS:
        s, i = run_experiment(here(folder))
        all_s.append(s)
        all_i.append(i)
    if len(EXPERIMENT_FOLDERS) > 1:
        pooled = here(POOLED_OUTPUT) if POOLED_OUTPUT else os.path.join(here(EXPERIMENT_FOLDERS[0]),
                                                                          OUTPUT_SUBFOLDER, "pooled")
        os.makedirs(pooled, exist_ok=True)
        segs, imgs = pd.concat(all_s, ignore_index=True), pd.concat(all_i, ignore_index=True)
        segs.to_csv(os.path.join(pooled, "neurite_segments.csv"), index=False)
        imgs.to_csv(os.path.join(pooled, "neurite_images.csv"), index=False)
        write_summary_and_plots(segs, imgs, pooled)
        print(f"pooled results: {pooled}")


if __name__ == "__main__":
    main()
