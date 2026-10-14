#!/usr/bin/env python3
"""
Nuclear transport (NCT) analysis of photoactivatable reporter time-lapses.

HOW TO USE
  1. Edit the PARAMETERS block below (file names, folders, treatments, settings).
  2. Run the file:  python nct_analysis.py   (or press Run in your editor).
  Requires: numpy, pandas, scipy, matplotlib.

INPUT FILES (CSV) in each data folder
  RESULTS_FILE  Fiji/ImageJ "Measure" results for every nucleus ROI in every frame
                (Label like "<image>_561.tif:<slice> - ID <roi>:<slice>")
  BG_FILE       Background intensity per image and frame (Label like "<image>_561.tif", Slice 0-based)
  META_FILE     Image -> condition table (image, reporter, treatment, time, pyruvate, replicate;
                image = name without the _XY<n> field suffix). An image may appear on several rows
                with different treatments; its cells are then used in each of those treatments.

PIPELINE
  Fit stage (per cell, independent of treatment; cached in <output>/_cache/):
  1. Parse labels and attach per-frame background; reporter per image comes from META_FILE.
  2. Keep ROIs tracked through every frame with area >= MIN_AREA_PX.
  3. Normalise: norm(t) = (Mean(t) - bg(t)) / (mean(Mean, frames 1..A) - bg(t)), A = activation frame.
  4. QC: "bright" cells have baseline - bg >= MIN_BASELINE grey levels.
  5. Photobleaching: one rate b per experiment (folder) x reporter, from the pre-activation frames ("preact"),
     fitted across cells ("fit"), or none. Trace is multiplied by exp(b * t).
  6. Per-cell fit from activation (t0 = (A-1) * interval):
       y(t) = y_inf + (1 - y_inf) * exp(-k (t - t0)),  k <= 1/interval
     -> amplitude = y_inf - 1, half-time ln2/k, initial flux = k * amplitude (recommended rate),
        plus a model-free initial slope, plateau-subtracted log fit, gain at a fixed time and AUC.
  Assign stage: attach treatments (and replicates) from META_FILE. REASSIGN_ONLY = True reruns from here.
  Report stage, per plot set: per reporter x treatment summaries, statistics (Kruskal-Wallis;
     Mann-Whitney vs the control treatment, Holm) and figures.
"""

# =============================================================================
# PARAMETERS — edit these
# =============================================================================

# --- Which data to analyse ---------------------------------------------------
# Single run: leave BATCH_FOLDERS empty and set DATA_FOLDER.
# Batch: list one folder per experiment (the folder name is the experiment name). The three input
# files are read from <folder>/<BATCH_DATA_SUBFOLDER> ("" = the folder itself).
DATA_FOLDER = "inputs"
BATCH_FOLDERS = [
    "/Volumes/holtl02lab/holtl02labspace/Holt_Lab_Members/Yi-Shuan/microscope/2026/07-09/20260925_MN_110,111",
    "/Volumes/holtl02lab/holtl02labspace/Holt_Lab_Members/Yi-Shuan/microscope/2026/07-09/20260930_MN_110,111",
    "/Volumes/holtl02lab/holtl02labspace/Holt_Lab_Members/Yi-Shuan/microscope/2026/10-12/20261010_MN_110,111",
]
BATCH_DATA_SUBFOLDER = "analysis"
COMBINE_BATCH = True         # True: also run one pooled analysis of all BATCH_FOLDERS together
POOLED_ONLY = False          # True: run only the pooled analysis (the folders are still read, but their
                             # own outputs are not redone); combine with REASSIGN_ONLY for a fast rerun

# File names inside each folder. Wildcards are allowed (e.g. "*combined_Results*.csv"),
# but each pattern must match exactly one file.
RESULTS_FILE = "combined_Results.csv"
BG_FILE = "mCherry_BGInt_Results.csv"
META_FILE = "fileid_treatment.csv"

# Results go to <data folder>/<OUTPUT_SUBFOLDER>; the pooled batch run goes to POOLED_OUTPUT
# (a relative POOLED_OUTPUT is placed next to this script).
OUTPUT_SUBFOLDER = "nct_results"
POOLED_OUTPUT = "nct_results_pooled"

# --- Rerun only the treatment assignment ---------------------------------------
# True: skip loading the measurements and fitting; reuse the per-cell fits saved by the last full run
# (<output>/_cache/) and only re-read META_FILE, then redo summaries, statistics and plots.
# Use it after fixing treatments or replicates. Changing a reporter, the measurements or any
# fitting/bleaching/QC setting needs a full run (False).
REASSIGN_ONLY = False

# --- Treatments and plot sets --------------------------------------------------
# Text replacements applied to the treatment column when it is read (regular expressions, case-insensitive).
# The default turns "glu" (not followed by a letter) into "gluc": "2.5mM glu" -> "2.5mM gluc";
# "gluc", "glucose" and "glutamine" are left alone.
TREATMENT_REPLACE = {r"glu(?![a-z])": "gluc"}
# Cells are grouped by the full text of the `treatment` column.
# Plot sets: each set gets its own summaries, statistics and plots (in a subfolder named after the set
# when there is more than one set). A set is
#   None                     -> every treatment
#   "text"                   -> every treatment whose name contains "text" (case-insensitive)
#   ["name 1", "name 2", ...] -> exactly these treatments, in this order
PLOT_SETS = {
    "all": None,
    # "2DG": "2DG",
    # "gluc": ["2.5mM gluc", "0.15mM gluc", "0mM gluc+FCCP+2DG"],
    # "MMS": "MMS",
}
# Order of treatments on the x axis (for sets given as None or "text"). Treatments listed here come
# first, in this order; the rest follow, sorted by the first number in their name, then alphabetically.
TREATMENT_ORDER = []
# Control (reference) treatment per plot set, for the statistics and "% of control".
# A set without an entry uses its first treatment.
CONTROLS = {
    # "gluc": "2.5mM gluc",
}

# --- Imaging -----------------------------------------------------------------
INTERVAL_S = 30.0            # seconds between frames
ACTIVATION_FRAME = 4         # last pre-activation frame (1-based); activation time t0 = (frame-1)*interval
DIRECTIONS = {110: "import", 111: "export"}   # reporters not listed are guessed from the data
EXCLUDE_REPORTERS = [185]    # reporters left out of fitting, statistics and plots

# --- Cell selection ----------------------------------------------------------
MIN_AREA_PX = 1000           # minimum ROI area in every frame
MIN_BASELINE = 20            # minimum baseline - background (grey levels) to count as "bright"

# --- Photobleaching ----------------------------------------------------------
BLEACH_METHOD = "preact"     # "preact" (pre-activation frames), "fit" (shared fit) or "none"
BLEACH_FIRST_FRAME = 2       # first frame used by "preact" (frame 1 can carry a one-off drop)

# --- Fitting and QC ----------------------------------------------------------
SLOPE_FRAMES = 3             # frames after activation used for the model-free initial slope
GAIN_TIME_S = 600.0          # time (s) at which gain is reported
MIN_R2 = 0.8                 # well-fit: R^2 of the exponential fit
MIN_AMPLITUDE = 0.05         # responder: |amplitude| in the expected direction
LNZ_WINDOW_S = {"import": 510.0, "export": 210.0}   # plateau-subtracted log fit window after t0
WRONG_DIR_WINDOW_S = 300.0   # wrong-direction QC: window after activation (s) ...
WRONG_DIR_MAX = 0.10         # ... cells whose median signal there moved >= this far the wrong way
                             # (import falling / export rising, as a fraction of baseline) are dropped
                             # and listed in wrong_direction_cells.csv; None turns the QC off

# --- Rate plots (initial flux and model-free slope) ---------------------------
RATE_UNIT = "per_min"        # "per_min": min^-1 (e.g. 0.072); "pct": % of baseline per min (e.g. 7.2)
EXPORT_POSITIVE = True       # plot export rates with the sign reversed, so faster export = larger positive value
RATE_YSCALE = "linear"       # "linear" or "log" (log leaves out cells with a rate <= 0; their count is noted)
                             # CSVs and statistics keep the signed % of baseline per min values

# --- Timecourse plots -------------------------------------------------------
TIMECOURSE_ERRORS = ["sem", "sd"]   # error band(s) to plot; one file per entry ("sem" and/or "sd")

# --- Statistics on plots -----------------------------------------------------
SHOW_STARS = True            # stars from the tests vs control, Holm-adjusted
STATS_UNIT = "auto"          # "cells": Mann-Whitney / Kruskal-Wallis on cells;
                             # "replicate": t-tests / ANOVA on experiment+replicate medians (paired by
                             # experiment+replicate when every unit has the control);
                             # "auto": replicate when every treatment has >= 3 units (pooled runs), else cells
STAR_LEVELS = [(0.001, "***"), (0.01, "**"), (0.05, "*")]

# --- Dot plots (superplots) ---------------------------------------------------
CELL_DOT_COLOR = "#b0b0b0"   # single cells in the background

# --- Export vs import rate scatter -------------------------------------------
RATE_SCATTER_METRICS = ["k_per_min", "flux_pct_per_min", "slope_pct_per_min"]   # one plot each:
                                    # k (well-fit cells), initial flux and model-free slope (bright cells; export sign flipped)
RATE_SCATTER_ERROR = "auto"         # "units": SD of experiment+replicate medians; "cells": SD across cells;
                                    # "auto": units when every treatment has >= 2 units, otherwise cells

# --- Per-cell plots ----------------------------------------------------------
SAVE_CELL_PLOTS = False      # one PNG per cell: trace + fitted exponential (slow for many cells)
CELL_PLOTS_BRIGHT_ONLY = True
CELL_PLOT_DPI = 100

# =============================================================================
# END OF PARAMETERS
# =============================================================================

import glob
import os
import re
import warnings
from types import SimpleNamespace

import numpy as np
import pandas as pd
from scipy import stats
from scipy.optimize import curve_fit

# numpy < 2.0 has trapz but not trapezoid
_trapezoid = getattr(np, "trapezoid", None) or np.trapz

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

warnings.filterwarnings("ignore")

DOSE_COLORS = ["#86b6ef", "#5598e7", "#2a78d6", "#1c5cab", "#104281"]
REPORTER_COLORS = {"import": "#2a78d6", "export": "#c2410c"}
# per-treatment shades of the reporter colour (light -> dark along the treatment order) for the dot plots
SHADES = {"import": ["#6fa5e8", "#4a8fe0", "#2a78d6", "#1c5cab", "#104281"],
          "export": ["#ec9460", "#df6a2c", "#c2410c", "#923109", "#632106"]}


def ramp(colors, n):
    if n <= len(colors):
        return [colors[round(i * (len(colors) - 1) / max(n - 1, 1))] for i in range(n)] if n > 1 else [colors[2]]
    from matplotlib.colors import LinearSegmentedColormap
    cmap = LinearSegmentedColormap.from_list("r", colors)
    return [cmap(i / (n - 1)) for i in range(n)]


def shades(direction, n):
    return ramp(SHADES[direction], n)


def settings():
    return SimpleNamespace(
        interval=INTERVAL_S, activation_frame=ACTIVATION_FRAME, min_area=MIN_AREA_PX,
        min_baseline=MIN_BASELINE, bleach=BLEACH_METHOD, bleach_first_frame=BLEACH_FIRST_FRAME,
        slope_frames=SLOPE_FRAMES, gain_time=GAIN_TIME_S, min_r2=MIN_R2, min_amp=MIN_AMPLITUDE)


def safe_name(s):
    return re.sub(r"[^\w.+-]+", "_", str(s)).strip("_")


# --------------------------------------------------------------------------- loading
def find_file(folder, pattern):
    hits = sorted(glob.glob(os.path.join(folder, pattern)))
    if len(hits) != 1:
        raise FileNotFoundError(f"'{pattern}' matched {len(hits)} files in {folder}: {hits}")
    return hits[0]


def load_meta(path, source):
    """Treatment table: one row per image x treatment (an image may have several treatments)."""
    meta = pd.read_csv(path).rename(columns={"image": "image_base", "time": "treat_time"})
    meta["image_base"] = meta["image_base"].astype(str).str.strip()
    meta["treatment"] = meta["treatment"].astype(str).str.strip()
    for pattern, new in TREATMENT_REPLACE.items():
        meta["treatment"] = meta["treatment"].str.replace(pattern, new, case=False, regex=True)
    if "replicate" not in meta:
        meta["replicate"] = "A"
    meta["replicate"] = meta["replicate"].fillna("A").astype(str).str.strip()
    dup = meta.duplicated(["image_base", "treatment"])
    if dup.any():
        print(f"  Note: {dup.sum()} duplicate image + treatment rows in the meta file were ignored")
        meta = meta[~dup]
    reps = meta.groupby("image_base").reporter.nunique()
    if (reps > 1).any():
        raise ValueError("An image has more than one reporter in the meta file: "
                         + ", ".join(reps[reps > 1].index[:5]))
    meta["source"] = source
    return meta


def load_measurements(results_path, bg_path, meta):
    res = pd.read_csv(results_path)
    res = res.rename(columns={res.columns[0]: "row_id"})
    lab = res["Label"].str.extract(r"^(?P<Image>.+?)_\d+\.tif:(?P<Slice>\d+) - ID (?P<roi>\d+):\d+$")
    if lab.isna().any().any():
        bad = res.loc[lab.isna().any(axis=1), "Label"].head(3).tolist()
        raise ValueError(f"Unparsed labels, e.g. {bad}")
    res = pd.concat([res, lab], axis=1)
    res["Slice"] = res["Slice"].astype(int)
    res["image_base"] = res["Image"].str.replace(r"_XY\d+$", "", regex=True)
    res["position"] = res["Image"].str.extract(r"_XY(\d+)$")[0].astype(int)

    reporter = meta.groupby("image_base").reporter.first()
    res["reporter"] = res["image_base"].map(reporter)
    missing = res.loc[res.reporter.isna(), "image_base"].unique()
    if len(missing):
        print(f"  Note: {len(missing)} image(s) are not in the meta file and were left out, e.g. {missing[0]}")
        res = res[res.reporter.notna()]
    res["reporter"] = res["reporter"].astype(int)
    excluded = res.reporter.isin(EXCLUDE_REPORTERS)
    if excluded.any():
        print(f"  Reporter(s) {sorted(int(r) for r in res[excluded].reporter.unique())} left out (EXCLUDE_REPORTERS)")
        res = res[~excluded]

    bg = pd.read_csv(bg_path)
    bg["Image"] = bg["Label"].str.replace(r"_\d+\.tif$", "", regex=True)
    bg = bg.assign(Slice=bg["Slice"].astype(int) + 1)[["Image", "Slice", "Mean"]].rename(columns={"Mean": "bg_mean"})
    res = res.merge(bg, on=["Image", "Slice"], how="left", validate="many_to_one")
    if res["bg_mean"].isna().any():
        raise ValueError("Background missing for some image/frame pairs")
    return res


def experiment(folder, data_folder=None):
    """(experiment name, folder with the input files)."""
    return os.path.basename(os.path.normpath(folder)), data_folder or folder


def load_folder(folder, data_folder=None):
    source, data_folder = experiment(folder, data_folder)
    meta = load_meta(find_file(data_folder, META_FILE), source)
    raw = load_measurements(find_file(data_folder, RESULTS_FILE), find_file(data_folder, BG_FILE), meta)
    raw["source"] = source
    return raw, meta


def filter_rois(df, args):
    n_frames = df.groupby("Image")["Slice"].transform("max")
    df = df.assign(n_frames=n_frames)
    g = df.groupby(["Image", "roi"])
    keep = (g["Slice"].transform("size") == df["n_frames"]) & (g["Area"].transform("min") >= args.min_area)
    out = df[keep].copy()
    out["cell_key"] = out["reporter"].astype(str) + "_" + out["Image"] + "_" + out["roi"]
    return out, df.groupby(["Image", "roi"]).ngroups, out["cell_key"].nunique()


def normalise(df, args):
    A = args.activation_frame
    df = df.sort_values(["cell_key", "Slice"]).copy()
    df["time"] = (df["Slice"] - 1) * args.interval
    df["Mean_SN"] = df["Mean"] - df["bg_mean"]
    df["Median_SN"] = df["Median"] - df["bg_mean"]
    pre = df[df["Slice"] <= A].groupby("cell_key")
    df["pre_mean"] = df["cell_key"].map(pre["Mean"].mean())
    df["pre_median"] = df["cell_key"].map(pre["Median"].mean())
    df["pre_bg"] = df["cell_key"].map(pre["bg_mean"].mean())
    df["norm_mean"] = df["Mean_SN"] / (df["pre_mean"] - df["bg_mean"])
    df["norm_median"] = df["Median_SN"] / (df["pre_median"] - df["bg_mean"])
    df["baseline_bgsub"] = df["pre_mean"] - df["pre_bg"]
    df["bright"] = df["baseline_bgsub"] >= args.min_baseline
    return df


# --------------------------------------------------------------------------- fitting
def relax(t, yinf, k, t0):
    return yinf + (1 - yinf) * np.exp(-k * (t - t0))


def fit_relax(t, y, t0, kmax):
    m = t >= t0
    try:
        p, _ = curve_fit(lambda tt, a, b: relax(tt, a, b, t0), t[m], y[m],
                         p0=(np.mean(y[-3:]), 0.01), bounds=([0, 1e-5], [5, kmax]), maxfev=20000)
        pred = relax(t[m], *p, t0)
        ss = ((y[m] - y[m].mean()) ** 2).sum()
        r2 = 1 - ((y[m] - pred) ** 2).sum() / ss if ss > 0 else np.nan
        return p[0], p[1], r2
    except Exception:
        return np.nan, np.nan, np.nan


def bleach_rate_preact(df, args):
    """Per-reporter b (1/s) = -(mean per-cell slope of ln(norm)) over pre-activation frames."""
    out = {}
    sub = df[df.bright & df.Slice.between(args.bleach_first_frame, args.activation_frame)]
    for rep, g in sub.groupby("reporter"):
        slopes = g.groupby("cell_key").apply(
            lambda c: np.polyfit(c.time, np.log(c.norm_mean.clip(lower=1e-3)), 1)[0])
        b, se = -slopes.mean(), slopes.std() / np.sqrt(len(slopes))
        out[rep] = dict(b_raw=b, b=max(b, 0.0), ci_lo=b - 1.96 * se, ci_hi=b + 1.96 * se, n=len(slopes))
    return out


def bleach_rate_fit(df, args):
    """Per-reporter b shared by all bright cells, chosen by minimising the summed SSE of
    y = A exp(-b t) S(t) over a grid (A, y_inf, k per cell)."""
    t0, kmax = (args.activation_frame - 1) * args.interval, 1 / args.interval
    grid = np.linspace(0, 3e-4, 16)

    def model(b):
        def f(t, A, yinf, k):
            S = np.where(t < t0, 1.0, yinf + (1 - yinf) * np.exp(-k * np.clip(t - t0, 0, None)))
            return A * np.exp(-b * t) * S
        return f

    out = {}
    for rep, g in df[df.bright].groupby("reporter"):
        cells = [(c.time.values, c.norm_mean.values) for _, c in g.groupby("cell_key")]
        sse = []
        for b in grid:
            f, tot = model(b), 0.0
            for t, y in cells:
                try:
                    p, _ = curve_fit(f, t, y, p0=(1, np.mean(y[-3:]), 0.01),
                                     bounds=([0.5, 0, 1e-5], [2, 5, kmax]), maxfev=5000)
                    tot += ((y - f(t, *p)) ** 2).sum()
                except Exception:
                    pass
            sse.append(tot)
        sse = np.array(sse)
        i = int(np.argmin(sse))
        b = grid[i]
        if 0 < i < len(grid) - 1:
            a, bb, _ = np.polyfit(grid[i - 1:i + 2], sse[i - 1:i + 2], 2)
            if a > 0:
                b = -bb / (2 * a)
        out[rep] = dict(b_raw=b, b=max(b, 0.0), ci_lo=np.nan, ci_hi=np.nan, n=len(cells))
    return out


def per_cell(df, args, directions):
    """Per-cell fits and metrics; one row per cell, independent of treatment."""
    t0 = (args.activation_frame - 1) * args.interval
    kmax = 1 / args.interval
    t_slope_end = t0 + args.slope_frames * args.interval
    rows = []
    for key, c in df.groupby("cell_key"):
        t, y, yc = c.time.values, c.norm_mean.values, c.norm_corr.values
        rep = c.reporter.iloc[0]
        sign = 1 if directions[rep] == "import" else -1
        yinf, k, r2 = fit_relax(t, yc, t0, kmax)
        yinf0, k0, _ = fit_relax(t, y, t0, kmax)             # uncorrected, for comparison
        amp = yinf - 1
        w = (t >= t0) & (t <= t_slope_end)
        slope = np.polyfit(t[w], yc[w], 1)[0] * 60 * 100
        ln_slope = -np.polyfit(t[w], np.log(np.clip(yc[w], 1e-3, None)), 1)[0] * 60
        post = t >= t0
        gain_idx = np.argmin(np.abs(t - args.gain_time))
        rows.append(dict(
            cell_key=key, source=c.source.iloc[0], reporter=rep, direction=directions[rep],
            image=c.Image.iloc[0], image_base=c.image_base.iloc[0], position=c.position.iloc[0], roi=c.roi.iloc[0],
            area_px=c.Area.mean(), baseline_bgsub=c.baseline_bgsub.iloc[0], bright=bool(c.bright.iloc[0]),
            y_inf=yinf, amplitude=amp, k_per_min=k * 60, t_half_s=np.log(2) / k if k > 0 else np.nan, r2=r2,
            flux_pct_per_min=k * amp * 6000,
            flux_uncorrected_pct_per_min=k0 * (yinf0 - 1) * 6000,
            slope_pct_per_min=slope, ln_early_rate_per_min=ln_slope,
            gain_at_t=yc[gain_idx] - 1, auc_frac_min=_trapezoid(yc[post] - 1, t[post]) / 60,
            k_at_bound=bool(k <= 2e-4 or k >= 0.95 * kmax),
            responder=bool(sign * amp >= args.min_amp),
        ))
    r = pd.DataFrame(rows)
    r["good_fit"] = r.bright & r.responder & (r.r2 >= args.min_r2) & ~r.k_at_bound

    # plateau-subtracted log fit (well-fit cells): ln((y - y_inf)/(1 - y_inf)) = -k (t - t0)
    lnz = {}
    by_key = dict(tuple(df.groupby("cell_key")))
    for _, row in r[r.good_fit].iterrows():
        c = by_key[row.cell_key]
        t, yc = c.time.values, c.norm_corr.values
        z = (yc - row.y_inf) / (1 - row.y_inf)
        m = (t >= t0) & (t <= t0 + LNZ_WINDOW_S[row.direction]) & (z > 0.05)
        if m.sum() >= 3:
            lnz[row.cell_key] = -np.polyfit(t[m], np.log(z[m]), 1)[0] * 60
    r["k_lnz_per_min"] = r.cell_key.map(lnz)
    return r


TRACE_COLS = ["source", "row_id", "cell_key", "Image", "image_base", "roi", "Slice", "time", "reporter",
              "position", "Area", "Mean", "StdDev", "Min", "Max", "Median", "IntDen", "bg_mean", "Mean_SN",
              "Median_SN", "pre_mean", "pre_median", "baseline_bgsub", "bright", "norm_mean", "norm_median",
              "norm_corr"]


def fit_stage(raw, out, args):
    """Fit every cell once (independent of treatment) and cache the result in <out>/_cache/."""
    df, n_rois, n_kept = filter_rois(raw, args)
    df = normalise(df, args)
    print(f"ROIs: {n_rois} total, {n_kept} kept (full track, area >= {args.min_area:g} px), "
          f"{df[df.bright].cell_key.nunique()} bright")

    reps = sorted(df.reporter.unique())
    directions = dict(DIRECTIONS)
    t0 = (args.activation_frame - 1) * args.interval
    for rep in reps:
        if rep not in directions:
            late = df[(df.reporter == rep) & (df.time > t0)].norm_mean.median()
            directions[rep] = "import" if late > 1 else "export"
            print(f"Reporter {rep} not in DIRECTIONS; guessed '{directions[rep]}' from the data")

    # one bleaching rate per experiment (folder) x reporter
    bleach = {}
    for src, g in df.groupby("source"):
        if args.bleach == "preact":
            est = bleach_rate_preact(g, args)
        elif args.bleach == "fit":
            est = bleach_rate_fit(g, args)
        elif args.bleach == "none":
            est = {r: dict(b_raw=0.0, b=0.0, ci_lo=np.nan, ci_hi=np.nan, n=0) for r in g.reporter.unique()}
        else:
            raise ValueError("BLEACH_METHOD must be 'preact', 'fit' or 'none'")
        bleach.update({(src, r): v for r, v in est.items()})
    df["b_per_s"] = [bleach.get(k, {"b": 0.0})["b"] for k in zip(df.source, df.reporter)]
    t_ref = (args.activation_frame + 1) / 2 * args.interval - args.interval   # mid-baseline time
    df["norm_corr"] = df.norm_mean * np.exp(df.b_per_s * (df.time - t_ref))
    bleach_tab = pd.DataFrame([dict(source=src, reporter=r, direction=directions[r], method=args.bleach,
                                    b_pct_per_min_estimate=v["b_raw"] * 6000, b_pct_per_min_applied=v["b"] * 6000,
                                    ci_lo=v["ci_lo"] * 6000, ci_hi=v["ci_hi"] * 6000, n_cells=v["n"])
                               for (src, r), v in bleach.items()])

    fits = per_cell(df, args, directions)
    traces = df[[k for k in TRACE_COLS if k in df]]

    cache = os.path.join(out, "_cache")
    os.makedirs(cache, exist_ok=True)
    fits.to_csv(os.path.join(cache, "cell_fits.csv"), index=False)
    traces.to_csv(os.path.join(cache, "traces.csv.gz"), index=False)
    bleach_tab.to_csv(os.path.join(out, "bleach_rates.csv"), index=False)
    return fits, traces, bleach_tab


def load_cache(out):
    cache = os.path.join(out, "_cache")
    try:
        fits = pd.read_csv(os.path.join(cache, "cell_fits.csv"), dtype={"roi": str})
        traces = pd.read_csv(os.path.join(cache, "traces.csv.gz"), dtype={"roi": str})
        bleach_tab = pd.read_csv(os.path.join(out, "bleach_rates.csv"))
    except FileNotFoundError:
        raise FileNotFoundError(f"No saved fits in {cache}; run once with REASSIGN_ONLY = False first")
    return fits, traces, bleach_tab


# --------------------------------------------------------------------------- assign
def first_number(name):
    m = re.search(r"[\d.]+", str(name))
    try:
        return float(m.group()) if m else np.inf
    except ValueError:
        return np.inf


def order_treatments(names, explicit=None):
    names = list(dict.fromkeys(names))
    if explicit is not None:
        return [n for n in explicit if n in names]
    head = [n for n in TREATMENT_ORDER if n in names]
    rest = sorted([n for n in names if n not in head], key=lambda n: (first_number(n), n))
    return head + rest


def assign(fits, traces, meta):
    """Attach treatments: one row per cell x treatment (cells of an image with several treatments
    are used in each of them)."""
    m = meta[["source", "image_base", "reporter", "treatment", "replicate"]
             + [c for c in ["treat_time", "pyruvate"] if c in meta]]
    check = fits[["source", "image_base", "reporter"]].drop_duplicates().merge(
        m[["source", "image_base", "reporter"]].drop_duplicates(), on=["source", "image_base"],
        how="left", suffixes=("", "_meta"))
    changed = check[check.reporter_meta.notna() & (check.reporter != check.reporter_meta)]
    if len(changed):
        raise ValueError("The reporter of some images changed in the meta file (e.g. "
                         f"{changed.image_base.iloc[0]}); run with REASSIGN_ONLY = False")
    unassigned = check[check.reporter_meta.isna()].image_base.unique()
    if len(unassigned):
        print(f"  Note: {len(unassigned)} image(s) have no treatment in the meta file and were left out, "
              f"e.g. {unassigned[0]}")
    m = m.drop(columns="reporter")
    cells = fits.merge(m, on=["source", "image_base"], how="inner")
    cells = cells[~cells.reporter.isin(EXCLUDE_REPORTERS)]
    cells["cell_id"] = cells.cell_key + " | " + cells.treatment
    cells["unit"] = cells.source + "_" + cells.replicate
    cells["dose_mM"] = cells.treatment.str.extract(r"([\d.]+)\s*mM")[0].astype(float)
    tr = traces.merge(m[["source", "image_base", "treatment", "replicate"]], on=["source", "image_base"], how="inner")
    tr = tr[~tr.reporter.isin(EXCLUDE_REPORTERS)]
    tr["cell_id"] = tr.cell_key + " | " + tr.treatment
    n_multi = cells.groupby("cell_key").treatment.nunique().gt(1).sum()
    if n_multi:
        print(f"  {n_multi} cells belong to more than one treatment and are used in each")
    return cells, tr


def wrong_direction_qc(cells, tr, args, out):
    """Drop bright cells whose signal clearly moved the wrong way after activation (import falling,
    export rising): median of the bleach-corrected trace over WRONG_DIR_WINDOW_S after activation."""
    if WRONG_DIR_MAX is None:
        return cells, tr
    t0 = (args.activation_frame - 1) * args.interval
    w = tr[(tr.time > t0) & (tr.time <= t0 + WRONG_DIR_WINDOW_S)].drop_duplicates(["cell_key", "Slice"])
    early = w.groupby("cell_key").norm_corr.median() - 1
    sign = np.where(cells.direction == "import", 1, -1)
    cells = cells.assign(early_change=cells.cell_key.map(early).values)
    bad = cells.bright & (-sign * cells.early_change >= WRONG_DIR_MAX)
    flagged = cells[bad]
    flagged[["cell_id", "source", "replicate", "reporter", "direction", "treatment", "image", "roi",
             "early_change", "amplitude", "flux_pct_per_min", "r2"]].to_csv(
        os.path.join(out, "wrong_direction_cells.csv"), index=False)
    if len(flagged):
        print(f"  Wrong-direction QC: dropped {flagged.cell_key.nunique()} bright cells "
              f"(listed in wrong_direction_cells.csv)")
    return cells[~bad], tr[~tr.cell_id.isin(flagged.cell_id)]


def select(names, spec):
    """Treatments in a plot set, in plotting order."""
    if spec is None:
        return order_treatments(names)
    if isinstance(spec, str):
        return order_treatments([n for n in names if spec.lower() in n.lower()])
    return order_treatments(names, explicit=list(spec))


# --------------------------------------------------------------------------- stats
STAT_METRICS = ["flux_pct_per_min", "slope_pct_per_min", "t_half_s", "k_per_min",
                "gain_at_t", "auc_frac_min"]
GOODFIT_METRICS = {"t_half_s", "k_per_min"}    # tested on well-fit cells only


def holm(ps):
    ps = np.asarray(ps, float)
    order, m, run = np.argsort(ps), len(ps), 0.0
    adj = np.empty(m)
    for i, j in enumerate(order):
        run = max(run, (m - i) * ps[j])
        adj[j] = min(1.0, run)
    return adj


def stars(p):
    if p is None or np.isnan(p):
        return ""
    for lvl, s in STAR_LEVELS:
        if p < lvl:
            return s
    return "ns"


def replicate_tests(units, rep, col, cells_used, conds, control):
    """Tests on experiment+replicate medians. units: DataFrame with treatment, unit and col."""
    out = []
    groups = {c: units[units.treatment == c].set_index("unit")[col].dropna() for c in conds}
    groups = {c: g for c, g in groups.items() if len(g)}
    if control not in groups or len(groups) < 2 or any(len(g) < 2 for g in groups.values()):
        return out
    comp = [c for c in conds if c != control and c in groups]
    ps, names, ns = [], [], []
    for c in comp:
        shared = groups[control].index.intersection(groups[c].index)
        if len(shared) >= 2 and len(shared) == len(groups[c]):
            ps.append(stats.ttest_rel(groups[control][shared], groups[c][shared]).pvalue)
            names.append(f"Paired t-test vs {control}")
        else:
            ps.append(stats.ttest_ind(groups[control], groups[c], equal_var=False).pvalue)
            names.append(f"Welch t-test vs {control}")
    ps = [p if np.isfinite(p) else 1.0 for p in ps]
    for c, p, pa, name in zip(comp, ps, holm(ps), names):
        out.append(dict(reporter=rep, metric=col, cells=cells_used, level="replicates", test=name, treatment=c,
                        n_control=len(groups[control]), n_treatment=len(groups[c]), p=p, p_holm=pa, stars=stars(pa)))
    an = stats.f_oneway(*groups.values()).pvalue
    out.append(dict(reporter=rep, metric=col, cells=cells_used, level="replicates", test="ANOVA (replicate medians)",
                    treatment=np.nan, p=an, stars=stars(an)))
    return out


def stats_level(cells, conds):
    if STATS_UNIT in ("cells", "replicate"):
        return "replicates" if STATS_UNIT == "replicate" else "cells"
    n = cells[cells.bright].groupby(["reporter", "treatment"]).unit.nunique()
    n = n[n.index.get_level_values("treatment").isin(conds)]
    return "replicates" if len(n) and n.min() >= 3 else "cells"


def summarise(cells, conds, control):
    b = cells[cells.bright]
    gf = cells[cells.good_fit]
    others = [c for c in conds if c != control]
    rows, tests = [], []
    for rep, x in b.groupby("reporter"):
        ref_flux = x[x.treatment == control].flux_pct_per_min.median()
        for cond in conds:
            g = x[x.treatment == cond]
            if g.empty:
                continue
            f = gf[(gf.reporter == rep) & (gf.treatment == cond)]
            rows.append(dict(
                reporter=rep, direction=g.direction.iloc[0], treatment=cond, is_control=cond == control,
                n_bright=len(g), n_fields=g.image.nunique(), n_good_fit=len(f), pct_good_fit=100 * len(f) / len(g),
                pct_responders=100 * g.responder.mean(),
                flux_median=g.flux_pct_per_min.median(), flux_pct_of_control=100 * g.flux_pct_per_min.median() / ref_flux,
                flux_uncorrected_median=g.flux_uncorrected_pct_per_min.median(),
                slope_median=g.slope_pct_per_min.median(), gain_median=g.gain_at_t.median(),
                auc_median=g.auc_frac_min.median(),
                k_per_min_goodfit=f.k_per_min.median(), t_half_s_goodfit=f.t_half_s.median(),
                amplitude_goodfit=f.amplitude.median(), k_lnz_per_min_goodfit=f.k_lnz_per_min.median()))
        for col in STAT_METRICS:
            xx = gf[gf.reporter == rep] if col in GOODFIT_METRICS else x
            cells_used = "well-fit" if col in GOODFIT_METRICS else "bright"
            present = [c for c in conds if (xx.treatment == c).any()]
            groups = {c: xx[xx.treatment == c][col].dropna() for c in present}
            if control not in groups or len(present) < 2 or any(len(g) < 2 for g in groups.values()):
                continue
            kw = stats.kruskal(*groups.values()).pvalue
            fields = xx.groupby(["treatment", "image"])[col].median().reset_index()
            fgroups = [g[col].dropna() for _, g in fields.groupby("treatment")]
            kwf = stats.kruskal(*fgroups).pvalue if all(len(g) >= 2 for g in fgroups) else np.nan
            ucol = xx.groupby(["treatment", "unit"])[col].median().reset_index()
            tests += replicate_tests(ucol, rep, col, cells_used, conds, control)
            comp = [c for c in others if c in groups]
            ps = [stats.mannwhitneyu(groups[control], groups[c]).pvalue for c in comp]
            for c, p, pa in zip(comp, ps, holm(ps)):
                tests.append(dict(reporter=rep, metric=col, cells=cells_used, level="cells", test=f"Mann-Whitney vs {control}",
                                  treatment=c, n_control=len(groups[control]), n_treatment=len(groups[c]),
                                  p=p, p_holm=pa, stars=stars(pa)))
            tests.append(dict(reporter=rep, metric=col, cells=cells_used, level="cells", test="Kruskal-Wallis (cells)",
                              treatment=np.nan, p=kw, stars=stars(kw)))
            tests.append(dict(reporter=rep, metric=col, cells=cells_used, level="cells", test="Kruskal-Wallis (field medians)",
                              treatment=np.nan, p=kwf, stars=stars(kwf)))
        # fraction of bright cells that are well-fit: Fisher's exact test vs control, chi-square across treatments
        present = [c for c in conds if (x.treatment == c).any()]
        if control in present and len(present) > 1:
            table = {c: [int(x[x.treatment == c].good_fit.sum()), int((~x[x.treatment == c].good_fit).sum())]
                     for c in present}
            comp = [c for c in others if c in table]
            ps = [stats.fisher_exact([table[control], table[c]])[1] for c in comp]
            for c, p, pa in zip(comp, ps, holm(ps)):
                tests.append(dict(reporter=rep, metric="pct_good_fit", cells="bright", level="cells",
                                  test=f"Fisher exact vs {control}", treatment=c, n_control=sum(table[control]),
                                  n_treatment=sum(table[c]), p=p, p_holm=pa, stars=stars(pa)))
            pu = x.groupby(["treatment", "unit"]).good_fit.mean().mul(100).rename("pct_good_fit").reset_index()
            tests += replicate_tests(pu, rep, "pct_good_fit", "bright", conds, control)
            tab = list(table.values())
            chi = stats.chi2_contingency(tab)[1] if min(map(sum, zip(*tab))) > 0 else np.nan
            tests.append(dict(reporter=rep, metric="pct_good_fit", cells="bright", level="cells", test="Chi-square (treatments)",
                              treatment=np.nan, p=chi, stars=stars(chi)))
    tests = pd.DataFrame(tests)
    if tests.empty:
        tests = pd.DataFrame(columns=["reporter", "metric", "cells", "level", "test", "treatment", "n_control",
                                      "n_treatment", "p", "p_holm", "stars"])
    return pd.DataFrame(rows), tests


def timecourse_tables(df, conds):
    """Mean / SD / SEM / n per reporter x treatment x time, bright cells, bleach-corrected and uncorrected."""
    b = df[df.bright]
    long = (b.groupby(["reporter", "treatment", "time"])
             .agg(n=("norm_corr", "size"), mean_corr=("norm_corr", "mean"), sd_corr=("norm_corr", "std"),
                  mean_uncorr=("norm_mean", "mean"), sd_uncorr=("norm_mean", "std"))
             .reset_index())
    long["sem_corr"] = long.sd_corr / np.sqrt(long.n)
    long["sem_uncorr"] = long.sd_uncorr / np.sqrt(long.n)
    long["order"] = long.treatment.map({c: i for i, c in enumerate(conds)})
    long = long.sort_values(["reporter", "order", "time"])[
        ["reporter", "treatment", "time", "n", "mean_corr", "sd_corr", "sem_corr",
         "mean_uncorr", "sd_uncorr", "sem_uncorr"]]
    wide = {}
    for rep, g in long.groupby("reporter"):
        w = g.pivot(index="time", columns="treatment", values=["mean_corr", "sd_corr", "sem_corr", "n"])
        w.columns = [f"{c}_{v.replace('_corr', '')}" for v, c in w.columns]
        order = [f"{c}_{v}" for c in conds if (g.treatment == c).any() for v in ["mean", "sd", "sem", "n"]]
        wide[rep] = w[order].reset_index()
    return long, wide


# --------------------------------------------------------------------------- figures
def style(ax):
    for sp in ["top", "right"]:
        ax.spines[sp].set_visible(False)
    ax.grid(axis="y", color="#eee")
    ax.set_axisbelow(True)


def cond_colors(conds):
    return dict(zip(conds, ramp(DOSE_COLORS, len(conds))))


def set_cond_ticks(ax, labels):
    long = max((len(s.split("\n")[0]) for s in labels), default=0) > 9
    ax.set_xticks(range(len(labels)), labels, fontsize=8, rotation=30 if long else 0,
                  ha="right" if long else "center")
    ax.set_xlabel("treatment")


def save(fig, path):
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)


def fig_timecourse(df, conds, rep, direction, args, path, err="sem"):
    t0 = (args.activation_frame - 1) * args.interval
    fig, ax = plt.subplots(figsize=(6.4, 4.3))
    x = df[(df.reporter == rep) & df.bright]
    cols = cond_colors(conds)
    for cond in conds:
        g = x[x.treatment == cond]
        if g.empty:
            continue
        m = g.groupby("time").norm_corr.agg(["mean", {"sd": "std"}.get(err, err)])
        e = m.iloc[:, 1]
        ax.fill_between(m.index, m["mean"] - e, m["mean"] + e, color=cols[cond], alpha=.18, lw=0)
        ax.plot(m.index, m["mean"], color=cols[cond], lw=2, label=f"{cond} (n={g.cell_id.nunique()})")
    ax.axhline(1, color="#999", lw=.8, ls="--")
    ax.axvspan(0, t0, color="#eee", zorder=0)
    ax.set_title(f"Reporter {rep} ({direction}): mean ± {err.upper()}, bright cells\n(grey = pre-activation)",
                 loc="left", fontsize=10)
    ax.set_xlabel("time (s)")
    ax.set_ylabel("nuclear signal / baseline")
    ax.legend(frameon=False, fontsize=8, loc="upper left", bbox_to_anchor=(1.01, 1))
    style(ax)
    save(fig, path)


UNIT_MARKERS = ["o", "s", "^", "D", "v", "P", "X", "<", ">"]


def strip(ax, x, col, conds, colors, fmt="{:.1f}", clip=None):
    """Superplot: every cell as a small grey dot, an unfilled box per treatment, and the median of
    each experiment + replicate as a large dot (marker shape = experiment + replicate)."""
    rng = np.random.default_rng(1)
    units = sorted(x.unit.unique())
    markers = {u: UNIT_MARKERS[i % len(UNIT_MARKERS)] for i, u in enumerate(units)}
    for i, cond in enumerate(conds):
        color = colors[i]
        g = x[(x.treatment == cond) & x[col].notna()]
        if len(g) == 0:
            continue
        v = g[col].values
        vv = np.clip(v, *clip) if clip else v
        ax.scatter(i + rng.uniform(-.3, .3, len(v)), vv, s=5, color=CELL_DOT_COLOR, alpha=.5, lw=0, zorder=1)
        ax.boxplot(vv, positions=[i], widths=.55, showfliers=False, patch_artist=False, zorder=2,
                   boxprops=dict(color=color, lw=1.5), whiskerprops=dict(color=color, lw=1.5),
                   capprops=dict(color=color, lw=1.5), medianprops=dict(color=color, lw=1.5))
        meds = g.groupby("unit")[col].median()
        off = np.linspace(-.12, .12, len(meds)) if len(meds) > 1 else [0]
        for o, (u, m) in zip(off, meds.items()):
            m = np.clip(m, *clip) if clip else m
            ax.scatter(i + o, m, s=55, marker=markers[u], color=color, edgecolor="white", lw=.8,
                       alpha=.9, zorder=3)
    if len(units) > 1:
        from matplotlib.lines import Line2D
        ax.legend([Line2D([], [], ls="", marker=markers[u], color="#555", ms=6) for u in units], units,
                  frameon=False, fontsize=7, title="experiment_replicate median", title_fontsize=7,
                  loc="upper left", bbox_to_anchor=(1.01, 1))
    ax.set_xlim(-.6, len(conds) - .4)
    set_cond_ticks(ax, [f"{c}\n(n={x[x.treatment == c][col].notna().sum()})" for c in conds])
    style(ax)


STARS_LEVEL = "cells"        # set per plot set by report() from STATS_UNIT


def add_stars(ax, tests, rep, metric, conds, control):
    if not SHOW_STARS or tests.empty:
        return
    t = tests[(tests.reporter == rep) & (tests.metric == metric) & (tests.level == STARS_LEVEL)]
    if t.empty:
        return
    lo, hi = ax.get_ylim()
    y = hi * (lo / hi) ** 0.03 if ax.get_yscale() == "log" else hi - 0.03 * (hi - lo)
    for i, cond in enumerate(conds):
        s = t[t.treatment == cond]
        if len(s):
            ax.text(i, y, s.stars.iloc[0], ha="center", va="top", fontsize=10, color="#222")
    pair = t[t.p_holm.notna()].test.iloc[0].split(" vs ")[0] if t.p_holm.notna().any() else "test"
    kw = t[t.test.isin(["Kruskal-Wallis (cells)", "Chi-square (treatments)", "ANOVA (replicate medians)"])]
    unit = "experiment+replicate medians" if STARS_LEVEL == "replicates" else "cells"
    note = f"Stars: {pair} vs {control} on {unit}, Holm-adjusted"
    if len(kw):
        note += f"   {kw.test.iloc[0].split(' (')[0]} p = {kw.p.iloc[0]:.2g}"
    ax.annotate(note, xy=(0, 0), xycoords="figure fraction", xytext=(6, 4), textcoords="offset points",
                fontsize=7, color="#666")


RATE_COLS = {"flux_pct_per_min", "slope_pct_per_min"}


def rate_scale(direction):
    """Factor that converts a stored rate (signed % of baseline per min) to the plotted value."""
    sign = -1 if (EXPORT_POSITIVE and direction == "export") else 1
    return sign * (0.01 if RATE_UNIT == "per_min" else 1)


def rate_label(direction):
    unit = "min⁻¹" if RATE_UNIT == "per_min" else "% of baseline per min"
    if EXPORT_POSITIVE and direction == "export":
        return f"export rate, sign reversed ({unit})"
    return f"{direction} rate ({unit})"


def fig_metric(cells, tests, conds, control, rep, direction, metric, title, ylabel, path, goodfit=False,
               fmt="{:.1f}", clip=None, zero_line=True):
    x = cells[cells.good_fit if goodfit else cells.bright]
    x = x[x.reporter == rep]
    log = False
    if metric in RATE_COLS:
        x = x.assign(**{metric: x[metric] * rate_scale(direction)})
        ylabel = rate_label(direction)
        if RATE_YSCALE == "log":
            log, n_drop = True, int((x[metric] <= 0).sum())
            x = x[x[metric] > 0]
            zero_line = False
            title += f"\n(log scale; {n_drop} cells with rate ≤ 0 not shown)"
    fig, ax = plt.subplots(figsize=(max(6, 1 + 1.0 * len(conds)), 5))
    if log:
        ax.set_yscale("log")
    strip(ax, x, metric, conds, shades(direction, len(conds)), fmt=fmt, clip=clip)
    lo, hi = ax.get_ylim()
    ax.set_ylim(lo, hi * (hi / lo) ** 0.12 if log else hi + 0.12 * (hi - lo))   # headroom for stars
    if zero_line:
        ax.axhline(0, color="#999", lw=.8, ls="--")
    ax.set_title(title, loc="left", fontsize=10)
    ax.set_ylabel(ylabel)
    fig.tight_layout(rect=(0, 0.03, 1, 1))
    add_stars(ax, tests, rep, metric, conds, control)
    fig.savefig(path, dpi=150)
    plt.close(fig)


def fig_pct_goodfit(cells, tests, conds, control, rep, direction, args, path):
    """% of bright cells that pass the well-fit criteria, per treatment."""
    x = cells[cells.bright & (cells.reporter == rep)]
    k = np.array([x[x.treatment == c].good_fit.sum() for c in conds])
    n = np.array([(x.treatment == c).sum() for c in conds])
    pct = np.where(n > 0, 100 * k / np.maximum(n, 1), np.nan)
    fig, ax = plt.subplots(figsize=(max(6, 1 + 1.0 * len(conds)), 5))
    ax.bar(range(len(conds)), pct, 0.6, color=REPORTER_COLORS[direction])
    for i in range(len(conds)):
        if n[i]:
            ax.text(i, pct[i] + 1.5, f"{pct[i]:.0f}%\n({k[i]}/{n[i]})", ha="center", va="bottom", fontsize=8)
    set_cond_ticks(ax, conds)
    ax.set_ylabel("well-fit cells (% of bright cells)")
    ax.set_ylim(0, 125)
    ax.set_yticks(range(0, 101, 20))
    ax.set_title(f"{rep} ({direction}): well-fit cells\n(R² ≥ {args.min_r2:g}, |amplitude| ≥ {args.min_amp:g} in the "
                 f"expected direction, k not at a bound)", loc="left", fontsize=10)
    style(ax)
    fig.tight_layout(rect=(0, 0.03, 1, 1))
    add_stars(ax, tests, rep, "pct_good_fit", conds, control)
    fig.savefig(path, dpi=150)
    plt.close(fig)


def fig_logfits(df, cells, conds, rep, direction, args, path):
    t0 = (args.activation_frame - 1) * args.interval
    gf = cells[cells.good_fit].set_index("cell_id")
    fig, axs = plt.subplots(2, 1, figsize=(6.4, 8.4))
    x = df[(df.reporter == rep) & df.bright & (df.time >= t0)].copy()
    cols = cond_colors(conds)
    x["ln"] = np.log(x.norm_corr.clip(lower=1e-3))
    ax = axs[0]
    for cond in conds:
        g = x[x.treatment == cond]
        if g.empty:
            continue
        m = g.groupby("time").ln.agg(["mean", "sem"])
        ax.errorbar(m.index, m["mean"], yerr=m["sem"], fmt="o", ms=3, color=cols[cond], lw=1, label=cond)
        e = m[m.index <= t0 + args.slope_frames * args.interval]
        p = np.polyfit(e.index, e["mean"], 1)
        tt = np.linspace(t0, m.index.max(), 30)
        ax.plot(tt, np.polyval(p, tt), color=cols[cond], lw=1, ls="--")
    ax.axhline(0, color="#999", lw=.8)
    ax.set_title(f"{rep} ({direction}): ln(norm); dashed = early-window fit", loc="left", fontsize=10)
    ax.legend(frameon=False, fontsize=8, loc="upper left", bbox_to_anchor=(1.01, 1))
    style(ax)
    ax = axs[1]
    y = x[x.cell_id.isin(gf.index) & (x.time <= t0 + LNZ_WINDOW_S[direction])].copy()
    y["yinf"] = y.cell_id.map(gf.y_inf)
    y["z"] = (y.norm_corr - y.yinf) / (1 - y.yinf)
    y = y[y.z > 0.05]
    y["lnz"] = np.log(y.z)
    for cond in conds:
        g = y[y.treatment == cond]
        m = g.groupby("time").lnz.agg(["mean", "sem", "size"])
        m = m[m["size"] >= 5]
        if len(m) < 2:
            continue
        ax.errorbar(m.index, m["mean"], yerr=m["sem"], fmt="o", ms=3, color=cols[cond], lw=1,
                    label=f"{cond} (n={g.cell_id.nunique()})")
        p = np.polyfit(m.index, m["mean"], 1)
        ax.plot(m.index, np.polyval(p, m.index), color=cols[cond], lw=1.2)
    ax.axhline(0, color="#999", lw=.8)
    ax.set_title(f"{rep}: ln((norm − plateau)/(1 − plateau)), well-fit cells; slope = −k", loc="left", fontsize=10)
    ax.set_xlabel("time (s)")
    ax.legend(frameon=False, fontsize=8, loc="upper left", bbox_to_anchor=(1.01, 1))
    style(ax)
    save(fig, path)


def fig_bleach(bleach_tab, cells, conds, control, reps, directions, args, path):
    fig, axs = plt.subplots(1, 1 + len(reps), figsize=(5 + max(5.5, 1 + .9 * len(conds)) * len(reps), 4.6),
                            squeeze=False)
    ax = axs[0][0]
    bt = bleach_tab.sort_values(["reporter", "source"]).reset_index(drop=True)
    multi = bt.source.nunique() > 1
    for i, r in bt.iterrows():
        c = REPORTER_COLORS[directions[r.reporter]]
        b = r.b_pct_per_min_estimate
        yerr = None if np.isnan(r.ci_lo) else [[b - r.ci_lo], [r.ci_hi - b]]
        ax.errorbar(i, b, yerr=yerr, fmt="o", color=c, ms=8, capsize=4)
        ax.text(i + .08, b, f"{b:.2f}", fontsize=8, va="center")
    ax.axhline(0, color="#999", lw=.8)
    ax.set_xticks(range(len(bt)), [f"{r.reporter}\n{r.source}" if multi else f"{r.reporter} ({directions[r.reporter]})"
                                   for r in bt.itertuples()], fontsize=8 if multi else 10)
    ax.set_xlim(-.5, len(bt) - .5)
    ax.set_ylabel("bleaching rate b (% per min)")
    ax.set_title(f"Bleaching rate per {'experiment × ' if multi else ''}reporter ({args.bleach})", loc="left", fontsize=10)
    style(ax)
    b = cells[cells.bright]
    for ax, rep in zip(axs[0][1:], reps):
        x = b[b.reporter == rep]
        w = .38
        for k, (col, lab, alpha) in enumerate([("flux_uncorrected_pct_per_min", "uncorrected", .35),
                                               ("flux_pct_per_min", "bleach-corrected", .9)]):
            med = x.groupby("treatment")[col].median().reindex(conds)
            pct = med / med.get(control, np.nan) * 100
            ax.bar(np.arange(len(conds)) + (k - .5) * w, pct, w * .95, color=REPORTER_COLORS[directions[rep]],
                   alpha=alpha, label=lab)
            for i, v in enumerate(pct):
                if np.isfinite(v):
                    ax.text(i + (k - .5) * w, max(v, 0) + 1.5, f"{v:.0f}", ha="center", fontsize=7.5)
        set_cond_ticks(ax, conds)
        ax.set_ylim(min(0, ax.get_ylim()[0]), max(125, ax.get_ylim()[1]))
        ax.set_title(f"{rep}: median initial flux, % of {control}", loc="left", fontsize=10)
        ax.legend(frameon=False, fontsize=8)
        style(ax)
    save(fig, path)


RATE_NAMES = {"k_per_min": ("k", "min⁻¹", "k"),
              "flux_pct_per_min": ("initial flux", "% of baseline per min", "flux"),
              "slope_pct_per_min": ("model-free slope", "% of baseline per min", "slope")}


def fig_rate_scatter(cells, conds, directions, metric, path):
    """Export rate (x) vs import rate (y) per treatment. Faint dots: median of each experiment + replicate
    (import and export paired by experiment, replicate and treatment). Large dot: mean of those medians
    (or the cell median when error = cells) with ± SD error bars."""
    imp = [r for r, d in directions.items() if d == "import" and r in set(cells.reporter)]
    exp = [r for r, d in directions.items() if d == "export" and r in set(cells.reporter)]
    if not imp or not exp:
        return False
    goodfit = metric in GOODFIT_METRICS
    c = cells[cells.good_fit if goodfit else cells.bright].copy()
    # flip export's sign so both axes read as positive rates (taking abs() would turn a blocked
    # import cell's slightly negative value into a positive rate and inflate the median)
    if metric in RATE_COLS:
        c["v"] = c[metric] * np.where(c.direction == "export", -1, 1) * (0.01 if RATE_UNIT == "per_min" else 1)
    else:
        c["v"] = c[metric]
    xi, xe = c[c.reporter == imp[0]], c[c.reporter == exp[0]]
    ui = xi.groupby(["treatment", "unit"]).v.median()
    ue = xe.groupby(["treatment", "unit"]).v.median()
    units = pd.concat([ue.rename("x"), ui.rename("y")], axis=1).dropna().reset_index()
    both = [cd for cd in conds if (xi.treatment == cd).any() and (xe.treatment == cd).any()]
    if not both:
        return False
    err = RATE_SCATTER_ERROR
    if err == "auto":
        err = "units" if all((units.treatment == cd).sum() >= 2 for cd in both) else "cells"
    cols = cond_colors(conds)
    name, unit_label, _ = RATE_NAMES[metric]
    if metric in RATE_COLS and RATE_UNIT == "per_min":
        unit_label = "min⁻¹"
    fig, ax = plt.subplots(figsize=(6.8, 4.6))
    for cd in both:
        u = units[units.treatment == cd]
        ax.scatter(u.x, u.y, s=40, color=cols[cd], alpha=.35, lw=0, zorder=2)
        if err == "units":
            mx, my, sx, sy = u.x.mean(), u.y.mean(), u.x.std(), u.y.std()
        else:
            ve, vi = xe[xe.treatment == cd].v, xi[xi.treatment == cd].v
            mx, my, sx, sy = ve.median(), vi.median(), ve.std(), vi.std()
        ax.errorbar(mx, my, xerr=sx, yerr=sy, fmt="o", ms=8, color=cols[cd], ecolor=cols[cd], elinewidth=1.3,
                    capsize=4, zorder=3, label=cd)
    ax.set_xlim(left=min(0, ax.get_xlim()[0]))      # include 0 without hiding error bars below it
    ax.set_ylim(bottom=min(0, ax.get_ylim()[0]))
    ax.set_xlabel(f"export rate, {name}{', sign reversed' if metric in RATE_COLS else ''} ({unit_label})")
    ax.set_ylabel(f"import rate, {name} ({unit_label})")
    centre = "mean of experiment+replicate medians" if err == "units" else "median of cells"
    ax.set_title(f"Export vs import rate: {centre} ± SD\n({imp[0]} import, {exp[0]} export; "
                 f"{'well-fit' if goodfit else 'bright'} cells; faint dots = experiment+replicate medians)",
                 loc="left", fontsize=9)
    ax.legend(frameon=False, fontsize=8, loc="upper left", bbox_to_anchor=(1.01, 1))
    style(ax)
    ax.grid(axis="x", color="#eee")
    save(fig, path)
    return True


def cell_plots(traces, cells, args, outdir):
    """One PNG per cell: normalised and bleach-corrected trace with the fitted exponential.
    Saved under <reporter>_<treatment>/; a cell in several treatments is saved in each folder."""
    t0 = (args.activation_frame - 1) * args.interval
    info = cells.drop_duplicates("cell_key").set_index("cell_key")
    folders = cells.groupby("cell_key").apply(
        lambda g: [f"{r}_{safe_name(t)}" for r, t in zip(g.reporter, g.treatment)])
    if CELL_PLOTS_BRIGHT_ONLY:
        traces = traces[traces.bright]
    n = 0
    for key, c in traces.drop_duplicates(["cell_key", "Slice"]).groupby("cell_key"):
        if key not in info.index:
            continue
        r = info.loc[key]
        fig, ax = plt.subplots(figsize=(5, 3.6))
        ax.axvspan(0, t0, color="#eee", zorder=0)
        ax.plot(c.time, c.norm_mean, "o", ms=3, color="#aaa", label="normalised")
        ax.plot(c.time, c.norm_corr, "o", ms=3.5, color=REPORTER_COLORS[r.direction], label="bleach-corrected")
        if np.isfinite(r.k_per_min):
            tt = np.linspace(t0, c.time.max(), 200)
            ax.plot(tt, relax(tt, r.y_inf, r.k_per_min / 60, t0), color="#222", lw=1.5, label="exponential fit")
        ax.axhline(1, color="#999", lw=.8, ls="--")
        verdict = ("well-fit" if r.good_fit else "dim" if not r.bright else
                   "non-responder" if not r.responder else "k at bound" if r.k_at_bound else "low R²")
        treat = ", ".join(cells[cells.cell_key == key].treatment)
        ax.set_title(f"{os.path.basename(r.image)} ROI {r.roi}  ({r.reporter}; {treat})\n" +
                     (f"flux {r.flux_pct_per_min / 100:.3f} min⁻¹" if RATE_UNIT == "per_min" else f"flux {r.flux_pct_per_min:.1f} %/min") + f" · amp {r.amplitude:.2f} · "
                     f"t½ {r.t_half_s:.0f} s · R² {r.r2:.2f} · {verdict}", loc="left", fontsize=8)
        ax.set_xlabel("time (s)")
        ax.set_ylabel("nuclear signal / baseline")
        ax.legend(frameon=False, fontsize=7)
        style(ax)
        fig.tight_layout()
        for sub in folders[key]:
            os.makedirs(os.path.join(outdir, sub), exist_ok=True)
            fig.savefig(os.path.join(outdir, sub, f"{r.image.replace('/', '_')}_ROI{r.roi}.png"), dpi=CELL_PLOT_DPI)
        plt.close(fig)
        n += 1
    return n


# --------------------------------------------------------------------------- report
def report(cells, tr, bleach_tab, conds, control, directions, args, out):
    """Summaries, statistics and plots for one plot set."""
    os.makedirs(out, exist_ok=True)
    cells = cells[cells.treatment.isin(conds)]
    tr = tr[tr.treatment.isin(conds)]
    global STARS_LEVEL
    reps = sorted(cells.reporter.unique())
    summary, tests = summarise(cells, conds, control)
    STARS_LEVEL = stats_level(cells, conds)
    print(f"  Statistics on plots: {'experiment+replicate medians' if STARS_LEVEL == 'replicates' else 'cells'}")
    tc_long, tc_wide = timecourse_tables(tr, conds)

    keep = ["source", "row_id", "cell_id", "Image", "roi", "Slice", "time", "reporter", "treatment", "replicate",
            "position", "Area", "Mean", "StdDev", "Min", "Max", "Median", "IntDen", "bg_mean", "Mean_SN",
            "Median_SN", "pre_mean", "pre_median", "baseline_bgsub", "bright", "norm_mean", "norm_median",
            "norm_corr"]
    tr[[k for k in keep if k in tr]].to_csv(os.path.join(out, "timecourses_long.csv"), index=False)
    tc_long.to_csv(os.path.join(out, "timecourse_mean_by_treatment.csv"), index=False)
    first = ["cell_id", "source", "replicate", "unit", "reporter", "direction", "treatment", "dose_mM", "image",
             "position", "roi"]
    cells[first + [c for c in cells if c not in first]].to_csv(os.path.join(out, "rates_per_cell.csv"), index=False)
    rep_cols = ["flux_pct_per_min", "slope_pct_per_min", "gain_at_t", "auc_frac_min"]
    fit_cols = ["k_per_min", "t_half_s", "amplitude"]
    b_ = cells[cells.bright].groupby(["reporter", "treatment", "source", "replicate"])
    f_ = cells[cells.good_fit].groupby(["reporter", "treatment", "source", "replicate"])
    pd.concat([b_.size().rename("n_bright"), b_[rep_cols].median(),
               f_.size().rename("n_good_fit"), f_[fit_cols].median().add_suffix("_goodfit")], axis=1
              ).reset_index().to_csv(os.path.join(out, "replicate_medians.csv"), index=False)
    summary.to_csv(os.path.join(out, "summary_by_treatment.csv"), index=False)
    tests.to_csv(os.path.join(out, "statistics.csv"), index=False)

    for rep in reps:
        d = directions[rep]
        tc_wide[rep].to_csv(os.path.join(out, f"timecourse_{rep}_wide.csv"), index=False)
        for err in TIMECOURSE_ERRORS:
            suffix = "" if err == "sem" else f"_{err}"
            fig_timecourse(tr, conds, rep, d, args, os.path.join(out, f"timecourse_{rep}{suffix}.png"), err=err)
        ref = cells[cells.bright & (cells.reporter == rep) & (cells.treatment == control)].flux_pct_per_min.median()
        fig_metric(cells, tests, conds, control, rep, d, "flux_pct_per_min",
                   f"{rep} ({d}): initial flux k × amplitude, all bright cells\n({control} median "
                   f"{ref * rate_scale(d):.3g})",
                   "% of baseline per min", os.path.join(out, f"flux_{rep}.png"))
        fig_metric(cells, tests, conds, control, rep, d, "slope_pct_per_min",
                   f"{rep} ({d}): model-free initial slope, all bright cells\n(linear fit, frames "
                   f"{ACTIVATION_FRAME}–{ACTIVATION_FRAME + SLOPE_FRAMES})", "% of baseline per min",
                   os.path.join(out, f"slope_{rep}.png"))
        fig_metric(cells, tests, conds, control, rep, d, "t_half_s",
                   f"{rep} ({d}): half-time of fitted exponential, well-fit cells", "t½ (s)",
                   os.path.join(out, f"half_time_{rep}.png"), goodfit=True, fmt="{:.0f}", clip=(0, 600),
                   zero_line=False)
        fig_metric(cells, tests, conds, control, rep, d, "k_per_min",
                   f"{rep} ({d}): rate constant k of fitted exponential, well-fit cells", "k (min⁻¹)",
                   os.path.join(out, f"k_{rep}.png"), goodfit=True, fmt="{:.2f}", zero_line=False)
        fig_pct_goodfit(cells, tests, conds, control, rep, d, args, os.path.join(out, f"pct_well_fit_{rep}.png"))
        fig_logfits(tr, cells, conds, rep, d, args, os.path.join(out, f"log_fits_{rep}.png"))
    fig_bleach(bleach_tab, cells, conds, control, reps, directions, args, os.path.join(out, "bleach_correction.png"))
    for m in RATE_SCATTER_METRICS:
        fig_rate_scatter(cells, conds, directions, m, os.path.join(out, f"export_vs_import_{RATE_NAMES[m][2]}.png"))

    pd.set_option("display.width", 200)
    print(summary[["reporter", "direction", "treatment", "n_bright", "pct_responders", "flux_median",
                   "flux_pct_of_control", "slope_median", "t_half_s_goodfit"]].round(2).to_string(index=False))


def analyse(raw, meta, out, args):
    """raw: measurements (None when REASSIGN_ONLY); meta: treatment table(s) with a source column."""
    os.makedirs(out, exist_ok=True)
    if raw is None:
        fits, traces, bleach_tab = load_cache(out)
        print(f"Reusing saved fits ({fits.cell_key.nunique()} cells); re-reading treatments only")
    else:
        fits, traces, bleach_tab = fit_stage(raw, out, args)
    directions = dict(zip(fits.reporter, fits.direction))
    print("Bleaching (%/min):", {f"{r.source} {r.reporter}": round(r.b_pct_per_min_applied, 3)
                                 for r in bleach_tab.itertuples()})

    cells, tr = assign(fits, traces, meta)
    cells, tr = wrong_direction_qc(cells, tr, args, out)
    names = list(cells.treatment.unique())
    for set_name, spec in PLOT_SETS.items():
        conds = select(names, spec)
        if not conds:
            print(f"Plot set '{set_name}': no matching treatments, skipped")
            continue
        control = CONTROLS.get(set_name, conds[0])
        if control not in conds:
            print(f"Plot set '{set_name}': control '{control}' not found; using '{conds[0]}'")
            control = conds[0]
        set_out = out if len(PLOT_SETS) == 1 else os.path.join(out, safe_name(set_name))
        print(f"--- plot set '{set_name}': {', '.join(conds)} (control: {control})")
        report(cells, tr, bleach_tab, conds, control, directions, args, set_out)

    if SAVE_CELL_PLOTS:
        n = cell_plots(tr, cells, args, os.path.join(out, "cell_plots"))
        print(f"Saved {n} per-cell plots")
    print(f"Outputs written to {os.path.abspath(out)}\n")


def main():
    args = settings()
    if not BATCH_FOLDERS:
        jobs = [(DATA_FOLDER, DATA_FOLDER)]
    else:
        jobs = [(f, os.path.join(f, BATCH_DATA_SUBFOLDER)) for f in BATCH_FOLDERS]
    pooled_only = POOLED_ONLY and bool(BATCH_FOLDERS)
    if pooled_only and not COMBINE_BATCH:
        print("POOLED_ONLY is set but COMBINE_BATCH is False; running each folder only")
    loaded = []
    for folder, data_folder in jobs:
        print(f"=== {folder}")
        out = os.path.join(data_folder, OUTPUT_SUBFOLDER)
        try:
            if REASSIGN_ONLY:
                source, _ = experiment(folder, data_folder)
                raw, meta = None, load_meta(find_file(data_folder, META_FILE), source)
            else:
                raw, meta = load_folder(folder, data_folder)
            if not (pooled_only and COMBINE_BATCH):
                analyse(raw, meta, out, args)
            loaded.append((raw, meta))
        except Exception as e:
            if not BATCH_FOLDERS:
                raise
            print(f"Skipped {folder}: {e}\n")
    if BATCH_FOLDERS and COMBINE_BATCH and len(loaded) > 1:
        print("=== pooled")
        out = POOLED_OUTPUT if os.path.isabs(POOLED_OUTPUT) else \
            os.path.join(os.path.dirname(os.path.abspath(__file__)), POOLED_OUTPUT)
        meta = pd.concat([m for _, m in loaded], ignore_index=True)
        raw = None
        if not REASSIGN_ONLY:
            pooled = []
            for r, _ in loaded:      # keep image names unique across folders
                r = r.copy()
                r["Image"] = r["source"] + "/" + r["Image"]
                pooled.append(r)
            raw = pd.concat(pooled, ignore_index=True)
        analyse(raw, meta, out, args)


if __name__ == "__main__":
    main()
