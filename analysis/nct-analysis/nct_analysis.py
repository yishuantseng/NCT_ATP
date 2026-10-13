#!/usr/bin/env python3
"""
Nuclear transport (NCT) analysis of photoactivatable reporter time-lapses.

HOW TO USE
  1. Edit the PARAMETERS block below (file names, folders, settings).
  2. Run the file:  python nct_analysis.py   (or press Run in your editor).
  Requires: numpy, pandas, scipy, matplotlib.

INPUT FILES (CSV) in each data folder
  RESULTS_FILE  Fiji/ImageJ "Measure" results for every nucleus ROI in every frame
                (Label like "<image>_561.tif:<slice> - ID <roi>:<slice>")
  BG_FILE       Background intensity per image and frame (Label like "<image>_561.tif", Slice 0-based)
  META_FILE     Image -> condition table (image, reporter, treatment, time, pyruvate, replicate;
                image = name without the _XY<n> field suffix)

PIPELINE
  1. Parse labels, attach conditions and per-frame background.
  2. Keep ROIs tracked through every frame with area >= MIN_AREA_PX.
  3. Normalise: norm(t) = (Mean(t) - bg(t)) / (mean(Mean, frames 1..A) - bg(t)), A = activation frame.
  4. QC: "bright" cells have baseline - bg >= MIN_BASELINE grey levels.
  5. Photobleaching: one rate b per experiment (folder) x reporter, from the pre-activation frames ("preact"),
     fitted across cells ("fit"), or none. Trace is multiplied by exp(b * t).
  6. Per-cell fit from activation (t0 = (A-1) * interval):
       y(t) = y_inf + (1 - y_inf) * exp(-k (t - t0)),  k <= 1/interval
     -> amplitude = y_inf - 1, half-time ln2/k, initial flux = k * amplitude (recommended rate),
        plus a model-free initial slope, plateau-subtracted log fit, gain at a fixed time and AUC.
  7. Per reporter x dose summaries, statistics (Kruskal-Wallis; Mann-Whitney vs control, Holm)
     and figures, with significance stars from the Holm-adjusted Mann-Whitney tests.
"""

# =============================================================================
# PARAMETERS — edit these
# =============================================================================

# --- Which data to analyse ---------------------------------------------------
# Single run: leave BATCH_FOLDERS empty and set DATA_FOLDER.
# Batch: list several folders; each must contain the three input files named below.
DATA_FOLDER = "inputs"
BATCH_FOLDERS = [
    # "/path/to/experiment_1",
    # "/path/to/experiment_2",
]
COMBINE_BATCH = False        # True: also run one pooled analysis of all BATCH_FOLDERS together

# File names inside each folder. Wildcards are allowed (e.g. "*combined_Results*.csv"),
# but each pattern must match exactly one file.
RESULTS_FILE = "combined_Results.csv"
BG_FILE = "mCherry_BGInt_Results.csv"
META_FILE = "fileid_treatment.csv"

# Results go to <folder>/<OUTPUT_SUBFOLDER>; the pooled batch run goes to POOLED_OUTPUT.
OUTPUT_SUBFOLDER = "nct_results"
POOLED_OUTPUT = "nct_results_pooled"

# --- Imaging -----------------------------------------------------------------
INTERVAL_S = 30.0            # seconds between frames
ACTIVATION_FRAME = 4         # last pre-activation frame (1-based); activation time t0 = (frame-1)*interval
DIRECTIONS = {110: "import", 111: "export"}   # reporters not listed are guessed from the data

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

# --- Timecourse plots -------------------------------------------------------
TIMECOURSE_ERRORS = ["sem", "sd"]   # error band(s) to plot; one file per entry ("sem" and/or "sd")

# --- Statistics on plots -----------------------------------------------------
SHOW_STARS = True            # stars from Mann-Whitney vs control, Holm-adjusted
STAR_LEVELS = [(0.001, "***"), (0.01, "**"), (0.05, "*")]

# --- Dot plots (superplots) ---------------------------------------------------
CELL_DOT_COLOR = "#b0b0b0"   # single cells in the background

# --- Per-cell plots ----------------------------------------------------------
SAVE_CELL_PLOTS = False      # one PNG per cell: trace + fitted exponential (slow for many cells)
CELL_PLOTS_BRIGHT_ONLY = True
CELL_PLOT_DPI = 100

# =============================================================================
# END OF PARAMETERS
# =============================================================================

import glob
import os
import warnings
from types import SimpleNamespace

import numpy as np
import pandas as pd
from scipy import stats
from scipy.optimize import curve_fit

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

warnings.filterwarnings("ignore")

DOSE_COLORS = ["#86b6ef", "#5598e7", "#2a78d6", "#1c5cab", "#104281"]
REPORTER_COLORS = {"import": "#2a78d6", "export": "#c2410c"}
# per-dose shades of the reporter colour (light -> dark = low -> high dose) for the dot plots
SHADES = {"import": ["#6fa5e8", "#4a8fe0", "#2a78d6", "#1c5cab", "#104281"],
          "export": ["#ec9460", "#df6a2c", "#c2410c", "#923109", "#632106"]}


def shades(direction, n):
    ramp = SHADES[direction]
    if n <= len(ramp):
        return [ramp[round(i * (len(ramp) - 1) / max(n - 1, 1))] for i in range(n)] if n > 1 else [ramp[2]]
    from matplotlib.colors import LinearSegmentedColormap
    cmap = LinearSegmentedColormap.from_list("r", ramp)
    return [cmap(i / (n - 1)) for i in range(n)]


def settings():
    return SimpleNamespace(
        interval=INTERVAL_S, activation_frame=ACTIVATION_FRAME, min_area=MIN_AREA_PX,
        min_baseline=MIN_BASELINE, bleach=BLEACH_METHOD, bleach_first_frame=BLEACH_FIRST_FRAME,
        slope_frames=SLOPE_FRAMES, gain_time=GAIN_TIME_S, min_r2=MIN_R2, min_amp=MIN_AMPLITUDE)


# --------------------------------------------------------------------------- loading
def find_file(folder, pattern):
    hits = sorted(glob.glob(os.path.join(folder, pattern)))
    if len(hits) != 1:
        raise FileNotFoundError(f"'{pattern}' matched {len(hits)} files in {folder}: {hits}")
    return hits[0]


def load(results_path, bg_path, meta_path):
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

    meta = pd.read_csv(meta_path).rename(columns={"image": "image_base", "time": "treat_time"})
    res = res.merge(meta, on="image_base", how="left", validate="many_to_one")
    if res["treatment"].isna().any():
        raise ValueError("Images missing from the meta file: "
                         + ", ".join(res.loc[res.treatment.isna(), "image_base"].unique()[:5]))
    res["dose_mM"] = res["treatment"].str.extract(r"([\d.]+)\s*mM")[0].astype(float)
    if "replicate" not in res:
        res["replicate"] = "A"
    res["replicate"] = res["replicate"].fillna("A").astype(str)

    bg = pd.read_csv(bg_path)
    bg["Image"] = bg["Label"].str.replace(r"_\d+\.tif$", "", regex=True)
    bg = bg.assign(Slice=bg["Slice"].astype(int) + 1)[["Image", "Slice", "Mean"]].rename(columns={"Mean": "bg_mean"})
    res = res.merge(bg, on=["Image", "Slice"], how="left", validate="many_to_one")
    if res["bg_mean"].isna().any():
        raise ValueError("Background missing for some image/frame pairs")
    return res


def load_folder(folder):
    raw = load(find_file(folder, RESULTS_FILE), find_file(folder, BG_FILE), find_file(folder, META_FILE))
    raw["source"] = os.path.basename(os.path.normpath(folder))
    return raw


def filter_rois(df, args):
    n_frames = df.groupby("Image")["Slice"].transform("max")
    df = df.assign(n_frames=n_frames)
    g = df.groupby(["Image", "roi"])
    keep = (g["Slice"].transform("size") == df["n_frames"]) & (g["Area"].transform("min") >= args.min_area)
    out = df[keep].copy()
    out["cell_id"] = out["reporter"].astype(str) + "_" + out["Image"] + "_" + out["roi"]
    return out, df.groupby(["Image", "roi"]).ngroups, out["cell_id"].nunique()


def normalise(df, args):
    A = args.activation_frame
    df = df.sort_values(["cell_id", "Slice"]).copy()
    df["time"] = (df["Slice"] - 1) * args.interval
    df["Mean_SN"] = df["Mean"] - df["bg_mean"]
    df["Median_SN"] = df["Median"] - df["bg_mean"]
    pre = df[df["Slice"] <= A].groupby("cell_id")
    df["pre_mean"] = df["cell_id"].map(pre["Mean"].mean())
    df["pre_median"] = df["cell_id"].map(pre["Median"].mean())
    df["pre_bg"] = df["cell_id"].map(pre["bg_mean"].mean())
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
        slopes = g.groupby("cell_id").apply(
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
        cells = [(c.time.values, c.norm_mean.values) for _, c in g.groupby("cell_id")]
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
    t0 = (args.activation_frame - 1) * args.interval
    kmax = 1 / args.interval
    t_slope_end = t0 + args.slope_frames * args.interval
    rows = []
    for cid, c in df.groupby("cell_id"):
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
            cell_id=cid, source=c.source.iloc[0], replicate=c.replicate.iloc[0],
            unit=f"{c.source.iloc[0]}_{c.replicate.iloc[0]}", reporter=rep, direction=directions[rep],
            treatment=c.treatment.iloc[0], dose_mM=c.dose_mM.iloc[0], image=c.Image.iloc[0],
            position=c.position.iloc[0], roi=c.roi.iloc[0],
            area_px=c.Area.mean(), baseline_bgsub=c.baseline_bgsub.iloc[0], bright=bool(c.bright.iloc[0]),
            y_inf=yinf, amplitude=amp, k_per_min=k * 60, t_half_s=np.log(2) / k if k > 0 else np.nan, r2=r2,
            flux_pct_per_min=k * amp * 6000,
            flux_uncorrected_pct_per_min=k0 * (yinf0 - 1) * 6000,
            slope_pct_per_min=slope, ln_early_rate_per_min=ln_slope,
            gain_at_t=yc[gain_idx] - 1, auc_frac_min=np.trapezoid(yc[post] - 1, t[post]) / 60,
            k_at_bound=bool(k <= 2e-4 or k >= 0.95 * kmax),
            responder=bool(sign * amp >= args.min_amp),
        ))
    r = pd.DataFrame(rows)
    r["good_fit"] = r.bright & r.responder & (r.r2 >= args.min_r2) & ~r.k_at_bound

    # plateau-subtracted log fit (well-fit cells): ln((y - y_inf)/(1 - y_inf)) = -k (t - t0)
    lnz = {}
    for _, row in r[r.good_fit].iterrows():
        c = df[df.cell_id == row.cell_id]
        t, yc = c.time.values, c.norm_corr.values
        z = (yc - row.y_inf) / (1 - row.y_inf)
        m = (t >= t0) & (t <= t0 + LNZ_WINDOW_S[row.direction]) & (z > 0.05)
        if m.sum() >= 3:
            lnz[row.cell_id] = -np.polyfit(t[m], np.log(z[m]), 1)[0] * 60
    r["k_lnz_per_min"] = r.cell_id.map(lnz)
    return r


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


def summarise(cells):
    b = cells[cells.bright]
    gf = cells[cells.good_fit]
    rows, tests = [], []
    for rep, x in b.groupby("reporter"):
        ref_flux = x[x.dose_mM == x.dose_mM.min()].flux_pct_per_min.median()
        for dose, g in x.groupby("dose_mM"):
            f = gf[(gf.reporter == rep) & (gf.dose_mM == dose)]
            rows.append(dict(
                reporter=rep, direction=g.direction.iloc[0], dose_mM=dose, n_bright=len(g),
                n_fields=g.image.nunique(), n_good_fit=len(f), pct_good_fit=100 * len(f) / len(g), pct_responders=100 * g.responder.mean(),
                flux_median=g.flux_pct_per_min.median(), flux_pct_of_control=100 * g.flux_pct_per_min.median() / ref_flux,
                flux_uncorrected_median=g.flux_uncorrected_pct_per_min.median(),
                slope_median=g.slope_pct_per_min.median(), gain_median=g.gain_at_t.median(),
                auc_median=g.auc_frac_min.median(),
                k_per_min_goodfit=f.k_per_min.median(), t_half_s_goodfit=f.t_half_s.median(),
                amplitude_goodfit=f.amplitude.median(), k_lnz_per_min_goodfit=f.k_lnz_per_min.median()))
        for col in STAT_METRICS:
            xx = gf[gf.reporter == rep] if col in GOODFIT_METRICS else x
            cells_used = "well-fit" if col in GOODFIT_METRICS else "bright"
            doses = sorted(xx.dose_mM.unique())
            groups = [xx[xx.dose_mM == d][col].dropna() for d in doses]
            if len(doses) < 2 or any(len(g) < 2 for g in groups):
                continue
            kw = stats.kruskal(*groups).pvalue
            fields = xx.groupby(["dose_mM", "image"])[col].median().reset_index()
            fgroups = [g[col].dropna() for _, g in fields.groupby("dose_mM")]
            kwf = stats.kruskal(*fgroups).pvalue if all(len(g) >= 2 for g in fgroups) else np.nan
            ctrl = groups[0]
            ps = [stats.mannwhitneyu(ctrl, g).pvalue for g in groups[1:]]
            for d, g, p, pa in zip(doses[1:], groups[1:], ps, holm(ps)):
                tests.append(dict(reporter=rep, metric=col, cells=cells_used, test=f"Mann-Whitney vs {doses[0]:g} mM",
                                  dose_mM=d, n_control=len(ctrl), n_dose=len(g), p=p, p_holm=pa, stars=stars(pa)))
            tests.append(dict(reporter=rep, metric=col, cells=cells_used, test="Kruskal-Wallis (cells)",
                              dose_mM=np.nan, p=kw, stars=stars(kw)))
            tests.append(dict(reporter=rep, metric=col, cells=cells_used, test="Kruskal-Wallis (field medians)",
                              dose_mM=np.nan, p=kwf, stars=stars(kwf)))
        # fraction of bright cells that are well-fit: Fisher's exact test vs control, chi-square across doses
        doses = sorted(x.dose_mM.unique())
        table = [[int(x[x.dose_mM == d].good_fit.sum()), int((~x[x.dose_mM == d].good_fit).sum())] for d in doses]
        if len(doses) > 1:
            ps = [stats.fisher_exact([table[0], t])[1] for t in table[1:]]
            for d, t, p, pa in zip(doses[1:], table[1:], ps, holm(ps)):
                tests.append(dict(reporter=rep, metric="pct_good_fit", cells="bright",
                                  test=f"Fisher exact vs {doses[0]:g} mM", dose_mM=d, n_control=sum(table[0]),
                                  n_dose=sum(t), p=p, p_holm=pa, stars=stars(pa)))
            chi = stats.chi2_contingency(table)[1] if min(map(sum, zip(*table))) > 0 else np.nan
            tests.append(dict(reporter=rep, metric="pct_good_fit", cells="bright", test="Chi-square (doses)",
                              dose_mM=np.nan, p=chi, stars=stars(chi)))
    return pd.DataFrame(rows), pd.DataFrame(tests)


def timecourse_tables(df):
    """Mean / SD / SEM / n per reporter x dose x time, bright cells, bleach-corrected and uncorrected."""
    b = df[df.bright]
    long = (b.groupby(["reporter", "dose_mM", "time"])
             .agg(n=("norm_corr", "size"), mean_corr=("norm_corr", "mean"), sd_corr=("norm_corr", "std"),
                  mean_uncorr=("norm_mean", "mean"), sd_uncorr=("norm_mean", "std"))
             .reset_index())
    long["sem_corr"] = long.sd_corr / np.sqrt(long.n)
    long["sem_uncorr"] = long.sd_uncorr / np.sqrt(long.n)
    long = long[["reporter", "dose_mM", "time", "n", "mean_corr", "sd_corr", "sem_corr",
                 "mean_uncorr", "sd_uncorr", "sem_uncorr"]]
    wide = {}
    for rep, g in long.groupby("reporter"):
        w = g.pivot(index="time", columns="dose_mM", values=["mean_corr", "sd_corr", "sem_corr", "n"])
        w.columns = [f"{d:g}mM_{v.replace('_corr', '')}" for v, d in w.columns]
        order = [f"{d:g}mM_{v}" for d in sorted(g.dose_mM.unique()) for v in ["mean", "sd", "sem", "n"]]
        wide[rep] = w[order].reset_index()
    return long, wide


# --------------------------------------------------------------------------- figures
def style(ax):
    for sp in ["top", "right"]:
        ax.spines[sp].set_visible(False)
    ax.grid(axis="y", color="#eee")
    ax.set_axisbelow(True)


def dose_colors(doses):
    return {d: DOSE_COLORS[min(i, len(DOSE_COLORS) - 1)] for i, d in enumerate(sorted(doses))}


def save(fig, path):
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)


def fig_timecourse(df, rep, direction, args, path, err="sem"):
    t0 = (args.activation_frame - 1) * args.interval
    fig, ax = plt.subplots(figsize=(6, 4.3))
    x = df[(df.reporter == rep) & df.bright]
    cols = dose_colors(x.dose_mM.unique())
    for dose, g in x.groupby("dose_mM"):
        m = g.groupby("time").norm_corr.agg(["mean", {"sd": "std"}.get(err, err)])
        e = m.iloc[:, 1]
        ax.fill_between(m.index, m["mean"] - e, m["mean"] + e, color=cols[dose], alpha=.18, lw=0)
        ax.plot(m.index, m["mean"], color=cols[dose], lw=2, label=f"{dose:g} mM (n={g.cell_id.nunique()})")
    ax.axhline(1, color="#999", lw=.8, ls="--")
    ax.axvspan(0, t0, color="#eee", zorder=0)
    ax.set_title(f"Reporter {rep} ({direction}): mean ± {err.upper()}, bright cells\n(grey = pre-activation)",
                 loc="left", fontsize=10)
    ax.set_xlabel("time (s)")
    ax.set_ylabel("nuclear signal / baseline")
    ax.legend(frameon=False, fontsize=8)
    style(ax)
    save(fig, path)


UNIT_MARKERS = ["o", "s", "^", "D", "v", "P", "X", "<", ">"]


def strip(ax, x, col, doses, colors, fmt="{:.1f}", clip=None):
    """Superplot: every cell as a small grey dot, an unfilled box per dose, and the median of
    each experiment + replicate as a large dot (marker shape = experiment + replicate)."""
    rng = np.random.default_rng(1)
    units = sorted(x.unit.unique())
    markers = {u: UNIT_MARKERS[i % len(UNIT_MARKERS)] for i, u in enumerate(units)}
    for i, d in enumerate(doses):
        color = colors[i]
        g = x[(x.dose_mM == d) & x[col].notna()]
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
    ax.set_xticks(range(len(doses)), [f"{d:g}\n(n={x[x.dose_mM == d][col].notna().sum()})" for d in doses],
                  fontsize=8)
    ax.set_xlabel("2DG (mM)")
    style(ax)


def add_stars(ax, tests, rep, metric, doses):
    if not SHOW_STARS or tests.empty:
        return
    t = tests[(tests.reporter == rep) & (tests.metric == metric)]
    if t.empty:
        return
    lo, hi = ax.get_ylim()
    for i, d in enumerate(doses):
        s = t[t.dose_mM == d]
        if len(s):
            ax.text(i, hi - 0.03 * (hi - lo), s.stars.iloc[0], ha="center", va="top", fontsize=10, color="#222")
    pair = t[t.p_holm.notna()].test.iloc[0].split(" vs ")[0]
    kw = t[t.test.isin(["Kruskal-Wallis (cells)", "Chi-square (doses)"])]
    note = f"Stars: {pair} vs {doses[0]:g} mM, Holm-adjusted"
    if len(kw):
        note += f"   {kw.test.iloc[0].split(' (')[0]} p = {kw.p.iloc[0]:.2g}"
    ax.annotate(note, xy=(0, 0), xycoords="figure fraction", xytext=(6, 4), textcoords="offset points",
                fontsize=7, color="#666")


def fig_metric(cells, tests, rep, direction, metric, title, ylabel, path, goodfit=False,
               fmt="{:.1f}", clip=None, zero_line=True, robust_ylim=False):
    x = cells[cells.good_fit if goodfit else cells.bright]
    x = x[x.reporter == rep]
    doses = sorted(cells[cells.reporter == rep].dose_mM.unique())
    fig, ax = plt.subplots(figsize=(6, 4.8))
    strip(ax, x, metric, doses, shades(direction, len(doses)), fmt=fmt, clip=clip)
    if robust_ylim and x[metric].notna().any():
        lo, hi = np.percentile(x[metric].dropna(), [1, 99])
        ax.set_ylim(min(lo, 0) * 1.1, max(hi, 0) * 1.1)
    lo, hi = ax.get_ylim()
    ax.set_ylim(lo, hi + 0.12 * (hi - lo))            # headroom for stars
    if zero_line:
        ax.axhline(0, color="#999", lw=.8, ls="--")
    ax.set_title(title, loc="left", fontsize=10)
    ax.set_ylabel(ylabel)
    fig.tight_layout(rect=(0, 0.03, 1, 1))
    add_stars(ax, tests, rep, metric, doses)
    fig.savefig(path, dpi=150)
    plt.close(fig)


def fig_pct_goodfit(cells, tests, rep, direction, args, path):
    """% of bright cells that pass the well-fit criteria, per dose."""
    x = cells[cells.bright & (cells.reporter == rep)]
    doses = sorted(x.dose_mM.unique())
    k = np.array([x[x.dose_mM == d].good_fit.sum() for d in doses])
    n = np.array([(x.dose_mM == d).sum() for d in doses])
    pct = 100 * k / n
    fig, ax = plt.subplots(figsize=(6, 4.8))
    ax.bar(range(len(doses)), pct, 0.6, color=REPORTER_COLORS[direction])
    for i in range(len(doses)):
        ax.text(i, pct[i] + 1.5, f"{pct[i]:.0f}%\n({k[i]}/{n[i]})", ha="center", va="bottom", fontsize=8)
    ax.set_xticks(range(len(doses)), [f"{d:g}" for d in doses])
    ax.set_xlabel("2DG (mM)")
    ax.set_ylabel("well-fit cells (% of bright cells)")
    ax.set_ylim(0, 125)
    ax.set_yticks(range(0, 101, 20))
    ax.set_title(f"{rep} ({direction}): well-fit cells\n(R² ≥ {args.min_r2:g}, |amplitude| ≥ {args.min_amp:g} in the "
                 f"expected direction, k not at a bound)", loc="left", fontsize=10)
    style(ax)
    fig.tight_layout(rect=(0, 0.03, 1, 1))
    add_stars(ax, tests, rep, "pct_good_fit", doses)
    fig.savefig(path, dpi=150)
    plt.close(fig)


def fig_logfits(df, cells, rep, direction, args, path):
    t0 = (args.activation_frame - 1) * args.interval
    gf = cells[cells.good_fit].set_index("cell_id")
    fig, axs = plt.subplots(2, 1, figsize=(6, 8.4))
    x = df[(df.reporter == rep) & df.bright & (df.time >= t0)].copy()
    cols = dose_colors(x.dose_mM.unique())
    x["ln"] = np.log(x.norm_corr.clip(lower=1e-3))
    ax = axs[0]
    for dose, g in x.groupby("dose_mM"):
        m = g.groupby("time").ln.agg(["mean", "sem"])
        ax.errorbar(m.index, m["mean"], yerr=m["sem"], fmt="o", ms=3, color=cols[dose], lw=1, label=f"{dose:g} mM")
        e = m[m.index <= t0 + args.slope_frames * args.interval]
        p = np.polyfit(e.index, e["mean"], 1)
        tt = np.linspace(t0, m.index.max(), 30)
        ax.plot(tt, np.polyval(p, tt), color=cols[dose], lw=1, ls="--")
    ax.axhline(0, color="#999", lw=.8)
    ax.set_title(f"{rep} ({direction}): ln(norm); dashed = early-window fit", loc="left", fontsize=10)
    ax.legend(frameon=False, fontsize=8)
    style(ax)
    ax = axs[1]
    y = x[x.cell_id.isin(gf.index) & (x.time <= t0 + LNZ_WINDOW_S[direction])].copy()
    y["yinf"] = y.cell_id.map(gf.y_inf)
    y["z"] = (y.norm_corr - y.yinf) / (1 - y.yinf)
    y = y[y.z > 0.05]
    y["lnz"] = np.log(y.z)
    for dose, g in y.groupby("dose_mM"):
        m = g.groupby("time").lnz.agg(["mean", "sem", "size"])
        m = m[m["size"] >= 5]
        if len(m) < 2:
            continue
        ax.errorbar(m.index, m["mean"], yerr=m["sem"], fmt="o", ms=3, color=cols[dose], lw=1,
                    label=f"{dose:g} mM (n={g.cell_id.nunique()})")
        p = np.polyfit(m.index, m["mean"], 1)
        ax.plot(m.index, np.polyval(p, m.index), color=cols[dose], lw=1.2)
    ax.axhline(0, color="#999", lw=.8)
    ax.set_title(f"{rep}: ln((norm − plateau)/(1 − plateau)), well-fit cells; slope = −k", loc="left", fontsize=10)
    ax.set_xlabel("time (s)")
    ax.legend(frameon=False, fontsize=8)
    style(ax)
    save(fig, path)


def fig_bleach(bleach, cells, reps, directions, args, path):
    fig, axs = plt.subplots(1, 1 + len(reps), figsize=(5 + 5.5 * len(reps), 4.4), squeeze=False)
    ax = axs[0][0]
    keys = sorted(bleach, key=lambda k: (k[1], k[0]))
    multi = len({k[0] for k in keys}) > 1
    for i, (src, rep) in enumerate(keys):
        bb = bleach[(src, rep)]
        c = REPORTER_COLORS[directions[rep]]
        yerr = None if np.isnan(bb["ci_lo"]) else [[bb["b_raw"] * 6000 - bb["ci_lo"] * 6000],
                                                     [bb["ci_hi"] * 6000 - bb["b_raw"] * 6000]]
        ax.errorbar(i, bb["b_raw"] * 6000, yerr=yerr, fmt="o", color=c, ms=8, capsize=4)
        ax.text(i + .08, bb["b_raw"] * 6000, f"{bb['b_raw'] * 6000:.2f}", fontsize=8, va="center")
    ax.axhline(0, color="#999", lw=.8)
    ax.set_xticks(range(len(keys)), [f"{r}\n{src}" if multi else f"{r} ({directions[r]})" for src, r in keys],
                  fontsize=8 if multi else 10)
    ax.set_xlim(-.5, len(keys) - .5)
    ax.set_ylabel("bleaching rate b (% per min)")
    ax.set_title(f"Bleaching rate per {'experiment × ' if multi else ''}reporter ({args.bleach})", loc="left", fontsize=10)
    style(ax)
    b = cells[cells.bright]
    for ax, rep in zip(axs[0][1:], reps):
        x = b[b.reporter == rep]
        doses = sorted(x.dose_mM.unique())
        w = .38
        for k, (col, lab, alpha) in enumerate([("flux_uncorrected_pct_per_min", "uncorrected", .35),
                                               ("flux_pct_per_min", "bleach-corrected", .9)]):
            med = x.groupby("dose_mM")[col].median().reindex(doses)
            pct = med / med.iloc[0] * 100
            ax.bar(np.arange(len(doses)) + (k - .5) * w, pct, w * .95, color=REPORTER_COLORS[directions[rep]],
                   alpha=alpha, label=lab)
            for i, v in enumerate(pct):
                ax.text(i + (k - .5) * w, v + 1.5, f"{v:.0f}", ha="center", fontsize=7.5)
        ax.set_xticks(range(len(doses)), [f"{d:g}" for d in doses])
        ax.set_xlabel("2DG (mM)")
        ax.set_ylim(0, 125)
        ax.set_title(f"{rep}: median initial flux, % of control", loc="left", fontsize=10)
        ax.legend(frameon=False, fontsize=8)
        style(ax)
    save(fig, path)


def cell_plots(df, cells, args, outdir):
    """One PNG per cell: normalised and bleach-corrected trace with the fitted exponential."""
    t0 = (args.activation_frame - 1) * args.interval
    info = cells.set_index("cell_id")
    if CELL_PLOTS_BRIGHT_ONLY:
        df = df[df.bright]
    n = 0
    for cid, c in df.groupby("cell_id"):
        r = info.loc[cid]
        sub = os.path.join(outdir, f"{r.reporter}_{r.dose_mM:g}mM")
        os.makedirs(sub, exist_ok=True)
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
        ax.set_title(f"{os.path.basename(r.image)} ROI {r.roi}  ({r.reporter}, {r.dose_mM:g} mM)\n"
                     f"flux {r.flux_pct_per_min:.1f} %/min · amp {r.amplitude:.2f} · "
                     f"t½ {r.t_half_s:.0f} s · R² {r.r2:.2f} · {verdict}", loc="left", fontsize=8)
        ax.set_xlabel("time (s)")
        ax.set_ylabel("nuclear signal / baseline")
        ax.legend(frameon=False, fontsize=7)
        style(ax)
        fig.tight_layout()
        fig.savefig(os.path.join(sub, f"{r.image.replace('/', '_')}_ROI{r.roi}.png"), dpi=CELL_PLOT_DPI)
        plt.close(fig)
        n += 1
    return n


# --------------------------------------------------------------------------- run
def analyse(raw, out, args):
    os.makedirs(out, exist_ok=True)
    df, n_rois, n_kept = filter_rois(raw, args)
    df = normalise(df, args)
    print(f"ROIs: {n_rois} total, {n_kept} kept (full track, area >= {args.min_area:g} px), "
          f"{df[df.bright].cell_id.nunique()} bright")

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
    pd.DataFrame([dict(source=src, reporter=r, method=args.bleach, b_pct_per_min_estimate=v["b_raw"] * 6000,
                       b_pct_per_min_applied=v["b"] * 6000, ci_lo=v["ci_lo"] * 6000, ci_hi=v["ci_hi"] * 6000,
                       n_cells=v["n"]) for (src, r), v in bleach.items()]).to_csv(
        os.path.join(out, "bleach_rates.csv"), index=False)

    cells = per_cell(df, args, directions)
    summary, tests = summarise(cells)
    tc_long, tc_wide = timecourse_tables(df)

    keep_cols = ["source", "row_id", "cell_id", "Image", "roi", "Slice", "time", "reporter", "treatment",
                 "dose_mM", "treat_time", "pyruvate", "replicate", "position", "Area", "Mean", "StdDev", "Min",
                 "Max", "Median", "IntDen", "bg_mean", "Mean_SN", "Median_SN", "pre_mean", "pre_median",
                 "baseline_bgsub", "bright", "norm_mean", "norm_median", "norm_corr"]
    df[[k for k in keep_cols if k in df]].to_csv(os.path.join(out, "timecourses_long.csv"), index=False)
    tc_long.to_csv(os.path.join(out, "timecourse_mean_by_dose.csv"), index=False)
    cells.to_csv(os.path.join(out, "rates_per_cell.csv"), index=False)
    rep_cols = ["flux_pct_per_min", "slope_pct_per_min", "gain_at_t", "auc_frac_min"]
    fit_cols = ["k_per_min", "t_half_s", "amplitude"]
    b_ = cells[cells.bright].groupby(["reporter", "dose_mM", "source", "replicate"])
    f_ = cells[cells.good_fit].groupby(["reporter", "dose_mM", "source", "replicate"])
    pd.concat([b_.size().rename("n_bright"), b_[rep_cols].median(),
               f_.size().rename("n_good_fit"), f_[fit_cols].median().add_suffix("_goodfit")], axis=1
              ).reset_index().to_csv(os.path.join(out, "replicate_medians.csv"), index=False)
    summary.to_csv(os.path.join(out, "summary_by_dose.csv"), index=False)
    tests.to_csv(os.path.join(out, "statistics.csv"), index=False)

    for rep in reps:
        d = directions[rep]
        tc_wide[rep].to_csv(os.path.join(out, f"timecourse_{rep}_wide.csv"), index=False)
        for err in TIMECOURSE_ERRORS:
            suffix = "" if err == "sem" else f"_{err}"
            fig_timecourse(df, rep, d, args, os.path.join(out, f"timecourse_{rep}{suffix}.png"), err=err)
        ref = cells[cells.bright & (cells.reporter == rep)]
        ref = ref[ref.dose_mM == ref.dose_mM.min()].flux_pct_per_min.median()
        fig_metric(cells, tests, rep, d, "flux_pct_per_min",
                   f"{rep} ({d}): initial flux k × amplitude, all bright cells\n(control median {ref:.1f})",
                   "% of baseline per min", os.path.join(out, f"flux_{rep}.png"))
        fig_metric(cells, tests, rep, d, "slope_pct_per_min",
                   f"{rep} ({d}): model-free initial slope, all bright cells\n(linear fit, frames "
                   f"{ACTIVATION_FRAME}–{ACTIVATION_FRAME + SLOPE_FRAMES})", "% of baseline per min",
                   os.path.join(out, f"slope_{rep}.png"))
        fig_metric(cells, tests, rep, d, "t_half_s",
                   f"{rep} ({d}): half-time of fitted exponential, well-fit cells", "t½ (s)",
                   os.path.join(out, f"half_time_{rep}.png"), goodfit=True, fmt="{:.0f}", clip=(0, 600),
                   zero_line=False)
        fig_metric(cells, tests, rep, d, "k_per_min",
                   f"{rep} ({d}): rate constant k of fitted exponential, well-fit cells", "k (per min)",
                   os.path.join(out, f"k_{rep}.png"), goodfit=True, fmt="{:.2f}", zero_line=False)
        fig_pct_goodfit(cells, tests, rep, d, args, os.path.join(out, f"pct_well_fit_{rep}.png"))
        fig_logfits(df, cells, rep, d, args, os.path.join(out, f"log_fits_{rep}.png"))
    fig_bleach(bleach, cells, reps, directions, args, os.path.join(out, "bleach_correction.png"))

    if SAVE_CELL_PLOTS:
        n = cell_plots(df, cells, args, os.path.join(out, "cell_plots"))
        print(f"Saved {n} per-cell plots")

    pd.set_option("display.width", 200)
    print("Bleaching (%/min):", {f"{src} {r}": round(v["b"] * 6000, 3) for (src, r), v in bleach.items()})
    print(summary[["reporter", "direction", "dose_mM", "n_bright", "pct_responders", "flux_median",
                   "flux_pct_of_control", "slope_median", "t_half_s_goodfit"]].round(2).to_string(index=False))
    print(f"Outputs written to {os.path.abspath(out)}\n")


def main():
    args = settings()
    if not BATCH_FOLDERS:
        print(f"=== {DATA_FOLDER}")
        analyse(load_folder(DATA_FOLDER), os.path.join(DATA_FOLDER, OUTPUT_SUBFOLDER), args)
        return
    loaded = []
    for folder in BATCH_FOLDERS:
        print(f"=== {folder}")
        try:
            raw = load_folder(folder)
            analyse(raw, os.path.join(folder, OUTPUT_SUBFOLDER), args)
            loaded.append(raw)
        except Exception as e:
            print(f"Skipped {folder}: {e}\n")
    if COMBINE_BATCH and len(loaded) > 1:
        print("=== pooled")
        pooled = []
        for raw in loaded:      # keep image names unique across folders
            raw = raw.copy()
            raw["Image"] = raw["source"] + "/" + raw["Image"]
            pooled.append(raw)
        analyse(pd.concat(pooled, ignore_index=True), POOLED_OUTPUT, args)


if __name__ == "__main__":
    main()
