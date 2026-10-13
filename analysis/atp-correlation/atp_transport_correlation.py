#!/usr/bin/env python3
"""
Correlate cellular ATP/ADP (ratiometric sensor, 488/405 excitation) with nuclear import and export
rates from nct_analysis.py.

HOW TO USE
  1. Edit the PARAMETERS block below (file names, folders, settings).
  2. Run the file:  python atp_transport_correlation.py   (or press Run in your editor).
  Requires: numpy, pandas, matplotlib.

INPUT FILES (single run: ATP_FILE + RATES_FOLDER; batch: PAIRS_FILE or BATCH_RUNS pair them per experiment)
  ATP_FILE      Per-cell ATP/ADP measurements (one row per cell; e.g. all_images_measurements.csv)
                with a treatment column ("0mM 2DG", ...) and a ratio column (ATP_COLUMN).
  RATES_FOLDER  An nct_analysis.py output folder (nct_results/) with rates_per_cell.csv.

HOW THEY ARE MATCHED
  MATCH_MODE = "condition" (default)
    ATP and transport were measured in different cells and dishes, so they are paired by condition
    (2DG dose). For each dose: median ATP/ADP and median rate, with bootstrap 95% CIs over cells.
    Across doses: Pearson r, Spearman rho and a linear fit of rate vs ATP/ADP.
    p-values are exact permutation tests over the doses. Only doses present in both files are used;
    with 4 shared doses the smallest possible two-sided p is 2/24 = 0.083, with 5 it is 2/120 = 0.017. r_CI is a bootstrap CI that resamples cells within each
    dose, i.e. how stable r is given cell-to-cell variability (not dose-to-dose).
  MATCH_MODE = "cell"
    Only when the same cells carry both measurements. Rows are joined on ATP_CELL_KEY / RATES_CELL_KEY
    and correlated per cell (per reporter, and per dose if PER_DOSE_CELL_STATS).

OUTPUTS (in OUTPUT_FOLDER, or OUTPUT_FOLDER/<batch name> for each pair in batch mode)
  atp_by_dose.csv          ATP/ADP per dose: n cells, n images, median, 95% CI, mean, SEM
  rates_vs_atp_by_dose.csv One row per reporter x dose: ATP/ADP and every rate metric (median, CI)
  correlations.csv         One row per reporter x rate metric: r, rho, permutation p, slope, r CI
  atp_by_dose.png          ATP/ADP per dose (cells + median)
  rate_vs_atp.png          Rate vs ATP/ADP per reporter (dose medians with 95% CI, linear fit)
  rate_vs_atp_all_metrics.png  Every metric in RATE_METRICS vs ATP/ADP
  (cell mode) cell_matched.csv, cell_correlations.csv, rate_vs_atp_cells.png
"""

# =============================================================================
# PARAMETERS — edit these
# =============================================================================

# --- Which data to analyse ---------------------------------------------------
# Single run: leave PAIRS_FILE empty and BATCH_RUNS empty, and set ATP_FILE and RATES_FOLDER.
ATP_FILE = "inputs/all_images_measurements.csv"
RATES_FOLDER = "../nct-analysis/inputs/nct_results"

# Batch, option 1: a pairs table (same format as folder_pairs.csv used by plot_ratio_vs_k.py).
#   batch,imaging_folders,kinetics_folders,include
#   20260925_MN,20260925_MN_92,"20260925_MN_110,111",yes
# imaging_folders hold the ATP file (ATP_FILE_NAME); kinetics_folders are the NCT data folders, read
# from their NCT_SUBFOLDER (nct_results/) written by nct_analysis.py. Folder names are relative to
# PAIRS_ROOT; several folders in one cell are separated by ";" and pooled. include = "no" skips a row.
PAIRS_FILE = ""                               # e.g. "folder_pairs.csv"; "" = not used
PAIRS_ROOT = ""                               # folder holding the experiment folders ("" = the pairs file's folder)
NCT_SUBFOLDER = "nct_results"

# Batch, option 2: list the pairs here. Same idea; "atp" and "nct" are paths (absolute, or relative
# to this script). "atp" may be the CSV or the folder holding it; "nct" may be the nct_results folder
# or the data folder containing it. Use a list for several folders. "name" = output subfolder.
BATCH_RUNS = [
    # {"name": "20260925", "nct": "/path/to/20260925_MN_110,111", "atp": "/path/to/20260925_MN_92"},
]
ATP_FILE_NAME = "all_images_measurements.csv"  # file looked for when an ATP path is a folder
RATES_FILE = "rates_per_cell.csv"             # file looked for inside each NCT results folder
OUTPUT_FOLDER = "atp_results"                 # batch: one subfolder per pair

# --- ATP/ADP measurement -----------------------------------------------------
ATP_COLUMN = "whole_cell_SN_ratio_488_405"   # also: bright_SN_ratio_488_405, nucleus_SN_ratio_488_405, *_median_*
ATP_LABEL = "ATP/ADP (488/405)"
ATP_TREATMENT_COLUMN = "treatment"            # values like "50mM 2DG"; dose is read from the number
ATP_IMAGE_COLUMN = "image_name"               # used to count images per dose
ATP_FILTERS = {                               # keep rows where column == value (empty dict = keep all)
    "cell_type": "MN",
    "used_whole_cell_fallback": False,
}
ATP_MIN_AREA_PX = 500                         # cell_area_px limits (drops debris and merged cells)
ATP_MAX_AREA_PX = 20000
ATP_AREA_COLUMN = "cell_area_px"

# --- Transport rates ---------------------------------------------------------
RATES_BRIGHT_ONLY = True                      # use only "bright" cells (as in nct_analysis summaries)
# metric column -> (label, use well-fit cells only)
RATE_METRICS = {
    "flux_pct_per_min": ("Initial flux (% baseline/min)", False),
    "slope_pct_per_min": ("Initial slope (% baseline/min)", False),
    "gain_at_t": ("Gain at 600 s (norm - 1)", False),
    "auc_frac_min": ("AUC (fraction x min)", False),
    "k_per_min": ("k (1/min), good fits", True),
    "t_half_s": ("t1/2 (s), good fits", True),
}
MAIN_METRIC = "flux_pct_per_min"              # metric shown in rate_vs_atp.png
ABSOLUTE_EXPORT = True                        # plot/correlate |rate| for export (export flux is negative)

# --- Matching ----------------------------------------------------------------
MATCH_MODE = "condition"                      # "condition" (by 2DG dose) or "cell"
ATP_CELL_KEY = "cell_id"                      # cell mode only: column in ATP_FILE ...
RATES_CELL_KEY = "cell_id"                    # ... and in rates_per_cell.csv with the same values
PER_DOSE_CELL_STATS = True                    # cell mode: also correlate within each dose

# --- Statistics --------------------------------------------------------------
N_BOOT = 2000                                 # bootstrap resamples for CIs
RANDOM_SEED = 1

# --- Figures -----------------------------------------------------------------
FIG_DPI = 200

# =============================================================================
# END OF PARAMETERS
# =============================================================================

import itertools
import os
import re
import warnings

import numpy as np
import pandas as pd

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

warnings.filterwarnings("ignore")

DOSE_COLORS = ["#86b6ef", "#5598e7", "#2a78d6", "#1c5cab", "#104281"]
REPORTER_COLORS = {"import": "#2a78d6", "export": "#c2410c"}
HERE = os.path.dirname(os.path.abspath(__file__))


def here(path):
    return path if os.path.isabs(path) else os.path.join(HERE, path)


# --------------------------------------------------------------------------- statistics
def rank(x):
    return pd.Series(x).rank().to_numpy()


def pearson(x, y):
    x, y = np.asarray(x, float), np.asarray(y, float)
    if len(x) < 3 or np.std(x) == 0 or np.std(y) == 0:
        return np.nan
    return float(np.corrcoef(x, y)[0, 1])


def spearman(x, y):
    return pearson(rank(x), rank(y))


def perm_p(x, y, fn):
    """Two-sided permutation p for correlation fn(x, y); exact for n <= 8, else 20000 shuffles."""
    x, y = np.asarray(x, float), np.asarray(y, float)
    r0 = fn(x, y)
    if np.isnan(r0):
        return np.nan
    n = len(x)
    if n <= 8:
        perms = (np.array(p) for p in itertools.permutations(range(n)))
        rs = np.array([fn(x, y[p]) for p in perms])
    else:
        rng = np.random.default_rng(RANDOM_SEED)
        rs = np.array([fn(x, rng.permutation(y)) for _ in range(20000)])
    return float(np.mean(np.abs(rs) >= abs(r0) - 1e-12))


def boot_median(v, rng):
    v = np.asarray(v, float)
    v = v[np.isfinite(v)]
    if len(v) == 0:
        return np.nan, np.nan, np.nan, np.empty(0)
    meds = np.median(rng.choice(v, (N_BOOT, len(v))), axis=1)
    return float(np.median(v)), *np.percentile(meds, [2.5, 97.5]), meds


def dose_from_treatment(s):
    m = re.search(r"([\d.]+)\s*mM", str(s))
    return float(m.group(1)) if m else np.nan


# --------------------------------------------------------------------------- loading
def as_list(x):
    return list(x) if isinstance(x, (list, tuple)) else [x]


def resolve_atp(path):
    path = here(path)
    return os.path.join(path, ATP_FILE_NAME) if os.path.isdir(path) else path


def resolve_rates(folder):
    folder = here(folder)
    if not os.path.exists(os.path.join(folder, RATES_FILE)) and os.path.isdir(os.path.join(folder, NCT_SUBFOLDER)):
        folder = os.path.join(folder, NCT_SUBFOLDER)
    return os.path.join(folder, RATES_FILE)


def read_pairs():
    """[(name, [atp paths], [nct folders])] from PAIRS_FILE."""
    path = here(PAIRS_FILE)
    root = here(PAIRS_ROOT) if PAIRS_ROOT else os.path.dirname(path)
    df = pd.read_csv(path, dtype=str).fillna("")
    runs = []
    for _, row in df.iterrows():
        if row.get("include", "yes").strip().lower() in ("no", "n", "false", "0"):
            print(f"skipping pair {row['batch']} (include = {row['include']})")
            continue
        split = lambda col: [os.path.join(root, f.strip()) for f in row[col].split(";") if f.strip()]
        atp, nct = split("imaging_folders"), split("kinetics_folders")
        if not atp or not nct:
            raise SystemExit(f"{path}: pair {row['batch']} needs imaging and kinetics folders")
        runs.append((row["batch"].strip() or os.path.basename(atp[0]), atp, nct))
    return runs


def load_atp(atp_paths):
    files = [resolve_atp(p) for p in as_list(atp_paths)]
    for f in files:
        print(f"ATP file: {f}")
    atp = pd.concat([pd.read_csv(f) for f in files], ignore_index=True)
    n0 = len(atp)
    for col, val in ATP_FILTERS.items():
        atp = atp[atp[col] == val]
    if ATP_AREA_COLUMN in atp:
        atp = atp[atp[ATP_AREA_COLUMN].between(ATP_MIN_AREA_PX, ATP_MAX_AREA_PX)]
    atp = atp[np.isfinite(atp[ATP_COLUMN]) & (atp[ATP_COLUMN] > 0)].copy()
    atp["dose_mM"] = atp[ATP_TREATMENT_COLUMN].map(dose_from_treatment)
    print(f"ATP: kept {len(atp)} of {n0} cells after filters")
    return atp


def load_rates(folders):
    files = [resolve_rates(f) for f in as_list(folders)]
    for f in files:
        print(f"Rates file: {f}")
    rates = pd.concat([pd.read_csv(f) for f in files], ignore_index=True)
    if RATES_BRIGHT_ONLY:
        rates = rates[rates["bright"].astype(str) == "True"]
    rates = rates.copy()
    if ABSOLUTE_EXPORT:
        exp = rates["direction"] == "export"
        for m in ("flux_pct_per_min", "slope_pct_per_min", "gain_at_t", "auc_frac_min"):
            if m in rates:
                rates.loc[exp, m] = rates.loc[exp, m].abs()
    return rates


def metric_values(df, metric):
    _, goodfit_only = RATE_METRICS[metric]
    if goodfit_only:
        df = df[df["good_fit"].astype(str) == "True"]
    return df[metric].to_numpy(float)


# --------------------------------------------------------------------------- condition mode
def atp_by_dose(atp, rng):
    rows, boots = [], {}
    for dose, g in atp.groupby("dose_mM"):
        med, lo, hi, meds = boot_median(g[ATP_COLUMN], rng)
        boots[dose] = meds
        v = g[ATP_COLUMN]
        rows.append(dict(dose_mM=dose, n_cells=len(g), n_images=g[ATP_IMAGE_COLUMN].nunique(),
                         atp_median=med, atp_ci_low=lo, atp_ci_high=hi,
                         atp_mean=v.mean(), atp_sem=v.sem()))
    return pd.DataFrame(rows), boots


def condition_mode(atp, rates, out, rng):
    atp_tab, atp_boot = atp_by_dose(atp, rng)
    atp_tab.to_csv(os.path.join(out, "atp_by_dose.csv"), index=False)

    rows, corr_rows = [], []
    for (rep, direction), rg in rates.groupby(["reporter", "direction"]):
        doses = sorted(set(rg["dose_mM"]) & set(atp_tab["dose_mM"]))
        for metric in RATE_METRICS:
            if metric not in rg:
                continue
            med, rboot = [], []
            for d in doses:
                m, lo, hi, meds = boot_median(metric_values(rg[rg["dose_mM"] == d], metric), rng)
                med.append(m)
                rboot.append(meds)
                a = atp_tab.set_index("dose_mM").loc[d]
                rows.append(dict(reporter=rep, direction=direction, metric=metric, dose_mM=d,
                                 atp_median=a.atp_median, atp_ci_low=a.atp_ci_low,
                                 atp_ci_high=a.atp_ci_high, n_atp_cells=a.n_cells,
                                 rate_median=m, rate_ci_low=lo, rate_ci_high=hi,
                                 n_rate_cells=int(np.isfinite(metric_values(
                                     rg[rg["dose_mM"] == d], metric)).sum())))
            x = atp_tab.set_index("dose_mM").loc[doses, "atp_median"].to_numpy()
            y = np.array(med)
            ok = np.isfinite(x) & np.isfinite(y)
            x, y = x[ok], y[ok]
            slope, icpt = (np.polyfit(x, y, 1) if len(x) >= 2 else (np.nan, np.nan))
            # bootstrap r: resample cells within each dose (ATP and rate independently)
            bd = [d for d, k in zip(doses, ok) if k]
            br = [b for b, k in zip(rboot, ok) if k]
            if len(bd) >= 3 and all(len(b) for b in br):
                bx = np.column_stack([atp_boot[d] for d in bd])
                by = np.column_stack(br)
                rs = np.array([pearson(bx[i], by[i]) for i in range(N_BOOT)])
                r_lo, r_hi = np.nanpercentile(rs, [2.5, 97.5])
            else:
                r_lo = r_hi = np.nan
            corr_rows.append(dict(reporter=rep, direction=direction, metric=metric,
                                  n_doses=len(x), pearson_r=pearson(x, y),
                                  pearson_perm_p=perm_p(x, y, pearson),
                                  pearson_r_boot_ci_low=r_lo, pearson_r_boot_ci_high=r_hi,
                                  spearman_rho=spearman(x, y),
                                  spearman_perm_p=perm_p(x, y, spearman),
                                  slope_rate_per_atp=slope, intercept=icpt))
    tab = pd.DataFrame(rows)
    corr = pd.DataFrame(corr_rows)
    tab.to_csv(os.path.join(out, "rates_vs_atp_by_dose.csv"), index=False)
    corr.to_csv(os.path.join(out, "correlations.csv"), index=False)

    plot_atp(atp, atp_tab, out)
    plot_rate_vs_atp(tab, corr, [MAIN_METRIC], os.path.join(out, "rate_vs_atp.png"))
    plot_rate_vs_atp(tab, corr, list(RATE_METRICS), os.path.join(out, "rate_vs_atp_all_metrics.png"))
    print(corr[corr["metric"] == MAIN_METRIC][["direction", "pearson_r", "pearson_perm_p",
                                                "spearman_rho", "spearman_perm_p"]].to_string(index=False))


# --------------------------------------------------------------------------- cell mode
def cell_mode(atp, rates, out):
    m = rates.merge(atp[[ATP_CELL_KEY, ATP_COLUMN]], left_on=RATES_CELL_KEY,
                    right_on=ATP_CELL_KEY, how="inner")
    if m.empty:
        raise ValueError(f"No cells matched between {ATP_CELL_KEY} (ATP) and {RATES_CELL_KEY} (rates). "
                         "Use MATCH_MODE = 'condition' when ATP and transport are different cells.")
    m.to_csv(os.path.join(out, "cell_matched.csv"), index=False)
    rows = []
    groups = [("all", None)] + ([(d, d) for d in sorted(m["dose_mM"].unique())] if PER_DOSE_CELL_STATS else [])
    for (rep, direction), rg in m.groupby(["reporter", "direction"]):
        for metric in RATE_METRICS:
            for name, d in groups:
                g = rg if d is None else rg[rg["dose_mM"] == d]
                _, goodfit_only = RATE_METRICS[metric]
                if goodfit_only:
                    g = g[g["good_fit"].astype(str) == "True"]
                g = g[np.isfinite(g[metric]) & np.isfinite(g[ATP_COLUMN])]
                x, y = g[ATP_COLUMN].to_numpy(), g[metric].to_numpy()
                rows.append(dict(reporter=rep, direction=direction, metric=metric, dose_mM=name,
                                 n_cells=len(g), pearson_r=pearson(x, y),
                                 pearson_perm_p=perm_p(x, y, pearson),
                                 spearman_rho=spearman(x, y),
                                 spearman_perm_p=perm_p(x, y, spearman)))
    pd.DataFrame(rows).to_csv(os.path.join(out, "cell_correlations.csv"), index=False)

    fig, axes = plt.subplots(1, m["direction"].nunique(), figsize=(5 * m["direction"].nunique(), 4), squeeze=False)
    doses = sorted(m["dose_mM"].unique())
    for ax, (direction, g) in zip(axes[0], m.groupby("direction")):
        for i, d in enumerate(doses):
            gd = g[g["dose_mM"] == d]
            ax.scatter(gd[ATP_COLUMN], gd[MAIN_METRIC], s=10, alpha=0.7,
                       color=DOSE_COLORS[i % len(DOSE_COLORS)], label=f"{d:g} mM")
        ax.set_xlabel(ATP_LABEL)
        ax.set_ylabel(RATE_METRICS[MAIN_METRIC][0])
        ax.set_title(direction)
        ax.legend(frameon=False, fontsize=8, title="2DG")
    fig.tight_layout()
    fig.savefig(os.path.join(out, "rate_vs_atp_cells.png"), dpi=FIG_DPI)
    plt.close(fig)


# --------------------------------------------------------------------------- figures
def style(ax):
    ax.spines[["top", "right"]].set_visible(False)


def plot_atp(atp, atp_tab, out):
    fig, ax = plt.subplots(figsize=(4.5, 4))
    rng = np.random.default_rng(RANDOM_SEED)
    for i, r in atp_tab.reset_index(drop=True).iterrows():
        v = atp.loc[atp["dose_mM"] == r.dose_mM, ATP_COLUMN]
        ax.scatter(i + rng.uniform(-0.2, 0.2, len(v)), v, s=6, alpha=0.35,
                   color=DOSE_COLORS[i % len(DOSE_COLORS)], linewidths=0)
        ax.errorbar(i, r.atp_median, yerr=[[r.atp_median - r.atp_ci_low], [r.atp_ci_high - r.atp_median]],
                    fmt="o", color="black", ms=5, capsize=4, zorder=3)
        ax.text(i, v.max(), f"n={int(r.n_cells)}", ha="center",
                va="bottom", fontsize=7, color="0.4")
    ax.set_xticks(range(len(atp_tab)), [f"{d:g}" for d in atp_tab["dose_mM"]])
    ax.set_xlabel("2DG (mM)")
    ax.set_ylabel(ATP_LABEL)
    ax.set_title("ATP/ADP by dose (median, 95% CI)", fontsize=10)
    style(ax)
    fig.tight_layout()
    fig.savefig(os.path.join(out, "atp_by_dose.png"), dpi=FIG_DPI)
    plt.close(fig)


def plot_rate_vs_atp(tab, corr, metrics, path):
    directions = list(dict.fromkeys(tab["direction"]))
    metrics = [m for m in metrics if m in set(tab["metric"])]
    fig, axes = plt.subplots(len(metrics), len(directions), squeeze=False,
                             figsize=(4.6 * len(directions), 3.8 * len(metrics)))
    for row, metric in enumerate(metrics):
        for col, direction in enumerate(directions):
            ax = axes[row, col]
            t = tab[(tab["metric"] == metric) & (tab["direction"] == direction)].sort_values("dose_mM")
            c = corr[(corr["metric"] == metric) & (corr["direction"] == direction)].iloc[0]
            color = REPORTER_COLORS.get(direction, "0.3")
            for i, r in enumerate(t.itertuples()):
                ax.errorbar(r.atp_median, r.rate_median,
                            xerr=[[r.atp_median - r.atp_ci_low], [r.atp_ci_high - r.atp_median]],
                            yerr=[[r.rate_median - r.rate_ci_low], [r.rate_ci_high - r.rate_median]],
                            fmt="o", ms=6, color=DOSE_COLORS[i % len(DOSE_COLORS)],
                            ecolor="0.6", capsize=3, zorder=3)
                ax.annotate(f"{r.dose_mM:g} mM", (r.atp_median, r.rate_median), xytext=(5, 4),
                            textcoords="offset points", fontsize=7, color="0.35")
            if np.isfinite(c.slope_rate_per_atp):
                xs = np.linspace(t["atp_ci_low"].min(), t["atp_ci_high"].max(), 50)
                ax.plot(xs, c.slope_rate_per_atp * xs + c.intercept, color=color, lw=1.2, alpha=0.8)
            ax.set_title(f"{direction} — r = {c.pearson_r:.2f} (p = {c.pearson_perm_p:.3f}), "
                         f"rho = {c.spearman_rho:.2f}", fontsize=9)
            ax.set_xlabel(ATP_LABEL)
            label = RATE_METRICS[metric][0]
            if ABSOLUTE_EXPORT and direction == "export" and metric in (
                    "flux_pct_per_min", "slope_pct_per_min", "gain_at_t", "auc_frac_min"):
                label = f"|{label}|"
            ax.set_ylabel(label)
            style(ax)
    fig.suptitle("Transport rate vs ATP/ADP (dose medians, 95% bootstrap CI; permutation p over doses)",
                 fontsize=9)
    fig.tight_layout()
    fig.savefig(path, dpi=FIG_DPI)
    plt.close(fig)


# --------------------------------------------------------------------------- main
def run(atp_paths, rates_folders, out):
    os.makedirs(out, exist_ok=True)
    rng = np.random.default_rng(RANDOM_SEED)
    atp = load_atp(atp_paths)
    rates = load_rates(rates_folders)
    print(f"Rates: {len(rates)} cells")
    if MATCH_MODE == "cell":
        cell_mode(atp, rates, out)
    else:
        condition_mode(atp, rates, out, rng)
    print(f"Wrote results to {out}")


if __name__ == "__main__":
    runs = read_pairs() if PAIRS_FILE else [
        (r.get("name") or f"run{i}", r["atp"], r["nct"]) for i, r in enumerate(BATCH_RUNS, 1)]
    if runs:
        for name, atp_paths, nct_folders in runs:
            print(f"\n=== {name} ===")
            run(atp_paths, nct_folders, os.path.join(here(OUTPUT_FOLDER), name))
    else:
        run(ATP_FILE, RATES_FOLDER, here(OUTPUT_FOLDER))
