"""Pair imaging and kinetics experiments from the same cell batch and plot
the fitted rate constant against the normalised 488/405 ATP/ADP ratio.

Folder pairing
--------------
Each folder is one experiment, named <date>_<cell>_<reporter ids>, e.g.
    20260925_MN_92          -> imaging   (all_images_measurements.csv)
    20260925_MN_110,111     -> kinetics  (MN_no-glu_<id>_30s_fit-k_md_SN_drop.csv)
Pairs are read from folder_pairs.csv (--pairs), one row per pair:
    batch,imaging_folders,kinetics_folders,include
    20260925_MN,20260925_MN_92,"20260925_MN_110,111",yes
Several folders in one cell are separated by ";". Set include to "no" to
skip a pair. To add a new experiment, add a row.
If the pairs file does not exist, folders sharing the same <date>_<cell>
prefix are paired automatically; --make-pairs writes those to the file so
they can be edited.

Conditions
----------
A condition is treatment_treat_time_cell_type, e.g. "50mM 2DG_1h_MN".
For the kinetics files the condition is parsed from cell_id
("50mM 2DG_1h_yes_A_MN_0_1_0001" -> treatment, treat_time, pyruvate, ...).
Rows without cell_type information are assigned "MN".

Normalisation
-------------
Within each batch and cell type, every cell's ATP/ADP ratio is divided by
the median ratio of the 0mM 2DG control, so the control median is 1.
The kinetics values are not normalised.

y-axis metrics (x-axis is always the normalised ATP/ADP ratio)
--------------
    k_vs_ratio      median adj_k_0.80_cond (--k-col)
    rate_vs_ratio   median (1-b)*k, using only fits with R-sq > 0.8.
                    "1-b" is the column in the kinetics file (= |1-b|).

Import vs export (kinetics / NCT folders only, no imaging)
----------------
x = export reporter 111, y = import reporter 110, one point per condition,
for k and for (1-b)*k. Colour = treatment (light -> dark blue with dose),
marker and line style = batch. Uses every condition measured by both reporters.

Error bars
----------
Every graph is drawn twice:
    _iqr    interquartile range (25th-75th percentile) of the cells
    _ci95   95% bootstrap confidence interval of the median
            (5000 resamples, fixed seed, so reruns give the same result)

Output (--out-dir, default ./output)
------
    <metric>/<metric>_<ratio>_all_<err>.png            all pairs and reporters
    <metric>/<metric>_<ratio>_reporter<id>_<err>.png   one reporter, all pairs
    <metric>/<metric>_<ratio>_<batch>_<err>.png        one pair, both reporters
    import_vs_export/import_vs_export_<k|rate>_<all|batch>_<err>.png
    import_vs_export/import_vs_export_summary.csv
    ratio_vs_k_summary.csv             medians, IQRs and n per condition
    imaging_normalised_<batch>.csv     per-cell ratios incl. normalised columns

Usage
-----
    python plot_ratio_vs_k.py                       # both ratio metrics
    python plot_ratio_vs_k.py --ratio median
    python plot_ratio_vs_k.py --base-dir . --k-col adj_k_0.80_cond --out-dir output
    python plot_ratio_vs_k.py --make-pairs          # (re)write folder_pairs.csv
"""

import argparse
import re
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd

PAIRS_FILE = "folder_pairs.csv"
PAIRS_COLUMNS = ["batch", "imaging_folders", "kinetics_folders", "include"]
IMAGING_FILE = "all_images_measurements.csv"
KINETICS_GLOB = "*_fit-k_md_SN_drop.csv"
DEFAULT_CELL_TYPE = "MN"
CONTROL_TREATMENT = "0mM 2DG"
RATIO_COLS = {
    "median": "bright_median_SN_ratio_488_405",
    "mean": "bright_mean_SN_ratio_488_405",
}
# Some exports name the mean-based ratio without "_mean".
RATIO_FALLBACKS = {"bright_mean_SN_ratio_488_405": "bright_SN_ratio_488_405"}
GROUP_COLS = ["treatment", "treat_time", "cell_type"]
FOLDER_RE = re.compile(r"^(?P<batch>\d{8}_[^_]+)_(?P<ids>[\d,]+)$")
# Categorical slots 1..8 (fixed order, never cycled).
SERIES_COLORS = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100",
                 "#e87ba4", "#008300", "#8a5cd6", "#6b6b6b"]
RSQ_THRESHOLD = 0.8
RATE_COL = f"(1-b)*k_Rsq>{RSQ_THRESHOLD}"
BATCH_MARKERS = ["o", "s", "^", "D", "v", "P"]
BATCH_LINESTYLES = ["-", "--", ":", "-."]
# Sequential blue ramp (steps 250-700) for treatment dose, light -> dark.
DOSE_RAMP = ["#86b6ef", "#6da7ec", "#5598e7", "#3987e5", "#2a78d6", "#256abf",
             "#1c5cab", "#184f95", "#104281", "#0d366b"]
IMPORT_REPORTER = "110"
EXPORT_REPORTER = "111"
N_BOOT = 5000
BOOT_SEED = 0
# Error-bar styles: name -> (low suffix, high suffix, footnote text).
ERROR_BARS = {
    "iqr": ("q25", "q75", "bars: interquartile range (25th-75th percentile)"),
    "ci95": ("ci_lo", "ci_hi", "bars: 95% bootstrap confidence interval of the median"),
}


# ---------------------------------------------------------------- pairing
def find_batches(base_dir):
    """Return {batch: {"imaging": [folders], "kinetics": [folders]}}."""
    batches = {}
    for folder in sorted(p for p in base_dir.iterdir() if p.is_dir()):
        m = FOLDER_RE.match(folder.name)
        if not m:
            continue
        entry = batches.setdefault(m["batch"], {"imaging": [], "kinetics": []})
        if (folder / IMAGING_FILE).exists():
            entry["imaging"].append(folder)
        if any(folder.glob(KINETICS_GLOB)):
            entry["kinetics"].append(folder)
    return {b: e for b, e in batches.items() if e["imaging"] and e["kinetics"]}


def write_pairs_file(batches, path):
    rows = [{"batch": b,
             "imaging_folders": ";".join(f.name for f in e["imaging"]),
             "kinetics_folders": ";".join(f.name for f in e["kinetics"]),
             "include": "yes"} for b, e in batches.items()]
    pd.DataFrame(rows, columns=PAIRS_COLUMNS).to_csv(path, index=False)
    print(f"wrote {len(rows)} pair(s) to {path}")


def read_pairs_file(path, base_dir):
    """Return {batch: {"imaging": [folders], "kinetics": [folders]}} from the pairs file."""
    df = pd.read_csv(path, dtype=str).fillna("")
    missing = set(PAIRS_COLUMNS[:3]) - set(df.columns)
    if missing:
        raise SystemExit(f"{path} is missing column(s) {sorted(missing)}")
    batches = {}
    for _, row in df.iterrows():
        if row.get("include", "yes").strip().lower() in ("no", "n", "false", "0"):
            print(f"skipping pair {row['batch']} (include = {row['include']})")
            continue
        entry = {}
        for key, col in (("imaging", "imaging_folders"), ("kinetics", "kinetics_folders")):
            folders = [base_dir / f.strip() for f in row[col].split(";") if f.strip()]
            for f in folders:
                if not f.is_dir():
                    raise SystemExit(f"{path}: folder not found for {row['batch']}: {f}")
            entry[key] = folders
        if not entry["imaging"] or not entry["kinetics"]:
            raise SystemExit(f"{path}: pair {row['batch']} needs imaging and kinetics folders")
        batch = row["batch"].strip() or entry["imaging"][0].name
        batches[batch] = entry
    return batches


# ---------------------------------------------------------------- loading
def fill_cell_type(df):
    if "cell_type" not in df.columns:
        df["cell_type"] = DEFAULT_CELL_TYPE
    df["cell_type"] = df["cell_type"].fillna(DEFAULT_CELL_TYPE)
    df.loc[df["cell_type"].astype(str).str.strip() == "", "cell_type"] = DEFAULT_CELL_TYPE
    return df


def load_imaging(folder):
    df = pd.read_csv(folder / IMAGING_FILE)
    for col, alt in RATIO_FALLBACKS.items():
        if col not in df.columns and alt in df.columns:
            df[col] = df[alt]
    df["source"] = folder.name
    return fill_cell_type(df)


def parse_cell_id(cell_id):
    """'50mM 2DG_1h_yes_A_MN_0_1_0001' -> treatment, treat_time, pyruvate, cell_type."""
    parts = str(cell_id).split("_")
    return pd.Series({
        "treatment": parts[0] if len(parts) > 0 else None,
        "treat_time": parts[1] if len(parts) > 1 else None,
        "pyruvate": parts[2] if len(parts) > 2 else None,
        "cell_type_from_id": parts[4] if len(parts) > 4 else None,
    })


def load_kinetics(path):
    df = pd.read_csv(path)
    df = pd.concat([df, df["cell_id"].apply(parse_cell_id)], axis=1)
    # Prefer an explicit cell_type column, then "cell", then the cell_id token.
    if "cell_type" not in df.columns:
        df["cell_type"] = df["cell"] if "cell" in df.columns else None
    df["cell_type"] = df["cell_type"].fillna(df["cell_type_from_id"])
    if "reporter" not in df.columns:
        m = re.search(r"_(\d+)_\d+s_", path.name)
        df["reporter"] = m.group(1) if m else path.stem
    df["reporter"] = df["reporter"].astype(str)
    df["source"] = path.name
    # Rate = (1-b)*k, kept only for fits with R-sq above the threshold.
    # The file's "1-b" column holds |1-b|; fall back to computing it.
    one_minus_b = df["1-b"] if "1-b" in df.columns else (1 - df["b"]).abs()
    df[RATE_COL] = (one_minus_b * df["k"]).where(df["R-sq"] > RSQ_THRESHOLD)
    return fill_cell_type(df)


# ---------------------------------------------------------------- normalise
def normalise_to_control(df, ratio_col):
    """Divide ratio_col by the 0mM 2DG median of the same cell type."""
    control = (df[df["treatment"] == CONTROL_TREATMENT]
               .groupby("cell_type")[ratio_col].median())
    missing = set(df["cell_type"]) - set(control.index)
    if missing:
        print(f"  no {CONTROL_TREATMENT} control for cell type(s) {sorted(missing)}; "
              f"their {ratio_col} is left unnormalised (NaN)")
    norm_col = f"{ratio_col}_norm"
    df[norm_col] = df[ratio_col] / df["cell_type"].map(control)
    return norm_col


# ---------------------------------------------------------------- summary
def q25(s):
    return s.quantile(0.25)


def q75(s):
    return s.quantile(0.75)


def bootstrap_median_ci(values, n_boot=N_BOOT, level=0.95, seed=BOOT_SEED):
    """Percentile bootstrap confidence interval of the median."""
    values = np.asarray(values, dtype=float)
    if len(values) < 2:
        return np.nan, np.nan
    rng = np.random.default_rng(seed)
    medians = np.median(rng.choice(values, size=(n_boot, len(values))), axis=1)
    alpha = (1 - level) / 2
    return tuple(np.quantile(medians, [alpha, 1 - alpha]))


def ci_lo(s):
    return bootstrap_median_ci(s)[0]


def ci_hi(s):
    return bootstrap_median_ci(s)[1]


def summarise(df, value_col, prefix):
    stats = ["median", q25, q75, ci_lo, ci_hi, "count"]
    g = df.groupby(GROUP_COLS, dropna=False)[value_col].agg(stats)
    g.columns = [f"{prefix}_{c}" for c in ["median", "q25", "q75", "ci_lo", "ci_hi", "n"]]
    return g


def dose_key(treatment):
    m = re.match(r"([\d.]+)", str(treatment))
    return float(m.group(1)) if m else float("inf")


def build_summary(imaging, kinetics, ratio_col, norm_col, y_cols):
    """Join per-condition imaging medians with kinetics medians.

    y_cols maps a prefix used in the output columns to a kinetics column.
    """
    img = summarise(imaging.dropna(subset=[norm_col]), norm_col, "ratio_norm")
    raw = summarise(imaging.dropna(subset=[ratio_col]), ratio_col, "ratio_raw")
    img = img.join(raw[["ratio_raw_median"]])
    rows = []
    for reporter, kdf in kinetics.groupby("reporter"):
        kin = pd.concat([summarise(kdf.dropna(subset=[col]), col, prefix)
                         for prefix, col in y_cols.items()], axis=1)
        merged = img.join(kin, how="inner").reset_index()
        merged.insert(0, "reporter", reporter)
        rows.append(merged)
        unmatched = sorted(set(kin.index) - set(img.index))
        if unmatched:
            print(f"  reporter {reporter}: no imaging match for "
                  + ", ".join("_".join(map(str, u)) for u in unmatched))
    out = pd.concat(rows, ignore_index=True)
    out.insert(1, "condition", out[GROUP_COLS].astype(str).agg("_".join, axis=1))
    out["ratio_col"] = ratio_col
    for prefix, col in y_cols.items():
        out[f"{prefix}_col"] = col
    out["_dose"] = out["treatment"].map(dose_key)
    return out.sort_values(["reporter", "cell_type", "_dose"]).drop(columns="_dose")


# ---------------------------------------------------------------- plotting
def series_styles(summary):
    """Fixed colour per (batch, reporter) and marker per batch, so a series
    looks the same in every graph."""
    batches = list(dict.fromkeys(summary["batch"]))
    keys = list(dict.fromkeys(zip(summary["batch"], summary["reporter"])))
    return {key: (SERIES_COLORS[i % len(SERIES_COLORS)],
                  BATCH_MARKERS[batches.index(key[0]) % len(BATCH_MARKERS)])
            for i, key in enumerate(keys)}


def plot_summary(summary, y_prefix, ratio_col, y_label, title, out_png, styles, err):
    """One point per (batch, reporter, condition).

    x = normalised ATP/ADP ratio, y = the kinetics metric named by y_prefix.
    """
    fig, ax = new_axes()
    ax.axvline(1, color="#999999", linewidth=1, linestyle="--", zorder=1)

    for (batch, reporter), sub in summary.groupby(["batch", "reporter"], sort=False):
        sub = sub.dropna(subset=[f"{y_prefix}_median"])
        sub = sub.sort_values("treatment", key=lambda s: s.map(dose_key))
        color, marker = styles[(batch, reporter)]
        x, y = sub["ratio_norm_median"], sub[f"{y_prefix}_median"]
        lo, hi, _ = ERROR_BARS[err]
        xerr = [x - sub[f"ratio_norm_{lo}"], sub[f"ratio_norm_{hi}"] - x]
        yerr = [y - sub[f"{y_prefix}_{lo}"], sub[f"{y_prefix}_{hi}"] - y]
        ax.errorbar(x, y, xerr=xerr, yerr=yerr,
                    fmt="none", ecolor=color, elinewidth=1, alpha=0.35, capsize=0)
        ax.plot(x, y, color=color, linewidth=1, alpha=0.5, zorder=2)
        ax.scatter(x, y, s=64, color=color,
                   marker=marker, edgecolor="white", linewidth=1.5, zorder=3,
                   label=f"{batch} · reporter {reporter}")
        for xi, yi, t in zip(x, y, sub["treatment"]):
            ax.annotate(str(t).replace(" 2DG", ""), (xi, yi),
                        xytext=(5, 4), textcoords="offset points",
                        fontsize=7, color="#555555")

    ax.set_xlabel(f"median {ratio_col}\nnormalised to {CONTROL_TREATMENT} median")
    ax.set_ylabel(y_label)
    ax.set_title(title, fontsize=11, loc="left")
    ax.legend(frameon=False, fontsize=8)
    finish(fig, "points: median per condition (labels = 2DG dose); "
           f"{ERROR_BARS[err][2]}; dashed line: control ratio = 1", out_png)


def new_axes():
    fig, ax = plt.subplots(figsize=(8, 6))
    for spine in ("top", "right"):
        ax.spines[spine].set_visible(False)
    ax.grid(True, color="#e6e6e6", linewidth=0.8)
    ax.set_axisbelow(True)
    return fig, ax


def finish(fig, footnote, out_png):
    fig.text(0.01, 0.01, footnote, fontsize=7, color="#777777")
    fig.tight_layout(rect=(0, 0.02, 1, 1))
    fig.savefig(out_png, dpi=200)
    plt.close(fig)


def plot_all(combined, ratio_name, y_metrics, out_dir):
    """Write graphs for every pair together, each reporter, and each pair."""
    styles = series_styles(combined)
    ratio_col = RATIO_COLS[ratio_name]
    batches = list(dict.fromkeys(combined["batch"]))
    scopes = [("all", f"All pairs: {', '.join(batches)}", combined)]
    for reporter, sub in combined.groupby("reporter"):
        scopes.append((f"reporter{reporter}",
                       f"Reporter {reporter}: {', '.join(batches)}", sub))
    for batch, sub in combined.groupby("batch", sort=False):
        scopes.append((batch, f"{batch}: {sub['pair'].iloc[0]}", sub))

    for y_prefix, (folder, y_label) in y_metrics.items():
        sub_dir = out_dir / folder
        sub_dir.mkdir(parents=True, exist_ok=True)
        for scope, title, data in scopes:
            for err in ERROR_BARS:
                out_png = sub_dir / f"{folder}_{ratio_name}_{scope}_{err}.png"
                plot_summary(data, y_prefix, ratio_col, y_label, title, out_png,
                             styles, err)
                print(f"saved {out_png}")


# ------------------------------------------------- import vs export (NCT only)
def summarise_kinetics(kinetics, y_cols):
    """Per-condition summary of each kinetics metric, one row per reporter."""
    rows = []
    for reporter, kdf in kinetics.groupby("reporter"):
        kin = pd.concat([summarise(kdf.dropna(subset=[col]), col, prefix)
                         for prefix, col in y_cols.items()], axis=1).reset_index()
        kin.insert(0, "reporter", reporter)
        rows.append(kin)
    return pd.concat(rows, ignore_index=True)


def import_export_table(kin_summary, y_cols):
    """Wide table: one row per (batch, condition), export (111) and import (110)
    columns side by side. Only conditions measured by both reporters are kept."""
    keys = ["batch"] + GROUP_COLS
    stat_cols = [c for c in kin_summary.columns
                 if any(c.startswith(f"{p}_") for p in y_cols)]
    exp = kin_summary[kin_summary["reporter"] == EXPORT_REPORTER][keys + stat_cols]
    imp = kin_summary[kin_summary["reporter"] == IMPORT_REPORTER][keys + stat_cols]
    wide = exp.merge(imp, on=keys, suffixes=("_export", "_import"))
    wide.insert(1, "condition", wide[GROUP_COLS].astype(str).agg("_".join, axis=1))
    wide["_dose"] = wide["treatment"].map(dose_key)
    return wide.sort_values(["batch", "cell_type", "_dose"]).drop(columns="_dose")


def treatment_colors(treatments):
    """Sequential blue ramp (light -> dark) by dose, so colour follows the
    treatment and stays the same in every import vs export graph."""
    ordered = sorted(set(treatments), key=dose_key)
    n = len(ordered)
    idx = [round(i * (len(DOSE_RAMP) - 1) / max(n - 1, 1)) for i in range(n)]
    return {t: DOSE_RAMP[j] for t, j in zip(ordered, idx)}


def plot_import_export(wide, y_prefix, metric_label, title, out_png, styles, err):
    """x = export reporter (111), y = import reporter (110); one point per condition.

    Colour = treatment, marker / line style = batch.
    """
    colors, batch_styles = styles
    fig, ax = new_axes()
    lo, hi, err_text = ERROR_BARS[err]
    xc, yc = f"{y_prefix}_median_export", f"{y_prefix}_median_import"
    for batch, sub in wide.groupby("batch", sort=False):
        sub = sub.dropna(subset=[xc, yc])
        sub = sub.sort_values("treatment", key=lambda s: s.map(dose_key))
        marker, linestyle = batch_styles[batch]
        ax.plot(sub[xc], sub[yc], color="#9a9a9a", linewidth=1, linestyle=linestyle,
                alpha=0.7, zorder=2)
        for _, r in sub.iterrows():
            color = colors[r["treatment"]]
            x, y = r[xc], r[yc]
            ax.errorbar(x, y,
                        xerr=[[x - r[f"{y_prefix}_{lo}_export"]],
                              [r[f"{y_prefix}_{hi}_export"] - x]],
                        yerr=[[y - r[f"{y_prefix}_{lo}_import"]],
                              [r[f"{y_prefix}_{hi}_import"] - y]],
                        fmt="none", ecolor=color, elinewidth=1, alpha=0.5, capsize=0)
            ax.scatter(x, y, s=70, color=color, marker=marker, edgecolor="white",
                       linewidth=1.5, zorder=3)
            ax.annotate(str(r["treatment"]).replace(" 2DG", ""), (x, y),
                        xytext=(5, 4), textcoords="offset points",
                        fontsize=7, color="#555555")

    present = sorted(set(wide["treatment"]), key=dose_key)
    handles = [Line2D([], [], linestyle="", marker="o", markersize=8,
                      markerfacecolor=colors[t], markeredgecolor="white", label=t)
               for t in present]
    handles += [Line2D([], [], color="#9a9a9a", linestyle=ls, marker=m, markersize=7,
                       markerfacecolor="#9a9a9a", markeredgecolor="white", label=b)
                for b, (m, ls) in batch_styles.items() if b in set(wide["batch"])]
    ax.legend(handles=handles, frameon=False, fontsize=8,
              loc="upper left", bbox_to_anchor=(1.01, 1), borderaxespad=0)
    ax.set_xlabel(f"export: reporter {EXPORT_REPORTER}, median {metric_label}")
    ax.set_ylabel(f"import: reporter {IMPORT_REPORTER}, median {metric_label}")
    ax.set_title(title, fontsize=11, loc="left")
    finish(fig, "points: median per condition, colour = treatment, "
           f"marker = batch; {err_text}", out_png)


def plot_all_import_export(kin_summary, y_cols, out_dir):
    wide = import_export_table(kin_summary, y_cols)
    if wide.empty:
        print(f"no conditions with both reporters {IMPORT_REPORTER} and "
              f"{EXPORT_REPORTER}; skipping import vs export graphs")
        return
    batches = list(dict.fromkeys(wide["batch"]))
    batch_styles = {b: (BATCH_MARKERS[i % len(BATCH_MARKERS)],
                        BATCH_LINESTYLES[i % len(BATCH_LINESTYLES)])
                    for i, b in enumerate(batches)}
    styles = (treatment_colors(wide["treatment"]), batch_styles)
    labels = {"k": y_cols["k"], "rate": f"(1-b)*k (R-sq > {RSQ_THRESHOLD})"}
    folder = "import_vs_export"
    sub_dir = out_dir / folder
    sub_dir.mkdir(parents=True, exist_ok=True)
    scopes = [("all", f"All batches: {', '.join(batches)}", wide)]
    scopes += [(b, b, sub) for b, sub in wide.groupby("batch", sort=False)]
    for y_prefix in y_cols:
        for scope, title, data in scopes:
            for err in ERROR_BARS:
                out_png = sub_dir / f"{folder}_{y_prefix}_{scope}_{err}.png"
                plot_import_export(data, y_prefix, labels[y_prefix], title,
                                   out_png, styles, err)
                print(f"saved {out_png}")
    out_csv = sub_dir / "import_vs_export_summary.csv"
    wide.to_csv(out_csv, index=False)
    print(f"saved {out_csv}")


# ---------------------------------------------------------------- main
def process_batch(batch, folders, ratios, y_cols, out_dir):
    """Return ({ratio name: summary DataFrame}, kinetics summary) for one pair."""
    pair = (f"{', '.join(f.name for f in folders['imaging'])} + "
            f"{', '.join(f.name for f in folders['kinetics'])}")
    print(f"Batch {batch}: {pair}")
    imaging = pd.concat([load_imaging(f) for f in folders["imaging"]], ignore_index=True)
    kin_files = sorted(p for f in folders["kinetics"] for p in f.glob(KINETICS_GLOB))
    for p in kin_files:
        print(f"  kinetics file: {p.parent.name}/{p.name}")
    kinetics = pd.concat([load_kinetics(p) for p in kin_files], ignore_index=True)
    for col in y_cols.values():
        if col not in kinetics.columns:
            raise SystemExit(f"Column {col!r} not found in kinetics files")

    results = {}
    for name in ratios:
        ratio_col = RATIO_COLS[name]
        if ratio_col not in imaging.columns:
            print(f"  skipping {ratio_col}: not in imaging data")
            continue
        norm_col = normalise_to_control(imaging, ratio_col)
        summary = build_summary(imaging, kinetics, ratio_col, norm_col, y_cols)
        summary.insert(0, "batch", batch)
        summary.insert(1, "pair", pair)
        results[name] = summary

    norm_cols = [c for c in imaging.columns if c.endswith("_norm")]
    keep = ["source", "image_name", "cell_id", "treatment", "treat_time", "pyruvate",
            "cell_type"] + [RATIO_COLS[n] for n in results] + norm_cols
    out_csv = out_dir / f"imaging_normalised_{batch}.csv"
    imaging[[c for c in keep if c in imaging.columns]].to_csv(out_csv, index=False)
    print(f"  saved {out_csv}")
    kin_summary = summarise_kinetics(kinetics, y_cols)
    kin_summary.insert(0, "batch", batch)
    return results, kin_summary


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--base-dir", type=Path, default=Path(__file__).resolve().parent,
                    help="directory containing the experiment folders")
    ap.add_argument("--ratio", choices=["median", "mean", "both"], default="both",
                    help="bright_median_SN_ratio_488_405 and/or bright_mean_SN_ratio_488_405")
    ap.add_argument("--k-col", default="adj_k_0.80_cond",
                    help="rate-constant column in the kinetics files")
    ap.add_argument("--pairs", type=Path, default=None,
                    help=f"folder pairs file (default: <base-dir>/{PAIRS_FILE})")
    ap.add_argument("--make-pairs", action="store_true",
                    help="auto-detect folder pairs, write them to the pairs file and exit")
    ap.add_argument("--out-dir", type=Path, default=None,
                    help="output directory (default: <base-dir>/output)")
    args = ap.parse_args()

    ratios = ["median", "mean"] if args.ratio == "both" else [args.ratio]
    out_dir = args.out_dir or args.base_dir / "output"
    out_dir.mkdir(parents=True, exist_ok=True)
    pairs_path = args.pairs or args.base_dir / PAIRS_FILE
    if args.make_pairs:
        write_pairs_file(find_batches(args.base_dir), pairs_path)
        return
    if pairs_path.exists():
        print(f"reading folder pairs from {pairs_path}")
        batches = read_pairs_file(pairs_path, args.base_dir)
    else:
        print(f"{pairs_path} not found; pairing folders by <date>_<cell> prefix")
        batches = find_batches(args.base_dir)
    if not batches:
        raise SystemExit(f"No imaging/kinetics folder pairs found in {args.base_dir}")

    # y-axis metrics: summary prefix -> kinetics column / (output folder, axis label)
    y_cols = {"k": args.k_col, "rate": RATE_COL}
    y_metrics = {
        "k": ("k_vs_ratio", f"median {args.k_col}"),
        "rate": ("rate_vs_ratio",
                 f"median (1-b)*k\n(fits with R-sq > {RSQ_THRESHOLD} only)"),
    }

    per_ratio = {name: [] for name in ratios}
    kin_summaries = []
    for batch, folders in batches.items():
        results, kin_summary = process_batch(batch, folders, ratios, y_cols, out_dir)
        kin_summaries.append(kin_summary)
        for name, summary in results.items():
            per_ratio[name].append(summary)

    plot_all_import_export(pd.concat(kin_summaries, ignore_index=True), y_cols, out_dir)

    all_summaries = []
    for name, summaries in per_ratio.items():
        if not summaries:
            continue
        combined = pd.concat(summaries, ignore_index=True)
        all_summaries.append(combined)
        plot_all(combined, name, y_metrics, out_dir)

    if all_summaries:
        out_csv = out_dir / "ratio_vs_k_summary.csv"
        pd.concat(all_summaries, ignore_index=True).to_csv(out_csv, index=False)
        print(f"saved {out_csv}")


if __name__ == "__main__":
    main()
