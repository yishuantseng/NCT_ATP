"""Pair imaging and kinetics experiments from the same cell batch and plot
the median 488/405 ATP ratio against the median fitted rate constant.

Folder pairing
--------------
Each folder is one experiment, named <date>_<cell>_<reporter ids>, e.g.
    20260925_MN_92          -> imaging   (all_images_measurements.csv)
    20260925_MN_110,111     -> kinetics  (MN_no-glu_<id>_30s_fit-k_md_SN_drop.csv)
Folders sharing the same <date>_<cell> prefix are the same batch of cells
and are paired.

Conditions
----------
A condition is treatment_treat_time_cell_type, e.g. "50mM 2DG_1h_MN".
For the kinetics files the condition is parsed from cell_id
("50mM 2DG_1h_yes_A_MN_0_1_0001" -> treatment, treat_time, pyruvate, ...).
Rows without cell_type information are assigned "MN".

Output (written into each imaging folder by default)
------
    ratio_vs_k_<ratio>_<batch>.png   one point per condition per kinetics reporter
    ratio_vs_k_summary_<batch>.csv   medians, IQRs and n per condition

Usage
-----
    python plot_ratio_vs_k.py                       # both ratio metrics
    python plot_ratio_vs_k.py --ratio median
    python plot_ratio_vs_k.py --base-dir . --k-col adj_k_0.80_cond
"""

import argparse
import re
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd

IMAGING_FILE = "all_images_measurements.csv"
KINETICS_GLOB = "*_fit-k_md_SN_drop.csv"
DEFAULT_CELL_TYPE = "MN"
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
    return fill_cell_type(df)


# ---------------------------------------------------------------- summary
def q25(s):
    return s.quantile(0.25)


def q75(s):
    return s.quantile(0.75)


def summarise(df, value_col, prefix):
    g = df.groupby(GROUP_COLS, dropna=False)[value_col].agg(["median", q25, q75, "count"])
    g.columns = [f"{prefix}_{c}" for c in ["median", "q25", "q75", "n"]]
    return g


def condition_label(row):
    return "_".join(str(row[c]) for c in GROUP_COLS)


def dose_key(treatment):
    m = re.match(r"([\d.]+)", str(treatment))
    return float(m.group(1)) if m else float("inf")


def build_summary(imaging, kinetics, ratio_col, k_col):
    img = summarise(imaging.dropna(subset=[ratio_col]), ratio_col, "ratio")
    rows = []
    for reporter, kdf in kinetics.groupby("reporter"):
        kin = summarise(kdf.dropna(subset=[k_col]), k_col, "k")
        merged = img.join(kin, how="inner").reset_index()
        merged.insert(0, "reporter", reporter)
        rows.append(merged)
        unmatched = sorted(set(kin.index) - set(img.index))
        if unmatched:
            print(f"  reporter {reporter}: no imaging match for "
                  + ", ".join("_".join(map(str, u)) for u in unmatched))
    out = pd.concat(rows, ignore_index=True)
    out.insert(1, "condition", out.apply(condition_label, axis=1))
    out["ratio_col"] = ratio_col
    out["k_col"] = k_col
    return out.sort_values(["reporter", "cell_type", "treatment"],
                           key=lambda s: s.map(dose_key) if s.name == "treatment" else s)


# ---------------------------------------------------------------- plotting
def plot_summary(summary, ratio_col, k_col, title, out_png):
    fig, ax = plt.subplots(figsize=(7, 5.5))
    for spine in ("top", "right"):
        ax.spines[spine].set_visible(False)
    ax.grid(True, color="#e6e6e6", linewidth=0.8)
    ax.set_axisbelow(True)

    for i, (reporter, sub) in enumerate(summary.groupby("reporter")):
        color = SERIES_COLORS[i % len(SERIES_COLORS)]
        xerr = [sub["k_median"] - sub["k_q25"], sub["k_q75"] - sub["k_median"]]
        yerr = [sub["ratio_median"] - sub["ratio_q25"], sub["ratio_q75"] - sub["ratio_median"]]
        ax.errorbar(sub["k_median"], sub["ratio_median"], xerr=xerr, yerr=yerr,
                    fmt="none", ecolor=color, elinewidth=1, alpha=0.45, capsize=0)
        ax.scatter(sub["k_median"], sub["ratio_median"], s=64, color=color,
                   edgecolor="white", linewidth=1.5, zorder=3,
                   label=f"reporter {reporter}")
        for _, r in sub.iterrows():
            ax.annotate(r["condition"], (r["k_median"], r["ratio_median"]),
                        xytext=(6, 4), textcoords="offset points",
                        fontsize=8, color="#444444")

    ax.set_xlabel(f"median {k_col}")
    ax.set_ylabel(f"median {ratio_col}")
    ax.set_title(title, fontsize=11, loc="left")
    ax.legend(frameon=False, fontsize=9)
    fig.text(0.01, 0.01, "points: median per condition; bars: interquartile range",
             fontsize=7, color="#777777")
    fig.tight_layout()
    fig.savefig(out_png, dpi=200)
    plt.close(fig)


# ---------------------------------------------------------------- main
def run_batch(batch, folders, ratios, k_col, out_dir):
    print(f"Batch {batch}: imaging={[f.name for f in folders['imaging']]} "
          f"kinetics={[f.name for f in folders['kinetics']]}")
    imaging = pd.concat([load_imaging(f) for f in folders["imaging"]], ignore_index=True)
    kin_files = sorted(p for f in folders["kinetics"] for p in f.glob(KINETICS_GLOB))
    for p in kin_files:
        print(f"  kinetics file: {p.parent.name}/{p.name}")
    kinetics = pd.concat([load_kinetics(p) for p in kin_files], ignore_index=True)
    if k_col not in kinetics.columns:
        raise SystemExit(f"Column {k_col!r} not found in kinetics files")

    out_dir = out_dir or folders["imaging"][0]
    out_dir.mkdir(parents=True, exist_ok=True)
    summaries = []
    for name in ratios:
        ratio_col = RATIO_COLS[name]
        if ratio_col not in imaging.columns:
            print(f"  skipping {ratio_col}: not in imaging data")
            continue
        summary = build_summary(imaging, kinetics, ratio_col, k_col)
        summaries.append(summary)
        title = (f"{batch}: {', '.join(f.name for f in folders['imaging'])} vs "
                 f"{', '.join(f.name for f in folders['kinetics'])}")
        out_png = out_dir / f"ratio_vs_k_{name}_{batch}.png"
        plot_summary(summary, ratio_col, k_col, title, out_png)
        print(f"  saved {out_png}")

    if summaries:
        out_csv = out_dir / f"ratio_vs_k_summary_{batch}.csv"
        pd.concat(summaries, ignore_index=True).to_csv(out_csv, index=False)
        print(f"  saved {out_csv}")


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--base-dir", type=Path, default=Path(__file__).resolve().parent,
                    help="directory containing the experiment folders")
    ap.add_argument("--ratio", choices=["median", "mean", "both"], default="both",
                    help="bright_median_SN_ratio_488_405 and/or bright_mean_SN_ratio_488_405")
    ap.add_argument("--k-col", default="adj_k_0.80_cond",
                    help="rate-constant column in the kinetics files")
    ap.add_argument("--out-dir", type=Path, default=None,
                    help="output directory (default: the imaging folder)")
    args = ap.parse_args()

    ratios = ["median", "mean"] if args.ratio == "both" else [args.ratio]
    batches = find_batches(args.base_dir)
    if not batches:
        raise SystemExit(f"No imaging/kinetics folder pairs found in {args.base_dir}")
    for batch, folders in batches.items():
        run_batch(batch, folders, ratios, args.k_col, args.out_dir)


if __name__ == "__main__":
    main()
