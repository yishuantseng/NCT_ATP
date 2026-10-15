"""Nucleus clustering vs nucleocytoplasmic transport rate.

For every nucleus in the Cellpose nuclear masks (one label per nucleus, label = Fiji ROI id) this
measures how clustered it is with its neighbours, then joins that to the per-cell transport rates
from nct_analysis.py (rates_per_cell.csv) and asks whether clustered cells transport faster or slower,
within each reporter x treatment.

Clustering metrics (all from one mask frame, pre-activation by default; pixels):
  n_neighbors    other nuclei whose centroid lies within RADIUS_PX of this one's centroid
  mean_dist_px   mean centroid distance to the K_NEAREST nearest nuclei ("average distance to
                 surrounding nuclei")
  cluster_position  isolated (cluster_size 1), middle (neighbours within RADIUS_PX on all sides: the
                 largest empty angle around the nucleus is < EDGE_ANGLE_DEG) or edge (otherwise)
  nn_dist_px     centroid distance to the nearest other nucleus
  nn_gap_px      edge-to-edge gap to the nearest other nucleus (0-1 = touching)
  cluster_size   size of the group of nuclei linked by gaps <= CLUSTER_GAP_PX (1 = isolated)
  clustered      cluster_size >= 2
  n_nuclei_field nuclei in the whole field (a density confound to keep in mind)
  near_edge      centroid within RADIUS_PX of the image border (neighbour count is an underestimate)

Each EXPERIMENTS entry is one experiment: its mask folder, its rates CSV and an output folder,
plus a pooled folder (POOLED_OUT) if there is more than one experiment.
"""
import os
import re

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker
import numpy as np
import pandas as pd
import tifffile
from scipy import ndimage, stats
from scipy.sparse.csgraph import connected_components

# ---------------------------------------------------------------- parameters
EXPERIMENTS = [
    # name: experiment label (used when the rates file has no source column); masks: folder of nuclear label stacks; rates: per-cell rates CSV, either rates_per_cell.csv from
    # nct_analysis.py or rates_nuc_vs_nc_per_cell.csv from nc_rates.py; out: output folder
    dict(name="20261010_MN_110,111",
         masks="/Volumes/holtl02lab/holtl02labspace/Holt_Lab_Members/Yi-Shuan/microscope/2026/10-12/20261010_MN_110,111/combine/DNA/erode",
         mask_suffix="_640_mask.tif",
         rates="/Volumes/holtl02lab/holtl02labspace/Holt_Lab_Members/Yi-Shuan/microscope/2026/10-12/20261010_MN_110,111/analysis/nct_results/rates_per_cell.csv",
         out="/Volumes/holtl02lab/holtl02labspace/Holt_Lab_Members/Yi-Shuan/microscope/2026/10-12/20261010_MN_110,111/analysis/nct_results/nucleus_clustering"),
    dict(name="20260925_MN_110,111",
         masks="/Volumes/holtl02lab/holtl02labspace/Holt_Lab_Members/Yi-Shuan/microscope/2026/07-09/20260925_MN_110,111/combine/DNA",
         mask_suffix="_640_cp_masks.tif",
         rates="/Volumes/holtl02lab/holtl02labspace/Holt_Lab_Members/Yi-Shuan/microscope/2026/07-09/20260925_MN_110,111/combine/nc_ratio/rates_nuc_vs_nc_per_cell.csv",
         out="/Volumes/holtl02lab/holtl02labspace/Holt_Lab_Members/Yi-Shuan/microscope/2026/07-09/20260925_MN_110,111/combine/nc_ratio/nucleus_clustering"),
]
POOLED_OUT = None            # folder for pooled outputs when there are several experiments (None = skip)
MASK_FRAME = 0               # 0-based frame of the mask stack to measure (0 = first, pre-activation)
MIN_NUC_AREA = 300           # px; smaller labels are fragments and are not counted as neighbours
RADIUS_PX = 100              # neighbourhood radius for n_neighbors (~2 nuclear diameters)
CLUSTER_GAP_PX = 10          # nuclei closer than this (edge to edge) belong to the same cluster
K_NEAREST = 3                # mean_dist_px = mean distance to this many nearest nuclei
EDGE_ANGLE_DEG = 180         # clustered nucleus is "edge" if its neighbours leave an empty sector >= this
EXCLUDE_REPORTERS = [185]
RATE_METRICS = ["flux_pct_per_min", "slope_pct_per_min", "k_per_min"]   # k uses well-fit cells only
CLUSTER_METRICS = ["mean_dist_px", "n_neighbors", "nn_gap_px", "nn_dist_px", "cluster_size", "n_nuclei_field"]
SCATTER_METRICS = ["mean_dist_px", "n_neighbors"]
POSITIONS = ["isolated", "edge", "middle"]
NC_RATES_READOUT = "nuclear"   # nc_rates.py files only: "nuclear" or "N/C"
NC_RATES_TREATMENT = "{dose}mM 2DG"   # nc_rates.py files only: treatment label from dose_mM
MIN_R2 = 0.5                 # nc_rates.py files only: well-fit = r2 >= MIN_R2 and k below its cap
K_CAP_PER_MIN = 2.0          # nc_rates.py files only: k cap (1/30 s)
POOL_TREATMENTS = True       # also analyse all treatments together (per reporter) as an extra "all" group
ALL_LABEL = "all"
DEFAULT_REPLICATE = "A"
RANDOM_SEED = 0
# ----------------------------------------------------------------------------

KEY_RE = re.compile(r"com_(\d+)_.*?_(\d{4})_XY(\d+)")


def field_key(name):
    m = KEY_RE.search(name)
    return (int(m.group(1)), m.group(2), int(m.group(3))) if m else None


def measure_field(mask):
    """Clustering metrics for every nucleus in one label image."""
    ids, areas = np.unique(mask[mask > 0], return_counts=True)
    ids = ids[areas >= MIN_NUC_AREA]
    areas = areas[areas >= MIN_NUC_AREA]
    if len(ids) == 0:
        return pd.DataFrame()
    cents = np.array(ndimage.center_of_mass(mask > 0, mask, ids))
    d = np.hypot(*(cents[:, None, :] - cents[None, :, :]).transpose(2, 0, 1))
    np.fill_diagonal(d, np.inf)
    # edge-to-edge gaps: distance map from each nucleus, minimum over every other nucleus
    gap = np.full((len(ids), len(ids)), np.inf)
    for i, lab in enumerate(ids):
        dt = ndimage.distance_transform_edt(mask != lab)
        gap[i] = ndimage.minimum(dt, labels=mask, index=ids)
    np.fill_diagonal(gap, np.inf)
    # clusters = connected components of the "gap <= CLUSTER_GAP_PX" graph
    _, comp = connected_components(gap <= CLUSTER_GAP_PX, directed=False)
    size = np.bincount(comp)[comp]
    # edge vs middle of a cluster: largest empty angle between neighbours within RADIUS_PX
    position = []
    for i in range(len(ids)):
        nb = np.where(d[i] <= RADIUS_PX)[0]
        if size[i] < 2:
            position.append("isolated")
            continue
        ang = np.sort(np.degrees(np.arctan2(cents[nb, 0] - cents[i, 0], cents[nb, 1] - cents[i, 1])))
        gaps = np.diff(np.r_[ang, ang[0] + 360])
        position.append("middle" if len(nb) >= 3 and gaps.max() < EDGE_ANGLE_DEG else "edge")
    k = min(K_NEAREST, len(ids) - 1)
    mean_dist = np.sort(d, 1)[:, :k].mean(1) if k else np.full(len(ids), np.nan)
    h, w = mask.shape
    edge = np.minimum.reduce([cents[:, 0], cents[:, 1], h - 1 - cents[:, 0], w - 1 - cents[:, 1]])
    one = len(ids) > 1
    return pd.DataFrame(dict(
        roi=ids.astype(int), nuc_area_px=areas, cy=cents[:, 0], cx=cents[:, 1],
        n_neighbors=(d <= RADIUS_PX).sum(1), mean_dist_px=mean_dist,
        nn_dist_px=d.min(1) if one else np.nan, nn_gap_px=gap.min(1) if one else np.nan,
        cluster_size=size, clustered=size >= 2, cluster_position=position, n_nuclei_field=len(ids), near_edge=edge < RADIUS_PX))


def measure_folder(folder, suffix):
    rows = []
    for f in sorted(os.listdir(folder)):
        if not f.endswith(suffix) or field_key(f) is None:
            continue
        mask = tifffile.imread(os.path.join(folder, f), key=MASK_FRAME)
        t = measure_field(mask)
        rep, fid, pos = field_key(f)
        rows.append(t.assign(mask_file=f, reporter=rep, file_id=fid, position=pos))
        print(f"  {f}: {len(t)} nuclei")
    return pd.concat(rows, ignore_index=True)


def load_rates(path, name):
    r = pd.read_csv(path, dtype={"roi": str, "file_id": str})
    if "readout" in r:   # nc_rates.py output: one row per cell and readout, bright cells only
        r = r[(r.readout == NC_RATES_READOUT) & ~r.wrong_direction].copy()
        r["image"] = r.cell_key.str.rsplit("_", n=1).str[0]
        r["roi"] = r.cell_key.str.rsplit("_", n=1).str[1]
        r["treatment"] = [NC_RATES_TREATMENT.format(dose=f"{d:g}") for d in r.dose_mM]
        r["bright"] = True
        r["good_fit"] = (r.r2 >= MIN_R2) & (r.k_per_min < 0.999 * K_CAP_PER_MIN)
        r = r.drop(columns=["file_id", "position"])
    r = r[~r.reporter.isin(EXCLUDE_REPORTERS)].copy()
    r["roi"] = r.roi.astype(int)
    keys = r.image.map(field_key)
    r["file_id"] = [k[1] for k in keys]
    r["position"] = [k[2] for k in keys]
    if "source" not in r:
        r["source"] = name
    if "replicate" not in r:
        r["replicate"] = DEFAULT_REPLICATE
    if "unit" not in r:
        r["unit"] = r.source.astype(str) + "_" + r.replicate.astype(str)
    return r


def correlations(t):
    out = []
    for (rep, direc, trt), g in t.groupby(["reporter", "direction", "treatment"]):
        for rm in RATE_METRICS:
            cells = g[g.good_fit] if rm == "k_per_min" else g[g.bright]
            for cm in CLUSTER_METRICS:
                x = cells[[cm, rm]].dropna()
                if len(x) >= 5:
                    rho, p = stats.spearmanr(x[cm], x[rm])
                    out.append(dict(reporter=rep, direction=direc, treatment=trt, rate=rm, cluster_metric=cm,
                                    test="Spearman", n=len(x), stat=rho, p=p))
                if trt == ALL_LABEL:   # also with each cell's rate centred on its own treatment median
                    x = cells[[cm, rm, "orig_treatment"]].dropna()
                    y = x[rm] - x.groupby("orig_treatment")[rm].transform("median")
                    if len(x) >= 5:
                        rho, p = stats.spearmanr(x[cm], y)
                        out.append(dict(reporter=rep, direction=direc, treatment=trt, rate=rm, cluster_metric=cm,
                                        test="Spearman, rate centred per treatment", n=len(x), stat=rho, p=p))
            grp = {pos: cells.loc[cells.cluster_position == pos, rm].dropna() for pos in POSITIONS}
            grp = {k: v for k, v in grp.items() if len(v) >= 3}
            med = {f"median_{k}": v.median() for k, v in grp.items()}
            nn = {f"n_{k}": len(v) for k, v in grp.items()}
            if len(grp) >= 2:
                out.append(dict(reporter=rep, direction=direc, treatment=trt, rate=rm,
                                cluster_metric="cluster_position", test="Kruskal-Wallis (positions)",
                                n=sum(map(len, grp.values())), stat=stats.kruskal(*grp.values()).statistic,
                                p=stats.kruskal(*grp.values()).pvalue, **nn, **med))
            if "edge" in grp and "middle" in grp:
                out.append(dict(reporter=rep, direction=direc, treatment=trt, rate=rm,
                                cluster_metric="cluster_position", test="Mann-Whitney middle vs edge",
                                n=len(grp["edge"]) + len(grp["middle"]),
                                stat=np.nan, p=stats.mannwhitneyu(grp["middle"], grp["edge"]).pvalue, **nn, **med))
    return pd.DataFrame(out)


def treatment_order(t):
    trts = list(t.sort_values("dose_mM").treatment.unique()) if "dose_mM" in t else sorted(t.treatment.unique())
    return [x for x in trts if x != ALL_LABEL] + ([ALL_LABEL] if ALL_LABEL in trts else [])


def plot_scatter(t, corr, rm, cm, path):
    reps = sorted(t.reporter.unique())
    trts = treatment_order(t)
    fig, axes = plt.subplots(len(reps), len(trts), figsize=(3.2 * len(trts), 3 * len(reps)), squeeze=False)
    rng = np.random.default_rng(RANDOM_SEED)
    for i, rep in enumerate(reps):
        for j, trt in enumerate(trts):
            ax = axes[i, j]
            g = t[(t.reporter == rep) & (t.treatment == trt)]
            g = g[g.good_fit] if rm == "k_per_min" else g[g.bright]
            discrete = cm in ("n_neighbors", "cluster_size")
            x = g[cm] + (rng.uniform(-0.2, 0.2, len(g)) if discrete else 0)
            if trt == ALL_LABEL:   # pooled panel: colour cells by their own treatment
                for ti, ot in enumerate(trts[:-1]):
                    m = (g.orig_treatment == ot).values
                    ax.scatter(x[m], g[rm][m], s=10, alpha=0.5, color=f"C{ti}", lw=0, label=ot)
                ax.legend(fontsize=6, frameon=False)
            else:
                ax.scatter(x, g[rm], s=10, alpha=0.4, color="0.3", lw=0)
            if discrete:
                for k, gk in g.groupby(cm):
                    ax.plot([k - 0.3, k + 0.3], [gk[rm].median()] * 2, color="C3", lw=2)
            c = corr[(corr.reporter == rep) & (corr.treatment == trt) & (corr.rate == rm)
                     & (corr.cluster_metric == cm) & (corr.test == "Spearman")]
            txt = f"rho={c.stat.iloc[0]:.2f}, p={c.p.iloc[0]:.2g}, n={c.n.iloc[0]}" if len(c) else ""
            ax.set_title(f"{rep} {g.direction.iloc[0] if len(g) else ''} | {trt}\n{txt}", fontsize=9)
            if cm.endswith("_px"):
                ax.set_xscale("log")
                ax.xaxis.set_major_formatter(matplotlib.ticker.ScalarFormatter())
                ax.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
            ax.set_xlabel(f"{cm} (within {RADIUS_PX} px)" if cm == "n_neighbors" else cm)
            ax.set_ylabel(rm)
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)


def plot_superplot(t, rm, path):
    """Rate by position in the cluster (isolated / edge / middle) within each treatment: faint per-cell
    dots, unfilled boxes, big dots for per experiment+replicate medians."""
    reps = sorted(t.reporter.unique())
    trts = treatment_order(t)
    fig, axes = plt.subplots(1, len(reps), figsize=(1.2 * len(POSITIONS) * len(trts) + 1, 3.8), squeeze=False)
    rng = np.random.default_rng(RANDOM_SEED)
    units = sorted(t.unit.unique())
    for ax, rep in zip(axes[0], reps):
        g = t[t.reporter == rep]
        g = g[g.good_fit] if rm == "k_per_min" else g[g.bright]
        pos, labels = 0, []
        for trt in trts:
            for ci, cl in enumerate(POSITIONS):
                v = g[(g.treatment == trt) & (g.cluster_position == cl)]
                ax.scatter(pos + rng.uniform(-0.15, 0.15, len(v)), v[rm], s=8, alpha=0.3,
                           color=f"C{ci}", lw=0)
                if len(v):
                    ax.boxplot(v[rm], positions=[pos], widths=0.5, showfliers=False,
                               medianprops=dict(color="k"), boxprops=dict(color="k"))
                for u, vu in v.groupby("unit"):
                    ax.scatter(pos, vu[rm].median(), s=70, edgecolor="k", zorder=3,
                               color=f"C{len(POSITIONS) + units.index(u)}", label=u)
                labels.append(f"{trt}\n{cl}\nn={len(v)}")
                pos += 1
        ax.set_xticks(range(pos))
        ax.set_xticklabels(labels, fontsize=7)
        ax.set_title(f"{rep} ({g.direction.iloc[0]})")
        ax.set_ylabel(rm)
    h, l = axes[0, 0].get_legend_handles_labels()
    if h:
        uniq = dict(zip(l, h))
        axes[0, -1].legend(uniq.values(), uniq.keys(), fontsize=7, frameon=False)
    fig.suptitle(f"cluster = nuclei within {CLUSTER_GAP_PX} px (edge to edge); middle = neighbours within "
                 f"{RADIUS_PX} px on all sides (no empty sector >= {EDGE_ANGLE_DEG} deg)", fontsize=8)
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)


def analyse(t, out):
    os.makedirs(out, exist_ok=True)
    t.to_csv(os.path.join(out, "nucleus_clustering_per_cell.csv"), index=False)
    t = t.assign(orig_treatment=t.treatment)
    if POOL_TREATMENTS:
        t = pd.concat([t, t.assign(treatment=ALL_LABEL, dose_mM=np.inf)], ignore_index=True)
    corr = correlations(t)
    corr.to_csv(os.path.join(out, "clustering_vs_rate_stats.csv"), index=False)
    for rm in RATE_METRICS:
        for cm in SCATTER_METRICS:
            plot_scatter(t, corr, rm, cm, os.path.join(out, f"{cm}_vs_{rm}.png"))
        plot_superplot(t, rm, os.path.join(out, f"cluster_position_vs_{rm}.png"))
    s = corr[(corr.rate == "flux_pct_per_min") & corr.cluster_metric.isin(SCATTER_METRICS + ["cluster_position"])]
    print(s[["reporter", "treatment", "cluster_metric", "test", "n", "stat", "p"]][s.treatment == ALL_LABEL].round(3).to_string(index=False))
    return corr


def main():
    allt = []
    for e in EXPERIMENTS:
        print(f"Experiment: {e['rates']}")
        nuc = measure_folder(e["masks"], e.get("mask_suffix", "_640_cp_masks.tif"))
        t = load_rates(e["rates"], e["name"]).merge(nuc, on=["reporter", "file_id", "position", "roi"], how="left")
        miss = t.n_neighbors.isna().sum()
        if miss:
            print(f"  {miss} rate cells had no matching nucleus label (dropped)")
        t = t.dropna(subset=["n_neighbors"]).astype({"clustered": bool, "near_edge": bool})
        analyse(t, e["out"])
        allt.append(t)
    if POOLED_OUT and len(allt) > 1:
        analyse(pd.concat(allt, ignore_index=True), POOLED_OUT)


if __name__ == "__main__":
    main()
