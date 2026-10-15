"""Nucleus clustering vs nucleocytoplasmic transport rate.

For every nucleus in the Cellpose nuclear masks (one label per nucleus, label = Fiji ROI id) this
measures how clustered it is with its neighbours, then joins that to the per-cell transport rates
from nct_analysis.py (rates_per_cell.csv) and asks whether clustered cells transport faster or slower,
within each reporter x treatment.

Clustering metrics (all from one mask frame, pre-activation by default; pixels):
  n_neighbors    other nuclei whose centroid lies within RADIUS_PX of this one's centroid
  nn_dist_px     centroid distance to the nearest other nucleus
  nn_gap_px      edge-to-edge gap to the nearest other nucleus (0-1 = touching)
  cluster_size   size of the group of nuclei linked by gaps <= CLUSTER_GAP_PX (1 = isolated)
  clustered      cluster_size >= 2
  n_nuclei_field nuclei in the whole field (a density confound to keep in mind)
  near_edge      centroid within RADIUS_PX of the image border (neighbour count is an underestimate)

Each EXPERIMENTS entry is one experiment: its mask folder and its nct_results folder.
Outputs go to <nct_results>/nucleus_clustering/ for each experiment, plus a pooled folder if
there is more than one experiment.
"""
import os
import re

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import tifffile
from scipy import ndimage, stats
from scipy.sparse.csgraph import connected_components

# ---------------------------------------------------------------- parameters
EXPERIMENTS = [
    dict(masks="/Volumes/holtl02lab/holtl02labspace/Holt_Lab_Members/Yi-Shuan/microscope/2026/10-12/20261010_MN_110,111/combine/DNA/erode",
         nct_results="/Volumes/holtl02lab/holtl02labspace/Holt_Lab_Members/Yi-Shuan/microscope/2026/10-12/20261010_MN_110,111/analysis/nct_results"),
]
POOLED_OUT = None            # folder for pooled outputs when there are several experiments (None = skip)
MASK_SUFFIX = "_640_mask.tif"
MASK_FRAME = 0               # 0-based frame of the mask stack to measure (0 = first, pre-activation)
MIN_NUC_AREA = 300           # px; smaller labels are fragments and are not counted as neighbours
RADIUS_PX = 100              # neighbourhood radius for n_neighbors (~2 nuclear diameters)
CLUSTER_GAP_PX = 10          # nuclei closer than this (edge to edge) belong to the same cluster
EXCLUDE_REPORTERS = [185]
RATE_METRICS = ["flux_pct_per_min", "slope_pct_per_min", "k_per_min"]   # k uses well-fit cells only
CLUSTER_METRICS = ["n_neighbors", "nn_gap_px", "nn_dist_px", "cluster_size", "n_nuclei_field"]
SCATTER_METRIC = "n_neighbors"
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
    h, w = mask.shape
    edge = np.minimum.reduce([cents[:, 0], cents[:, 1], h - 1 - cents[:, 0], w - 1 - cents[:, 1]])
    one = len(ids) > 1
    return pd.DataFrame(dict(
        roi=ids.astype(int), nuc_area_px=areas, cy=cents[:, 0], cx=cents[:, 1],
        n_neighbors=(d <= RADIUS_PX).sum(1),
        nn_dist_px=d.min(1) if one else np.nan, nn_gap_px=gap.min(1) if one else np.nan,
        cluster_size=size, clustered=size >= 2, n_nuclei_field=len(ids), near_edge=edge < RADIUS_PX))


def measure_folder(folder):
    rows = []
    for f in sorted(os.listdir(folder)):
        if not f.endswith(MASK_SUFFIX) or field_key(f) is None:
            continue
        mask = tifffile.imread(os.path.join(folder, f), key=MASK_FRAME)
        t = measure_field(mask)
        rep, fid, pos = field_key(f)
        rows.append(t.assign(mask_file=f, reporter=rep, file_id=fid, position=pos))
        print(f"  {f}: {len(t)} nuclei")
    return pd.concat(rows, ignore_index=True)


def load_rates(path):
    r = pd.read_csv(path, dtype={"roi": str})
    r = r[~r.reporter.isin(EXCLUDE_REPORTERS)].copy()
    r["roi"] = r.roi.astype(int)
    keys = r.image.map(field_key)
    r["file_id"] = [k[1] for k in keys]
    r["position"] = [k[2] for k in keys]
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
            for cm in CLUSTER_METRICS + ["clustered"]:
                x = cells[[cm, rm]].dropna()
                if len(x) < 5:
                    continue
                if cm == "clustered":
                    a, b = x[x[cm]][rm], x[~x[cm]][rm]
                    if min(len(a), len(b)) < 3:
                        continue
                    p = stats.mannwhitneyu(a, b).pvalue
                    out.append(dict(reporter=rep, direction=direc, treatment=trt, rate=rm, cluster_metric=cm,
                                    test="Mann-Whitney clustered vs isolated", n=len(x), n_clustered=len(a),
                                    median_clustered=a.median(), median_isolated=b.median(), stat=np.nan, p=p))
                else:
                    rho, p = stats.spearmanr(x[cm], x[rm])
                    out.append(dict(reporter=rep, direction=direc, treatment=trt, rate=rm, cluster_metric=cm,
                                    test="Spearman", n=len(x), stat=rho, p=p))
    return pd.DataFrame(out)


def treatment_order(t):
    return list(t.sort_values("dose_mM").treatment.unique()) if "dose_mM" in t else sorted(t.treatment.unique())


def plot_scatter(t, corr, rm, path):
    reps = sorted(t.reporter.unique())
    trts = treatment_order(t)
    fig, axes = plt.subplots(len(reps), len(trts), figsize=(3.2 * len(trts), 3 * len(reps)), squeeze=False)
    rng = np.random.default_rng(RANDOM_SEED)
    for i, rep in enumerate(reps):
        for j, trt in enumerate(trts):
            ax = axes[i, j]
            g = t[(t.reporter == rep) & (t.treatment == trt)]
            g = g[g.good_fit] if rm == "k_per_min" else g[g.bright]
            x = g[SCATTER_METRIC] + rng.uniform(-0.2, 0.2, len(g))
            ax.scatter(x, g[rm], s=10, alpha=0.4, color="0.3", lw=0)
            for k, gk in g.groupby(SCATTER_METRIC):
                ax.plot([k - 0.3, k + 0.3], [gk[rm].median()] * 2, color="C3", lw=2)
            c = corr[(corr.reporter == rep) & (corr.treatment == trt) & (corr.rate == rm)
                     & (corr.cluster_metric == SCATTER_METRIC)]
            txt = f"rho={c.stat.iloc[0]:.2f}, p={c.p.iloc[0]:.2g}, n={c.n.iloc[0]}" if len(c) else ""
            ax.set_title(f"{rep} {g.direction.iloc[0] if len(g) else ''} | {trt}\n{txt}", fontsize=9)
            ax.set_xlabel(f"{SCATTER_METRIC} (within {RADIUS_PX} px)")
            ax.set_ylabel(rm)
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)


def plot_superplot(t, rm, path):
    """Isolated vs clustered nuclei within each treatment: faint per-cell dots, unfilled boxes,
    big dots for per experiment+replicate medians."""
    reps = sorted(t.reporter.unique())
    trts = treatment_order(t)
    fig, axes = plt.subplots(1, len(reps), figsize=(1.6 * 2 * len(trts) + 1, 3.6), squeeze=False)
    rng = np.random.default_rng(RANDOM_SEED)
    units = sorted(t.unit.unique())
    for ax, rep in zip(axes[0], reps):
        g = t[t.reporter == rep]
        g = g[g.good_fit] if rm == "k_per_min" else g[g.bright]
        pos, labels = 0, []
        for trt in trts:
            for cl in (False, True):
                v = g[(g.treatment == trt) & (g.clustered == cl)]
                ax.scatter(pos + rng.uniform(-0.15, 0.15, len(v)), v[rm], s=8, alpha=0.3,
                           color="C1" if cl else "C0", lw=0)
                if len(v):
                    ax.boxplot(v[rm], positions=[pos], widths=0.5, showfliers=False,
                               medianprops=dict(color="k"), boxprops=dict(color="k"))
                for u, vu in v.groupby("unit"):
                    ax.scatter(pos, vu[rm].median(), s=70, edgecolor="k", zorder=3,
                               color=f"C{2 + units.index(u)}", label=u)
                labels.append(f"{trt}\n{'clustered' if cl else 'isolated'}\nn={len(v)}")
                pos += 1
        ax.set_xticks(range(pos))
        ax.set_xticklabels(labels, fontsize=7)
        ax.set_title(f"{rep} ({g.direction.iloc[0]})")
        ax.set_ylabel(rm)
    h, l = axes[0, 0].get_legend_handles_labels()
    if h:
        uniq = dict(zip(l, h))
        axes[0, -1].legend(uniq.values(), uniq.keys(), fontsize=7, frameon=False)
    fig.suptitle(f"clustered = within {CLUSTER_GAP_PX} px of another nucleus", fontsize=9)
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)


def analyse(t, out):
    os.makedirs(out, exist_ok=True)
    t.to_csv(os.path.join(out, "nucleus_clustering_per_cell.csv"), index=False)
    corr = correlations(t)
    corr.to_csv(os.path.join(out, "clustering_vs_rate_stats.csv"), index=False)
    for rm in RATE_METRICS:
        plot_scatter(t, corr, rm, os.path.join(out, f"{SCATTER_METRIC}_vs_{rm}.png"))
        plot_superplot(t, rm, os.path.join(out, f"clustered_vs_isolated_{rm}.png"))
    s = corr[(corr.rate == "flux_pct_per_min") & corr.cluster_metric.isin([SCATTER_METRIC, "clustered"])]
    print(s[["reporter", "treatment", "cluster_metric", "n", "stat", "median_clustered", "median_isolated", "p"]]
          .to_string(index=False) if "median_clustered" in s else s.to_string(index=False))
    return corr


def main():
    allt = []
    for e in EXPERIMENTS:
        print(f"Experiment: {e['nct_results']}")
        nuc = measure_folder(e["masks"])
        rates = load_rates(os.path.join(e["nct_results"], "rates_per_cell.csv"))
        t = rates.merge(nuc, on=["reporter", "file_id", "position", "roi"], how="left")
        miss = t.n_neighbors.isna().sum()
        if miss:
            print(f"  {miss} rate cells had no matching nucleus label (dropped)")
        t = t.dropna(subset=["n_neighbors"])
        analyse(t, os.path.join(e["nct_results"], "nucleus_clustering"))
        allt.append(t)
    if POOLED_OUT and len(allt) > 1:
        analyse(pd.concat(allt, ignore_index=True), POOLED_OUT)


if __name__ == "__main__":
    main()
