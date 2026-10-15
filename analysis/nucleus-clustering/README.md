# nucleus-clustering

`nucleus_clustering.py` asks whether cells whose nuclei sit in clusters transport faster or slower.

It reads the Cellpose nuclear label masks (one label per nucleus; label = Fiji ROI id) and, for one
frame (`MASK_FRAME`, default the first, pre-activation), measures per nucleus:

| column | meaning |
|---|---|
| `mean_dist_px` | mean centroid distance to the `K_NEAREST` (default 3) nearest nuclei |
| `n_neighbors` | other nuclei with a centroid within `RADIUS_PX` (default 100 px, about 2 nuclear diameters) |
| `nn_dist_px` | centroid distance to the nearest nucleus |
| `nn_gap_px` | edge-to-edge gap to the nearest nucleus |
| `cluster_size` | nuclei in the group linked by gaps <= `CLUSTER_GAP_PX` (default 10 px); 1 = isolated |
| `clustered` | `cluster_size >= 2` |
| `cluster_position` | `isolated` (cluster_size 1), `middle` (neighbours within `RADIUS_PX` on all sides: no empty sector >= `EDGE_ANGLE_DEG`, default 180°, and at least 3 neighbours) or `edge` |
| `n_nuclei_field` | nuclei in the field (density confound) |
| `near_edge` | centroid within `RADIUS_PX` of the image border (neighbour count is an underestimate) |

Labels smaller than `MIN_NUC_AREA` px are fragments and are ignored. The metrics are joined to
the rates by reporter, file id, XY position and ROI. The rates file can be `rates_per_cell.csv` from
`nct_analysis.py` or `rates_nuc_vs_nc_per_cell.csv` from `nc_rates.py` (readout chosen by
`NC_RATES_READOUT`; wrong-direction cells dropped; well-fit = r2 >= `MIN_R2` and k below its cap).

Stats per reporter x treatment: Spearman of each rate against each clustering metric, Kruskal-Wallis
across isolated/edge/middle, and Mann-Whitney middle vs edge. With `POOL_TREATMENTS = True` every
reporter also gets an `all` group (all treatments together; the pooled scatter colours cells by
treatment), plus a Spearman with each cell's rate centred on its own treatment median, so a dose
effect on both crowding and rate can't drive the pooled correlation. Flux and slope use bright cells; k uses well-fit cells.

Each `EXPERIMENTS` entry is one experiment (name, mask folder + suffix, rates CSV, output folder).
Outputs per experiment:
- `nucleus_clustering_per_cell.csv`: rates plus clustering metrics per cell
- `clustering_vs_rate_stats.csv`
- `mean_dist_px_vs_<rate>.png`: rate vs average distance to the nearest nuclei (log x)
- `n_neighbors_vs_<rate>.png`: rate vs neighbour count, red bars = medians
- `cluster_position_vs_<rate>.png`: isolated / edge / middle superplot, faint per-cell dots, unfilled boxes, big dots = per experiment+replicate medians

Needs numpy, pandas, scipy, matplotlib, tifffile.
