# nucleus-clustering

`nucleus_clustering.py` asks whether cells whose nuclei sit in clusters transport faster or slower.

It reads the Cellpose nuclear label masks (one label per nucleus; label = Fiji ROI id) and, for one
frame (`MASK_FRAME`, default the first, pre-activation), measures per nucleus:

| column | meaning |
|---|---|
| `n_neighbors` | other nuclei with a centroid within `RADIUS_PX` (default 100 px, about 2 nuclear diameters) |
| `nn_dist_px` | centroid distance to the nearest nucleus |
| `nn_gap_px` | edge-to-edge gap to the nearest nucleus |
| `cluster_size` | nuclei in the group linked by gaps <= `CLUSTER_GAP_PX` (default 10 px); 1 = isolated |
| `clustered` | `cluster_size >= 2` |
| `n_nuclei_field` | nuclei in the field (density confound) |
| `near_edge` | centroid within `RADIUS_PX` of the image border (neighbour count is an underestimate) |

Labels smaller than `MIN_NUC_AREA` px are fragments and are ignored. The metrics are joined to
`rates_per_cell.csv` from `nct_analysis.py` by reporter, file id, XY position and ROI.

Stats per reporter x treatment: Spearman of each rate against each clustering metric, and
Mann-Whitney of clustered vs isolated cells. Flux and slope use bright cells; k uses well-fit cells.

Each `EXPERIMENTS` entry is one experiment (mask folder + nct_results folder). Outputs go to
`<nct_results>/nucleus_clustering/`:
- `nucleus_clustering_per_cell.csv`: rates plus clustering metrics per cell
- `clustering_vs_rate_stats.csv`
- `n_neighbors_vs_<rate>.png`: rate vs neighbour count, red bars = medians
- `clustered_vs_isolated_<rate>.png`: superplot, faint per-cell dots, unfilled boxes, big dots = per experiment+replicate medians

Needs numpy, pandas, scipy, matplotlib, tifffile.
