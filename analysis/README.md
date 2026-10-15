# analysis

Analysis scripts for the MN no-glucose / 2DG experiments.

| Folder | Script | What it does |
|---|---|---|
| `nct-analysis/` | `nct_analysis.py` | Nuclear import/export rates (bleach correction, fits, per-dose superplots, stats) |
| `cyto-masks/` | `cyto_mask.py` | Cytosol masks for the 561 channel from Cellpose nuclear masks |
| `atp-correlation/` | `atp_transport_correlation.py` | ATP/ADP ratio vs transport rates by dose |
| `nucleus-clustering/` | `nucleus_clustering.py` | Nuclei clustering (neighbours, gaps, cluster size) vs transport rates |

Each script has its parameters in a block at the top of the file. Data and outputs are not stored
here (see `.gitignore`); keep them next to the scripts locally.
