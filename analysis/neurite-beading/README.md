# neurite-beading

`neurite_beading.py` quantifies beading (focal swellings / varicosities) of neurites in vLH92
motor-neuron `neurite_*` fields (no somata, dense networks). First used on `20261015_MN_92`.

## Run

Edit the PARAMETERS block at the top (experiment folders, file names, thresholds), then run
`python neurite_beading.py`. Needs numpy, pandas, scipy, scikit-image, tifffile, matplotlib.
Each experiment folder needs `tif/` images and `analysis/fileid_treatment.csv`
(columns `image`, `treatment`, `treat_time`, optional `recovery time` and `replicate`).
Outputs go to `<experiment>/analysis/neurite_beading/`.

## Method

1. The 488 and 405 (`Dapi-ex_Green-em`) images are background-subtracted and summed. vLH92 fills
   the cytosol, so the sum follows neurite volume and does not depend on the ATP/ADP ratio.
2. Neurite mask and segments as in `neurite_atp.py` (ridge filter, threshold, skeleton cut at
   crossings). Blobs wider than ~5.5 µm (cell bodies, clumps) are removed.
3. Beads are parts of a neurite at least 1.5 µm wide (a grey opening with a 1.5 µm disk removes
   anything narrower: thin shafts, crossings of thin neurites, short bright spots) and at least
   1.5x brighter than the shaft around them (median along the skeleton in rings around the bead).
   They must be at most 8 µm long and at most 3x longer than wide (else a thick stretch). Where
   4 or more neurite arms meet (a crossing of thick neurites), the spot must be 2.5x the shaft.
4. Large, nearly circular blobs (aspect < 1.25, ≥ 2 µm) are taken as out-of-focus dead cells or
   debris and are not counted; they are listed in `beading_blobs.csv` and drawn blue in overlays.

The criteria were tuned on beads, dead cells and non-beads marked by hand in 20261015_MN_92.

## Metrics

| column | meaning |
|---|---|
| `beads_per_100um` | beads per 100 µm of neurite (skeleton) length; main readout |
| `bead_area_fraction` | fraction of the neurite mask area inside beads |
| `intensity_cv` | SD / mean of the intensity along the skeleton (threshold-free unevenness) |

Image-level values pool the whole field (total beads / total length). Segment-level values are in
`beading_segments.csv`; every bead is in `beading_beads.csv`. `stats_vs_control.csv` compares each
condition with the time-matched `2.5mM gluc` control (Mann–Whitney U, image and segment level).

## Plots

`superplot_images.png` and `superplot_segments.png`: faint dots per image / segment, unfilled
boxplots, big dots for each experiment + replicate median (replicate defaults to "A").
`overlays/` shows the mask and the detected beads for each field, for checking the detection.
