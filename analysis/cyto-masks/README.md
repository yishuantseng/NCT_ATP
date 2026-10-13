# cyto-masks

`cyto_mask.py` builds per-frame cytosol masks for the 561 nm (PA) channel from Cellpose nuclear masks.

For each `<name>_561.tif` in a folder it reads `<name>_640_cp_masks.tif` as nuclei and writes:

- `<name>_561_cyto_mask.tif`: cytosol labels (same IDs as the nuclei), only for cells that express the reporter
- `<name>_561_cyto_mask_nopuncta.tif`: the same masks with very bright puncta removed

```
python cyto_mask.py <folder> [--overwrite]
```

Default parameters (`PARAMS`) were tuned on hand-drawn ROIs. Requires numpy, scipy, scikit-image and tifffile.
