import numpy as np
from scipy import ndimage as ndi
from skimage import filters, morphology, segmentation

def soma_frame(img, nuc, expressing, frac=0.5, max_dist=40, ring=(2, 8), sigma=1.5, ring_stat=50, bg=None, open_r=2):
    sm = filters.gaussian(img.astype(np.float32), sigma, preserve_range=True)
    if bg is None: bg = np.percentile(sm, 20)
    nucmask = nuc > 0
    keep = [l for l in np.unique(nuc) if l in expressing]
    fg = np.zeros(nuc.shape, bool)
    for l in keep:
        n = nuc == l
        d = ndi.distance_transform_edt(~n)
        r = (d > ring[0]) & (d <= ring[1]) & ~nucmask
        if r.sum() < 10: continue
        thr = bg + frac * (np.percentile(sm[r], ring_stat) - bg)
        cand = ((sm > thr) & (d <= max_dist) & ~nucmask) | n
        cc, _ = ndi.label(cand)
        own = np.unique(cc[n]); fg |= np.isin(cc, own[own > 0])
    seeds = np.where(np.isin(nuc, keep), nuc, 0).astype(np.int32)
    dist = ndi.distance_transform_edt(seeds == 0)
    cells = segmentation.watershed(dist, markers=seeds, mask=fg)
    out = np.zeros(nuc.shape, np.uint16)
    for l in keep:
        m = morphology.opening(cells == l, morphology.disk(open_r)) | (nuc == l)
        cc, k = ndi.label(m)
        if k > 1:
            own = np.unique(cc[nuc == l]); m = np.isin(cc, own[own > 0])
        out[m & (out == 0)] = l
    return out

def expressing_cells(stack, nucs, min_nuc=25):
    bg = np.median([np.percentile(f, 20) for f in stack])
    vals = {}
    for t in range(len(stack)):
        for l in np.unique(nucs[t])[1:]:
            vals.setdefault(l, []).append(stack[t][nucs[t] == l].mean() - bg)
    return {l for l, v in vals.items() if np.median(v) > min_nuc}

# Defaults tuned on the user's hand-drawn ROIs (XY1 frame 1, 7 cells; checked against a second drawing at frame 33).
PARAMS = dict(frac=0.45, max_dist=25, ring_stat=75, open_r=5)

def cyto_stack(stack, nucs, rim=2, **kw):
    """Cytosol labels (same IDs as the nuclear labels) for every frame.
    Only cells whose nucleus is clearly above background in 561 (expressing cells) get a mask.
    Soma = connected pixels around the nucleus brighter than bg + frac*(perinuclear 75th pct - bg),
    within max_dist px; touching somas are split by distance to each nucleus. Cytosol = soma minus
    the nucleus dilated by `rim` px."""
    p = {**PARAMS, **kw}
    ex = expressing_cells(stack, nucs)
    bg = np.median([np.percentile(filters.gaussian(f.astype(np.float32), 1.5, preserve_range=True), 20) for f in stack])
    out = np.zeros(stack.shape, np.uint16)
    for t in range(stack.shape[0]):
        soma = soma_frame(stack[t], nucs[t], ex, bg=bg, **p)
        nd = morphology.dilation(nucs[t] > 0, morphology.disk(rim))
        soma[nd] = 0
        out[t] = soma
    return out

def puncta_mask(img, cyto, bg, k=0.5, radius=3, min_th=25, grow=1):
    """Very bright puncta inside each cell's cytosol: pixels on small spots whose white top-hat
    (disk `radius`, after a 1 px Gaussian) exceeds max(min_th, k * the cell's baseline), where the
    baseline is the 30th percentile of that cell's cytosol above bg. Grown by `grow` px to drop the halo."""
    sm = filters.gaussian(img.astype(np.float32), 1, preserve_range=True)
    th = morphology.white_tophat(sm, morphology.disk(radius))
    out = np.zeros(img.shape, bool)
    for l in np.unique(cyto)[1:]:
        m = cyto == l
        base = max(np.percentile(sm[m], 30) - bg, 5)
        out |= m & (th > max(min_th, k * base))
    if grow: out = morphology.dilation(out, morphology.disk(grow))
    return out

def remove_puncta(stack, cyto, **kw):
    bg = np.median([np.percentile(filters.gaussian(f.astype(np.float32), 1.5, preserve_range=True), 20) for f in stack])
    out = cyto.copy()
    for t in range(stack.shape[0]):
        out[t][puncta_mask(stack[t], cyto[t], bg, **kw)] = 0
    return out

if __name__ == "__main__":
    # usage: python cyto_mask.py <folder> [--overwrite]
    # writes <name>_561_cyto_mask.tif next to each <name>_561.tif, using <name>_640_cp_masks.tif as nuclei,
    # and <name>_561_cyto_mask_nopuncta.tif with the very bright puncta removed
    import sys, glob, os, tifffile
    folder = sys.argv[1]; overwrite = "--overwrite" in sys.argv
    for f561 in sorted(glob.glob(os.path.join(folder, "*_561.tif"))):
        base = f561[:-len("_561.tif")]
        out = base + "_561_cyto_mask.tif"
        if os.path.exists(out) and not overwrite:
            print("skip, exists:", out); continue
        stack = tifffile.imread(f561)
        masks = cyto_stack(stack, tifffile.imread(base + "_640_cp_masks.tif"))
        tifffile.imwrite(out, masks[:, None, None], imagej=True, metadata={"axes": "TZCYX"}, compression="zlib")
        clean = remove_puncta(stack, masks)
        out2 = base + "_561_cyto_mask_nopuncta.tif"
        tifffile.imwrite(out2, clean[:, None, None], imagej=True, metadata={"axes": "TZCYX"}, compression="zlib")
        print("wrote", out, "and", out2)
