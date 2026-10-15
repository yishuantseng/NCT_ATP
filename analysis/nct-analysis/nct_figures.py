"""
Fine-tune figures from an existing nct_analysis.py run (usually the pooled one).

Run nct_analysis.py first. This script only reads its outputs (rates_per_cell.csv and
timecourses_long.csv) and redraws the figures listed in FIGURES, each with its own treatments,
order, control, axis limits, size and labels. Statistics are recomputed for the treatments in
each figure, so stars always compare against that figure's control.

Edit the parameters below and run:  python nct_figures.py
It needs nct_analysis.py in the same folder (plot style and statistics are shared).
"""

# =============================================================================
# PARAMETERS
# =============================================================================
# Folder with the analysis outputs (relative paths are next to this script).
RESULTS_FOLDER = "nct_results_pooled"
OUTPUT_FOLDER = "figures"            # created inside RESULTS_FOLDER
FORMATS = ["png", "pdf"]             # one file per format
DPI = 300

# Defaults used by every figure unless the figure sets its own value.
DEFAULTS = dict(
    stars=True,             # stars vs control (Holm-adjusted)
    stats_unit="auto",      # "cells", "replicate" or "auto" (replicate when every group has >= 3 units)
    figsize=None,           # (width, height) in inches; None = automatic
    ylim=None,              # (min, max); None = automatic. Use None inside, e.g. (0, None), to fix one end
    xlim=None,
    title=None,             # None = automatic title; "" = no title
    ylabel=None,            # None = automatic
    xlabel=None,
    labels=None,            # x tick labels replacing the treatment names, e.g. {"0mM 2DG": "0"}
    show_n=True,            # "(n=…)" under each x tick
    error="sem",            # timecourse: "sem" or "sd"
)

# One entry per figure. Keys:
#   name        file name (without extension)
#   kind        "metric"       dot plot (superplot) of one per-cell value
#               "timecourse"   mean trace ± error per treatment
#               "rate_scatter" export rate (x) vs import rate (y) per treatment
#               "pct_well_fit" % of bright cells that are well fit
#   metric      for "metric" / "rate_scatter": "flux", "slope", "k" or "half_time"
#   reporter    for "metric" / "timecourse" / "pct_well_fit": e.g. 110
#   treatments  list of treatments, in plot order (exact names as in rates_per_cell.csv)
#   control     reference for the stars (default: first treatment)
#   plus any key from DEFAULTS
FIGURES = [
    dict(name="flux_110_2DG", kind="metric", metric="flux", reporter=110,
         treatments=["0mM 2DG", "2mM 2DG", "6mM 2DG", "10mM 2DG", "50mM 2DG"], ylim=(-0.05, 0.25)),
    dict(name="flux_111_2DG", kind="metric", metric="flux", reporter=111,
         treatments=["0mM 2DG", "2mM 2DG", "6mM 2DG", "10mM 2DG", "50mM 2DG"], ylim=(-0.2, 1.2)),
    dict(name="timecourse_110_2DG", kind="timecourse", reporter=110,
         treatments=["0mM 2DG", "10mM 2DG", "50mM 2DG"], ylim=(0.9, 1.4)),
    dict(name="export_vs_import_flux", kind="rate_scatter", metric="flux",
         treatments=["0mM 2DG", "2mM 2DG", "6mM 2DG", "10mM 2DG", "50mM 2DG"], xlim=(0, None), ylim=(0, None)),
    dict(name="well_fit_110", kind="pct_well_fit", reporter=110,
         treatments=["0mM 2DG", "2mM 2DG", "6mM 2DG", "10mM 2DG", "50mM 2DG"]),
]

# =============================================================================
# END OF PARAMETERS
# =============================================================================

import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import nct_analysis as nct      # noqa: E402  (shared style, colours and statistics)
import matplotlib.pyplot as plt  # noqa: E402

METRICS = {   # name: (column, cells, default y label, number format)
    "flux": ("flux_pct_per_min", "bright", None, "{:.2f}"),
    "slope": ("slope_pct_per_min", "bright", None, "{:.2f}"),
    "k": ("k_per_min", "well-fit", "k (min⁻¹)", "{:.2f}"),
    "half_time": ("t_half_s", "well-fit", "t½ (s)", "{:.0f}"),
}
METRIC_TITLES = {"flux": "initial flux k × amplitude", "slope": "model-free initial slope",
                 "k": "rate constant k", "half_time": "half-time"}


def here(path):
    return path if os.path.isabs(path) else os.path.join(os.path.dirname(os.path.abspath(__file__)), path)


def opt(fig, key):
    return fig.get(key, DEFAULTS.get(key))


def check_treatments(spec, available):
    conds = spec.get("treatments") or sorted(available, key=nct.first_number)
    missing = [c for c in conds if c not in available]
    if missing:
        raise ValueError(f"{spec['name']}: treatments not in the results: {missing}. "
                         f"Available: {sorted(available)}")
    control = spec.get("control", conds[0])
    if control not in conds:
        raise ValueError(f"{spec['name']}: control '{control}' is not in its treatments")
    return conds, control


def stats_for(cells, conds, control, spec):
    """Tests vs control for this figure's treatments; sets the level the stars use."""
    nct.STATS_UNIT = opt(spec, "stats_unit")
    _, tests = nct.summarise(cells, conds, control)
    nct.STARS_LEVEL = nct.stats_level(cells, conds)
    nct.SHOW_STARS = bool(opt(spec, "stars"))
    return tests


def tick_labels(spec, conds, counts=None):
    names = opt(spec, "labels") or {}
    out = [names.get(c, c) for c in conds]
    if counts is not None and opt(spec, "show_n"):
        out = [f"{o}\n(n={n})" for o, n in zip(out, counts)]
    return out


def finish(ax, spec, title, ylabel, xlabel=None):
    t = opt(spec, "title")
    ax.set_title(title if t is None else t, loc="left", fontsize=10)
    y = opt(spec, "ylabel")
    ax.set_ylabel(ylabel if y is None else y)
    x = opt(spec, "xlabel")
    if x is not None or xlabel is not None:
        ax.set_xlabel(xlabel if x is None else x)
    for key, setter in [("ylim", ax.set_ylim), ("xlim", ax.set_xlim)]:
        lim = opt(spec, key)
        if lim is not None:
            setter(*lim)


def save(fig, out, name):
    for fmt in FORMATS:
        fig.savefig(os.path.join(out, f"{name}.{fmt}"), dpi=DPI)
    plt.close(fig)


def fig_metric(spec, cells, out):
    rep = spec["reporter"]
    col, which, ylabel, fmt = METRICS[spec["metric"]]
    x = cells[(cells.reporter == rep) & (cells.good_fit if which == "well-fit" else cells.bright)]
    conds, control = check_treatments(spec, set(x.treatment))
    x = x[x.treatment.isin(conds)]
    direction = x.direction.iloc[0]
    tests = stats_for(cells[cells.reporter == rep], conds, control, spec)
    if col in nct.RATE_COLS:
        x = x.assign(**{col: x[col] * nct.rate_scale(direction)})
        ylabel = nct.rate_label(direction)
    fig, ax = plt.subplots(figsize=opt(spec, "figsize") or (max(6, 1 + 1.0 * len(conds)), 4.8))
    nct.strip(ax, x, col, conds, nct.shades(direction, len(conds)), fmt=fmt)
    ax.set_xticks(range(len(conds)), tick_labels(spec, conds, [(x.treatment == c).sum() for c in conds]),
                  fontsize=8)
    if col in nct.RATE_COLS:
        ax.axhline(0, color="#999", lw=.8, ls="--")
    finish(ax, spec, f"{rep} ({direction}): {METRIC_TITLES[spec['metric']]}, {which} cells", ylabel)
    if opt(spec, "ylim") is None:
        lo, hi = ax.get_ylim()
        ax.set_ylim(lo, hi + 0.12 * (hi - lo))           # headroom for stars
    fig.tight_layout(rect=(0, 0.03, 1, 1))
    nct.add_stars(ax, tests, rep, col, conds, control)
    save(fig, out, spec["name"])


def fig_pct_well_fit(spec, cells, out):
    rep = spec["reporter"]
    x = cells[cells.bright & (cells.reporter == rep)]
    conds, control = check_treatments(spec, set(x.treatment))
    direction = x.direction.iloc[0]
    tests = stats_for(cells[cells.reporter == rep], conds, control, spec)
    k = np.array([x[x.treatment == c].good_fit.sum() for c in conds])
    n = np.array([(x.treatment == c).sum() for c in conds])
    pct = 100 * k / np.maximum(n, 1)
    fig, ax = plt.subplots(figsize=opt(spec, "figsize") or (max(6, 1 + 1.0 * len(conds)), 4.8))
    ax.bar(range(len(conds)), pct, 0.6, color=nct.shades(direction, len(conds)))
    for i in range(len(conds)):
        ax.text(i, pct[i] + 1.5, f"{pct[i]:.0f}%\n({k[i]}/{n[i]})", ha="center", va="bottom", fontsize=8)
    ax.set_xticks(range(len(conds)), tick_labels(spec, conds), fontsize=8)
    ax.set_ylim(0, 125)
    ax.set_yticks(range(0, 101, 20))
    nct.style(ax)
    finish(ax, spec, f"{rep} ({direction}): well-fit cells", "well-fit cells (% of bright cells)")
    fig.tight_layout(rect=(0, 0.03, 1, 1))
    nct.add_stars(ax, tests, rep, "pct_good_fit", conds, control)
    save(fig, out, spec["name"])


def fig_timecourse(spec, traces, out):
    rep = spec["reporter"]
    x = traces[(traces.reporter == rep) & traces.bright]
    conds, _ = check_treatments(spec, set(x.treatment))
    err = opt(spec, "error")
    cols = nct.cond_colors(conds)
    names = opt(spec, "labels") or {}
    t0 = (nct.ACTIVATION_FRAME - 1) * nct.INTERVAL_S
    fig, ax = plt.subplots(figsize=opt(spec, "figsize") or (6.4, 4.3))
    for cond in conds:
        g = x[x.treatment == cond]
        m = g.groupby("time").norm_corr.agg(["mean", "std" if err == "sd" else "sem"])
        e = m.iloc[:, 1]
        ax.fill_between(m.index, m["mean"] - e, m["mean"] + e, color=cols[cond], alpha=.18, lw=0)
        ax.plot(m.index, m["mean"], color=cols[cond], lw=2,
                label=f"{names.get(cond, cond)} (n={g.cell_id.nunique()})")
    ax.axhline(1, color="#999", lw=.8, ls="--")
    ax.axvspan(0, t0, color="#eee", zorder=0)
    ax.legend(frameon=False, fontsize=8, loc="upper left", bbox_to_anchor=(1.01, 1))
    nct.style(ax)
    direction = nct.DIRECTIONS.get(rep, "")
    finish(ax, spec, f"Reporter {rep} ({direction}): mean ± {err.upper()}, bright cells",
           "nuclear signal / baseline", "time (s)")
    fig.tight_layout()
    save(fig, out, spec["name"])


def fig_rate_scatter(spec, cells, out):
    col, which, _, _ = METRICS[spec["metric"]]
    c = cells[cells.good_fit if which == "well-fit" else cells.bright].copy()
    imp = [r for r in sorted(c.reporter.unique()) if nct.DIRECTIONS.get(r) == "import"]
    exp = [r for r in sorted(c.reporter.unique()) if nct.DIRECTIONS.get(r) == "export"]
    if not imp or not exp:
        raise ValueError(f"{spec['name']}: needs one import and one export reporter")
    both = set(c[c.reporter == imp[0]].treatment) & set(c[c.reporter == exp[0]].treatment)
    conds, _ = check_treatments(spec, both)
    if col in nct.RATE_COLS:
        c["v"] = c[col] * np.where(c.direction == "export", -1, 1) * (0.01 if nct.RATE_UNIT == "per_min" else 1)
        unit = "min⁻¹" if nct.RATE_UNIT == "per_min" else "% of baseline per min"
    else:
        c["v"] = c[col]
        unit = "min⁻¹" if spec["metric"] == "k" else "s"
    xi, xe = c[c.reporter == imp[0]], c[c.reporter == exp[0]]
    units = pd.concat([xe.groupby(["treatment", "unit"]).v.median().rename("x"),
                       xi.groupby(["treatment", "unit"]).v.median().rename("y")], axis=1).dropna().reset_index()
    by_units = all((units.treatment == cd).sum() >= 2 for cd in conds)
    cols = nct.cond_colors(conds)
    names = opt(spec, "labels") or {}
    fig, ax = plt.subplots(figsize=opt(spec, "figsize") or (6.4, 4.6))
    for cd in conds:
        u = units[units.treatment == cd]
        ax.scatter(u.x, u.y, s=40, color=cols[cd], alpha=.35, lw=0, zorder=2)
        if by_units:
            mx, my, sx, sy = u.x.mean(), u.y.mean(), u.x.std(), u.y.std()
        else:
            ve, vi = xe[xe.treatment == cd].v, xi[xi.treatment == cd].v
            mx, my, sx, sy = ve.median(), vi.median(), ve.std(), vi.std()
        ax.errorbar(mx, my, xerr=sx, yerr=sy, fmt="o", ms=8, color=cols[cd], elinewidth=1.3, capsize=4,
                    zorder=3, label=names.get(cd, cd))
    ax.legend(frameon=False, fontsize=8, loc="upper left", bbox_to_anchor=(1.01, 1))
    nct.style(ax)
    ax.grid(axis="x", color="#eee")
    name = METRIC_TITLES[spec["metric"]]
    rev = ", sign reversed" if col in nct.RATE_COLS else ""
    centre = "mean of experiment+replicate medians" if by_units else "median of cells"
    finish(ax, spec, f"Export vs import: {centre} ± SD", f"import rate, {name} ({unit})",
           f"export rate, {name}{rev} ({unit})")
    fig.tight_layout()
    save(fig, out, spec["name"])


def main():
    folder = here(RESULTS_FOLDER)
    out = os.path.join(folder, OUTPUT_FOLDER)
    os.makedirs(out, exist_ok=True)
    cells = pd.read_csv(os.path.join(folder, "rates_per_cell.csv"), dtype={"roi": str, "replicate": str})
    traces = None
    for spec in FIGURES:
        kind = spec["kind"]
        try:
            if kind == "metric":
                fig_metric(spec, cells, out)
            elif kind == "pct_well_fit":
                fig_pct_well_fit(spec, cells, out)
            elif kind == "rate_scatter":
                fig_rate_scatter(spec, cells, out)
            elif kind == "timecourse":
                if traces is None:
                    traces = pd.read_csv(os.path.join(folder, "timecourses_long.csv"),
                                         usecols=["cell_id", "time", "reporter", "treatment", "bright", "norm_corr"])
                fig_timecourse(spec, traces, out)
            else:
                raise ValueError(f"unknown kind '{kind}'")
            print(f"  {spec['name']}")
        except Exception as e:
            print(f"  {spec.get('name', '?')}: skipped ({e})")
    print(f"Figures written to {out}")


if __name__ == "__main__":
    main()
