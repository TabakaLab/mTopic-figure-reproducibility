import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import mtopic
import numpy as np
import pandas as pd
import os

from matplotlib.colors import to_rgb
from matplotlib.ticker import MultipleLocator

INPUT_PATHS = {
    "P22H3K4me3": "../data/P22MouseBrainH3K4me3_trained.h5mu",
    "P22H3K27me3": "../data/P22MouseBrainH3K27me3_trained.h5mu",
    "P22H3K27ac": "../data/P22MouseBrainH3K27ac_trained.h5mu",
}

# panels h, i: the benchmark caches, copied over as they are
SCORE_PATHS = {
    "frame": "../data/benchmark_null_hungarian_5sites_frame.csv",
    "points": "../data/benchmark_null_hungarian_5sites_points.csv",
}

# panel g: one file per condition, written by prepare_E3g_cache.py
G_MAP_PATH = "../data/E3g_{mode}_t{topic}_cont{cont}.csv"

OUT = "../figures/Figure_E3"
os.makedirs(OUT, exist_ok=True)


DEFAULT_PALETTE = {
    "topic_1": "#f4f391",
    "topic_2": "#cad9db",
    "topic_3": "#f52128",
    "topic_4": "#b8cfad",
    "topic_5": "#d74f2d",
    "topic_6": "#623dbc",
    "topic_7": "#4df5a8",
    "topic_8": "#2a3b32",
    "topic_9": "#ad1568",
    "topic_10": "#ee73c1",
    "topic_11": "#c0900a",
    "topic_12": "#07e92e",
    "topic_13": "#6eb7a1",
    "topic_14": "#f7b6d2",
    "topic_15": "#7f7f7f",
    "topic_16": "#c7c7c7",
    "topic_17": "#bcbd21",
    "topic_18": "#dbdb8d",
    "topic_19": "#1abecf",
    "topic_20": "#9edbe5",
    "topic_21": "#3a3b79",
    "topic_22": "#5354a3",
    "topic_23": "#6c6ecf",
    "topic_24": "#9d9ede",
    "topic_25": "#647a39",
    "topic_26": "#8ca352",
    "topic_27": "#1263de",
    "topic_28": "#53e463",
    "topic_29": "#2c7cb3",
    "topic_30": "#c47014",
    "topic_31": "#a181f7",
    "topic_32": "#9e5a19",
    "topic_33": "#201e83",
    "topic_34": "#145351",
    "topic_35": "#7c63ea",
    "topic_36": "#684562",
    "topic_37": "#6c72f2",
    "topic_38": "#f3c541",
    "topic_39": "#2a00d4",
    "topic_40": "#f4c127",
    "topic_41": "#56579b",
    "topic_42": "#63d54b",
    "topic_43": "#8aa802",
    "topic_44": "#252fba",
    "topic_45": "#e90157",
    "topic_46": "#615cd2",
    "topic_47": "#d88583",
    "topic_48": "#c3b306",
    "topic_49": "#93c428",
    "topic_50": "#c71581",
}


# ---------------------------------------------------------------------------
# E3g-i settings
# ---------------------------------------------------------------------------

MODE_LABEL = {"frame": "Frame spike-in", "points": "Points spike-in"}

# contamination in the cache, label, shade toward white. As in the benchmark
# script, the shades run lightest-first over ascending contamination, so 100%
# spike-in (contamination 0) is the light one and 50% the dark one.
LEVELS = [(0.0, "100%", 0.45), (0.5, "50%", 0.0)]

DATASET_COLOR = {
    "P22ATAC": "#E9C46A",
    "P22H3K4me3": "#8FBC5A",
    "P22H3K27me3": "#7B93AB",
    "P22H3K27ac": "#B822B0",
}
DATASET_LABEL = {
    "P22ATAC": "ATAC-RNA",
    "P22H3K4me3": "H3K4me3-RNA",
    "P22H3K27me3": "H3K27me3-RNA",
    "P22H3K27ac": "H3K27ac-RNA",
}

MOD_COLOR = {"rna": "#921F18", "chrom": "#215C37"}
MOD_LABEL = {"rna": "RNA", "chrom": "Chromatin"}

# Set outright rather than fitted, so the same panel is comparable between
# runs. A point outside is clipped and reported rather than silently dropped.
Y_LIM = {
    "ari": (0.65, 0.93),
    "auroc": (0.69, 1.0),
    "beta_cos": (0.9, 1.0),
    "beta_rbo": (0.1, 1.0),
}
Y_STEP = {"ari": 0.1, "auroc": 0.1, "beta_cos": 0.05, "beta_rbo": 0.2}

POINT_SIZE = 13
BOX_EDGE = "black"
BOX_LW = 1.0
BOX_MEDIAN_LW = 1.3

G_TOPICS = [10, 25]
G_DATASET_LABEL = "ATAC-RNA"
G_CMAP = "gnuplot"
G_SPOT_SIZE = 6.6


def make_panel(sample, panel_label):
    mdata = mtopic.read.h5mu(INPUT_PATHS[sample])

    mtopic.pl.dominant_topics(
        mdata,
        x="coords",
        s=28,
        figsize=(8, 4),
        legend=True,
        markerscale=3,
        legend_ncol=3,
        palette=DEFAULT_PALETTE,
        save=f"{OUT}/Figure_{panel_label}_1.png",
    )

    palette = DEFAULT_PALETTE
    figsize = (18, 8)

    topic_order = list(mdata.obsm["topics"].columns)
    gamma = mdata.obsm["topics"][topic_order]

    dominant_topic = gamma.idxmax(axis=1)
    dominant_proportion = gamma.max(axis=1)
    topic_rank = {t: i for i, t in enumerate(topic_order)}

    cell_order = (
        pd.DataFrame(
            {
                "topic": dominant_topic,
                "topic_rank": dominant_topic.map(topic_rank),
                "proportion": dominant_proportion,
            }
        )
        .sort_values(["topic_rank", "proportion"], ascending=[True, False])
        .index
    )
    gamma_sorted = gamma.loc[cell_order].T.loc[topic_order]

    N = gamma_sorted.shape[1]
    K = gamma_sorted.shape[0]

    sorted_topics = dominant_topic.loc[cell_order].values
    color_bar = np.array(
        [plt.matplotlib.colors.to_rgb(palette.get(t, "#bebebe")) for t in sorted_topics]
    )

    fig, (ax_bar, ax_heat) = plt.subplots(
        2,
        1,
        figsize=figsize,
        gridspec_kw={"height_ratios": [0.04, 1], "hspace": 0.01},
    )
    ax_bar.imshow(
        color_bar[np.newaxis, :, :],
        aspect="auto",
        interpolation="none",
        extent=[0, N, 0, 1],
    )
    ax_bar.set_xlim(0, N)
    ax_bar.axis("off")

    im = ax_heat.imshow(
        gamma_sorted.values,
        aspect="auto",
        cmap="gnuplot",
        vmin=0,
        vmax=1,
        interpolation="none",
        extent=[0, N, K, 0],
    )
    ax_heat.set_xlim(0, N)
    ax_heat.set_xticks([])
    ax_heat.set_yticks([])

    handles = [mpatches.Patch(color=palette[t], label=t) for t in topic_order]
    fig.legend(
        handles=handles,
        loc="lower center",
        ncol=min(10, len(topic_order)),
        fontsize=7,
        frameon=False,
        bbox_to_anchor=(0.5, -0.02),
    )
    plt.savefig(f"{OUT}/Figure_{panel_label}_2.png", bbox_inches="tight", dpi=300)
    plt.close(fig)


# ---------------------------------------------------------------------------
# E3g-i helpers
# ---------------------------------------------------------------------------

def _shade(color, frac):
    rgb = np.array(to_rgb(color))
    return tuple(rgb * (1 - frac) + frac)


def _dot(color, label):
    return plt.Line2D([], [], marker="o", ls="", color=color, markersize=6,
                      label=label)


def _g_map(mode, topic, cont):
    path = G_MAP_PATH.format(mode=mode, topic=topic, cont=f"{cont:g}")
    if not os.path.exists(path):
        print(f"  missing {path}")
        return None
    return pd.read_csv(path)


def make_panel_g(panel_label="E3g"):
    """Topic distributions only: the spike-in region is not marked."""
    modes = list(MODE_LABEL)
    rows = [(t, c, lab) for t in G_TOPICS for c, lab, _ in LEVELS]

    fig, axes = plt.subplots(
        len(rows), len(modes),
        figsize=(2.3 * len(modes), 2.3 * len(rows)),
        squeeze=False,
    )
    fig.subplots_adjust(left=0.14, right=0.88, top=0.93, bottom=0.10,
                        wspace=0.02, hspace=0.02)

    im = None
    for ri, (topic, cont, lab) in enumerate(rows):
        for mi, mode in enumerate(modes):
            ax = axes[ri][mi]
            ax.set_aspect("equal")
            ax.axis("off")

            d = _g_map(mode, topic, cont)
            if d is None:
                continue

            im = ax.scatter(d["x"], d["y"], c=d["theta"], s=G_SPOT_SIZE,
                            marker="s", cmap=G_CMAP, vmin=0, vmax=1,
                            edgecolors="none")

            if ri == 0:
                ax.set_title(MODE_LABEL[mode], fontsize=9, pad=6)

        pos = axes[ri][-1].get_position()
        fig.text(pos.x1 + 0.012, (pos.y0 + pos.y1) / 2, lab,
                 va="center", ha="left", fontsize=9)

    # topic label spanning the two levels of each topic
    for ti, topic in enumerate(G_TOPICS):
        top = axes[ti * len(LEVELS)][0].get_position()
        bot = axes[ti * len(LEVELS) + len(LEVELS) - 1][0].get_position()
        fig.text(top.x0 - 0.02, (top.y1 + bot.y0) / 2,
                 f"{G_DATASET_LABEL} Topic {topic}",
                 rotation=90, va="center", ha="right", fontsize=9)

    top = axes[0][-1].get_position()
    bot = axes[-1][-1].get_position()
    fig.text(top.x1 + 0.055, (top.y1 + bot.y0) / 2, "Spike-in level",
             rotation=90, va="center", ha="left", fontsize=9)

    if im is not None:
        cax = fig.add_axes([0.16, 0.045, 0.16, 0.018])
        cb = fig.colorbar(im, cax=cax, orientation="horizontal")
        cb.set_ticks([0, 1])
        cb.ax.tick_params(labelsize=8, length=2)
        cb.set_label("Topic-spot distribution", fontsize=8, labelpad=4)

    fig.patch.set_facecolor("white")
    fig.savefig(f"{OUT}/Figure_{panel_label}.png", dpi=300, facecolor="white")
    plt.close(fig)


def load_scores():
    frames = []
    for mode, path in SCORE_PATHS.items():
        if not os.path.exists(path):
            print(f"  missing {path}")
            continue
        df = pd.read_csv(path)
        df["mode"] = mode
        frames.append(df)
    if not frames:
        raise SystemExit("no score CSVs in ../data")

    res = pd.concat(frames, ignore_index=True)
    res = res[~res["baseline"].astype(bool)]
    print(res.groupby(["mode", "dataset"])["model"].count().to_string())
    return res


def vals(col, cont=None, mod=None):
    name = col if mod is None else f"{col}_{mod}"

    def _v(d, ds):
        s = d[d["dataset"] == ds]
        if cont is not None:
            s = s[np.isclose(s["contamination"].astype(float), cont)]
        if name not in s.columns:
            return np.array([])
        v = s[name].to_numpy(float)
        return v[np.isfinite(v)]

    return _v


def by_dataset(shade):
    return lambda ds: _shade(DATASET_COLOR.get(ds, "0.5"), shade)


def by_modality(mod, shade):
    return lambda _ds: _shade(MOD_COLOR[mod], shade)


def _strip(ax, d, datasets, series, y_lim, y_step, tag):
    rng = np.random.default_rng(0)
    w = 0.8 / max(len(series), 1)

    for gi, ds in enumerate(datasets):
        for si, (label, colour_of, values_of) in enumerate(series):
            v = values_of(d, ds)
            if not len(v):
                continue
            x = gi + (si - (len(series) - 1) / 2) * w
            ax.scatter(x + rng.uniform(-w * 0.28, w * 0.28, len(v)), v,
                       s=POINT_SIZE, color=colour_of(ds), edgecolors="none",
                       zorder=2)
            ax.boxplot([v], positions=[x], widths=w * 0.62, showfliers=False,
                       patch_artist=True, zorder=5,
                       medianprops=dict(color=BOX_EDGE, lw=BOX_MEDIAN_LW),
                       whiskerprops=dict(color=BOX_EDGE, lw=BOX_LW),
                       capprops=dict(color=BOX_EDGE, lw=BOX_LW),
                       boxprops=dict(facecolor="none", edgecolor=BOX_EDGE,
                                     lw=BOX_LW))
            out = v[(v < y_lim[0]) | (v > y_lim[1])]
            if out.size:
                print(f"  {tag} / {ds} / {label}: {out.size} point(s) outside "
                      f"{y_lim} NOT DRAWN: "
                      + ", ".join(f"{x_:.3f}" for x_ in np.sort(out)))

    ax.set_xticks(np.arange(len(datasets)))
    ax.set_xlim(-0.5, len(datasets) - 0.5)
    ax.set_ylim(*y_lim)
    ax.yaxis.set_major_locator(MultipleLocator(y_step))
    ax.grid(axis="y", alpha=0.25)
    ax.set_axisbelow(True)
    ax.tick_params(labelsize=8)


def facet_figure(res, rows, series_for, panel_label, legends, figsize):
    """rows: (name, ylabel, column). series_for(column) -> strip series."""
    modes = [m for m in MODE_LABEL if m in set(res["mode"])]
    datasets = [d for d in DATASET_COLOR if d in set(res["dataset"])]

    fig, axes = plt.subplots(len(rows), len(modes), figsize=figsize,
                             sharey="row", sharex="col", squeeze=False)
    fig.subplots_adjust(left=0.13, right=0.98, top=0.86, bottom=0.26,
                        wspace=0.06, hspace=0.10)

    for ri, (name, ylab, col) in enumerate(rows):
        for mi, mode in enumerate(modes):
            ax = axes[ri][mi]
            _strip(ax, res[res["mode"] == mode], datasets, series_for(col),
                   Y_LIM[name], Y_STEP[name], f"{name}/{mode}")
            if ri == 0:
                ax.set_title(MODE_LABEL[mode], fontsize=9)
            if ri == len(rows) - 1:
                ax.set_xticklabels([DATASET_LABEL.get(d, d) for d in datasets],
                                   rotation=45, ha="right", fontsize=8)
            else:
                ax.set_xticklabels([])
            if mi == 0:
                ax.set_ylabel(ylab, fontsize=9)

    for handles, title, anchor in legends:
        leg = fig.legend(handles=handles, title=title, frameon=False,
                         fontsize=8, title_fontsize=8, ncol=len(handles),
                         loc="upper center", bbox_to_anchor=anchor)
        leg._legend_box.align = "left"

    fig.patch.set_facecolor("white")
    fig.savefig(f"{OUT}/Figure_{panel_label}.png", dpi=300, facecolor="white")
    plt.close(fig)


def make_panel_h(res, panel_label="E3h"):
    def series_for(col):
        return [(lab, by_dataset(sh), vals(col, cont))
                for cont, lab, sh in LEVELS]

    level_handles = [_dot(_shade("0.1", sh), lab) for _c, lab, sh in LEVELS]

    facet_figure(
        res,
        rows=[("ari", "ARI", "ari_outside"),
              ("auroc", "AUROC", "auroc_inside")],
        series_for=series_for,
        panel_label=panel_label,
        legends=[(level_handles, "Spike-in level", (0.5, 1.0))],
        figsize=(6.2, 5.0),
    )


def make_panel_i(res, panel_label="E3i"):
    def series_for(col):
        return [(f"{MOD_LABEL[mod]}, {lab}", by_modality(mod, sh),
                 vals(col, cont, mod))
                for mod in ("rna", "chrom")
                for cont, lab, sh in LEVELS]

    level_handles = [_dot(_shade("0.1", sh), lab) for _c, lab, sh in LEVELS]
    mod_handles = [_dot(MOD_COLOR[m], MOD_LABEL[m]) for m in ("rna", "chrom")]

    facet_figure(
        res,
        rows=[("beta_cos", "Cosine similarity\n(all features)", "beta_cos"),
              ("beta_rbo", "Rank-biased overlap\n(leading features)",
               "beta_rbo")],
        series_for=series_for,
        panel_label=panel_label,
        legends=[(level_handles, "Spike-in level", (0.42, 1.0)),
                 (mod_handles, None, (0.82, 1.0))],
        figsize=(6.6, 5.0),
    )


# ═══════════════════════════════════════════════════════════════════════════
# Figure E3a
# ═══════════════════════════════════════════════════════════════════════════
make_panel("P22H3K4me3", "E3a")

# ═══════════════════════════════════════════════════════════════════════════
# Figure E3b
# ═══════════════════════════════════════════════════════════════════════════
make_panel("P22H3K27me3", "E3b")

# ═══════════════════════════════════════════════════════════════════════════
# Figure E3c
# ═══════════════════════════════════════════════════════════════════════════
make_panel("P22H3K27ac", "E3c")

# ═══════════════════════════════════════════════════════════════════════════
# Figure E3g
# ═══════════════════════════════════════════════════════════════════════════
make_panel_g("E3g")

# ═══════════════════════════════════════════════════════════════════════════
# Figure E3h
# ═══════════════════════════════════════════════════════════════════════════
scores = load_scores()
make_panel_h(scores, "E3h")

# ═══════════════════════════════════════════════════════════════════════════
# Figure E3i
# ═══════════════════════════════════════════════════════════════════════════
make_panel_i(scores, "E3i")