from matplotlib.colors import LinearSegmentedColormap as LSC
from matplotlib.offsetbox import AnnotationBbox, HPacker, TextArea
from matplotlib.ticker import FixedLocator, FuncFormatter, NullLocator
import matplotlib.pyplot as plt
import mtopic
import numpy as np
import pandas as pd
import os
import re

INPUT_PATHS = {
    "trained_model": "../data/HumanPBMC_trained.h5mu",
    "peak_gene_3h": "../data/PBMC_Cyto_PeakGeneAssoc.csv",
    "scores_3i": "../data/PRF1_scores.csv",
    "colmeta_3i": "../data/PRF1_colmeta.csv",
    "rowmeta_3i": "../data/PRF1_rowmeta.csv",
}

plt.rcParams["font.family"] = "Arial"
plt.rcParams["pdf.fonttype"] = 42

OUT = "../figures/Figure_3"
os.makedirs(OUT, exist_ok=True)

path_trained = INPUT_PATHS["trained_model"]
mdata = mtopic.read.h5mu(path_trained)
TOPIC_COLOR = mdata.uns["TOPIC_COLOR"]
TOPIC_CELLTYPE = mdata.uns["TOPIC_CELLTYPE"]
CELLTYPE_COLOR = mdata.uns["CELLTYPE_COLOR"]


# ═══════════════════════════════════════════════════════════════════════════
# Figure 3a
# ═══════════════════════════════════════════════════════════════════════════
mtopic.pl.dominant_topics(
    mdata,
    x="umap",
    s=28,
    figsize=(8, 4),
    legend=True,
    palette=TOPIC_COLOR,
    legend_ncol=2,
    markerscale=2,
    annotation=TOPIC_CELLTYPE,
    save=f"{OUT}/Figure_3a.png",
)


# ═══════════════════════════════════════════════════════════════════════════
# Figure 3b
# ═══════════════════════════════════════════════════════════════════════════
mdata = mtopic.read.h5mu(path_trained)

CELLTYPE_ORDER = [
    "CD4 Naive-1",
    "CD4 Naive-2",
    "CD4 Naive-3",
    "CD4 Naive-4",
    "CD4 Naive-5",
    "CD4 Recently activated",
    "CD4 Activated",
    "CD4 Effector",
    "CD4 Memory",
    "CD4 Tregs double positive",
    "CD4 Treg",
    "CD4 Tfh",
    "CD4 Th1",
    "CD4 Th2 Th17",
    "CD8 Naive-1",
    "CD8 Naive-2",
    "CD8 Effector Memory",
    "CD8 Effector",
    "NK",
    "gdT",
    "MAIT",
    "Naive B",
    "Memory B",
    "pDCs",
    "cDCs",
    "Activated DCs",
    "HSCs",
]
EXCLUDE_KEYWORDS = ["Doublet"]

palette = TOPIC_COLOR
figsize = (10, 20)


def topic_sort_key(t):
    ann = TOPIC_CELLTYPE.get(t, t)
    return CELLTYPE_ORDER.index(ann) if ann in CELLTYPE_ORDER else len(CELLTYPE_ORDER)


topic_order = sorted(mdata.obsm["topics"].columns, key=topic_sort_key)
topic_order = [
    t
    for t in topic_order
    if not any(kw in TOPIC_CELLTYPE.get(t, "") for kw in EXCLUDE_KEYWORDS)
]

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

data = gamma_sorted.values.T
color_bar_v = color_bar[:, np.newaxis, :]

fig, (ax_bar, ax_heat) = plt.subplots(
    1,
    2,
    figsize=(8, 18),
    gridspec_kw={"width_ratios": [0.04, 1], "wspace": 0.01},
)

ax_bar.imshow(
    color_bar_v[::-1], aspect="auto", interpolation="none", extent=[0, 1, N, 0]
)
ax_bar.set_ylim(0, N)
ax_bar.axis("off")

im = ax_heat.imshow(
    data[::-1],
    aspect="auto",
    cmap="gnuplot",
    vmin=0,
    vmax=1,
    interpolation="none",
    extent=[0, K, N, 0],
)
ax_heat.set_ylim(0, N)
ax_heat.set_yticks([])
ax_heat.set_xticks(np.arange(K) + 0.5)
ax_heat.set_xticklabels(
    [TOPIC_CELLTYPE.get(t, t) for t in topic_order],
    fontsize=7,
    rotation=90,
)

fig.colorbar(im, ax=ax_heat, fraction=0.02, pad=0.01)

plt.savefig(f"{OUT}/Figure_3b.png", bbox_inches="tight", dpi=300)
plt.close()


# ═══════════════════════════════════════════════════════════════════════════
# Figure 3f
# ═══════════════════════════════════════════════════════════════════════════
cmap_3f = LSC.from_list(
    "blue_green",
    [
        (0.0, "#ffffff"),
        (0.05, "#ffffff"),
        (0.3, "#78c679"),
        (0.5, "#005a32"),
        (0.8, "#000000"),
        (1, "#000000"),
    ],
)

mdata = mtopic.read.h5mu(path_trained)

celltype_order = [
    "CD4 Naive-1",
    "CD4 Naive-2",
    "CD4 Naive-3",
    "CD4 Naive-4",
    "CD4 Naive-5",
    "CD4 Recently activated",
    "CD4 Activated",
    "CD4 Effector",
    "CD4 Memory",
    "CD4 Tregs double positive",
    "CD4 Treg",
    "CD4 Tfh",
    "CD4 Th1",
    "CD4 Th2 Th17",
    "CD8 Naive-1",
    "CD8 Naive-2",
    "CD8 Effector Memory",
    "CD8 Effector",
    "NK",
    "gdT",
    "MAIT",
    "Naive B",
    "Memory B",
    "pDCs",
    "cDCs",
    "Activated DCs",
    "HSCs",
]
celltype_to_topic = {v: k for k, v in TOPIC_CELLTYPE.items()}
topic_order = [
    celltype_to_topic[ct] for ct in celltype_order if ct in celltype_to_topic
]

n_top = {"atac": 500}
mod = "atac"

df_sig = mdata.mod[mod].varm["signatures"]
df = df_sig.div(df_sig.sum(axis=0), axis=1)
df = df[topic_order]

top_features = pd.Index(
    pd.concat([df.nlargest(n_top[mod], col).index.to_series() for col in df.columns])
)
df_plot = df.loc[top_features]

fig, ax = plt.subplots(figsize=(10, 5))
im = ax.imshow(df_plot.values, aspect="auto", cmap=cmap_3f, interpolation="none")

ax.set_xticks(range(df_plot.shape[1]))
ax.set_xticklabels([TOPIC_CELLTYPE.get(t, t) for t in topic_order], rotation=90)
ax.set_yticks([])
ax.set_xlabel("Topics")
ax.set_ylabel("Peaks")

cbar = fig.colorbar(im, ax=ax, label="Probability")
cbar.set_ticks([im.norm.vmin, im.norm.vmax])
plt.tight_layout()
plt.savefig(f"{OUT}/Figure_3f.png", dpi=300, bbox_inches="tight")
plt.close()


# ═══════════════════════════════════════════════════════════════════════════
# Figure 3h
# ═══════════════════════════════════════════════════════════════════════════
PT = 72.27 / 25.4
SIZE = 12
LABEL_SIZE = 16
POINT_D = 1.25 * (0.5 * PT + 0.5 * 96 / 25.4 / 2)
LINE_W = 0.5 * PT
TICK_LEN = 2.75

df_3h = pd.read_csv(INPUT_PATHS["peak_gene_3h"])

top_hhi = df_3h[
    df_3h["HHI"] == df_3h.groupby("Avg_log2FC_Type")["HHI"].transform("max")
].pipe(lambda d: d[d["GA"] == d.groupby("Avg_log2FC_Type")["GA"].transform("max")])

PANELS_3H = {
    "Avg_log2FC_NK cells": "NK",
    "Avg_log2FC_MAIT": "MAIT",
    "Avg_log2FC_CD8 Effector Memory": "CD8 Effector Memory",
    "Avg_log2FC_gdT cells": "gdT",
}
groups = list(PANELS_3H)
missing = [g for g in groups if g not in set(df_3h["Avg_log2FC_Type"])]
if missing:
    raise ValueError(
        f"{missing} not in Avg_log2FC_Type; available: "
        f"{sorted(df_3h['Avg_log2FC_Type'].dropna().unique())}"
    )


def log_breaks(lo, hi):
    lo_e, hi_e = int(np.floor(np.log10(lo))), int(np.ceil(np.log10(hi)))
    for mult in ([1], [1, 3], [1, 2, 5], list(range(1, 10))):
        br = [m * 10.0**e for e in range(lo_e, hi_e + 1) for m in mult]
        br = [b for b in br if lo <= b <= hi]
        if len(br) >= 3:
            return br
    return br


def fmt_g(x, _):
    return f"{x:g}"


def split_label(label):
    m = re.match(r"^atac:(chr\w+?)[-_:](\d+)[-_]\d+_rna:(.+)$", str(label))
    if not m:
        return None, str(label)
    chrom, start, gene = m.groups()
    return gene, f"{chrom}:{start}"


def gene_peak_box(label, fontsize):
    gene, peak = split_label(label)
    if gene is None:
        return TextArea(peak, textprops=dict(fontsize=fontsize))
    return HPacker(
        children=[
            TextArea(gene, textprops=dict(fontsize=fontsize, fontstyle="italic")),
            TextArea(f"–{peak}", textprops=dict(fontsize=fontsize)),
        ],
        pad=0,
        sep=0,
        align="baseline",
    )


FIG_W, FIG_H = 3.6, 3.9
AX_L, AX_B, AX_S = 0.8, 0.7, 2.6
px_per_pt = 300 / 72

TOTAL_W = FIG_W * len(groups)
fig = plt.figure(figsize=(TOTAL_W, FIG_H), dpi=300)

for i, grp in enumerate(groups):
    ax = fig.add_axes(
        [(i * FIG_W + AX_L) / TOTAL_W, AX_B / FIG_H, AX_S / TOTAL_W, AX_S / FIG_H]
    )

    d = df_3h[(df_3h["Avg_log2FC_Type"] == grp) & (df_3h["GA"] > 0)]

    ax.scatter(
        d["GA"], d["HHI"], s=POINT_D**2, color="black", linewidths=0, clip_on=False
    )

    ax.set_xscale("log")
    lx = np.log10(d["GA"])
    pad = 0.05 * (lx.max() - lx.min())
    xlo, xhi = 10 ** (lx.min() - pad), 10 ** (lx.max() + pad)
    ax.set_xlim(xlo, xhi)
    ax.xaxis.set_major_locator(FixedLocator(log_breaks(xlo, xhi)))
    ax.xaxis.set_minor_locator(NullLocator())
    ax.xaxis.set_major_formatter(FuncFormatter(fmt_g))

    ax.set_ylim(0, 1.05)
    ax.set_yticks([0, 0.25, 0.5, 0.75, 1])
    ax.yaxis.set_major_formatter(FuncFormatter(fmt_g))

    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_linewidth(LINE_W)
    ax.tick_params(
        labelsize=SIZE, colors="black", width=LINE_W, length=TICK_LEN, direction="out"
    )
    ax.set_title(PANELS_3H[grp], fontsize=SIZE, pad=2.2 * SIZE)
    ax.set_xlabel("Global association score (log10)", fontsize=SIZE)
    ax.set_ylabel("HHI", fontsize=SIZE)

    fig.canvas.draw()
    a_box = ax.get_window_extent()

    for _, r in top_hhi[top_hhi["Avg_log2FC_Type"] == grp].iterrows():
        ab = AnnotationBbox(
            gene_peak_box(r["label"], LABEL_SIZE),
            (r["GA"], r["HHI"]),
            xybox=(0.5, 1.0),
            boxcoords="axes fraction",
            box_alignment=(1, 0),
            frameon=False,
            pad=0,
            annotation_clip=False,
            arrowprops=dict(
                arrowstyle="-", color="black", lw=1, shrinkA=2, shrinkB=2
            ),
        )
        ax.add_artist(ab)
        fig.canvas.draw()

        w = ab.offsetbox.get_window_extent().width
        px, py = ax.transData.transform((r["GA"], r["HHI"]))
        margin = 0.3 * SIZE * px_per_pt
        right = min(px + w / 2, a_box.x1)
        if right - w < a_box.x0 + margin:
            right = a_box.x0 + margin + w
        bottom = min(py + 1.8 * SIZE * px_per_pt, a_box.y1 + 0.3 * SIZE * px_per_pt)
        ab.xybox = tuple(ax.transAxes.inverted().transform((right, bottom)))

fig.savefig(f"{OUT}/Figure_3h.png", dpi=300)
plt.close(fig)


# ═══════════════════════════════════════════════════════════════════════════
# Figure 3i
# ═══════════════════════════════════════════════════════════════════════════
AS = pd.read_csv(INPUT_PATHS["scores_3i"], index_col=0)
col_meta = pd.read_csv(INPUT_PATHS["colmeta_3i"])
row_meta = pd.read_csv(INPUT_PATHS["rowmeta_3i"])

n_rows, n_cols = AS.shape
TOP_PEAKS = 1
COLOR_PEAKS = [
    "#CD0000", "#000080", "#00FFFF", "#006400",
    "#551A8B", "#CD853F", "#642915", "#58508d",
]
if TOP_PEAKS <= len(COLOR_PEAKS):
    peak_colors = COLOR_PEAKS[:TOP_PEAKS]
else:
    peak_colors = [
        LSC.from_list("peaks", COLOR_PEAKS)(v) for v in np.linspace(0, 1, TOP_PEAKS)
    ]
peak_colors = list(peak_colors) + ["black"] * (n_rows - len(peak_colors))

cmap_3i = LSC.from_list(
    "atac",
    ["#F2F2F2", "#A1D99B", "#74C476", "#41AB5D", "#238B45", "#005A32", "#000000"],
    N=100,
)

R_COLORS = {
    "magenta": "#FF00FF",
    "plum1": "#FFBBFF",
    "purple": "#A020F0",
    "purple4": "#551A8B",
}


def r_color(c):
    c = str(c).strip()
    if c.lower() in R_COLORS:
        return R_COLORS[c.lower()]
    try:
        return plt.matplotlib.colors.to_hex(c)
    except ValueError:
        raise ValueError(f"Unknown R colour {c!r}; add its hex code to R_COLORS")


TEXT_3I = 10
BORDER_3I = 0.8
PATH_C_3I = "#999999"
PATH_W_3I = 0.8
DOT_3I = 32
TOP_DOT_3I = 60
HHI_TICKS = [0, 0.75]
GA_TICKS = [0.01, 0.05]

HM_L, HM_B = 1.5, 2.3
HM_W, HM_H = 4.57, 3.57
TOP_GAP, TOP_H = 0.25, 0.71
R_GAP, R_W = 0.2, 0.68
CB_W, CB_H, CB_GAP = 1.5, 0.14, 0.32

FIG_W_3I = HM_L + HM_W + 2 * (R_GAP + R_W) + 0.8
FIG_H_3I = HM_B + HM_H + TOP_GAP + TOP_H + 0.4


def axes_in(fig, left, bottom, width, height):
    return fig.add_axes(
        [left / FIG_W_3I, bottom / FIG_H_3I, width / FIG_W_3I, height / FIG_H_3I]
    )


def style_box(ax):
    for sp in ax.spines.values():
        sp.set_visible(True)
        sp.set_linewidth(BORDER_3I)
        sp.set_color("black")
    ax.grid(False)
    ax.tick_params(
        labelsize=TEXT_3I, colors="black", width=BORDER_3I, length=TICK_LEN
    )


def row_positions(meta):
    if "Rows" in meta and meta["Rows"].astype(str).isin(AS.index.astype(str)).all():
        lookup = {str(r): i for i, r in enumerate(AS.index)}
        return meta["Rows"].astype(str).map(lookup).to_numpy()
    if "Rows" in meta and pd.api.types.is_numeric_dtype(meta["Rows"]):
        return meta["Rows"].to_numpy() - 1
    return np.arange(len(meta))


fig = plt.figure(figsize=(FIG_W_3I, FIG_H_3I))

ax_hm = axes_in(fig, HM_L, HM_B, HM_W, HM_H)
vmin, vmax = AS.values.min(), AS.values.max()
im = ax_hm.imshow(
    AS.values,
    cmap=cmap_3i,
    vmin=vmin,
    vmax=vmax,
    aspect="auto",
    interpolation="nearest",
    origin="lower",
    extent=[-0.5, n_cols - 0.5, -0.5, n_rows - 0.5],
)
style_box(ax_hm)
ax_hm.set_xticks(range(n_cols))
ax_hm.set_xticklabels(AS.columns, rotation=90, ha="center", va="top")
ax_hm.set_yticks(range(n_rows))
ax_hm.set_yticklabels(AS.index)
for lab, col in zip(ax_hm.get_yticklabels(), peak_colors[::-1]):
    lab.set_color(col)

top_y = HM_B + HM_H + TOP_GAP
ax_top = axes_in(fig, HM_L, top_y, HM_W, TOP_H)
ax_top.sharex(ax_hm)
col_x = col_meta["Rank"].to_numpy() - 1
score = col_meta["Score"].to_numpy()
ax_top.plot(col_x, score, color=PATH_C_3I, lw=PATH_W_3I, zorder=1)
ax_top.scatter(
    col_x,
    score,
    s=TOP_DOT_3I,
    c=[r_color(c) for c in col_meta["color"]],
    linewidths=0,
    zorder=2,
)
ax_top.set_yscale("log")
ly = np.log10(score)
pad_y = 0.18 * (ly.max() - ly.min())
ax_top.set_ylim(10 ** (ly.min() - pad_y), 10 ** (ly.max() + pad_y))
ax_top.yaxis.set_major_locator(FixedLocator([score.min(), score.max()]))
ax_top.yaxis.set_minor_locator(NullLocator())
ax_top.set_yticklabels(["Min", "Max"])
ax_top.set_ylabel("RNA score\n(log10)", fontsize=TEXT_3I)
style_box(ax_top)
ax_top.tick_params(labelbottom=False, bottom=False)

row_y = row_positions(row_meta)
for k, (col, title, ticks) in enumerate(
    [
        ("HHI", "HHI", HHI_TICKS),
        ("ModelScores", "Global\nassociation\nscore", GA_TICKS),
    ]
):
    ax_r = axes_in(fig, HM_L + HM_W + R_GAP + k * (R_W + R_GAP), HM_B, R_W, HM_H)
    ax_r.sharey(ax_hm)
    vals = row_meta[col].to_numpy()
    ax_r.plot(vals, row_y, color=PATH_C_3I, lw=PATH_W_3I, zorder=1)
    ax_r.scatter(vals, row_y, s=DOT_3I, color="black", linewidths=0, zorder=2)
    lo, hi = min(vals.min(), min(ticks)), max(vals.max(), max(ticks))
    span = hi - lo
    ax_r.set_xlim(lo - 0.15 * span, hi + 0.15 * span)
    ax_r.xaxis.set_major_locator(FixedLocator(ticks))
    ax_r.xaxis.set_major_formatter(FuncFormatter(fmt_g))
    ax_r.set_xlabel(title, fontsize=TEXT_3I)
    style_box(ax_r)
    ax_r.tick_params(labelleft=False, left=False)

fig.canvas.draw()
labels_bottom = fig.dpi_scale_trans.inverted().transform(
    (0, ax_hm.xaxis.get_tightbbox().y0)
)[1]
ax_cb = axes_in(
    fig, HM_L + (HM_W - CB_W) / 2, labels_bottom - CB_GAP - CB_H, CB_W, CB_H
)
cb = fig.colorbar(im, cax=ax_cb, orientation="horizontal", ticks=[vmin, vmax])
cb.ax.set_xticklabels(["Min", "Max"], fontsize=TEXT_3I)
cb.outline.set_linewidth(BORDER_3I)
cb.ax.tick_params(length=0, pad=2)
cb.ax.set_title("ATAC score", fontsize=TEXT_3I, pad=4)

fig.savefig(f"{OUT}/Figure_3i.png", dpi=300, bbox_inches="tight")
plt.close(fig)