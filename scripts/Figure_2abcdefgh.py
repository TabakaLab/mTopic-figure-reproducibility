from collections import defaultdict
import re

import matplotlib.colors as mcolors
from matplotlib.collections import PatchCollection
from matplotlib.colors import LinearSegmentedColormap
import matplotlib.gridspec as gridspec
from matplotlib.patches import ConnectionPatch, Rectangle, Wedge
import matplotlib.pyplot as plt
import mtopic
import numpy as np
import pandas as pd
import scanpy as sc
from scipy.sparse import issparse
from scipy.stats import gmean
import os

INPUT_PATHS = {
    "trained_model": "../data/HumanTonsil_trained.h5mu",
}

plt.rcParams["font.family"] = "Arial"
plt.rcParams["pdf.fonttype"] = 42

OUT = "../figures/Figure_2"
os.makedirs(OUT, exist_ok=True)

CMAP = LinearSegmentedColormap.from_list(
    "topic_cmap",
    ["#000000", "#6600ff", "#cc0066", "#ff4400", "#ffaa00", "#ffff00"],
)

path_trained = INPUT_PATHS["trained_model"]
mdata = mtopic.read.h5mu(path_trained)
TOPIC_COLOR = mdata.uns["TOPIC_COLOR"]
TOPIC_CELLTYPE = mdata.uns["TOPIC_CELLTYPE"]
CELLTYPE_COLOR = mdata.uns["CELLTYPE_COLOR"]

FONT = 10
TITLE_FONT = 11
BG_KWARGS = dict(s=2, c="lightgrey", linewidths=0, rasterized=True)


def pretty(name):
    name = name.replace("Myeloid cells", "Myeloid")
    return re.sub(r"-(\d+)$", r" \1", name)


def spot_map(ax, coords, values, cmap, vmin, vmax, title, italic=False):
    ax.scatter(coords[:, 0], coords[:, 1], **BG_KWARGS, zorder=1)
    sc_map = ax.scatter(
        coords[:, 0],
        coords[:, 1],
        c=values,
        s=12,
        marker="h",
        cmap=cmap,
        linewidths=0,
        rasterized=True,
        zorder=2,
        vmin=vmin,
        vmax=vmax,
    )
    ax.set_title(
        title, fontsize=TITLE_FONT, pad=4, fontstyle="italic" if italic else "normal"
    )
    ax.axis("off")
    ax.set_aspect("equal")
    return sc_map


def unit_colorbar(fig, mappable, axes, label, ticklabels=("0", "1"), **kwargs):
    cbar = fig.colorbar(mappable, ax=axes, ticks=[0, 1], **kwargs)
    cbar.ax.set_yticklabels(list(ticklabels), fontsize=FONT)
    cbar.set_label(label, fontsize=FONT)
    return cbar


def normalised_coords(md):
    xy = md.obsm["coords"].copy()
    xy = xy - xy.min()
    xy = xy / xy.max().max()
    return np.asarray(xy)


def draw_pies(ax, topics_df, xy, palette, radius, xrange=(0, 1), yrange=(0, 1)):
    mask = (
        (xy[:, 0] >= xrange[0])
        & (xy[:, 0] <= xrange[1])
        & (xy[:, 1] >= yrange[0])
        & (xy[:, 1] <= yrange[1])
    )
    props = topics_df.values[mask]
    pts = xy[mask]
    cumsum = np.cumsum(props, axis=1) * 360
    starts = np.hstack([np.zeros((len(props), 1)), cumsum[:, :-1]])
    for k, topic in enumerate(topics_df.columns):
        patches = [
            Wedge(tuple(pts[i]), radius, 90 - cumsum[i, k], 90 - starts[i, k])
            for i in range(len(props))
        ]
        ax.add_collection(
            PatchCollection(
                patches, facecolor=palette[topic], edgecolor="none", rasterized=True
            )
        )
    ax.set(
        aspect="equal",
        xlim=[xrange[0] - 0.01, xrange[1] + 0.01],
        ylim=[yrange[0] - 0.01, yrange[1] + 0.01],
    )


def pie_with_zoom(
    md, palette, radius, xrange, yrange, title=None, legend_ax=None, frame=True
):
    topics_df = md.obsm["topics"]
    xy = normalised_coords(md)

    fig = plt.figure(figsize=(11, 6.5 if legend_ax else 5))
    gs = gridspec.GridSpec(
        2 if legend_ax else 1,
        2,
        figure=fig,
        height_ratios=[1, 0.42] if legend_ax else [1],
        wspace=0.25,
        hspace=0.05,
    )
    ax_full = fig.add_subplot(gs[0, 0])
    ax_zoom = fig.add_subplot(gs[0, 1])

    draw_pies(ax_full, topics_df, xy, palette, radius)
    ax_full.axis("off")
    if title:
        ax_full.set_title(title, fontsize=TITLE_FONT)

    draw_pies(ax_zoom, topics_df, xy, palette, radius, xrange, yrange)
    if not frame:
        ax_zoom.axis("off")
        if legend_ax:
            legend_ax(fig.add_subplot(gs[1, :]))
        return fig

    ax_zoom.set_xticks([])
    ax_zoom.set_yticks([])
    for sp in ax_zoom.spines.values():
        sp.set_linewidth(1.2)

    ax_full.add_patch(
        Rectangle(
            (xrange[0], yrange[0]),
            xrange[1] - xrange[0],
            yrange[1] - yrange[0],
            fill=False,
            lw=1.2,
            zorder=5,
        )
    )
    start = (xrange[1], (yrange[0] + yrange[1]) / 2)
    fig.add_artist(
        ConnectionPatch(
            xyA=start,
            coordsA=ax_full.transData,
            xyB=(0, 0.5),
            coordsB=ax_zoom.transAxes,
            lw=1.2,
            color="black",
        )
    )
    ax_full.plot(*start, "o", color="black", ms=4, zorder=6)
    ax_zoom.plot(0, 0.5, "o", color="black", ms=4, transform=ax_zoom.transAxes, clip_on=False)

    if legend_ax:
        legend_ax(fig.add_subplot(gs[1, :]))
    return fig


# ═══════════════════════════════════════════════════════════════════════════
# Figure 2a
# ═══════════════════════════════════════════════════════════════════════════
GROUP_ORDER_2a = [
    "Light zone",
    "Dark zone",
    "Mantle zone",
    "Plasma cells",
    "T cell zone",
    "Myeloid",
    "MRC",
    "Epithelial apical",
    "Epithelial transitional",
    "Epithelial basal",
]

groups_2a = defaultdict(list)
for topic, celltype in TOPIC_CELLTYPE.items():
    if celltype == "Unknown" or topic not in TOPIC_COLOR:
        continue
    name = pretty(celltype)
    m = re.match(r"^(.*?)(?: (\d+))?$", name)
    groups_2a[m.group(1)].append((int(m.group(2) or 0), TOPIC_COLOR[topic]))

legend_rows = []
for base in sorted(
    groups_2a,
    key=lambda b: GROUP_ORDER_2a.index(b) if b in GROUP_ORDER_2a else len(GROUP_ORDER_2a),
):
    members = sorted(groups_2a[base])
    nums = [n for n, _ in members if n]
    label = f"{base} {min(nums)}–{max(nums)}" if len(nums) > 1 else base
    legend_rows.append((label, [c for _, c in members]))


def topic_legend(ax):
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")
    ax.text(0.5, 0.88, "Topics", ha="center", va="center", fontsize=FONT, weight="bold")
    n_col = int(np.ceil(len(legend_rows) / 2))
    dy = 0.6 / max(n_col - 1, 1)
    for j, (label, colors) in enumerate(legend_rows):
        col, row = divmod(j, n_col)
        tx = 0.33 + 0.3 * col
        ty = 0.72 - row * dy
        ax.text(tx, ty, label, ha="left", va="center", fontsize=FONT)
        for k, c in enumerate(reversed(colors)):
            ax.scatter(tx - 0.012 - k * 0.014, ty, s=45, color=c, linewidths=0)


fig = pie_with_zoom(
    mdata,
    TOPIC_COLOR,
    radius=0.0074,
    xrange=[0.37, 0.52],
    yrange=[0.45, 0.63],
    legend_ax=topic_legend,
    frame=False,
)
fig.savefig(f"{OUT}/Figure_2a.png", bbox_inches="tight", dpi=300)
plt.close(fig)


# ═══════════════════════════════════════════════════════════════════════════
# Figure 2b
# ═══════════════════════════════════════════════════════════════════════════
mdata = mtopic.read.h5mu(path_trained)

ZONE_ORDER = [
    "Light zone",
    "Dark zone",
    "Mantle zone",
    "Plasma cells",
    "T cell zone",
    "Myeloid",
    "MRC",
    "Epithelial",
]
ZONE_TITLES = {"Epithelial": "Epithelial cells"}

topics_df = mdata.obsm["topics"]
coords = np.array(mdata.obsm["coords"])

base_to_topics = defaultdict(list)
tcell_to_topics = defaultdict(list)
for topic, celltype in TOPIC_CELLTYPE.items():
    if celltype == "Unknown" or topic not in topics_df.columns:
        continue
    base = next((z for z in ZONE_ORDER if celltype.startswith(z)), celltype)
    base_to_topics[base].append(topic)
    if base == "T cell zone":
        tcell_to_topics[celltype].append(topic)

zones = [
    (ZONE_TITLES.get(z, z), topics_df[base_to_topics[z]].sum(axis=1).values)
    for z in ZONE_ORDER
    if z in base_to_topics
]
tcell_order = ["T cell zone-1", "T cell zone-2", "T cell zone-3", "T cell zone-4"]
tcells = [
    (pretty(ct), topics_df[tcell_to_topics[ct]].sum(axis=1).values)
    for ct in tcell_order
    if ct in tcell_to_topics
]

prot_adata = mdata.mod["prot"]
protein_cmap = mcolors.LinearSegmentedColormap.from_list(
    "protein_level",
    [
        (0.0, "#E4E4E4"),
        (0.5, "#E4E4E4"),
        (0.7, "#CDAC3B"),
        (0.8, "#967402"),
        (1.0, "#000000"),
    ],
)
proteins_to_plot = ["prot:CD8A", "prot:PTPRC-1", "prot:CCR7", "prot:VIM"]
protein_labels = ["CD8A", "CD45RO", "CCR7", "VIM"]

raw = prot_adata.layers["counts"]
raw = (raw.toarray() if issparse(raw) else np.asarray(raw)).astype(float)
gm = gmean(raw + 1, axis=1, keepdims=True)
clr = np.log((raw + 1) / gm)

proteins = []
for protein, label in zip(proteins_to_plot, protein_labels):
    expr = clr[:, list(prot_adata.var_names).index(protein)]
    lo, hi = expr.min(), expr.max()
    proteins.append((label, (expr - lo) / (hi - lo if hi > lo else 1)))

PANEL_W, GAP_X, TITLE_H, BLOCK_GAP = 2.6, 0.15, 0.3, 0.25
CB_W, CB_GAP, MARGIN = 0.12, 0.2, 0.1
span = np.ptp(coords, axis=0)
panel_h = PANEL_W * span[1] / span[0]
row_h = panel_h + TITLE_H
blocks = [
    (zones, CMAP, "Cumulative\ntopic-spot distribution", ("0", "1")),
    (tcells, CMAP, "Topic-spot distribution", ("0", "1")),
    (proteins, protein_cmap, "Protein level", ("0", "Max")),
]
n_rows_b = [int(np.ceil(len(items) / 4)) for items, *_ in blocks]
fig_w = MARGIN + 4 * PANEL_W + 3 * GAP_X + CB_GAP + CB_W + 1.0
fig_h = MARGIN * 2 + sum(n_rows_b) * row_h + BLOCK_GAP * (len(blocks) - 1)
fig = plt.figure(figsize=(fig_w, fig_h))


def add_axes_in(left, bottom, width, height):
    return fig.add_axes([left / fig_w, bottom / fig_h, width / fig_w, height / fig_h])


top = fig_h - MARGIN
for (items, cmap, label, ticklabels), nr in zip(blocks, n_rows_b):
    block_top = top
    for j, (title, v) in enumerate(items):
        r, c = divmod(j, 4)
        ax = add_axes_in(
            MARGIN + c * (PANEL_W + GAP_X),
            top - (r + 1) * row_h,
            PANEL_W,
            panel_h,
        )
        mappable = spot_map(ax, coords, v, cmap, 0, 1, title)
    block_h = nr * row_h
    cb_h = min(0.6 * block_h, 2.0)
    cax = add_axes_in(
        MARGIN + 4 * PANEL_W + 3 * GAP_X + CB_GAP,
        block_top - block_h / 2 - cb_h / 2 - TITLE_H / 2,
        CB_W,
        cb_h,
    )
    cbar = fig.colorbar(mappable, cax=cax, ticks=[0, 1])
    cbar.ax.set_yticklabels(list(ticklabels), fontsize=FONT)
    cbar.set_label(label, fontsize=FONT)
    top -= block_h + BLOCK_GAP

fig.savefig(f"{OUT}/Figure_2b.png", bbox_inches="tight", dpi=300)
plt.close(fig)


# ═══════════════════════════════════════════════════════════════════════════
# Figure 2c
# ═══════════════════════════════════════════════════════════════════════════
mdata = mtopic.read.h5mu(path_trained)

ZONE_ORDER_2c = [
    "Light zone",
    "Dark zone",
    "Mantle zone",
    "Plasma cells",
    "T cell zone",
    "Myeloid",
    "MRC",
    "Epithelial apical",
    "Epithelial transitional",
    "Epithelial basal",
]
GROUP_LABELS_2c = {"Epithelial transitional": "Epithelial trans."}


def topic_sort_key(t):
    ann = TOPIC_CELLTYPE.get(t, t)
    for i, zone in enumerate(ZONE_ORDER_2c):
        if ann.startswith(zone):
            suffix = ann.replace(zone, "").strip("-").strip()
            num = int(suffix) if suffix.isdigit() else 0
            return (i, num)
    return (len(ZONE_ORDER_2c), 0)


def zone_of(ann):
    return next((z for z in ZONE_ORDER_2c if ann.startswith(z)), ann)


topic_order = sorted(mdata.obsm["topics"].columns, key=topic_sort_key)
topic_order = [t for t in topic_order if TOPIC_CELLTYPE.get(t) != "Unknown"]

gamma = mdata.obsm["topics"][topic_order]
dominant_topic = gamma.idxmax(axis=1)
dominant_celltype = dominant_topic.map(TOPIC_CELLTYPE)
dominant_proportion = gamma.max(axis=1)

celltype_order = {
    ann: i for i, ann in enumerate(dict.fromkeys(TOPIC_CELLTYPE[t] for t in topic_order))
}
cell_order = (
    pd.DataFrame(
        {
            "celltype": dominant_celltype,
            "celltype_order": dominant_celltype.map(celltype_order),
            "topic": dominant_topic,
            "proportion": dominant_proportion,
        }
    )
    .sort_values(
        ["celltype_order", "topic", "proportion"], ascending=[True, True, False]
    )
    .index
)
gamma_sorted = gamma.loc[cell_order].T.loc[topic_order]

N = gamma_sorted.shape[1]
K = gamma_sorted.shape[0]

sorted_celltypes = dominant_celltype.loc[cell_order].values
color_bar = np.array(
    [mcolors.to_rgb(CELLTYPE_COLOR.get(ct, "#bebebe")) for ct in sorted_celltypes]
)

fig = plt.figure(figsize=(12, 4))
gs = gridspec.GridSpec(
    2,
    2,
    figure=fig,
    height_ratios=[0.04, 1],
    width_ratios=[1, 0.015],
    hspace=0.02,
    wspace=0.03,
)
ax_bar = fig.add_subplot(gs[0, 0])
ax_heat = fig.add_subplot(gs[1, 0])
ax_cbar = fig.add_subplot(gs[1, 1])

ax_bar.imshow(
    color_bar[np.newaxis, :, :],
    aspect="auto",
    interpolation="none",
    extent=[0, N, 0, 1],
)
ax_bar.set_xlim(0, N)
ax_bar.axis("off")

zones_sorted = np.array([zone_of(ct) for ct in sorted_celltypes])
change = np.flatnonzero(np.r_[True, zones_sorted[1:] != zones_sorted[:-1], True])
for s, e in zip(change[:-1], change[1:]):
    z = zones_sorted[s]
    ax_bar.text(
        (s + e) / 2,
        1.3,
        GROUP_LABELS_2c.get(z, z),
        rotation=45,
        ha="left",
        va="bottom",
        fontsize=FONT,
        transform=ax_bar.get_xaxis_transform(),
    )

im = ax_heat.imshow(
    gamma_sorted.values,
    aspect="auto",
    cmap=CMAP,
    vmin=0,
    vmax=1,
    interpolation="none",
    extent=[0, N, K, 0],
)
ax_heat.set_xlim(0, N)
ax_heat.set_xticks([])
ax_heat.set_yticks([])
ax_heat.set_xlabel("Spots", fontsize=FONT)
ax_heat.set_ylabel("Topics", fontsize=FONT)

cbar = fig.colorbar(im, cax=ax_cbar, ticks=[0, 1])
cbar.ax.tick_params(labelsize=FONT)
cbar.set_label("Topic-spot distribution", fontsize=FONT)

fig.savefig(f"{OUT}/Figure_2c.png", bbox_inches="tight", dpi=300)
plt.close(fig)


# ═══════════════════════════════════════════════════════════════════════════
# Figure 2d
# ═══════════════════════════════════════════════════════════════════════════
mdata = mtopic.read.h5mu(path_trained)
topics_df = mdata.obsm["topics"]
coords = np.array(mdata.obsm["coords"])

only_celltypes = [
    "Epithelial apical-1",
    "Epithelial apical-2",
    "Epithelial transitional-1",
    "Epithelial transitional-2",
    "Epithelial basal-1",
    "Epithelial basal-2",
    "Myeloid cells-1",
    "Myeloid cells-2",
    "Myeloid cells-3",
]

celltype_to_topics = defaultdict(list)
for topic, celltype in TOPIC_CELLTYPE.items():
    if celltype in only_celltypes and topic in topics_df.columns:
        celltype_to_topics[celltype].append(topic)

ordered_items = [
    (pretty(ct), topics_df[celltype_to_topics[ct]].sum(axis=1).values)
    for ct in only_celltypes
    if ct in celltype_to_topics
]

n = len(ordered_items)
ncols = 3
nrows = int(np.ceil(n / ncols))
fig, axes = plt.subplots(
    nrows, ncols, figsize=(3 * ncols, 3.2 * nrows), constrained_layout=True
)
axes_flat = np.atleast_1d(axes).flatten()
for ax, (title, weights) in zip(axes_flat, ordered_items):
    sc_d = spot_map(ax, coords, weights, CMAP, 0, 1, title)
for ax in axes_flat[n:]:
    ax.axis("off")
unit_colorbar(
    fig, sc_d, axes_flat.tolist(), "Topic-spot distribution", shrink=0.5, pad=0.02
)
fig.savefig(f"{OUT}/Figure_2d.png", bbox_inches="tight", dpi=300)
plt.close(fig)


# ═══════════════════════════════════════════════════════════════════════════
# Figure 2e
# ═══════════════════════════════════════════════════════════════════════════
GENES = ["KRT6B", "KRT4", "KRT18", "KRT19", "KRT15", "EGFR", "IL1B", "S100A2", "CXCL14"]
modality = "rna"

CMAP_2e = LinearSegmentedColormap.from_list(
    "cmap",
    [
        (0.00, "#eaeaea"),
        (0.30, "#eaeaea"),
        (0.60, "#DE2B25"),
        (0.90, "#99000D"),
        (1.00, "#000000"),
    ],
    N=256,
)

mdata = mtopic.read.h5mu(path_trained)
coords = np.array(mdata.obsm["coords"])
rna = mdata.mod[modality].copy()
rna.X = rna.layers["counts"].copy()
sc.pp.normalize_total(rna, target_sum=1e4)
sc.pp.log1p(rna)

strip_prefix = any(v.startswith("rna:") for v in rna.var_names)

panels = []
for g in GENES:
    var_name = f"rna:{g}" if strip_prefix else g
    if var_name not in rna.var_names:
        print(f"skipping {g}: not in var_names")
        continue
    x = rna[:, var_name].X
    vals = np.asarray(x.toarray() if issparse(x) else x).flatten()
    panels.append((g, vals))

n = len(panels)
ncols = 3
nrows = int(np.ceil(n / ncols))
fig, axes = plt.subplots(
    nrows, ncols, figsize=(3.8 * ncols, 3.2 * nrows), constrained_layout=True
)
axes_flat = np.atleast_1d(axes).flatten()
for ax, (g, vals) in zip(axes_flat, panels):
    sc_plot = spot_map(
        ax, coords, vals, CMAP_2e, 0, np.quantile(vals, 0.99), g, italic=True
    )
    cbar = fig.colorbar(sc_plot, ax=ax, fraction=0.04, pad=0.02, shrink=0.6)
    cbar.outline.set_visible(False)
    cbar.ax.tick_params(labelsize=FONT - 1)
for ax in axes_flat[n:]:
    ax.axis("off")
fig.savefig(f"{OUT}/Figure_2e.png", bbox_inches="tight", dpi=300)
plt.close(fig)


# ═══════════════════════════════════════════════════════════════════════════
# Figure 2f
# ═══════════════════════════════════════════════════════════════════════════
modality = "rna"
vmin, vmax = -2.0, 2.0
cmap_2f = LinearSegmentedColormap.from_list(
    "custom_bwr",
    ["#08306B", "#1562A9", "#FFFFFF", "#DE2B25", "#99000D"],
    N=256,
)

CELL_TYPE_ORDER = [
    "Epithelial apical-1",
    "Epithelial apical-2",
    "Epithelial transitional-1",
    "Epithelial transitional-2",
    "Epithelial basal-1",
    "Epithelial basal-2",
    "Myeloid cells-1",
    "Myeloid cells-2",
    "Myeloid cells-3",
]
FEATURE_GROUPS = [
    (
        "Epithelial apical",
        [
            "rna:KRT6B",
            "rna:CRYAB",
            "rna:IL36G",
            "rna:IL22RA1",
            "rna:KRT23",
            "rna:KRT4",
            "rna:CEACAM6",
            "rna:AREG",
            "rna:IL1A",
            "rna:IL23A",
        ],
    ),
    (
        "Epithelial transitional",
        [
            "rna:DSC2",
            "rna:CLDN10",
            "rna:FXYD3",
            "rna:PTGS1",
            "rna:KRT7",
            "rna:KRT8",
            "rna:KRT18",
            "rna:KRT19",
            "rna:CD274",
            "rna:EPCAM",
            "rna:CCL22",
        ],
    ),
    (
        "Epithelial basal",
        [
            "rna:CLEC5A",
            "rna:KRT15",
            "rna:NGFR",
            "rna:FGF1",
            "rna:EGFR",
        ],
    ),
    (
        "Myeloid\nepithelial region",
        [
            "rna:CXCL8",
            "rna:IL1B",
            "rna:CXCR1",
            "rna:PTGS2",
            "rna:VCAN",
            "rna:IL1R2",
            "rna:CXCL10",
            "rna:S100A2",
            "rna:LAMP3",
            "rna:CD14",
        ],
    ),
    (
        "Myeloid\nstromal region",
        [
            "rna:CXCL14",
            "rna:PTGDS",
            "rna:SELENOP",
            "rna:APOE",
        ],
    ),
]

mdata = mtopic.read.h5mu(path_trained)

sig = mdata.mod[modality].varm["signatures"]
var_names = mdata.mod[modality].var_names
if isinstance(sig, pd.DataFrame):
    signatures = sig.copy()
    if not signatures.index.equals(var_names):
        signatures.index = var_names
else:
    arr = np.asarray(sig)
    signatures = pd.DataFrame(
        arr,
        index=var_names,
        columns=[f"topic_{i + 1}" for i in range(arr.shape[1])],
    )
signatures = signatures.rename(columns=TOPIC_CELLTYPE)

ordered_features = []
group_lengths = []
group_labels = []
for label, genes in FEATURE_GROUPS:
    present = [g for g in genes if g in signatures.index]
    ordered_features.extend(present)
    group_lengths.append(len(present))
    group_labels.append(label)

missing_cell_types = [c for c in CELL_TYPE_ORDER if c not in signatures.columns]
if missing_cell_types:
    print(f"warning: cell types missing from signatures: {missing_cell_types}")

mat = signatures.loc[
    ordered_features,
    [c for c in CELL_TYPE_ORDER if c in signatures.columns],
]

mat_t = mat.T
mu_ = mat_t.mean(axis=0)
sd = mat_t.std(axis=0, ddof=1).replace(0, np.nan)
mat_scaled = (mat_t - mu_) / sd

n_rows, n_cols = mat_scaled.shape
fig, ax = plt.subplots(figsize=(0.32 * n_cols + 3.5, 0.24 * n_rows + 1.4))

im = ax.imshow(
    mat_scaled.values,
    aspect="auto",
    cmap=cmap_2f,
    vmin=vmin,
    vmax=vmax,
    interpolation="nearest",
)

ax.set_yticks(np.arange(n_rows))
ax.set_yticklabels([pretty(c) for c in mat_scaled.index], fontsize=FONT)
ax.yaxis.tick_right()

ax.set_xticks(np.arange(n_cols))
ax.set_xticklabels(
    [g.replace("rna:", "") for g in mat_scaled.columns],
    rotation=90,
    fontsize=FONT,
    fontstyle="italic",
)

boundary = 0
for L in group_lengths[:-1]:
    boundary += L
    if L > 0:
        ax.axvline(boundary - 0.5, color="black", lw=1.0)

starts = np.cumsum([0] + group_lengths[:-1])
centers = starts + np.array(group_lengths) / 2.0 - 0.5
trans = ax.get_xaxis_transform()
for c, lab, L in zip(centers, group_labels, group_lengths):
    if L == 0:
        continue
    ax.text(
        c,
        1.02,
        lab,
        transform=trans,
        ha="center",
        va="bottom",
        fontsize=FONT,
        clip_on=False,
    )

for spine in ax.spines.values():
    spine.set_linewidth(0.8)

fig.canvas.draw()
inv = fig.transFigure.inverted()
labels_right = inv.transform((ax.yaxis.get_tightbbox().x1, 0))[0]
ax_pos = ax.get_position()
cax = fig.add_axes(
    [labels_right + 0.03, ax_pos.y0 + 0.2 * ax_pos.height, 0.012, 0.6 * ax_pos.height]
)
cbar = fig.colorbar(im, cax=cax, ticks=[-2, 0, 2])
cbar.ax.tick_params(labelsize=FONT)
cbar.set_label("Relative RNA score", fontsize=FONT)

fig.savefig(f"{OUT}/Figure_2f.png", bbox_inches="tight", dpi=300)
plt.close(fig)


# ═══════════════════════════════════════════════════════════════════════════
# Figure 2g
# ═══════════════════════════════════════════════════════════════════════════
mdata = mtopic.read.h5mu(path_trained)

palette = {
    topic: (
        TOPIC_COLOR[topic]
        if TOPIC_CELLTYPE[topic].startswith(
            ("Epithelial", "Myeloid cells-1", "Myeloid cells-2")
        )
        else "#bebebe"
    )
    for topic in TOPIC_COLOR
}

fig = pie_with_zoom(
    mdata,
    palette,
    radius=0.0074,
    xrange=[0.4, 0.7],
    yrange=[0.29, 0.6],
    title="Epithelial region",
)
fig.savefig(f"{OUT}/Figure_2g.png", bbox_inches="tight", dpi=300)
plt.close(fig)


# ═══════════════════════════════════════════════════════════════════════════
# Figure 2h
# ═══════════════════════════════════════════════════════════════════════════
n_subsample = 100
threshold = 0.5
random_seed = 17

TYPES_IN_EPI = [
    "Epithelial apical-1",
    "Epithelial apical-2",
    "Epithelial transitional-1",
    "Epithelial transitional-2",
    "Epithelial basal-1",
    "Epithelial basal-2",
    "Myeloid cells-1",
    "Myeloid cells-2",
]
CELLTYPE_ORDER_2h = TYPES_IN_EPI.copy()
APICAL1 = "Epithelial apical-1"
MYELOID2 = "Myeloid cells-2"

mdata = mtopic.read.h5mu(path_trained)
gamma = mdata.obsm["topics"]
if not isinstance(gamma, pd.DataFrame):
    n_topics = np.asarray(gamma).shape[1]
    gamma = pd.DataFrame(
        np.asarray(gamma),
        index=mdata.obs_names,
        columns=[f"topic_{i + 1}" for i in range(n_topics)],
    )
gamma = gamma.rename(columns=TOPIC_CELLTYPE)

missing_cell_types = [c for c in TYPES_IN_EPI if c not in gamma.columns]
if missing_cell_types:
    raise KeyError(f"Cell types missing from gamma columns: {missing_cell_types}")

missing_colors = [c for c in TYPES_IN_EPI if c not in CELLTYPE_COLOR]
if missing_colors:
    raise KeyError(f"Missing colors in CELLTYPE_COLOR: {missing_colors}")

COLORS = {c: CELLTYPE_COLOR[c] for c in TYPES_IN_EPI}

region_gamma = gamma.loc[:, TYPES_IN_EPI].copy()
region_gamma = region_gamma.sort_values(MYELOID2, ascending=True, kind="mergesort")
region_gamma = region_gamma.sort_values(APICAL1, ascending=False, kind="mergesort")
region_spots = region_gamma.index[region_gamma.sum(axis=1) > threshold]

if len(region_spots) == 0:
    raise ValueError("No spots cleared the threshold; try lowering it.")

gamma_epi_wide = gamma.loc[region_spots, TYPES_IN_EPI].copy()
dominant = gamma_epi_wide.idxmax(axis=1)

celltype_rank = {cell_type: i for i, cell_type in enumerate(CELLTYPE_ORDER_2h)}
dominant_rank = dominant.map(celltype_rank)
if dominant_rank.isna().any():
    bad = dominant[dominant_rank.isna()].unique().tolist()
    raise ValueError(f"Dominant cell types missing from CELLTYPE_ORDER: {bad}")

spot_order = dominant_rank.sort_values(kind="mergesort").index

rng = np.random.default_rng(random_seed)
sample_pool = np.array(region_spots)
if len(sample_pool) > n_subsample:
    spots_keep = set(rng.choice(sample_pool, size=n_subsample, replace=False))
else:
    spots_keep = set(sample_pool)
final_spot_order = [spot for spot in spot_order if spot in spots_keep]

region = gamma_epi_wide.loc[final_spot_order].copy()
row_sum = region.sum(axis=1).replace(0, np.nan)
region = region.div(row_sum, axis=0).fillna(0.0)

LEGEND_ORDER = CELLTYPE_ORDER_2h
STACK_ORDER = list(reversed(CELLTYPE_ORDER_2h))

H_W, H_H = 5.2, 2.6
fig = plt.figure(figsize=(H_W + 2.6, H_H + 0.7))
ax = fig.add_axes([0.6 / (H_W + 2.6), 0.45 / (H_H + 0.7), H_W / (H_W + 2.6), H_H / (H_H + 0.7)])
x = np.arange(len(region))
bottom = np.zeros(len(region))

for cell_type in STACK_ORDER:
    vals = region[cell_type].to_numpy()
    ax.bar(
        x,
        vals,
        bottom=bottom,
        width=0.9,
        color=COLORS[cell_type],
        edgecolor="none",
    )
    bottom += vals

ax.set_xlim(-0.5, len(region) - 0.5)
ax.set_ylim(0.0, 1.0)
ax.set_xticks([])
ax.set_yticks([0, 1])
ax.tick_params(labelsize=FONT)
ax.set_xlabel("Epithelial region spots", fontsize=FONT)
ax.set_ylabel("Topic proportions", fontsize=FONT)
ax.set_facecolor("white")
fig.patch.set_facecolor("white")
for spine in ax.spines.values():
    spine.set_visible(True)

legend_handles = [
    plt.Rectangle((0, 0), 1, 1, color=COLORS[cell_type]) for cell_type in LEGEND_ORDER
]
ax.legend(
    legend_handles,
    [pretty(c) for c in LEGEND_ORDER],
    loc="center left",
    bbox_to_anchor=(1.01, 0.5),
    frameon=False,
    fontsize=FONT,
)
fig.savefig(f"{OUT}/Figure_2h.png", bbox_inches="tight", dpi=300)
plt.close(fig)