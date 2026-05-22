from collections import defaultdict
import re

import matplotlib.colors as mcolors
from matplotlib.colors import LinearSegmentedColormap
import matplotlib.gridspec as gridspec
import matplotlib.patches as mpatches
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


# ═══════════════════════════════════════════════════════════════════════════
# Figure 2a
# ═══════════════════════════════════════════════════════════════════════════

mtopic.pl.scatter_pie(
    mdata,
    x="coords",
    radius=0.0074,
    palette=TOPIC_COLOR,
    annotation=TOPIC_CELLTYPE,
    figsize=(8, 4),
    legend_ncol=2,
    legend=True,
    save=f"{OUT}/Figure_2a_1.png",
)

mtopic.pl.scatter_pie(
    mdata,
    x="coords",
    radius=0.0074,
    xrange=[0.37, 0.52],
    yrange=[0.45, 0.63],
    palette=TOPIC_COLOR,
    annotation=TOPIC_CELLTYPE,
    legend=True,
    figsize=(8, 4),
    legend_ncol=2,
    save=f"{OUT}/Figure_2a_2.png",
)


# ═══════════════════════════════════════════════════════════════════════════
# Figure 2b
# ═══════════════════════════════════════════════════════════════════════════
mdata = mtopic.read.h5mu(path_trained)

ZONE_ORDER = [
    "Light zone",
    "Dark zone",
    "Mantle zone",
    "T cell zone",
    "Epithelial",
    "Myeloid",
    "Plasma cells",
    "MRC",
]

topics_df = mdata.obsm["topics"]
coords = np.array(mdata.obsm["coords"])

base_to_topics = defaultdict(list)
for topic, celltype in TOPIC_CELLTYPE.items():
    if celltype == "Unknown":
        continue
    base = next((z for z in ZONE_ORDER if celltype.startswith(z)), celltype)
    base_to_topics[base].append(topic)

summed = {}
for base, topic_list in base_to_topics.items():
    valid = [t for t in topic_list if t in topics_df.columns]
    summed[base] = topics_df[valid].sum(axis=1).values

ordered_items = sorted(
    summed.items(),
    key=lambda x: ZONE_ORDER.index(x[0]) if x[0] in ZONE_ORDER else len(ZONE_ORDER),
)

bg_kwargs = dict(s=2, c="lightgrey", linewidths=0, rasterized=True)

n_a = len(ordered_items)
fig, axes = plt.subplots(1, n_a, figsize=(3 * n_a, 3.2), constrained_layout=True)
axes = np.atleast_1d(axes)
for ax, (base, weights) in zip(axes, ordered_items):
    ax.scatter(coords[:, 0], coords[:, 1], **bg_kwargs, zorder=1)
    sc_a = ax.scatter(
        coords[:, 0],
        coords[:, 1],
        c=weights,
        s=12,
        marker="h",
        cmap=CMAP,
        linewidths=0,
        rasterized=True,
        zorder=2,
        vmin=0,
        vmax=1,
    )
    ax.set_title(base, fontsize=8, pad=4)
    ax.axis("off")
    ax.set_aspect("equal")
cbar = fig.colorbar(sc_a, ax=axes, shrink=0.7, fraction=0.02, pad=0.01)
cbar.set_label("Topic weight", fontsize=8)
cbar.ax.tick_params(labelsize=7)
plt.savefig(f"{OUT}/Figure_2b_1.png", bbox_inches="tight", dpi=300)
plt.close()

individual_celltypes = ["T cell zone"]
only_celltypes = ["T cell zone-1", "T cell zone-2", "T cell zone-3", "T cell zone-4"]

base_to_topics_b = defaultdict(list)
for topic, celltype in TOPIC_CELLTYPE.items():
    if celltype == "Unknown":
        continue
    base = next((z for z in ZONE_ORDER if celltype.startswith(z)), celltype)
    if base in individual_celltypes:
        base_to_topics_b[celltype].append(topic)
    else:
        base_to_topics_b[base].append(topic)

summed_b = {}
for base, topic_list in base_to_topics_b.items():
    valid = [t for t in topic_list if t in topics_df.columns]
    summed_b[base] = topics_df[valid].sum(axis=1).values
summed_b = {k: v for k, v in summed_b.items() if k in only_celltypes}

ordered_items_b = sorted(
    summed_b.items(),
    key=lambda x: (
        (
            ZONE_ORDER.index(next((z for z in ZONE_ORDER if x[0].startswith(z)), x[0]))
            if any(x[0].startswith(z) for z in ZONE_ORDER)
            else len(ZONE_ORDER)
        ),
        int(re.search(r"-(\d+)$", x[0]).group(1)) if re.search(r"-(\d+)$", x[0]) else 0,
    ),
)

n_b = len(ordered_items_b)
fig, axes = plt.subplots(1, n_b, figsize=(3 * n_b, 3.2), constrained_layout=True)
axes = np.atleast_1d(axes)
for ax, (base, weights) in zip(axes, ordered_items_b):
    ax.scatter(coords[:, 0], coords[:, 1], **bg_kwargs, zorder=1)
    sc_b = ax.scatter(
        coords[:, 0],
        coords[:, 1],
        c=weights,
        s=12,
        marker="h",
        cmap=CMAP,
        linewidths=0,
        rasterized=True,
        zorder=2,
        vmin=0,
        vmax=1,
    )
    ax.set_title(base, fontsize=8, pad=4)
    ax.axis("off")
    ax.set_aspect("equal")
cbar = fig.colorbar(sc_b, ax=axes, shrink=0.7, fraction=0.02, pad=0.01)
cbar.set_label("Topic weight", fontsize=8)
cbar.ax.tick_params(labelsize=7)
plt.savefig(f"{OUT}/Figure_2b_2.png", bbox_inches="tight", dpi=300)
plt.close()

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
labels = ["CD8A", "CD45RO", "CCR7", "VIM"]

raw = prot_adata.layers["counts"].toarray().astype(float)
gm = gmean(raw + 1, axis=1, keepdims=True)
clr = np.log((raw + 1) / gm)

n_c = len(proteins_to_plot)
fig, axes = plt.subplots(1, n_c, figsize=(3 * n_c, 3.2), constrained_layout=True)
axes = np.atleast_1d(axes)
for ax, protein, label in zip(axes, proteins_to_plot, labels):
    feat_idx = list(prot_adata.var_names).index(protein)
    expr = clr[:, feat_idx]
    vmin_val = expr.min()
    vmax_val = expr.max() if expr.max() > vmin_val else vmin_val + 1
    expr_n = (expr - vmin_val) / (vmax_val - vmin_val)

    ax.scatter(coords[:, 0], coords[:, 1], **bg_kwargs, zorder=1)
    sc_c = ax.scatter(
        coords[:, 0],
        coords[:, 1],
        c=expr_n,
        s=12,
        marker="h",
        cmap=protein_cmap,
        linewidths=0,
        rasterized=True,
        zorder=2,
        vmin=0,
        vmax=1,
    )
    ax.set_title(label, fontsize=12, pad=4)
    ax.axis("off")
    ax.set_aspect("equal")
cbar = fig.colorbar(sc_c, ax=axes, shrink=0.7, fraction=0.02, pad=0.01)
cbar.set_label("CLR (min–max scaled)", fontsize=8)
cbar.ax.tick_params(labelsize=7)
plt.savefig(f"{OUT}/Figure_2b_3.png", bbox_inches="tight", dpi=300)
plt.close()


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


def topic_sort_key(t):
    ann = TOPIC_CELLTYPE.get(t, t)
    for i, zone in enumerate(ZONE_ORDER_2c):
        if ann.startswith(zone):
            suffix = ann.replace(zone, "").strip("-").strip()
            num = int(suffix) if suffix.isdigit() else 0
            return (i, num)
    return (len(ZONE_ORDER_2c), 0)


topic_order = sorted(mdata.obsm["topics"].columns, key=topic_sort_key)
topic_order = [t for t in topic_order if TOPIC_CELLTYPE.get(t) != "Unknown"]

gamma = mdata.obsm["topics"][topic_order]
topic_annotation = TOPIC_CELLTYPE
celltype_color = CELLTYPE_COLOR
figsize = (18, 8)

dominant_topic = gamma.idxmax(axis=1)
dominant_celltype = dominant_topic.map(topic_annotation)
dominant_proportion = gamma.max(axis=1)

celltype_order = {
    ann: i
    for i, ann in enumerate(dict.fromkeys(topic_annotation[t] for t in topic_order))
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
    [
        plt.matplotlib.colors.to_rgb(celltype_color.get(ct, "#bebebe"))
        for ct in sorted_celltypes
    ]
)

fig = plt.figure(figsize=figsize)
gs = gridspec.GridSpec(
    2,
    2,
    figure=fig,
    height_ratios=[0.03, 1],
    width_ratios=[1, 0.02],
    hspace=0.01,
    wspace=0.02,
)
ax_bar = fig.add_subplot(gs[0, 0])
ax_heat = fig.add_subplot(gs[1, 0])
ax_cbar = fig.add_subplot(gs[:, 1])

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
    cmap=CMAP,
    vmin=0,
    vmax=1,
    interpolation="none",
    extent=[0, N, K, 0],
)
ax_heat.set_xlim(0, N)
ax_heat.set_yticks(np.arange(K) + 0.5)
ax_heat.set_yticklabels(
    [topic_annotation.get(t, t) for t in topic_order],
    fontsize=7,
)
ax_heat.set_xticks([])

fig.colorbar(im, cax=ax_cbar)
ax_cbar.set_ylabel("Topic-spot distribution")
ax_cbar.yaxis.set_label_position("right")
ax_cbar.yaxis.tick_right()

handles = [
    mpatches.Patch(color=color, label=ct) for ct, color in celltype_color.items()
]
fig.legend(
    handles=handles,
    loc="lower center",
    ncol=min(6, len(celltype_color)),
    fontsize=7,
    frameon=False,
    bbox_to_anchor=(0.5, -0.02),
)
plt.savefig(f"{OUT}/Figure_2c.png", bbox_inches="tight", dpi=300)
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
    if celltype in only_celltypes:
        celltype_to_topics[celltype].append(topic)

summed = {}
for celltype, topic_list in celltype_to_topics.items():
    valid = [t for t in topic_list if t in topics_df.columns]
    if valid:
        summed[celltype] = topics_df[valid].sum(axis=1).values

ordered_items = [(ct, summed[ct]) for ct in only_celltypes if ct in summed]

bg_kwargs = dict(s=2, c="lightgrey", linewidths=0, rasterized=True)

n = len(ordered_items)
ncols = 3
nrows = int(np.ceil(n / ncols))
fig, axes = plt.subplots(nrows, ncols, figsize=(3 * ncols, 3.2 * nrows), constrained_layout=True)
axes_flat = np.atleast_1d(axes).flatten()
for ax, (celltype, weights) in zip(axes_flat, ordered_items):
    ax.scatter(coords[:, 0], coords[:, 1], **bg_kwargs, zorder=1)
    sc_d = ax.scatter(
        coords[:, 0],
        coords[:, 1],
        c=weights,
        s=12,
        marker="h",
        cmap=CMAP,
        linewidths=0,
        rasterized=True,
        zorder=2,
        vmin=0,
        vmax=1,
    )
    ax.set_title(celltype, fontsize=10, pad=4)
    ax.axis("off")
    ax.set_aspect("equal")
for ax in axes_flat[n:]:
    ax.axis("off")
cbar = fig.colorbar(sc_d, ax=axes_flat.tolist(), shrink=0.6, fraction=0.02, pad=0.01)
cbar.set_label("Topic weight", fontsize=9)
cbar.ax.tick_params(labelsize=8)
plt.savefig(f"{OUT}/Figure_2d.png", bbox_inches="tight", dpi=300)
plt.close()


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

bg_kwargs = dict(s=2, c="lightgrey", linewidths=0, rasterized=True)
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
fig, axes = plt.subplots(nrows, ncols, figsize=(3.8 * ncols, 3.2 * nrows), constrained_layout=True)
axes_flat = np.atleast_1d(axes).flatten()
for ax, (g, vals) in zip(axes_flat, panels):
    ax.scatter(coords[:, 0], coords[:, 1], **bg_kwargs, zorder=1)
    sc_plot = ax.scatter(
        coords[:, 0],
        coords[:, 1],
        c=vals,
        s=12,
        marker="h",
        cmap=CMAP_2e,
        linewidths=0,
        rasterized=True,
        zorder=2,
        vmin=0,
        vmax=np.quantile(vals, 0.99),
    )
    ax.set_title(g, fontsize=12)
    ax.axis("off")
    ax.set_aspect("equal")
    cbar = fig.colorbar(sc_plot, ax=ax, fraction=0.04, pad=0.02, shrink=0.6)
    cbar.outline.set_visible(False)
    cbar.ax.tick_params(labelsize=7)
for ax in axes_flat[n:]:
    ax.axis("off")
plt.savefig(f"{OUT}/Figure_2e.png", bbox_inches="tight", dpi=300)
plt.close()


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
        "epithelial apical",
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
        "epithelial transitional",
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
        "epithelial basal",
        [
            "rna:CLEC5A",
            "rna:KRT15",
            "rna:NGFR",
            "rna:FGF1",
            "rna:EGFR",
        ],
    ),
    (
        "myeloid epi",
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
        "myeloid stroma",
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
fig_w = max(12.0, 0.32 * n_cols + 5.0)
fig_h = 0.45 * n_rows + 2.0
fig, ax = plt.subplots(figsize=(fig_w, fig_h))

im = ax.imshow(
    mat_scaled.values,
    aspect="auto",
    cmap=cmap_2f,
    vmin=vmin,
    vmax=vmax,
    interpolation="nearest",
)

ax.set_yticks(np.arange(n_rows))
ax.set_yticklabels(mat_scaled.index, fontsize=9)
ax.yaxis.tick_right()
ax.yaxis.set_label_position("right")

ax.set_xticks(np.arange(n_cols))
ax.set_xticklabels([f"rna:{g}" for g in mat_scaled.columns], rotation=90, fontsize=8)

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
        -0.32,
        lab,
        transform=trans,
        ha="center",
        va="top",
        fontsize=10,
        clip_on=False,
    )

for spine in ax.spines.values():
    spine.set_linewidth(0.8)

fig.subplots_adjust(right=0.78)
cax = fig.add_axes([0.90, 0.30, 0.012, 0.45])
cbar = fig.colorbar(im, cax=cax)
cbar.set_label("Relative RNA score", fontsize=9)

fig.savefig(f"{OUT}/Figure_2f.png", bbox_inches="tight")
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

mtopic.pl.scatter_pie(
    mdata,
    x="coords",
    radius=0.0074,
    palette=palette,
    annotation=TOPIC_CELLTYPE,
    figsize=(8, 4),
    legend_ncol=2,
    legend=True,
    save=f"{OUT}/Figure_2g_1.png",
)

mtopic.pl.scatter_pie(
    mdata,
    x="coords",
    radius=0.0074,
    xrange=[0.4, 0.7],
    yrange=[0.29, 0.6],
    palette=palette,
    annotation=TOPIC_CELLTYPE,
    legend=True,
    figsize=(8, 4),
    legend_ncol=2,
    save=f"{OUT}/Figure_2g_2.png",
)


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

fig, ax = plt.subplots(figsize=(12, 3))
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
        label=cell_type,
    )
    bottom += vals

ax.set_xlim(-0.5, len(region) - 0.5)
ax.set_ylim(0.0, 1.0)
ax.set_xticks([])
ax.set_xlabel("Epithelial region spots")
ax.set_ylabel("Topic proportions")
ax.set_facecolor("white")
fig.patch.set_facecolor("white")
for spine in ax.spines.values():
    spine.set_visible(True)

legend_handles = [
    plt.Rectangle((0, 0), 1, 1, color=COLORS[cell_type]) for cell_type in LEGEND_ORDER
]
ax.legend(
    legend_handles,
    LEGEND_ORDER,
    title="Topic",
    loc="center left",
    bbox_to_anchor=(1.01, 0.5),
    frameon=False,
    fontsize=9,
)
fig.tight_layout()
fig.savefig(f"{OUT}/Figure_2h.png", bbox_inches="tight", dpi=300)
plt.close(fig)
