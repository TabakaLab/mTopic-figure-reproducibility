from matplotlib.colors import LinearSegmentedColormap as LSC
import matplotlib.gridspec as gridspec
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import mtopic
import numpy as np
import pandas as pd
import os

INPUT_PATHS = {
    "trained_model": "../data/P22MouseBrainATAC_trained.h5mu",
    "raw_data": "../data/P22MouseBrainATAC_filtered.h5mu",
}
OUT = "../figures/Figure_1"
os.makedirs(OUT, exist_ok=True)

mdata = mtopic.read.h5mu(INPUT_PATHS["trained_model"])
P22ATAC_TOPIC_CELLTYPE = mdata.uns["TOPIC_CELLTYPE"]
P22ATAC_CELLTYPE_COLOR = mdata.uns["CELLTYPE_COLOR"]

# ═══════════════════════════════════════════════════════════════════════════
# Figure 1b
# ═══════════════════════════════════════════════════════════════════════════
mtopic.pl.dominant_topics(
    mdata,
    x="coords",
    s=20,
    figsize=(8, 4),
    legend=True,
    palette=mdata.uns["TOPIC_COLOR"],
    annotation=mdata.uns["TOPIC_CELLTYPE"],
    save=f"{OUT}/Figure_1b.png",
    legend_ncol=3,
    fontsize=6,
    markerscale=2,
)


# ═══════════════════════════════════════════════════════════════════════════
# Figure 1c
# ═══════════════════════════════════════════════════════════════════════════
mdata = mtopic.read.h5mu(INPUT_PATHS["trained_model"])
mtopic.tl.zscores(
    mdata,
    raw_data_path=INPUT_PATHS["raw_data"],
    mod="rna",
    n_top=20,
)
mtopic.pl.corr_heatmap(
    arr1=mdata.obsm["topics"],
    label1="Topic-spot distribution",
    arr2=mdata.mod["rna"].obsm["zscores"],
    label2="RNA signatures",
    cmap=LSC.from_list(
        "continuous_cmap",
        [
            (0, "#08306B"),
            (0.3, "#1562A9"),
            (0.5, "#FFFFFF"),
            (0.7, "#DE2B25"),
            (1, "#99000D"),
        ],
    ),
    save=f"{OUT}/Figure_1c.png",
)


# ═══════════════════════════════════════════════════════════════════════════
# Figure 1d
# ═══════════════════════════════════════════════════════════════════════════
mdata = mtopic.read.h5mu(INPUT_PATHS["trained_model"])
mtopic.tl.zscores(
    mdata,
    raw_data_path=INPUT_PATHS["raw_data"],
    mod="atac",
    n_top=100,
)
mtopic.pl.corr_heatmap(
    arr1=mdata.obsm["topics"],
    label1="Topic-spot distribution",
    arr2=mdata.mod["atac"].obsm["zscores"],
    label2="ATAC signatures",
    cmap=LSC.from_list(
        "continuous_cmap",
        [
            (0, "#08306B"),
            (0.3, "#1562A9"),
            (0.5, "#FFFFFF"),
            (0.7, "#78C679"),
            (1, "#005A32"),
        ],
    ),
    save=f"{OUT}/Figure_1d.png",
)


# ═══════════════════════════════════════════════════════════════════════════
# Figure 1e
# ═══════════════════════════════════════════════════════════════════════════
mdata = mtopic.read.h5mu(INPUT_PATHS["trained_model"])

CELLTYPE_ORDER = [
    "Excitatory neurons-1",
    "Excitatory neurons-2",
    "Excitatory neurons-3",
    "Excitatory neurons-4",
    "Excitatory neurons-5",
    "Excitatory neurons-6",
    "Excitatory neurons-7",
    "Excitatory neurons-8",
    "Excitatory neurons-9",
    "Excitatory neurons-10",
    "Excitatory neurons-11",
    "Excitatory neurons-12",
    "Excitatory neurons-13",
    "Excitatory neurons-14",
    "Excitatory neurons-15",
    "Excitatory neurons-16",
    "Excitatory neurons-17",
    "Excitatory neurons-18",
    "Inhibitory neurons-3",
    "Inhibitory neurons-4",
    "Inhibitory neurons-5",
    "Inhibitory neurons-6",
    "Inhibitory neurons-7",
    "Inhibitory neurons-8",
    "Medium spiny neurons-1",
    "Medium spiny neurons-2",
    "Medium spiny neurons-3",
    "Medium spiny neurons-4",
    "Medium spiny neurons-5",
    "Medium spiny neurons-6",
    "Medium spiny neurons-7",
    "Medium spiny neurons-8",
    "Oligodendrocytes",
    "Inhibitory neurons-2",
    "Inhibitory neurons-1",
    "Vascular-1",
    "Neuronal progenitors",
    "Ependymal cells-1",
    "Astrocytes",
    "Ependymal cells-2",
    "Microglia-1",
    "Microglia-2",
    "Perivascular macrophages-1",
    "Perivascular macrophages-2",
    "Perivascular macrophages-3",
    "Perivascular macrophages-4",
    "Perivascular macrophages-5",
    "OPC",
    "Pericytes",
    "Vascular-2",
]
celltype_rank = {ct: i for i, ct in enumerate(CELLTYPE_ORDER)}


def topic_sort_key(t):
    ann = P22ATAC_TOPIC_CELLTYPE.get(t, t)
    return (celltype_rank.get(ann, len(CELLTYPE_ORDER)), ann, t)


topic_order = sorted(mdata.obsm["topics"].columns, key=topic_sort_key)
topic_order = [t for t in topic_order if P22ATAC_TOPIC_CELLTYPE.get(t) != "Unknown"]

gamma = mdata.obsm["topics"][topic_order]
topic_annotation = P22ATAC_TOPIC_CELLTYPE
celltype_color = P22ATAC_CELLTYPE_COLOR
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
        ["celltype_order", "topic", "proportion"],
        ascending=[True, True, False],
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
    cmap=LSC.from_list(
        "topic_cmap",
        ["#000000", "#6600ff", "#cc0066", "#ff4400", "#ffaa00", "#ffff00"],
    ),
    vmin=0,
    vmax=1,
    interpolation="none",
    extent=[0, N, K, 0],
)
ax_heat.set_xlim(0, N)
ax_heat.set_xticks([])
ax_heat.set_yticks([])

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
    bbox_to_anchor=(0.5, -0.12),
)

plt.subplots_adjust(bottom=0.15)
plt.savefig(f"{OUT}/Figure_1e.png", bbox_inches="tight", dpi=300)
