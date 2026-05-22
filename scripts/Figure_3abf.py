from matplotlib.colors import LinearSegmentedColormap as LSC
import matplotlib.pyplot as plt
import mtopic
import numpy as np
import pandas as pd
import os

INPUT_PATHS = {"trained_model": "../data/HumanPBMC_trained.h5mu"}
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
