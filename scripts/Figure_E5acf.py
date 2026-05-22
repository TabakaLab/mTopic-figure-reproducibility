import matplotlib.pyplot as plt
import mtopic
import muon as mu
import numpy as np
import pandas as pd
import os

mu.set_options(pull_on_update=False)

INPUT_PATHS = {
    "P22ATAC": "../data/P22MouseBrainATAC_trained.h5mu",
    "P22H3K27me3": "../data/P22MouseBrainH3K27me3_trained.h5mu",
    "umi_stats": "../data/P22MouseBrain_UMI_stats.csv",
    "atlas_celltypes": "../data/P22MouseBrain_AMBA_celltype_detection.csv",
}

OUT = "../figures/Figure_E5"
os.makedirs(OUT, exist_ok=True)


def plot_selected_topics(sample, path, topics, save_path):
    mdata = mtopic.read.h5mu(path)
    topics_df = mdata.obsm["topics"]
    coords = np.array(mdata.obsm["coords"])

    ncols = len(topics)
    nrows = 1
    fig, axes = plt.subplots(
        nrows, ncols, figsize=(ncols * 1.5, nrows * 1.7), constrained_layout=True
    )
    axes = np.atleast_1d(axes).flatten()

    bg_kwargs = dict(s=1, c="lightgrey", linewidths=0, rasterized=True)

    for ax, topic in zip(axes, topics):
        col = f"topic_{topic}"
        ax.scatter(coords[:, 0], coords[:, 1], **bg_kwargs, zorder=1)
        ax.scatter(
            coords[:, 0],
            coords[:, 1],
            c=topics_df[col].values,
            s=0.9,
            marker="s",
            cmap="gnuplot",
            linewidths=0,
            rasterized=True,
            zorder=2,
            vmin=0,
            vmax=1,
        )
        ax.axis("off")
        ax.set_aspect("equal")
        ax.set_title(f"{sample}\nTopic {topic}", fontsize=12)

    plt.savefig(save_path, bbox_inches="tight", dpi=300)
    plt.close()


# ═══════════════════════════════════════════════════════════════════════════
# Figure E5a
# ═══════════════════════════════════════════════════════════════════════════

df_umi = pd.read_csv(INPUT_PATHS["umi_stats"])
df_atlas = pd.read_csv(INPUT_PATHS["atlas_celltypes"])

MTOPIC_COLOR = "#b322ad"

DATASET_LABELS = {
    "P22_ATAC_RNA": "ATAC",
    "P22_H3K4me3_RNA": "H3K4me3",
    "P22_H3K27me3_RNA": "H3K27me3",
    "P22_H3K27ac_RNA": "H3K27ac",
}

merged = df_atlas.merge(df_umi[["dataset", "median_UMI"]], on="dataset", how="left")
mt = merged[merged["method"] == "mTopic"].sort_values("median_UMI")

fig, ax = plt.subplots(figsize=(5, 3))
ax.plot(
    mt["median_UMI"],
    mt["detected_types"],
    color=MTOPIC_COLOR,
    linewidth=1.2,
    linestyle="-",
    alpha=0.7,
    zorder=2,
)
ax.scatter(
    mt["median_UMI"],
    mt["detected_types"],
    color=MTOPIC_COLOR,
    s=70,
    edgecolor="white",
    linewidth=1.0,
    zorder=3,
)

for _, row in mt.iterrows():
    ax.annotate(
        DATASET_LABELS.get(row["dataset"], row["dataset"]),
        xy=(row["median_UMI"], row["detected_types"]),
        xytext=(6, 4),
        textcoords="offset points",
        fontsize=8,
        color="#444444",
    )

ax.set_xlabel("Median UMI per spot")
ax.set_ylabel("Detected cell types")
ax.grid(alpha=0.3, zorder=1)

fig.tight_layout()

fig.savefig(f"{OUT}/Figure_E5a.png", dpi=300, bbox_inches="tight")
plt.close(fig)

# ═══════════════════════════════════════════════════════════════════════════
# Figure E5c
# ═══════════════════════════════════════════════════════════════════════════
plot_selected_topics(
    "P22ATAC",
    INPUT_PATHS["P22ATAC"],
    topics=[10, 37, 12, 45],
    save_path=f"{OUT}/Figure_E5c_1.png",
)
plot_selected_topics(
    "P22H3K27me3",
    INPUT_PATHS["P22H3K27me3"],
    topics=[19, 2, 46, 48, 32, 1, 33, 3],
    save_path=f"{OUT}/Figure_E5c_2.png",
)


# ═══════════════════════════════════════════════════════════════════════════
# Figure E5f
# ═══════════════════════════════════════════════════════════════════════════
plot_selected_topics(
    "P22ATAC",
    INPUT_PATHS["P22ATAC"],
    topics=[31],
    save_path=f"{OUT}/Figure_E5f_1.png",
)
plot_selected_topics(
    "P22H3K27me3",
    INPUT_PATHS["P22H3K27me3"],
    topics=[14, 45],
    save_path=f"{OUT}/Figure_E5f_2.png",
)
