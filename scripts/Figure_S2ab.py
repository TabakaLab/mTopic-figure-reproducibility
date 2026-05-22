import matplotlib.pyplot as plt
import mtopic
import numpy as np
import os

INPUT_PATHS = {
    "P22H3K27me3": "../data/P22MouseBrainH3K27me3_trained.h5mu",
    "P22H3K27ac": "../data/P22MouseBrainH3K27ac_trained.h5mu",
}

OUT = "../figures/Figure_S2"
os.makedirs(OUT, exist_ok=True)

def all_topics_grid(sample, save_path):
    mdata = mtopic.read.h5mu(INPUT_PATHS[sample])
    topics_df = mdata.obsm["topics"]
    coords = np.array(mdata.obsm["coords"])

    topics = topics_df.columns.tolist()

    ncols = 10
    nrows = int(np.ceil(len(topics) / ncols))

    fig, axes = plt.subplots(
        nrows, ncols, figsize=(ncols * 1.5, nrows * 1.7), constrained_layout=True
    )
    axes = np.atleast_2d(axes).flatten()

    bg_kwargs = dict(s=1, c="lightgrey", linewidths=0, rasterized=True)

    for ax, topic in zip(axes, topics):
        ax.scatter(coords[:, 0], coords[:, 1], **bg_kwargs, zorder=1)
        ax.scatter(
            coords[:, 0],
            coords[:, 1],
            c=topics_df[topic].values,
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
        ax.set_title(f"Topic {topic.split('_')[1]}", fontsize=12)

    for ax in axes[len(topics) :]:
        ax.axis("off")

    plt.savefig(save_path, bbox_inches="tight", dpi=300)
    plt.close()


# ═══════════════════════════════════════════════════════════════════════════
# Figure S2a
# ═══════════════════════════════════════════════════════════════════════════
all_topics_grid("P22H3K27me3", f"{OUT}/Figure_S2a.png")


# ═══════════════════════════════════════════════════════════════════════════
# Figure S2b
# ═══════════════════════════════════════════════════════════════════════════
all_topics_grid("P22H3K27ac", f"{OUT}/Figure_S2b.png")
