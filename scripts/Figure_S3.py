import matplotlib.pyplot as plt
import mudata as md
import numpy as np
import os

INPUT_PATHS = {
    "mtopic": "../data/HumanTonsil_trained.h5mu",
}

OUT = "../figures/Figure_S3"
os.makedirs(OUT, exist_ok=True)


# ═══════════════════════════════════════════════════════════════════════════
# Figure S3
# ═══════════════════════════════════════════════════════════════════════════
mdata = md.read_h5mu(INPUT_PATHS["mtopic"])
topics_df = mdata.obsm["topics"]
coords = np.array(mdata.obsm["coords"])

topics = [f"topic_{i}" for i in range(1, 27) if f"topic_{i}" in topics_df.columns]

ncols = 7
nrows = int(np.ceil(len(topics) / ncols))

fig, axes = plt.subplots(
    nrows, ncols, figsize=(ncols * 2, nrows * 2.3), constrained_layout=True
)
axes = np.atleast_2d(axes).flatten()

bg_kwargs = dict(s=1, c="lightgrey", linewidths=0, rasterized=True)

for ax, topic in zip(axes, topics):
    ax.scatter(coords[:, 0], coords[:, 1], **bg_kwargs, zorder=1)
    ax.scatter(
        coords[:, 0],
        coords[:, 1],
        c=topics_df[topic].values,
        s=6,
        marker="h",
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

plt.savefig(f"{OUT}/Figure_S3.png", bbox_inches="tight", dpi=300)
plt.close()
