from matplotlib.collections import PatchCollection
from matplotlib.patches import Wedge
import matplotlib.pyplot as plt
import mtopic
import mudata as md
import numpy as np
import pandas as pd
import scanpy as sc
from scipy.optimize import linear_sum_assignment
import os

INPUT_PATHS = {
    "mtopic": "../data/HumanTonsil_trained.h5mu",
    "spatialglue": "../data/HumanTonsil_SpatialGlue.h5ad",
    "miso": "../data/HumanTonsil_MISO.h5ad",
    "coral": "../data/HumanTonsil_CORAL.h5ad",
}

OUT = "../figures/Figure_S3"
os.makedirs(OUT, exist_ok=True)


# ═══════════════════════════════════════════════════════════════════════════
# Figure S3a
# ═══════════════════════════════════════════════════════════════════════════
mdata = mtopic.read.h5mu(INPUT_PATHS["mtopic"])
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

plt.savefig(f"{OUT}/Figure_S3a.png", bbox_inches="tight", dpi=300)
plt.close()


# ═══════════════════════════════════════════════════════════════════════════
# Figure S3b
# ═══════════════════════════════════════════════════════════════════════════
def add_pie_scatter(
    ax, xs, ys, radius, proportions=None, colors=None, na_color="#cccccc"
):
    xs = np.asarray(xs, dtype=float)
    ys = np.asarray(ys, dtype=float)
    n = len(xs)

    if proportions is None:
        spot_colors = colors if colors is not None else [na_color] * n
        patches = [Wedge((x, y), radius, 0, 360) for x, y in zip(xs, ys)]
        ax.add_collection(
            PatchCollection(patches, facecolors=spot_colors, edgecolors="none")
        )
        return

    props = np.asarray(proportions, dtype=float)
    patches, patch_colors = [], []
    cum = np.zeros(n)
    for seg_idx in range(props.shape[1]):
        fracs = props[:, seg_idx]
        theta1 = 90.0 - cum * 360.0
        cum += fracs
        theta2 = 90.0 - cum * 360.0
        for i in np.where(fracs > 0)[0]:
            patches.append(Wedge((xs[i], ys[i]), radius, theta2[i], theta1[i]))
            patch_colors.append(colors[seg_idx])
    ax.add_collection(
        PatchCollection(patches, facecolors=patch_colors, edgecolors="none")
    )


mdata = md.read_h5mu(INPUT_PATHS["mtopic"])
mtopic_dom = np.argmax(mdata.obsm["topics"].values, axis=1)

sg = sc.read_h5ad(INPUT_PATHS["spatialglue"])
miso = sc.read_h5ad(INPUT_PATHS["miso"])
coral = sc.read_h5ad(INPUT_PATHS["coral"])

mdata.obs["spatialglue"] = sg.obs["SpatialGlue26"]
mdata.obs["miso"] = miso.obs["miso26"]
mdata.obs["coral"] = coral.obs["cluster"]

palette26 = [
    "#00cd00",
    "#caf0f8",
    "#fff8dc",
    "#ff69b4",
    "#cd853f",
    "#00ff7f",
    "#ffff00",
    "#00a1d8",
    "#cc0000",
    "#bebebe",
    "#dea873",
    "#ffc0cb",
    "#20b2aa",
    "#fff39c",
    "#8b2252",
    "#006400",
    "#ffe600",
    "#ffa5d2",
    "#8e24aa",
    "#03ff03",
    "#00b4d8",
    "#0000cd",
    "#00cd99",
    "#000000",
    "#cdba96",
    "#ff0000",
]

coords = np.array(mdata.obsm["coords"])
xs, ys = coords[:, 0], coords[:, 1]
radius = 0.0075 * max(xs.max() - xs.min(), ys.max() - ys.min())

topics_df = mdata.obsm["topics"]
topic_colors = palette26
topic_props = topics_df.values


def correlated_cmap_bijective(obs_labels, mtopic_dom, palette26):
    labels = pd.Series(obs_labels)
    unique_clusters = sorted(labels.dropna().unique(), key=lambda x: int(x))
    n_clusters = len(unique_clusters)

    sim = np.zeros((n_clusters, len(palette26)), dtype=int)
    for i, c in enumerate(unique_clusters):
        mask = (labels == c).values
        for t in mtopic_dom[mask].astype(int):
            sim[i, t] += 1

    row_ind, col_ind = linear_sum_assignment(-sim)
    return {c: palette26[col_ind[i]] for i, c in enumerate(unique_clusters)}


methods = ["mtopic", "spatialglue", "miso", "coral"]
col_titles = ["mTopic", "SpatialGlue", "MISO", "Coral"]

color_maps = {
    m: correlated_cmap_bijective(mdata.obs[m].values, mtopic_dom, palette26)
    for m in methods[1:]
}

fig, axes = plt.subplots(1, 4, figsize=(12, 4))

for col, (method, title) in enumerate(zip(methods, col_titles)):
    ax = axes[col]
    if method == "mtopic":
        add_pie_scatter(
            ax, xs, ys, radius, proportions=topic_props, colors=topic_colors
        )
    else:
        cmap = color_maps[method]
        spot_colors = [
            cmap.get(v, "#cccccc") if pd.notna(v) else "#cccccc"
            for v in mdata.obs[method]
        ]
        add_pie_scatter(ax, xs, ys, radius, colors=spot_colors)

    ax.set_xlim(xs.min() - radius, xs.max() + radius)
    ax.set_ylim(ys.min() - radius, ys.max() + radius)
    ax.set_aspect("equal")
    ax.set_title(title, fontsize=15)
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_facecolor("white")
    for spine in ax.spines.values():
        spine.set_visible(False)


plt.tight_layout(rect=[0, 0.07, 1, 1])
plt.savefig(f"{OUT}/Figure_S3b.png", bbox_inches="tight", dpi=300)
plt.close()
