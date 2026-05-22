import math
import pickle

import matplotlib.colors as mcolors
from matplotlib.colors import LinearSegmentedColormap, Normalize
import matplotlib.pyplot as plt
import mtopic
import muon as mu
import numpy as np
import pandas as pd
from scipy.stats import gmean, zscore
import os

mu.set_options(pull_on_update=False)

INPUT_PATHS = {
    "trained_model": "../data/HumanTonsil_trained.h5mu",
    "islands_results": "../data/HumanTonsil_follicle_islands.pkl",
}

OUT = "../figures/Figure_E6"
os.makedirs(OUT, exist_ok=True)

mdata = mu.read_h5mu(INPUT_PATHS["trained_model"])
rna_adata = mdata.mod["rna"]
prot_adata = mdata.mod["prot"]

TOPIC_CELLTYPE = mdata.uns["TOPIC_CELLTYPE"]
CELLTYPE_COLOR = mdata.uns["CELLTYPE_COLOR"]
topic_map = TOPIC_CELLTYPE
color_map = CELLTYPE_COLOR

with open(INPUT_PATHS["islands_results"], "rb") as f:
    R = pickle.load(f)

label_order = [
    "Light zone-1",
    "Light zone-2",
    "Light zone-3",
    "Dark zone",
    "Mantle zone-1",
    "Mantle zone-2",
    "Mantle zone-3",
    "Mantle zone-4",
    "Mantle zone-5",
    "Mantle zone-6",
    "Plasma cells",
    "T cell zone-1",
    "T cell zone-2",
    "T cell zone-3",
    "T cell zone-4",
    "Myeloid cells-1",
    "Myeloid cells-2",
    "Myeloid cells-3",
    "MRC",
    "Epithelial apical-1",
    "Epithelial apical-2",
    "Epithelial transitional-1",
    "Epithelial transitional-2",
    "Epithelial basal-1",
    "Epithelial basal-2",
]

gamma = pd.DataFrame(
    mdata.obsm["topics"],
    index=mdata.obs_names,
    columns=[f"topic_{i+1}" for i in range(mdata.obsm["topics"].shape[1])],
)
mdata.obs["max_topic"] = gamma.idxmax(axis=1)
mdata.obs["label"] = mdata.obs["max_topic"].map(topic_map)
mdata.obs["label"] = pd.Categorical(
    mdata.obs["label"],
    categories=label_order,
    ordered=True,
)

cell_order = mdata.obs.sort_values("label").index.tolist()


def hex_to_rgb(hex_color):
    hex_color = hex_color.lstrip("#")
    return [int(hex_color[i : i + 2], 16) / 255 for i in (0, 2, 4)]


# ═══════════════════════════════════════════════════════════════════════
# Figure E6a
# ═══════════════════════════════════════════════════════════════════════
raw_rna = rna_adata.layers["counts"].toarray().astype(float)
raw_rna = raw_rna / raw_rna.sum(axis=1, keepdims=True) * 1e4
log_rna = np.log1p(raw_rna)
scaled_rna = zscore(log_rna, axis=0)

rna_signatures = pd.DataFrame(
    rna_adata.varm["signatures"],
    index=rna_adata.var_names,
    columns=[f"topic_{i+1}" for i in range(rna_adata.varm["signatures"].shape[1])],
)
rna_signatures.columns = [topic_map.get(c, c) for c in rna_signatures.columns]

rna_var = list(rna_adata.var_names)
rna_labels = mdata.obs.loc[list(rna_adata.obs_names), "label"].values

rna_cmap = mcolors.LinearSegmentedColormap.from_list(
    "rna_heatmap",
    [
        "#022860",
        "#08306B",
        "#1562A9",
        "#FFFFFF",
        "#FFFFFF",
        "#EA6661",
        "#920713",
        "#4E0A10",
    ],
)

rna_features = []
for label in label_order:
    if label in rna_signatures.columns:
        topn = rna_signatures[label].sort_values(ascending=False).head(5).index.tolist()
        rna_features.extend(topn)
rna_features = list(dict.fromkeys(rna_features))
rna_features = [f for f in rna_features if f in rna_var]
feat_idx = [rna_var.index(f) for f in rna_features]

feat_mat = scaled_rna[:, feat_idx]
present_labels = [l for l in label_order if (rna_labels == l).sum() > 0]
mat_rna = np.vstack(
    [feat_mat[rna_labels == label].mean(axis=0) for label in present_labels]
)
mat_rna = np.clip(mat_rna, -2.5, 2.5)

fig_a, axes_a = plt.subplots(
    1,
    2,
    figsize=(14, 7),
    gridspec_kw={"width_ratios": [0.03, 1]},
    constrained_layout=True,
)
axes_a[0].imshow(
    np.array([[hex_to_rgb(color_map[l])] for l in present_labels]),
    aspect="auto",
)
axes_a[0].axis("off")

im_a = axes_a[1].imshow(
    mat_rna,
    aspect="auto",
    cmap=rna_cmap,
    vmin=-2.5,
    vmax=2.5,
    interpolation="none",
    rasterized=True,
)
axes_a[1].set_xticks(np.arange(len(rna_features)))
axes_a[1].set_xticklabels(rna_features, rotation=90, fontsize=7)
axes_a[1].set_yticks(np.arange(len(present_labels)))
axes_a[1].set_yticklabels(present_labels, fontsize=8)

cbar_a = fig_a.colorbar(im_a, ax=axes_a[1], shrink=0.4, pad=0.01)
cbar_a.set_label("Mean relative gene expression", fontsize=10)
cbar_a.set_ticks([-2, 0, 2])
plt.savefig(f"{OUT}/Figure_E6a.png", bbox_inches="tight", dpi=300)
plt.close(fig_a)


# ═══════════════════════════════════════════════════════════════════════════
# Figure E6b
# ═══════════════════════════════════════════════════════════════════════════
raw_prot = prot_adata.layers["counts"].toarray().astype(float)
gm = gmean(raw_prot + 1, axis=1, keepdims=True)
clr_prot = np.log((raw_prot + 1) / gm)
scaled_prot = zscore(clr_prot, axis=0)

prot_var = list(prot_adata.var_names)

prot_cells = [c for c in cell_order if c in prot_adata.obs_names]
ct_for_prot = mdata.obs.loc[prot_cells, "label"]

prot_df = pd.DataFrame(scaled_prot, index=prot_adata.obs_names, columns=prot_var)
prot_df = prot_df.loc[prot_cells]
prot_df["celltype"] = ct_for_prot.values
ct_mean = prot_df.groupby("celltype")[prot_var].mean()

ct_order = [c for c in label_order if c in ct_mean.index]
ct_mean = ct_mean.loc[ct_order]

peak_ct = ct_mean.idxmax(axis=0)
peak_ct = peak_ct.iloc[
    pd.Categorical(peak_ct, categories=ct_order, ordered=True).argsort()
]
prot_order = peak_ct.index.tolist()

mat_prot = ct_mean.loc[ct_order, prot_order].values
mat_prot = np.clip(mat_prot, -2.5, 2.5)

prot_labels = [
    p.replace("prot:", "").replace("PTPRC-1", "CD45RO").replace("PTPRC", "CD45RA")
    for p in prot_order
]

prot_cmap = mcolors.LinearSegmentedColormap.from_list(
    "prot_heatmap",
    ["#08306B", "#1562A9", "#FFFFFF", "#CDAC3B", "#967402"],
)

fig_b, axes_b = plt.subplots(
    1,
    2,
    figsize=(14, max(3, 0.4 * len(ct_order) + 1)),
    gridspec_kw={"width_ratios": [0.03, 1]},
    constrained_layout=True,
)
axes_b[0].imshow(
    np.array([[hex_to_rgb(color_map[c])] for c in ct_order]),
    aspect="auto",
)
axes_b[0].axis("off")

im_b = axes_b[1].imshow(
    mat_prot,
    aspect="auto",
    cmap=prot_cmap,
    vmin=-2.5,
    vmax=2.5,
    interpolation="none",
    rasterized=True,
)
axes_b[1].set_xticks(np.arange(len(prot_labels)))
axes_b[1].set_xticklabels(prot_labels, rotation=90, fontsize=8)
axes_b[1].set_yticks(np.arange(len(ct_order)))
axes_b[1].set_yticklabels(ct_order, fontsize=9)

cbar_b = fig_b.colorbar(im_b, ax=axes_b[1], shrink=0.6, pad=0.01)
cbar_b.set_label("Relative protein level", fontsize=10)
cbar_b.set_ticks([-2, 0, 2])
plt.savefig(f"{OUT}/Figure_E6b.png", bbox_inches="tight", dpi=300)
plt.close(fig_b)


# ═══════════════════════════════════════════════════════════════════════════
# Figure E6c
# ═══════════════════════════════════════════════════════════════════════════

mdata_full = mtopic.read.h5mu(INPUT_PATHS["trained_model"])
TOPIC_COLOR = mdata_full.uns["TOPIC_COLOR"]

mtopic.pl.scatter_pie(
    mdata_full,
    x="coords",
    radius=0.0074,
    xrange=[0.55, 0.69],
    yrange=[0.43, 0.64],
    palette=TOPIC_COLOR,
    annotation=TOPIC_CELLTYPE,
    legend=True,
    figsize=(8, 4),
    legend_ncol=2,
    save=f"{OUT}/Figure_E6c_1.png",
)

coords = R["coords"]
islands = R["initial_island_dict"]
labelled = True

fig, ax = plt.subplots(figsize=(6.5, 6.5))
ax.scatter(coords[:, 0], coords[:, 1], color="#bebebe", s=35, edgecolor="none")
for i, (name, ids) in enumerate(islands.items()):
    ax.scatter(coords[ids, 0], coords[ids, 1], color="#eaeaea", s=35, edgecolor="none")
    if labelled:
        ax.text(
            coords[ids, 0].mean(),
            coords[ids, 1].mean(),
            str(i + 1),
            fontsize=15,
            ha="center",
            va="center",
            color="black",
        )
ax.axis("off")
ax.set_aspect("equal")
fig.tight_layout()
fig.savefig(f"{OUT}/Figure_E6c_2.png", dpi=300, bbox_inches="tight")
plt.close()


# ═══════════════════════════════════════════════════════════════════════════
# Figure E6d
# ═══════════════════════════════════════════════════════════════════════════
ISLANDS_KEEP = {f"island_{i}" for i in (9, 20, 22, 24, 29, 38)}

n_sectors = R["config"]["n_sectors"]
theta = np.linspace(0, 2 * np.pi, n_sectors, endpoint=False)
inner_radius = 0.25

plot_idx = 1
for name, merged_sorted, island_max in R["all_plot_data"]:
    if name not in ISLANDS_KEEP:
        continue
    fig, ax = plt.subplots(figsize=(4, 4), subplot_kw=dict(polar=True))
    for t in theta:
        ax.plot([t, t], [inner_radius, 1.0], color="#bebebe", linewidth=1)
    for label, values, color in merged_sorted:
        v_closed = np.append(values, values[0])
        t_closed = np.append(theta, theta[0])
        v_scaled = (v_closed / island_max) * (1 - inner_radius) + inner_radius
        ax.plot(t_closed, v_scaled, color=color, linewidth=2)
    ax.set_xticks([])
    ax.set_yticks([])
    ax.grid(False)
    ax.set_title(name, fontsize=10)
    fig.tight_layout()
    fig.savefig(f"{OUT}/Figure_E6d_{plot_idx}.png", dpi=300, bbox_inches="tight")
    plt.close()
    plot_idx += 1


# ═══════════════════════════════════════════════════════════════════════════
# Figure E6e
# ═══════════════════════════════════════════════════════════════════════════
coords = R["coords"]
clusters = R["clusters"]
correlation_df = R["correlation_df"]
island_to_cluster = R["island_to_cluster"]
n_clusters = R["config"]["n_clusters"]
palette_e6e = plt.colormaps["tab10"].resampled(n_clusters)

custom_bwr = LinearSegmentedColormap.from_list(
    "custom_bwr",
    ["#08306B", "#1562A9", "#FFFFFF", "#FFFFFF", "#DE2B25", "#99000D"],
    N=256,
)

labelled = False
fig, ax = plt.subplots(figsize=(6.5, 6.5))
ax.scatter(coords[:, 0], coords[:, 1], color="#bebebe", s=35, edgecolor="none")
for idx, (name, ids) in enumerate(R["initial_island_dict"].items()):
    ax.scatter(
        coords[ids, 0],
        coords[ids, 1],
        s=35,
        color=palette_e6e(island_to_cluster[name] - 1),
        edgecolor="none",
    )
    if labelled:
        ax.text(
            coords[ids, 0].mean(),
            coords[ids, 1].mean(),
            str(idx + 1),
            fontsize=15,
            ha="center",
            va="center",
        )
ax.axis("off")
ax.set_aspect("equal")
fig.tight_layout()
fig.savefig(f"{OUT}/Figure_E6e_1.png", dpi=300, bbox_inches="tight")
plt.close(fig)

ordered_idx = clusters.sort_values().index
ordered_corr = correlation_df.loc[ordered_idx, ordered_idx]
cluster_colors = clusters.loc[ordered_idx].map(lambda x: palette_e6e(x - 1))

fig, ax = plt.subplots(figsize=(10, 10))
im = ax.imshow(ordered_corr, cmap=custom_bwr, vmin=-1, vmax=1)
for i, c in enumerate(cluster_colors):
    ax.add_patch(plt.Rectangle((i - 0.5, -1.5), 1, 0.5, color=c, clip_on=False))
    ax.add_patch(plt.Rectangle((-1.5, i - 0.5), 0.5, 1, color=c, clip_on=False))
ax.set_xticks([])
ax.set_yticks([])
cbar = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
cbar.set_label("Correlation")
fig.tight_layout()
fig.savefig(f"{OUT}/Figure_E6e_2.png", dpi=300, bbox_inches="tight")
plt.close(fig)


# ═══════════════════════════════════════════════════════════════════════════
# Figure E6f
# ═══════════════════════════════════════════════════════════════════════════
cluster_pca = R["cluster_pca"]["unit_unmerged"]
spoke_color_map = R["topic_to_color"]

custom_cmap_e6f = LinearSegmentedColormap.from_list(
    "custom_cmap",
    ["#08306B", "#1562A9", "#FFFFFF", "#DE2B25", "#99000D"],
    N=256,
)
c_inner, c_outer = 0.4, 1.0
ring_inner, ring_outer = 1.02, 1.06


def draw_radar(ax, pc1, labels, var_exp, title, markersize=8, lw=1.0):
    n = len(labels)
    t_spokes = np.linspace(0, 2 * np.pi, n, endpoint=False)
    a_step = 2 * np.pi / n
    abs_max = np.max(np.abs(pc1)) + 1e-8
    norm_ = Normalize(vmin=-abs_max, vmax=abs_max)
    colors = [custom_cmap_e6f(norm_(v)) for v in pc1]
    radii = np.abs(pc1)
    r = (radii / (radii.max() + 1e-6)) * (c_outer - c_inner) + c_inner
    t_closed = np.append(t_spokes, t_spokes[0])
    r_closed = np.append(r, r[0])
    ax.plot(t_closed, r_closed, color="black", linewidth=lw)
    for t_, r_, col in zip(t_spokes, r, colors):
        ax.plot(
            [t_],
            [r_],
            "o",
            color=col,
            markersize=markersize,
            markeredgecolor="black",
            markeredgewidth=0.8,
            zorder=10,
        )
    for t_ in t_spokes:
        ax.plot([t_, t_], [0, c_outer], linestyle="--", linewidth=0.4, color="gray")
    for radius in [c_inner, c_outer]:
        ax.plot(
            np.linspace(0, 2 * np.pi, 200),
            [radius] * 200,
            linestyle="--",
            linewidth=0.4,
            color="gray",
        )
    for i, label in enumerate(labels):
        arc = np.linspace(t_spokes[i] - a_step / 2, t_spokes[i] + a_step / 2, 100)
        ax.fill_between(
            arc,
            np.full_like(arc, ring_inner),
            np.full_like(arc, ring_outer),
            color=spoke_color_map.get(label, "#CCCCCC"),
            linewidth=0,
            zorder=1,
        )
    ax.spines["polar"].set_visible(False)
    ax.set_theta_direction(-1)
    ax.set_theta_zero_location("N")
    ax.set_xticks([])
    ax.set_yticks([])
    ax.grid(False)
    ax.set_title(f"{title}\nPC1={var_exp:.0%}", fontsize=8)


n = len(cluster_pca)
ncols = math.ceil(math.sqrt(n))
nrows = math.ceil(n / ncols)
fig, axes = plt.subplots(
    nrows, ncols, figsize=(4 * ncols, 4 * nrows), subplot_kw=dict(polar=True)
)
axes_flat = np.array(axes).flatten()
for ax_i, cid in enumerate(sorted(cluster_pca)):
    info = cluster_pca[cid]
    draw_radar(
        axes_flat[ax_i],
        info["pcs"]["PC1"].values,
        info["labels"],
        var_exp=info["var_exp"],
        title=f"cluster {cid}",
        markersize=6,
        lw=0.8,
    )
for ax in axes_flat[n:]:
    ax.set_visible(False)
fig.tight_layout()
fig.savefig(f"{OUT}/Figure_E6f.png", dpi=300, bbox_inches="tight")
plt.close()
