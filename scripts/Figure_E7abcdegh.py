import pickle

from matplotlib.colors import LinearSegmentedColormap as LSC, TwoSlopeNorm
from matplotlib.legend_handler import HandlerPatch
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
import mtopic
import numpy as np
import pandas as pd
from scipy.sparse import issparse
from sklearn.neighbors import NearestNeighbors
import os

INPUT_PATHS = {
    "trained_model": "../data/HumanPBMC_trained.h5mu",
    "trained_model_rna": "../data/HumanPBMC_RNA_trained.h5mu",
    "trained_model_rna_atac": "../data/HumanPBMC_RNA-ATAC_trained.h5mu",
    "neighborhood_overlap": "../data/HumanPBMC_neighborhood_overlap_nn200.pkl",
    "signature_dominance": "../data/HumanPBMC_signatures_diagonal_dominance_ratio.pkl",
    "resolution_sweep": "../data/Leiden_modality_test.csv",
}

OUT = "../figures/Figure_E7"
os.makedirs(OUT, exist_ok=True)

path_trained = INPUT_PATHS["trained_model"]
mdata = mtopic.read.h5mu(path_trained)
TOPIC_COLOR = mdata.uns["TOPIC_COLOR"]
TOPIC_CELLTYPE = mdata.uns["TOPIC_CELLTYPE"]
CELLTYPE_COLOR = mdata.uns["CELLTYPE_COLOR"]


method_col = {
    "cobolt": "#e1bb3e",
    "matilda": "#e35436",
    "mofa": "#7991a8",
    "mojitoo": "#7eba53",
    "mtopic": "#b322ad",
    "wnn": "#fec59a",
}

CMAP_PBMC = {
    "atac": LSC.from_list(
        "blue_green",
        [
            (0.0, "#ffffff"),
            (0.05, "#ffffff"),
            (0.4, "#78c679"),
            (0.7, "#005a32"),
            (1, "#000000"),
        ],
    ),
    "rna": LSC.from_list(
        "blue_red",
        [
            (0.0, "#ffffff"),
            (0.01, "#ffffff"),
            (0.3, "#de2b25"),
            (0.7, "#99000d"),
            (1, "#000000"),
        ],
    ),
    "prot": LSC.from_list(
        "prot",
        [(0.0, "#ffffff"), (0.4, "#cdac3b"), (0.7, "#967402"), (1, "#000000")],
    ),
}


class HandlerCircle(HandlerPatch):
    def create_artists(
        self, legend, orig_handle, xdescent, ydescent, width, height, fontsize, trans
    ):
        center = width / 2, height / 2
        p = mpatches.Circle(
            center,
            radius=width / 3,
            color=orig_handle.get_facecolor(),
            alpha=orig_handle.get_alpha(),
            edgecolor=None,
            linewidth=0,
        )
        self.update_prop(p, orig_handle, legend)
        p.set_transform(trans)
        return [p]


def hex_to_rgb(h):
    h = h.lstrip("#")
    return [int(h[i : i + 2], 16) / 255 for i in (0, 2, 4)]


# ═══════════════════════════════════════════════════════════════════════════
# Figure E7a
# ═══════════════════════════════════════════════════════════════════════════
with open(INPUT_PATHS["neighborhood_overlap"], "rb") as f:
    d = pickle.load(f)

d_sum = {key: np.sum([d[key][mod] for mod in d[key]]) for key in d}
d_plot = dict()
for method in ["cobolt", "matilda", "mofa", "mojitoo", "mtopic", "wnn"]:
    d_method = {key: d_sum[key] for key in d_sum if method in key}
    if d_method:
        max_key = max(d_method, key=d_method.get)
        d_plot[max_key] = d[max_key]

fig, ax = plt.subplots(figsize=(6, 5))

c = []
for key in d_plot:
    for method in method_col:
        if method in key:
            c.append(method_col[method])
            break

total_scores = {
    key: sum(d_plot[key][mod] for mod in ["atac", "rna", "prot"]) for key in d_plot
}
sorted_keys = sorted(total_scores, key=total_scores.get, reverse=True)
d_plot = {key: d_plot[key] for key in sorted_keys}
c = [next(method_col[m] for m in method_col if m in key) for key in sorted_keys]

scores = {
    mod: np.asarray([d_plot[key][mod] for key in d_plot])
    for mod in ["atac", "rna", "prot"]
}
totals = scores["atac"] + scores["rna"] + scores["prot"]
norm_factor = np.max(totals)
for mod in scores:
    scores[mod] = scores[mod] / norm_factor

ax.bar(range(len(d_plot)), scores["atac"], color=c, alpha=1)
ax.bar(range(len(d_plot)), scores["rna"], bottom=scores["atac"], color=c, alpha=0.7)
ax.bar(
    range(len(d_plot)),
    scores["prot"],
    bottom=scores["atac"] + scores["rna"],
    color=c,
    alpha=0.4,
)

circle1 = mpatches.Circle(
    (0, 0),
    radius=0.3,
    color="black",
    alpha=0.7,
    edgecolor=None,
    linewidth=0,
    label="ATAC",
)
circle2 = mpatches.Circle(
    (0, 0),
    radius=0.3,
    color="black",
    alpha=0.4,
    edgecolor=None,
    linewidth=0,
    label="RNA",
)
circle3 = mpatches.Circle(
    (0, 0),
    radius=0.3,
    color="black",
    alpha=0.1,
    edgecolor=None,
    linewidth=0,
    label="PROT",
)
ax.legend(
    handles=[circle3, circle2, circle1],
    loc="center left",
    bbox_to_anchor=(1, 0.5),
    title="Modality",
    frameon=False,
    handler_map={mpatches.Circle: HandlerCircle()},
    handletextpad=0.1,
    labelspacing=0.6,
)
ax.set_xticks(range(len(d_plot)))
ax.set_xticklabels(list(d_plot.keys()), rotation=90)
ax.set_ylabel("Normalized cumulative\nneighborhood overlap score")
plt.savefig(f"{OUT}/Figure_E7a.png", dpi=300, bbox_inches="tight")
plt.close()


# ═══════════════════════════════════════════════════════════════════════════
# Figure E7b
# ═══════════════════════════════════════════════════════════════════════════
with open(INPUT_PATHS["signature_dominance"], "rb") as f:
    raw_scores = pickle.load(f)
    raw_scores = dict(sorted(raw_scores.items()))

fig, ax = plt.subplots(figsize=(3, 6))

total_scores = {
    key: sum(raw_scores[key][mod] for mod in ["atac", "rna", "prot"])
    for key in raw_scores
}
sorted_keys = sorted(total_scores, key=total_scores.get, reverse=True)
c = [next(method_col[m] for m in method_col if m in key) for key in sorted_keys]
methods = sorted_keys
scores = {
    mod: np.asarray([raw_scores[key][mod] for key in sorted_keys])
    for mod in ["atac", "rna", "prot"]
}

ax.bar(range(len(scores["atac"])), scores["atac"], color=c, alpha=1)
ax.bar(
    range(len(scores["rna"])), scores["rna"], bottom=scores["atac"], color=c, alpha=0.7
)
ax.bar(
    range(len(scores["prot"])),
    scores["prot"],
    bottom=scores["atac"] + scores["rna"],
    color=c,
    alpha=0.4,
)

circle1 = mpatches.Circle(
    (0, 0),
    radius=0.3,
    color="black",
    alpha=0.7,
    edgecolor=None,
    linewidth=0,
    label="ATAC",
)
circle2 = mpatches.Circle(
    (0, 0),
    radius=0.3,
    color="black",
    alpha=0.4,
    edgecolor=None,
    linewidth=0,
    label="RNA",
)
circle3 = mpatches.Circle(
    (0, 0),
    radius=0.3,
    color="black",
    alpha=0.1,
    edgecolor=None,
    linewidth=0,
    label="Protein",
)
ax.legend(
    handles=[circle3, circle2, circle1],
    loc="center left",
    bbox_to_anchor=(1, 0.5),
    title="Modality",
    frameon=False,
    handler_map={mpatches.Circle: HandlerCircle()},
    handletextpad=0.1,
    labelspacing=0.6,
)
ax.set_xticks(range(len(methods)))
ax.set_xticklabels(methods, rotation=90)
ax.set_ylabel("Signature specificity score")
ax.set_xlim(-1, len(methods))
plt.savefig(f"{OUT}/Figure_E7b.png", dpi=300, bbox_inches="tight")
plt.close()


# ═══════════════════════════════════════════════════════════════════════════
# Figure E7c
# ═══════════════════════════════════════════════════════════════════════════
mdata_1 = mtopic.read.h5mu(INPUT_PATHS["trained_model_rna"])
mdata_2 = mtopic.read.h5mu(INPUT_PATHS["trained_model_rna_atac"])
mdata_3 = mdata

mtopic.pl.dominant_topics(
    mdata_1,
    x="umap",
    topics="topics_mm",
    palette=TOPIC_COLOR,
    annotation=TOPIC_CELLTYPE,
    markerscale=2,
    title="RNA",
    figsize=(8, 4),
    legend=True,
    s=50,
    legend_ncol=2,
    save=f"{OUT}/Figure_E7c_1.png",
)
mtopic.pl.dominant_topics(
    mdata_2,
    x="umap",
    topics="topics_mm",
    palette=TOPIC_COLOR,
    annotation=TOPIC_CELLTYPE,
    markerscale=2,
    title="RNA+ATAC",
    figsize=(8, 4),
    legend=True,
    s=50,
    legend_ncol=2,
    save=f"{OUT}/Figure_E7c_2.png",
)
mtopic.pl.dominant_topics(
    mdata_3,
    x="umap",
    palette=TOPIC_COLOR,
    annotation=TOPIC_CELLTYPE,
    markerscale=2,
    title="RNA+ATAC+PROT",
    figsize=(8, 4),
    legend=True,
    s=50,
    legend_ncol=2,
    save=f"{OUT}/Figure_E7c_3.png",
)


# ═══════════════════════════════════════════════════════════════════════════
# Figure E7d
# ═══════════════════════════════════════════════════════════════════════════
K_NN = 15

PROTEINS_LABEL = [
    "CD4-1",
    "CD4-2",
    "CD45RA",
    "CD38-1",
    "CD279",
    "CD45RO",
    "CD8",
    "CD57_Recombinant",
    "CD314",
    "CD56(NCAM)",
    "TCR_V_7.2",
    "CD19",
    "IgD",
    "CD123",
    "CD141",
    "CD11c",
    "CD86",
    "CD41",
    "CD278",
    "CD195",
    "CD69",
    "CD71",
    "CD25",
]

mdata_rna = mtopic.read.h5mu(INPUT_PATHS["trained_model_rna"])
mdata_rna_atac = mtopic.read.h5mu(INPUT_PATHS["trained_model_rna_atac"])
mdata_full_e7d = mtopic.read.h5mu(INPUT_PATHS["trained_model"])

CONDITIONS = [
    ("RNA + ATAC + Protein", mdata_full_e7d, "topics"),
    ("RNA + ATAC", mdata_rna_atac, "topics"),
    ("RNA", mdata_rna, "topics"),
]

prot = mdata_full_e7d.mod["prot"]
counts = prot.layers["counts"]
counts = counts.toarray() if issparse(counts) else np.asarray(counts)
counts = counts.astype(float)
log1p_ = np.log1p(counts)
clr = log1p_ - log1p_.mean(axis=1, keepdims=True)
prot_df = pd.DataFrame(
    clr, index=prot.obs_names, columns=[v.replace("prot:", "") for v in prot.var_names]
)


def knn_smoothness_per_protein(mdata, topics_key, prot_df, k=K_NN):
    common = mdata.obs_names.intersection(prot_df.index)
    X = np.asarray(mdata.obsm[topics_key])[mdata.obs_names.isin(common)]
    expr = prot_df.loc[mdata.obs_names[mdata.obs_names.isin(common)]].values
    idx = (
        NearestNeighbors(n_neighbors=k + 1)
        .fit(X)
        .kneighbors(X, return_distance=False)[:, 1:]
    )
    local_sq = ((expr[idx] - expr[:, None, :]) ** 2).mean(axis=1)
    local_var = local_sq.mean(axis=0)
    global_var = expr.var(axis=0, ddof=0)
    global_var = np.where(global_var > 0, global_var, np.nan)
    purity = 1 - local_var / (2 * global_var)
    return pd.Series(purity, index=prot_df.columns)


purity = pd.DataFrame(
    {name: knn_smoothness_per_protein(m, key, prot_df) for name, m, key in CONDITIONS}
)

gain = purity["RNA + ATAC + Protein"] - purity["RNA"]
purity = purity.loc[gain.sort_values(ascending=False).index]

purity_T = purity.T
col_mean = purity_T.mean(axis=0)
centered = purity_T.subtract(col_mean, axis=1)
row_mean = centered.mean(axis=1)

DIVERGING = LSC.from_list(
    "blue_white_red",
    [
        (0.0, "#153655"),
        (0.2, "#3579b8"),
        (0.5, "#f5f5f5"),
        (0.8, "#c12f2f"),
        (1.0, "#8a0202"),
    ],
)
c_max = max(float(np.nanmax(np.abs(centered.values))), 1e-3)
c_norm = TwoSlopeNorm(vcenter=0.0, vmin=-c_max, vmax=c_max)

n_rows, n_cols = centered.shape

fig, (ax, ax_mean) = plt.subplots(
    1,
    2,
    figsize=(10, 2),
    gridspec_kw={"width_ratios": [4, 1], "wspace": 0.15},
    sharey=True,
)
im = ax.imshow(
    centered.values,
    aspect="auto",
    cmap=DIVERGING,
    norm=c_norm,
    interpolation="nearest",
    interpolation_stage="rgba",
)

proteins_in_order = list(centered.columns)
tick_positions = [i for i, p in enumerate(proteins_in_order) if p in PROTEINS_LABEL]
tick_labels = [proteins_in_order[i] for i in tick_positions]
ax.set_xticks(tick_positions)
ax.set_xticklabels(tick_labels, rotation=90, fontsize=7)
ax.set_yticks(range(n_rows))
ax.set_yticklabels(centered.index)
ax.tick_params(top=False, bottom=False, left=False)
ax.set_title("Protein neighborhood coherence (PNC)")
for spine in ax.spines.values():
    spine.set_visible(False)

y = np.arange(n_rows)
ax_mean.scatter(
    row_mean.values,
    y,
    s=200,
    c=row_mean.values,
    cmap=DIVERGING,
    norm=c_norm,
    edgecolor="black",
    linewidth=1,
    zorder=2,
)
ax_mean.invert_yaxis()
ax_mean.set_xlabel("Mean relative\nPNC")
ax_mean.tick_params(left=False, labelleft=False)
ax_mean.spines["top"].set_visible(False)
ax_mean.spines["right"].set_visible(False)
ax_mean.spines["left"].set_visible(False)
ax_mean.set_xticks([-0.03, 0, 0.03])
ax_mean.set_xlim(-0.03, 0.03)
ax_mean.grid(True, which="both", color="#e0e0e0", lw=2, zorder=0)
ax_mean.set_axisbelow(True)

fig.tight_layout()
mean_pos = ax_mean.get_position()
cax = fig.add_axes([mean_pos.x1 + 0.04, mean_pos.y0, 0.015, mean_pos.height])
cbar = fig.colorbar(im, cax=cax)
cbar.set_label("Relative PNC")
cbar.outline.set_visible(False)

fig.savefig(f"{OUT}/Figure_E7d.png", bbox_inches="tight", dpi=300)
plt.close()


# ═══════════════════════════════════════════════════════════════════════════
# Figure E7e
# ═══════════════════════════════════════════════════════════════════════════
df_sweep = pd.read_csv(INPUT_PATHS["resolution_sweep"])
pbmc_sweep = df_sweep[df_sweep["dataset"].str.startswith("PBMC")].copy()

PAIRS = [
    ("PBMC_rna-atac", "RNA", "ATAC", "#b1b1b1"),
    ("PBMC_atac-prot", "ATAC", "Protein", "#dedede"),
    ("PBMC_rna-prot", "RNA", "Protein", "#8b8b8b"),
]
MODALITY_COLOR = {
    "ATAC": "#719C71",
    "RNA": "#D6483B",
    "Protein": "#DABF69",
}

fig, ax = plt.subplots(figsize=(6, 3))
seen_mm = set()
for dataset, label_a, label_b, base_color in PAIRS:
    sub = pbmc_sweep[pbmc_sweep["dataset"] == dataset].sort_values("res")
    if label_a not in seen_mm:
        ax.plot(
            sub["res"],
            sub["nmi_mm_a"],
            "-o",
            color=MODALITY_COLOR[label_a],
            label=f"MTM vs {label_a}",
            lw=2,
            ms=5,
        )
        seen_mm.add(label_a)
    if label_b not in seen_mm:
        ax.plot(
            sub["res"],
            sub["nmi_mm_b"],
            "-o",
            color=MODALITY_COLOR[label_b],
            label=f"MTM vs {label_b}",
            lw=2,
            ms=5,
        )
        seen_mm.add(label_b)
    ax.plot(
        sub["res"],
        sub["nmi_a_b"],
        "-o",
        color=base_color,
        label=f"{label_a} vs {label_b}",
        lw=1.8,
        ms=4,
    )

ax.set_xlabel("Clustering resolution")
ax.set_ylabel("NMI")
ax.set_xticks([0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0])
ax.set_ylim(0.3, 0.7)
ax.yaxis.set_major_locator(MultipleLocator(0.1))
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
ax.grid(True, axis="y", alpha=0.5)
ax.legend(frameon=False, fontsize=8, loc="center left", bbox_to_anchor=(1.02, 0.5))
fig.tight_layout()
fig.savefig(f"{OUT}/Figure_E7e.png", bbox_inches="tight", dpi=300)
plt.close()


# ═══════════════════════════════════════════════════════════════════════════
# Figure E7g
# ═══════════════════════════════════════════════════════════════════════════
CELLTYPES_KEEP_G = [
    "CD4 Treg",
    "CD4 Th2 Th17",
    "CD4 Tfh",
    "CD4 Th1",
    "CD4 Effector",
    "CD4 Activated",
    "CD4 Tregs double positive",
    "CD4 Memory",
    "CD4 Recently activated",
]
CELLTYPES_HIGHLIGHT = [
    "CD4 Treg",
    "CD4 Th2 Th17",
    "CD4 Th1",
    "CD4 Tfh",
    "CD4 Tregs double positive",
]
GRAY = "#bebebe"

mdata = mtopic.read.h5mu(path_trained)
topics = mdata.obsm["topics"]
dominant_topic = topics.idxmax(axis=1)
mdata.obs["celltype"] = dominant_topic.map(TOPIC_CELLTYPE).astype("category")

mdata = mdata[mdata.obs["celltype"].isin(CELLTYPES_KEEP_G)].copy()
umap = np.asarray(mdata.obsm["umap"])
celltype = mdata.obs["celltype"].astype(str).values
colors = np.array(
    [CELLTYPE_COLOR[c] if c in CELLTYPES_HIGHLIGHT else GRAY for c in celltype]
)
order = np.argsort(np.isin(celltype, CELLTYPES_HIGHLIGHT))

fig, ax = plt.subplots(figsize=(4, 4))
ax.scatter(umap[order, 0], umap[order, 1], c=colors[order], s=20, edgecolor="none")
ax.set_xticks([])
ax.set_yticks([])
ax.set_aspect("equal")
ax.set_xlim((-6.5, -0.5))
ax.set_ylim((-2, 4))
ax.axis("off")
fig.tight_layout()
plt.savefig(f"{OUT}/Figure_E7g.png", dpi=300)
plt.close()


# ═══════════════════════════════════════════════════════════════════════════
# Figure E7h
# ═══════════════════════════════════════════════════════════════════════════
CELLTYPES_KEEP_H = [
    "CD4 Tfh",
    "CD4 Th1",
    "CD4 Th2 Th17",
    "CD4 Treg",
    "CD4 Tregs double positive",
]
proteins_h = ["prot:CD278", "prot:CD195", "prot:CD69", "prot:CD71", "prot:CD25"]

mdata_pkg = mtopic.read.h5mu(INPUT_PATHS["trained_model"])

mod = "prot"
df_sig = mdata_pkg.mod[mod].varm["signatures"]
df = df_sig.div(df_sig.sum(axis=0), axis=1)
proteins_h = [p for p in proteins_h if p in df.index]

celltype_to_topic = {v: k for k, v in TOPIC_CELLTYPE.items()}
topic_order = [
    celltype_to_topic[ct] for ct in CELLTYPES_KEEP_H if ct in celltype_to_topic
]

df_plot = df.loc[proteins_h, topic_order].T
df_rank = df_plot.rank(ascending=False, axis=0)
n_rows, n_cols = df_plot.shape
max_size, min_size = 200, 5
df_size = min_size + (max_size - min_size) * (1 - (df_rank - 1) / (n_rows - 1))

fig, ax = plt.subplots(figsize=(3, 1.8))
rows, cols = np.meshgrid(range(n_rows), range(n_cols), indexing="ij")
sc = ax.scatter(
    cols.ravel(),
    rows.ravel(),
    c=df_plot.values.ravel(),
    s=df_size.values.ravel(),
    cmap=LSC.from_list(
        "prot",
        [
            (0.0, "#ffffff"),
            (0.4, "#cdac3b"),
            (0.5, "#967402"),
            (0.9, "#000000"),
            (1, "#000000"),
        ],
    ),
    vmin=0,
    linewidths=0,
)

fig.colorbar(sc, ax=ax, label="Probability", shrink=0.4)

ax.set_xticks(range(n_cols))
ax.set_xticklabels([p.replace("prot:", "") for p in df_plot.columns], rotation=90)
ax.set_yticks(range(n_rows))
ax.set_yticklabels([TOPIC_CELLTYPE.get(t, t) for t in df_plot.index])
ax.set_xlim(-0.5, n_cols - 0.5)
ax.set_ylim(-0.5, n_rows - 0.5)
ax.invert_yaxis()
plt.savefig(f"{OUT}/Figure_E7h_3.png", dpi=300, bbox_inches="tight")
plt.close()
