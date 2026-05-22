from matplotlib.colors import LinearSegmentedColormap as LSC
from matplotlib.lines import Line2D
import matplotlib.pyplot as plt
import mtopic
import numpy as np
import pandas as pd
from scipy.stats import gmean, zscore
import os

INPUT_PATHS = {
    "trained_model": "../data/HumanPBMC_trained.h5mu",
}

OUT = "../figures/Figure_E8"
os.makedirs(OUT, exist_ok=True)

mdata = mtopic.read.h5mu(INPUT_PATHS["trained_model"])
TOPIC_CELLTYPE = mdata.uns["TOPIC_CELLTYPE"]

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


# ═══════════════════════════════════════════════════════════════════════════
# Figure E8a
# ═══════════════════════════════════════════════════════════════════════════
cmap_e8a = LSC.from_list(
    "blue_red",
    [
        (0.0, "#ffffff"),
        (0.01, "#ffffff"),
        (0.3, "#de2b25"),
        (0.7, "#99000d"),
        (1, "#000000"),
    ],
)

celltype_to_topic = {v: k for k, v in TOPIC_CELLTYPE.items()}
topic_order = [
    celltype_to_topic[ct] for ct in CELLTYPE_ORDER if ct in celltype_to_topic
]

n_top = {"rna": 50}
mod = "rna"

df_sig = mdata.mod[mod].varm["signatures"]
df = df_sig.div(df_sig.sum(axis=0), axis=1)
df = df[topic_order]

top_features = pd.Index(
    pd.concat([df.nlargest(n_top[mod], col).index.to_series() for col in df.columns])
)
df_plot = df.loc[top_features]

fig, ax = plt.subplots(figsize=(10, 5))
im = ax.imshow(df_plot.values, aspect="auto", cmap=cmap_e8a, interpolation="none")

ax.set_xticks(range(df_plot.shape[1]))
ax.set_xticklabels([TOPIC_CELLTYPE.get(t, t) for t in topic_order], rotation=90)
ax.set_yticks([])
ax.set_ylabel("RNA")

cbar = fig.colorbar(im, ax=ax, label="Probability")
cbar.set_ticks([im.norm.vmin, im.norm.vmax])
plt.tight_layout()
plt.savefig(f"{OUT}/Figure_E8a.png", dpi=300, bbox_inches="tight")
plt.close()


# ═══════════════════════════════════════════════════════════════════════════
# Figure E8d
# ═══════════════════════════════════════════════════════════════════════════
mdata_pkg = mdata
prot_adata = mdata_pkg.mod["prot"]

raw_prot = prot_adata.layers["counts"].toarray().astype(float)
gm = gmean(raw_prot + 1, axis=1, keepdims=True)
clr_prot = np.log((raw_prot + 1) / gm)
scaled_prot = zscore(clr_prot, axis=0)
scaled_prot = np.clip(scaled_prot, -2, 2)

gamma = mdata_pkg.obsm["topics"]
cell_topics = gamma.idxmax(axis=1)

proteins = [
    "prot:CD4-1",
    "prot:CD4-2",
    "prot:CD45RA",
    "prot:CD38-1",
    "prot:CD279",
    "prot:CD45RO",
    "prot:CD8",
    "prot:CD57_Recombinant",
    "prot:CD314",
    "prot:CD56(NCAM)",
    "prot:TCR_V_7.2",
    "prot:CD19",
    "prot:IgD",
    "prot:CD123",
    "prot:CD141",
    "prot:CD11c",
    "prot:CD86",
    "prot:CD41",
]
var_names = list(prot_adata.var_names)
prot_idx = [var_names.index(p) for p in proteins if p in var_names]
prot_labels = [p.replace("prot:", "") for p in proteins if p in var_names]

celltype_to_topic = {v: k for k, v in TOPIC_CELLTYPE.items()}
topic_order = [
    celltype_to_topic[ct] for ct in CELLTYPE_ORDER if ct in celltype_to_topic
]

expr = scaled_prot[:, prot_idx]
counts = raw_prot[:, prot_idx]

mean_expr = pd.DataFrame(index=topic_order, columns=prot_labels, dtype=float)
pct_expr = pd.DataFrame(index=topic_order, columns=prot_labels, dtype=float)

for topic in topic_order:
    mask = (cell_topics == topic).values
    if mask.sum() == 0:
        mean_expr.loc[topic] = 0
        pct_expr.loc[topic] = 0
    else:
        mean_expr.loc[topic] = expr[mask].mean(axis=0)
        pct_expr.loc[topic] = (counts[mask] > 0).mean(axis=0) * 100

CMAP_E8D = LSC.from_list(
    "bwg",
    ["#153068", "#274992", "#FFFFFF", "#CAB157", "#967926"],
)

n_rows, n_cols = len(prot_labels), len(topic_order)
max_size, min_size = 200, 2

fig, ax = plt.subplots(figsize=(n_cols * 0.5 + 2, n_rows * 0.4 + 1))
rows, cols = np.meshgrid(range(n_rows), range(n_cols), indexing="ij")

sc = ax.scatter(
    cols.ravel(),
    rows.ravel(),
    c=mean_expr.T.values.ravel(),
    s=((pct_expr.T.values / 100) * max_size + min_size).ravel(),
    cmap=CMAP_E8D,
    vmin=-2,
    vmax=2,
    linewidths=0,
)
ax.set_xticks(range(n_cols))
ax.set_xticklabels([TOPIC_CELLTYPE.get(t, t) for t in topic_order], rotation=90)
ax.set_yticks(range(n_rows))
ax.set_yticklabels(prot_labels)
ax.set_xlim(-0.5, n_cols - 0.5)
ax.set_ylim(-0.5, n_rows - 0.5)
ax.invert_yaxis()

fig.colorbar(sc, ax=ax, label="Relative protein level", shrink=0.4)

for pct in [0, 25, 50, 75, 100]:
    size = (pct / 100) * max_size + min_size
    ax.add_line(
        Line2D(
            [],
            [],
            marker="o",
            color="none",
            markerfacecolor="grey",
            markeredgewidth=0,
            markersize=np.sqrt(size),
            label=str(pct),
        )
    )
ax.legend(
    title="Percent expressed", frameon=False, loc="upper left", bbox_to_anchor=(1.15, 1)
)
ax.set_title("Protein expression per topic")
plt.savefig(f"{OUT}/Figure_E8d.png", dpi=300, bbox_inches="tight")
plt.close()


# ═══════════════════════════════════════════════════════════════════════════
# Figure E8e
# ═══════════════════════════════════════════════════════════════════════════
mdata_pkg = mdata

mod = "prot"
df_sig = mdata_pkg.mod[mod].varm["signatures"]
df = df_sig.div(df_sig.sum(axis=0), axis=1)

proteins = [
    "prot:CD4-1",
    "prot:CD4-2",
    "prot:CD45RA",
    "prot:CD38-1",
    "prot:CD279",
    "prot:CD45RO",
    "prot:CD8",
    "prot:CD57_Recombinant",
    "prot:CD314",
    "prot:CD56(NCAM)",
    "prot:TCR_V_7.2",
    "prot:CD19",
    "prot:IgD",
    "prot:CD123",
    "prot:CD141",
    "prot:CD11c",
    "prot:CD86",
    "prot:CD41",
]
proteins = [p for p in proteins if p in df.index]

celltype_to_topic = {v: k for k, v in TOPIC_CELLTYPE.items()}
topic_order = [
    celltype_to_topic[ct] for ct in CELLTYPE_ORDER if ct in celltype_to_topic
]
df_plot = df.loc[proteins, topic_order]

df_rank = df_plot.rank(ascending=False, axis=0)
n_rows, n_cols = df_plot.shape
max_size, min_size = 200, 5
df_size = min_size + (max_size - min_size) * (1 - (df_rank - 1) / (n_rows - 1))

fig, ax = plt.subplots(figsize=(n_cols * 0.5 + 2, n_rows * 0.4 + 1))
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
    vmax=np.percentile(df_plot.values, 99),
    linewidths=0,
)

ax.set_xticks(range(n_cols))
ax.set_xticklabels([TOPIC_CELLTYPE.get(t, t) for t in topic_order], rotation=90)
ax.set_yticks(range(n_rows))
ax.set_yticklabels([p.replace("prot:", "") for p in df_plot.index])
ax.set_xlim(-0.5, n_cols - 0.5)
ax.set_ylim(-0.5, n_rows - 0.5)
ax.invert_yaxis()

fig.colorbar(sc, ax=ax, label="Probability", shrink=0.4)

legend_elements = [
    Line2D(
        [0],
        [0],
        marker="o",
        color="none",
        markerfacecolor="grey",
        markeredgewidth=0,
        markersize=np.sqrt(
            min_size + (max_size - min_size) * (1 - (rank - 1) / (n_rows - 1))
        ),
        label=f"Rank {rank}",
    )
    for rank in [1, 5, 10, 20]
]
ax.legend(
    handles=legend_elements,
    title="Rank",
    frameon=False,
    loc="upper left",
    bbox_to_anchor=(1.15, 1),
)
ax.set_title("Protein signatures per topic")
plt.savefig(f"{OUT}/Figure_E8e.png", dpi=300, bbox_inches="tight")
plt.close()
