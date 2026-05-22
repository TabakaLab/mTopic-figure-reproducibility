import matplotlib.pyplot as plt
import mudata as mu
import pandas as pd
import seaborn as sns
import os

INPUT_PATHS = {
    "cossim_E9b": "../data/HumanPBMC_cossim_per_gene.csv",
    "cossim_E9c": "../data/HumanPBMC_cossim_per_celltype.csv",
    "trained_model": "../data/HumanPBMC_trained.h5mu",
}

OUT = "../figures/Figure_E9"
os.makedirs(OUT, exist_ok=True)

# ═══════════════════════════════════════════════════════════════════════════
# Figure E9b
# ═══════════════════════════════════════════════════════════════════════════

df = pd.read_csv(INPUT_PATHS["cossim_E9b"], index_col=0)

fig, ax = plt.subplots(figsize=(2.2, 4.5))

sns.histplot(
    df["cosine_similarity"], bins=30, kde=False, ax=ax, color="seagreen", stat="percent"
)

ax.grid(False)
ax.set_xlabel("Number of cells")
ax.set_ylabel("Percentage of genes (%)")

fig.savefig(f"{OUT}/Figure_E9b.png", dpi=300, bbox_inches="tight")


# ═══════════════════════════════════════════════════════════════════════════
# Figure E9c
# ═══════════════════════════════════════════════════════════════════════════

mdata = mu.read_h5mu(INPUT_PATHS["trained_model"])
topics = mdata.obsm["topics"]
mdata.obs["celltype"] = topics.idxmax(axis=1).map(mdata.uns["TOPIC_CELLTYPE"]).values

counts = (
    mdata.obs["celltype"].value_counts().rename_axis("celltype").to_frame("n_cells")
)
cossim = pd.read_csv(INPUT_PATHS["cossim_E9c"], index_col=0)["cosine_similarity"]

common = counts.index.intersection(cossim.index)
counts = counts.loc[common]
cossim = cossim.loc[common]

colors = [mdata.uns["CELLTYPE_COLOR"][ct] for ct in counts.index]

fig, ax = plt.subplots(figsize=(3, 4.5))
ax.scatter(counts["n_cells"], cossim, c=colors, s=50, zorder=2)

ax.axvline(50, color="gray", linestyle="--", linewidth=1, zorder=1)

ax.set_xlabel("Number of cells")
ax.set_ylabel("Cosine similarity")

fig.savefig(f"{OUT}/Figure_E9c.png", dpi=300, bbox_inches="tight")
