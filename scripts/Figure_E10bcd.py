import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter, MultipleLocator
import mudata as mu
import numpy as np
import pandas as pd
import rdata
import seaborn as sns
import xarray as xr
import os

INPUT_PATHS = {
    "cossim_E10b": "../data/HumanPBMC_cossim_per_gene.csv",
    "cossim_E10c": "../data/HumanPBMC_cossim_per_celltype.csv",
    "trained_model": "../data/HumanPBMC_trained.h5mu",
    "grn_E10d": "../data/PBMC_cyto_GRN.RDS",
    "rte_E10d": "../data/RTE.RDS",
}

plt.rcParams["font.family"] = "Arial"
plt.rcParams["pdf.fonttype"] = 42

OUT = "../figures/Figure_E10"
os.makedirs(OUT, exist_ok=True)

# ═══════════════════════════════════════════════════════════════════════════
# Figure E10b
# ═══════════════════════════════════════════════════════════════════════════

df = pd.read_csv(INPUT_PATHS["cossim_E10b"], index_col=0)

fig, ax = plt.subplots(figsize=(2.2, 4.5))

sns.histplot(
    df["cosine_similarity"], bins=30, kde=False, ax=ax, color="seagreen", stat="percent"
)

ax.grid(False)
ax.set_xlabel("Number of cells")
ax.set_ylabel("Percentage of genes (%)")

fig.savefig(f"{OUT}/Figure_E10b.png", dpi=300, bbox_inches="tight")


# ═══════════════════════════════════════════════════════════════════════════
# Figure E10c
# ═══════════════════════════════════════════════════════════════════════════

mdata = mu.read_h5mu(INPUT_PATHS["trained_model"])
topics = mdata.obsm["topics"]
mdata.obs["celltype"] = topics.idxmax(axis=1).map(mdata.uns["TOPIC_CELLTYPE"]).values

counts = (
    mdata.obs["celltype"].value_counts().rename_axis("celltype").to_frame("n_cells")
)
cossim = pd.read_csv(INPUT_PATHS["cossim_E10c"], index_col=0)["cosine_similarity"]

common = counts.index.intersection(cossim.index)
counts = counts.loc[common]
cossim = cossim.loc[common]

colors = [mdata.uns["CELLTYPE_COLOR"][ct] for ct in counts.index]

fig, ax = plt.subplots(figsize=(3, 4.5))
ax.scatter(counts["n_cells"], cossim, c=colors, s=50, zorder=2)

ax.axvline(50, color="gray", linestyle="--", linewidth=1, zorder=1)

ax.set_xlabel("Number of cells")
ax.set_ylabel("Cosine similarity")

fig.savefig(f"{OUT}/Figure_E10c.png", dpi=300, bbox_inches="tight")


# ═══════════════════════════════════════════════════════════════════════════
# Figure E10d
# ═══════════════════════════════════════════════════════════════════════════

TOP_N = 2000
N_TOP_LABELS = 10
EXTRA_LABELS = {
    "TBX21": ["SAMD3", "LINC00299", "KLRG1", "GZMA"],
    "EOMES": ["PTPN22", "PTPRM", "RIN3", "TRGC2"],
    "RORA": [],
}
LABEL_SIZE = 24


def to_frame(x):
    if isinstance(x, pd.DataFrame):
        return x
    if isinstance(x, xr.DataArray):
        return pd.DataFrame(
            x.values,
            index=x.coords[x.dims[0]].values,
            columns=x.coords[x.dims[1]].values,
        )
    raise TypeError(f"Unexpected GRN type: {type(x)}")


def spread(target, step, lo, hi):
    n = len(target)
    offset = np.arange(n) * step
    blocks = []
    for v in -(target + offset):
        blocks.append([v, 1])
        while len(blocks) > 1 and blocks[-2][0] / blocks[-2][1] > blocks[-1][0] / blocks[-1][1]:
            total, count = blocks.pop()
            blocks[-1][0] += total
            blocks[-1][1] += count
    z = -np.concatenate([np.full(c, t / c) for t, c in blocks])
    z = np.clip(z, lo + (n - 1) * step, hi)
    return z - offset


def place_labels(ax, df, df_lab, fontsize, gap=0.15, line_spacing=1.25, isolation=1.5):
    df_lab = df_lab.sort_values("Score", ascending=False)
    n = len(df_lab)
    to_axes = ax.transAxes.inverted()
    pts = to_axes.transform(ax.transData.transform(df_lab[["Rank", "Score"]].values))
    curve = to_axes.transform(ax.transData.transform(df[["Rank", "Score"]].values))

    ax_h = ax.get_window_extent().height
    line_h = fontsize * line_spacing * ax.figure.dpi / 72 / ax_h

    x_col = min(pts[:, 0].max() + gap, 0.7)
    under = curve[curve[:, 0] >= x_col - 0.02, 1]
    lo = max(0.03, under.max() + line_h if len(under) else 0.03)
    hi = 0.97
    step = min(line_h, (hi - lo) / max(n - 1, 1))

    gaps = np.full(n, np.inf)
    dy = -np.diff(pts[:, 1])
    gaps[:-1] = np.minimum(gaps[:-1], dy)
    gaps[1:] = np.minimum(gaps[1:], dy)
    y = spread(pts[:, 1], step, lo, hi)
    isolated = (gaps >= isolation * step) & (np.abs(y - pts[:, 1]) < 0.75 * step)

    dx = fontsize * 0.4 * ax.figure.dpi / 72 / ax.get_window_extent().width

    for gene, (px, py), yi, iso in zip(df_lab.index, pts, y, isolated):
        if iso:
            ax.annotate(
                gene,
                xy=(px, py),
                xycoords="axes fraction",
                xytext=(px + dx, yi),
                textcoords="axes fraction",
                ha="left",
                va="center",
                fontsize=fontsize,
                fontstyle="italic",
            )
        else:
            ax.annotate(
                gene,
                xy=(px, py),
                xycoords="axes fraction",
                xytext=(x_col, yi),
                textcoords="axes fraction",
                ha="left",
                va="center",
                fontsize=fontsize,
                fontstyle="italic",
                arrowprops=dict(
                    arrowstyle="-", color="black", lw=1,
                    shrinkA=2, shrinkB=1, relpos=(0, 0.5),
                ),
            )


def gene_names(x):
    if isinstance(x, xr.DataArray):
        return list(x.coords[x.dims[0]].values)
    if isinstance(x, pd.Series):
        return list(x.index)
    if isinstance(x, dict):
        return list(x.keys())
    return list(x)


tf_x_gene = to_frame(rdata.read_rds(INPUT_PATHS["grn_E10d"]))
rte = rdata.read_rds(INPUT_PATHS["rte_E10d"])

panels = [
    ("TBX21", "TBX21--MA0690.2", "#00FF00", 0.0005),
    ("EOMES", "EOMES--MA0800.1", "#FF0000", 0.0005),
    ("RORA", "RORA--MA0072.1", "#0000FF", 0.005),
]

PANEL_W, PANEL_H = 5, 6
fig = plt.figure(figsize=(PANEL_W * len(panels), PANEL_H))
labels = []

for i, (tf, motif, color, y_step) in enumerate(panels):
    ax = fig.add_axes([
        (i * PANEL_W + 0.22 * PANEL_W) / (PANEL_W * len(panels)),
        0.12,
        0.73 * PANEL_W / (PANEL_W * len(panels)),
        0.80,
    ])

    scores = tf_x_gene.loc[motif].astype(float).sort_values(ascending=False)
    df = pd.DataFrame(
        {"Rank": np.arange(1, len(scores) + 1), "Score": scores.values},
        index=scores.index,
    ).iloc[:TOP_N]

    rte_all = gene_names(rte[tf])
    df_rte = df.loc[[g for g in rte_all if g in df.index]]

    lab_genes = rte_all[:N_TOP_LABELS] + EXTRA_LABELS[tf]
    df_lab = df.loc[[g for g in dict.fromkeys(lab_genes) if g in df.index]]

    ax.scatter(df["Rank"], df["Score"], s=8, color="black", linewidths=0, zorder=1)
    ax.scatter(df_rte["Rank"], df_rte["Score"], s=10, color=color, linewidths=0, zorder=2)

    ax.set_xticks([1, 1000, 2000], ["1", "1,000", "2,000"])
    ax.yaxis.set_major_locator(MultipleLocator(y_step))
    ax.yaxis.set_major_formatter(FormatStrFormatter("%g"))
    ax.grid(False)
    ax.spines[["top", "right"]].set_visible(False)
    ax.tick_params(labelsize=LABEL_SIZE)
    ax.set_title(tf, fontsize=LABEL_SIZE)

    if i == 0:
        ax.set_ylabel("TF-gene association score", fontsize=LABEL_SIZE)
    if i == 1:
        ax.set_xlabel("Rank-sorted genes", fontsize=LABEL_SIZE)

    labels.append((ax, df, df_lab))

fig.canvas.draw()

for ax, df, df_lab in labels:
    place_labels(ax, df, df_lab, LABEL_SIZE)

fig.savefig(f"{OUT}/Figure_E10d.png", dpi=300, bbox_inches="tight")
plt.close(fig)