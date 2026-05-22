import pickle

from matplotlib.lines import Line2D
from matplotlib.patches import Patch
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
import numpy as np
import pandas as pd
from scipy.spatial import ConvexHull
import os

INPUT_PATHS = {
    "spatial_benchmarks_csv": "../data/Spatial_benchmarks.csv",
    "p22_atac_res_csv": "../data/P22MouseBrainATAC_metrics_res_comparison.csv",
    "p22_h3k4me3_res_csv": "../data/P22MouseBrainH3K4me3_metrics_res_comparison.csv",
    "p22_h3k27ac_res_csv": "../data/P22MouseBrainH3K27ac_metrics_res_comparison.csv",
    "p22_h3k27me3_res_csv": "../data/P22MouseBrainH3K27me3_metrics_res_comparison.csv",
    "atlas_celltypes": "../data/P22MouseBrain_AMBA_celltype_detection.csv",
    "p22_atac_diag": "../data/P22MouseBrainATAC_signatures_diagonal_dominance_ratio.pkl",
    "p22_h3k4me3_diag": "../data/P22MouseBrainH3K4me3_signatures_diagonal_dominance_ratio.pkl",
    "p22_h3k27me3_diag": "../data/P22MouseBrainH3K27me3_signatures_diagonal_dominance_ratio.pkl",
    "p22_h3k27ac_diag": "../data/P22MouseBrainH3K27ac_signatures_diagonal_dominance_ratio.pkl",
}

OUT = "../figures/Figure_E4"
os.makedirs(OUT, exist_ok=True)


# ═══════════════════════════════════════════════════════════════════════════
# Figure E4a
# ═══════════════════════════════════════════════════════════════════════════
df_e4a = pd.read_csv(INPUT_PATHS["spatial_benchmarks_csv"])

METHOD_COLOR = {
    "mtopic": "#b322ad",
    "spatialglue": "#6acbd0",
    "miso": "#be8046",
    "coral": "#1f47a3",
}
DATASET_MARKER = {
    "Human tonsil": "X",
    "P22 ATAC": "P",
    "P22 H3K4me3": "*",
    "P22 H3K27ac": "v",
    "P22 H3K27me3": "<",
}
METHOD_ORDER = ["spatialglue", "miso", "coral", "mtopic"]
DATASET_ORDER = [
    "P22 ATAC",
    "P22 H3K4me3",
    "P22 H3K27me3",
    "P22 H3K27ac",
    "Human tonsil",
]

legend_methods = [
    Line2D([0], [0], marker="s", color=c, lw=0, label=m, ms=10)
    for m, c in METHOD_COLOR.items()
]
legend_datasets = [
    Line2D([0], [0], marker=mk, color="gray", lw=0, label=d, ms=8)
    for d, mk in DATASET_MARKER.items()
]

for i, x_metric in enumerate(["CHAOS", "PAS"], start=1):
    fig, ax = plt.subplots(figsize=(5.5, 4.5))
    for method in METHOD_ORDER:
        sub = df_e4a[df_e4a["method"] == method]
        pts = sub[[x_metric, "LISI"]].values
        if len(pts) >= 3:
            hull = ConvexHull(pts)
            verts = pts[hull.vertices]
            ax.fill(
                verts[:, 0],
                verts[:, 1],
                facecolor=METHOD_COLOR[method],
                edgecolor="none",
                alpha=0.18,
                zorder=1,
            )
            closed = np.vstack([verts, verts[:1]])
            ax.plot(
                closed[:, 0],
                closed[:, 1],
                color=METHOD_COLOR[method],
                lw=1,
                alpha=0.7,
                zorder=1,
            )
        elif len(pts) == 2:
            ax.plot(
                pts[:, 0],
                pts[:, 1],
                color=METHOD_COLOR[method],
                alpha=0.4,
                lw=8,
                zorder=1,
                solid_capstyle="round",
            )
        for _, r in sub.iterrows():
            ax.scatter(
                r[x_metric],
                r["LISI"],
                color=METHOD_COLOR[method],
                marker=DATASET_MARKER[r["dataset"]],
                s=120,
                edgecolor="black",
                linewidths=0.5,
                zorder=2,
            )

    ax.invert_xaxis()
    ax.invert_yaxis()
    ax.set_xlabel(x_metric)
    ax.set_ylabel("LISI")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.legend(
        handles=legend_methods + legend_datasets,
        loc="lower left",
        bbox_to_anchor=(1.02, 0),
        frameon=False,
        fontsize=8,
    )
    fig.tight_layout()
    fig.savefig(f"{OUT}/Figure_E4a_{i}.png", bbox_inches="tight", dpi=300)
    plt.close()


# ═══════════════════════════════════════════════════════════════════════════
# Figure E4b
# ═══════════════════════════════════════════════════════════════════════════
DATASETS_E4B = {
    "P22 ATAC": INPUT_PATHS["p22_atac_res_csv"],
    "P22 H3K4me3": INPUT_PATHS["p22_h3k4me3_res_csv"],
    "P22 H3K27ac": INPUT_PATHS["p22_h3k27ac_res_csv"],
    "P22 H3K27me3": INPUT_PATHS["p22_h3k27me3_res_csv"],
}
RESOLUTION_MARKER = {
    "20 clusters": "d",
    "30 clusters": "^",
    "40 clusters": "s",
    "50 clusters": "p",
}

frames = []
for name, path in DATASETS_E4B.items():
    df = pd.read_csv(path).rename(columns={"dataset": "resolution"})
    df["dataset_name"] = name
    frames.append(df)
full = pd.concat(frames, ignore_index=True)

agg_mean = (
    full.groupby(["method", "resolution"])[["CHAOS", "PAS", "LISI"]]
    .mean()
    .reset_index()
)

legend_methods_b = [
    Line2D([0], [0], marker="o", color=c, lw=0, label=m, ms=8)
    for m, c in METHOD_COLOR.items()
]
legend_res_b = [
    Line2D([0], [0], marker=m, color="gray", lw=0, label=k, ms=8)
    for k, m in RESOLUTION_MARKER.items()
]

for i, x_metric in enumerate(["CHAOS", "PAS"], start=1):
    fig, ax = plt.subplots(figsize=(5.5, 4.5))
    for method in METHOD_ORDER:
        sub_all = full[full["method"] == method]
        sub_mean = agg_mean[agg_mean["method"] == method].sort_values("resolution")

        for _, r in sub_all.iterrows():
            ax.scatter(
                r[x_metric],
                r["LISI"],
                color=METHOD_COLOR[method],
                marker=RESOLUTION_MARKER[r["resolution"]],
                s=100,
                alpha=0.2,
                edgecolor="none",
                zorder=1,
            )
        for _, r in sub_mean.iterrows():
            ax.scatter(
                r[x_metric],
                r["LISI"],
                color=METHOD_COLOR[method],
                marker=RESOLUTION_MARKER[r["resolution"]],
                s=150,
                edgecolor="black",
                linewidths=0.5,
                zorder=3,
            )

    ax.invert_xaxis()
    ax.invert_yaxis()
    ax.set_xlabel(x_metric)
    ax.set_ylabel("LISI")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.legend(
        handles=legend_methods_b + legend_res_b,
        loc="lower left",
        bbox_to_anchor=(1.02, 0),
        frameon=False,
        fontsize=8,
    )
    fig.tight_layout()
    fig.savefig(f"{OUT}/Figure_E4b_{i}.png", bbox_inches="tight", dpi=300)
    plt.close()


# ═══════════════════════════════════════════════════════════════════════════
# Figure E4c
# ═══════════════════════════════════════════════════════════════════════════

df = pd.read_csv(INPUT_PATHS["atlas_celltypes"])

DATASET_ORDER = [
    "P22_ATAC_RNA",
    "P22_H3K4me3_RNA",
    "P22_H3K27me3_RNA",
    "P22_H3K27ac_RNA",
]
DATASET_LABELS = {
    "P22_ATAC_RNA": "ATAC",
    "P22_H3K4me3_RNA": "H3K4me3",
    "P22_H3K27me3_RNA": "H3K27me3",
    "P22_H3K27ac_RNA": "H3K27ac",
}
METHOD_COLORS = {
    "mTopic": "#b322ad",
    "SpatialGlue": "#6acbd0",
    "CORAL": "#1f47a3",
    "MISO": "#be8046",
}
METHOD_ORDER = ["mTopic", "SpatialGlue", "CORAL", "MISO"]


def _ordered_present(values, preferred_order):
    present = list(pd.Series(values).dropna().unique())
    ordered = [x for x in preferred_order if x in present]
    ordered += [x for x in sorted(present) if x not in ordered]
    return ordered


datasets = _ordered_present(df["dataset"], DATASET_ORDER)
methods = _ordered_present(df["method"], METHOD_ORDER)

mat = (
    df.pivot_table(
        index="method", columns="dataset", values="detected_types", aggfunc="first"
    )
    .reindex(index=methods, columns=datasets)
    .fillna(0)
)

n_datasets = len(datasets)
group_gap = 0.9
bar_width = 0.9

ymax = 1.1 * (mat.values.max())
text_pad = 0.02 * ymax

fig, ax = plt.subplots(figsize=(8.5, 5.2))

all_positions = []
all_ticklabels = []

for mi, method in enumerate(methods):
    group_start = mi * (n_datasets + group_gap)
    positions = group_start + np.arange(n_datasets)

    values = mat.loc[method, datasets].astype(float).values

    bars = ax.bar(
        positions,
        values,
        width=bar_width,
        color=METHOD_COLORS[method],
        edgecolor="none",
        zorder=3,
    )

    for b, v in zip(bars, values):
        ax.text(
            b.get_x() + b.get_width() / 2,
            v + text_pad,
            str(int(v)),
            ha="center",
            va="bottom",
            fontsize=9,
        )

    all_positions.extend(positions)
    all_ticklabels.extend([DATASET_LABELS.get(d, d) for d in datasets])


ax.set_xticks(all_positions)
ax.set_xticklabels(all_ticklabels, rotation=45, ha="right", fontsize=9)

ax.set_ylabel("Detected cell types")
ax.set_ylim(0, ymax)
ax.grid(alpha=0.3, axis="y", zorder=1)

if all_positions:
    ax.set_xlim(min(all_positions) - 0.8, max(all_positions) + 0.8)

legend_handles = [
    Patch(facecolor=METHOD_COLORS[method], edgecolor="none", label=method)
    for method in methods
]

ax.legend(
    handles=legend_handles,
    frameon=False,
    ncol=len(methods),
    loc="upper center",
    bbox_to_anchor=(0.5, 1.18),
)

fig.subplots_adjust(bottom=0.28, top=0.80)

fig.savefig(f"{OUT}/Figure_E4c.png", bbox_inches="tight", dpi=300)
plt.close()


# ═══════════════════════════════════════════════════════════════════════════
# Figure E4d
# ═══════════════════════════════════════════════════════════════════════════

dsets = ["p22_atac", "p22_h3k4me3", "p22_h3k27me3", "p22_h3k27ac"]
method_col = {"mofa": "#7991a8", "mtopic": "#b322ad"}

scores = {
    "mofa": [
        pickle.load(open(INPUT_PATHS[f"{dset}_diag"], "rb"))["mofa"] for dset in dsets
    ],
    "mtopic": [
        pickle.load(open(INPUT_PATHS[f"{dset}_diag"], "rb"))["mtopic"] for dset in dsets
    ],
}
scores = {
    key: [
        {k if k == "rna" else "peaks": v for k, v in instance.items()}
        for instance in value
    ]
    for key, value in scores.items()
}
modalities = ["peaks", "rna"]

method_totals = {
    method: sum(
        instance.get(mod, 0) for instance in instances for mod in ["peaks", "rna"]
    )
    for method, instances in scores.items()
}
sorted_methods = sorted(method_totals, key=method_totals.get, reverse=True)
scores = {method: scores[method] for method in sorted_methods}
methods = list(scores.keys())
c = [method_col[method] for method in methods]

bar_width = 0.2
x_positions = np.arange(len(methods))
num_instances = max(len(scores[m]) for m in methods)
x_offsets = np.linspace(-bar_width * 1.5, bar_width * 1.5, num_instances)

method_sums = [
    sum(instance.get(mod, 0) for instance in scores[method] for mod in modalities)
    for method in methods
]
x_margin = bar_width
xlim = (
    x_positions[0] + x_offsets[0] - x_margin,
    x_positions[-1] + x_offsets[-1] + x_margin,
)

fig, ax_line = plt.subplots(figsize=(6, 2.2))
for y in method_sums:
    ax_line.axhline(y=y, color="gray", linestyle="--", linewidth=1, alpha=0.5)
ax_line.scatter(x_positions, method_sums, color=c, s=100, zorder=3)
ax_line.set_ylabel("Cumulative\nscore", fontsize=12)
ax_line.spines[["right", "top"]].set_visible(False)
ax_line.set_yticks([0, 1, 2, 3])
ax_line.set_xticks(x_positions)
ax_line.set_xticklabels([m.upper() for m in methods], fontsize=10)
ax_line.set_ylim([0, max(method_sums) + 0.25])
ax_line.set_xlim(xlim)
for spine in ax_line.spines.values():
    spine.set_linewidth(1.5)
plt.tight_layout()
plt.savefig(f"{OUT}/Figure_E4d_1.png", dpi=300, bbox_inches="tight")
plt.close()

fig, ax_bar = plt.subplots(figsize=(6, 5))
xtick_positions = []
xtick_labels = []

for i, method in enumerate(methods):
    for j in range(len(scores[method])):
        bottom = 0
        label = f"{method.upper()} - {dsets[j].replace('p22_', '')}"
        x_pos = x_positions[i] + x_offsets[j]
        for mod, alpha in zip(["peaks", "rna"], [1, 0.7]):
            value = scores[method][j].get(mod, 0)
            ax_bar.bar(
                x_pos,
                value,
                width=bar_width * 0.9,
                bottom=bottom,
                color=c[i],
                alpha=alpha,
            )
            bottom += value
        xtick_positions.append(x_pos)
        xtick_labels.append(label)

ax_bar.set_xticks(xtick_positions)
ax_bar.set_xticklabels(xtick_labels, rotation=90, fontsize=9)
ax_bar.set_ylabel("Signature\nspecificity score", fontsize=10)
ax_bar.set_ylim([0, 0.8])
ax_bar.set_xlim(xlim)
ax_bar.spines[["right", "top"]].set_visible(False)
ax_bar.tick_params(axis="x", which="both", bottom=False, top=False)
ax_bar.yaxis.set_major_locator(MaxNLocator(nbins=5))
for spine in ax_bar.spines.values():
    spine.set_linewidth(1.5)
plt.xticks(fontsize=9)
plt.yticks(fontsize=10)
plt.tight_layout()
plt.savefig(f"{OUT}/Figure_E4d_2.png", dpi=300, bbox_inches="tight")
plt.close()
