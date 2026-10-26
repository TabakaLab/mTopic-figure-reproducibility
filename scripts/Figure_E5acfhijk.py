import matplotlib.pyplot as plt
import mtopic
import muon as mu
import numpy as np
import pandas as pd
import os

from matplotlib.cm import ScalarMappable
from matplotlib.colors import LinearSegmentedColormap, Normalize, to_rgb
from matplotlib.ticker import FormatStrFormatter, MultipleLocator
from scipy.stats import gaussian_kde

mu.set_options(pull_on_update=False)

INPUT_PATHS = {
    "P22ATAC": "../data/P22MouseBrainATAC_trained.h5mu",
    "P22H3K27me3": "../data/P22MouseBrainH3K27me3_trained.h5mu",
    "umi_stats": "../data/P22MouseBrain_UMI_stats.csv",
    "atlas_celltypes": "../data/P22MouseBrain_AMBA_celltype_detection.csv",
    "radius_sweep_metrics": "../data/radius_sweep_metrics.tsv",
    "corr_cache": "../data/topic_signature_corr_H3K27me3_cache.npz",
    "corr_similarity": "../data/topic_signature_corr_H3K27me3_similarity_all.csv",
    "rbo_per_topic": "../data/P22H3K27me3_rbo_rearranges_per_topic.csv",
    "rbo_preserved": "../data/P22H3K27me3_rbo_rearranges_preserved.csv",
    "theta_MTM": "../data/theta_H3K27me3_MTM.csv",
    "theta_r0p06": "../data/theta_H3K27me3_sMTM_r0p06.csv",
    "theta_r0p8": "../data/theta_H3K27me3_sMTM_r0p8.csv",
    "gini_moran_shift": "../data/gini_moran_H3K27me3_shift.csv",
}

OUT = "../figures/Figure_E5"
os.makedirs(OUT, exist_ok=True)


_TOPIC_CACHE = {}


def _topics_and_coords(path):
    if path not in _TOPIC_CACHE:
        mdata = mtopic.read.h5mu(path)
        _TOPIC_CACHE[path] = (mdata.obsm["topics"], np.array(mdata.obsm["coords"]))
    return _TOPIC_CACHE[path]


def plot_selected_topics(groups, save_path, panel_w=1.5, panel_h=1.7, gap=0.25):
    widths = []
    for gi, (_, _, topics) in enumerate(groups):
        if gi:
            widths.append(gap)
        widths += [1.0] * len(topics)

    fig = plt.figure(
        figsize=(sum(widths) * panel_w, panel_h), constrained_layout=True
    )
    gs = fig.add_gridspec(1, len(widths), width_ratios=widths)

    bg_kwargs = dict(s=1, c="lightgrey", linewidths=0, rasterized=True)

    col_i = 0
    for gi, (sample, path, topics) in enumerate(groups):
        if gi:
            col_i += 1
        topics_df, coords = _topics_and_coords(path)
        for topic in topics:
            ax = fig.add_subplot(gs[0, col_i])
            col_i += 1
            ax.scatter(coords[:, 0], coords[:, 1], **bg_kwargs, zorder=1)
            ax.scatter(
                coords[:, 0],
                coords[:, 1],
                c=topics_df[f"topic_{topic}"].values,
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

    fig.savefig(save_path, bbox_inches="tight", dpi=300)
    plt.close(fig)

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
    [
        ("P22ATAC", INPUT_PATHS["P22ATAC"], [10, 37, 12, 45]),
        ("P22H3K27me3", INPUT_PATHS["P22H3K27me3"], [19, 2, 46, 48, 32, 1, 33, 3]),
    ],
    save_path=f"{OUT}/Figure_E5c.png",
)


# ═══════════════════════════════════════════════════════════════════════════
# Figure E5f
# ═══════════════════════════════════════════════════════════════════════════
plot_selected_topics(
    [
        ("P22ATAC", INPUT_PATHS["P22ATAC"], [31]),
        ("P22H3K27me3", INPUT_PATHS["P22H3K27me3"], [14, 45]),
    ],
    save_path=f"{OUT}/Figure_E5f.png",
)

# ═══════════════════════════════════════════════════════════════════════════
# Figure E5h
# ═══════════════════════════════════════════════════════════════════════════

DATASET_E5H = "H3K27me3"
SOURCES_E5H = ["nonspatial", "paper", "sweep"]

METRICS_E5H = [
    ("pas", "PAS", "#b322ad"),
    ("chaos", "CHAOS", "#1f47a3"),
    ("lisi", "LISI", "#be8046"),
]

RADII_E5H = [0.2, 0.4, 0.6, 0.8]
PAPER_RADIUS = 0.06
NONSPATIAL_RADIUS = 0.0

res_h = pd.read_csv(INPUT_PATHS["radius_sweep_metrics"], sep="\t")
d_h = (
    res_h[
        (res_h["dataset"] == DATASET_E5H) & (res_h["source"].isin(SOURCES_E5H))
    ]
    .sort_values("radius_num")
    .copy()
)

base_h = d_h[d_h["source"] == "nonspatial"].iloc[0]
for m, _, _ in METRICS_E5H:
    d_h[f"{m}_rel"] = d_h[m] / base_h[m]

xticks_h = [NONSPATIAL_RADIUS, PAPER_RADIUS] + RADII_E5H

fig, ax = plt.subplots(figsize=(9.0, 5))
for m, name, color in METRICS_E5H:
    ax.plot(
        d_h["radius_num"],
        d_h[f"{m}_rel"],
        marker="o",
        ms=6,
        lw=1.6,
        color=color,
        label=name,
    )

ax.axhline(1.0, color="0.4", lw=1.0, zorder=0)
ax.set_xlabel("Spatial radius (r)")
ax.set_ylabel("Normalized metric score")
ax.set_xticks(
    xticks_h,
    ["0" if v == NONSPATIAL_RADIUS else f"{v:g}" for v in xticks_h],
    fontsize=8,
)
ax.grid(alpha=0.25)
ax.set_axisbelow(True)
ax.legend(
    fontsize=9,
    frameon=False,
    loc="center left",
    bbox_to_anchor=(1.02, 0.5),
    borderaxespad=0,
)

fig.tight_layout()

fig.savefig(f"{OUT}/Figure_E5h.png", dpi=300, bbox_inches="tight")
plt.close(fig)

# ═══════════════════════════════════════════════════════════════════════════
# Figure E5i
# ═══════════════════════════════════════════════════════════════════════════

PANELS_E5I = [
    ("MTM vs MTM", "MTM"),
    ("r = 0.06", "sMTM r=0.06"),
    ("r = 0.8", "sMTM r=0.8"),
]
MODS_E5I = ["rna", "h3k27me3"]
MOD_LABEL_E5I = {"rna": "RNA", "h3k27me3": "H3K27me3"}
MOD_COLOR_E5I = {"rna": "#99000D", "h3k27me3": "#005A32"}
MOD_MARKER_E5I = {"rna": "o", "h3k27me3": "s"}
MOD_CMAP_E5I = {
    "rna": LinearSegmentedColormap.from_list(
        "rna_cmap",
        [
            (0.0, "#08306B"),
            (0.3, "#1562A9"),
            (0.45, "#FFFFFF"),
            (0.55, "#FFFFFF"),
            (0.7, "#DE2B25"),
            (1.0, "#99000D"),
        ],
    ),
    "h3k27me3": LinearSegmentedColormap.from_list(
        "h3k27me3_cmap",
        [
            (0.0, "#08306B"),
            (0.3, "#1562A9"),
            (0.5, "#FFFFFF"),
            (0.7, "#78C679"),
            (1.0, "#005A32"),
        ],
    ),
}
CBAR_LABEL_E5I = {"rna": "Pearson\ncorrelation\nRNA", "h3k27me3": "H3K27me3"}
AXIS_NOTE = "(Topic 1\u219250)"
ROW_LABEL = "MTM"

VARIANT_E5I = "centered"
METRIC_E5I = "full_pearson"

RADIUS_OF = {
    "r = 0.06": 0.06,
    "r = 0.2": 0.2,
    "r = 0.4": 0.4,
    "r = 0.6": 0.6,
    "r = 0.8": 0.8,
}

CORR_LIM = 1.0
PANEL_SIZE = 1.5
RIGHT_WIDTH = 3.2
SIM_BOX_ASPECT_E5I = 1.27

cache = np.load(INPUT_PATHS["corr_cache"])

sim = pd.read_csv(INPUT_PATHS["corr_similarity"])
sim = sim[sim["variant"] == VARIANT_E5I].copy()
sim["radius"] = sim["model"].map(RADIUS_OF)
sim = sim.dropna(subset=["radius"]).sort_values("radius")
radii = sorted(sim["radius"].unique())

n_rows, n_cols = len(PANELS_E5I), len(MODS_E5I)

fig = plt.figure(
    figsize=(
        PANEL_SIZE * n_cols + RIGHT_WIDTH + 1.0,
        PANEL_SIZE * n_rows + 1.0,
    )
)
outer = fig.add_gridspec(
    1, 2, width_ratios=[PANEL_SIZE * n_cols, RIGHT_WIDTH], wspace=0.35
)
left = outer[0].subgridspec(n_rows, n_cols, hspace=0.08, wspace=0.08)
right = outer[1].subgridspec(3, 1, height_ratios=[1.0, 0.35, 1.4])
cbars = right[0].subgridspec(
    1, 5, width_ratios=[1.0, 0.12, 1.1, 0.12, 0.6], wspace=0
)

for i, (label, title) in enumerate(PANELS_E5I):
    for j, mod in enumerate(MODS_E5I):
        labels = list(cache[f"{mod}__{VARIANT_E5I}__labels"])
        C = cache[f"{mod}__{VARIANT_E5I}__C"][labels.index(label)]

        ax = fig.add_subplot(left[i, j])
        ax.imshow(
            C,
            cmap=MOD_CMAP_E5I[mod],
            vmin=-CORR_LIM,
            vmax=CORR_LIM,
            aspect="equal",
            interpolation="nearest",
            rasterized=True,
        )
        ax.set_xticks([])
        ax.set_yticks([])
        for sp in ax.spines.values():
            sp.set_linewidth(0.6)
            sp.set_color("0.6")
        if i == 0:
            ax.set_title(MOD_LABEL_E5I[mod], fontsize=11)
        if j == 0:
            ax.set_ylabel(f"{title}\n{AXIS_NOTE}", fontsize=10)

ax_lab = fig.add_subplot(left[n_rows - 1, :], frameon=False)
ax_lab.set_xticks([])
ax_lab.set_yticks([])
ax_lab.patch.set_visible(False)
ax_lab.set_xlabel(f"{ROW_LABEL}\n{AXIS_NOTE}", fontsize=10)

for k, mod in enumerate(MODS_E5I):
    cax = fig.add_subplot(cbars[0, 1 + 2 * k])
    cb = fig.colorbar(
        ScalarMappable(Normalize(-CORR_LIM, CORR_LIM), MOD_CMAP_E5I[mod]),
        cax=cax,
        ticks=[-CORR_LIM, 0, CORR_LIM],
    )
    cb.ax.yaxis.set_major_formatter(FormatStrFormatter("%g"))
    cb.ax.tick_params(labelsize=9)
    cb.outline.set_linewidth(0.6)
    cb.set_label(CBAR_LABEL_E5I[mod], fontsize=10)
    cax.yaxis.set_label_position("left")

ax_sim = fig.add_subplot(right[2])
ax_sim.set_box_aspect(SIM_BOX_ASPECT_E5I)
for mod in MODS_E5I:
    d = sim[sim["modality"] == mod]
    ax_sim.plot(
        d["radius"],
        d[METRIC_E5I],
        marker=MOD_MARKER_E5I[mod],
        ms=5 if mod == "rna" else 4.4,
        lw=1.4,
        color=MOD_COLOR_E5I[mod],
        label=MOD_LABEL_E5I[mod],
    )

ax_sim.set_ylim(0.90, 0.98)
ax_sim.yaxis.set_major_locator(MultipleLocator(0.02))
ax_sim.yaxis.set_minor_locator(plt.NullLocator())
ax_sim.yaxis.set_major_formatter(FormatStrFormatter("%.2f"))
ax_sim.set_xticks(radii, [f"{r:g}" for r in radii])
ax_sim.set_xlim(min(radii) - 0.04, max(radii) + 0.04)
ax_sim.set_xlabel("Spatial radius (r)", fontsize=10)
ax_sim.set_ylabel("Matrix correlation", fontsize=10)
ax_sim.tick_params(labelsize=9)
ax_sim.grid(axis="y", alpha=0.25, lw=0.6)
ax_sim.set_axisbelow(True)
ax_sim.legend(
    fontsize=9,
    frameon=False,
    loc="lower center",
    bbox_to_anchor=(0.5, 1.0),
    ncol=len(MODS_E5I),
    columnspacing=1.5,
    borderaxespad=0.3,
)

fig.savefig(f"{OUT}/Figure_E5i.png", dpi=300, bbox_inches="tight")
plt.close(fig)

# ═══════════════════════════════════════════════════════════════════════════
# Figure E5j
# ═══════════════════════════════════════════════════════════════════════════

MODS_E5J = ["rna", "h3k27me3"]
MOD_LABEL_E5J = {"rna": "RNA", "h3k27me3": "H3K27me3"}
MOD_COLOR_E5J = {"rna": "#921F18", "h3k27me3": "#215C37"}
MOD_MARKER_E5J = {"rna": "o", "h3k27me3": "s"}
MOD_AREA_E5J = {"rna": 1.0, "h3k27me3": np.pi / 4}

KDE_GRID = 200
KDE_BW = 0.40
KDE_LEVELS = 5
KDE_LOW_CUT = 0.12
POINT_SIZE_E5J = 16
POINT_SIZE_TOP_E5J = 30
JITTER_E5J = 0.12
FONT_E5J = 13

BOX_EDGE_E5J = "black"
BOX_LW_E5J = 1.3
BOX_MEDIAN_LW_E5J = 1.6
BOX_WHIS_E5J = 1.5
BOX_WIDTH_E5J = 2.8 * JITTER_E5J

res_j = pd.read_csv(INPUT_PATHS["rbo_per_topic"])
res_j = res_j[res_j["used"] & res_j["matched"]]
norm_j = pd.read_csv(INPUT_PATHS["rbo_preserved"])

xx, yy = np.mgrid[0 : 1 : complex(0, KDE_GRID), 0 : 1 : complex(0, KDE_GRID)]
grid_j = np.vstack([xx.ravel(), yy.ravel()])

fig = plt.figure(figsize=(4.4, 6.72))
outer = fig.add_gridspec(2, 1, height_ratios=[4.0, 1.1], hspace=0.38)

ax = fig.add_subplot(outer[0])
handles = []
for mod in MODS_E5J:
    sub = res_j[res_j["modality"] == mod]
    x = sub["cosine_full"].to_numpy(float)
    y = sub["rbo_main"].to_numpy(float)
    ok = np.isfinite(x) & np.isfinite(y)
    if ok.sum() < 5:
        print(f"too few points for {mod}")
        continue
    col = to_rgb(MOD_COLOR_E5J[mod])

    kde = gaussian_kde(np.vstack([x[ok], y[ok]]), bw_method=KDE_BW)
    zz = kde(grid_j).reshape(xx.shape)
    ax.contourf(
        xx,
        yy,
        zz,
        levels=np.linspace(KDE_LOW_CUT * zz.max(), zz.max(), KDE_LEVELS + 1),
        cmap=LinearSegmentedColormap.from_list(
            f"kde_{mod}", [(*col, 0.18), (*col, 0.80)]
        ),
        zorder=2,
    )
    ax.scatter(
        x[ok],
        y[ok],
        s=POINT_SIZE_TOP_E5J * MOD_AREA_E5J[mod],
        marker=MOD_MARKER_E5J[mod],
        color=col,
        edgecolors="none",
        zorder=3,
    )
    handles.append(
        plt.Line2D(
            [],
            [],
            marker=MOD_MARKER_E5J[mod],
            ls="",
            color=col,
            markeredgecolor="none",
            markersize=14 * np.sqrt(MOD_AREA_E5J[mod]),
            label=MOD_LABEL_E5J[mod],
        )
    )

ax.plot([0, 1], [0, 1], color="0.55", ls="--", lw=1, zorder=1)
ax.set_xlim(0, 1)
ax.set_ylim(0, 1)
ax.set_aspect("equal", adjustable="box")
ax.xaxis.set_major_formatter(FormatStrFormatter("%g"))
ax.yaxis.set_major_formatter(FormatStrFormatter("%g"))
ax.set_xlabel("Cosine similarity (all features)", fontsize=FONT_E5J)
ax.set_ylabel("RBO (leading features)", fontsize=FONT_E5J)
ax.tick_params(labelsize=FONT_E5J)
ax.grid(alpha=0.25, lw=0.6)
ax.set_axisbelow(True)
ax.legend(
    handles=handles,
    frameon=False,
    fontsize=FONT_E5J,
    loc="lower center",
    bbox_to_anchor=(0.5, 1.0),
    ncol=len(handles),
    handletextpad=0.3,
    columnspacing=1.5,
    borderaxespad=0.3,
)

ax_d = fig.add_subplot(outer[1])
rng_j = np.random.default_rng(0)
for gi, mod in enumerate(MODS_E5J):
    v = norm_j[norm_j["modality"] == mod]["difference"].dropna().to_numpy(float)
    if not len(v):
        continue
    print(f"E5j bottom / {MOD_LABEL_E5J[mod]}: n = {len(v)} topics")
    ax_d.scatter(
        v,
        gi + rng_j.uniform(-JITTER_E5J, JITTER_E5J, len(v)),
        s=POINT_SIZE_E5J * MOD_AREA_E5J[mod],
        marker=MOD_MARKER_E5J[mod],
        color=MOD_COLOR_E5J[mod],
        edgecolors="none",
        zorder=2,
    )
    horiz = (
        dict(orientation="horizontal")
        if tuple(int(p) for p in plt.matplotlib.__version__.split(".")[:2])
        >= (3, 10)
        else dict(vert=False)
    )
    ax_d.boxplot(
        [v],
        positions=[gi],
        **horiz,
        widths=BOX_WIDTH_E5J,
        whis=BOX_WHIS_E5J,
        showfliers=False,
        patch_artist=True,
        manage_ticks=False,
        zorder=5,
        medianprops=dict(color=BOX_EDGE_E5J, lw=BOX_MEDIAN_LW_E5J),
        whiskerprops=dict(color=BOX_EDGE_E5J, lw=BOX_LW_E5J),
        capprops=dict(color=BOX_EDGE_E5J, lw=BOX_LW_E5J),
        boxprops=dict(facecolor="none", edgecolor=BOX_EDGE_E5J, lw=BOX_LW_E5J),
    )

ax_d.axvline(0, color="0.55", ls="--", lw=1, zorder=1)
ax_d.set_yticks(range(len(MODS_E5J)), [MOD_LABEL_E5J[m] for m in MODS_E5J])
ax_d.set_ylim(len(MODS_E5J) - 0.5, -0.5)
ax_d.xaxis.set_major_formatter(FormatStrFormatter("%g"))
ax_d.set_xlabel("Differential preservation", fontsize=FONT_E5J)
ax_d.tick_params(labelsize=FONT_E5J)
ax_d.grid(axis="x", alpha=0.25)

fig.savefig(f"{OUT}/Figure_E5j.png", dpi=300, bbox_inches="tight")
plt.close(fig)

# ═══════════════════════════════════════════════════════════════════════════
# Figure E5k
# ═══════════════════════════════════════════════════════════════════════════

MAPS_E5K = [
    ("theta_MTM", "MTM"),
    ("theta_r0p06", "sMTM r=0.06"),
    ("theta_r0p8", "sMTM r=0.8"),
]
TOPICS_E5K = ["topic_28", "topic_1"]

TARGET_E5K = "r = 0.06"
TOP_N_BARS = 20

MAP_PANEL_SIZE = 2.6
BAR_WIDTH = 4.2
SPOT_SIZE = 3
VMAX_E5K = 1.0
CBAR_LABEL = "Topic-spot distribution"

theta = {label: pd.read_csv(INPUT_PATHS[key]) for key, label in MAPS_E5K}
map_labels = [label for _, label in MAPS_E5K]

n_rows, n_cols = len(TOPICS_E5K), len(map_labels)

fig = plt.figure(
    figsize=(
        MAP_PANEL_SIZE * n_cols + BAR_WIDTH + 1.2,
        max(MAP_PANEL_SIZE * n_rows, 0.24 * TOP_N_BARS + 1.6),
    )
)
outer = fig.add_gridspec(
    1, 2, width_ratios=[MAP_PANEL_SIZE * n_cols, BAR_WIDTH], wspace=0.28
)
grid = outer[0].subgridspec(n_rows, n_cols, hspace=0.06, wspace=0.06)

map_axes, last_sc = [], None
for i, topic in enumerate(TOPICS_E5K):
    for j, label in enumerate(map_labels):
        df = theta[label]
        v = df[topic].to_numpy(float)

        ax = fig.add_subplot(grid[i, j])
        o = np.argsort(v)
        sc = ax.scatter(
            df["x"].to_numpy(float)[o],
            df["y"].to_numpy(float)[o],
            c=v[o],
            marker="s",
            s=SPOT_SIZE,
            cmap="gnuplot",
            vmin=0,
            vmax=VMAX_E5K,
            linewidths=0,
            rasterized=True,
        )
        ax.set_aspect("equal")
        ax.set_xticks([])
        ax.set_yticks([])
        for sp in ax.spines.values():
            sp.set_visible(False)
        if i == 0:
            ax.set_title(label, fontsize=10)
        if j == 0:
            ax.set_ylabel(topic.replace("topic_", "Topic "), fontsize=10)

        map_axes.append(ax)
        last_sc = sc

cb = fig.colorbar(
    last_sc, ax=map_axes, fraction=0.03, pad=0.02, location="right"
)
cb.set_label(CBAR_LABEL, fontsize=9)
cb.set_ticks([0, VMAX_E5K])
cb.ax.tick_params(labelsize=8)

shift = pd.read_csv(INPUT_PATHS["gini_moran_shift"])
tgt = shift[shift["model"] == TARGET_E5K].nlargest(TOP_N_BARS, "cosine_dist")
bar = tgt.sort_values("cosine_dist")
y = np.arange(len(bar))

ax_bar = fig.add_subplot(outer[1])
ax_bar.barh(
    y,
    bar["cosine_dist"],
    color="#31688e",
    edgecolor="0.3",
    linewidth=0.3,
    height=0.85,
)
ax_bar.set_yticks(y, bar["topic"], fontsize=8)
ax_bar.set_ylim(-0.5, len(bar) - 0.5)
ax_bar.set_xlim(0, float(bar["cosine_dist"].max()) * 1.02)
ax_bar.set_xlabel("Cosine distance")
ax_bar.set_ylabel("Topic")
ax_bar.grid(axis="x", alpha=0.25)
ax_bar.set_axisbelow(True)

fig.savefig(f"{OUT}/Figure_E5k.png", dpi=300, bbox_inches="tight")
plt.close(fig)