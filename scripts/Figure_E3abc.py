import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import mtopic
import numpy as np
import pandas as pd
import os

INPUT_PATHS = {
    "P22H3K4me3": "../data/P22MouseBrainH3K4me3_trained.h5mu",
    "P22H3K27me3": "../data/P22MouseBrainH3K27me3_trained.h5mu",
    "P22H3K27ac": "../data/P22MouseBrainH3K27ac_trained.h5mu",
}

OUT = "../figures/Figure_E3"
os.makedirs(OUT, exist_ok=True)


DEFAULT_PALETTE = {
    "topic_1": "#1f77b4",
    "topic_2": "#aec7e8",
    "topic_3": "#ff7f0e",
    "topic_4": "#ffbb78",
    "topic_5": "#2ca02c",
    "topic_6": "#98df8a",
    "topic_7": "#d62728",
    "topic_8": "#ff9896",
    "topic_9": "#9467bd",
    "topic_10": "#c5b0d5",
    "topic_11": "#8c564b",
    "topic_12": "#c49c94",
    "topic_13": "#e377c2",
    "topic_14": "#f7b6d2",
    "topic_15": "#7f7f7f",
    "topic_16": "#c7c7c7",
    "topic_17": "#bcbd22",
    "topic_18": "#dbdb8d",
    "topic_19": "#17becf",
    "topic_20": "#9edae5",
    "topic_21": "#393b79",
    "topic_22": "#5254a3",
    "topic_23": "#6b6ecf",
    "topic_24": "#9c9ede",
    "topic_25": "#637939",
    "topic_26": "#8ca252",
    "topic_27": "#b5cf6b",
    "topic_28": "#cedb9c",
    "topic_29": "#8c6d31",
    "topic_30": "#bd9e39",
    "topic_31": "#e7ba52",
    "topic_32": "#e7cb94",
    "topic_33": "#843c39",
    "topic_34": "#ad494a",
    "topic_35": "#d6616b",
    "topic_36": "#e7969c",
    "topic_37": "#7b4173",
    "topic_38": "#a55194",
    "topic_39": "#ce6dbd",
    "topic_40": "#de9ed6",
    "topic_41": "#3182bd",
    "topic_42": "#6baed6",
    "topic_43": "#9ecae1",
    "topic_44": "#c6dbef",
    "topic_45": "#e6550d",
    "topic_46": "#fd8d3c",
    "topic_47": "#fdae6b",
    "topic_48": "#fdd0a2",
    "topic_49": "#31a354",
    "topic_50": "#74c476",
}


def make_panel(sample, panel_label):
    mdata = mtopic.read.h5mu(INPUT_PATHS[sample])

    mtopic.pl.dominant_topics(
        mdata,
        x="coords",
        s=28,
        figsize=(8, 4),
        legend=True,
        markerscale=3,
        legend_ncol=3,
        palette=DEFAULT_PALETTE,
        save=f"{OUT}/Figure_{panel_label}_1.png",
    )

    palette = DEFAULT_PALETTE
    figsize = (18, 8)

    topic_order = list(mdata.obsm["topics"].columns)
    gamma = mdata.obsm["topics"][topic_order]

    dominant_topic = gamma.idxmax(axis=1)
    dominant_proportion = gamma.max(axis=1)
    topic_rank = {t: i for i, t in enumerate(topic_order)}

    cell_order = (
        pd.DataFrame(
            {
                "topic": dominant_topic,
                "topic_rank": dominant_topic.map(topic_rank),
                "proportion": dominant_proportion,
            }
        )
        .sort_values(["topic_rank", "proportion"], ascending=[True, False])
        .index
    )
    gamma_sorted = gamma.loc[cell_order].T.loc[topic_order]

    N = gamma_sorted.shape[1]
    K = gamma_sorted.shape[0]

    sorted_topics = dominant_topic.loc[cell_order].values
    color_bar = np.array(
        [plt.matplotlib.colors.to_rgb(palette.get(t, "#bebebe")) for t in sorted_topics]
    )

    fig, (ax_bar, ax_heat) = plt.subplots(
        2,
        1,
        figsize=figsize,
        gridspec_kw={"height_ratios": [0.04, 1], "hspace": 0.01},
    )
    ax_bar.imshow(
        color_bar[np.newaxis, :, :],
        aspect="auto",
        interpolation="none",
        extent=[0, N, 0, 1],
    )
    ax_bar.set_xlim(0, N)
    ax_bar.axis("off")

    im = ax_heat.imshow(
        gamma_sorted.values,
        aspect="auto",
        cmap="gnuplot",
        vmin=0,
        vmax=1,
        interpolation="none",
        extent=[0, N, K, 0],
    )
    ax_heat.set_xlim(0, N)
    ax_heat.set_xticks([])
    ax_heat.set_yticks([])

    handles = [mpatches.Patch(color=palette[t], label=t) for t in topic_order]
    fig.legend(
        handles=handles,
        loc="lower center",
        ncol=min(10, len(topic_order)),
        fontsize=7,
        frameon=False,
        bbox_to_anchor=(0.5, -0.02),
    )
    plt.savefig(f"{OUT}/Figure_{panel_label}_2.png", bbox_inches="tight", dpi=300)
    plt.close(fig)


# ═══════════════════════════════════════════════════════════════════════════
# Figure E3a
# ═══════════════════════════════════════════════════════════════════════════
make_panel("P22H3K4me3", "E3a")

# ═══════════════════════════════════════════════════════════════════════════
# Figure E3b
# ═══════════════════════════════════════════════════════════════════════════
make_panel("P22H3K27me3", "E3b")

# ═══════════════════════════════════════════════════════════════════════════
# Figure E3c
# ═══════════════════════════════════════════════════════════════════════════
make_panel("P22H3K27ac", "E3c")
