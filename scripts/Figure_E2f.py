import mtopic
import os

INPUT_PATHS = {
    "trained_model": "../data/P22MouseBrainATAC_trained.h5mu",
}

OUT = "../figures/Figure_E2"
os.makedirs(OUT, exist_ok=True)


# ═══════════════════════════════════════════════════════════════════════════
# Figure E2f
# ═══════════════════════════════════════════════════════════════════════════
GRAYS = {
    "topic_1": "#e7e7e7",
    "topic_2": "#a4a4a4",
    "topic_3": "#999999",
    "topic_4": "#b9b9b9",
    "topic_5": "#b5b5b5",
    "topic_6": "#b2b2b2",
    "topic_7": "#a7a7a7",
    "topic_8": "#a3a3a3",
    "topic_9": "#ececec",
    "topic_10": "#dbdbdb",
    "topic_11": "#a1a1a1",
    "topic_12": "#e1e1e1",
    "topic_13": "#cccccc",
    "topic_14": "#9a9a9a",
    "topic_15": "#f2f2f2",
    "topic_16": "#eaeaea",
    "topic_17": "#b1b1b1",
    "topic_18": "#b3b3b3",
    "topic_19": "#d6d6d6",
    "topic_20": "#e6e6e6",
    "topic_21": "#dddddd",
    "topic_22": "#afafaf",
    "topic_23": "#ebebeb",
    "topic_24": "#cbcbcb",
    "topic_25": "#efefef",
    "topic_26": "#cfcfcf",
    "topic_27": "#f1f1f1",
    "topic_28": "#969696",
    "topic_29": "#aaaaaa",
    "topic_30": "#e8e8e8",
    "topic_31": "#c1c1c1",
    "topic_32": "#dadada",
    "topic_33": "#9f9f9f",
    "topic_34": "#ededed",
    "topic_35": "#c6c6c6",
    "topic_36": "#ababab",
    "topic_37": "#9c9c9c",
    "topic_38": "#9b9b9b",
    "topic_39": "#aeaeae",
    "topic_40": "#d0d0d0",
    "topic_41": "#acacac",
    "topic_42": "#d7d7d7",
    "topic_43": "#bcbcbc",
    "topic_44": "#a6a6a6",
    "topic_45": "#989898",
    "topic_46": "#c4c4c4",
    "topic_47": "#e3e3e3",
    "topic_48": "#b8b8b8",
    "topic_49": "#9d9d9d",
    "topic_50": "#cecece",
}

mdata = mtopic.read.h5mu(INPUT_PATHS["trained_model"])
palette = GRAYS.copy()
palette["topic_21"] = "#ffd500"
palette["topic_43"] = "#5c0099"

mtopic.pl.scatter_pie(
    mdata,
    radius=0.005,
    figsize=(10, 5),
    legend=True,
    legend_ncol=3,
    xrange=[0.11, 0.34],
    yrange=[0.22, 0.69],
    palette=palette,
    annotation=mdata.uns["TOPIC_CELLTYPE"],
    save=f"{OUT}/Figure_E2f_1.png",
)

mtopic.pl.scatter_pie(
    mdata,
    radius=0.005,
    figsize=(10, 5),
    legend=True,
    legend_ncol=3,
    xrange=[0.14, 0.24],
    yrange=[0.48, 0.61],
    palette=palette,
    annotation=mdata.uns["TOPIC_CELLTYPE"],
    save=f"{OUT}/Figure_E2f_2.png",
)

mtopic.pl.scatter_pie(
    mdata,
    radius=0.005,
    figsize=(8, 4),
    legend=True,
    legend_ncol=3,
    fontsize=6,
    legend_markersize=6,
    palette=palette,
    annotation=mdata.uns["TOPIC_CELLTYPE"],
    save=f"{OUT}/Figure_E2f_3.png",
)
