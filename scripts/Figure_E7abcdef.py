from matplotlib.collections import PatchCollection
from matplotlib.patches import Wedge
import matplotlib.pyplot as plt
import mudata as md
import mtopic as mt
import numpy as np
import pandas as pd
import scanpy as sc
import inspect
from scipy.optimize import linear_sum_assignment
import os

md.set_options(pull_on_update=False)

INPUT_PATHS = {
    "mtopic": "../data/HumanTonsil_trained.h5mu",
    "spatialglue": "../data/HumanTonsil_SpatialGlue.h5ad",
    "miso": "../data/HumanTonsil_MISO.h5ad",
    "coral": "../data/HumanTonsil_CORAL.h5ad",
    "cosmx_sc": "../data/C1_filtered_sc_K43_s1.h5mu",
    "cosmx_sq13": "../data/C1_filtered_sq13_K43_s1.h5mu",
    "cosmx_sq18": "../data/C1_filtered_sq18_K43_s1.h5mu",
    "cosmx_sq26": "../data/C1_filtered_sq26_K43_s1.h5mu",
    "cosmx_sq36": "../data/C1_filtered_sq36_K43_s1.h5mu",
    "cosmx_sq52": "../data/C1_filtered_sq52_K43_s1.h5mu",
    "benchmark_scores": "../data/benchmark_scores_ent400.csv",
    "benchmark_cloud": "../data/pred_vs_true_cloud_ent400_s1.npz",
    "benchmark_stats": "../data/pred_vs_true_stats_ent400.csv",
    "benchmark_per_topic": "../data/per_topic_r_ent400.csv",
}

OUT = "../figures/Figure_E7"
os.makedirs(OUT, exist_ok=True)

# ═══════════════════════════════════════════════════════════════════════════
# Topic annotations
# ═══════════════════════════════════════════════════════════════════════════

CORE_COLOR = {'Dark zone': '#0000cd', 'Light zone': '#ffff00', 'Mantle zone': '#65ff65', 'T cell zone': '#dea873', 'NK cells': '#834b14', 'Plasma cells': '#000000', 'Epithelial apical': '#ff69b4', 'Epithelial transitional': '#ff0000', 'Epithelial basal': '#8e24aa', 'Myeloid cells': '#00b4d8', 'MRC': '#f84866', 'Stromal': '#9e7d55', 'Endothelial': '#1a9fa1', 'Unknown': '#bebebe'}

N_SHADES = {'Dark zone': 3, 'Light zone': 6, 'Mantle zone': 5, 'T cell zone': 3, 'NK cells': 1, 'Plasma cells': 6, 'Epithelial apical': 4, 'Epithelial transitional': 2, 'Epithelial basal': 4, 'Myeloid cells': 9, 'MRC': 3, 'Stromal': 4, 'Endothelial': 2, 'Unknown': 2}

FAMILY_SHADES = {'Dark zone': ['#0000cd', '#0000f1', '#1515ff'], 'Light zone': ['#ffff00', '#ffff24', '#ffff47', '#ffff6b', '#ffff8f', '#ffffb3'], 'Mantle zone': ['#65ff65', '#89ff89', '#acffac', '#d0ffd0', '#d6ffd6'], 'T cell zone': ['#dea873', '#e5ba90', '#ecccad'], 'NK cells': ['#834b14'], 'Plasma cells': ['#000000', '#1a1a1a', '#333333', '#4d4d4d', '#666666', '#808080'], 'Epithelial apical': ['#ff69b4', '#ff8dc6', '#ffb0d8', '#ffd4ea'], 'Epithelial transitional': ['#ff0000', '#ff2424'], 'Epithelial basal': ['#8e24aa', '#a72ac7', '#b73fd6', '#c25ddd'], 'Myeloid cells': ['#00b4d8', '#00d2fc', '#20daff', '#44e0ff', '#68e6ff', '#8cecff', '#aff2ff', '#d3f8ff', '#d6f8ff'], 'MRC': ['#f84866', '#f96a83', '#fb8d9f'], 'Stromal': ['#9e7d55', '#ae8f69', '#bba080', '#c7b197'], 'Endothelial': ['#1a9fa1', '#1fbdc0'], 'Unknown': ['#bebebe', '#d0d0d0']}

TOPIC_CELLTYPE_sc = {
    "topic_1": "Epithelial apical-1",
    "topic_2": "Mantle zone-1",
    "topic_3": "Light zone-4",
    "topic_4": "Mantle zone-3",
    "topic_5": "Epithelial apical-4",
    "topic_6": "Light zone-6",
    "topic_7": "Endothelial-2",
    "topic_8": "Myeloid cells-3",
    "topic_9": "Myeloid cells-1",
    "topic_10": "MRC-1",
    "topic_11": "Myeloid cells-6",
    "topic_12": "Plasma cells-4",
    "topic_13": "Stromal-2",
    "topic_14": "Plasma cells-2",
    "topic_15": "T cell zone-1",
    "topic_16": "Light zone-2",
    "topic_17": "Dark zone-1",
    "topic_18": "Light zone-1",
    "topic_19": "Myeloid cells-5",
    "topic_20": "Endothelial-1",
    "topic_21": "Epithelial apical-3",
    "topic_22": "Myeloid cells-2",
    "topic_23": "Mantle zone-2",
    "topic_24": "Dark zone-2",
    "topic_25": "Dark zone-3",
    "topic_26": "Epithelial transitional-1",
    "topic_27": "Light zone-5",
    "topic_28": "Stromal-4",
    "topic_29": "Myeloid cells-7",
    "topic_30": "Unknown-1",
    "topic_31": "Light zone-3",
    "topic_32": "Stromal-3",
    "topic_33": "MRC-2",
    "topic_34": "Epithelial basal-2",
    "topic_35": "NK cells-1",
    "topic_36": "Epithelial basal-1",
    "topic_37": "Epithelial apical-2",
    "topic_38": "Myeloid cells-4",
    "topic_39": "T cell zone-2",
    "topic_40": "Plasma cells-5",
    "topic_41": "Plasma cells-1",
    "topic_42": "Stromal-1",
    "topic_43": "Plasma cells-3",
}

TOPIC_COLOR_sc = {
    "topic_1": "#ff69b4",
    "topic_2": "#65ff65",
    "topic_3": "#ffff6b",
    "topic_4": "#acffac",
    "topic_5": "#ffd4ea",
    "topic_6": "#ffffb3",
    "topic_7": "#1fbdc0",
    "topic_8": "#20daff",
    "topic_9": "#00b4d8",
    "topic_10": "#f84866",
    "topic_11": "#8cecff",
    "topic_12": "#4d4d4d",
    "topic_13": "#ae8f69",
    "topic_14": "#1a1a1a",
    "topic_15": "#dea873",
    "topic_16": "#ffff24",
    "topic_17": "#0000cd",
    "topic_18": "#ffff00",
    "topic_19": "#68e6ff",
    "topic_20": "#1a9fa1",
    "topic_21": "#ffb0d8",
    "topic_22": "#00d2fc",
    "topic_23": "#89ff89",
    "topic_24": "#0000f1",
    "topic_25": "#1515ff",
    "topic_26": "#ff0000",
    "topic_27": "#ffff8f",
    "topic_28": "#c7b197",
    "topic_29": "#aff2ff",
    "topic_30": "#bebebe",
    "topic_31": "#ffff47",
    "topic_32": "#bba080",
    "topic_33": "#f96a83",
    "topic_34": "#a72ac7",
    "topic_35": "#834b14",
    "topic_36": "#8e24aa",
    "topic_37": "#ff8dc6",
    "topic_38": "#44e0ff",
    "topic_39": "#e5ba90",
    "topic_40": "#666666",
    "topic_41": "#000000",
    "topic_42": "#9e7d55",
    "topic_43": "#333333",
}

TOPIC_CELLTYPE_sq13 = {
    "topic_1": "Epithelial apical-1",
    "topic_2": "Mantle zone-1",
    "topic_3": "Light zone-4",
    "topic_4": "Myeloid cells-1",
    "topic_5": "Epithelial apical-3",
    "topic_6": "Light zone-3",
    "topic_7": "Dark zone-2",
    "topic_8": "Light zone-2",
    "topic_9": "Myeloid cells-3",
    "topic_10": "MRC-1",
    "topic_11": "Myeloid cells-4",
    "topic_12": "Plasma cells-5",
    "topic_13": "Stromal-2",
    "topic_14": "Plasma cells-1",
    "topic_15": "T cell zone-1",
    "topic_16": "T cell zone-2",
    "topic_17": "Dark zone-1",
    "topic_18": "Light zone-1",
    "topic_19": "T cell zone-3",
    "topic_20": "Endothelial-1",
    "topic_21": "Epithelial apical-2",
    "topic_22": "Myeloid cells-2",
    "topic_23": "Mantle zone-2",
    "topic_24": "Dark zone-3",
    "topic_25": "Myeloid cells-7",
    "topic_26": "Epithelial transitional-1",
    "topic_27": "Light zone-5",
    "topic_28": "Plasma cells-3",
    "topic_29": "Myeloid cells-9",
    "topic_30": "MRC-3",
    "topic_31": "Epithelial basal-2",
    "topic_32": "Stromal-3",
    "topic_33": "MRC-2",
    "topic_34": "Epithelial basal-3",
    "topic_35": "NK cells-1",
    "topic_36": "Epithelial basal-1",
    "topic_37": "Myeloid cells-5",
    "topic_38": "Myeloid cells-6",
    "topic_39": "Myeloid cells-8",
    "topic_40": "Plasma cells-6",
    "topic_41": "Plasma cells-4",
    "topic_42": "Stromal-1",
    "topic_43": "Plasma cells-2",
}

TOPIC_COLOR_sq13 = {
    "topic_1": "#ff69b4",
    "topic_2": "#65ff65",
    "topic_3": "#ffff6b",
    "topic_4": "#00b4d8",
    "topic_5": "#ffb0d8",
    "topic_6": "#ffff47",
    "topic_7": "#0000f1",
    "topic_8": "#ffff24",
    "topic_9": "#20daff",
    "topic_10": "#f84866",
    "topic_11": "#44e0ff",
    "topic_12": "#666666",
    "topic_13": "#ae8f69",
    "topic_14": "#000000",
    "topic_15": "#dea873",
    "topic_16": "#e5ba90",
    "topic_17": "#0000cd",
    "topic_18": "#ffff00",
    "topic_19": "#ecccad",
    "topic_20": "#1a9fa1",
    "topic_21": "#ff8dc6",
    "topic_22": "#00d2fc",
    "topic_23": "#89ff89",
    "topic_24": "#1515ff",
    "topic_25": "#aff2ff",
    "topic_26": "#ff0000",
    "topic_27": "#ffff8f",
    "topic_28": "#333333",
    "topic_29": "#d6f8ff",
    "topic_30": "#fb8d9f",
    "topic_31": "#a72ac7",
    "topic_32": "#bba080",
    "topic_33": "#f96a83",
    "topic_34": "#b73fd6",
    "topic_35": "#834b14",
    "topic_36": "#8e24aa",
    "topic_37": "#68e6ff",
    "topic_38": "#8cecff",
    "topic_39": "#d3f8ff",
    "topic_40": "#808080",
    "topic_41": "#4d4d4d",
    "topic_42": "#9e7d55",
    "topic_43": "#1a1a1a",
}

TOPIC_CELLTYPE_sq18 = {
    "topic_1": "Epithelial apical-1",
    "topic_2": "Mantle zone-3",
    "topic_3": "Light zone-2",
    "topic_4": "Mantle zone-2",
    "topic_5": "Plasma cells-3",
    "topic_6": "Light zone-5",
    "topic_7": "Dark zone-3",
    "topic_8": "Light zone-3",
    "topic_9": "Myeloid cells-3",
    "topic_10": "Myeloid cells-4",
    "topic_11": "Myeloid cells-2",
    "topic_12": "Plasma cells-6",
    "topic_13": "Stromal-1",
    "topic_14": "Plasma cells-2",
    "topic_15": "T cell zone-1",
    "topic_16": "T cell zone-2",
    "topic_17": "Dark zone-1",
    "topic_18": "Light zone-1",
    "topic_19": "T cell zone-3",
    "topic_20": "Endothelial-1",
    "topic_21": "Epithelial apical-2",
    "topic_22": "Myeloid cells-1",
    "topic_23": "Mantle zone-5",
    "topic_24": "Dark zone-2",
    "topic_25": "Myeloid cells-5",
    "topic_26": "Epithelial transitional-1",
    "topic_27": "Light zone-4",
    "topic_28": "Plasma cells-5",
    "topic_29": "Myeloid cells-7",
    "topic_30": "MRC-1",
    "topic_31": "Epithelial basal-3",
    "topic_32": "Stromal-2",
    "topic_33": "MRC-2",
    "topic_34": "Epithelial basal-4",
    "topic_35": "NK cells-1",
    "topic_36": "Epithelial basal-1",
    "topic_37": "Mantle zone-4",
    "topic_38": "Epithelial basal-2",
    "topic_39": "Myeloid cells-6",
    "topic_40": "Plasma cells-4",
    "topic_41": "Epithelial transitional-2",
    "topic_42": "Mantle zone-1",
    "topic_43": "Plasma cells-1",
}

TOPIC_COLOR_sq18 = {
    "topic_1": "#ff69b4",
    "topic_2": "#acffac",
    "topic_3": "#ffff24",
    "topic_4": "#89ff89",
    "topic_5": "#333333",
    "topic_6": "#ffff8f",
    "topic_7": "#1515ff",
    "topic_8": "#ffff47",
    "topic_9": "#20daff",
    "topic_10": "#44e0ff",
    "topic_11": "#00d2fc",
    "topic_12": "#808080",
    "topic_13": "#9e7d55",
    "topic_14": "#1a1a1a",
    "topic_15": "#dea873",
    "topic_16": "#e5ba90",
    "topic_17": "#0000cd",
    "topic_18": "#ffff00",
    "topic_19": "#ecccad",
    "topic_20": "#1a9fa1",
    "topic_21": "#ff8dc6",
    "topic_22": "#00b4d8",
    "topic_23": "#d6ffd6",
    "topic_24": "#0000f1",
    "topic_25": "#68e6ff",
    "topic_26": "#ff0000",
    "topic_27": "#ffff6b",
    "topic_28": "#666666",
    "topic_29": "#aff2ff",
    "topic_30": "#f84866",
    "topic_31": "#b73fd6",
    "topic_32": "#ae8f69",
    "topic_33": "#f96a83",
    "topic_34": "#c25ddd",
    "topic_35": "#834b14",
    "topic_36": "#8e24aa",
    "topic_37": "#d0ffd0",
    "topic_38": "#a72ac7",
    "topic_39": "#8cecff",
    "topic_40": "#4d4d4d",
    "topic_41": "#ff2424",
    "topic_42": "#65ff65",
    "topic_43": "#000000",
}

TOPIC_CELLTYPE_sq26 = {
    "topic_1": "Epithelial apical-2",
    "topic_2": "Mantle zone-3",
    "topic_3": "Light zone-2",
    "topic_4": "Mantle zone-2",
    "topic_5": "Plasma cells-3",
    "topic_6": "Light zone-5",
    "topic_7": "Dark zone-1",
    "topic_8": "Light zone-3",
    "topic_9": "Myeloid cells-4",
    "topic_10": "Myeloid cells-2",
    "topic_11": "Unknown-1",
    "topic_12": "Myeloid cells-7",
    "topic_13": "Stromal-2",
    "topic_14": "Plasma cells-2",
    "topic_15": "Stromal-1",
    "topic_16": "T cell zone-1",
    "topic_17": "Dark zone-2",
    "topic_18": "Light zone-1",
    "topic_19": "NK cells-1",
    "topic_20": "Unknown-2",
    "topic_21": "Epithelial apical-1",
    "topic_22": "Myeloid cells-1",
    "topic_23": "Mantle zone-4",
    "topic_24": "Dark zone-3",
    "topic_25": "Myeloid cells-6",
    "topic_26": "Epithelial transitional-1",
    "topic_27": "Light zone-4",
    "topic_28": "Plasma cells-5",
    "topic_29": "Myeloid cells-8",
    "topic_30": "MRC-1",
    "topic_31": "Epithelial basal-2",
    "topic_32": "MRC-3",
    "topic_33": "MRC-2",
    "topic_34": "Epithelial basal-1",
    "topic_35": "Myeloid cells-3",
    "topic_36": "Endothelial-1",
    "topic_37": "Epithelial transitional-2",
    "topic_38": "Plasma cells-4",
    "topic_39": "Myeloid cells-5",
    "topic_40": "Mantle zone-5",
    "topic_41": "Epithelial basal-3",
    "topic_42": "Mantle zone-1",
    "topic_43": "Plasma cells-1",
}

TOPIC_COLOR_sq26 = {
    "topic_1": "#ff8dc6",
    "topic_2": "#acffac",
    "topic_3": "#ffff24",
    "topic_4": "#89ff89",
    "topic_5": "#333333",
    "topic_6": "#ffff8f",
    "topic_7": "#0000cd",
    "topic_8": "#ffff47",
    "topic_9": "#44e0ff",
    "topic_10": "#00d2fc",
    "topic_11": "#bebebe",
    "topic_12": "#aff2ff",
    "topic_13": "#ae8f69",
    "topic_14": "#1a1a1a",
    "topic_15": "#9e7d55",
    "topic_16": "#dea873",
    "topic_17": "#0000f1",
    "topic_18": "#ffff00",
    "topic_19": "#834b14",
    "topic_20": "#d0d0d0",
    "topic_21": "#ff69b4",
    "topic_22": "#00b4d8",
    "topic_23": "#d0ffd0",
    "topic_24": "#1515ff",
    "topic_25": "#8cecff",
    "topic_26": "#ff0000",
    "topic_27": "#ffff6b",
    "topic_28": "#666666",
    "topic_29": "#d3f8ff",
    "topic_30": "#f84866",
    "topic_31": "#a72ac7",
    "topic_32": "#fb8d9f",
    "topic_33": "#f96a83",
    "topic_34": "#8e24aa",
    "topic_35": "#20daff",
    "topic_36": "#1a9fa1",
    "topic_37": "#ff2424",
    "topic_38": "#4d4d4d",
    "topic_39": "#68e6ff",
    "topic_40": "#d6ffd6",
    "topic_41": "#b73fd6",
    "topic_42": "#65ff65",
    "topic_43": "#000000",
}

TOPIC_CELLTYPE_sq36 = {
    "topic_1": "Epithelial apical-2",
    "topic_2": "Mantle zone-3",
    "topic_3": "Light zone-2",
    "topic_4": "Mantle zone-2",
    "topic_5": "Plasma cells-3",
    "topic_6": "Light zone-4",
    "topic_7": "Dark zone-1",
    "topic_8": "Light zone-3",
    "topic_9": "Myeloid cells-4",
    "topic_10": "MRC-1",
    "topic_11": "Unknown-2",
    "topic_12": "Myeloid cells-7",
    "topic_13": "Stromal-1",
    "topic_14": "Plasma cells-2",
    "topic_15": "T cell zone-1",
    "topic_16": "T cell zone-2",
    "topic_17": "Dark zone-3",
    "topic_18": "Light zone-1",
    "topic_19": "Mantle zone-4",
    "topic_20": "Unknown-1",
    "topic_21": "Epithelial apical-3",
    "topic_22": "Myeloid cells-2",
    "topic_23": "Mantle zone-5",
    "topic_24": "Dark zone-2",
    "topic_25": "Myeloid cells-3",
    "topic_26": "Epithelial transitional-2",
    "topic_27": "Light zone-5",
    "topic_28": "Plasma cells-5",
    "topic_29": "T cell zone-3",
    "topic_30": "MRC-2",
    "topic_31": "Myeloid cells-1",
    "topic_32": "NK cells-1",
    "topic_33": "MRC-3",
    "topic_34": "Epithelial apical-1",
    "topic_35": "Myeloid cells-6",
    "topic_36": "Endothelial-1",
    "topic_37": "Epithelial transitional-1",
    "topic_38": "Plasma cells-4",
    "topic_39": "Myeloid cells-5",
    "topic_40": "Stromal-2",
    "topic_41": "Epithelial basal-1",
    "topic_42": "Mantle zone-1",
    "topic_43": "Plasma cells-1",
}

TOPIC_COLOR_sq36 = {
    "topic_1": "#ff8dc6",
    "topic_2": "#acffac",
    "topic_3": "#ffff24",
    "topic_4": "#89ff89",
    "topic_5": "#333333",
    "topic_6": "#ffff6b",
    "topic_7": "#0000cd",
    "topic_8": "#ffff47",
    "topic_9": "#44e0ff",
    "topic_10": "#f84866",
    "topic_11": "#d0d0d0",
    "topic_12": "#aff2ff",
    "topic_13": "#9e7d55",
    "topic_14": "#1a1a1a",
    "topic_15": "#dea873",
    "topic_16": "#e5ba90",
    "topic_17": "#1515ff",
    "topic_18": "#ffff00",
    "topic_19": "#d0ffd0",
    "topic_20": "#bebebe",
    "topic_21": "#ffb0d8",
    "topic_22": "#00d2fc",
    "topic_23": "#d6ffd6",
    "topic_24": "#0000f1",
    "topic_25": "#20daff",
    "topic_26": "#ff2424",
    "topic_27": "#ffff8f",
    "topic_28": "#666666",
    "topic_29": "#ecccad",
    "topic_30": "#f96a83",
    "topic_31": "#00b4d8",
    "topic_32": "#834b14",
    "topic_33": "#fb8d9f",
    "topic_34": "#ff69b4",
    "topic_35": "#8cecff",
    "topic_36": "#1a9fa1",
    "topic_37": "#ff0000",
    "topic_38": "#4d4d4d",
    "topic_39": "#68e6ff",
    "topic_40": "#ae8f69",
    "topic_41": "#8e24aa",
    "topic_42": "#65ff65",
    "topic_43": "#000000",
}

TOPIC_CELLTYPE_sq52 = {
    "topic_1": "Epithelial apical-1",
    "topic_2": "Mantle zone-2",
    "topic_3": "Light zone-1",
    "topic_4": "Mantle zone-3",
    "topic_5": "Plasma cells-3",
    "topic_6": "Light zone-5",
    "topic_7": "Dark zone-1",
    "topic_8": "Light zone-4",
    "topic_9": "MRC-3",
    "topic_10": "MRC-2",
    "topic_11": "Unknown-2",
    "topic_12": "Plasma cells-4",
    "topic_13": "Stromal-1",
    "topic_14": "Plasma cells-2",
    "topic_15": "T cell zone-1",
    "topic_16": "T cell zone-2",
    "topic_17": "Dark zone-3",
    "topic_18": "Light zone-2",
    "topic_19": "Mantle zone-4",
    "topic_20": "Unknown-1",
    "topic_21": "Myeloid cells-1",
    "topic_22": "Myeloid cells-4",
    "topic_23": "Mantle zone-5",
    "topic_24": "Myeloid cells-2",
    "topic_25": "Myeloid cells-6",
    "topic_26": "Epithelial transitional-1",
    "topic_27": "Light zone-6",
    "topic_28": "Light zone-3",
    "topic_29": "T cell zone-3",
    "topic_30": "MRC-1",
    "topic_31": "Myeloid cells-3",
    "topic_32": "Dark zone-2",
    "topic_33": "Stromal-3",
    "topic_34": "Epithelial apical-2",
    "topic_35": "Myeloid cells-5",
    "topic_36": "Endothelial-1",
    "topic_37": "Epithelial transitional-2",
    "topic_38": "Plasma cells-5",
    "topic_39": "Epithelial basal-2",
    "topic_40": "Stromal-2",
    "topic_41": "Epithelial basal-1",
    "topic_42": "Mantle zone-1",
    "topic_43": "Plasma cells-1",
}

TOPIC_COLOR_sq52 = {
    "topic_1": "#ff69b4",
    "topic_2": "#89ff89",
    "topic_3": "#ffff00",
    "topic_4": "#acffac",
    "topic_5": "#333333",
    "topic_6": "#ffff8f",
    "topic_7": "#0000cd",
    "topic_8": "#ffff6b",
    "topic_9": "#fb8d9f",
    "topic_10": "#f96a83",
    "topic_11": "#d0d0d0",
    "topic_12": "#4d4d4d",
    "topic_13": "#9e7d55",
    "topic_14": "#1a1a1a",
    "topic_15": "#dea873",
    "topic_16": "#e5ba90",
    "topic_17": "#1515ff",
    "topic_18": "#ffff24",
    "topic_19": "#d0ffd0",
    "topic_20": "#bebebe",
    "topic_21": "#00b4d8",
    "topic_22": "#44e0ff",
    "topic_23": "#d6ffd6",
    "topic_24": "#00d2fc",
    "topic_25": "#8cecff",
    "topic_26": "#ff0000",
    "topic_27": "#ffffb3",
    "topic_28": "#ffff47",
    "topic_29": "#ecccad",
    "topic_30": "#f84866",
    "topic_31": "#20daff",
    "topic_32": "#0000f1",
    "topic_33": "#bba080",
    "topic_34": "#ff8dc6",
    "topic_35": "#68e6ff",
    "topic_36": "#1a9fa1",
    "topic_37": "#ff2424",
    "topic_38": "#666666",
    "topic_39": "#a72ac7",
    "topic_40": "#ae8f69",
    "topic_41": "#8e24aa",
    "topic_42": "#65ff65",
    "topic_43": "#000000",
}

TOPIC_CELLTYPE = {
    "sc": TOPIC_CELLTYPE_sc,
    "sq13": TOPIC_CELLTYPE_sq13,
    "sq18": TOPIC_CELLTYPE_sq18,
    "sq26": TOPIC_CELLTYPE_sq26,
    "sq36": TOPIC_CELLTYPE_sq36,
    "sq52": TOPIC_CELLTYPE_sq52,
}

TOPIC_COLOR = {
    "sc": TOPIC_COLOR_sc,
    "sq13": TOPIC_COLOR_sq13,
    "sq18": TOPIC_COLOR_sq18,
    "sq26": TOPIC_COLOR_sq26,
    "sq36": TOPIC_COLOR_sq36,
    "sq52": TOPIC_COLOR_sq52,
}

# ═══════════════════════════════════════════════════════════════════════════
# Figure E7a
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
plt.savefig(f"{OUT}/Figure_E7a.png", bbox_inches="tight", dpi=300)
plt.close()

# ═══════════════════════════════════════════════════════════════════════════
# Figure E7b
# ═══════════════════════════════════════════════════════════════════════════
mdata_sc = md.read_h5mu(INPUT_PATHS["cosmx_sc"])
mt.pl.dominant_topics(
    mdata_sc,
    x="coords",
    s=1,
    legend=False,
    figsize=(5, 5),
    annotation=TOPIC_CELLTYPE["sc"],
    palette=TOPIC_COLOR["sc"],
    transparent=False,
    save=f"{OUT}/Figure_E7b.png",
)
del mdata_sc

# ═══════════════════════════════════════════════════════════════════════════
# Figure E7c
# ═══════════════════════════════════════════════════════════════════════════
ZOOM = dict(xrange=[0.53, 0.70], yrange=[0.35, 0.50])

RADIUS = {
    "sq13": 0.001,
    "sq18": 0.0015,
    "sq26": 0.0021,
    "sq36": 0.0029,
    "sq52": 0.0043,
}

scores = pd.read_csv(INPUT_PATHS["benchmark_scores"])
cells_per_spot = (
    scores[scores["seed"] == 1]
    .set_index("pitch")
    .eval("cells_covered / n_spots")
)

scatter_pie_params = inspect.signature(mt.pl.scatter_pie).parameters

for i, (level, radius_c) in enumerate(RADIUS.items(), start=1):
    mdata_p = md.read_h5mu(INPUT_PATHS[f"cosmx_{level}"])
    out = f"{OUT}/Figure_E7c_{i}.png"
    pitch = int(level.replace("sq", ""))
    cps = f"{cells_per_spot.get(pitch, float('nan')):.1f}".rstrip("0").rstrip(".")
    title = f"Pooled resolution: {pitch} $\\mu$m"

    kwargs = dict(
        radius=radius_c,
        legend=False,
        **ZOOM,
        figsize=(5, 5),
        annotation=TOPIC_CELLTYPE[level],
        palette=TOPIC_COLOR[level],
        transparent=False,
    )

    if "title" in scatter_pie_params:
        mt.pl.scatter_pie(mdata_p, title=title, save=out, **kwargs)
    else:
        mt.pl.scatter_pie(mdata_p, **kwargs)
        fig = plt.gcf()
        fig.suptitle(title, fontsize=13)
        fig.savefig(out, dpi=300, bbox_inches="tight")
        plt.close(fig)

    del mdata_p

# ═══════════════════════════════════════════════════════════════════════════
# Figure E7d
# ═══════════════════════════════════════════════════════════════════════════
PITCHES = (13, 18, 26, 36, 52)
HEX = 55

stats = pd.read_csv(INPUT_PATHS["benchmark_stats"])
agg = stats.drop(columns="seed").groupby("pitch").agg(["mean", "std"])

z = np.load(INPUT_PATHS["benchmark_cloud"], allow_pickle=True)
cloud = {p: (z[f"t_{p}"], z[f"q_{p}"]) for p in PITCHES if f"t_{p}" in z}

fig, axes = plt.subplots(
    1, len(PITCHES), figsize=(3.4 * len(PITCHES), 3.9), sharex=True, sharey=True
)

for ax, p in zip(np.atleast_1d(axes), PITCHES):
    t, q = cloud[p]
    hb = ax.hexbin(
        t, q, gridsize=HEX, bins="log", cmap="magma_r", mincnt=1,
        extent=(0, 1, 0, 1),
    )
    ax.plot([0, 1], [0, 1], ":", color="tab:blue", lw=1)
    ax.text(
        0.04,
        0.96,
        f"r={agg[('r', 'mean')][p]:.3f}±{agg[('r', 'std')][p]:.3f}\n"
        f"RMSE={agg[('rmse', 'mean')][p]:.3f}±{agg[('rmse', 'std')][p]:.3f}",
        transform=ax.transAxes,
        va="top",
        ha="left",
        fontsize=9,
    )
    ax.set_title(f"{p} $\\mu$m", fontsize=11)
    ax.set_xlabel("Reference topic proportion")
    ax.set_aspect("equal")

axes[0].set_ylabel("Recovered topic proportion")
cb = fig.colorbar(hb, ax=axes, fraction=0.015)
cb.set_label("Spot-topic pairs")

fig.savefig(f"{OUT}/Figure_E7d.png", dpi=300, bbox_inches="tight")
plt.close(fig)

# ═══════════════════════════════════════════════════════════════════════════
# Figure E7e
# ═══════════════════════════════════════════════════════════════════════════
TREND_FRAC = 0.6


def trend(x, y, frac=TREND_FRAC, n=120):
    """Smooth recovery-vs-usage curve in log-usage space.

    LOWESS if statsmodels is available, otherwise a quadratic fit - the shape is
    what matters here, not the estimator."""
    lx = np.log10(x)
    grid = np.linspace(lx.min(), lx.max(), n)
    try:
        from statsmodels.nonparametric.smoothers_lowess import lowess

        sm = lowess(y, lx, frac=frac, xvals=grid)
    except Exception:
        sm = np.polyval(np.polyfit(lx, y, 2), grid)
    return 10 ** grid, sm


per = pd.read_csv(INPUT_PATHS["benchmark_per_topic"])

fig, ax = plt.subplots(figsize=(8.6, 5))

for p, c in zip(PITCHES, [f"C{i}" for i in range(len(PITCHES))]):
    s = per[per["pitch"] == p]
    ax.scatter(s["usage"], s["r"], s=40, color=c, linewidths=0,
               label=f"{p} $\\mu$m")

gx, gy = trend(per["usage"].to_numpy(), per["r"].to_numpy())
ax.plot(gx, gy, color="white", lw=7, solid_capstyle="round", zorder=4)
ax.plot(gx, gy, color="black", lw=4, solid_capstyle="round", zorder=5)

ax.axhline(0, color="grey", lw=0.8)
ax.set_xscale("log")
ax.set_xlabel("Mean topic proportion across cells")
ax.set_ylabel("Per-topic recovery (Pearson $r$)")
ax.set_ylim(0, 1.05)
ax.legend(
    fontsize=8,
    title="Pooling\nresolution",
    title_fontsize=8,
    loc="upper left",
    bbox_to_anchor=(1.02, 1),
    borderaxespad=0,
    frameon=False,
)

fig.tight_layout()

fig.savefig(f"{OUT}/Figure_E7e.png", dpi=300, bbox_inches="tight")
plt.close(fig)

# ═══════════════════════════════════════════════════════════════════════════
# Figure E7f
# ═══════════════════════════════════════════════════════════════════════════
PANELS = (
    ("cos", "Cosine similarity", "tab:blue"),
    ("top1", "Dominant topic agreement", "tab:green"),
)

NULLS = (
    ("perm", "Shuffled cell-spot assignment", ":"),
    ("permW", "Shuffled gene labels", "--"),
)

TRUTH = "soft"

num = scores.select_dtypes("number").drop(columns="seed")
mean = num.groupby("pitch").mean()
sd = num.groupby("pitch").std()
pitches = list(mean.index)

fig, axes = plt.subplots(1, 2, figsize=(11, 4))

for ax, (key, name, col) in zip(axes, PANELS):
    m, s = mean[f"{key}_{TRUTH}"], sd[f"{key}_{TRUTH}"]
    ax.plot(pitches, m, "o-", color=col, ms=4, label="mTopic")
    ax.fill_between(pitches, m - s, m + s, color=col, alpha=0.25)

    lo = m.min() - s.max()
    for suffix, null_label, ls in NULLS:
        y = mean[f"{key}_{TRUTH}_{suffix}"]
        ax.plot(pitches, y, ls, color=col, lw=1.4, label=null_label)
        lo = min(lo, y.min())

    ax.set_ylim(max(0.0, np.floor((lo - 0.05) * 10) / 10), 1.0)
    ax.set_xticks(pitches, [str(p) for p in pitches])
    ax.set_xlabel("Pooling resolution ($\\mu$m)")
    ax.set_ylabel(name)

handles, labels = axes[0].get_legend_handles_labels()
fig.legend(
    handles,
    labels,
    fontsize=9,
    frameon=False,
    ncol=3,
    loc="upper center",
    bbox_to_anchor=(0.5, 1.06),
)

fig.tight_layout()

fig.savefig(f"{OUT}/Figure_E7f.pdf", dpi=300, bbox_inches="tight")
plt.close(fig)
