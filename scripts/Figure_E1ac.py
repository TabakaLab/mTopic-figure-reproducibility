import json

import matplotlib.pyplot as plt
from matplotlib.ticker import NullFormatter
import numpy as np
import pandas as pd
import os

INPUT_PATHS = {
    "cpu_10k": "../data/Scalability_CPU_10k.csv",
    "gpu_10k": "../data/Scalability_GPU_10k.csv",
    "cpu_100k": "../data/Scalability_CPU_100k.csv",
    "gpu_100k": "../data/Scalability_GPU_100k.csv",
    "tsll_optimal_params": "../data/Topic_selection_curves.json",
    "tsll_p22_atac": "../data/P22MouseBrainATAC_TSLL_CV.csv",
    "tsll_p22_h3k4me3": "../data/P22MouseBrainH3K4me3_TSLL_CV.csv",
    "tsll_p22_h3k27ac": "../data/P22MouseBrainH3K27ac_TSLL_CV.csv",
    "tsll_p22_h3k27me3": "../data/P22MouseBrainH3K27me3_TSLL_CV.csv",
    "tsll_human_tonsil": "../data/HumanTonsil_TSLL_CV.csv",
    "tsll_pbmc": "../data/HumanPBMC_TSLL_CV.csv",
}

OUT = "../figures/Figure_E1"
os.makedirs(OUT, exist_ok=True)


# ═══════════════════════════════════════════════════════════════════════════
# Figure E1a
# ═══════════════════════════════════════════════════════════════════════════
mod_col = {2: "#1f77b4", 3: "#ff7f0e", 4: "#2ca02c"}


def plot_runtime(cpu_path, gpu_path, out_path, title):
    cpu = pd.read_csv(cpu_path)
    gpu = pd.read_csv(gpu_path)

    fig, (ax, ax2) = plt.subplots(
        2,
        1,
        figsize=(6.5, 9),
        gridspec_kw={"height_ratios": [2.5, 1]},
        sharex=True,
    )

    cell_counts = sorted(cpu["n_cells"].unique())
    xtick_labels = [f"{c // 1000}k" if c < 1_000_000 else "1M" for c in cell_counts]

    for n_mod in sorted(cpu["n_mod"].unique()):
        c = mod_col.get(n_mod, "gray")
        sub_cpu = cpu[cpu["n_mod"] == n_mod].sort_values("n_cells")
        sub_gpu = gpu[gpu["n_mod"] == n_mod].sort_values("n_cells")
        merged = sub_cpu.merge(sub_gpu, on="n_cells", suffixes=("_cpu", "_gpu"))
        speedup = merged["time_cpu"] / merged["time_gpu"]
        ax2.plot(
            merged["n_cells"],
            speedup,
            "-o",
            color=c,
            lw=2.5,
            ms=11,
            label=f"{n_mod} mod",
        )

    ax2.axhline(1, color="black", lw=0.8, ls=":", alpha=0.6)
    ax2.set_ylabel("Runtime ratio\n(CPU / GPU)")
    fig.suptitle(title)
    ax2.grid(axis="y", which="major", linestyle=":", linewidth=1.5, alpha=0.5)
    ax2.legend(fontsize=9, loc="upper right", frameon=False)
    ax2.spines[["right", "top"]].set_visible(False)
    ax2.set_xlabel("Number of cells (log scale)")
    ax2.set_yticks([1, 2, 3, 4])
    ax2.set_yticklabels(["1", "2", "3", "4"])
    ax2.set_ylim(1, 4.3)

    for n_mod in sorted(cpu["n_mod"].unique()):
        c = mod_col.get(n_mod, "gray")
        sub_cpu = cpu[cpu["n_mod"] == n_mod].sort_values("n_cells")
        sub_gpu = gpu[gpu["n_mod"] == n_mod].sort_values("n_cells")
        ax.plot(
            sub_cpu["n_cells"],
            sub_cpu["time"],
            "-o",
            color=c,
            lw=2.5,
            ms=11,
            label=f"CPU, {n_mod} mod",
        )
        ax.plot(
            sub_gpu["n_cells"],
            sub_gpu["time"],
            "--o",
            color=c,
            lw=1.8,
            ms=10,
            mfc="white",
            mew=1.8,
            label=f"GPU, {n_mod} mod",
        )

    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_ylabel("Runtime (log scale)")
    ax.set_xticks(cell_counts)
    ax.set_xticklabels(xtick_labels)
    ax.xaxis.set_minor_formatter(NullFormatter())

    y_major = [10, 30, 60, 300, 600, 1800, 3600, 6300]
    y_labels = ["10s", "30s", "1min", "5min", "10min", "30min", "1h", "1h45min"]
    ax.set_yticks(y_major)
    ax.set_yticklabels(y_labels)
    ax.yaxis.set_minor_locator(plt.NullLocator())

    ax.set_xlim(cell_counts[0] * 0.8, cell_counts[-1] * 1.2)
    ax.grid(axis="y", which="major", linestyle=":", linewidth=1.5, alpha=0.5)
    ax.legend(fontsize=8, loc="lower right", frameon=False)
    ax.spines[["right", "top"]].set_visible(False)
    ax.set_ylim(5, 7200)

    for a in (ax, ax2):
        for spine in a.spines.values():
            spine.set_linewidth(1.2)

    plt.tight_layout()
    plt.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close()


plot_runtime(
    INPUT_PATHS["cpu_10k"],
    INPUT_PATHS["gpu_10k"],
    f"{OUT}/Figure_E1a_1.png",
    "10k features",
)
plot_runtime(
    INPUT_PATHS["cpu_100k"],
    INPUT_PATHS["gpu_100k"],
    f"{OUT}/Figure_E1a_2.png",
    "100k features",
)

# ═══════════════════════════════════════════════════════════════════════════
# Figure E1c
# ═══════════════════════════════════════════════════════════════════════════
with open(INPUT_PATHS["tsll_optimal_params"]) as f:
    saved_params = json.load(f)

p22_datasets = ["P22_ATAC", "P22_H3K4me3", "P22_H3K27ac", "P22_H3K27me3"]
single_datasets = ["Human_tonsil", "PBMC"]

x = np.array(
    [5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 90, 95, 100]
)

model_name = "Rational (a+bx)/(1+cx+dx²)"


def model_fn(x_, p):
    return (p[0] + p[1] * x_) / (1 + p[2] * x_ + p[3] * x_**2)


x_fit_norm = np.linspace(1e-6, 1, 300)
x_fit = x_fit_norm * 95 + 5


def norm01(arr):
    return (arr - arr.min()) / (arr.max() - arr.min())


def find_plateau(x_fit, y_fit, threshold_fraction=0.05):
    dy = np.gradient(y_fit, x_fit)
    threshold = threshold_fraction * dy.max()
    below = np.where(dy <= threshold)[0]
    if len(below) == 0:
        return x_fit[-1]
    return x_fit[below[0]]


def get_curve(dataset):
    meta = saved_params[dataset][model_name]
    params = meta["params"]
    y_fit = model_fn(x_fit_norm, params) * meta["y_range"] + meta["y_min"]
    df = pd.read_csv(INPUT_PATHS[f"tsll_{dataset.lower()}"])
    df = df[df.CV != "total"]
    y_data = (
        df.groupby("K")["total"].mean().reset_index().sort_values("K")["total"].values
    )
    return norm01(y_fit), norm01(y_data)


threshold_fraction = 0.05

all_fit, all_data, p22_K = [], [], []
for ds in p22_datasets:
    yf, yd = get_curve(ds)
    all_fit.append(yf)
    all_data.append(yd)
    p22_K.append(find_plateau(x_fit, yf, threshold_fraction))

avg_fit = np.array(all_fit).mean(axis=0)
avg_data = np.array(all_data).mean(axis=0)
plateau_x_p22 = float(np.mean(p22_K))
plateau_x_p22_std = float(np.std(p22_K, ddof=1))

panels = [("P22 mouse brain", avg_fit, avg_data, plateau_x_p22, plateau_x_p22_std)]
for ds in single_datasets:
    yf, yd = get_curve(ds)
    plateau_x = find_plateau(x_fit, yf, threshold_fraction)
    panels.append((ds.replace("_", " "), yf, yd, plateau_x, None))

for i, (title, y_fit, y_data, plateau_x, plateau_x_std) in enumerate(panels, start=1):
    plateau_y = np.interp(plateau_x, x_fit, y_fit)
    fig, ax = plt.subplots(figsize=(5.5, 5))
    ax.scatter(x, y_data, color="black", s=25, zorder=5, label="Data")
    ax.plot(x_fit, y_fit, color="black", linewidth=2, label="Fit")
    ax.axvline(plateau_x, color="black", linestyle=":", linewidth=1.4, zorder=4)
    if plateau_x_std is not None:
        ax.axvspan(
            plateau_x - plateau_x_std,
            plateau_x + plateau_x_std,
            alpha=0.1,
            color="black",
            zorder=3,
        )
    ax.scatter(
        [plateau_x], [plateau_y], color="white", s=100, zorder=6, edgecolor="black"
    )
    ax.text(
        plateau_x + 1.5, plateau_y - 0.1, f"K={plateau_x:.0f}", fontsize=9, va="center"
    )
    ax.set_title(title, fontsize=11)
    ax.set_xticks([5, 25, 50, 75, 100])
    ax.set_xlabel("Number of topics")
    ax.set_ylabel("Normalized log-likelihood")
    ax.legend(fontsize=8, frameon=False)
    fig.tight_layout()
    fig.savefig(f"{OUT}/Figure_E1c_{i}.png", bbox_inches="tight")
    plt.close()
