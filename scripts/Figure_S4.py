from matplotlib.gridspec import GridSpec
import matplotlib.pyplot as plt
import mtopic
import muon as mu
import numpy as np
import os

mu.set_options(pull_on_update=False)

INPUT_PATHS = {
    "pbmc_trained": "../data/HumanPBMC_trained.h5mu",
}
OUT = "../figures/Figure_S4"
os.makedirs(OUT, exist_ok=True)


# ═══════════════════════════════════════════════════════════════════════════
# Figure S4
# ═══════════════════════════════════════════════════════════════════════════
mdata = mtopic.read.h5mu(INPUT_PATHS["pbmc_trained"])

x_key = "umap"
topics_key = "topics"
cmap = "gnuplot"
marker = "o"
s_size = 1
fontsize = 12

n_topics = mdata.obsm[topics_key].shape[1]
nrow, ncol = 4, 8

fig = plt.figure(constrained_layout=True, figsize=(ncol * 2, nrow * 2.2))
gs = GridSpec(nrow, ncol, figure=fig)

for t in range(n_topics):
    ax = fig.add_subplot(gs[t // ncol, t % ncol])
    thr = np.percentile(mdata.obsm[topics_key].iloc[:, t], 97)
    vmax = np.percentile(mdata.obsm[topics_key].iloc[:, t], 99.9)
    mask = mdata.obsm[topics_key].iloc[:, t] >= thr

    ax.scatter(
        x=mdata.obsm[x_key].values[~mask, 0],
        y=mdata.obsm[x_key].values[~mask, 1],
        edgecolor="none",
        c=mdata.obsm[topics_key].values[~mask, t],
        cmap=cmap,
        s=s_size,
        vmin=0,
        vmax=vmax,
        marker=marker,
    )
    p = ax.scatter(
        x=mdata.obsm[x_key].values[mask, 0],
        y=mdata.obsm[x_key].values[mask, 1],
        edgecolor="none",
        c=mdata.obsm[topics_key].values[mask, t],
        cmap=cmap,
        s=s_size,
        vmin=0,
        vmax=vmax,
        marker=marker,
    )

    ax.set(aspect="equal")
    ax.title.set_size(fontsize)
    ax.set_title(f"Topic {t+1}")
    ax.set_xticks([])
    ax.set_yticks([])

plt.savefig(f"{OUT}/Figure_S4.png", dpi=300)
plt.close()
