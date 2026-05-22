# mTopic figure reproducibility

This repository contains the scripts used to reproduce the figures from the manuscript:

> **Identifying multimodal molecular programs with mTopic**

The repository is designed for figure reproduction using processed datasets and intermediate analysis files hosted on Zenodo.

## Repository layout

```text
.
├── data/       # processed datasets and intermediate files downloaded from Zenodo
├── figures/    # output directory created/populated by the figure scripts
└── scripts/    # Python and R scripts used to generate manuscript figures
```

## Data availability

The processed datasets and intermediate analysis files required for figure reproduction are available from Zenodo:

**Zenodo record:** <https://zenodo.org/records/20341450>

Download all files from the Zenodo record and place them directly in the `data/` directory at the root of this repository.

After downloading, the expected layout is:

```text
.
├── data/
│   ├── P22MouseBrainATAC_trained.h5mu
│   ├── HumanTonsil_trained.h5mu
│   ├── HumanPBMC_trained.h5mu
│   └── ...
├── figures/
└── scripts/
```

## Reproducing the figures

Run the figure scripts from inside the `scripts/` directory.

Each script creates a figure-specific subdirectory inside `figures/` and saves the corresponding plots there. For example, after running the Figure 1 scripts, the output may look like this:

```text
figures/
└── Figure_1/
    ├── Figure_1b.pdf
    ├── Figure_1c.pdf
    ├── Figure_1d.pdf
    └── Figure_1e.pdf
```

Extended Data and Supplementary figure scripts write to corresponding folders, for example:

```text
figures/
├── Figure_1/
├── Figure_2/
├── Figure_E1/
├── Figure_S1/
└── ...
```

## Software environment

### Python environment

The Python figure-generation scripts were run with **Python 3.11** and the following package versions:

```text
anndata: 0.11.4
h5py: 3.16.0
kneed: 0.8.5
matplotlib: 3.10.9
mtopic: 1.1
mudata: 0.3.8
muon: 0.1.7
natsort: 8.4.0
networkx: 3.6.1
numba: 0.65.1
numpy: 1.26.1
pandas: 2.2.2
pillow: 12.2.0
pynndescent: 0.6.0
scanpy: 1.11.5
scikit-learn: 1.3.0
scipy: 1.11.3
seaborn: 0.13.2
statsmodels: 0.14.6
torch: 2.6.0+cu124
tqdm: 4.67.3
umap-learn: 0.5.6
```

Additional Python standard-library modules used by the scripts include:

```text
os, re, json, math, pickle, collections
```

### R environment

The R figure-generation scripts were run with **R 4.5.2** and the following package versions:

```text
Seurat: 5.4.0
SeuratObject: 5.3.0
SeuratData: 0.2.2.9002
Signac: 1.16.0
data.table: 1.18.2.1
dplyr: 1.2.1
tidyr: 1.3.2
tibble: 3.3.1
purrr: 1.2.1
stringr: 1.6.0
ggplot2: 4.0.2
ggrepel: 0.9.8
ggforce: 0.5.0
ggplotify: 0.1.3
ggraph: 2.2.2
tidygraph: 1.3.1
gridExtra: 2.3
gtable: 0.3.6
patchwork: 1.3.2
cowplot: 1.2.0
pheatmap: 1.0.13
ComplexHeatmap: 2.26.1
circlize: 0.4.17
RColorBrewer: 1.1.3
viridis: 0.6.5
colorRamps: 2.3.4
scales: 1.4.0
igraph: 2.2.3
rgexf: 0.16.3
networkD3: 0.4.1
htmlwidgets: 1.6.4
webshot2: 0.1.2
reshape2: 1.4.5
jsonlite: 2.0.0
this.path: 2.8.0
MASS: 7.3.65
GenomicRanges: 1.62.1
Rsamtools: 2.26.0
clusterProfiler: 4.18.1
org.Hs.eg.db: 3.22.0
seqLogo: 1.76.0
MuDataSeurat: 0.0.0.9000
TFBSTools: 1.48.0
ggalign: 1.2.0
ggpubr: 0.6.2
msigdbr: 26.1.0
```

Base R packages used by the scripts include:

```text
stats, utils, graphics, grDevices, grid
```
