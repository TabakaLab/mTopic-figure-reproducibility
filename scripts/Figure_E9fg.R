require(data.table)
library(ggplot2)
library(dplyr)
library(tidyr)
library(grid)
library(SeuratData)
library(MuDataSeurat)
library(viridis)
library(gridExtra)
library(Seurat)
library(Signac)
library(patchwork)
library(RColorBrewer)
library(this.path)

setwd(this.dir())

dir.create("../figures/Figure_E9", recursive = T)

RTE <- readRDS("../data/RTE.RDS")
data <- readRDS("../data/Seurat_PBMCs.RDS")
dictionary <- data.frame(fread("../data/cell_type_annotations_PBMCs.csv",header=T))

DefaultAssay(data) <- "rna"
Genes <- c("rna:TBX21", "rna:EOMES", "rna:RORA")
p <- FeaturePlot(
  data,
  features = Genes,
  order = TRUE,
  ncol = 1,
  pt.size = 0.1,
  keep.scale = "all"
) &
  scale_colour_gradientn(
    colours = c(
      rep("grey90",5),
      rev(brewer.pal(n = 10, name = "RdBu"))[7:10],
      "black"
    )
  ) &
  theme(legend.position = "right")

ggsave("../figures/Figure_E9/Figure_E9f.pdf", width = 5, height = 12, p)

celltypes <- as.vector(unique(data@meta.data$cell_type))
G1 <- names(unlist(RTE$TBX21))
G2 <- names(unlist(RTE$EOMES))
G3 <- names(unlist(RTE$RORA))
G1 <- list(paste0("rna:", G1))
G2 <- list(paste0("rna:", G2))
G3 <- list(paste0("rna:", G3))

data <- AddModuleScore(nbin = 20, object = data, features = G1, ctrl = 100, name = "G1")
data <- AddModuleScore(nbin = 20, object = data, features = G2, ctrl = 100, name = "G2")
data <- AddModuleScore(nbin = 20, object = data, features = G3, ctrl = 100, name = "G3")

df <- data[[]][, c("cell_type", "G11", "G21", "G31"), drop = FALSE]
df <- df %>%
  dplyr::rename(
    TBX21 = G11,
    EOMES = G21,
    RORA  = G31
  ) %>%
  mutate(cell_type = as.character(cell_type)) %>%
  filter(
    !is.na(cell_type),
    cell_type != "",
    toupper(cell_type) != "NA",
    !grepl("doublet", cell_type, ignore.case = TRUE)
  )

cell_order <- if (is.factor(data$cell_type)) {
  levels(data$cell_type)
} else {
  unique(as.character(data$cell_type))
}

cell_order <- cell_order[
  !is.na(cell_order) &
    cell_order != "" &
    toupper(cell_order) != "NA" &
    !grepl("doublet", cell_order, ignore.case = TRUE)
]

df$cell_type <- factor(df$cell_type, levels = cell_order)

## colors
if (!is.null(names(dictionary$colors)) && all(cell_order %in% names(dictionary$colors))) {
  cols_use <- dictionary$colors[cell_order]
} else {
  cols_use <- setNames(dictionary$colors[seq_along(cell_order)], cell_order)
}

## long format
df_long <- df %>%
  pivot_longer(
    cols = c(TBX21, EOMES, RORA),
    names_to = "signature",
    values_to = "score"
  ) %>%
  mutate(signature = factor(signature, levels = c("TBX21", "EOMES", "RORA")))

## labels inside each panel
lab_df <- data.frame(
  signature = factor(c("TBX21", "EOMES", "RORA"), levels = c("TBX21", "EOMES", "RORA")),
  cell_type = factor(rep(cell_order[1], 3), levels = cell_order),
  y = c(1.15, 1.15, 0.95),
  label = c("TBX21", "EOMES", "RORA")
)

p <- ggplot(df_long, aes(x = cell_type, y = score, fill = cell_type)) +
  geom_violin(
    width = 0.9,
    scale = "width",
    trim = TRUE,
    color = "#3a3a3a",
    linewidth = 0.8
  ) +
  facet_grid(
    rows = vars(signature),
    scales = "free_y",
    switch = "y"
  ) +
  geom_text(
    data = lab_df,
    aes(x = cell_type, y = y, label = label),
    inherit.aes = FALSE,
    hjust = 0,
    vjust = 1,
    size = 10
  ) +
  scale_fill_manual(values = cols_use) +
  labs(x = NULL, y = "RNA signature scores") +
  coord_cartesian(clip = "off") +
  theme_classic(base_size = 20) +
  theme(
    legend.position = "none",
    strip.background = element_blank(),
    strip.text = element_blank(),
    panel.spacing.y = unit(0, "mm"),   
    axis.line = element_line(color = "black", linewidth = 1.2),
    axis.ticks = element_line(color = "black", linewidth = 1),
    axis.text.x = element_text(
      angle = 90,
      vjust = 0.5,
      hjust = 1,
      color = "black",
      size = 16
    ),
    axis.text.y = element_text(color = "black", size = 16),
    axis.title.y = element_text(size = 22),
    panel.background = element_blank(),
    plot.background = element_blank(),
    plot.margin = margin(5, 20, 5, 40)
  )

ggsave("../figures/Figure_E9/Figure_E9g.pdf", width = 12, height = 8, p)

# > sessionInfo()
# R version 4.5.2 (2025-10-31)
# Running under: Ubuntu 24.04.4 LTS
# 
# Matrix products: default
# BLAS:   /usr/lib/x86_64-linux-gnu/openblas-pthread/libblas.so.3 
# LAPACK: /usr/lib/x86_64-linux-gnu/openblas-pthread/libopenblasp-r0.3.26.so;  LAPACK version 3.12.0
# 
# locale:
#   [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C               LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8     LC_MONETARY=en_US.UTF-8   
# [6] LC_MESSAGES=en_US.UTF-8    LC_PAPER=en_US.UTF-8       LC_NAME=C                  LC_ADDRESS=C               LC_TELEPHONE=C            
# [11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       
# 
# time zone: Europe/Warsaw
# tzcode source: system (glibc)
# 
# attached base packages:
#   [1] grid      stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
#   [1] RColorBrewer_1.1-3      patchwork_1.3.2         Signac_1.16.0           Seurat_5.4.0            SeuratObject_5.3.0     
# [6] sp_2.2-1                gridExtra_2.3           viridis_0.6.5           viridisLite_0.4.3       MuDataSeurat_0.0.0.9000
# [11] SeuratData_0.2.2.9002   tidyr_1.3.2             dplyr_1.2.1             data.table_1.18.2.1     ggpubr_0.6.2           
# [16] ggplot2_4.0.2          
# 
# loaded via a namespace (and not attached):
#   [1] rstudioapi_0.18.0      jsonlite_2.0.0         magrittr_2.0.5         spatstat.utils_3.2-2   farver_2.1.2           ragg_1.5.0            
# [7] vctrs_0.7.2            ROCR_1.0-12            Rsamtools_2.26.0       spatstat.explore_3.8-0 RcppRoll_0.3.1         rstatix_0.7.3         
# [13] htmltools_0.5.9        broom_1.0.12           Formula_1.2-5          sctransform_0.4.3      parallelly_1.46.1      KernSmooth_2.23-26    
# [19] htmlwidgets_1.6.4      ica_1.0-3              plyr_1.8.9             plotly_4.12.0          zoo_1.8-15             igraph_2.2.3          
# [25] mime_0.13              lifecycle_1.0.5        pkgconfig_2.0.3        Matrix_1.7-4           R6_2.6.1               fastmap_1.2.0         
# [31] fitdistrplus_1.2-6     future_1.70.0          shiny_1.13.0           digest_0.6.39          S4Vectors_0.48.0       tensor_1.5.1          
# [37] RSpectra_0.16-2        irlba_2.3.7            GenomicRanges_1.62.1   textshaping_1.0.5      labeling_0.4.3         progressr_0.19.0      
# [43] spatstat.sparse_3.1-0  httr_1.4.8             polyclip_1.10-7        abind_1.4-8            compiler_4.5.2         bit64_4.6.0-1         
# [49] withr_3.0.2            S7_0.2.1               backports_1.5.1        BiocParallel_1.44.0    carData_3.0-6          fastDummies_1.7.5     
# [55] ggsignif_0.6.4         MASS_7.3-65            rappdirs_0.3.4         tools_4.5.2            lmtest_0.9-40          otel_0.2.0            
# [61] httpuv_1.6.17          future.apply_1.20.2    goftest_1.2-3          glue_1.8.0             nlme_3.1-168           promises_1.5.0        
# [67] Rtsne_0.17             cluster_2.1.8.1        reshape2_1.4.5         generics_0.1.4         hdf5r_1.3.12           gtable_0.3.6          
# [73] spatstat.data_3.1-9    XVector_0.50.0         car_3.1-5              BiocGenerics_0.56.0    spatstat.geom_3.7-3    RcppAnnoy_0.0.23      
# [79] ggrepel_0.9.8          RANN_2.6.2             pillar_1.11.1          stringr_1.6.0          spam_2.11-3            RcppHNSW_0.6.0        
# [85] later_1.4.8            splines_4.5.2          lattice_0.22-7         survival_3.8-3         bit_4.6.0              deldir_2.0-4          
# [91] tidyselect_1.2.1       Biostrings_2.78.0      miniUI_0.1.2           pbapply_1.7-4          Seqinfo_1.0.0          IRanges_2.44.0        
# [97] scattermore_1.2        stats4_4.5.2           matrixStats_1.5.0      UCSC.utils_1.6.1       stringi_1.8.7          lazyeval_0.2.3        
# [103] codetools_0.2-20       tibble_3.3.1           cli_3.6.6              uwot_0.2.4             xtable_1.8-8           reticulate_1.46.0     
# [109] systemfonts_1.3.2      GenomeInfoDb_1.46.2    Rcpp_1.1.1             globals_0.19.1         spatstat.random_3.4-5  png_0.1-9             
# [115] spatstat.univar_3.1-7  parallel_4.5.2         dotCall64_1.2          bitops_1.0-9           listenv_0.10.1         scales_1.4.0          
# [121] ggridges_0.5.7         purrr_1.2.1            crayon_1.5.3           rlang_1.2.0            fastmatch_1.1-8        cowplot_1.2.0    