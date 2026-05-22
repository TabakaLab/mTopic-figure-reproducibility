library(Seurat)
library(RColorBrewer)
library(scales)
library(ggplot2)
library(TFBSTools)
library(seqLogo)
library(ggplotify)
library(ggpubr)
library(grid)
library(this.path)

setwd(this.dir())

dir.create("../figures/Figure_E2", recursive = T)
data <- readRDS("../data/Seurat_RNA_ATAC.RDS")
Genes <- c("Sox10","Mef2c","Sox2","Pbx3","Nfib","Sox4","Sox8","Pou2f1")
motifs <- readRDS("../data/JASPAR2022.RDS")
plot <- FeaturePlot(data, features = Genes, pt.size = 0.01, max.cutoff = 3,
                    keep.scale = "all", ncol = 4) & scale_colour_gradientn(
                    colours = c("grey90", rev(brewer.pal(n = 10, name = "RdBu"))[5:10],
                                "black", "black"), limits = c(0, 3), oob = squish) &
  theme(legend.position = "bottom")

ggsave("../figures/Figure_E2/Figure_E2e_1.pdf", width = 8, height = 5.5, plot)


motif_vec <- c("MA0442.2", "MA0497.1", "MA0143.4", "MA1114.1",
               "MA1643.1", "MA0867.2", "MA0868.2", "MA0785.1")
names(motif_vec) <- c("Sox10", "Mef2c", "Sox2", "Pbx3", 
                      "Nfib", "Sox4", "Sox8", "Pou2f1")

plotlist <- list()
for (gene_name in names(motif_vec)) {
  motif_id <- motif_vec[gene_name]
  pfm <- motifs[[motif_id]]
  icm <- toICM(pfm, pseudocounts = 0.1, schneider = TRUE)
  plotlist[[gene_name]] <- TFBSTools::seqLogo(icm, ic.scale = TRUE)
  p_logo <- as.ggplot(expression(TFBSTools::seqLogo(icm, ic.scale = TRUE))) +
    ggtitle(paste0(gene_name, " (", motif_id, ")")) +
    theme(plot.title = element_text(hjust = 0.5, size = 10))
  plotlist[[gene_name]] <- p_logo
}
plot2 <- ggpubr::ggarrange(plotlist = plotlist, ncol = 4, nrow = 2)

ggsave("../figures/Figure_E2/Figure_E2e_2.pdf", width = 20, height = 5, plot2)

# > sessionInfo()
# R version 4.5.2 (2025-10-31)
# Platform: x86_64-pc-linux-gnu
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
#   [1] stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
#   [1] ggplot2_4.0.2      scales_1.4.0       RColorBrewer_1.1-3 Seurat_5.4.0       SeuratObject_5.3.0 sp_2.2-1          
# 
# loaded via a namespace (and not attached):
#   [1] deldir_2.0-4           pbapply_1.7-4          gridExtra_2.3          rlang_1.2.0            magrittr_2.0.5         RcppAnnoy_0.0.23      
# [7] otel_0.2.0             matrixStats_1.5.0      ggridges_0.5.7         compiler_4.5.2         spatstat.geom_3.7-3    systemfonts_1.3.2     
# [13] png_0.1-9              vctrs_0.7.2            reshape2_1.4.5         stringr_1.6.0          pkgconfig_2.0.3        fastmap_1.2.0         
# [19] labeling_0.4.3         promises_1.5.0         ragg_1.5.0             purrr_1.2.1            jsonlite_2.0.0         goftest_1.2-3         
# [25] later_1.4.8            spatstat.utils_3.2-2   irlba_2.3.7            parallel_4.5.2         cluster_2.1.8.1        R6_2.6.1              
# [31] ica_1.0-3              stringi_1.8.7          spatstat.data_3.1-9    reticulate_1.46.0      parallelly_1.46.1      spatstat.univar_3.1-7 
# [37] lmtest_0.9-40          scattermore_1.2        Rcpp_1.1.1             tensor_1.5.1           future.apply_1.20.2    zoo_1.8-15            
# [43] sctransform_0.4.3      httpuv_1.6.17          Matrix_1.7-4           splines_4.5.2          igraph_2.2.3           tidyselect_1.2.1      
# [49] rstudioapi_0.18.0      abind_1.4-8            spatstat.random_3.4-5  codetools_0.2-20       miniUI_0.1.2           spatstat.explore_3.8-0
# [55] listenv_0.10.1         lattice_0.22-7         tibble_3.3.1           plyr_1.8.9             withr_3.0.2            shiny_1.13.0          
# [61] S7_0.2.1               ROCR_1.0-12            Rtsne_0.17             future_1.70.0          fastDummies_1.7.5      survival_3.8-3        
# [67] polyclip_1.10-7        fitdistrplus_1.2-6     pillar_1.11.1          KernSmooth_2.23-26     plotly_4.12.0          generics_0.1.4        
# [73] RcppHNSW_0.6.0         globals_0.19.1         xtable_1.8-8           glue_1.8.0             lazyeval_0.2.3         tools_4.5.2           
# [79] data.table_1.18.2.1    RSpectra_0.16-2        RANN_2.6.2             dotCall64_1.2          cowplot_1.2.0          grid_4.5.2            
# [85] tidyr_1.3.2            nlme_3.1-168           patchwork_1.3.2        cli_3.6.6              spatstat.sparse_3.1-0  textshaping_1.0.5     
# [91] spam_2.11-3            viridisLite_0.4.3      dplyr_1.2.1            uwot_0.2.4             gtable_0.3.6           digest_0.6.39         
# [97] progressr_0.19.0       ggrepel_0.9.8          htmlwidgets_1.6.4      farver_2.1.2           htmltools_0.5.9        lifecycle_1.0.5       
# [103] httr_1.4.8             mime_0.13              MASS_7.3-65  