library(SeuratData)
library(MuDataSeurat)
library(viridis)
library(gridExtra)
library(ggplot2)
library(RColorBrewer)
library(cowplot)
library(Seurat)
library(Signac)
library(this.path)

setwd(this.dir())

dir.create("../figures/Figure_E7", recursive = T)


data<-readRDS("../data/HumanPBMC_RNA_unfiltered.RDS")
annot<-read.csv("../data/HumanPBMC_annotations.csv")

data<-NormalizeData(data, normalization.method = "CLR", assay = "prot")
data<-NormalizeData(data, assay = "rna")

# ═══════════════════════════════════════════════════════════════════════════
# Figure E6f
# ═══════════════════════════════════════════════════════════════════════════

DefaultAssay(data)<-"prot"

prot<-c("prot:CD45RA", "prot:CD45RO")
plotList<-lapply(prot, function(x) {
  FeaturePlot(data, features=x, pt.size = 0.2, order = T)+
    scale_color_gradientn(colours = c('#E4E4E4', '#E4E4E4','#E4E4E4', '#E4E4E4','#CDAC3B', '#967402', '#000000'))+
    theme(axis.text=element_blank(),
          axis.ticks=element_blank(),
          panel.border = element_rect(colour = "black", fill = NA, size = 1))+
    coord_fixed()
})

ml <- marrangeGrob(plotList, layout_matrix = matrix(1:2, 1, 2, TRUE))
ggsave("../figures/Figure_E7/Figure_E7f.png", ml, width = 10, height = 5, bg = "white")

# ═══════════════════════════════════════════════════════════════════════════
# Figure E6h
# ═══════════════════════════════════════════════════════════════════════════

cell_types<-levels(data$cell_type)[grep("CD4", levels(data$cell_type))]
cell_types<-cell_types[-grep("Naive", cell_types)]

data_sub<-subset(data, subset = cell_type %in% cell_types)

annot2<-annot[annot$cell_type %in% cell_types,]
annot2$color[1:4]<-"grey80"
plot<-DimPlot(data, group.by = "cell_type", cols = annot2$color, pt.size = 2.5)+
  xlim(-7,-1)+
  ylim(-2.5, 4.5)+
  theme_nothing()

DefaultAssay(data_sub)<-"rna"
marker<-c("rna:BCL6", "rna:TBX21", "rna:GATA3", "rna:RORC", "rna:FOXP3")

plotList<-lapply(marker, function(x) {
  FeaturePlot(data_sub, features=x, order = T, pt.size = 0.5, max.cutoff = 'q99')+
    xlim(-7,-1)+
    ylim(-2, 4)+
    scale_color_gradientn(colours = c('#E4E4E4','#E4E4E4', '#DE2B25', '#99000D', '#000000'))+
    theme(axis.title=element_blank(),
          axis.text=element_blank(),
          axis.ticks=element_blank())+
    coord_fixed()
})

ml <- marrangeGrob(plotList, layout_matrix = matrix(1:5, 1, 5, TRUE))
ggsave("../figures/Figure_E7/Figure_E7h_1.png", ml, width = 15, height = 3, bg = "white")


DefaultAssay(data)<-"prot"
DefaultAssay(data_sub)<-"prot"


prot<-c("prot:CD278", "prot:CD195", "prot:CD69", "prot:CD71", "prot:CD25")

plotList<-lapply(prot, function(x) {
  FeaturePlot(data_sub, features=x, pt.size = 0.2, order = T)+
    xlim(-7,-1)+
    ylim(-2, 4)+
    scale_color_gradientn(colours = c('#E4E4E4', '#E4E4E4','#E4E4E4', '#E4E4E4','#CDAC3B', '#967402', '#000000'))+
    theme(axis.text=element_blank(),
          axis.ticks=element_blank(),
          panel.border = element_rect(colour = "black", fill = NA, size = 1))+
    coord_fixed()
})

ml <- marrangeGrob(plotList, layout_matrix = matrix(1:5, 1, 5, TRUE))
ggsave("../figures/Figure_E7/Figure_E7h_2.png", ml, width = 15, height = 3, bg = "white")

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
#   [1] Signac_1.16.0           Seurat_5.4.0            SeuratObject_5.3.0      sp_2.2-1                cowplot_1.2.0          
# [6] RColorBrewer_1.1-3      ggplot2_4.0.2           gridExtra_2.3           viridis_0.6.5           viridisLite_0.4.3      
# [11] MuDataSeurat_0.0.0.9000 SeuratData_0.2.2.9002  
# 
# loaded via a namespace (and not attached):
#   [1] rstudioapi_0.18.0      jsonlite_2.0.0         magrittr_2.0.5         spatstat.utils_3.2-2   farver_2.1.2           ragg_1.5.0            
# [7] vctrs_0.7.2            ROCR_1.0-12            spatstat.explore_3.8-0 Rsamtools_2.26.0       RcppRoll_0.3.1         htmltools_0.5.9       
# [13] sctransform_0.4.3      parallelly_1.46.1      KernSmooth_2.23-26     htmlwidgets_1.6.4      ica_1.0-3              plyr_1.8.9            
# [19] plotly_4.12.0          zoo_1.8-15             igraph_2.2.3           mime_0.13              lifecycle_1.0.5        pkgconfig_2.0.3       
# [25] Matrix_1.7-4           R6_2.6.1               fastmap_1.2.0          fitdistrplus_1.2-6     future_1.70.0          shiny_1.13.0          
# [31] digest_0.6.39          patchwork_1.3.2        S4Vectors_0.48.0       tensor_1.5.1           RSpectra_0.16-2        irlba_2.3.7           
# [37] textshaping_1.0.5      GenomicRanges_1.62.1   labeling_0.4.3         progressr_0.19.0       spatstat.sparse_3.1-0  httr_1.4.8            
# [43] polyclip_1.10-7        abind_1.4-8            compiler_4.5.2         bit64_4.6.0-1          withr_3.0.2            S7_0.2.1              
# [49] BiocParallel_1.44.0    fastDummies_1.7.5      MASS_7.3-65            rappdirs_0.3.4         tools_4.5.2            lmtest_0.9-40         
# [55] otel_0.2.0             httpuv_1.6.17          future.apply_1.20.2    goftest_1.2-3          glue_1.8.0             nlme_3.1-168          
# [61] promises_1.5.0         grid_4.5.2             Rtsne_0.17             cluster_2.1.8.1        reshape2_1.4.5         generics_0.1.4        
# [67] hdf5r_1.3.12           gtable_0.3.6           spatstat.data_3.1-9    tidyr_1.3.2            data.table_1.18.2.1    XVector_0.50.0        
# [73] BiocGenerics_0.56.0    spatstat.geom_3.7-3    RcppAnnoy_0.0.23       ggrepel_0.9.8          RANN_2.6.2             pillar_1.11.1         
# [79] stringr_1.6.0          spam_2.11-3            RcppHNSW_0.6.0         later_1.4.8            splines_4.5.2          dplyr_1.2.1           
# [85] lattice_0.22-7         survival_3.8-3         bit_4.6.0              deldir_2.0-4           tidyselect_1.2.1       Biostrings_2.78.0     
# [91] miniUI_0.1.2           pbapply_1.7-4          IRanges_2.44.0         Seqinfo_1.0.0          scattermore_1.2        stats4_4.5.2          
# [97] matrixStats_1.5.0      stringi_1.8.7          UCSC.utils_1.6.1       lazyeval_0.2.3         codetools_0.2-20       tibble_3.3.1          
# [103] cli_3.6.6              uwot_0.2.4             systemfonts_1.3.2      xtable_1.8-8           reticulate_1.46.0      Rcpp_1.1.1            
# [109] GenomeInfoDb_1.46.2    globals_0.19.1         spatstat.random_3.4-5  png_0.1-9              spatstat.univar_3.1-7  parallel_4.5.2        
# [115] dotCall64_1.2          bitops_1.0-9           listenv_0.10.1         scales_1.4.0           ggridges_0.5.7         purrr_1.2.1           
# [121] crayon_1.5.3           rlang_1.2.0            fastmatch_1.1-8      