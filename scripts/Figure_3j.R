library(SeuratData)
library(MuDataSeurat)
library(viridis)
library(gridExtra)
library(Seurat)
library(Signac)
library(stringr)
library(GenomicRanges)
library(data.table)
library(stringr)
library(ggplot2)
library(scales)
library(Rsamtools)
library(this.path)

setwd(this.dir())

dir.create("../figures/Figure_3", recursive = T)

data <- readRDS("../data/Seurat_PBMCs.RDS")
data[["peaks"]]@fragments[[1]]@path <- "../data/GSM5065524_LLL_CTRL_fragments.tsv.gz"
d <- subset(data, cell_type %in% setdiff(names(table(data@meta.data$cell_type)),
  c("NA doublets-1", "NA doublets-2", "NA doublets-3")))

Topic <- 15
TOI <- c(15, 20, 7, 28)
Gene <- "PRF1"
Reg <- "chr10-70551963-70652264"
DFLinks <- readRDS(file="../data/PRF_GlobalModel.RDS")
min.color <- 0

GlobalModel <- ggplot(data = DFLinks) +
  ggforce::geom_bezier(
    mapping = aes_string(x = "x", y = "y", group='group', color = "score"))

GlobalModel <- GlobalModel +
  geom_hline(yintercept = 0, color = 'grey') +
  scale_color_continuous(low = "grey90",  high = "black",
                         limits = c(min(DFLinks$score), max(DFLinks$score)),
                         n.breaks = 9)

gr <- unlist(strsplit(Reg, split="-"))

GlobalModel <- GlobalModel +
  theme_classic() +
  theme(axis.ticks.y = element_blank(),
        axis.text.y = element_blank()) +
  ylab("Links") +
  xlab(label = paste0(gr[1], " position (bp)")) +
  xlim(c(as.integer(gr[2]), as.integer(gr[3])))

atac_signatures <- "../data/PBMC_dogmaseq_signatures_atac.csv"
atac_signatures <- fread(atac_signatures)
tmp <- colnames(atac_signatures)
tmp[1] <- "V1"
colnames(atac_signatures) <- tmp

Elements <- readRDS(file="../data/PRF1_elementsIDs.RDS")
PromoterID <- Elements[[1]]
DistalElementsID <- Elements[[2]]

CreateLink<-function(peak_1, peak_2,gr, int_score){
  df <- data.frame(
    x = c((peak_1$start + peak_1$end)/2,
          ((peak_1$start + peak_1$end)/2 + (peak_2$start + peak_2$end)/2)/2,
          (peak_2$start + peak_2$end)/2),
    y = c(0, -int_score, 0),
    group = rep(x = gr, 3),
    score = rep(int_score, 3)
  )
  df <- df[order(df$x, decreasing = F),]
  return(df)
}

GP_assoc <- list()
for(i in 1:length(TOI)){
  ListLinks <- lapply(1:nrow(DistalElementsID), function(x) CreateLink(
    PromoterID, DistalElementsID[x,], x, DistalElementsID[x, paste0("topic_", TOI[i])]))
  DFLinks <- do.call(rbind, ListLinks)
  min.color <- 0
  
  p1 <- ggplot(data = DFLinks) + ggforce::geom_bezier(mapping = aes_string(
    x = "x", y = "y",group='group', color = "score"))
  
  p1 <- p1 + geom_hline(yintercept = 0, color = 'grey') +
    scale_color_continuous(low = "grey90", high = "black",
                           limits = c(min(DFLinks$score), max(DFLinks$score)),
                           n.breaks = 9)
  
  gr <- unlist(strsplit(Reg, split="-"))
  p1 <- p1 +
    theme_classic() +
    theme(axis.ticks.y = element_blank(),
          axis.text.y = element_blank()) +
    ylab("Links") +
    xlab(label = paste0(gr[1], " position (bp)")) +
    xlim(c(as.integer(gr[2]), as.integer(gr[3])))
  GP_assoc[[i]] <- p1
}

cell_type <- read.csv("../data/cell_type_annotations_PBMCs.csv", row.names = 1)
names(GP_assoc) <- cell_type$cell_type[match(paste0("topic_", TOI), cell_type$topic)]

ranges.show <- StringToGRanges("chr10-70572022-70572023")
Ext <- 1000
end(ranges.show) <- end(ranges.show) + Ext
start(ranges.show) <- start(ranges.show) - Ext
ranges.show$color <- "red"

cov_plot <- CoveragePlot(object = d,
                         assay = "peaks",
                         group.by = "cell_type",
                         window = 200,
                         scale.factor = 1e7,
                         downsample.rate = 1,
                         ymax=500,
                         assay.scale = "separate",
                         peaks.group.by = cell_type,	
                         show.bulk = F,
                         region = Reg,
                         annotation = F,
                         peaks = F,
                         region.highlight = ranges.show
) & scale_fill_manual(values = cell_type$colors)

cov_plot <- cov_plot + theme(panel.border = element_rect(
  color = "black", fill = NA, size = 0.4), text = element_text(size = 14)) +  
  scale_x_continuous(labels = label_number(scale = 1e-6, accuracy = 0.001, suffix = "")) +
  labs(x = gsub(pattern = "bp", replacement = "Mbp", cov_plot$labels$x))

bulk <- CoveragePlot(assay = "peaks",
                     group.by = "orig.ident",
                     window=200,#scale.factor=1e5,    
                     object = data,
                     region = Reg,
                     ymax=300,
                     annotation = TRUE,
                     peaks = FALSE
) & scale_fill_manual(values = "grey60")

CT <- CombineTracks(
  plotlist = list(cov_plot, bulk, GlobalModel, GP_assoc$`NK cells`),
  heights = c(15,2,2,2),
  widths = c(10, 1)
)

pdf(file = paste0("../figures/Figure_3/Figure_3j_1.pdf"), width = 7, height=7)
  CT
dev.off()

expr_plot <- ExpressionPlot(
  group.by = "cell_type", object = d,
  features = paste0("rna:", Gene),
  assay = "rna"
  ) & scale_fill_manual(values = cell_type$colors)

pdf(file = paste0("../figures/Figure_3/Figure_3j_2.pdf"), width = 1.5, height=6)
  expr_plot
dev.off()


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
#   [1] stats4    stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
#   [1] Rsamtools_2.26.0        Biostrings_2.78.0       XVector_0.50.0          scales_1.4.0            ggplot2_4.0.2          
# [6] data.table_1.18.2.1     GenomicRanges_1.62.1    Seqinfo_1.0.0           IRanges_2.44.0          S4Vectors_0.48.0       
# [11] BiocGenerics_0.56.0     generics_0.1.4          stringr_1.6.0           Signac_1.16.0           Seurat_5.4.0           
# [16] SeuratObject_5.3.0      sp_2.2-1                gridExtra_2.3           viridis_0.6.5           viridisLite_0.4.3      
# [21] MuDataSeurat_0.0.0.9000 SeuratData_0.2.2.9002  
# 
# loaded via a namespace (and not attached):
#   [1] RColorBrewer_1.1-3     rstudioapi_0.18.0      jsonlite_2.0.0         magrittr_2.0.5         spatstat.utils_3.2-2   farver_2.1.2          
# [7] vctrs_0.7.2            ROCR_1.0-12            spatstat.explore_3.8-0 RcppRoll_0.3.1         htmltools_0.5.9        sctransform_0.4.3     
# [13] parallelly_1.46.1      KernSmooth_2.23-26     htmlwidgets_1.6.4      ica_1.0-3              plyr_1.8.9             plotly_4.12.0         
# [19] zoo_1.8-15             igraph_2.2.3           mime_0.13              lifecycle_1.0.5        pkgconfig_2.0.3        Matrix_1.7-4          
# [25] R6_2.6.1               fastmap_1.2.0          fitdistrplus_1.2-6     future_1.70.0          shiny_1.13.0           digest_0.6.39         
# [31] patchwork_1.3.2        tensor_1.5.1           RSpectra_0.16-2        irlba_2.3.7            labeling_0.4.3         progressr_0.19.0      
# [37] spatstat.sparse_3.1-0  httr_1.4.8             polyclip_1.10-7        abind_1.4-8            compiler_4.5.2         withr_3.0.2           
# [43] bit64_4.6.0-1          S7_0.2.1               BiocParallel_1.44.0    fastDummies_1.7.5      ggforce_0.5.0          MASS_7.3-65           
# [49] rappdirs_0.3.4         tools_4.5.2            lmtest_0.9-40          otel_0.2.0             httpuv_1.6.17          future.apply_1.20.2   
# [55] goftest_1.2-3          glue_1.8.0             nlme_3.1-168           promises_1.5.0         grid_4.5.2             Rtsne_0.17            
# [61] cluster_2.1.8.1        reshape2_1.4.5         hdf5r_1.3.12           gtable_0.3.6           spatstat.data_3.1-9    tidyr_1.3.2           
# [67] spatstat.geom_3.7-3    RcppAnnoy_0.0.23       ggrepel_0.9.8          RANN_2.6.2             pillar_1.11.1          spam_2.11-3           
# [73] RcppHNSW_0.6.0         later_1.4.8            splines_4.5.2          tweenr_2.0.3           dplyr_1.2.1            lattice_0.22-7        
# [79] survival_3.8-3         bit_4.6.0              deldir_2.0-4           tidyselect_1.2.1       miniUI_0.1.2           pbapply_1.7-4         
# [85] scattermore_1.2        matrixStats_1.5.0      stringi_1.8.7          UCSC.utils_1.6.1       lazyeval_0.2.3         codetools_0.2-20      
# [91] tibble_3.3.1           cli_3.6.6              uwot_0.2.4             xtable_1.8-8           reticulate_1.46.0      Rcpp_1.1.1            
# [97] GenomeInfoDb_1.46.2    globals_0.19.1         spatstat.random_3.4-5  png_0.1-9              spatstat.univar_3.1-7  parallel_4.5.2        
# [103] dotCall64_1.2          bitops_1.0-9           listenv_0.10.1         ggridges_0.5.7         purrr_1.2.1            crayon_1.5.3          
# [109] rlang_1.2.0            fastmatch_1.1-8        cowplot_1.2.0    