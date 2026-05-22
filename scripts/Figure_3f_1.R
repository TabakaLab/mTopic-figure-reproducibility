# Second part of the figure is generated in python

library(grid)
library(viridis)
library(pheatmap)
library(RColorBrewer)
library(ggplot2)
library(this.path)

setwd(this.dir())

dir.create("../figures/Figure_3", recursive = T)

cell_type <- read.csv("../data/cell_type_annotations_PBMCs.csv", row.names = 1)
data <- readRDS("../data/TFMotifs_PBMCs.RDS")

row_anno <- data.frame(RowCat = rownames(data), row.names = rownames(data))
col_anno <- data.frame(ColCat = cell_type$cell_type, row.names = cell_type$cell_type)

unique_rows <- unique(row_anno$RowCat)
row_palette <- turbo(length(unique_rows))
names(row_palette) <- unique_rows

unique_cols <- unique(col_anno$ColCat)
col_palette <- cell_type$colors
names(col_palette) <- unique_cols

ann_colors <- list(RowCat = row_palette, ColCat = col_palette)

plot <- pheatmap(
  data,cluster_rows = F, cluster_cols = F, border_color = NA,
  annotation_row = row_anno, annotation_col = col_anno,
  annotation_colors = ann_colors, fontsize_row = 6, 
  annotation_names_row = F, fontsize = 6,
  color = colorRampPalette(colors = c(
    rev(brewer.pal(n = 9, name = 'Blues')),
    "white", "White", "White", "White", "White", "White", "White",
    brewer.pal(n = 8, name = 'Greens')))(100)
  )


labels <- c(
  "STAT1::STAT2--MA0517.1", "IRF8--MA0652.1", "IRF4--MA1419.1", 
  "ETV3--MA0763.1", "TCF7--MA0769.2", "ELK3--MA0759.2", "TCF7L2--MA0523.1",
  "RUNX2--MA0511.2", "CTCF--MA0139.1", "STAT3--MA0144.2", "POU2F1--MA0785.1",
  "TBX21--MA0690.2", "EOMES--MA0800.1", "RORA--MA0072.1", "NFKB1--MA0105.4",
  "NFKB2--MA0778.1", "RELB--MA1117.1", "FOS--MA1951.1", "JUN--MA0488.1",
  "BATF--MA1634.1", "BATF3--MA0835.2", "CEBPD--MA0836.2", "JUNB--MA0490.2",
  "ATF4--MA0833.2","ATF7--MA0834.1")

add.flag <- function(pheatmap, kept.labels, repel.degree) {
  heatmap <- pheatmap$gtable
  new.label <- heatmap$grobs[[which(heatmap$layout$name == "row_names")]] 
  new.label$label <- ifelse(new.label$label %in% kept.labels, new.label$label, "")

  repelled.y <- function(d, d.select, k = repel.degree){
    # d = vector of distances for labels
    # d.select = vector of T/F for which labels are significant
    
    # recursive function to get current label positions
    # (note the unit is "npc" for all components of each distance)
    strip.npc <- function(dd){
      if(!"unit.arithmetic" %in% class(dd)) {
        return(as.numeric(dd))
      }
      
      d1 <- strip.npc(dd$arg1)
      d2 <- strip.npc(dd$arg2)
      fn <- dd$fname
      return(lazyeval::lazy_eval(paste(d1, fn, d2)))
    }
    
    full.range <- sapply(seq_along(d), function(i) strip.npc(d[i]))
    selected.range <- sapply(seq_along(d[d.select]), function(i) strip.npc(d[d.select][i]))
    
    return(unit(seq(from = max(selected.range) + k*(max(full.range) - max(selected.range)),
                    to = min(selected.range) - k*(min(selected.range) - min(full.range)), 
                    length.out = sum(d.select)), "npc"))
  }
  
  new.y.positions <- repelled.y(new.label$y, d.select = new.label$label != "")
  new.flag <- segmentsGrob(x0 = new.label$x,
                           x1 = new.label$x + unit(0.15, "npc"),
                           y0 = new.label$y[new.label$label != ""],
                           y1 = new.y.positions)

  new.label$x <- new.label$x + unit(0.2, "npc")
  new.label$y[new.label$label != ""] <- new.y.positions
  heatmap <- gtable::gtable_add_grob(x = heatmap, grobs = new.flag, t = 4, l = 4)
  heatmap$grobs[[which(heatmap$layout$name == "row_names")]] <- new.label
  grid.newpage()
  grid.draw(heatmap)
  invisible(heatmap)
}

plot <- add.flag(plot, kept.labels = labels, repel.degree = 0)
ggsave("../figures/Figure_3/Figure_3f_1.pdf", width = 6, height = 4, plot)


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
#   [1] grid      stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
#   [1] ggplot2_4.0.2      RColorBrewer_1.1-3 pheatmap_1.0.13    viridis_0.6.5      viridisLite_0.4.3 
# 
# loaded via a namespace (and not attached):
#   [1] deldir_2.0-4           pbapply_1.7-4          gridExtra_2.3          rlang_1.2.0            magrittr_2.0.5         RcppAnnoy_0.0.23      
# [7] otel_0.2.0             matrixStats_1.5.0      ggridges_0.5.7         compiler_4.5.2         spatstat.geom_3.7-3    systemfonts_1.3.2     
# [13] png_0.1-9              vctrs_0.7.2            reshape2_1.4.5         stringr_1.6.0          pkgconfig_2.0.3        fastmap_1.2.0         
# [19] promises_1.5.0         ragg_1.5.0             purrr_1.2.1            jsonlite_2.0.0         goftest_1.2-3          later_1.4.8           
# [25] spatstat.utils_3.2-2   irlba_2.3.7            parallel_4.5.2         cluster_2.1.8.1        R6_2.6.1               ica_1.0-3             
# [31] stringi_1.8.7          spatstat.data_3.1-9    reticulate_1.46.0      parallelly_1.46.1      spatstat.univar_3.1-7  lmtest_0.9-40         
# [37] scattermore_1.2        Rcpp_1.1.1             tensor_1.5.1           future.apply_1.20.2    zoo_1.8-15             sctransform_0.4.3     
# [43] httpuv_1.6.17          Matrix_1.7-4           splines_4.5.2          igraph_2.2.3           tidyselect_1.2.1       rstudioapi_0.18.0     
# [49] abind_1.4-8            spatstat.random_3.4-5  codetools_0.2-20       miniUI_0.1.2           spatstat.explore_3.8-0 listenv_0.10.1        
# [55] lattice_0.22-7         tibble_3.3.1           plyr_1.8.9             withr_3.0.2            shiny_1.13.0           S7_0.2.1              
# [61] ROCR_1.0-12            Rtsne_0.17             future_1.70.0          fastDummies_1.7.5      survival_3.8-3         polyclip_1.10-7       
# [67] fitdistrplus_1.2-6     pillar_1.11.1          Seurat_5.4.0           KernSmooth_2.23-26     plotly_4.12.0          generics_0.1.4        
# [73] RcppHNSW_0.6.0         sp_2.2-1               scales_1.4.0           globals_0.19.1         xtable_1.8-8           glue_1.8.0            
# [79] lazyeval_0.2.3         tools_4.5.2            data.table_1.18.2.1    RSpectra_0.16-2        RANN_2.6.2             dotCall64_1.2         
# [85] cowplot_1.2.0          tidyr_1.3.2            nlme_3.1-168           patchwork_1.3.2        cli_3.6.6              spatstat.sparse_3.1-0 
# [91] textshaping_1.0.5      spam_2.11-3            dplyr_1.2.1            uwot_0.2.4             gtable_0.3.6           digest_0.6.39         
# [97] progressr_0.19.0       ggrepel_0.9.8          htmlwidgets_1.6.4      SeuratObject_5.3.0     farver_2.1.2           htmltools_0.5.9       
# [103] lifecycle_1.0.5        httr_1.4.8             mime_0.13              MASS_7.3-65       