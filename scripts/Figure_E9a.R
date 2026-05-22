library(pheatmap)
library(RColorBrewer)
library(data.table)
library(this.path)

setwd(this.dir())

dir.create("../figures/Figure_E9", recursive = T)

dictionary <- data.frame(fread("../data/cell_type_annotations_PBMCs.csv", header=T))
data <- readRDS(file="../data/PeakGene_associations.RDS")
row_anno <- data.frame(RowCat = rownames(data), row.names = rownames(data))
col_anno <- data.frame(ColCat = colnames(data), row.names = colnames(data))
unique_cols <- unique(col_anno$ColCat)
col_palette <- dictionary$colors
names(col_palette) <- unique_cols
ann_colors <- list(ColCat = col_palette)

col_palette <- setNames(dictionary$colors, dictionary$cell_type)
col_palette <- col_palette[unique_cols]
ann_colors <- list(ColCat = col_palette)

ph <- pheatmap(
  as.matrix(data),cluster_rows = F,cluster_cols = F, border_color = NA,
  show_rownames = F,
  annotation_col = col_anno,
  annotation_colors = ann_colors,
  fontsize_row = 6, annotation_names_row = F, fontsize = 6,
  color = colorRampPalette(colors = c(brewer.pal(n = 8, name = 'Greens'),"black"))(100))

png(file = paste0("../figures/Figure_E9/Figure_E9a.png"), width = 180, height=120, units="mm", res =300)
ph
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
#   [1] stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
#   [1] data.table_1.18.2.1 RColorBrewer_1.1-3  pheatmap_1.0.13    
# 
# loaded via a namespace (and not attached):
#   [1] scales_1.4.0      compiler_4.5.2    R6_2.6.1          cli_3.6.6         tools_4.5.2       glue_1.8.0        gtable_0.3.6     
# [8] rstudioapi_0.18.0 farver_2.1.2      grid_4.5.2        pkgload_1.5.1     lifecycle_1.0.5   rlang_1.2.0      