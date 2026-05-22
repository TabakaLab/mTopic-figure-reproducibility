library(RColorBrewer)
library(ggplot2)
library(ggalign)
library(data.table)
library(this.path)

setwd(this.dir())

dir.create("../figures/Figure_3", recursive = T)

AS <- data.frame(fread(file="../data/PRF1_scores.csv"))
rownames(AS) <- AS$V1
AS <- AS[,-1]
ColMeta <- fread(file="../data/PRF1_colmeta.csv")
RowMeta <- fread(file="../data/PRF1_rowmeta.csv")

min_val <- min(AS)
max_val <- max(AS)

cp <- colorRampPalette(colors = c("grey95", 
  brewer.pal(n=8, name='Greens')[4:8], "black"))(100)

color_peaks <- c("#CD0000", "#000080", "#00FFFF", "#006400", 
                 "purple4", "#CD853F", "#642915", "#58508d")

Top_peaks <- 1

vec_colors <- if(Top_peaks <= 8){
  color_peaks[1:Top_peaks]
} else {
  colorRampPalette(color_peaks)(Top_peaks)
}
if(length(vec_colors) < nrow(AS)){
  vec_colors <- c(vec_colors, rep("black", nrow(AS)-length(vec_colors)))
}

ASgg <- ggheatmap(AS) +
  theme(axis.text.y = element_text(angle = 0, hjust = 0, color = rev(vec_colors))) +
  scale_fill_gradientn(
    colours = cp, name = "ATAC score", limits = c(min_val, max_val),
    breaks = c(min_val, max_val), labels = c("Min", "Max"),
    guide = guide_colorbar(
      frame.colour = "black", frame.linewidth = 0.5, title.position = "left",
      label.position = "right", ticks.colour = "black", ticks.linewidth = 0.5,
      barwidth = unit(0.5, "cm"), barheight = unit(1.4, "cm"))
    ) +
  theme(
    panel.border = element_rect(fill = NA),
    axis.text.x = element_text(angle = 90, hjust = 1),
    legend.title = element_text(angle = 90),
    legend.position = c(-0.15, -0.08),
    legend.box.background = element_blank()
    ) + 
  anno_top() +
  ggalign(data = ColMeta, size = unit(1.5, "cm")) +
  geom_path(aes(y = Score,x = Rank), linewidth = 0.1) +
  geom_point(aes(y = Score, x = Rank, color = color), size = 1.5) + 
  scale_color_identity() +
  theme(panel.border = element_rect(fill = NA)) +
  labs(y = "RNA score") + 
  scale_y_log10() +
  anno_right() +
  ggalign(data=RowMeta, size = unit(1.5, "cm")) +
  geom_point(aes(y = Rows, x = HHI)) +
  geom_path(aes(y = Rows, x = HHI), linewidth=0.1) +
  labs(x = "HHI") +
  scale_x_continuous(expand = expansion(mult = c(0.1, 0.1))) +
  theme(
    panel.border = element_rect(fill = NA),
    axis.text.x = element_text(angle = 90, hjust = 1)
  ) +
  anno_right() +
  ggalign(data = RowMeta,size = unit(1.5, "cm")) +
  geom_point(aes(y = Rows, x = ModelScores))+
  geom_path(aes(y = Rows, x = ModelScores), linewidth = 0.1) +
  labs(x = "Global \n association \n score") +
  scale_x_continuous(expand = expansion(mult = c(0.1, 0.1)))+
  theme(
    panel.border = element_rect(fill = NA),
    axis.text.x = element_text(angle = 90, hjust = 1)
  )

pdf(file = paste0("../figures/Figure_3/Figure_3i.pdf"), width = 10, height = 8)
  ASgg
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
#   [1] ggalign_1.2.0       RColorBrewer_1.1-3  dplyr_1.2.1         data.table_1.18.2.1 ggrepel_0.9.8       ggplot2_4.0.2      
# 
# loaded via a namespace (and not attached):
#   [1] vctrs_0.7.2       cli_3.6.6         rlang_1.2.0       generics_0.1.4    textshaping_1.0.5 S7_0.2.1          glue_1.8.0       
# [8] labeling_0.4.3    ragg_1.5.0        scales_1.4.0      grid_4.5.2        tibble_3.3.1      lifecycle_1.0.5   compiler_4.5.2   
# [15] Rcpp_1.1.1        pkgconfig_2.0.3   rstudioapi_0.18.0 systemfonts_1.3.2 farver_2.1.2      R6_2.6.1          tidyselect_1.2.1 
# [22] pillar_1.11.1     magrittr_2.0.5    tools_4.5.2       withr_3.0.2       gtable_0.3.6    