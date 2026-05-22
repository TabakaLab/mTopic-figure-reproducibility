library(pheatmap)
library(data.table)
library(RColorBrewer)
library(grid)
library(ggplot2)
library(gtable)
library(ggplotify)
library(this.path)

setwd(this.dir())

dir.create("../figures/Figure_E2", recursive = T)

P22_TFs <- readRDS(file = "../data/TFMotifs_P22atac.RDS")
P22_TFs <- t(P22_TFs)
labels <- unique(fread("../data/P22_atac_TF_scores.csv")$Motif)

plot <- pheatmap(
  P22_TFs,
  legend = TRUE,
  cluster_cols = TRUE,
  fontsize_row = 9,
  fontsize_col = 8,
  angle_col = "45",
  color = colorRampPalette(
    colors = c(
      rev(brewer.pal(n = 9, name = "Blues")),
      "white", "white", "white",
      brewer.pal(n = 8, name = "Greens")
    )
  )(50)
)

add.flag.col <- function(pheatmap, kept.labels, repel.degree = 0) {
  heatmap <- pheatmap$gtable
  idx <- which(heatmap$layout$name == "col_names")
  new.label <- heatmap$grobs[[idx]]
  keep <- new.label$label %in% kept.labels
  new.label$label <- ifelse(keep, new.label$label, "")
  
  if (!any(keep)) {
    warning("No column labels matched kept.labels. Returning heatmap unchanged.")
    grid.newpage()
    grid.draw(heatmap)
    return(invisible(heatmap))
  }
  
  repelled.x <- function(d, d.select, k = repel.degree) {
    strip.npc <- function(dd) {
      if (!"unit.arithmetic" %in% class(dd)) {
        return(as.numeric(dd))
      }
      d1 <- strip.npc(dd$arg1)
      d2 <- strip.npc(dd$arg2)
      fn <- dd$fname
      if (fn == "+") return(d1 + d2)
      if (fn == "-") return(d1 - d2)
      if (fn == "*") return(d1 * d2)
      if (fn == "/") return(d1 / d2)
      stop("Unsupported unit arithmetic")
    }
    
    full.range <- sapply(seq_along(d), function(i) strip.npc(d[i]))
    selected.range <- sapply(seq_along(d[d.select]), function(i) strip.npc(d[d.select][i]))
    
    unit(
      seq(
        from = min(selected.range) - k * (min(selected.range) - min(full.range)),
        to   = max(selected.range) + k * (max(full.range) - max(selected.range)),
        length.out = sum(d.select)
      ),
      "npc"
    )
  }
  new.x.positions <- repelled.x(new.label$x, keep, repel.degree)
  new.flag <- segmentsGrob(
    x0 = new.label$x[keep],
    x1 = new.x.positions,
    y0 = new.label$y,
    y1 = new.label$y - unit(0.15, "npc")
  )
  new.label$x[keep] <- new.x.positions
  new.label$y <- new.label$y - unit(0.2, "npc")
  heatmap <- gtable_add_grob(
    x = heatmap,
    grobs = new.flag,
    t = 5,
    l = 3
  )
  heatmap$grobs[[idx]] <- new.label
  grid.newpage()
  grid.draw(heatmap)
  invisible(heatmap)
}

plot2 <- add.flag.col(plot, kept.labels = labels, repel.degree = 0)
plot2$layout$clip <- "off"
p_save <- ggplotify::as.ggplot(plot2, scale = 0.92) +
  coord_cartesian(clip = "off") +
  theme(plot.margin = margin(t = 20, r = 10, b = 70, l = 10, unit = "pt"))

ggsave("../figures/Figure_E2/Figure_E2c.pdf", plot = p_save, width = 7, height = 5, 
       units = "in", limitsize = FALSE)

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
#   [1] ggplotify_0.1.3     gtable_0.3.6        ggplot2_4.0.2       RColorBrewer_1.1-3  data.table_1.18.2.1 pheatmap_1.0.13    
# 
# loaded via a namespace (and not attached):
#   [1] vctrs_0.7.2        cli_3.6.6          rlang_1.2.0        gridGraphics_0.5-1 generics_0.1.4     textshaping_1.0.5  S7_0.2.1          
# [8] labeling_0.4.3     glue_1.8.0         ragg_1.5.0         scales_1.4.0       rappdirs_0.3.4     tibble_3.3.1       lifecycle_1.0.5   
# [15] compiler_4.5.2     dplyr_1.2.1        fs_2.0.1           pkgconfig_2.0.3    rstudioapi_0.18.0  systemfonts_1.3.2  farver_2.1.2      
# [22] digest_0.6.39      R6_2.6.1           tidyselect_1.2.1   pillar_1.11.1      magrittr_2.0.5     tools_4.5.2        withr_3.0.2       
# [29] yulab.utils_0.2.4 