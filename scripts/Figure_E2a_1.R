library(Seurat)
library(dplyr)
library(stringr)
library(purrr)
library(igraph)
library(ggplot2)
library(ggrepel)
library(tibble)
library(this.path)

setwd(this.dir())

dir.create("../figures/Figure_E2", recursive = T)

MouseAdoAtlas <- readRDS("../data/MouseAdolescentBrainAtlas.RDS")

merge_gap   <- 0.3   # larger = more label merging
point_size  <- 0.25
point_alpha <- 0.90
label_size  <- 3.0
seed_use    <- 10

set.seed(seed_use)
stopifnot("ClusterName" %in% colnames(MouseAdoAtlas@meta.data))
stopifnot("umap" %in% names(MouseAdoAtlas@reductions))

umap_mat <- Embeddings(MouseAdoAtlas, reduction = "umap")[, 1:2, drop = FALSE]

plot_df <- tibble(
  cell        = rownames(umap_mat),
  UMAP_1      = umap_mat[, 1],
  UMAP_2      = umap_mat[, 2],
  ClusterName = as.character(MouseAdoAtlas@meta.data[rownames(umap_mat), "ClusterName"])
)

## ----------------------------
## Helper functions
## ----------------------------
get_prefix <- function(x) {
  stringr::str_extract(x, "^[A-Za-z]+")
}

get_index <- function(x) {
  suppressWarnings(as.integer(stringr::str_extract(x, "[0-9]+$")))
}

make_range_label <- function(prefix, idx) {
  idx <- sort(unique(idx[!is.na(idx)]))
  if (length(idx) == 0) return(prefix)
  grp <- cumsum(c(1, diff(idx) != 1))
  runs <- split(idx, grp)
  run_labels <- vapply(
    runs,
    function(v) {
      if (length(v) == 1) {
        as.character(v)
      } else {
        paste0(min(v), "-", max(v))
      }
    },
    character(1)
  )
  paste0(prefix, paste(run_labels, collapse = ","))
}

merge_nearby_labels <- function(df_prefix, gap_cut = 0.35) {
  df_prefix <- df_prefix[order(df_prefix$idx), , drop = FALSE]
  
  if (nrow(df_prefix) == 1) {
    df_prefix$DisplayLabel <- make_range_label(df_prefix$LabelPrefix[[1]], df_prefix$idx)
    return(df_prefix)
  }
  
  coords <- as.matrix(df_prefix[, c("cx", "cy"), drop = FALSE])
  rownames(coords) <- df_prefix$ClusterName
  
  centroid_dist <- as.matrix(stats::dist(coords))
  radius_sum    <- outer(df_prefix$radius, df_prefix$radius, "+")
  
  ## negative gap = overlapping clusters
  gap <- centroid_dist - radius_sum
  
  ## only allow consecutive numbered labels to merge
  consecutive <- abs(outer(df_prefix$idx, df_prefix$idx, "-")) == 1
  
  adj <- (gap <= gap_cut) & consecutive
  diag(adj) <- FALSE
  
  g <- igraph::graph_from_adjacency_matrix(adj, mode = "undirected", diag = FALSE)
  comp <- igraph::components(g)$membership
  
  df_prefix$component <- comp[match(df_prefix$ClusterName, names(comp))]
  
  df_prefix %>%
    dplyr::group_by(component) %>%
    dplyr::mutate(DisplayLabel = make_range_label(LabelPrefix[[1]], idx)) %>%
    dplyr::ungroup() %>%
    dplyr::select(-component)
}

## ----------------------------
## Per-cluster centroid and radius
## ----------------------------
cluster_summary <- plot_df %>%
  dplyr::group_by(ClusterName) %>%
  dplyr::summarise(
    cx = median(UMAP_1),
    cy = median(UMAP_2),
    n  = dplyr::n(),
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    LabelPrefix = get_prefix(ClusterName),
    idx         = get_index(ClusterName)
  )

plot_df2 <- plot_df %>%
  dplyr::left_join(
    cluster_summary %>% dplyr::select(ClusterName, cx, cy),
    by = "ClusterName"
  ) %>%
  dplyr::mutate(
    dist_to_centroid = sqrt((UMAP_1 - cx)^2 + (UMAP_2 - cy)^2)
  )

cluster_radius <- plot_df2 %>%
  dplyr::group_by(ClusterName) %>%
  dplyr::summarise(
    radius = as.numeric(stats::quantile(dist_to_centroid, probs = 0.90)),
    .groups = "drop"
  )

cluster_summary <- cluster_summary %>%
  dplyr::left_join(cluster_radius, by = "ClusterName")

## ----------------------------
## Build merged display labels
## ----------------------------
numbered_clusters <- cluster_summary %>%
  dplyr::filter(!is.na(idx))

non_numbered_clusters <- cluster_summary %>%
  dplyr::filter(is.na(idx)) %>%
  dplyr::mutate(DisplayLabel = ClusterName)

merged_numbered <- split(numbered_clusters, numbered_clusters$LabelPrefix) %>%
  purrr::map_dfr(merge_nearby_labels, gap_cut = merge_gap)

cluster_labels <- dplyr::bind_rows(
  merged_numbered,
  non_numbered_clusters
) %>%
  dplyr::select(ClusterName, LabelPrefix, idx, DisplayLabel)

plot_df_final <- plot_df %>%
  dplyr::left_join(cluster_labels, by = "ClusterName")

## ----------------------------
## Label positions
## ----------------------------
label_df <- plot_df_final %>%
  dplyr::group_by(DisplayLabel) %>%
  dplyr::summarise(
    UMAP_1 = median(UMAP_1),
    UMAP_2 = median(UMAP_2),
    .groups = "drop"
  )

## ----------------------------
## User palette
## ----------------------------
base_palette <- unique(c(
  "#cdba96", "#000080", "#ffff00", "#9ecae1", "#189e7f", "#fff8dc",
  "#642915", "#0000ff", "#bebebe", "#4292c6", "#cd0000", "#58508d",
  "#006400", "#10d894", "#666666", "#800080", "#ff00ff", "#ffc0cb",
  "#ffd700", "#00ffff", "#cd853f", "#dda0dd", "#000000", "#551a8b",
  "#00cd00", "#f4a460", "#00ff7f", "#ffbbff"
))

display_groups <- sort(unique(plot_df_final$DisplayLabel))

if (length(display_groups) <= length(base_palette)) {
  group_cols <- base_palette[seq_along(display_groups)]
} else {
  group_cols <- grDevices::colorRampPalette(base_palette)(length(display_groups))
}
names(group_cols) <- display_groups

plot_df_final$DisplayLabel <- factor(plot_df_final$DisplayLabel, levels = display_groups)
label_df$DisplayLabel <- factor(label_df$DisplayLabel, levels = display_groups)

## ----------------------------
## Plot
## ----------------------------
p_pub <- ggplot(plot_df_final, aes(x = UMAP_1, y = UMAP_2, color = DisplayLabel)) +
  geom_point(size = point_size, alpha = point_alpha, stroke = 0) +
  ggrepel::geom_text_repel(
    data = label_df,
    aes(x = UMAP_1, y = UMAP_2, label = DisplayLabel),
    inherit.aes = FALSE,
    color = "black",
    size = label_size,
    fontface = "bold",
    seed = seed_use,
    box.padding = 0.25,
    point.padding = 0.15,
    min.segment.length = 0,
    max.overlaps = Inf,
    segment.size = 0.25,
    segment.alpha = 0.7
  ) +
  scale_color_manual(values = group_cols, drop = FALSE) +
  coord_equal() +
  labs(
    x = "UMAP_1",
    y = "UMAP_2",
    color = NULL
  ) +
  guides(
    color = guide_legend(
      override.aes = list(size = 3, alpha = 1),
      ncol = 1
    )
  ) +
  theme_classic(base_size = 11) +
  theme(
    legend.position = "none",
    legend.text = element_text(size = 9, color = "black"),
    axis.title = element_text(size = 12, color = "black"),
    axis.text  = element_text(size = 10, color = "black"),
    axis.line  = element_line(linewidth = 0.5, color = "black"),
    plot.margin = margin(4, 6, 4, 4, "mm")
  )

p_pub

ggsave(
  filename = "../figures/Figure_E2/Figure_E2a_1.pdf",
  plot = p_pub,
  width = 6.2,
  height = 5.0,
  units = "in",
  useDingbats = FALSE
)


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
#   [1] tibble_3.3.1       ggrepel_0.9.8      ggplot2_4.0.2      igraph_2.2.3       purrr_1.2.1        stringr_1.6.0      dplyr_1.2.1       
# [8] Seurat_5.4.0       SeuratObject_5.3.0 sp_2.2-1          
# 
# loaded via a namespace (and not attached):
#   [1] deldir_2.0-4           pbapply_1.7-4          gridExtra_2.3          rlang_1.2.0            magrittr_2.0.5         RcppAnnoy_0.0.23      
# [7] otel_0.2.0             matrixStats_1.5.0      ggridges_0.5.7         compiler_4.5.2         spatstat.geom_3.7-3    systemfonts_1.3.2     
# [13] png_0.1-9              vctrs_0.7.2            reshape2_1.4.5         pkgconfig_2.0.3        fastmap_1.2.0          labeling_0.4.3        
# [19] promises_1.5.0         ragg_1.5.0             jsonlite_2.0.0         goftest_1.2-3          later_1.4.8            spatstat.utils_3.2-2  
# [25] irlba_2.3.7            parallel_4.5.2         cluster_2.1.8.1        R6_2.6.1               ica_1.0-3              stringi_1.8.7         
# [31] RColorBrewer_1.1-3     spatstat.data_3.1-9    reticulate_1.46.0      parallelly_1.46.1      spatstat.univar_3.1-7  lmtest_0.9-40         
# [37] scattermore_1.2        Rcpp_1.1.1             tensor_1.5.1           future.apply_1.20.2    zoo_1.8-15             sctransform_0.4.3     
# [43] httpuv_1.6.17          Matrix_1.7-4           splines_4.5.2          tidyselect_1.2.1       rstudioapi_0.18.0      abind_1.4-8           
# [49] spatstat.random_3.4-5  codetools_0.2-20       miniUI_0.1.2           spatstat.explore_3.8-0 listenv_0.10.1         lattice_0.22-7        
# [55] plyr_1.8.9             withr_3.0.2            shiny_1.13.0           S7_0.2.1               ROCR_1.0-12            Rtsne_0.17            
# [61] future_1.70.0          fastDummies_1.7.5      survival_3.8-3         polyclip_1.10-7        fitdistrplus_1.2-6     pillar_1.11.1         
# [67] KernSmooth_2.23-26     plotly_4.12.0          generics_0.1.4         RcppHNSW_0.6.0         scales_1.4.0           globals_0.19.1        
# [73] xtable_1.8-8           glue_1.8.0             lazyeval_0.2.3         tools_4.5.2            data.table_1.18.2.1    RSpectra_0.16-2       
# [79] RANN_2.6.2             dotCall64_1.2          cowplot_1.2.0          grid_4.5.2             tidyr_1.3.2            nlme_3.1-168          
# [85] patchwork_1.3.2        cli_3.6.6              spatstat.sparse_3.1-0  textshaping_1.0.5      spam_2.11-3            viridisLite_0.4.3     
# [91] uwot_0.2.4             gtable_0.3.6           digest_0.6.39          progressr_0.19.0       htmlwidgets_1.6.4      farver_2.1.2          
# [97] htmltools_0.5.9        lifecycle_1.0.5        httr_1.4.8             mime_0.13              MASS_7.3-65   
