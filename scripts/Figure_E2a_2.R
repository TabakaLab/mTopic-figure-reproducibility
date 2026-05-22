library(Seurat)
library(dplyr)
library(stringr)
library(ggplot2)
library(ggrepel)
library(tibble)
library(this.path)

setwd(this.dir())

dir.create("../figures/Figure_E2", recursive = T)

MouseAdoAtlas <- readRDS("../data/MouseAdolescentBrainAtlas.RDS")

get_prefix <- function(x) {
  stringr::str_extract(x, "^[A-Za-z]+")
}

get_index <- function(x) {
  suppressWarnings(as.integer(stringr::str_extract(x, "[0-9]+$")))
}

zoom_umap_clusters <- function(
    seurat_obj,
    target_clusters,
    reduction = "umap",
    cluster_col = "ClusterName",
    x = 2.5,
    point_size = 0.35,
    target_point_size = NULL,
    point_alpha = 0.95,
    other_alpha = 0.35,
    label_size = 3.8,
    repel_force = 1,
    show_box_on_full_umap = FALSE,
    full_plot_color_all = TRUE,   # NEW: color all clusters on full plot
    show_legend = TRUE,
    color_other_cells = FALSE,
    other_cell_color = "grey85",
    palette = c(
      "#cdba96", "#000080", "#ffff00", "#9ecae1", "#189e7f", "#fff8dc",
      "#642915", "#0000ff", "#bebebe", "#4292c6", "#cd0000", "#58508d",
      "#006400", "#10d894", "#666666", "#800080", "#ff00ff", "#ffc0cb",
      "#bebebe", "#ffd700", "#bebebe", "#00ffff", "#cd853f", "#dda0dd",
      "#000000", "#551a8b", "#00cd00", "#f4a460", "#00ff7f", "#ffbbff"
    ),
    cluster_colors = NULL,   # optional named vector
    seed_use = 10
) {
  set.seed(seed_use)
  
  if (is.null(target_point_size)) {
    target_point_size <- point_size
  }
  
  stopifnot(cluster_col %in% colnames(seurat_obj@meta.data))
  stopifnot(reduction %in% names(seurat_obj@reductions))
  
  emb <- Seurat::Embeddings(seurat_obj, reduction = reduction)[, 1:2, drop = FALSE]
  
  df <- tibble::tibble(
    cell        = rownames(emb),
    UMAP_1      = emb[, 1],
    UMAP_2      = emb[, 2],
    ClusterName = as.character(seurat_obj@meta.data[rownames(emb), cluster_col])
  )
  
  missing_targets <- setdiff(target_clusters, unique(df$ClusterName))
  if (length(missing_targets) > 0) {
    stop(
      paste0(
        "These target clusters were not found in ", cluster_col, ": ",
        paste(missing_targets, collapse = ", ")
      )
    )
  }
  
  cluster_levels <- df %>%
    dplyr::distinct(ClusterName) %>%
    dplyr::mutate(
      prefix = get_prefix(ClusterName),
      idx    = get_index(ClusterName)
    ) %>%
    dplyr::arrange(prefix, is.na(idx), idx, ClusterName) %>%
    dplyr::pull(ClusterName)
  
  palette <- unique(palette)
  
  if (!is.null(cluster_colors)) {
    if (is.null(names(cluster_colors))) {
      if (length(cluster_colors) < length(cluster_levels)) {
        cluster_cols <- grDevices::colorRampPalette(palette)(length(cluster_levels))
      } else {
        cluster_cols <- cluster_colors[seq_along(cluster_levels)]
      }
      names(cluster_cols) <- cluster_levels
    } else {
      cluster_cols <- cluster_colors
      missing_color_clusters <- setdiff(cluster_levels, names(cluster_cols))
      
      if (length(missing_color_clusters) > 0) {
        used_cols <- unname(cluster_cols)
        available_palette <- setdiff(palette, used_cols)
        
        if (length(available_palette) < length(missing_color_clusters)) {
          fill_cols <- grDevices::colorRampPalette(palette)(length(missing_color_clusters))
        } else {
          fill_cols <- available_palette[seq_along(missing_color_clusters)]
        }
        
        names(fill_cols) <- missing_color_clusters
        cluster_cols <- c(cluster_cols, fill_cols)
      }
      
      cluster_cols <- cluster_cols[cluster_levels]
    }
  } else {
    if (length(cluster_levels) <= length(palette)) {
      cluster_cols <- palette[seq_along(cluster_levels)]
    } else {
      cluster_cols <- grDevices::colorRampPalette(palette)(length(cluster_levels))
    }
    names(cluster_cols) <- cluster_levels
  }
  
  target_df_all <- df %>%
    dplyr::filter(ClusterName %in% target_clusters)
  
  center_x <- median(target_df_all$UMAP_1)
  center_y <- median(target_df_all$UMAP_2)
  
  x_min <- center_x - x
  x_max <- center_x + x
  y_min <- center_y - x
  y_max <- center_y + x
  
  zoom_df <- df %>%
    dplyr::filter(
      UMAP_1 >= x_min, UMAP_1 <= x_max,
      UMAP_2 >= y_min, UMAP_2 <= y_max
    ) %>%
    dplyr::mutate(is_target = ClusterName %in% target_clusters)
  
  if (nrow(zoom_df) == 0) {
    stop("No cells found inside the zoom box. Try increasing x.")
  }
  
  target_zoom_df <- zoom_df %>%
    dplyr::filter(is_target)
  
  other_zoom_df <- zoom_df %>%
    dplyr::filter(!is_target)
  
  label_df <- target_zoom_df %>%
    dplyr::group_by(ClusterName) %>%
    dplyr::summarise(
      UMAP_1 = median(UMAP_1),
      UMAP_2 = median(UMAP_2),
      .groups = "drop"
    )
  
  target_clusters_present <- target_clusters[target_clusters %in% unique(target_zoom_df$ClusterName)]
  target_zoom_df$ClusterName <- factor(target_zoom_df$ClusterName, levels = target_clusters_present)
  label_df$ClusterName <- factor(label_df$ClusterName, levels = target_clusters_present)
  
  p_zoom <- ggplot()
  
  if (nrow(other_zoom_df) > 0) {
    if (color_other_cells) {
      other_clusters_visible <- cluster_levels[cluster_levels %in% unique(other_zoom_df$ClusterName)]
      other_zoom_df$ClusterName <- factor(other_zoom_df$ClusterName, levels = other_clusters_visible)
      
      p_zoom <- p_zoom +
        geom_point(
          data = other_zoom_df,
          aes(x = UMAP_1, y = UMAP_2, color = ClusterName),
          size = point_size,
          alpha = other_alpha,
          stroke = 0,
          show.legend = FALSE
        )
    } else {
      p_zoom <- p_zoom +
        geom_point(
          data = other_zoom_df,
          aes(x = UMAP_1, y = UMAP_2),
          color = other_cell_color,
          size = point_size,
          alpha = other_alpha,
          stroke = 0,
          show.legend = FALSE
        )
    }
  }
  
  p_zoom <- p_zoom +
    geom_point(
      data = target_zoom_df,
      aes(x = UMAP_1, y = UMAP_2, color = ClusterName),
      size = target_point_size,
      alpha = point_alpha,
      stroke = 0
    ) +
    ggrepel::geom_text_repel(
      data = label_df,
      aes(x = UMAP_1, y = UMAP_2, label = ClusterName),
      inherit.aes = FALSE,
      color = "black",
      size = label_size,
      fontface = "bold",
      seed = seed_use,
      force = repel_force,
      box.padding = 0.25,
      point.padding = 0.15,
      min.segment.length = 0,
      max.overlaps = Inf,
      segment.size = 0.25,
      segment.alpha = 0.7
    ) +
    scale_color_manual(
      values = cluster_cols,
      breaks = target_clusters_present,
      limits = names(cluster_cols),
      drop = FALSE
    ) +
    coord_cartesian(
      xlim = c(x_min, x_max),
      ylim = c(y_min, y_max),
      expand = FALSE
    ) +
    labs(
      x = "UMAP_1",
      y = "UMAP_2",
      title = paste0("Zoom: ", paste(target_clusters, collapse = ", ")),
      color = NULL
    ) +
    guides(
      color = guide_legend(
        override.aes = list(size = 3.5, alpha = 1)
      )
    ) +
    theme_classic(base_size = 11) +
    theme(
      legend.position = if (show_legend) "right" else "none",
      legend.text = element_text(size = 10, color = "black"),
      axis.title = element_text(size = 12, color = "black"),
      axis.text  = element_text(size = 10, color = "black"),
      axis.line  = element_line(linewidth = 0.5, color = "black"),
      plot.title = element_text(size = 12, face = "bold"),
      plot.margin = margin(4, 6, 4, 4, "mm")
    )
  
  p_full <- NULL
  if (show_box_on_full_umap) {
    df$ClusterName <- factor(df$ClusterName, levels = cluster_levels)
    
    if (full_plot_color_all) {
      p_full <- ggplot(df, aes(x = UMAP_1, y = UMAP_2, color = ClusterName)) +
        geom_point(size = 0.20, alpha = 0.85, stroke = 0) +
        geom_rect(
          inherit.aes = FALSE,
          aes(xmin = x_min, xmax = x_max, ymin = y_min, ymax = y_max),
          fill = NA,
          color = "black",
          linewidth = 0.5
        ) +
        scale_color_manual(
          values = cluster_cols,
          limits = names(cluster_cols),
          drop = FALSE
        ) +
        labs(
          x = "UMAP_1",
          y = "UMAP_2",
          title = "Full UMAP with zoom box",
          color = NULL
        ) +
        theme_classic(base_size = 11) +
        theme(
          legend.position = if (show_legend) "right" else "none",
          legend.text = element_text(size = 10, color = "black"),
          axis.title = element_text(size = 12, color = "black"),
          axis.text  = element_text(size = 10, color = "black"),
          axis.line  = element_line(linewidth = 0.5, color = "black"),
          plot.title = element_text(size = 12, face = "bold"),
          plot.margin = margin(4, 6, 4, 4, "mm")
        )
    } else {
      full_target_df <- df %>%
        dplyr::filter(ClusterName %in% target_clusters)
      
      p_full <- ggplot(df, aes(x = UMAP_1, y = UMAP_2)) +
        geom_point(color = "grey85", size = 0.15, alpha = 0.7, stroke = 0) +
        geom_point(
          data = full_target_df,
          aes(color = ClusterName),
          size = 0.35,
          alpha = 0.95,
          stroke = 0
        ) +
        geom_rect(
          inherit.aes = FALSE,
          aes(xmin = x_min, xmax = x_max, ymin = y_min, ymax = y_max),
          fill = NA,
          color = "black",
          linewidth = 0.5
        ) +
        scale_color_manual(
          values = cluster_cols,
          breaks = target_clusters,
          limits = names(cluster_cols),
          drop = FALSE
        ) +
        labs(
          x = "UMAP_1",
          y = "UMAP_2",
          title = "Full UMAP with zoom box",
          color = NULL
        ) +
        theme_classic(base_size = 11) +
        theme(
          legend.position = if (show_legend) "right" else "none",
          legend.text = element_text(size = 10, color = "black"),
          axis.title = element_text(size = 12, color = "black"),
          axis.text  = element_text(size = 10, color = "black"),
          axis.line  = element_line(linewidth = 0.5, color = "black"),
          plot.title = element_text(size = 12, face = "bold"),
          plot.margin = margin(4, 6, 4, 4, "mm")
        )
    }
  }
  
  return(list(
    plot = p_zoom,
    full_plot = p_full,
    zoom_data = zoom_df,
    target_data = target_zoom_df,
    other_data = other_zoom_df,
    label_data = label_df,
    center = c(UMAP_1 = center_x, UMAP_2 = center_y),
    box = c(xmin = x_min, xmax = x_max, ymin = y_min, ymax = y_max),
    colors = cluster_cols
  ))
}

my_cluster_colors <- c(
  TEGLU1  = "#cdba96",
  TEGLU2  = "#642915",
  TEGLU3  = "#ffff00",
  TEGLU4  = "#9ecae1",
  TEGLU5  = "#189e7f",
  TEGLU6  = "#fff8dc",
  TEGLU7  = "#642915",
  TEGLU8  = "#0000ff",
  TEGLU9  = "#bebebe",
  TEGLU10 = "#4292c6",
  TEGLU11 = "#cd0000",
  TEGLU12 = "#58508d",
  TEGLU13 = "#006400",
  TEGLU14 = "#10d894",
  TEGLU15 = "#666666",
  TEGLU16 = "#800080",
  TEGLU17 = "#ff00ff",
  TEGLU18 = "#ffc0cb",
  TEGLU19 = "#ffd700",
  TEGLU20 = "#00ffff",
  TEGLU21 = "#cd853f",
  TEGLU22 = "#dda0dd"
)


res_teglu23_custom <- zoom_umap_clusters(
  seurat_obj = MouseAdoAtlas,
  target_clusters = c("TEGLU4", paste0("TEGLU",5:11), paste0("TEGLU", 13:22)),
  x = 5,
  point_size = 0.7, 
  label_size = 4,
  repel_force = 100,
  cluster_colors = my_cluster_colors,
  show_box_on_full_umap = TRUE,
  full_plot_color_all = FALSE,
  show_legend = TRUE
)

res_teglu23_custom$plot
res_teglu23_custom$full_plot

ggsave(
  filename = "../figures/Figure_E2/Figure_E2a_2.pdf",
  plot = res_teglu23_custom$plot,
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
#   [1] tibble_3.3.1        ggrepel_0.9.8       stringr_1.6.0       Seurat_5.4.0        SeuratObject_5.3.0  sp_2.2-1            colorRamps_2.3.4   
# [8] viridis_0.6.5       viridisLite_0.4.3   dplyr_1.2.1         data.table_1.18.2.1 rgexf_0.16.3        ggraph_2.2.2        ggplot2_4.0.2      
# [15] tidygraph_1.3.1     igraph_2.2.3       
# 
# loaded via a namespace (and not attached):
#   [1] RColorBrewer_1.1-3     rstudioapi_0.18.0      jsonlite_2.0.0         magrittr_2.0.5         spatstat.utils_3.2-2   farver_2.1.2          
# [7] ragg_1.5.0             vctrs_0.7.2            ROCR_1.0-12            memoise_2.0.1          spatstat.explore_3.8-0 htmltools_0.5.9       
# [13] sctransform_0.4.3      parallelly_1.46.1      KernSmooth_2.23-26     htmlwidgets_1.6.4      ica_1.0-3              plyr_1.8.9            
# [19] plotly_4.12.0          zoo_1.8-15             cachem_1.1.0           mime_0.13              lifecycle_1.0.5        pkgconfig_2.0.3       
# [25] Matrix_1.7-4           R6_2.6.1               fastmap_1.2.0          fitdistrplus_1.2-6     future_1.70.0          shiny_1.13.0          
# [31] digest_0.6.39          patchwork_1.3.2        tensor_1.5.1           RSpectra_0.16-2        irlba_2.3.7            textshaping_1.0.5     
# [37] labeling_0.4.3         progressr_0.19.0       spatstat.sparse_3.1-0  httr_1.4.8             polyclip_1.10-7        abind_1.4-8           
# [43] compiler_4.5.2         withr_3.0.2            S7_0.2.1               fastDummies_1.7.5      ggforce_0.5.0          MASS_7.3-65           
# [49] tools_4.5.2            lmtest_0.9-40          otel_0.2.0             httpuv_1.6.17          future.apply_1.20.2    goftest_1.2-3         
# [55] glue_1.8.0             nlme_3.1-168           promises_1.5.0         grid_4.5.2             Rtsne_0.17             cluster_2.1.8.1       
# [61] reshape2_1.4.5         generics_0.1.4         gtable_0.3.6           spatstat.data_3.1-9    tidyr_1.3.2            spatstat.geom_3.7-3   
# [67] RcppAnnoy_0.0.23       RANN_2.6.2             pillar_1.11.1          spam_2.11-3            RcppHNSW_0.6.0         servr_0.32            
# [73] later_1.4.8            splines_4.5.2          tweenr_2.0.3           lattice_0.22-7         survival_3.8-3         deldir_2.0-4          
# [79] tidyselect_1.2.1       miniUI_0.1.2           pbapply_1.7-4          gridExtra_2.3          scattermore_1.2        xfun_0.57             
# [85] graphlayouts_1.2.2     matrixStats_1.5.0      stringi_1.8.7          lazyeval_0.2.3         codetools_0.2-20       cli_3.6.6             
# [91] uwot_0.2.4             xtable_1.8-8           reticulate_1.46.0      systemfonts_1.3.2      Rcpp_1.1.1             globals_0.19.1        
# [97] spatstat.random_3.4-5  png_0.1-9              XML_3.99-0.22          spatstat.univar_3.1-7  parallel_4.5.2         dotCall64_1.2         
# [103] listenv_0.10.1         scales_1.4.0           ggridges_0.5.7         purrr_1.2.1            rlang_1.2.0            cowplot_1.2.0   

