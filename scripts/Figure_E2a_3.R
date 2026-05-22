library(Seurat)
library(dplyr)
library(stringr)
library(ggplot2)
library(ggrepel)
library(tibble)
library(MASS)
library(grid)
library(viridis)
library(data.table)
library(patchwork)
library(this.path)

setwd(this.dir())

dir.create("../figures/Figure_E2", recursive = T)

MouseAdoAtlas <- readRDS("../data/MouseAdolescentBrainAtlas.RDS")

order_topic <- function(signatures, topic=1) {
  topic_col <- paste0("topic_", topic)
  ind <- order(signatures[[topic_col]], decreasing = TRUE)
  df <- data.frame(
    Feature = signatures$V1[ind], 
    Score = signatures[[topic_col]][ind]
  )
  return(df)
}

add_modality_signatures <- function(modality_prefix, seurat_obj, n_topics = 50, top_genes = 25) {
  cat(paste0("\n--- Processing Modality: ", toupper(modality_prefix), " ---\n"))
  
  # Construct file paths
  path_rna_sig <- paste0("../data/p22_", modality_prefix, "_signatures_rna.csv")
  
  # Check if file exists
  if(!file.exists(path_rna_sig)) stop("File not found: ", path_rna_sig)
  
  # Load RNA signatures
  rna_signatures <- fread(path_rna_sig)
  
  # Extract top genes for each topic
  cat("Extracting top genes and calculating module scores...\n")
  for(i in 1:n_topics) {
    topic_df <- order_topic(rna_signatures, i)
    top_features <- list(head(topic_df, n = top_genes)$Feature)
    
    # Calculate Module Score. Note: Seurat auto-appends '1' to the name (e.g., Topic_11)
    seurat_obj <- AddModuleScore(
      nbin = 20,
      object = seurat_obj,
      features = top_features,
      ctrl = 100,
      name = paste0(toupper(modality_prefix), "_Topic_", i, "_")
    )
  }
  return(seurat_obj)
}

get_prefix <- function(x) {
  stringr::str_extract(x, "^[A-Za-z]+")
}

get_index <- function(x) {
  suppressWarnings(as.integer(stringr::str_extract(x, "[0-9]+$")))
}

plot_signature_density_zoom <- function(
    seurat_obj,
    signature,
    target_clusters,
    reduction = "umap",
    cluster_col = "ClusterName",
    x = 2.5,
    point_size = 0.20,
    target_point_size = NULL,
    point_alpha = 0.95,
    other_alpha = 0.35,
    label_size = 3.8,
    repel_force = 1,
    show_box_on_full_umap = TRUE,
    full_plot_color_all = FALSE,     # full plot: color all clusters by cell type
    color_other_cells = FALSE,       # zoom plot: color non-target cells by cell type
    grey_all_cells = FALSE,          # if TRUE, all cells are grey; no cell-type color legend
    grey_color = "gray40",
    show_celltype_legend = TRUE,     # legend for target cell types
    show_density_legend = TRUE,      # legend/colorbar for signature density
    show_celltype_labels = TRUE,     # show target cell-type labels on UMAP
    palette = c(
      "#cdba96", "#000080", "#ffff00", "#9ecae1", "#189e7f", "#fff8dc",
      "#642915", "#0000ff", "#bebebe", "#4292c6", "#cd0000", "#58508d",
      "#006400", "#10d894", "#666666", "#800080", "#ff00ff", "#ffc0cb",
      "#bebebe", "#ffd700", "#bebebe", "#00ffff", "#cd853f", "#dda0dd",
      "#000000", "#551a8b", "#00cd00", "#f4a460", "#00ff7f", "#ffbbff"
    ),
    cluster_colors = NULL,           # optional named vector of cluster colors
    bandwidth = c(0.5, 0.5),
    quantile_cutoff = 0.995,
    use_positive_only = FALSE,       # TRUE => dens_df = df %>% filter(Score > 0)
    grid_n = 300,
    density_alpha_max = 0.95,
    density_alpha_min_display = 0.02,
    #density_colors = c("white", "#F6E58D", "#F39C12", "red", "#cd0000","brown4"),
    density_colors = inferno(10),
    density_on_top = TRUE,
    seed_use = 10
) {
  set.seed(seed_use)
  
  if (is.null(target_point_size)) {
    target_point_size <- point_size
  }
  
  stopifnot(cluster_col %in% colnames(seurat_obj@meta.data))
  stopifnot(reduction %in% names(seurat_obj@reductions))
  
  emb <- Seurat::Embeddings(seurat_obj, reduction = reduction)[, 1:2, drop = FALSE]
  vals <- Seurat::FetchData(seurat_obj, vars = signature)[, 1]
  
  df <- tibble::tibble(
    cell        = rownames(emb),
    UMAP_1      = emb[, 1],
    UMAP_2      = emb[, 2],
    ClusterName = as.character(seurat_obj@meta.data[rownames(emb), cluster_col]),
    Score       = as.numeric(vals)
  ) %>%
    dplyr::filter(is.finite(Score))
  
  missing_targets <- setdiff(target_clusters, unique(df$ClusterName))
  if (length(missing_targets) > 0) {
    stop(
      paste0(
        "These target clusters were not found in ", cluster_col, ": ",
        paste(missing_targets, collapse = ", ")
      )
    )
  }
  
  ## ------------------------------------------------------------------
  ## stable cluster color mapping
  ## ------------------------------------------------------------------
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
  
  ## ------------------------------------------------------------------
  ## density source cells
  ## ------------------------------------------------------------------
  if (use_positive_only) {
    thr <- 0
    dens_df <- df %>% dplyr::filter(Score > 0)
    density_subtitle <- NULL
  } else {
    thr <- stats::quantile(df$Score, quantile_cutoff, na.rm = TRUE)
    dens_df <- df %>% dplyr::filter(Score >= thr)
    density_subtitle <- NULL
  }
  
  if (nrow(dens_df) < 10) {
    stop("Too few cells selected for density. Lower quantile_cutoff or use use_positive_only = TRUE.")
  }
  
  ## ------------------------------------------------------------------
  ## zoom center and box
  ## ------------------------------------------------------------------
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
    dplyr::mutate(
      is_target = ClusterName %in% target_clusters,
      PlotColor = unname(cluster_cols[ClusterName])
    )
  
  if (nrow(zoom_df) == 0) {
    stop("No cells found inside the zoom box. Try increasing x.")
  }
  
  dens_zoom_df <- dens_df %>%
    dplyr::filter(
      UMAP_1 >= x_min, UMAP_1 <= x_max,
      UMAP_2 >= y_min, UMAP_2 <= y_max
    )
  
  target_zoom_df <- zoom_df %>% dplyr::filter(is_target)
  other_zoom_df  <- zoom_df %>% dplyr::filter(!is_target)
  
  target_clusters_present <- target_clusters[target_clusters %in% unique(target_zoom_df$ClusterName)]
  
  label_df_zoom <- target_zoom_df %>%
    dplyr::group_by(ClusterName) %>%
    dplyr::summarise(
      UMAP_1 = median(UMAP_1),
      UMAP_2 = median(UMAP_2),
      .groups = "drop"
    )
  
  label_df_zoom$LegendCluster <- factor(label_df_zoom$ClusterName, levels = target_clusters_present)
  target_zoom_df$LegendCluster <- factor(target_zoom_df$ClusterName, levels = target_clusters_present)
  
  ## ------------------------------------------------------------------
  ## helper to build KDE grid
  ## ------------------------------------------------------------------
  build_density_grid <- function(density_cells, xlim, ylim, bandwidth, grid_n, alpha_floor) {
    if (nrow(density_cells) < 10) {
      return(data.frame(
        UMAP_1 = numeric(0),
        UMAP_2 = numeric(0),
        density = numeric(0),
        density_scaled = numeric(0),
        density_alpha = numeric(0)
      ))
    }
    
    kd <- MASS::kde2d(
      x = density_cells$UMAP_1,
      y = density_cells$UMAP_2,
      n = grid_n,
      h = bandwidth,
      lims = c(xlim[1], xlim[2], ylim[1], ylim[2])
    )
    
    kd_df <- expand.grid(
      UMAP_1 = kd$x,
      UMAP_2 = kd$y
    )
    kd_df$density <- as.vector(kd$z)
    
    max_d <- max(kd_df$density, na.rm = TRUE)
    kd_df$density_scaled <- if (is.finite(max_d) && max_d > 0) kd_df$density / max_d else 0
    kd_df$density_alpha  <- ifelse(kd_df$density_scaled < alpha_floor, 0, kd_df$density_scaled)
    
    kd_df
  }
  
  full_density_df <- build_density_grid(
    density_cells = dens_df,
    xlim = range(df$UMAP_1, na.rm = TRUE),
    ylim = range(df$UMAP_2, na.rm = TRUE),
    bandwidth = bandwidth,
    grid_n = grid_n,
    alpha_floor = density_alpha_min_display
  )
  
  zoom_density_df <- build_density_grid(
    density_cells = dens_zoom_df,
    xlim = c(x_min, x_max),
    ylim = c(y_min, y_max),
    bandwidth = bandwidth,
    grid_n = grid_n,
    alpha_floor = density_alpha_min_display
  )
  
  zoom_density_layer <- geom_raster(
    data = zoom_density_df,
    aes(x = UMAP_1, y = UMAP_2, fill = density, alpha = density_alpha),
    inherit.aes = FALSE,
    interpolate = TRUE,
    show.legend = show_density_legend
  )
  
  full_density_layer <- geom_raster(
    data = full_density_df,
    aes(x = UMAP_1, y = UMAP_2, fill = density, alpha = density_alpha),
    inherit.aes = FALSE,
    interpolate = TRUE,
    show.legend = show_density_legend
  )
  
  ## ------------------------------------------------------------------
  ## helper: add point layers
  ## ------------------------------------------------------------------
  add_zoom_points <- function(p) {
    if (grey_all_cells) {
      p <- p +
        geom_point(
          data = zoom_df,
          aes(x = UMAP_1, y = UMAP_2),
          color = grey_color,
          size = point_size,
          alpha = other_alpha,
          stroke = 0,
          show.legend = FALSE
        )
    } else {
      if (nrow(other_zoom_df) > 0) {
        if (color_other_cells) {
          p <- p +
            geom_point(
              data = other_zoom_df,
              aes(x = UMAP_1, y = UMAP_2),
              color = other_zoom_df$PlotColor,
              size = point_size,
              alpha = other_alpha,
              stroke = 0,
              show.legend = FALSE
            )
        } else {
          p <- p +
            geom_point(
              data = other_zoom_df,
              aes(x = UMAP_1, y = UMAP_2),
              color = grey_color,
              size = point_size,
              alpha = other_alpha,
              stroke = 0,
              show.legend = FALSE
            )
        }
      }
      
      if (nrow(target_zoom_df) > 0) {
        p <- p +
          geom_point(
            data = target_zoom_df,
            aes(x = UMAP_1, y = UMAP_2, color = LegendCluster),
            size = target_point_size,
            alpha = point_alpha,
            stroke = 0,
            show.legend = show_celltype_legend
          )
      }
    }
    p
  }
  
  add_full_points <- function(p) {
    df_full <- df %>%
      dplyr::mutate(PlotColor = unname(cluster_cols[ClusterName]))
    full_target_df <- df_full %>%
      dplyr::filter(ClusterName %in% target_clusters) %>%
      dplyr::mutate(LegendCluster = factor(ClusterName, levels = target_clusters_present))
    
    if (grey_all_cells) {
      p <- p +
        geom_point(
          data = df_full,
          aes(x = UMAP_1, y = UMAP_2),
          color = grey_color,
          size = point_size,
          alpha = other_alpha,
          stroke = 0,
          show.legend = FALSE
        )
    } else {
      if (full_plot_color_all) {
        p <- p +
          geom_point(
            data = df_full,
            aes(x = UMAP_1, y = UMAP_2),
            color = df_full$PlotColor,
            size = point_size,
            alpha = other_alpha,
            stroke = 0,
            show.legend = FALSE
          )
      } else {
        p <- p +
          geom_point(
            data = df_full,
            aes(x = UMAP_1, y = UMAP_2),
            color = grey_color,
            size = point_size,
            alpha = other_alpha,
            stroke = 0,
            show.legend = FALSE
          )
      }
      
      if (nrow(full_target_df) > 0) {
        p <- p +
          geom_point(
            data = full_target_df,
            aes(x = UMAP_1, y = UMAP_2, color = LegendCluster),
            size = target_point_size,
            alpha = point_alpha,
            stroke = 0,
            show.legend = show_celltype_legend
          )
      }
    }
    p
  }
  
  ## ------------------------------------------------------------------
  ## zoom plot
  ## ------------------------------------------------------------------
  p_zoom <- ggplot()
  
  if (density_on_top) {
    p_zoom <- add_zoom_points(p_zoom)
    p_zoom <- p_zoom + zoom_density_layer
  } else {
    p_zoom <- p_zoom + zoom_density_layer
    p_zoom <- add_zoom_points(p_zoom)
  }
  
  if (show_celltype_labels && nrow(label_df_zoom) > 0) {
    p_zoom <- p_zoom +
      ggrepel::geom_text_repel(
        data = label_df_zoom,
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
      )
  }
  
  p_zoom <- p_zoom +
    scale_fill_gradientn(
      colours = density_colors,
      name = signature,
      guide = if (show_density_legend) {
        guide_colorbar(barheight = grid::unit(35, "mm"))
      } else {
        "none"
      }
    ) +
    scale_alpha(
      range = c(0, density_alpha_max),
      guide = "none"
    ) +
    coord_cartesian(
      xlim = c(x_min, x_max),
      ylim = c(y_min, y_max),
      expand = FALSE
    ) +
    labs(
      x = "UMAP_1",
      y = "UMAP_2",
      title = paste0(signature, " | Zoom: ", paste(target_clusters, collapse = ", ")),
      subtitle = density_subtitle
    ) +
    theme_classic(base_size = 11) +
    theme(
      axis.title = element_text(size = 12, color = "black"),
      axis.text  = element_text(size = 10, color = "black"),
      axis.line  = element_line(linewidth = 0.5, color = "black"),
      plot.title = element_text(size = 12, face = "bold"),
      plot.margin = margin(4, 6, 4, 4, "mm")
    )
  
  if (!grey_all_cells && show_celltype_legend && length(target_clusters_present) > 0) {
    p_zoom <- p_zoom +
      scale_color_manual(
        values = cluster_cols[target_clusters_present],
        breaks = target_clusters_present,
        drop = FALSE,
        guide = guide_legend(
          title = NULL,
          override.aes = list(shape = 16, size = 3.5, alpha = 1, stroke = 0)
        )
      )
  }
  
  ## ------------------------------------------------------------------
  ## full plot
  ## ------------------------------------------------------------------
  p_full <- NULL
  if (show_box_on_full_umap) {
    p_full <- ggplot()
    
    if (density_on_top) {
      p_full <- add_full_points(p_full)
      p_full <- p_full + full_density_layer
    } else {
      p_full <- p_full + full_density_layer
      p_full <- add_full_points(p_full)
    }
    
    p_full <- p_full +
      #     geom_rect(
      #        inherit.aes = FALSE,
      #        aes(xmin = x_min, xmax = x_max, ymin = y_min, ymax = y_max),
      #        fill = NA,
      #        color = "black",
      #        linewidth = 0.5
      #      ) +
      scale_fill_gradientn(
        colours = density_colors,
        name = signature,
        guide = if (show_density_legend) {
          guide_colorbar(barheight = grid::unit(35, "mm"))
        } else {
          "none"
        }
      ) +
      scale_alpha(
        range = c(0, density_alpha_max),
        guide = "none"
      ) +
      labs(
        x = "UMAP_1",
        y = "UMAP_2",
        title = signature,
        subtitle = density_subtitle
      ) +
      theme_classic(base_size = 11) +
      theme(
        axis.title = element_text(size = 12, color = "black"),
        axis.text  = element_text(size = 10, color = "black"),
        axis.line  = element_line(linewidth = 0.5, color = "black"),
        plot.title = element_text(size = 12, face = "bold"),
        plot.margin = margin(4, 6, 4, 4, "mm")
      )
    
    if (!grey_all_cells && show_celltype_legend && length(target_clusters_present) > 0) {
      p_full <- p_full +
        scale_color_manual(
          values = cluster_cols[target_clusters_present],
          breaks = target_clusters_present,
          drop = FALSE,
          guide = guide_legend(
            title = NULL,
            override.aes = list(shape = 16, size = 3.5, alpha = 1, stroke = 0)
          )
        )
    }
  }
  
  return(list(
    plot = p_zoom,
    full_plot = p_full,
    all_data = df,
    zoom_data = zoom_df,
    density_data_full = full_density_df,
    density_data_zoom = zoom_density_df,
    density_source_data = dens_df,
    label_data_zoom = label_df_zoom,
    threshold = thr,
    center = c(UMAP_1 = center_x, UMAP_2 = center_y),
    box = c(xmin = x_min, xmax = x_max, ymin = y_min, ymax = y_max),
    colors = cluster_cols[target_clusters_present]
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

#Select modality
modality<-"atac"
MouseAdoAtlas_atac <- add_modality_signatures(modality_prefix = modality, seurat_obj = MouseAdoAtlas, n_topics = 50, top_genes = 25)

quantile_topic <- list('10' = 0.99, '37' = 0.99, '4' = 0.99, '45' = 0.99,
                       '12' = 0.99, '19' = 0.99, '21' = 0.9995, '25' = 0.99,
                       '31' = 0.999, '33' = 0.999, '41' = 0.9995, '43' = 0.999,
                       '49' = 0.99)

list_of_figures <- list()
for (Topic in c(10, 12, 37, 45, 4, 41, 25, 33, 19, 21, 31, 43, 49)){
  res_sig <- plot_signature_density_zoom(
    seurat_obj = MouseAdoAtlas_atac,
    signature = paste0(toupper(modality), "_Topic_", Topic, "_1"),
    target_clusters = "TEGLU4",
    x = 5,
    point_size = 0.4,
    label_size = 4,
    cluster_colors = my_cluster_colors,
    show_box_on_full_umap = T,
    full_plot_color_all = FALSE,
    color_other_cells = FALSE,
    grey_all_cells = T,
    show_celltype_labels = F,
    show_celltype_legend = TRUE,
    show_density_legend = TRUE,
    bandwidth = c(4.5, 4.5),
    quantile_cutoff = quantile_topic[[as.character(Topic)]],
    use_positive_only = FALSE,
    density_on_top = TRUE,
    grid_n = 500,
    density_alpha_max = 1,
    density_alpha_min_display = 0.02
  )
  list_of_figures[[as.character(Topic)]] <- res_sig$full_plot
}

png("../figures/Figure_E2/Figure_E2a_3.png", height = 3500, width = 5000, res = 200)
wrap_plots(list_of_figures, ncol = 4)
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
#   [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C               LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8    
# [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8    LC_PAPER=en_US.UTF-8       LC_NAME=C                 
# [9] LC_ADDRESS=C               LC_TELEPHONE=C             LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       
# 
# time zone: Europe/Warsaw
# tzcode source: system (glibc)
# 
# attached base packages:
#   [1] grid      stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
#   [1] patchwork_1.3.2     data.table_1.18.2.1 viridis_0.6.5       viridisLite_0.4.3   MASS_7.3-65         tibble_3.3.1       
# [7] ggrepel_0.9.8       ggplot2_4.0.2       stringr_1.6.0       dplyr_1.2.1         Seurat_5.4.0        SeuratObject_5.3.0 
# [13] sp_2.2-1           
# 
# loaded via a namespace (and not attached):
#   [1] deldir_2.0-4           pbapply_1.7-4          gridExtra_2.3          rlang_1.2.0            magrittr_2.0.5        
# [6] RcppAnnoy_0.0.23       otel_0.2.0             matrixStats_1.5.0      ggridges_0.5.7         compiler_4.5.2        
# [11] spatstat.geom_3.7-3    png_0.1-9              vctrs_0.7.2            reshape2_1.4.5         pkgconfig_2.0.3       
# [16] fastmap_1.2.0          labeling_0.4.3         promises_1.5.0         purrr_1.2.1            jsonlite_2.0.0        
# [21] goftest_1.2-3          later_1.4.8            spatstat.utils_3.2-2   irlba_2.3.7            parallel_4.5.2        
# [26] cluster_2.1.8.1        R6_2.6.1               ica_1.0-3              stringi_1.8.7          RColorBrewer_1.1-3    
# [31] spatstat.data_3.1-9    reticulate_1.46.0      parallelly_1.46.1      spatstat.univar_3.1-7  lmtest_0.9-40         
# [36] scattermore_1.2        Rcpp_1.1.1             tensor_1.5.1           future.apply_1.20.2    zoo_1.8-15            
# [41] sctransform_0.4.3      httpuv_1.6.17          Matrix_1.7-4           splines_4.5.2          igraph_2.2.3          
# [46] tidyselect_1.2.1       rstudioapi_0.18.0      abind_1.4-8            spatstat.random_3.4-5  codetools_0.2-20      
# [51] miniUI_0.1.2           spatstat.explore_3.8-0 listenv_0.10.1         lattice_0.22-7         plyr_1.8.9            
# [56] withr_3.0.2            shiny_1.13.0           S7_0.2.1               ROCR_1.0-12            Rtsne_0.17            
# [61] future_1.70.0          fastDummies_1.7.5      survival_3.8-3         polyclip_1.10-7        fitdistrplus_1.2-6    
# [66] pillar_1.11.1          KernSmooth_2.23-26     plotly_4.12.0          generics_0.1.4         RcppHNSW_0.6.0        
# [71] scales_1.4.0           globals_0.19.1         xtable_1.8-8           glue_1.8.0             lazyeval_0.2.3        
# [76] tools_4.5.2            RSpectra_0.16-2        RANN_2.6.2             dotCall64_1.2          cowplot_1.2.0         
# [81] tidyr_1.3.2            nlme_3.1-168           cli_3.6.6              spatstat.sparse_3.1-0  spam_2.11-3           
# [86] uwot_0.2.4             gtable_0.3.6           digest_0.6.39          progressr_0.19.0       htmlwidgets_1.6.4     
# [91] farver_2.1.2           htmltools_0.5.9        lifecycle_1.0.5        httr_1.4.8             mime_0.13   