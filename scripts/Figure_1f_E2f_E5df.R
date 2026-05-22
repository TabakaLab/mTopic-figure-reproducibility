library(Seurat)
library(ggplot2)
library(data.table)
library(dplyr)
library(ggrepel)
library(reshape2)
library(this.path)

setwd(this.dir())

dir.create("../figures/Figure_1", recursive = T)
dir.create("../figures/Figure_E2", recursive = T)
dir.create("../figures/Figure_E5", recursive = T)

source("allen_sections_grid.R")
clear_allen_cache()

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

scale_values <- function(x) {
  (x - min(x, na.rm = TRUE)) / (max(x, na.rm = TRUE) - min(x, na.rm = TRUE))
}

plot_modality_signatures <- function(modality_prefix, seurat_obj, n_topics = 50, top_genes = 25) {
  cat(paste0("\n--- Processing Modality: ", toupper(modality_prefix), " ---\n"))

  path_rna_sig <- paste0("../data/p22_", modality_prefix, "_signatures_rna.csv")
  if(!file.exists(path_rna_sig)) stop("File not found: ", path_rna_sig)

  rna_signatures <- fread(path_rna_sig)
  cat("Extracting top genes and calculating module scores...\n")
  for(i in 1:n_topics) {
    topic_df <- order_topic(rna_signatures, i)
    top_features <- list(head(topic_df, n = top_genes)$Feature)

    seurat_obj <- AddModuleScore(
      nbin = 20,
      object = seurat_obj,
      features = top_features,
      ctrl = 100,
      name = paste0(toupper(modality_prefix), "_Topic_", i, "_")
    )
  }

  topic_cols <- paste0(toupper(modality_prefix), "_Topic_", 1:n_topics, "_1")
  tmp <- seurat_obj@meta.data[, c("ClusterName", topic_cols)] 
  tmp <- reshape2::melt(tmp, id.vars = "ClusterName", variable.name = "Topic", value.name = "Score")

  cat("Formatting data for plotting...\n")
  n_clusters <- length(unique(tmp$ClusterName))

  plot_df <- tmp %>%
    group_by(Topic, ClusterName) %>% 
    summarise(Score = mean(Score, na.rm = TRUE), .groups = 'drop') %>%
    group_by(Topic) %>% 
    arrange(desc(Score), .by_group = TRUE) %>%
    mutate(
      Rank = row_number(),
      Scaled = scale_values(Score)
    ) %>%
    ungroup()

  cat("Generating output plots...\n")
  plot_size <- 12
  if(modality_prefix == "h3k27me3") {
    out_pdf <- paste0("../figures/Figure_E5/Figure_E5df_Ranksorted_Signatures_", toupper(modality_prefix), "_Top", top_genes, ".pdf")
  } else {
    out_pdf <- paste0("../figures/Figure_1/Figure_1f_E2f_Ranksorted_Signatures_", toupper(modality_prefix), "_Top", top_genes, ".pdf")
  }
  #out_csv <- paste0("../figures/Figure_1/Figure_1f_Ranksorted_Signatures_", toupper(modality_prefix), "_Top", top_genes, ".csv")

  top_labels <- plot_df %>% filter(Rank <= 5)
  p_final <- ggplot(plot_df, aes(x = Rank, y = Scaled)) +
    geom_point(size = 0.2, color = 'grey30') +
    facet_wrap(~Topic, ncol = 10) +
    theme_classic(base_size = plot_size) +
    theme(
      axis.line = element_line(colour = "black", linewidth = 0.5),
      axis.text.x = element_blank(), 
      axis.ticks.x = element_blank(),
      strip.background = element_blank(),
      strip.text = element_text(face = "bold"),
      legend.position = "none"
    ) +
    labs(
      title = paste0(toupper(modality_prefix), " - Rank-sorted Signatures (Top ", top_genes, " Genes)"),
      x = "Cluster Rank",
      y = "Scaled Module Score"
    ) + 
    geom_text_repel(
      data = top_labels,
      aes(label = ClusterName),
      size = 3.5, 
      segment.size = 0.2, 
      colour = "black",
      max.overlaps = 50
    )

  ggsave(filename = out_pdf, plot = p_final, width = 35, height = 25, limitsize = FALSE)
  #fwrite(plot_df, file = out_csv, row.names = FALSE)
  return(seurat_obj)
}

# Run ATAC
MouseAdoAtlas_atac <- plot_modality_signatures(modality_prefix = "atac", seurat_obj = MouseAdoAtlas)

# # Run H3K4me3 (Uncomment to run)
# MouseAdoAtlas_h3k4 <- plot_modality_signatures(modality_prefix = "h3k4me3", seurat_obj = MouseAdoAtlas)
# 
# # Run H3K27me3 (Uncomment to run)
MouseAdoAtlas_h3k27me3 <- plot_modality_signatures(modality_prefix = "h3k27me3", seurat_obj = MouseAdoAtlas)
# 
# # Run H3K27ac (Uncomment to run)
# MouseAdoAtlas_h3k27ac <- plot_modality_signatures(modality_prefix = "h3k27ac", seurat_obj = MouseAdoAtlas)


get_signature <- function(modality_prefix, topic = 1, top_genes = 25) {
  cat(paste0("\n--- Processing Modality: ", toupper(modality_prefix),", Topic=",topic, " ---\n"))
  path_rna_sig <- paste0("../data/p22_", modality_prefix, "_signatures_rna.csv")

  if(!file.exists(path_rna_sig)) stop("File not found: ", path_rna_sig)
  rna_signatures <- fread(path_rna_sig)
  topic_df <- order_topic(rna_signatures, topic)
  top_features <- head(topic_df, n = top_genes)$Feature

  return(top_features)
}
order_topic <- function(signatures, topic=1) {
  topic_col <- paste0("topic_", topic)
  ind <- order(signatures[[topic_col]], decreasing = TRUE)
  df <- data.frame(
    Feature = signatures$V1[ind], 
    Score = signatures[[topic_col]][ind]
  )
  return(df)
}

Topics <- c(4,10,12,19,25,31,33,37,41,43,45,49)
GS <- c(10,10,5,10,50,25,10,10,10,25,25,50)

for (i in 1:length(Topics)){
  Genes <- get_signature(
    "atac", topic = Topics[i], top_genes = GS[i])
  plot_allen_sections_grid(
    list_genes = Genes,
    sections = 24,
    atlas_depth = 20,
    smooth_factor = 12,
    expr_interp = "nearest",
    boundary_lwd = 0.2,
    main_title = paste("Topic", Topics[i]),
    out_file = paste0("../figures/Figure_1/Figure_1f_topic_", Topics[i], ".png"),
    show_plot = FALSE
  )
}


Genes <- get_signature("atac", topic = 21, top_genes = 50)
plot_allen_sections_grid(
  list_genes = Genes,
  sections = 24,
  atlas_depth = 20,
  smooth_factor = 12,
  expr_interp = "nearest",
  boundary_lwd = 0.2,
  main_title = "Topic 21",
  show_plot = FALSE,
  out_file = paste0("../figures/Figure_E2/Figure_E2f_4.png")
)
