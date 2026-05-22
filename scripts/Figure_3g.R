library(igraph)
library(tidygraph)
library(ggraph)
library(rgexf)
library(data.table)
library(dplyr)
library(viridis)
library(colorRamps)
library(this.path)

setwd(this.dir())

dir.create("../figures/Figure_3", recursive = T)

set.seed(20)
graph_tbl <- readRDS(file="../data/PBMCs_GRNgraph.RDS")
cols <- as_tibble(graph_tbl, what = "edges") %>% pull(edge_color)

selected_nodes <- c(
  "STAT1::STAT2", "IRF8", "IRF4", "ETV3", "TCF7", "ELK3", "TCF7L2","RUNX2",
  "CTCF","STAT3", "POU2F1", "TBX21","EOMES","RORA", "NFKB1","NFKB2","RELB", 
  "FOS", "JUN","BATF","BATF3","CEBPD","JUNB","ATF4","ATF7")

bold_labels <- selected_nodes
gg <- ggraph(graph_tbl, layout = "fr") +
  geom_edge_bundle_force(edge_colour=rep(cols,each=100), 
                         n_cycle = 5,force=3, edge_width = 0.1, alpha=0.1) +
  scale_edge_colour_identity() +
  geom_node_point(aes(color = color),size = 1,alpha=0.2) +
  scale_colour_identity() +
  geom_node_text(aes(label = name), repel = F,size=3,check_overlap = F) +
  theme_graph() 


df <- data.frame(Name=gg$data$name, X=gg$data$x, Y=gg$data$y)
selected_df <- subset(df, Name %in% selected_nodes)

min_distance_to_selected <- function(x, y, selected) {
  distances <- sqrt((selected$X - x)^2 + (selected$Y - y)^2)
  min(distances)
}

# Calculate the minimum distance for each node to any selected node
df$min_dist <- apply(df[, c("X", "Y")], 1, function(coord) {
  min_distance_to_selected(coord[1], coord[2], selected_df)
})

plot(density(df$min_dist, a=0.1))
# Define a threshold for "far away" (adjust this based on your data's scale)
threshold <- 2
far_away_nodes <- subset(df, min_dist > threshold)

selected_nodes<-c(selected_nodes,far_away_nodes$Name)

graph_tbl <- graph_tbl %>%
  activate(nodes) %>%
  mutate(label = ifelse(name %in% selected_nodes, name, NA))

graph_tbl <- graph_tbl %>% 
  activate(nodes) %>% 
  mutate(fontface   = ifelse(label %in% bold_labels, "bold", "plain"),
         label_size = ifelse(label %in% bold_labels, 5, 4))

set.seed(20)
png(file = paste0("../figures/Figure_3/Figure_3g.png"), width = 300, height=200, units="mm", res=400)
ggraph(graph_tbl, layout = "fr") +
  geom_edge_bundle_force(edge_colour=rep(cols,each=50),n = 50,n_cycle = 10,
                         force=4.5,edge_width = 0.2,alpha=0.07) +
  scale_edge_colour_identity() +
  geom_node_point(aes(color = color),size = 1,alpha=0.7) +
  scale_colour_identity() +
  geom_node_text(aes(label = label,fontface = fontface,size = label_size), 
                 repel = T,check_overlap = T) +
  scale_size_identity() +
  theme_graph() + theme(legend.position="none")
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
#   [1] colorRamps_2.3.4    viridis_0.6.5       viridisLite_0.4.3   dplyr_1.2.1         data.table_1.18.2.1 rgexf_0.16.3        ggraph_2.2.2       
# [8] ggplot2_4.0.2       tidygraph_1.3.1     igraph_2.2.3       
# 
# loaded via a namespace (and not attached):
#   [1] deldir_2.0-4           pbapply_1.7-4          gridExtra_2.3          rlang_1.2.0            magrittr_2.0.5         RcppAnnoy_0.0.23      
# [7] otel_0.2.0             matrixStats_1.5.0      ggridges_0.5.7         compiler_4.5.2         spatstat.geom_3.7-3    systemfonts_1.3.2     
# [13] png_0.1-9              vctrs_0.7.2            reshape2_1.4.5         stringr_1.6.0          pkgconfig_2.0.3        fastmap_1.2.0         
# [19] labeling_0.4.3         promises_1.5.0         ragg_1.5.0             purrr_1.2.1            xfun_0.57              cachem_1.1.0          
# [25] jsonlite_2.0.0         goftest_1.2-3          later_1.4.8            tweenr_2.0.3           spatstat.utils_3.2-2   irlba_2.3.7           
# [31] parallel_4.5.2         cluster_2.1.8.1        R6_2.6.1               ica_1.0-3              stringi_1.8.7          RColorBrewer_1.1-3    
# [37] spatstat.data_3.1-9    reticulate_1.46.0      parallelly_1.46.1      spatstat.univar_3.1-7  lmtest_0.9-40          scattermore_1.2       
# [43] Rcpp_1.1.1             tensor_1.5.1           future.apply_1.20.2    zoo_1.8-15             sctransform_0.4.3      httpuv_1.6.17         
# [49] Matrix_1.7-4           splines_4.5.2          tidyselect_1.2.1       rstudioapi_0.18.0      abind_1.4-8            spatstat.random_3.4-5 
# [55] codetools_0.2-20       miniUI_0.1.2           spatstat.explore_3.8-0 listenv_0.10.1         lattice_0.22-7         tibble_3.3.1          
# [61] plyr_1.8.9             withr_3.0.2            shiny_1.13.0           S7_0.2.1               ROCR_1.0-12            Rtsne_0.17            
# [67] future_1.70.0          fastDummies_1.7.5      survival_3.8-3         polyclip_1.10-7        fitdistrplus_1.2-6     pillar_1.11.1         
# [73] Seurat_5.4.0           KernSmooth_2.23-26     plotly_4.12.0          generics_0.1.4         RcppHNSW_0.6.0         sp_2.2-1              
# [79] servr_0.32             scales_1.4.0           globals_0.19.1         xtable_1.8-8           glue_1.8.0             lazyeval_0.2.3        
# [85] tools_4.5.2            RSpectra_0.16-2        RANN_2.6.2             XML_3.99-0.22          graphlayouts_1.2.2     dotCall64_1.2         
# [91] cowplot_1.2.0          grid_4.5.2             tidyr_1.3.2            nlme_3.1-168           patchwork_1.3.2        ggforce_0.5.0         
# [97] cli_3.6.6              spatstat.sparse_3.1-0  textshaping_1.0.5      spam_2.11-3            uwot_0.2.4             gtable_0.3.6          
# [103] digest_0.6.39          progressr_0.19.0       ggrepel_0.9.8          htmlwidgets_1.6.4      SeuratObject_5.3.0     farver_2.1.2          
# [109] memoise_2.0.1          htmltools_0.5.9        lifecycle_1.0.5        httr_1.4.8             mime_0.13              MASS_7.3-65   