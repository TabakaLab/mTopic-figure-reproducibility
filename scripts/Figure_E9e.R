require(igraph)
require(tidygraph)
require(ggraph)
require(rgexf)
require(data.table)
library(dplyr)
require(viridis)
require(colorRamps)
library(this.path)

setwd(this.dir())

dir.create("../figures/Figure_E9", recursive = T)

graph_tbl <- readRDS(file="../data/PBMCs_GRN_RTE.RDS")
node_colors <- c("darkred","green","blue")
graph_tbl <- graph_tbl %>%
  activate(edges) %>%
  mutate(edge_color = node_colors[from])

graph_tbl <- graph_tbl %>% filter(!(from == to))

cols <- as_tibble(graph_tbl, what = "edges") %>% 
  pull(edge_color)

set.seed(3)
gg <- ggraph(graph_tbl,'fr') +
  geom_edge_bundle_force(edge_colour = rep(cols, each=20), n_cycle = 10, force=4, 
                         edge_width = 0.4, n = 20, alpha=0.4)+
  scale_edge_colour_identity() +
  geom_node_point(size = 1, aes(color = color)) + 
  scale_color_manual(values = c("blue" = "blue", "red" = "red"))+
  geom_node_text(aes(label = name), repel = T, size=4) +
  theme_graph() + 
  theme(legend.position="none")

png(file = paste0("../figures/Figure_E9/Figure_E9e.png"), width = 300, height=200, units="mm", res =400)
gg
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
#   [1] generics_0.1.4     tidyr_1.3.2        magrittr_2.0.5     grid_4.5.2         RColorBrewer_1.1-3 fastmap_1.2.0      jsonlite_2.0.0    
# [8] ggrepel_0.9.8      gridExtra_2.3      promises_1.5.0     purrr_1.2.1        scales_1.4.0       tweenr_2.0.3       XML_3.99-0.22     
# [15] cli_3.6.6          rlang_1.2.0        graphlayouts_1.2.2 polyclip_1.10-7    withr_3.0.2        cachem_1.1.0       otel_0.2.0        
# [22] tools_4.5.2        memoise_2.0.1      httpuv_1.6.17      vctrs_0.7.2        R6_2.6.1           lifecycle_1.0.5    MASS_7.3-65       
# [29] pkgconfig_2.0.3    pillar_1.11.1      later_1.4.8        gtable_0.3.6       glue_1.8.0         Rcpp_1.1.1         servr_0.32        
# [36] ggforce_0.5.0      xfun_0.57          tibble_3.3.1       tidyselect_1.2.1   rstudioapi_0.18.0  farver_2.1.2       labeling_0.4.3    
# [43] compiler_4.5.2     S7_0.2.1   
