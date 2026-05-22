library(networkD3)
library(htmlwidgets)
library(webshot2)
library(reshape2)
library(this.path)

setwd(this.dir())

dir.create("../figures/Figure_E5", recursive = T)

sim_matrix <- readRDS("../data/SimMatrix_ATAC_H3K27me3_rbo.RDS")
sim_matrix <- t(sim_matrix)

row_ids <- c(1, 2, 3, 19, 32, 33, 46, 48, 14, 45)
col_ids <- c(10, 12, 37, 45, 31)

sim_matrix <- sim_matrix[row_ids, col_ids, drop = FALSE]
nodes <- data.frame("name" = c(rownames(sim_matrix), colnames(sim_matrix))) # Node 6

rownames(sim_matrix) <- 0:(nrow(sim_matrix)-1)
colnames(sim_matrix) <- (nrow(sim_matrix)):(nrow(sim_matrix)+ncol(sim_matrix)-1)

keep_top <- function(row) {
  cutoff <- sort(row, decreasing = TRUE)[1]
  row[row < cutoff] <- 0
  return(row)
}

sim_matrix <- t(apply(sim_matrix, 1, keep_top))
tmp <- melt(sim_matrix,value.name = "Score")
links <- tmp
names(links) <- c("target", "source", "value")
links <- links[which(links$value!=0),]
sn <- sankeyNetwork(
  Links = links, Nodes = nodes,width = 400, height = 800, Source = "source",
  Target = "target", Value = "value", NodeID = "name", fontSize = 40,
  nodeWidth = 20, iterations = 0)

saveWidget(sn, "../figures/Figure_E5/Figure_E5bf.html", selfcontained = TRUE)
webshot2::webshot("../figures/Figure_E5/Figure_E5bf.html", 
                  "../figures/Figure_E5/Figure_E5bf.pdf")


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
#   [1] reshape2_1.4.5    webshot2_0.1.2    htmlwidgets_1.6.4 networkD3_0.4.1  
# 
# loaded via a namespace (and not attached):
#   [1] knitr_1.51        cli_3.6.6         data.tree_1.2.0   xfun_0.57         rlang_1.2.0       stringi_1.8.7     otel_0.2.0       
# [8] processx_3.8.7    promises_1.5.0    jsonlite_2.0.0    glue_1.8.0        plyr_1.8.9        htmltools_0.5.9   ps_1.9.2         
# [15] chromote_0.5.1    rmarkdown_2.31    evaluate_1.0.5    fastmap_1.2.0     yaml_2.3.12       lifecycle_1.0.5   stringr_1.6.0    
# [22] compiler_4.5.2    igraph_2.2.3      Rcpp_1.1.1        pkgconfig_2.0.3   websocket_1.4.4   rstudioapi_0.18.0 later_1.4.8      
# [29] digest_0.6.39     R6_2.6.1          magrittr_2.0.5    tools_4.5.2  