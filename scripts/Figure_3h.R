library(ggplot2)
library(ggrepel)
library(data.table)
library(dplyr)
library(this.path)

setwd(this.dir())

dir.create("../figures/Figure_3", recursive = T)

data <- fread(file="../data/PBMC_Cyto_PeakGeneAssoc.csv")

top_HHI <- data %>%
  group_by(Avg_log2FC_Type) %>%
  slice_max(order_by = HHI, n = 1) %>%
  slice_max(order_by = GA, n = 1) %>%
  ungroup()

size <- 12
RTE_GA <- ggplot(data, aes(x = GA, y = HHI)) +
  geom_point(size=0.5) +
  geom_text_repel(data = top_HHI, aes(label = label), box.padding = 0.5, size=2) +
  facet_wrap(~Avg_log2FC_Type, scales = "free", ncol = 4) +
  theme_minimal() +
  theme(panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour="black", size = 0.5),
        text=element_text(size=size),
        axis.ticks = element_line(size=0.5,colour='black'),
        axis.text.x = element_text(size=size, colour='black', hjust = 1),
        axis.text.y = element_text(size=size, colour='black'),
        plot.title=element_text(size=size, colour='black'),
        axis.title.x=element_text(size=size, colour='black'),
        axis.title.y=element_text(size=size, colour='black'),
        legend.position="none")+
  scale_x_log10() + 
  labs(x = "Global association score, (log10)", y = "HHI") +
  scale_y_continuous(limits = c(0, 1)) 

ggsave("../figures/Figure_3/Figure_3h.pdf", width = 12, height = 4, RTE_GA)


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
#   [1] dplyr_1.2.1         data.table_1.18.2.1 ggrepel_0.9.8       ggplot2_4.0.2      
# 
# loaded via a namespace (and not attached):
#   [1] labeling_0.4.3     RColorBrewer_1.1-3 R6_2.6.1           tidyselect_1.2.1   farver_2.1.2       magrittr_2.0.5     gtable_0.3.6      
# [8] glue_1.8.0         tibble_3.3.1       pkgconfig_2.0.3    generics_0.1.4     lifecycle_1.0.5    cli_3.6.6          S7_0.2.1          
# [15] scales_1.4.0       grid_4.5.2         vctrs_0.7.2        textshaping_1.0.5  withr_3.0.2        systemfonts_1.3.2  compiler_4.5.2    
# [22] rstudioapi_0.18.0  tools_4.5.2        ragg_1.5.0         pillar_1.11.1      Rcpp_1.1.1         rlang_1.2.0    