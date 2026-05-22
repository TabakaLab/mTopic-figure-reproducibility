library(ggplot2)
library(ggpubr)
library(this.path)

setwd(this.dir())

dir.create("../figures/Figure_E9", recursive = T)

data <- read.csv("../data/PBMC_GO.csv", sep=";")
data_clean <- data[, c("Gene.Set.Name", "FDR.q.value")]

data_clean$log10_FDR <- -log10(data_clean$FDR.q.value)

data_clean$Gene.Set.Name <- gsub("GOBP_", "", data_clean$Gene.Set.Name)  # Remove "GOBP" prefix
data_clean$Gene.Set.Name <- gsub("_", " ", data_clean$Gene.Set.Name)

size <- 12
cyto <- ggplot(data_clean, aes(x = reorder(Gene.Set.Name, log10_FDR), 
                               y = log10_FDR)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  labs(x = "Gene Set Name", y = "-log10(FDR)") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
  theme(panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black",size = 0.5),
        text=element_text(size=size),
        axis.ticks = element_line(size = 0.5,colour='black'),
        axis.text.x = element_text(size=size,colour='black', hjust = 1),
        axis.text.y = element_text(size=size,colour='black'),
        plot.title=element_text(size=size,colour='black'),
        axis.title.x=element_text(size=size,colour='black'),
        axis.title.y=element_text(size=size,colour='black'),
        legend.position="none")

ggsave("../figures/Figure_E9/Figure_E9h.pdf", width = 9, height = 3, cyto)


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
#   [1] ggpubr_0.6.2  ggplot2_4.0.2
# 
# loaded via a namespace (and not attached):
#   [1] vctrs_0.7.2        cli_3.6.6          rlang_1.2.0        Formula_1.2-5      purrr_1.2.1        car_3.1-5          generics_0.1.4    
# [8] textshaping_1.0.5  S7_0.2.1           labeling_0.4.3     glue_1.8.0         backports_1.5.1    ragg_1.5.0         scales_1.4.0      
# [15] grid_4.5.2         abind_1.4-8        carData_3.0-6      tibble_3.3.1       rstatix_0.7.3      lifecycle_1.0.5    ggsignif_0.6.4    
# [22] compiler_4.5.2     dplyr_1.2.1        RColorBrewer_1.1-3 pkgconfig_2.0.3    tidyr_1.3.2        rstudioapi_0.18.0  systemfonts_1.3.2 
# [29] farver_2.1.2       R6_2.6.1           tidyselect_1.2.1   pillar_1.11.1      magrittr_2.0.5     tools_4.5.2        withr_3.0.2       
# [36] gtable_0.3.6       broom_1.0.12 
