library(data.table)
library(ggplot2)
library(ggrepel)
library(this.path)

setwd(this.dir())

dir.create("../figures/Figure_E2", recursive = T)

plot_df <- fread("../data/TF_rms.csv")
plot <- ggplot(plot_df, aes(x = n_cells, y = value)) +
  geom_point(size = 1, alpha = 0.8) +
  geom_text_repel(
    data = plot_df[label_topic == TRUE],
    aes(label = topic_clean),
    size = 3,
    max.overlaps = Inf,
    box.padding = 0.3,
    point.padding = 0.2,
    show.legend = FALSE
  ) +
  scale_x_log10(breaks = c(20, 30, 50, 100, 300, 1000, 3000)) +
  labs(x = "Spots per topic", y = "RMS", title = "Labels: topics") +
  theme_bw(base_size = 12)
ggsave("../figures/Figure_E2/Figure_E2b.pdf", width = 3, height = 3, plot)

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
#   [1] ggrepel_0.9.8       here_1.0.2          ggplot2_4.0.2       data.table_1.18.2.1
# 
# loaded via a namespace (and not attached):
#   [1] vctrs_0.7.2        cli_3.6.6          rlang_1.2.0        generics_0.1.4     textshaping_1.0.5  S7_0.2.1           glue_1.8.0        
# [8] labeling_0.4.3     rprojroot_2.1.1    ragg_1.5.0         scales_1.4.0       grid_4.5.2         tibble_3.3.1       lifecycle_1.0.5   
# [15] compiler_4.5.2     dplyr_1.2.1        RColorBrewer_1.1-3 Rcpp_1.1.1         pkgconfig_2.0.3    rstudioapi_0.18.0  systemfonts_1.3.2 
# [22] farver_2.1.2       R6_2.6.1           tidyselect_1.2.1   pillar_1.11.1      magrittr_2.0.5     tools_4.5.2        withr_3.0.2       
# [29] gtable_0.3.6