library(data.table)
library(ggplot2)
library(dplyr)
library(this.path)

setwd(this.dir())

dir.create("../figures/Figure_E2", recursive = T)

df_colors <- data.frame(
  Topic = c("Topic_4","Topic_10","Topic_19","Topic_25","Topic_31",
            "Topic_33","Topic_37","Topic_41","Topic_43","Topic_45","Topic_49"),
  Color = c("#ffeef6","#00cd00","#000000","#f4a460","#003f5c",
            "#fff8dc","#c9f4c9","#ee3a8c","#bcb8f4","#006400","#8a508f")
)
topic_levels <- df_colors$Topic
topic_cols   <- setNames(df_colors$Color, df_colors$Topic)
df <- fread("../data/P22_atac_TF_scores.csv")
df$Topic <- factor(df$Topic, levels = topic_levels)
df <- df %>% filter(!is.na(Topic)) %>% arrange(Topic, desc(CombinedScore))  
df <- df %>% mutate(x_order = factor(seq_len(n()), levels = seq_len(n())))

topic_tab <- table(df$Topic)
topic_breaks <- cumsum(as.numeric(topic_tab)) + 0.5
topic_breaks <- topic_breaks[-length(topic_breaks)]
x_labels <- setNames(df$Motif, df$x_order)

plot <- ggplot(df, aes(x = x_order, y = CombinedScore, fill = Topic)) +
  geom_col(width = 0.9) +
  geom_vline(xintercept = topic_breaks, linetype = 2, colour = "grey70") +
  scale_x_discrete(labels = x_labels) +
  scale_fill_manual(
    values = topic_cols,
    breaks = topic_levels,
    limits = topic_levels,
    drop = FALSE
  ) +
  labs(x = "Motif", y = "TF Score") +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 10),
    panel.grid.major.x = element_blank(),
    plot.margin = margin(80, 10, 20, 40)
  )

ggsave("../figures/Figure_E2/Figure_E2d.pdf", width = 10, height = 6, plot)

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
#   [1] dplyr_1.2.1         ggplot2_4.0.2       data.table_1.18.2.1
# 
# loaded via a namespace (and not attached):
#   [1] labeling_0.4.3     RColorBrewer_1.1-3 R6_2.6.1           tidyselect_1.2.1   farver_2.1.2       magrittr_2.0.5     gtable_0.3.6      
# [8] glue_1.8.0         tibble_3.3.1       pkgconfig_2.0.3    generics_0.1.4     lifecycle_1.0.5    cli_3.6.6          S7_0.2.1          
# [15] scales_1.4.0       grid_4.5.2         vctrs_0.7.2        textshaping_1.0.5  withr_3.0.2        systemfonts_1.3.2  compiler_4.5.2    
# [22] rstudioapi_0.18.0  tools_4.5.2        ragg_1.5.0         pillar_1.11.1      rlang_1.2.0     