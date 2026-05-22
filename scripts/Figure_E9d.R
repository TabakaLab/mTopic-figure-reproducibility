require(data.table)
library(ggrepel)
library(gridExtra)
library(this.path)

setwd(this.dir())

dir.create("../figures/Figure_E9", recursive = T)

TFxGene <- readRDS(file="../data/PBMC_cyto_GRN.RDS")
RTE <- readRDS("../data/RTE.RDS")

genes <- TFxGene[which(rownames(TFxGene)=="EOMES--MA0800.1"),]
genes <- unlist(as.vector(genes))
genes <- sort(genes,decreasing = T)
df <- data.frame(Rank=1:length(genes), Genes=names(genes), Score=genes)
df <- df[1:2000,]
df_RTE <- df[names(RTE$EOMES),]

n <- 10
p <- 5
df_RTE_l <- df_RTE[1:n,]
tmp <- df_RTE[n+1:nrow(df_RTE),]
tmp <- tmp[sample(1:nrow(tmp),p),]
df_RTE_l <- rbind(df_RTE_l,tmp)

size <- 22
p1 <- ggplot(df, aes(x=Rank, y=Score)) +
  geom_point(size=0.2)+
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
        legend.position="none")+
  scale_x_continuous(breaks=seq(0,10000,1000)) +
  geom_point(data=df_RTE, size=0.2,color="red") +
  geom_text_repel(force = 1, max.overlaps = Inf,
                  data = df_RTE_l,
                  aes(label = Genes),
                  size=5,segment.size = 0.2,colour="black")+
  labs(title="EOMES")

genes <- TFxGene[which(rownames(TFxGene)=="TBX21--MA0690.2"),]
genes <- unlist(as.vector(genes))
genes <- sort(genes,decreasing = T)
df <- data.frame(Rank=1:length(genes), Genes=names(genes), Score=genes)
df <- df[1:2000,]
df_RTE <- df[names(RTE$TBX21),]

n <- 10
p <- 5
df_RTE_l <- df_RTE[1:n,]
tmp <- df_RTE[n+1:nrow(df_RTE),]
tmp <- tmp[sample(1:nrow(tmp),p),]
df_RTE_l <- rbind(df_RTE_l,tmp)

p2 <- ggplot(df, aes(x=Rank, y=Score)) +
  geom_point(size=0.2)+
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
        legend.position="none")+
  scale_x_continuous(breaks=seq(0,10000,1000)) +
  geom_point(data=df_RTE,size=0.2,color="green") +
  geom_text_repel(force = 1,max.overlaps = Inf,
                  data = df_RTE_l,
                  aes(label = Genes),
                  #sapply(strsplit(as.character(tmp$Var1[tmp$Rank %in% 1:5 ]),"--"), `[`, 1)),
                  size=5,segment.size = 0.2,colour="black") +
  labs(title="TBX21")



genes <- TFxGene[which(rownames(TFxGene)=="RORA--MA0072.1"),]
genes <- unlist(as.vector(genes))
genes <- sort(genes,decreasing = T)
df <- data.frame(Rank=1:length(genes), Genes=names(genes), Score=genes)
df <- df[1:2000,]
df_RTE <- df[names(RTE$RORA),]

n <- 10
p <- 5
df_RTE_l <- df_RTE[1:n,]
tmp <- df_RTE[n+1:nrow(df_RTE),]
tmp <- tmp[sample(1:nrow(tmp),p),]
df_RTE_l <- rbind(df_RTE_l,tmp)

p3 <- ggplot(df, aes(x=Rank, y=Score)) +
  geom_point(size=0.2) +
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
        legend.position="none") +
  scale_x_continuous(breaks=seq(0,10000,1000)) +
  geom_point(data=df_RTE,size=0.2,color="blue") +
  geom_text_repel(force = 1,max.overlaps = Inf,
                  data = df_RTE_l,
                  aes(label = Genes),
                  #sapply(strsplit(as.character(tmp$Var1[tmp$Rank %in% 1:5 ]),"--"), `[`, 1)),
                  size=5, segment.size = 0.2, colour="black") +
  labs(title="RORA")


pdf(file = paste0("../figures/Figure_E9/Figure_E9d.pdf"), width = 20, height=8)
grid.arrange(p2, p1, p3, ncol=3, nrow =1)
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
#   [1] gridExtra_2.3       ggrepel_0.9.8       ggplot2_4.0.2       data.table_1.18.2.1
# 
# loaded via a namespace (and not attached):
#   [1] labeling_0.4.3     RColorBrewer_1.1-3 R6_2.6.1           tidyselect_1.2.1   farver_2.1.2       magrittr_2.0.5     gtable_0.3.6      
# [8] glue_1.8.0         tibble_3.3.1       pkgconfig_2.0.3    generics_0.1.4     dplyr_1.2.1        lifecycle_1.0.5    cli_3.6.6         
# [15] S7_0.2.1           scales_1.4.0       grid_4.5.2         vctrs_0.7.2        withr_3.0.2        compiler_4.5.2     rstudioapi_0.18.0 
# [22] tools_4.5.2        pillar_1.11.1      Rcpp_1.1.1         rlang_1.2.0      