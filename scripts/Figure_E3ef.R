require(data.table)
require(stringr)
require(dplyr)
require(ggplot2)
library(gridExtra)
library(this.path)

setwd(this.dir())

dir.create("../figures/Figure_E3", recursive = T)

PromList_atac<-readRDS("../data/atac_Pr_list.RDS")
DistList_atac<-readRDS("../data/atac_Di_list.RDS")

PromList_h3K27ac<-readRDS("../data/H3K27ac_Pr_list.RDS")
DistList_h3K27ac<-readRDS("../data/H3K27ac_Di_list.RDS")

PromList_h3k27me3<-readRDS("../data/H3k27me3_Pr_list.RDS")
DistList_h3k27me3<-readRDS("../data/H3k27me3_Di_list.RDS")

PromList_h3K4me3<-readRDS("../data/H3K4me3_Pr_list.RDS")
DistList_h3K4me3<-readRDS("../data/H3K4me3_Di_list.RDS")

GetRE<-function(nTopics=50,REList,nTop=1000){
  newREList<-list()
  for(nt in 1:nTopics){
    tmp<-rbind(REList[[nt]])
    tmp<-tmp %>%
      group_by(GeneID) %>%
      mutate(Sum_Chr_score = sum(Chr_score, na.rm = TRUE)) %>%
      ungroup()
    tmp<-tmp[!duplicated(tmp$GeneID),]
    tmp<-tmp[order(tmp$RNA_score,decreasing = T),]
    newREList[[nt]]<-tmp[1:nTop,]
  }
  return(newREList) 
}

nTop<-1000

P_atac_REList<-GetRE(nTopics=50,PromList_atac,nTop=nTop)
D_atac_REList<-GetRE(nTopics=50,DistList_atac,nTop=nTop)

P_h3K27ac_REList<-GetRE(nTopics=50,PromList_h3K27ac,nTop=nTop)
D_h3K27ac_REList<-GetRE(nTopics=50,DistList_h3K27ac,nTop=nTop)


P_h3K4me3_REList<-GetRE(nTopics=50,PromList_h3K4me3,nTop=nTop)
D_h3K4me3_REList<-GetRE(nTopics=50,DistList_h3K4me3,nTop=nTop)

P_h3k27me3_REList<-GetRE(nTopics=50,PromList_h3k27me3,nTop=nTop)
D_h3k27me3_REList<-GetRE(nTopics=50,DistList_h3k27me3,nTop=nTop)

Prom_ATAC<-do.call(rbind,P_atac_REList)
Prom_h3K27ac<-do.call(rbind,P_h3K27ac_REList)
Prom_h3k27me3<-do.call(rbind,P_h3k27me3_REList)
Prom_h3K4me3<-do.call(rbind,P_h3K4me3_REList)

Prom_ATAC$RNA_score<-Prom_ATAC$RNA_score/max(Prom_ATAC$RNA_score)
Prom_ATAC$Sum_Chr_score<-Prom_ATAC$Sum_Chr_score/max(Prom_ATAC$Sum_Chr_score)
Prom_h3K27ac$RNA_score<-Prom_h3K27ac$RNA_score/max(Prom_h3K27ac$RNA_score)
Prom_h3K27ac$Sum_Chr_score<-Prom_h3K27ac$Sum_Chr_score/max(Prom_h3K27ac$Sum_Chr_score)
Prom_h3k27me3$RNA_score<-Prom_h3k27me3$RNA_score/max(Prom_h3k27me3$RNA_score)
Prom_h3k27me3$Sum_Chr_score<-Prom_h3k27me3$Sum_Chr_score/max(Prom_h3k27me3$Sum_Chr_score)
Prom_h3K4me3$RNA_score<-Prom_h3K4me3$RNA_score/max(Prom_h3K4me3$RNA_score)
Prom_h3K4me3$Sum_Chr_score<-Prom_h3K4me3$Sum_Chr_score/max(Prom_h3K4me3$Sum_Chr_score)

Dist_ATAC<-do.call(rbind,D_atac_REList)
Dist_h3K27ac<-do.call(rbind,D_h3K27ac_REList)
Dist_h3k27me3<-do.call(rbind,D_h3k27me3_REList)
Dist_h3K4me3<-do.call(rbind,D_h3K4me3_REList)

Dist_ATAC$RNA_score<-Dist_ATAC$RNA_score/max(Dist_ATAC$RNA_score)
Dist_ATAC$Sum_Chr_score<-Dist_ATAC$Sum_Chr_score/max(Dist_ATAC$Sum_Chr_score)
Dist_h3K27ac$RNA_score<-Dist_h3K27ac$RNA_score/max(Dist_h3K27ac$RNA_score)
Dist_h3K27ac$Sum_Chr_score<-Dist_h3K27ac$Sum_Chr_score/max(Dist_h3K27ac$Sum_Chr_score)
Dist_h3k27me3$RNA_score<-Dist_h3k27me3$RNA_score/max(Dist_h3k27me3$RNA_score)
Dist_h3k27me3$Sum_Chr_score<-Dist_h3k27me3$Sum_Chr_score/max(Dist_h3k27me3$Sum_Chr_score)
Dist_h3K4me3$RNA_score<-Dist_h3K4me3$RNA_score/max(Dist_h3K4me3$RNA_score)
Dist_h3K4me3$Sum_Chr_score<-Dist_h3K4me3$Sum_Chr_score/max(Dist_h3K4me3$Sum_Chr_score)

Prom<-rbind(cbind(Prom_ATAC,Type="ATAC"),
            cbind(Prom_h3K27ac,Type="H3K27ac"),
            cbind(Prom_h3K4me3,Type="H3K4me3"),
            cbind(Prom_h3k27me3,Type="H3K27me3")
)
Dist<-rbind(cbind(Dist_ATAC,Type="ATAC"),
            cbind(Dist_h3K27ac,Type="H3K27ac"),
            cbind(Dist_h3K4me3,Type="H3K4me3"),
            cbind(Dist_h3k27me3,Type="H3K27me3")
)

Prom <- Prom[setdiff(1:nrow(Prom),which(Prom$RNA_score==0.01 & Prom$Chr_score==0.01)),]
Dist <- Dist[setdiff(1:nrow(Dist),which(Dist$RNA_score==0.01 & Dist$Chr_score==0.01)),]

cols <- c(
  "H3K27me3" = "#7991A8",
  "ATAC"     = "#E1BB3E",
  "H3K27ac"  = "#B322AD",
  "H3K4me3"  = "#7EBA53"
)
size <- 18
png(file = paste0("../figures/Figure_E3/Figure_E3e_1.png"),width = 140,height=120,units="mm",res=300)
ggplot(Prom, aes(x = RNA_score, y = Sum_Chr_score, color = Type, fill = Type)) +
  geom_point(size = 0.2, alpha = 1) +
  scale_color_manual(values = cols) +
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
  xlab("RNA topic score") + 
  ylab("Chromatin topic score")
dev.off()

png(file = paste0("../figures/Figure_E3/Figure_E3e_2.png"),width = 200,height=120,units="mm",res =300)
ggplot(Dist, aes(x = RNA_score, y = Sum_Chr_score, color = Type, fill = Type)) +
  geom_point(size = 0.1, alpha = 1) +
  scale_color_manual(values = cols) +
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
        axis.title.y=element_text(size=size,colour='black'))+
  theme(legend.title = element_blank())+
  xlab("RNA topic score") + 
  ylab("Chromatin cumulative \n topic score")+guides(color = guide_legend(override.aes = list(size = 6)))
dev.off()



# User-specified RNA_score range
min_val <- 0.05
max_val <- 0.1 # adjust as needed

# Filter within the specified range
Prom_subset <- Prom %>%
  filter(RNA_score >= min_val, RNA_score <= max_val, Sum_Chr_score > 0)

# Plot density of Sum_Chr_score for this RNA_score window
#png(file = paste0("ECDF_promoters",,"_Chr_RNA.png"),width = 140,height=120,units="mm",res =300)
p1 <- ggplot(Prom_subset, aes(x = Sum_Chr_score, color = Type)) +
  stat_ecdf(geom = "step", size = 1) +
  scale_color_manual(values = cols) +
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
  xlab("Chromatin topic score") + 
  ylab("Ecdf") +
  xlim(0, 0.5)

min_val <- 0.15
max_val <- 0.2 # adjust as needed

# Filter within the specified range
Prom_subset <- Prom %>%
  filter(RNA_score >= min_val, RNA_score <= max_val, Sum_Chr_score > 0)

# Plot density of Sum_Chr_score for this RNA_score window
p2 <- ggplot(Prom_subset, aes(x = Sum_Chr_score, color = Type)) +
  stat_ecdf(geom = "step", size = 1) +
  scale_color_manual(values = cols) +
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
  xlab("Chromatin topic score") + 
  ylab("Ecdf") +
  xlim(0, 0.5)

# User-specified RNA_score range
min_val <- 0.05
max_val <- 0.1 # adjust as needed

# Filter within the specified range
Prom_subset <- Dist %>%
  filter(RNA_score >= min_val, RNA_score <= max_val, Sum_Chr_score > 0)

# Plot density of Sum_Chr_score for this RNA_score window
d1 <- ggplot(Prom_subset, aes(x = Sum_Chr_score, color = Type)) +
  stat_ecdf(geom = "step", size = 1) +
  # geom_density(alpha = 0.6) +
  scale_color_manual(values = cols) +
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
  xlab("Chromatin topic score") + 
  ylab("Ecdf") +
  xlim(0, 0.3)

min_val <- 0.15
max_val <- 0.2 # adjust as needed

# Filter within the specified range
Prom_subset <- Dist %>%
  filter(RNA_score >= min_val, RNA_score <= max_val, Sum_Chr_score > 0)

# Plot density of Sum_Chr_score for this RNA_score window
d2 <- ggplot(Prom_subset, aes(x = Sum_Chr_score, color = Type)) +
  stat_ecdf(geom = "step", size = 1) +
  #geom_density(alpha = 0.6) +
  scale_color_manual(values = cols) +
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
  xlab("Chromatin topic score") + 
  ylab("Ecdf") +
  xlim(0, 0.3)

pdf(file = paste0("../figures/Figure_E3/Figure_E3f.pdf"), width = 13, height=4)
grid.arrange(p1, p2, d1, d2, ncol=4, nrow =1)
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
#   [1] gridExtra_2.3       ggplot2_4.0.2       dplyr_1.2.1         stringr_1.6.0       data.table_1.18.2.1
# 
# loaded via a namespace (and not attached):
#   [1] labeling_0.4.3     RColorBrewer_1.1-3 R6_2.6.1           tidyselect_1.2.1   farver_2.1.2       magrittr_2.0.5     gtable_0.3.6      
# [8] glue_1.8.0         tibble_3.3.1       pkgconfig_2.0.3    generics_0.1.4     lifecycle_1.0.5    cli_3.6.6          S7_0.2.1          
# [15] scales_1.4.0       grid_4.5.2         vctrs_0.7.2        withr_3.0.2        compiler_4.5.2     rstudioapi_0.18.0  tools_4.5.2       
# [22] pillar_1.11.1      rlang_1.2.0        stringi_1.8.7  