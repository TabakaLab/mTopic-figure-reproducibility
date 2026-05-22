require(data.table)
require(ggplot2)
library(gridExtra)
library(this.path)

setwd(this.dir())

dir.create("../figures/Figure_E3", recursive = T)

#ATAC
topics_atac<-fread("../data/p22_atac_topics.csv")
rna_signatures_atac<-fread("../data/p22_atac_signatures_rna.csv")
atac_signatures<-fread("../data/p22_atac_signatures_atac.csv")
dominant_topics_atac<-fread("../data/p22_atac_dominant_topics.csv")

#H3K4me3
topics_H3K4me3<-fread("../data/p22_h3k4me3_topics.csv")
rna_signatures_H3K4me3<-fread("../data/p22_h3k4me3_signatures_rna.csv")
H3K4me3_signatures<-fread("../data/p22_h3k4me3_signatures_modality.csv")
dominant_topics_H3K4me3<-fread("../data/p22_h3k4me3_dominant_topics.csv")

#H3K27me3
topics_h3k27me3<-fread("../data/p22_h3k27me3_topics.csv")
rna_signatures_h3k27me3<-fread("../data/p22_h3k27me3_signatures_rna.csv")
h3k27me3_signatures<-fread("../data/p22_h3k27me3_signatures_modality.csv")
dominant_topics_h3k27me3<-fread("../data/p22_h3k27me3_dominant_topics.csv")

#H3K27ac
topics_h3k27ac<-fread("../data/p22_h3k27ac_topics.csv")
rna_signatures_h3k27ac<-fread("../data/p22_h3k27ac_signatures_rna.csv")
h3k27ac_signatures<-fread("../data/p22_h3k27ac_signatures_modality.csv")
dominant_topics_h3k27ac<-fread("../data/p22_h3k27ac_dominant_topics.csv")

################################################################################

title<-"Corpus collosum"
rna_atac<-rna_signatures_atac[,"topic_31"]
df_1<-data.frame(row.names = rna_signatures_atac$V1,TopicScore=rna_atac)

rna_H3K4me3<-rna_signatures_H3K4me3[,"topic_28"]
df_2<-data.frame(row.names = rna_signatures_H3K4me3$V1,TopicScore=rna_H3K4me3)

df<-cbind(df_1[intersect(rownames(df_1),rownames(df_2)),],
          df_2[intersect(rownames(df_1),rownames(df_2)),])
pc1<-cor(df[,1],df[,2],method = "kendall")
colnames(df)<-c("M1","M2")

size=18

df[,1]<-df[,1]/sum(df[,1])
df[,2]<-df[,2]/sum(df[,2])

p1<-ggplot(df, aes(x=M1, y=M2))+
  geom_point(size=1.0)+
  theme(panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black",size = 0.5),
        text=element_text(size=size),
        axis.ticks = element_line(size = 0.5,colour='black'),
        axis.text.x = element_text(size=size,colour='black', hjust = 0.5),
        axis.text.y = element_text(size=size,colour='black'),
        plot.title=element_text(size=size,colour='black'),
        axis.title.x=element_text(size=size,colour='black'),
        axis.title.y=element_text(size=size,colour='black'),
        legend.position="none")+
  labs(title=bquote(.(title) ~ tau == .(sprintf("%.2f", pc1))),
    x = "RNA probabilities, ATAC+RNA",
    y = "RNA probabilities, H3K4me3+RNA"
  )

########################################################################

rna_h3k27me3<-rna_signatures_h3k27me3[,"topic_14"]
df_2<-data.frame(row.names = rna_signatures_h3k27me3$V1,TopicScore=rna_h3k27me3)
df<-cbind(df_1[intersect(rownames(df_1),rownames(df_2)),],
          df_2[intersect(rownames(df_1),rownames(df_2)),])
pc2<-cor(df[,1],df[,2],method = "kendall")
colnames(df)<-c("M1","M2")
df[,1]<-df[,1]/sum(df[,1])
df[,2]<-df[,2]/sum(df[,2])
p2<-ggplot(df, aes(x=M1, y=M2))+
  geom_point(size=1.0)+
  #facet_wrap(~Var2,  ncol=10)+
  theme(panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black",size = 0.5),
        text=element_text(size=size),
        axis.ticks = element_line(size = 0.5,colour='black'),
        axis.text.x = element_text(size=size,colour='black', hjust = 0.5),
        axis.text.y = element_text(size=size,colour='black'),
        plot.title=element_text(size=size,colour='black'),
        axis.title.x=element_text(size=size,colour='black'),
        axis.title.y=element_text(size=size,colour='black'),
        legend.position="none")+
  labs(title=bquote(.(title) ~ tau == .(sprintf("%.2f", pc2))),
    x = "RNA probabilities, ATAC+RNA",
    y = "RNA probabilities, H3K27me3+RNA"
  )

################################################

rna_h3k27ac<-rna_signatures_h3k27ac[,"topic_22"]
df_2<-data.frame(row.names = rna_signatures_h3k27ac$V1,TopicScore=rna_h3k27ac)
df<-cbind(df_1[intersect(rownames(df_1),rownames(df_2)),],
          df_2[intersect(rownames(df_1),rownames(df_2)),])
pc3<-cor(df[,1],df[,2],method = "kendall")
colnames(df)<-c("M1","M2")
df[,1]<-df[,1]/sum(df[,1])
df[,2]<-df[,2]/sum(df[,2])
p3<-ggplot(df, aes(x=M1, y=M2))+
  geom_point(size=1.0)+
  #facet_wrap(~Var2,  ncol=10)+
  theme(panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black",size = 0.5),
        text=element_text(size=size),
        axis.ticks = element_line(size = 0.5,colour='black'),
        axis.text.x = element_text(size=size,colour='black', hjust = 0.5),
        axis.text.y = element_text(size=size,colour='black'),
        plot.title=element_text(size=size,colour='black'),
        axis.title.x=element_text(size=size,colour='black'),
        axis.title.y=element_text(size=size,colour='black'),
        legend.position="none")+
  labs(title=bquote(.(title) ~ tau == .(sprintf("%.2f", pc3))),
    x = "RNA probabilities, ATAC+RNA",
    y = "RNA probabilities, H3K27ac+RNA"
  )

########################################################################

title<-"Lateral ventricle"
rna_atac<-rna_signatures_atac[,"topic_43"]
df_1<-data.frame(row.names = rna_signatures_atac$V1,TopicScore=rna_atac)
rna_H3K4me3<-rna_signatures_H3K4me3[,"topic_25"]
df_2<-data.frame(row.names = rna_signatures_H3K4me3$V1,TopicScore=rna_H3K4me3)
df<-cbind(df_1[intersect(rownames(df_1),rownames(df_2)),],
          df_2[intersect(rownames(df_1),rownames(df_2)),])
pc4<-cor(df[,1],df[,2],method = "kendall")
colnames(df)<-c("M1","M2")
df[,1]<-df[,1]/sum(df[,1])
df[,2]<-df[,2]/sum(df[,2])
p4<-ggplot(df, aes(x=M1, y=M2))+
  geom_point(size=1.0)+
  #facet_wrap(~Var2,  ncol=10)+
  theme(panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black",size = 0.5),
        text=element_text(size=size),
        axis.ticks = element_line(size = 0.5,colour='black'),
        axis.text.x = element_text(size=size,colour='black', hjust = 0.5),
        axis.text.y = element_text(size=size,colour='black'),
        plot.title=element_text(size=size,colour='black'),
        axis.title.x=element_text(size=size,colour='black'),
        axis.title.y=element_text(size=size,colour='black'),
        legend.position="none")+
  labs(title=bquote(.(title) ~ tau == .(sprintf("%.2f", pc4))),
    x = "RNA probabilities, ATAC+RNA",
    y = "RNA probabilities, H3K4me3+RNA"
  )

#######################################

rna_h3k27me3<-rna_signatures_h3k27me3[,"topic_43"]
df_2<-data.frame(row.names = rna_signatures_h3k27me3$V1,TopicScore=rna_h3k27me3)
df<-cbind(df_1[intersect(rownames(df_1),rownames(df_2)),],df_2[intersect(rownames(df_1),rownames(df_2)),])
pc5<-cor(df[,1],df[,2],method = "kendall")
colnames(df)<-c("M1","M2")
df[,1]<-df[,1]/sum(df[,1])
df[,2]<-df[,2]/sum(df[,2])
p5<-ggplot(df, aes(x=M1, y=M2))+
  geom_point(size=1.0)+
  #facet_wrap(~Var2,  ncol=10)+
  theme(panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black",size = 0.5),
        text=element_text(size=size),
        axis.ticks = element_line(size = 0.5,colour='black'),
        axis.text.x = element_text(size=size,colour='black', hjust = 0.5),
        axis.text.y = element_text(size=size,colour='black'),
        plot.title=element_text(size=size,colour='black'),
        axis.title.x=element_text(size=size,colour='black'),
        axis.title.y=element_text(size=size,colour='black'),
        legend.position="none")+
  labs(title=bquote(.(title) ~ tau == .(sprintf("%.2f", pc5))),
    x = "RNA probabilities, ATAC+RNA",
    y = "RNA probabilities, H3K27me3+RNA"
  )


rna_h3k27ac<-rna_signatures_h3k27ac[,"topic_6"]
df_2<-data.frame(row.names = rna_signatures_h3k27ac$V1,TopicScore=rna_h3k27ac)
df<-cbind(df_1[intersect(rownames(df_1),rownames(df_2)),],
          df_2[intersect(rownames(df_1),rownames(df_2)),])
pc6<-cor(df[,1],df[,2],method = "pearson")
colnames(df)<-c("M1","M2")
df[,1]<-df[,1]/sum(df[,1])
df[,2]<-df[,2]/sum(df[,2])
p6<-ggplot(df, aes(x=M1, y=M2))+
  geom_point(size=1.0)+
  #facet_wrap(~Var2,  ncol=10)+
  theme(panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black",size = 0.5),
        text=element_text(size=size),
        axis.ticks = element_line(size = 0.5,colour='black'),
        axis.text.x = element_text(size=size,colour='black', hjust = 0.5),
        axis.text.y = element_text(size=size,colour='black'),
        plot.title=element_text(size=size,colour='black'),
        axis.title.x=element_text(size=size,colour='black'),
        axis.title.y=element_text(size=size,colour='black'),
        legend.position="none")+
  labs(title=bquote(.(title) ~ tau == .(sprintf("%.2f", pc6))),
    x = "RNA probabilities, ATAC+RNA",
    y = "RNA probabilities, H3K27ac+RNA"
  )

pdf(file = paste0("../figures/Figure_E3/Figure_E3d.pdf"), width = 10, height=15)
grid.arrange(p4, p1, p5, p2, p6, p3, ncol=2, nrow =3)
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
#   [1] gridExtra_2.3       ggplot2_4.0.2       data.table_1.18.2.1
# 
# loaded via a namespace (and not attached):
#   [1] labeling_0.4.3     RColorBrewer_1.1-3 R6_2.6.1           tidyselect_1.2.1   farver_2.1.2       magrittr_2.0.5     gtable_0.3.6      
# [8] glue_1.8.0         tibble_3.3.1       pkgconfig_2.0.3    generics_0.1.4     dplyr_1.2.1        lifecycle_1.0.5    cli_3.6.6         
# [15] S7_0.2.1           scales_1.4.0       grid_4.5.2         vctrs_0.7.2        withr_3.0.2        compiler_4.5.2     rstudioapi_0.18.0 
# [22] tools_4.5.2        pillar_1.11.1      rlang_1.2.0     