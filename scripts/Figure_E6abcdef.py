require(data.table)
library(ggrepel)
library(gridExtra)
library(this.path)

setwd(this.dir())
set.seed(1)

dir.create("../figures/Figure_E10", recursive = T, showWarnings = FALSE)

TFxGene <- readRDS(file="../data/PBMC_cyto_GRN.RDS")
RTE <- readRDS("../data/RTE.RDS")

size <- 7
label_size <- 7

genes <- TFxGene[which(rownames(TFxGene)=="EOMES--MA0800.1"),]
genes <- unlist(as.vector(genes))
genes <- sort(genes,decreasing = T)
df <- data.frame(Rank=1:length(genes), Genes=names(genes), Score=genes)
df <- df[1:2000,]
df_RTE <- df[names(RTE$EOMES),]

n <- 10
p <- 5
df_RTE_l <- df_RTE[1:n,]
tmp <- df_RTE[(n+1):nrow(df_RTE),]
tmp <- tmp[sample(1:nrow(tmp),p),]
df_RTE_l <- rbind(df_RTE_l,tmp)

p1 <- ggplot(df, aes(x=Rank, y=Score)) +
  geom_point(size=0.2)+
  theme(panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black",linewidth = 0.5),
        text=element_text(size=size),
        axis.ticks = element_line(linewidth = 0.5,colour='black'),
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
                  size=label_size,segment.size = 0.2,colour="black")+
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
tmp <- df_RTE[(n+1):nrow(df_RTE),]
tmp <- tmp[sample(1:nrow(tmp),p),]
df_RTE_l <- rbind(df_RTE_l,tmp)

p2 <- ggplot(df, aes(x=Rank, y=Score)) +
  geom_point(size=0.2)+
  theme(panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black",linewidth = 0.5),
        text=element_text(size=size),
        axis.ticks = element_line(linewidth = 0.5,colour='black'),
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
                  size=label_size,segment.size = 0.2,colour="black") +
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
tmp <- df_RTE[(n+1):nrow(df_RTE),]
tmp <- tmp[sample(1:nrow(tmp),p),]
df_RTE_l <- rbind(df_RTE_l,tmp)

p3 <- ggplot(df, aes(x=Rank, y=Score)) +
  geom_point(size=0.2) +
  theme(panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black",linewidth = 0.5),
        text=element_text(size=size),
        axis.ticks = element_line(linewidth = 0.5,colour='black'),
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
                  size=label_size, segment.size = 0.2, colour="black") +
  labs(title="RORA")

pdf(file = "../figures/Figure_E10/Figure_E10d.pdf", width = 20, height=8)
grid.arrange(p2, p1, p3, ncol=3, nrow =1)
dev.off()