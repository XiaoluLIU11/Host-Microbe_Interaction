### for RNA-seq 

rm(list=ls())
library(ggplot2)
library(dplyr)
library(factoextra)
library(pheatmap)
library(RColorBrewer)
library(reshape2)
library(data.table)
library(VennDiagram)
library(xlsx)
library(stringr)

set.seed(12334)
#source("D:/works/pg_omics/all_batch_readcounts/annotation.R")
#source("C:/Users/amin/Nutstore/1/拉曼数据分析平台/R code pipeline/commonAPI.R")
name.rank <- c("Wt_1","Wt_2","Wt_3","Ss_1","Ss_2","Ss_3","Pg_1","Pg_2","Pg_3","Sn_1","Sn_2","Sn_3")
pca_group <- factor(c(rep("Wt",3),rep("Ss",3),rep("Pg",3),rep("Sn",3)),levels = c("Wt","Ss","Pg","Sn"))


rna.m <- read.table("D:\\works\\pg_omics\\Rnaseq/DEGsList/Union_for_cluster_genename_sample.xls",sep = "\t",header = T,row.names = 1)
rna.1 <- rna.m %>% t() %>% .[c(1:12),]

row.names(rna.1) %>% 
  gsub("F2","",.) %>% 
  gsub("Sa","Ss",.) %>% 
  gsub("C","Wt",.) ->row.names(rna.1)

rna.2 <- rna.1[match(name.rank,row.names(rna.1)),]
rna.3 <- apply(rna.2, 2, as.numeric)

maxs <- apply(rna.3, 2, max)
mins <- apply(rna.3, 2, min)

rna.ok <- scale(rna.3,center = mins, scale = maxs - mins)
row.names(rna.ok) <- row.names(rna.2)

rna.pca_dat <- prcomp(rna.ok, scale = TRUE) 
pca_Variance <- round((rna.pca_dat$sdev^2/sum(rna.pca_dat$sdev^2)) * 100,2)
rna.p <-
  fviz_pca_ind(rna.pca_dat, repel = T,geom.ind = c('point','text'),
               col.ind= pca_group, mean.point=F,
               addEllipses = T, legend.title="Groups",
               ellipse.type="confidence", ellipse.level=0.9,
               palette = c("#BDBDBDFF", "#42B540FF", "#A50026","#4575B4"))+ 
  theme_classic()+
  labs(x = paste0("PC1 (",pca_Variance[1],"% explained var.)"), 
       y = paste0("PC2 (",pca_Variance[2],"% explained var.)" ))+
  theme_test()+
  ggtitle('RNA-Seq PCA plot')
#ggsave(plot = rna.p,filename = "D:/works/pg_omics/res/rna.pca.png",device = "png",height = 10,width = 12,units = "cm",dpi = 600,limitsize = F)
#output_ppt_m(rna.p,name = "D:/works/pg_omics/plot.pptx")

####Pheatmap
rna.m <- read.table("D:\\works\\pg_omics\\Rnaseq/DEGsList/Union_for_cluster_genename_sample.xls",
                    sep = "\t",header = T,row.names = 1) %>%  .[,c(7:9,4:6,1:3,10:12)]

colnames(rna.m) <- c("Wt_1","Wt_2","Wt_3",
                     "Ss_1","Ss_2","Ss_3",
                     "Pg_1","Pg_2","Pg_3",
                     "Sn_1","Sn_2","Sn_3")

#tiff(filename = 'D:/works/pg_omics/res/rna.dis.png',width = 6, height = 5, units = 'cm', res = 600, pointsize = 5)
t(rna.m[,1:12]) %>% dist(.,method = "canberra") %>% hclust() %>% plot()
#dev.off()
#output_ppt_m(t(rna.m[,1:12]) %>% dist(.,method = "canberra") %>% hclust() %>% plot(),name = "D:/works/pg_omics/plot.pptx")


rna.group <- read.table("D:\\works\\pg_omics\\Rnaseq/DEGsList/Union_for_cluster_genename_group.xls",
                        sep = "\t",header = T,row.names = 1) %>% .[,1:4]
colnames(rna.group) <- c("Pg","Ss","Wt","Sn")
rna.group <- rna.group [,c(3,2,1,4)]

pheatmap::pheatmap(rna.group,scale = "row",show_rownames = F,cutree_rows = 3,cluster_cols = F,
                   treeheight_row = 0, treeheight_col = 2,
                   colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(100),
                   height = 8.26,width = 8.08)-> degs.phm #,filename = "D:/works/pg_omics/res/DEGs.heatmap.png",show_colnames = F) 
#output_ppt_m(degs.phm,name = "D:/works/pg_omics/plot.pptx")

row_clusters <- cutree(degs.phm$tree_row, k = 3) 
row_clusters %>% data.frame() -> genes_cluster
rna.group[row.names(genes_cluster) %in% row.names(rna.group),] -> degs.log.re 
cbind(degs.log.re,genes_cluster) -> trends.df
trends.df <- trends.df[,c(3,2,1,4,5)]
colnames(trends.df)[5] <- "cluster"

trends.df %>% dplyr::filter(cluster == 1) %>% .[,-5]-> c1
c1$gene <- row.names(c1)
c1_new = melt(c1)
c1_new$variable <- factor(c1_new$variable,levels = c("Wt","Ss","Pg","Sn"))

ggplot(c1_new,aes(variable, value, group=gene)) +# geom_line(color="gray90",size=0.8) + 
  geom_hline(yintercept =0,linetype=2) +
  stat_summary(aes(group=1),fun.y=mean, geom="line", size=1.2, color="#c51b7d") + 
  stat_summary(aes(group=1),fun.data = mean_se, geom="errorbar", size=0.2, color="black",width = .4) + 
  theme_bw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text = element_text(size=8, face = "bold"),
        strip.text = element_text(size = 8, face = "bold")) -> c1trends
#output_ppt_m(c1trends,name = "D:/works/pg_omics/plot.pptx")


### add venn
color <- brewer.pal(4, "PiYG")
#color <- c("#F1B6DA","#D01C8B","#4DAC26","#B8E186")

PgvsWt_up <- fread("D:/works/pg_omics/20220531 转录组分析结果/Batch3/results/7.DiffExprAnalysis/7.1.DEGsList/PgvsC/PgvsC.DEG_up.xls") %>% .$Gene_id
PgvsWt_dn <- fread("D:/works/pg_omics/20220531 转录组分析结果/Batch3/results/7.DiffExprAnalysis/7.1.DEGsList/PgvsC/PgvsC.DEG_down.xls") %>% .$Gene_id

SsvsWt_up <- fread("D:/works/pg_omics/20220531 转录组分析结果/Batch3/results/7.DiffExprAnalysis/7.1.DEGsList/SavsC/SavsC.DEG_up.xls") %>% .$Gene_id
SsvsWt_dn <- fread("D:/works/pg_omics/20220531 转录组分析结果/Batch3/results/7.DiffExprAnalysis/7.1.DEGsList/SavsC/SavsC.DEG_down.xls") %>% .$Gene_id

SnvsPg_up <- fread("D:/works/pg_omics/20220531 转录组分析结果/Batch3/results/7.DiffExprAnalysis/7.1.DEGsList/SnvsPg/SnvsPg.DEG_up.xls") %>% .$Gene_id
SnvsPg_dn <- fread("D:/works/pg_omics/20220531 转录组分析结果/Batch3/results/7.DiffExprAnalysis/7.1.DEGsList/SnvsPg/SnvsPg.DEG_down.xls") %>% .$Gene_id

venn.p <- venn.diagram(
  x = list(PgvsWt_up,SnvsPg_up,SnvsPg_dn, PgvsWt_dn),
  category.names = c("PgvsWt_up" , "SnvsPg_up","SnvsPg_dn","PgvsWt_dn"),
  filename =NULL ,output=F, #'venn.png'
  resolution = 400,fill = color,
  lwd = 2
);grid.draw(venn.p)

## function
#enrichplotdata <- fread("D:\\works\\pg_omics\\res\\batch_care_degs\\enrichplot.txt")
enrichplotdata <- read.xlsx("D:\\works\\pg_omics\\res\\batch_care_degs\\General GO enrichment.xlsx",sheetIndex = 2)
enrichplotdata$Group <- factor(enrichplotdata$Group,c("Ss_vs_Wt","Pg_vs_Wt","Sn_vs_Pg"))
enrichplotdata$GeneRatio <- enrichplotdata$Hits/enrichplotdata$Total
plot_data <- enrichplotdata


color_map <- c("Up" = "red", "Down" = "blue")
enrichplot <- ggplot(plot_data, aes(x = Term, y = Group, size = -log(pvalue), color = Trends, alpha = GeneRatio)) +
  geom_point() +
  scale_color_manual(values = color_map) +  # 设置颜色映射
  scale_alpha_continuous(range = c(0.2, 1)) +
  scale_size_continuous(range = c(1, 9)) +
  scale_x_discrete(labels = function(x) str_wrap(x, width = 20)) +
  labs(x = "", y = "Group", size = "-log(pvalue)") +
  theme_classic() +
  coord_flip()+
  guides(color = guide_legend(override.aes = list(alpha = 1))) +  # 控制颜色的深浅程度
  theme(axis.title.x = element_blank(),  # 隐藏X轴标题
        axis.text.x = element_text(angle = 90, hjust = 1))+  # 调整X轴文字角度和对齐方式
  scale_x_discrete(position = "bottom", labels=function(x) str_wrap(x, width=45))
