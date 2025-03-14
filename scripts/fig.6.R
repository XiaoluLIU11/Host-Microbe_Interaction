rm(list = ls())
library(tidyr)
library(ggplot2)
library(ggpubr)
library(KEGGREST)
library(xlsx)
library(factoextra)
library(MetaboAnalystR)
library(stringr)
library(pheatmap)
library(data.table)
set.seed(12334)

###PCA

meta.m <- read.xlsx("D:/works/pg_omics/new_metabolome/metabolome/2.Basic_Analysis/Difference_analysis/sigMetabolitesSummary.xlsx",sheetIndex = 1,header = T,row.names =1)
meta.m <- meta.m[,c(11:22)] %>% t() %>% data.frame()
row.names(meta.m) %>% 
  gsub("C","Wt_",.) %>%
  gsub("Sa","Ss_",.) %>%
  gsub("Pg","Pg_",.) %>%
  gsub("Sn","Sn_",.) -> row.names(meta.m)

maxs <- apply(meta.m, 2, max)
mins <- apply(meta.m, 2, min)

meta.ok <- scale(meta.m,center = mins, scale = maxs - mins)
row.names(meta.ok) <- row.names(meta.m)

meta.pca_dat <- prcomp(meta.ok, scale = TRUE) 
pca_Variance <- round((meta.pca_dat$sdev^2/sum(meta.pca_dat$sdev^2)) * 100,2)
meta.p <- fviz_pca_ind(meta.pca_dat, repel = T,geom.ind = c('point','text'),
                       col.ind=pca_group, mean.point=F, 
                       addEllipses = T, legend.title="Groups", 
                       ellipse.type="confidence", ellipse.level=0.9, 
                       palette = c("#F8766D", "#00BA38", "#619CFF","black"))+ 
  theme_classic()+
  labs(x = paste0("PC1 (",pca_Variance[1],"% explained var.)"), 
       y = paste0("PC2 (",pca_Variance[2],"% explained var.)" ))+
  theme_test()+
  ggtitle('Metabolome PCA plot')


meta.raw <- read.xlsx("D:/works/pg_omics\\ALL_sample_data.xlsx",sheetIndex = 1)
meta.count <- meta.raw[,c(1,12:23)]
row.names(meta.count) <- meta.count$Index
meta.count <- meta.count[,-1]
colnames(meta.count) <- c("Wt_1","Wt_2","Wt_3","Ss_1","Ss_2","Ss_3","Pg_1","Pg_2","Pg_3","Sn_1","Sn_2","Sn_3")


t(meta.count)  %>% dist(.,method = "canberra") %>% hclust() %>% plot()



dems.box.df <- read.xlsx("D:/works/pg_omics/new_metabolome/metabolome/2.Basic_Analysis/Difference_analysis/sigMetabolitesSummary.xlsx",sheetIndex = 2)
dems.box.df %>%
  pivot_longer(cols = c(Wt, Ss, Pg, Sn), names_to = "Groups", values_to = "Value") -> dems.box.long
dems.box.long$Groups <- factor(dems.box.long$Groups,levels = c("Wt","Ss","Pg","Sn"))
dems.box.long$Value <- log(dems.box.long$Value)

p <- ggplot(dems.box.long, aes(x = str_wrap(Class.I, width = 20), y = Value, fill = Groups)) +
  geom_boxplot(outlier.shape = NA) +
  scale_fill_manual(values = c("#F8766D", "#00BA38", "#619CFF","grey"))+
  #scale_fill_manual(values = cal_palette("chaparral1"))+ #
  theme_classic()+
  ylim(c(3,32))+
  labs(title = "Pairwise Boxplot", x = "Class I", y = "Groups") +
  labs(x = "Metabolites Classes", y = "Log(Content)")+
  stat_compare_means(label = "p.signif", method = "kruskal.test", ref.group = NULL)+
  theme(axis.text.x = element_text(size = 8,angle = 60, hjust = 1,color = "black"))+
  ggtitle("Differential Expression Metabolites Content") 


###amino acid
dems.box.long %>% dplyr::filter(Class.I == "Amino acid and Its metabolites") -> aa.box.long

mean(dplyr::filter(aa.box.long,Groups == "Wt")$Value)
sd(dplyr::filter(aa.box.long,Groups == "Wt")$Value)

mean(dplyr::filter(aa.box.long,Groups == "Pg")$Value)
sd(dplyr::filter(aa.box.long,Groups == "Pg")$Value)

mean(dplyr::filter(aa.box.long,Groups == "Sn")$Value)
sd(dplyr::filter(aa.box.long,Groups == "Sn")$Value)

p <- ggplot(aa.box.long, aes(x = Groups, y = Value, fill = Groups)) +
  geom_boxplot(outlier.shape = NA) +
  scale_fill_manual(values = c("#F8766D", "#00BA38", "#619CFF","grey"))+
  theme_classic()+
  #ylim(c(3,32))+
  labs(title = paste0(aa.box.long$Class.I[1]), x = "Class I", y = "Groups") +
  labs(x = NULL, y = "Log(Content)")+
  stat_compare_means(comparisons = list( c("Wt", "Pg"), c("Wt", "Sn"),c("Ss", "Pg"), c("Pg", "Sn")), method = "wilcox.test")

# Fatty acyls
dems.box.long %>% dplyr::filter(Class.I == "Fatty acyls") -> aa.box.long
p <- ggplot(aa.box.long, aes(x = Groups, y = Value, fill = Groups)) +
  geom_boxplot(outlier.shape = NA) +
  scale_fill_manual(values = c("#F8766D", "#00BA38", "#619CFF","grey"))+
  theme_classic()+
  #ylim(c(3,32))+
  labs(title = paste0(aa.box.long$Class.I[1]), x = "Class I", y = "Groups") +
  labs(x = NULL, y = "Log(Content)")+
  stat_compare_means(comparisons = list(c("Wt", "Ss"), c("Wt", "Pg"), c("Wt", "Sn")), method = "wilcox.test")


###Glycerol phospholipids 

dems.box.long %>% dplyr::filter(Class.I == "Glycerol phospholipids") -> aa.box.long
p <- ggplot(aa.box.long, aes(x = Groups, y = Value, fill = Groups)) +
  geom_boxplot(outlier.shape = NA) +
  scale_fill_manual(values = c("#F8766D", "#00BA38", "#619CFF","grey"))+
  theme_classic()+
  #ylim(c(3,32))+
  labs(title = paste0(aa.box.long$Class.I[1]), x = "Class I", y = "Groups") +
  labs(x = NULL, y = "Log(Content)")+
  stat_compare_means(comparisons = list( c("Wt", "Sn"), c("Pg", "Sn")), method = "wilcox.test")


###Nucleotide And Its metabolites

dems.box.long %>% dplyr::filter(Class.I == "Nucleotide And Its metabolites") -> aa.box.long
p <- ggplot(aa.box.long, aes(x = Groups, y = Value, fill = Groups)) +
  geom_boxplot(outlier.shape = NA) +
  scale_fill_manual(values = c("#F8766D", "#00BA38", "#619CFF","grey"))+
  theme_classic()+
  #ylim(c(3,32))+
  labs(title = paste0(aa.box.long$Class.I[1]), x = "Class I", y = "Groups") +
  labs(x = NULL, y = "Log(Content)")+
  stat_compare_means(comparisons = list( c("Wt", "Sn"),c("Wt", "Pg"), c("Pg", "Ss"), c("Pg", "Sn")), method = "wilcox.test")




dems_dat <- read.xlsx(file = "D:/works/pg_omics/new_metabolome/metabolome/2.Basic_Analysis/Difference_analysis/sigMetabolitesSummary.xlsx",sheetIndex = 1)
colnames(dems_dat)[grep(colnames(dems_dat),pattern = "Type")]  -> group_names

trends = "up" # or down

for (i in 1:length(group_names)){
  dems_dat %>% dplyr::select(group_names[i]) %>%
    data.frame(dems_dat$HMDB,.) %>% .[.[,2] == trends,] %>%
    dplyr::filter(.,dems_dat.HMDB != "-") %>% 
    .[,1] -> dems
  
  mSet<-InitDataObjects("conc", "pathora", FALSE)
  mSet<-Setup.MapData(mSet, dems)
  mSet<-CrossReferencing(mSet, "hmdb")
  mSet<-SetKEGG.PathLib(mSet, "hsa", "current")
  mSet<-SetMetabolomeFilter(mSet, F)
  mSet<-CalculateOraScore(mSet, "rbc", "hyperg")
  mSet$api$ora.results %>% data.frame() -> msea_res
  rownames(mSet$analSet$ora.mat) -> msea_res$kegg_id
  #write.csv(x = msea_res,file = paste0("D:/works/pg_omics/res/dems_analysis2/",group_names[i],"_",trends,"_kegg_enrichment_res.csv"))
}

## enrichment plot
enrichplotdata <- read.xlsx("D:/works/pg_omics/res/dems_analysis2/kegg_enrichment_integration.xlsx",sheetIndex = 1)
enrichplotdata$Group <- factor(enrichplotdata$Group,c("Pg_vs_Wt","Sn_vs_Pg","Ss_vs_Wt"))
plot_data <- enrichplotdata
plot_data$GeneRatio <- plot_data$Hits/plot_data$Total
plot_data <- plot_data %>% dplyr::filter(FDR < 0.05)


color_map <- c("Up" = "red", "Down" = "blue")

enrichplot <- ggplot(plot_data, aes(x = Terms, y = Group, size = -log(FDR), color = Trends,alpha = GeneRatio)) +
  geom_point() +
  scale_color_manual(values = color_map) + 
  scale_alpha_continuous(range = c(0.2, 1)) +
  scale_size_continuous(range = c(1, 9)) +
  scale_x_discrete(labels = function(x) str_wrap(x, width = 20)) +
  labs(x = "", y = "Group", size = "-log(pvalue)") +
  theme_classic() +
  coord_flip()+
  guides(color = guide_legend(override.aes = list(alpha = 1))) + 
  theme(axis.title.x = element_blank(),  
        axis.text.x = element_text(angle = 90, hjust = 1))+ 
  scale_x_discrete(position = "bottom", labels=function(x) str_wrap(x, width=45))


### DEMs kegg heatmap
arg.m <- read.table("D:/works/pg_omics/res/dems_analysis2/arginine.txt",sep = "\t",row.names = 1,header = T)
arg.m <- arg.m[,c(4,8,12,16)]
pheatmap(arg.m,scale = "row",color = colorRampPalette(c("navyblue","white","firebrick3"))(200),
         cellwidth = 25,cellheight = 25,cutree_rows = 2,treeheight_col = 10) -> arg.phm


arac.m <- read.table("D:/works/pg_omics/res/dems_analysis2/arachidonic.txt",sep = "\t",row.names = 1,header = T,encoding = "UTF-8")
pheatmap(arac.m,scale = "row",color = colorRampPalette(c("navyblue","white","firebrick3"))(200),
         cellwidth = 25,cellheight = 25,treeheight_col = 10) -> arac.phm


### arg + gene expression level to heatmap
all_genes <- fread("D:\\works\\pg_omics\\20220531 转录组分析结果\\Batch3\\results\\5.GeneQuantification\\fpkm_genename.csv")
arg_bio_path <- keggGet("hsa00220")
arg_bio_path[[1]]$GENE -> gene.info 
genes <- unlist(lapply(gene.info,function(x) strsplit(x,";")))
gene.symbol <- genes[1:length(genes)%%3 == 2]

arg_bio_gene_df <- all_genes[all_genes$`Gene name` %in% gene.symbol,]
arg_bio_phm <- data.frame(row.names = arg_bio_gene_df$`Gene name`,arg_bio_gene_df[,c(2:13)])
arg_bio_phm$Pg <- rowMeans(arg_bio_phm[,1:3])
arg_bio_phm$Ss <- rowMeans(arg_bio_phm[,4:6])
arg_bio_phm$Wt <- rowMeans(arg_bio_phm[,7:9])
arg_bio_phm$Sn <- rowMeans(arg_bio_phm[,10:12])
arg_bio_phm <- arg_bio_phm[c(2,5,8,11,12,17),]
arg_bio_phm <- arg_bio_phm[,c(7:9,15,4:6,14,1:3,13,10:12,16)]
colnames(arg_bio_phm) <- colnames(arg.m)
rbind(arg.m,arg_bio_phm)

arg.all.phm <- pheatmap(rbind(arg.m,arg_bio_phm)[,c(4,8,12,16)],scale = "row",cluster_cols = F,color = colorRampPalette(c("navyblue","white","firebrick3"))(200),cellwidth = 25,cellheight = 25,treeheight_col = 10)