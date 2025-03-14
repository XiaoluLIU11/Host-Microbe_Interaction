rm(list = ls())
library(ggplot2)
library(WGCNA)
library(Seurat)
library(xlsx)
library(dplyr)
library(stringr)
library(data.table)
library(MetaboAnalystR)
library(RColorBrewer)


set.seed(12334)
name.rank <- c("Wt_1","Wt_2","Wt_3","Ss_1","Ss_2","Ss_3","Pg_1","Pg_2","Pg_3","Sn_1","Sn_2","Sn_3")
#mixomics wgcna analysis

#rna

readcount <- read.csv("D:/works/pg_omics/all_batch_readcounts/batch3/readcount.csv",header = T,row.names = 1)
readcount <- readcount[,c(7:9,4:6,1:3,10:12)]
readcount <- readcount[!grepl(row.names(readcount),pattern = "Novel*"),]
colnames(readcount) <- c("Wt_1","Wt_2","Wt_3","Ss_1","Ss_2","Ss_3","Pg_1","Pg_2","Pg_3","Sn_1","Sn_2","Sn_3")


#metabolome

meta.raw <- read.xlsx("D:/works/pg_omics/new_metabolome/metabolome/1.Data_Assess/all_group/ALL_sample_data.xlsx",sheetIndex = 1)
meta.count <- meta.raw[,c(1,12:23)]
row.names(meta.count) <- meta.count$Index
meta.count <- meta.count[,-1]

colnames(meta.count) <- c("Wt_1","Wt_2","Wt_3","Ss_1","Ss_2","Ss_3","Pg_1","Pg_2","Pg_3","Sn_1","Sn_2","Sn_3")

#raman
raman.m <- read.table("D:\\works\\pg_omics\\raman/alldata_pre.txt",header = T)
raman.meta <- read.table("D:\\works\\pg_omics\\raman/alldata_pre_meta.txt",header = T)
raman.meta$group_and_rep <- paste0(raman.meta$group_strain,"_",raman.meta$group_rep)

raman.data <- raman.meta[,c(1,11)] %>% dplyr::left_join(raman.m,.,by = "filename")
groups <- raman.meta$group_and_rep %>% unique()


mean.l <- list()
for (i in groups){
  raman.i <- raman.data[raman.data$group_and_rep == i,]
  mean.i <- apply(raman.i[,c(-1,-353)],2,mean)
  mean.l[[i]] <-  mean.i
}
do.call(rbind,mean.l) %>% data.frame() -> raman.m2

row.names(raman.m2) %>% 
  gsub("HGE_","Wt_",.) %>% 
  gsub("HGE-","",.) %>% 
  gsub("Pg-","",.) %>% 
  gsub("F2","",.)-> row.names(raman.m2)


raman.m3 <- raman.m2[match(name.rank,row.names(raman.m2)),]
raman.m4 <- apply(raman.m3, 2, as.numeric)
raman.ok <- raman.m4

maxs <- apply(raman.m4, 2, max)
mins <- apply(raman.m4, 2, min)


raman.ok <- scale(raman.m4,center = mins, scale = maxs - mins)
row.names(raman.ok) <- row.names(raman.m3)

raman.datrait <- t(raman.ok)
colnames(raman.datrait) <- name.rank

### raman marker

raman.seurat <- t(raman.m)
colnames(raman.seurat) <- raman.seurat[1,]
raman.seurat <- raman.seurat[-1,]
mode(raman.seurat) <- "numeric"

make_meta_table = function(x.n.mat.processed.1){
  meta_raw = colnames(raman.seurat)
  meta_raw = gsub('HGE_cell.*', 'Wt_', meta_raw)
  meta_raw = gsub('HGE-Pg_cell.*', 'Pg_', meta_raw)
  meta_raw = gsub('HGE-Pg-SnF2_cell.*', 'Sn_', meta_raw)
  meta_raw = gsub('HGE-Ss_cell.*', 'Ss_', meta_raw)
  matrix(unlist(str_split(meta_raw, "_")), 
         nrow=ncol(x.n.mat.processed.1), byrow=TRUE) -> mat
  mat = data.frame(mat)
  mat$ID = colnames(x.n.mat.processed.1)
  return(mat)
}

make_meta_table(raman.seurat) -> meta.table;meta.table

x <- list(spc=raman.seurat, meta = meta.table)


x$meta$id = names(x$spc)

dat = CreateSeuratObject(counts = x$spc)

dat[['drug']] = x$meta$X1

run_classic = function(dat.rep2){
  #dat.rep2 <- NormalizeData(object = dat.rep2)
  #dat.rep2 <- FindVariableFeatures(dat.rep2, selection.method = "vst", nfeatures = 400)
  # Identify the 10 most highly variable genes
  #top10 <- head(VariableFeatures(dat.rep2), 10)
  #VariableFeaturePlot(dat.rep2)
  all.genes <- rownames(dat.rep2)
  dat.rep2 <- ScaleData(dat.rep2, features = all.genes)
  # PCA
  dat.rep2 <- RunPCA(dat.rep2, features = all.genes)
  dat.rep2 <- RunUMAP(dat.rep2, dims = 1:10)
  dat.rep2 <- RunTSNE(dat.rep2, dims = 1:10,check_duplicates = FALSE)
  return(dat.rep2)
}

run_classic(dat) -> dat

Idents(dat) <- dat$drug

cluster2.markers <- FindAllMarkers(dat,min.pct = 0, logfc.threshold = 0,only.pos = T)

#pg_markers
pg_markers <- cluster2.markers[cluster2.markers$cluster == "Pg",]
pg_up_markers <- pg_markers %>% dplyr::filter(.,avg_log2FC > 0) %>% dplyr::filter(.,p_val_adj < 0.05)


sn_markers <- cluster2.markers[cluster2.markers$cluster == "Sn",]
sn_up_markers <- sn_markers %>% dplyr::filter(.,avg_log2FC > 0) %>% dplyr::filter(.,p_val_adj < 0.05)





##only consider Sn vs. Pg
dat.forsn <- subset(dat,drug == "Sn" | drug == "Pg")

cluster3.markers <- FindAllMarkers(dat.forsn,min.pct = 0, logfc.threshold = 0)

cluster3.markers %>% dplyr::filter(cluster == "Sn") %>% dplyr::filter(.,p_val_adj < 0.05) %>% dplyr::filter(avg_log2FC > 0) -> sn_up_info

sn_up_info$gene %>% gsub("X","",.) %>% gsub("[.].*","",.) -> sn_up_all


# using the key wavenumbers
c(1166.26,1398.21,1425.9,1563.03,1627.4) -> pg_find
c(1004.45,1029.82,1251.33,1304.08,1657.72) -> sn_find

#data preprocess
merge.expr <- rbind(readcount,meta.count)
m.mad <- apply(merge.expr, 1, mad)
dataExprVar <- merge.expr [which(m.mad > max(quantile(m.mad, probs=seq(0, 1, 0.25))[1],0)),];dim(dataExprVar)
dataExpr <- as.data.frame(t(dataExprVar))

gsg = goodSamplesGenes(dataExpr, verbose = 3)

if (!gsg$allOK){
  # Optionally, print the gene and sample names that were removed:
  if (sum(!gsg$goodGenes)>0)
    printFlush(paste("Removing genes:",
                     paste(names(dataExpr)[!gsg$goodGenes], collapse = ",")));
  if (sum(!gsg$goodSamples)>0)
    printFlush(paste("Removing samples:",
                     paste(rownames(dataExpr)[!gsg$goodSamples], collapse = ",")));
  # Remove the offending genes and samples from the data:
  dataExpr = dataExpr[gsg$goodSamples, gsg$goodGenes]
}

sampleTree = hclust(dist(dataExpr), method = "average")
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="")

#### raman intensity for dataTrait

row.names(raman.datrait) <- row.names(raman.datrait) %>% gsub("X","",.) 


care.raman_pg <- raman.datrait[row.names(raman.datrait) %in% pg_find,]
care.raman_sn <- raman.datrait[row.names(raman.datrait) %in% sn_find,]

datTraits <- t(rbind(care.raman_pg,care.raman_sn))

##
enableWGCNAThreads(nThreads = 2)
powers = c(c(1:20))
sft = pickSoftThreshold(dataExpr, powerVector=powers,
                        networkType="unsigned", verbose=5)

plot.new()
par(mfrow = c(1,2)) 
cex1 = 0.9 
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"))

text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");

abline(h=0.75,col="black",lty=4,lwd=1)


plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red") 


power = 11
adjacency = adjacency(dataExpr, power = power)
TOM = TOMsimilarity(adjacency)
dissTOM = 1-TOM


geneTree = hclust(as.dist(dissTOM), method = "average")

plot(geneTree, xlab="", sub="", main = "Gene clustering on TOM-based dissimilarity",
     labels = FALSE, hang = 0.04)

minModuleSize = 30      #模块基因数目
dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM,
                            deepSplit = 2, pamRespectsDendro = FALSE,
                            minClusterSize = minModuleSize);
table(dynamicMods)


dynamicColors = labels2colors(dynamicMods)
table(dynamicColors)
plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Gene dendrogram and module colors")


MEList = moduleEigengenes(dataExpr, colors = dynamicColors)
MEs = MEList$eigengenes
MEDiss = 1-cor(MEs);
METree = hclust(as.dist(MEDiss), method = "average")
plot(METree, main = "Clustering of module eigengenes",xlab = "", sub = "") 


MEDissThres = 0.3 #剪切高度可修改
abline(h=MEDissThres, col = "red")


merge = mergeCloseModules(dataExpr, dynamicColors, cutHeight = MEDissThres, verbose = 3)
mergedColors = merge$colors
mergedMEs = merge$newMEs

plotDendroAndColors(geneTree, cbind(dynamicColors,mergedColors),c("Dynamic Tree Cut ","mergedColors"),
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,cex.colorLabels = 1,
                    main = "Gene dendrogram and module colors")
moduleColors = mergedColors
table(moduleColors)
table(moduleColors) %>% data.frame() -> pie_chart_df

# pie plot
pie_chart_df$moduleColors <- reorder(pie_chart_df$moduleColors, -pie_chart_df$Freq)

pie.p <- ggplot(pie_chart_df, aes(x = "", y = Freq, fill = moduleColors)) +
  geom_bar(stat = "identity", width = 1) +
  coord_polar(theta = "y") +
  theme_void() +
  scale_fill_manual(values = pie_chart_df$moduleColors)

colorOrder = c("grey", standardColors(50))
moduleLabels = match(moduleColors, colorOrder)-1
MEs = mergedMEs

plotEigengeneNetworks(MEs,"Eigengene adjacency heatmap", marDendro = c(3,3,2,4), marHeatmap = c(3,4,2,2), plotDendrograms = T,)


nGenes = ncol(dataExpr)
nSamples = nrow(dataExpr)
moduleTraitCor = cor(MEs, datTraits, use = "p")
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)

textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                   signif(moduleTraitPvalue, 1), ")", sep = "")

par(mar = c(6, 8.5, 3, 3));


labeled.Heatmap <- labeledHeatmap(Matrix = moduleTraitCor,
                                  xLabels = colnames(datTraits),
                                  yLabels = names(MEs),
                                  ySymbols = names(MEs),
                                  colorLabels = FALSE,
                                  colors = blueWhiteRed(50),
                                  textMatrix = textMatrix,
                                  setStdMargins = FALSE,
                                  cex.text = 1,
                                  zlim = c(-1,1),
                                  main = paste("Module-trait relationships"))


module_cor_m <- list(moduleTraitCor,moduleTraitPvalue)

bk <- c(seq(-1,-0.01,by=0.01),seq(0,1,by=0.01))

module_cor_m_c <- module_cor_m[[1]]
module_cor_m_c[module_cor_m_c < 0.5 & module_cor_m_c > -0.5] <- module_cor_m_c[module_cor_m_c < 0.5 & module_cor_m_c > -0.5]/2
pheatmap(module_cor_m_c,cluster_rows = T,cluster_cols = F,border_color = 'black',cellwidth = 20,cellheight = 20,
         col =  rev(colorRampPalette(brewer.pal(9, "RdBu")) (240)))

### enrichment analysis
gene_anno <- fread("D:\\works\\pg_omics\\all_batch_readcounts\\batch3/fpkm_genename.xls")
probes = names(dataExpr)
probes2gene_anno = match(probes,gene_anno$geneID)

probes2meta_anno = match(probes,meta.raw$Index)

geneInfo0 = data.frame(gene_id = colnames(dataExpr),
                       geneanno = gene_anno$`Gene name`[probes2gene_anno],
                       metaanno = meta.raw$HMDB[probes2meta_anno],
                       moduleColor = moduleColors)


module_names <- c("white","orangered3","chocolate3")
for (m in module_names){
  module_members <- geneInfo0 %>% dplyr::filter(moduleColor == m)
  g_list <- module_members$geneanno %>% na.omit()
  g_list <- g_list[g_list != "-"]
  
  ### candidate modules kegg enrichment for joint-analysis
  mSet<-InitDataObjects("conc", "pathinteg", FALSE)
  mSet<-SetOrganism(mSet, "hsa")
  gene_list <- paste(g_list, collapse = "\r")
  mSet<-PerformGeneMapping(mSet, gene_list, "hsa", "symbol");
  
  m_list <- module_members$metaanno %>% na.omit()
  m_list <- m_list[m_list != "-"]  
  metabolite_list <- paste(m_list, collapse = "\r")
  mSet<-PerformCmpdMapping(mSet, metabolite_list, "hsa", "hmdb");
  mSet<-CreateMappingResultTable(mSet)
  
  mSet<-PrepareIntegData(mSet)
  mSet<-PerformIntegPathwayAnalysis(mSet, "dc", "hyper", "all", "query")
  mSet<-CreateIntegMatchingTable(mSet);
  #mSet<-SaveTransformedData(mSet)
  
  mSet$dataSet$path.m %>% data.frame() -> intergration_analysis_res
  #write.csv(intergration_analysis_res,paste0("D:/works/pg_omics/res/new_res_",m,"_kegg_joint_res.csv"))
  
}


wgcna_kegg_res <- read.table("D:/works/pg_omics/res/dems_analysis2/joint-kegg_enrichment.txt",sep = "\t", header = T)
kegg.plot <- ggplot(data = wgcna_kegg_res,aes(x = reorder(Pathway,Rank,decreasing = T), y = Hits,color = FDR))+ 
  geom_point(aes(size = GeneRatio))+
  scale_color_gradient(low="#E3170D",high = "#1E90FF")+
  #scale_shape_manual(values = c(21,16))+
  #scale_size_area(limits = c(1,15))+
  theme_classic()+
  coord_flip()+
  labs(x=NULL,y = " Member counts")+
  scale_x_discrete( labels=function(x) str_wrap(x, width=20))+ #position = "top",
  theme(axis.line = element_line(colour = "black"), 
        axis.text = element_text(color = "black",size = 9), #,angle = 90, hjust = 1
        legend.text = element_text(size = 9),
        legend.title=element_text(size=9),
        panel.grid.major.y = element_blank(),
        axis.title.x = element_text(size = 11),
        axis.title.y = element_text(size = 11),
        strip.text.x = element_text(size=8, color = "white"),
        strip.background = element_rect(fill="black"))+
  ggtitle("Candidate Modules Joint-KEGG enrichment")

### figE
nico.df <- read.table("D:/works/pg_omics/res/Nicotinamide_heatmap.txt",sep = "\t",header = T,row.names = 1)
nico.rna.pheatmap <- pheatmap(nico.df[,-1:-2],scale = "row",cluster_rows = F,cluster_cols = F,
                              cellwidth = 10,cellheight = 10,gaps_row = seq(1:10),border_color = "black",
                              color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(100),
                              filename = "D:/works/pg_omics/res/Nicotinamide_heatmap.png")


