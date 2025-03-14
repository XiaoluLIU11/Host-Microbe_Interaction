library("RamanEx")
library(ggbiplot)
library("MASS")
library("lattice")
library('optparse')
library("pROC")
library("R.utils")
library(xlsx)
library("scales")
library("pls")
func_lib()

rootpath <- "/Users/chenrongze/Myprogrames/2_PFSup/郑晓姗/Ramanome (P&G)/all_data"
outpath <- func_initialization(rootpath, output_folder = "output_0115")

##################### 1. get spc, meta, info files
folderpath <- rootpath
Name_group <- c("group_strain", "B", "group_cell","D","group_rep")
read_allsinspc(folderpath, Name_group, outpath,instrument ="CAST-R")

##################### 2. read and reorganize
filepath <- paste(folderpath,"/output_1101/files",sep = "")
all_data_list <- read_3data(filepath)

spc_reor_filt <- dcast(all_data_list[[1]],filename ~ wavenumber)
meta_filt <- filter(all_data_list[[2]],filename %in% spc_reor_filt$filename)
data_info <- all_data_list[[3]]

for(i in levels(factor(meta_filt$max_wave)))
{
  a <- filter(meta_filt,max_wave == i)$Group
  a <- a[!duplicated(a)]
  print(a)
}
spc_reor_filt_1 <- filter(spc_reor_filt,filename %in% filter(meta_filt,max_wave == "3200.75")$filename)
spc_reor_filt_2 <- filter(spc_reor_filt,filename %in% filter(meta_filt,max_wave == "3200.87")$filename)
spc_reor_filt_3 <- filter(spc_reor_filt,filename %in% filter(meta_filt,max_wave == "3201.15")$filename)
spc_reor_filt_1 <- as.data.frame(spc_reor_filt_1)
spc_reor_filt_2 <- as.data.frame(spc_reor_filt_2)
spc_reor_filt_3 <- as.data.frame(spc_reor_filt_3)

spc_reor_filt_1 <- spc_reor_filt_1 %>% dplyr::select(colnames(spc_reor_filt_1[1,!is.na(spc_reor_filt_1[1,])]))
spc_reor_filt_2 <- spc_reor_filt_2 %>% dplyr::select(colnames(spc_reor_filt_2[1,!is.na(spc_reor_filt_2[1,])]))
spc_reor_filt_3 <- spc_reor_filt_3 %>% dplyr::select(colnames(spc_reor_filt_3[1,!is.na(spc_reor_filt_3[1,])]))

colnames(spc_reor_filt_1) <- colnames(spc_reor_filt_2) <- colnames(spc_reor_filt_3)
spc_reor_filt_all <- rbind(spc_reor_filt_1,spc_reor_filt_2,spc_reor_filt_3)

##################### 3. data clean, QC, pretreatment
hs_baseline <- pre_baseline(spc_reor_filt_all,1)
hs_pre <- pre_norm(hs_baseline,method = "max")
hs_pre_data <- data.table(hs_pre)
meta_data <- as.data.frame(meta_filt)
fwrite(hs_pre_data, paste(filepath, "/alldata_pre", ".txt", sep = ""), 
       row.names = F, col.names = T, quote = F,sep = "\t")
fwrite(meta_data, paste(filepath, "/alldata_pre_meta", ".txt", sep = ""), 
       row.names = F, col.names = T, quote = F,sep = "\t")

##################### 5. 全谱分析、鉴定及差异判断RICH
hs_pre <- read.table(paste(filepath,"/","alldata_pre.txt", sep=""),header = T,sep="\t",fileEncoding = "utf-8")
colnames(hs_pre) <- gsub("X","",colnames(hs_pre))
hs_meta <- fread(paste(filepath, "/alldata_pre_meta.txt", sep = ""), header = TRUE, sep = "\t")

spc_melt_summary <- plot_func_spcgroup(hs_pre,hs_meta,labels = c("group_strain")) 
melt_summary <- spc_melt_summary %>% filter(wavenumber > 600,wavenumber < 3200)
melt_summary$group_strain <- factor(melt_summary$group_strain,levels = c("HGE-Pg","HGE-Pg-SnF2","HGE-Ss","HGE"))

levels <- levels(factor(melt_summary$group_strain))
n <- length(levels) * 0.1
for(i in levels)
{
  inds <- which(melt_summary$group_strain == i )
  melt_summary[inds,]$value <- melt_summary[inds,]$value + n
  n <- n - 0.1
}
fwrite(melt_summary, paste(filepath, "/Fig2_A ", ".txt", sep = ""), 
       row.names = F, col.names = T, quote = F,sep = "\t")

melt_summary_pg <- melt_summary %>% filter(group_strain == "HGE-Pg") %>% dplyr::select(c("wavenumber","value"))
fwrite(melt_summary_pg, paste(filepath, "/Fig2_A-pg ", ".txt", sep = ""), 
       row.names = F, col.names = T, quote = F,sep = "\t")
#### Fig2.B 差谱600~1800
spc_melt_summary <- func_diff(hs_pre,hs_meta)
spc_melt_summary$group_strain[which(spc_melt_summary$group_strain == "HGE")] <- "Control"
spc_melt_summary$group_strain <- factor(spc_melt_summary$group_strain,levels = c("HGE-Pg","HGE-Pg-SnF2","HGE-Ss","Control"))
levels <- levels(factor(spc_melt_summary$group_strain))

n <- length(levels)* 0.004
for(i in levels)
{
  n <- n - 0.004
  inds <- which(spc_melt_summary$group_strain == i )
  spc_melt_summary[inds,]$value <- spc_melt_summary[inds,]$value + n
}
fwrite(spc_melt_summary, paste(filepath, "/Fig2_B.txt", sep = ""), 
       row.names = F, col.names = T, quote = F,sep = "\t")

spc_melt_summary_pg <- spc_melt_summary %>% filter(group_strain == "HGE-Pg") %>% dplyr::select(c("wavenumber","value"))
fwrite(spc_melt_summary_pg, paste(filepath, "/Fig2_B-pg ", ".txt", sep = ""), 
       row.names = F, col.names = T, quote = F,sep = "\t")
#### Fig2.C 差谱600~1800
hs_ori <- spc_waveumber(spc_reor_filt_all,600,1800)
hs_baseline <- pre_baseline(hs_ori,1)
hs_pre <- pre_norm(hs_baseline,method = "sum")
hs_meta <- meta_filt

RICH_data <- data.table(left_join(select(meta_data,c("group_strain","filename")),hs_pre,by= "filename"))
RICH_data$RICH_group <- RICH_data$group_strain
RICH_data$RICH_group[which(RICH_data$RICH_group == "HGE")] <- "Control"
#RICH_data <- RICH_data[-c(2,68,138,204),]

RICH_data <- melt(RICH_data[,-1],id.vars = c("filename","RICH_group"),variable.name = "wavenumber",value.name = "value")
RICH_data$wavenumber <- as.numeric(as.character(RICH_data$wavenumber))

RICH_peaks_result <- cal_RICH_tt(RICH_data,Pvalue = 0.0001)
RICH_peaks_i <- RICH_peaks_result[[1]]
p_value_Data <- RICH_peaks_result[[2]]

###Fig2.C heatmap
wavenumber <- as.numeric(as.character(levels(factor(RICH_data$wavenumber))))

RICH_data$wavenumber <- ceiling(as.numeric(as.character(RICH_data$wavenumber)))
RICH_data$value <- RICH_data$value + min(RICH_data$value)
RICH_data$value <- RICH_data$value/max(RICH_data$value)

wavenumber_ins <- c()
for(i in seq(600,1700,100))
{
  wavenumber_ins_i <- RICH_data$wavenumber[which.min(abs(RICH_data$wavenumber-i))]
  RICH_data$wavenumber[which(RICH_data$wavenumber == wavenumber_ins_i)] <- i
}
fwrite(RICH_data, paste(filepath, "/Fig2_C.txt", sep = ""), 
       row.names = F, col.names = T, quote = F,sep = "\t")


plot_tile <- ggplot(all_data_table, aes(x=factor(wavenumber), y=1, fill=p_value)) + 
  geom_tile() + theme_bw() + coord_equal(ratio = 1) +
  scale_x_discrete(breaks = seq(600,1700,100)) +
  scale_fill_gradientn(colours =myPal) +
  theme(panel.grid = element_blank(),
        panel.border = element_blank(),
        legend.title = element_blank(),
        legend.key.width = unit(1,"pt"),
        legend.key.height = unit(10, "pt"),
        axis.ticks.y = element_blank(),
        axis.title = element_blank(),
        axis.text.x = element_text(size = 10,angle = 45,hjust = 1,vjust =1), 
        axis.text.y = element_blank(),        
        axis.text = element_text())
plot_tile

got_palette <- c("#1A5878", "#C44237", "#AD8941", "#E99093","#50594B", 
                 "#8968CD", "#9ACD32", "#9e4f49", "#adbd95","#c7dbd5")
plot <- ggplot(p_value_Data_i,aes(x=wavenumber,size = abs(p_value)))+
  geom_point(alpha=0.7) + scale_size(range=c(1,5)) + #facet_grid(.~CD_label) + 
  ylab("Assignment")+
  labs(color = "Communty",size = "Number of peaks")+
  scale_fill_manual(na.value = "white",values = got_palette,aesthetics = c("fill","color")) +
  scale_y_discrete(position = "right") +
  scale_x_discrete(position = "top") +
  default_theme()+
  theme(
    legend.position = "right",
    legend.direction = "vertical",
    strip.text = element_blank(),
    axis.title.x = element_blank(),
    axis.text.x.top = element_text(size = 15,angle = 45,
                                   hjust = 1,vjust = 1,family = "Arial"))
plot

plot_pvalue_2 <- ggplot(data = p_value_Data_i,group = 1) +
  geom_line(aes(x = wavenumber, y = diff),size = 0.8) +
  theme_bw() + ylab("Normalized Intensity (a.u.)") +
  xlab(expression(paste('Wavenumber (cm'^{-1},')',sep = ""))) +
  geom_vline(xintercept = RICH_peaks_i, colour = "red",alpha = 0.2,size = 0.5)+
  scale_x_continuous(breaks = seq(600,1800,100)) +
  default_theme()+
  theme(
    axis.text.x = element_text(size = 15,angle = 45,hjust = 1,vjust = 1),
    axis.ticks.y = element_blank(),
    axis.text.y = element_blank())

plot_pvalue <- plot_pvalue_1 / plot_pvalue_2
plot_pvalue
ggsave(filename=paste(outpath,"/Fig3_A.png", sep=""),plot=plot_pvalue,
       limitsize=T,width = 8,height = 8)#,


folderpath <- "/Users/chenrongze/Myprogrames/2_PFSup/郑晓姗/Ramanome (P&G)/all_data"
filepath <- paste(folderpath,"/output_1101/files",sep = "")

hs_pre <- fread(paste(filepath, "/alldata_pre", ".txt", sep = ""),header = TRUE,sep = "\t")
hs_meta <- fread(paste(filepath, "/alldata_pre_meta", ".txt", sep = ""),header = TRUE,sep = "\t")

RICH_data <- data.table(left_join(dplyr::select(hs_meta,c("group_strain","filename")),hs_pre,by= "filename"))
RICH_data$RICH_group <- RICH_data$group_strain
RICH_data$RICH_group[which(RICH_data$RICH_group == "HGE")] <- "Control"

RICH_data <- melt(RICH_data[,-1],id.vars = c("filename","RICH_group"),variable.name = "wavenumber",value.name = "value")
RICH_data$wavenumber <- as.numeric(as.character(RICH_data$wavenumber))

RICH_peaks_result <- cal_RICH_tt(RICH_data,Pvalue = 0.0001)
RICH_peaks_i <- RICH_peaks_result[[1]]
p_value_Data <- RICH_peaks_result[[2]]
p_value_Data$RICH_group <- factor(p_value_Data$RICH_group,levels = c("HGE-Pg","HGE-Pg-SnF2","HGE-Ss","Control"))
RICH_peaks_i <- sort(as.numeric(RICH_peaks_i[!duplicated(RICH_peaks_i)]))

control_mean <- RICH_data %>%
  filter(RICH_group == "Control") %>% 
  dplyr::group_by( wavenumber) %>% 
  dplyr::summarise(Value = mean(value))

melt_summary_a <- RICH_data %>%
  dplyr::group_by(filename)  %>%
  dplyr::mutate(diff = value - control_mean$Value)

RICH_slect_spc <- filter(melt_summary_a[,-4],wavenumber %in% RICH_peaks_i)
RICH_slect_spc <- dcast(RICH_slect_spc,filename+RICH_group ~ wavenumber)

weights_diff <- p_value_Data %>% filter(RICH_group == "HGE-Pg",wavenumber %in% RICH_peaks_i) %>% dplyr::select(diff) #\- 0.008
weights_diff <- weights_diff /sum(weights_diff^2)
fwrite(weights_diff, paste(filepath, "/weights_diff.txt", sep = ""), 
       row.names = F, col.names = T, quote = F,sep = "\t")
fwrite(as.data.table(RICH_peaks_i), paste(filepath, "/RICH_peaks_i.txt", sep = ""), 
       row.names = F, col.names = T, quote = F,sep = "\t")

RICH_value <- as.matrix(RICH_slect_spc[,3:ncol(RICH_slect_spc)]) %*% as.matrix(weights_diff)
RICH_value <- data.table(RICH_group = RICH_slect_spc$RICH_group,RICH_value)
#RICH_value <- RICH_value %>% filter(RICH_group != "HGE-Ss")
color_cluster <- c("#BDBDBDFF","#A50026","#4575B4","#42B540FF")
RICH_value$RICH_group <- factor(RICH_value$RICH_group,levels = c("Control","HGE-Pg","HGE-Pg-SnF2","HGE-Ss"))
fwrite(RICH_value, paste(filepath, "/Fig3_B-1.txt", sep = ""), 
       row.names = F, col.names = T, quote = F,sep = "\t")

RICH_value <- fread(paste(filepath, "/Fig3_B-1", ".txt", sep = ""),header = TRUE,sep = "\t")


levels_RICH_group <- levels(factor(RICH_value$RICH_group))
p_value <- matrix(0,nrow = length(levels_RICH_group),ncol = length(levels_RICH_group))
for( i in 1:length(levels_RICH_group)){
  for( j in 1:length(levels_RICH_group)){
    if (levels_RICH_group[i] == levels_RICH_group[j]){
      p_value[i,j] <- 0
    }else{
      RICH_value_i <- RICH_value %>% filter(RICH_group %in% c(levels_RICH_group[i],levels_RICH_group[j]))
      p_value[i,j] <- kruskal.test(diff~RICH_group,RICH_value_i)$p.value
    }
  }
}
p_value <- data.frame(p_value)
colnames(p_value) <- rownames(p_value) <- levels_RICH_group
fwrite(p_value, paste(filepath, "/Fig3_B-1-pvalue.txt", sep = ""), 
       row.names = T, col.names = T, quote = F,sep = "\t")

plot<-ggplot(RICH_value,aes(x=RICH_group,y=diff,group= RICH_group,fill = RICH_group))+
  theme_bw() + geom_point(aes(group=RICH_group,color = RICH_group),position=position_jitter(0.1))+
  geom_boxplot(aes(group=RICH_group),width = 0.4,alpha=0.2,outlier.alpha = 0)+ 
  scale_fill_manual(values = color_cluster,aesthetics = c("fill","color"))+
  geom_hline(aes(yintercept=0), colour="red", linetype="dashed",size = 1)+
  geom_hline(aes(yintercept=1), colour="red", linetype="dashed",size = 1)+
  theme(
    panel.grid = element_blank(),
    panel.grid.minor = element_blank(),
    legend.title = element_markdown(size = 15),
    legend.text = element_markdown(size = 15),
    legend.position = "none",
    legend.background = element_blank(),
    text = element_text(color="black"),
    axis.title.y = element_text(size = 20, angle = 90),
    axis.text.x = element_text(size = 15,angle = 45,vjust = 1,hjust = 1),
    axis.text.y = element_text(size = 15),
    strip.text = element_text(size = 20),
    axis.ticks = element_line(size = 1),
    axis.ticks.length = unit(0.4,"lines"),
    axis.title= element_text(size=20),
    axis.title.x= element_blank())+
  ylab('Stress Index')
plot 
ggsave(filename=paste(outpath,"//","Fig3_B-1.svg", sep=""),plot=plot,width = 8,height = 8)#,

 
levels_RICH_group <- levels(factor(RICH_value$RICH_group))
Number_batch <- list(1:22,23:44,45:66)
RICH_value_HI <- matrix(0,nrow = 3,ncol = length(levels_RICH_group))
for( i in 1:length(levels_RICH_group)){
  for( j in 1:3){
      RICH_value_i <- filter(RICH_value,RICH_group == levels_RICH_group[i])[Number_batch[[j]],]
      RICH_value_HI[j,i] <- sd(RICH_value_i$diff)#/mean(RICH_value_i$diff)
  }
}
RICH_value_HI <- data.frame(RICH_value_HI)
colnames(RICH_value_HI) <- levels_RICH_group
RICH_value_HI <- melt(RICH_value_HI)
RICH_value_HI <- Rmisc::summarySE(RICH_value_HI, measurevar = "value", groupvars = c("variable"))
plot<-ggplot(RICH_value_HI,aes(x=variable,y = value,fill = variable))+
  geom_col() +   geom_errorbar(aes(ymin = value - sd, ymax = value + sd),width = 0.2) +
  theme_bw() + 
  scale_fill_manual(values = color_cluster,aesthetics = c("fill","color"))+
  theme(
    panel.grid = element_blank(),
    panel.grid.minor = element_blank(),
    legend.title = element_markdown(size = 15),
    legend.text = element_markdown(size = 15),
    legend.position = "none",
    legend.background = element_blank(),
    text = element_text(color="black"),
    axis.title.y = element_text(size = 20, angle = 90),
    axis.text.x = element_text(size = 15,angle = 45,vjust = 1,hjust = 1),
    axis.text.y = element_text(size = 15),
    strip.text = element_text(size = 20),
    axis.ticks = element_line(size = 1),
    axis.ticks.length = unit(0.4,"lines"),
    axis.title= element_text(size=20),
    axis.title.x= element_blank())+
  ylab('Standard deviation')
plot 
ggsave(filename=paste(outpath,"//","Fig4-si-Standard deviation.pdf", sep=""),plot=plot,width = 6,height = 6)#,



RICH_value_summary <- Rmisc::summarySE(RICH_value, measurevar = "diff", groupvars = c("RICH_group"))
RICH_value_summary_melt <- melt(RICH_value_summary[,c(1,3,4)], id.vars = c("RICH_group"),
                                variable.name = "label", value.name = "value")
RICH_value_summary_melt$label <- gsub("diff","Mean",RICH_value_summary_melt$label)
RICH_value_summary_melt$label <- gsub("sd","Sd",RICH_value_summary_melt$label)
color_cluster <- c("#BDBDBDFF","#A50026","#4575B4","#42B540FF")
RICH_value_summary_melt$RICH_group <- factor(RICH_value_summary_melt$RICH_group,levels = c("Control","HGE-Pg","HGE-Pg-SnF2","HGE-Ss"))

plot<-ggplot(RICH_value_summary_melt,aes(x=RICH_group,y = value,fill = RICH_group))+
  geom_col() +   facet_grid(. ~ label) + 
  theme_bw() + 
  scale_fill_manual(values = color_cluster,aesthetics = c("fill","color"))+
  theme(
    panel.grid = element_blank(),
    panel.grid.minor = element_blank(),
    legend.title = element_markdown(size = 15),
    legend.text = element_markdown(size = 15),
    legend.position = "none",
    legend.background = element_blank(),
    text = element_text(color="black"),
    axis.title.y = element_text(size = 20, angle = 90),
    axis.text.x = element_text(size = 15,angle = 45,vjust = 1,hjust = 1),
    axis.text.y = element_text(size = 15),
    strip.text = element_text(size = 20),
    axis.ticks = element_line(size = 1),
    axis.ticks.length = unit(0.4,"lines"),
    axis.title= element_text(size=20),
    axis.title.x= element_blank())+
  ylab('Stress Index')
plot 
ggsave(filename=paste(outpath,"//","Fig3_C-1.svg", sep=""),plot=plot,width = 8,height = 8)#,


plot<-ggplot(RICH_value,aes(x=diff,group= RICH_group,fill = RICH_group))+
  geom_density(aes(y = ..ndensity..,color = RICH_group),position = "identity",alpha = 0.9) + 
  geom_histogram(aes(y = ..ncount..),position = "identity",bins = 100,alpha = 0.7) + 
  scale_fill_manual(values = color_cluster,aesthetics = c("fill","color"))+
  facet_grid(. ~ RICH_group) +  theme_bw() + 
  geom_vline(aes(xintercept=0), colour="red", linetype="dashed",size = 1)+
  geom_vline(aes(xintercept=1), colour="red", linetype="dashed",size = 1)+
  theme(
    panel.grid = element_blank(),
    panel.grid.minor = element_blank(),
    legend.title = element_markdown(size = 15),
    legend.text = element_markdown(size = 15),
    legend.position = "none",
    legend.background = element_blank(),
    text = element_text(color="black"),
    axis.title.y = element_text(size = 20, angle = 90),
    axis.text.x = element_text(size = 15,angle = 45,vjust = 1,hjust = 1),
    axis.text.y = element_text(size = 15),
    strip.text = element_text(size = 20),
    axis.ticks = element_line(size = 1),
    axis.ticks.length = unit(0.4,"lines"),
    axis.title= element_text(size=20),
    axis.title.x= element_blank())+
  ylab('Stress Index')
plot 
ggsave(filename=paste(outpath,"//","Fig3_D-1.svg", sep=""),plot=plot,width = 20,height = 6)#,

##########Pg Stress Index
RICH_data <- data.table(left_join(select(hs_meta,c("group_strain","filename")),hs_pre,by= "filename"))
RICH_data$RICH_group <- RICH_data$group_strain
RICH_data$RICH_group[which(RICH_data$RICH_group == "HGE")] <- "Control"
#RICH_data <- RICH_data[-c(2,68,138,204),]

RICH_data <- melt(RICH_data[,-1],id.vars = c("filename","RICH_group"),variable.name = "wavenumber",value.name = "value")
RICH_data$wavenumber <- as.numeric(as.character(RICH_data$wavenumber))

RICH_peaks_result <- cal_RICH_tt(RICH_data,Pvalue = 0.0001)
RICH_peaks_i <- RICH_peaks_result[[1]]
p_value_Data <- RICH_peaks_result[[2]]
p_value_Data$RICH_group <- factor(p_value_Data$RICH_group,levels = c("HGE-Pg","HGE-Pg-SnF2","HGE-Ss","Control"))

n <- length(levels(factor(p_value_Data$RICH_group))) * 0.004
for(i in levels(factor(p_value_Data$RICH_group)))
{
  n <- n - 0.004
  inds <- which(p_value_Data$RICH_group == i )
  p_value_Data[inds,]$diff <- p_value_Data[inds,]$diff + n
}
Color_table_1 <- c("#F8766D","#7CAE00","#00BFC4","gray")
RICH_peaks_i <- sort(as.numeric(RICH_peaks_i[!duplicated(RICH_peaks_i)]))
plot_pvalue <- ggplot(data = p_value_Data,aes(x = wavenumber, y = diff,fill = RICH_group),group = 1) +
  geom_ribbon(aes(ymin = diff - sd, ymax = diff + sd),alpha = 0.3,) +
  geom_line(aes(color = RICH_group),size = 0.8) + 
  theme_bw() + ylab("Normalized Intensity (a.u.)") + 
  xlab(expression(paste('Wavenumber (cm'^{-1},')',sep = ""))) +
  scale_fill_manual(values = Color_table_1,aesthetics = c("fill","color"))+
  geom_hline(yintercept = 0.004,linetype = "dashed", colour = "gray",size = 0.8)+
  geom_hline(yintercept = 0,linetype = "dashed", colour = "gray",size = 0.8)+
  geom_hline(yintercept = 0.008,linetype = "dashed", colour = "gray",size = 0.8)+
  geom_hline(yintercept = 0.012,linetype = "dashed", colour = "gray",size = 0.8)+
  geom_vline(xintercept = RICH_peaks_i, colour = "red",alpha = 0.2,size = 0.25)+
  theme(
    legend.title = element_blank(),
    legend.text = element_markdown(size = 15),
    #legend.position = "none",
    legend.background = element_blank(),
    text = element_text(color = "black"),
    axis.title.y = element_text(size = 15, angle = 90),
    axis.text.x = element_text(size = 10,angle = 45,hjust = 1,vjust = 1),
    axis.text.y = element_blank(),
    strip.background = element_rect(fill = "white"),
    strip.text = element_text(size = 20),
    axis.ticks = element_line(size = 1),
    axis.ticks.y = element_blank(),
    axis.ticks.length = unit(0.4,"lines"),
    axis.title = element_text(size = 20))
plot_pvalue
ggsave(filename=paste(outpath,"//","all_group-diff.png", sep=""),plot=plot_pvalue)#,

a <- 0.01
n <- length(levels(factor(melt_summary_a$RICH_group))) * a
for(i in levels(factor(melt_summary_a$RICH_group)))
{
  n <- n - a
  inds <- which(melt_summary_a$RICH_group == i )
  melt_summary_a[inds,]$diff <- melt_summary_a[inds,]$diff + n
}

Color_table_2 <- c("gray","#F8766D","#7CAE00","#00BFC4")
plot_pvalue <- ggplot(data = melt_summary_a,aes(x = wavenumber, y = diff,group = filename)) +
  geom_line(aes(color = RICH_group),size = 0.8) + 
  theme_bw() + ylab("Normalized Intensity (a.u.)") + 
  xlab(expression(paste('Wavenumber (cm'^{-1},')',sep = ""))) +
  scale_fill_manual(values = Color_table_2,aesthetics = c("fill","color"))+
  theme(
    legend.title = element_blank(),
    legend.text = element_markdown(size = 15),
    #legend.position = "none",
    legend.background = element_blank(),
    text = element_text(color = "black"),
    axis.title.y = element_text(size = 15, angle = 90),
    axis.text.x = element_text(size = 10,angle = 45,hjust = 1,vjust = 1),
    axis.text.y = element_blank(),
    strip.background = element_rect(fill = "white"),
    strip.text = element_text(size = 20),
    axis.ticks = element_line(size = 1),
    axis.ticks.y = element_blank(),
    axis.ticks.length = unit(0.4,"lines"),
    axis.title = element_text(size = 20))
plot_pvalue
ggsave(filename=paste(outpath,"//","all_file-diff.png", sep=""),plot=plot_pvalue)#,

control_mean <- RICH_data[-c(2,68,138,204),] %>%
  filter(RICH_group == "Control") %>% 
  summarySE(measurevar = "value",groupvars = c("RICH_group", "wavenumber")) %>% 
  select(wavenumber,value)

melt_summary_a <- RICH_data %>%
  group_by(filename)  %>%
  mutate(diff = value - control_mean$value)

RICH_slect_spc <- filter(melt_summary_a[,-4],wavenumber %in% RICH_peaks_i)
RICH_slect_spc <- dcast(RICH_slect_spc,filename+RICH_group ~ wavenumber)
RICH_slect_spc <- RICH_slect_spc[-c(2,68,138,204),]

weights_diff <- p_value_Data %>% filter(RICH_group == "HGE-Pg",wavenumber %in% RICH_peaks_i) %>% select(diff) #- 0.008
weights_diff <- weights_diff /sum(weights_diff^2)

RICH_value <- as.matrix(RICH_slect_spc[,3:ncol(RICH_slect_spc)]) %*% as.matrix(weights_diff)
RICH_value <- data.table(RICH_group = RICH_slect_spc$RICH_group,RICH_value)
Color_table_2 <- c("gray","#F8766D","#7CAE00","#00BFC4")
plot<-ggplot(RICH_value,aes(x=RICH_group,y=diff,group= RICH_group,fill = RICH_group))+
  theme_bw() + geom_point(aes(group=RICH_group,color = RICH_group),position=position_jitter(0.1))+
  geom_boxplot(aes(group=RICH_group),alpha=0.2)+ 
  scale_fill_manual(values = Color_table_2,aesthetics = c("fill","color"))+
  geom_hline(aes(yintercept=0), colour="red", linetype="dashed",size = 1)+
  geom_hline(aes(yintercept=1), colour="red", linetype="dashed",size = 1)+
  theme(
    panel.grid = element_blank(),
    panel.grid.minor = element_blank(),
    legend.title = element_markdown(size = 15),
    legend.text = element_markdown(size = 15),
    legend.position = "none",
    legend.background = element_blank(),
    text = element_text(color="black"),
    axis.title.y = element_text(size = 20, angle = 90),
    axis.text.x = element_text(size = 15,angle = 45,vjust = 1,hjust = 1),
    axis.text.y = element_text(size = 15),
    strip.text = element_text(size = 20),
    axis.ticks = element_line(size = 1),
    axis.ticks.length = unit(0.4,"lines"),
    axis.title= element_text(size=20),
    axis.title.x= element_blank())+
  ylab('Pg Stress Index')
plot 
ggsave(filename=paste(outpath,"//","all_RICH_value.png", sep=""),plot=plot,width = 7,height = 7)#,

melt_summary_b <- dcast(melt_summary_a[,-4],filename+RICH_group ~ wavenumber)
svm_data <- org_model(melt_summary_b)
training <- svm_data[[1]][-c(1,63,119,180),]
test <- svm_data[[2]]


#svm：由于是分类问题，此处我们选择C-classification
spc_svm<- svm(x = training[,-(1:2)] , y = factor(training$RICH_group),type = 'C',kernel = 'radial')
svm.pred <- predict(spc_svm, test[,-(1:2)])

svm.pred_matrix <- as.data.frame(table(svm.pred,test$RICH_group))
svm.pred_matrix <- svm.pred_matrix %>%
  mutate(x = factor(Var2), # alphabetical order by default
         y = factor(svm.pred, levels = rev(unique(svm.pred))))

text_color <- svm.pred_matrix$Freq
text_color[which(text_color > 1/2 * max(text_color))] <- "white"
  text_color[which(text_color != "white")] <- "black"
    
  plot_tile <- ggplot(svm.pred_matrix, aes(x=x, y=y, fill=Freq)) + 
    geom_tile(color = "black") + theme_bw() + coord_equal() +
    scale_fill_distiller(palette="Blues", direction=1) +
    geom_text(aes(label=Freq), color=text_color) +# printing values
    theme(panel.grid = element_blank(),
          # panel.border = element_blank(),
          legend.title = element_blank(),
          axis.ticks = element_blank(),
          axis.title = element_blank(),
          axis.text.x = element_text(angle = 45,vjust = 1,hjust = 1),
          axis.text = element_text() #family = "Times"
    )
  plot_tile
  ggsave(filename=paste(outpath,"//","confusion matrix tile_final.png", sep=""),plot=plot_tile, limitsize=T,width = 5,height = 7)
  


Group <- meta_data$Group
tb_pre_mean <- cal_spc_mean(hs_pre_data,Group,save = F)
tb_pre_peak <- cal_findpeak(tb_pre_mean[,-c(2:3)],step = 10^-5)

tb_pre_mean <- tb_pre_mean[,-c(2:3)]
colnames(tb_pre_mean) <- c("Group",colnames(tb_pre_mean[,-1]))

colnames(hs_pre_data) <- c("Group",colnames(hs_pre_data[,-1]))
hs_pre_data <- hs_pre_data[,-344]

peak_all_i <- org_peak_all(tb_pre_peak)
spc_peak <- org_findpeak(hs_pre_data,peak_all_i)
spc_peak$Group <- meta_data$Group

plot_trend <- plot_peaktrend(spc_peak,meta_data)
plot_trend

plot_meanpeak <- plot_findpeak(tb_pre_mean,peak_all_i)
plot_meanpeak



hs_pre_1$group_B[which(hs_pre_1$group_B != "0")] <- "Grug"
levels(factor(hs_pre_1$group_B))

a <- hist(cluster.meansd.dataframe_melt$value,freq = T,
          breaks = seq(min(cluster.meansd.dataframe_melt$value),
                       max(cluster.meansd.dataframe_melt$value),length.out = 100))
a$breaks[which(find_peak_window(a$counts) == TRUE)] 

cluster.meansd.dataframe_melt <- spc_melt(hs_pre_1,"Group",from = 500,to = 1800)
wavenum_range <- levels(factor(cluster.meansd.dataframe_melt$wavenumber))
peak_bywindow <- c()
for(i in wavenum_range)
{
  melt_i <- filter(cluster.meansd.dataframe_melt,wavenumber == i)
  a <- hist(melt_i$value,freq = T,breaks = seq(min(melt_i$value),max(melt_i$value),length.out = 100))
  peak_bywindow <- cbind(peak_bywindow,find_peak_window(a$counts))
}
peak_bywindow <- data.frame(peak_bywindow)
colnames(peak_bywindow) <- wavenum_range
peak_nums <-c()
for(i in 1:ncol(peak_bywindow))
{
  peak_nums <- c(peak_nums,length(which(peak_bywindow[,i] == TRUE)))
}
peak_nums <- cbind(peak_nums,wavenum_range)
peak_nums <- data.frame(peak_nums)
peak_nums <- filter(peak_nums,peak_nums > 1)

for(i in peak_nums$wavenum_range)
{
  melt_i <- filter(cluster.meansd.dataframe_melt,wavenumber == i)
  plot<-ggplot(melt_i,aes(x=factor(wavenumber),y=value))+
    geom_point(aes(group=wavenumber),position=position_jitter(0.2))+
    geom_boxplot(aes(group=wavenumber),alpha=0.2,outlier.alpha = 0)+ 
    theme_bw()+ facet_wrap(Group~.,ncol = 8)+
    theme(
      legend.title = element_markdown(size = 15),
      legend.text = element_markdown(size = 15),
      legend.position = "none",
      legend.background = element_blank(),
      text = element_text(color="black"),
      axis.title.y = element_text(size = 20, angle = 90),
      axis.text.x = element_text(size = 15,angle = 45,vjust = 1,hjust = 1),
      axis.text.y = element_text(size = 15),
      strip.text = element_text(size = 20),
      strip.background = element_rect(fill = "white"),
      axis.ticks = element_line(size = 1),
      axis.ticks.length = unit(0.4,"lines"),
      axis.title= element_text(size=20))+
    xlab('Metabolic activity (CD-ratio)')
  plot
  ggsave(
    filename = paste(outpath,"/",i, "_bygroup.png", sep = ""),
    plot = plot,
    width = 20, height = 5, dpi = dpi
  )
}


melt_i <- filter(cluster.meansd.dataframe_melt,wavenumber == 992.419)
plot<-ggplot(melt_i,aes(x=factor(wavenumber),y=value))+
  geom_point(aes(group=wavenumber),position=position_jitter(0.2))+
  geom_boxplot(aes(group=wavenumber),alpha=0.2,outlier.alpha = 0)+ 
  theme_bw()+ facet_wrap(Group~.,ncol = 6)+
  theme(
    legend.title = element_markdown(size = 15),
    legend.text = element_markdown(size = 15),
    legend.position = "none",
    legend.background = element_blank(),
    text = element_text(color="black"),
    axis.title.y = element_text(size = 20, angle = 90),
    axis.text.x = element_text(size = 15,angle = 45,vjust = 1,hjust = 1),
    axis.text.y = element_text(size = 15),
    strip.text = element_text(size = 20),
    strip.background = element_rect(fill = "white"),
    axis.ticks = element_line(size = 1),
    axis.ticks.length = unit(0.4,"lines"),
    axis.title= element_text(size=20))+
  xlab('Metabolic activity (CD-ratio)')
plot

hist_CD_ratio <- ggplot(cluster.meansd.dataframe_melt,aes(x=value,fill=factor(wavenumber)))+
  geom_density(alpha = 0.4)+ theme_bw() + facet_wrap(.~wavenumber)+
  theme_bw()+
  theme(
    legend.title = element_markdown(size = 15),
    legend.text = element_markdown(size = 10),
    legend.position = "none",
    legend.background = element_blank(),
    text = element_text(color="black"),
    axis.title.y = element_text(size = 10, angle = 90),
    axis.text.x = element_text(size = 10,angle = 45,vjust = 1,hjust = 1),
    axis.text.y = element_text(size = 10),
    strip.text = element_text(size = 5),
    strip.background = element_rect(fill = "white"),
    axis.ticks = element_line(size = 1),
    axis.ticks.length = unit(0.4,"lines"),
    axis.title= element_text(size=10))+
  xlab('Metabolic activity (CD-ratio)')
hist_CD_ratio

ggsave(
  filename = paste(outpath, "plot_single_bygroup.png", sep = "/"),
  plot = plot_single_bygroup,
  width = 10, height = 10, dpi = dpi
)
##################### 3. mean spectra
meta_group <- c("group_concentration", "group_time")
spc_mean_df <- cal_mean(hs_pre, meta_group, outpath = outpath)

##################### 4. plot mean spectra
plot_mean_concentration <- plot_Meanspc_bygroup(hs_pre[,,900~1700],0.002)
plot_mean_concentration
ggsave(
  filename = paste(outpath, "plot_mean_concentration.pdf", sep = "/"),
  plot = plot_hyperspec,
  width = 10, height = 5, dpi = dpi
)






###################### 2. pretreatment
hs_pre <- pre_allpipeline(hs_ori)

wavenum <- as.numeric(as.character(colnames(spc_mean_df[,-1])))
a <- spc_mean_df[1,-1]
hs_pre_pg <- new("hyperSpec",spc = a,wavelength = wavenum)

a_peak <- cal_findpeak(hs_pre_pg)
for(i in a_peak)
{
  a_peak[which(a_peak == i)] <- wavenum[which.min(abs(wavenum - i))]
}

a_peak_1<-a_peak
a_peak_2<-a_peak
a_peak_3<-a_peak
a_peak_4<-a_peak

all_peak <- c(a_peak_1,a_peak_2,a_peak_3,a_peak_4)
all_peak <- all_peak[which(!duplicated(all_peak))]
all_peak <- all_peak[which(all_peak < 1750 | all_peak > 2600)]
write.table(all_peak,
            file = paste(outpath, "/", "all_peaks.txt", sep = ""),
            row.names = FALSE,
            col.names = FALSE
)
###################### 3. mean spectra
meta_group <- "group_strain"
spc_mean_df <- cal_mean(hs_pre, meta_group, outpath = outpath)

spc_mean_df_melt <- melt(spc_mean_df,id.vars = "group_strain",variable.name = "wavenumber",value.name = "intensity")
spc_mean_df_melt <- filter(spc_mean_df_melt,group_strain == "HGE-Ss")
spc_mean_df_melt$intensity_sel <- spc_mean_df_melt$intensity
waves_sel_df <- a_peak
if (is.na(sum(waves_sel_df))) {
  next
}
logic <- spc_mean_df_melt$wavenumber %in% waves_sel_df
spc_mean_df_melt$intensity_sel[logic == FALSE] <- NA

spc_mean_df_melt$wavenumber <- as.numeric(as.character(spc_mean_df_melt$wavenumber))
plot_sel <- ggplot(spc_mean_df_melt, aes(x = wavenumber, y = intensity)) +
  geom_line(aes(group = group_strain)) +  
  geom_point(aes(x = wavenumber, y = intensity_sel),
             colour = "grey90",
             fill = "red",
             shape = 21,
             alpha = 0.8,
             size = 2.5
  ) +
  #ggtitle(paste(group_name, " = ", level, "\nmark = ", mark, sep = "")) +
  xlab(expression(paste("Wavenumber (cm"^"-1", ")", sep = ""))) +
  ylab("Normalized Intensity (a.u.)") +
  scale_x_continuous(breaks = seq(0, 4000, 200)) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    text = element_text(color = "black"),
    plot.title = element_text(hjust = 0.5, size = 15, face = "bold"),
    axis.title.x = element_text(size = 15),
    axis.title.y = element_text(size = 15),
    axis.text.x = element_text(size = 10),
    axis.text.y = element_text(size = 10),
    axis.ticks.x = element_line(size = 1),
    axis.ticks.y = element_line(size = 1)
  )
plot_sel
##################### 4. plot mean spectra
# hs_pre$Group <- factor(hs_pre$Group,levels = c("HGE-Pg","HGE-Pg-SnF2","HGE-Ss","HGE"))
hs_pre$Group <- hs_pre$group_strain
plot_mean_concentration <- plot_Meanspc_bygroup(hs_pre,0.005)
plot_mean_concentration


ggsave(filename = paste(outpath, "plot_mean_concentration-2-LZD.png", sep = "/"),
       plot = plot_mean_concentration,width = 20, height = 5, dpi = dpi)

##################### 4. PCA / tSNE
tsne_data <- hs_pre$spc
tsne_label <- paste(hs_pre$group_strain,hs_pre$group_rep,sep = "_")
tSNE_result_data <- cal_tSNE(tsne_data,tsne_label)
tSNE_result_data$label <- factor(hs_pre$group_strain)
tSNE_result_data$rep <- paste("rep",hs_pre$group_rep,sep = "_")
tSNE_result_plot <- tsne_plot(tSNE_result_data)
tSNE_result_plot
ggsave(
  filename = paste(outpath, "tSNE_mean_strain.png", sep = "/"),
  plot = tSNE_result_plot,
  width = 7, height = 5, dpi = dpi
)
