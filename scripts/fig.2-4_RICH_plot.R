library("RamanEx")
library("ggbreak")
func_lib()

cluster <- c("PG-红色","SnF2-浅蓝色","Ss-绿色","HGE-灰色") 

rootpath <- "/Users/chenrongze/Myprogrames/2_PFSup/郑晓姗/Ramanome (P&G)/all_data"
outpath <- func_initialization(rootpath, output_folder = "output_0115")
dpi <- 800

##################### 1. get spc, meta, info files
filepath <- paste(rootpath,"/output_1101/files",sep = "") 
Name_group <- c("group_strain", "B", "group_cell","D","group_rep")


###Fig2.A 平均光谱
melt_summary <- fread(paste(filepath,"/","Fig2_A.txt", sep=""),header = T,sep="\t")

color_cluster <- c("#A50026","#4575B4","#42B540FF","#BDBDBDFF")
melt_summary$group_strain <- factor(melt_summary$group_strain,levels = c("HGE-Pg","HGE-Pg-SnF2","HGE-Ss","HGE"))
plot_hyperspec <- ggplot(data = melt_summary,aes(x = wavenumber, y = value, group = factor(group_strain))) +
  geom_ribbon(aes(ymin = value - sd, ymax = value + sd, fill = factor(group_strain)),alpha = 0.3) +
  geom_line(aes(color = factor(group_strain)),size = 0.8) + 
  scale_x_continuous(breaks = seq(600,3200,200))+ default_theme()+
  scale_fill_manual(values = color_cluster,aesthetics = c("fill","color")) +
  scale_x_break(c(1750,2800),space = 0.3,scales = 0.3) +
  theme(
    legend.position = "none",
    axis.text.x.top = element_blank(),
    axis.ticks.x.top = element_blank(),
    axis.text.x = element_text(size = 15,angle = 45,hjust = 1,vjust = 1),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank()) +
  ylab("Normalized Intensity (a.u.)") +
  xlab(expression(paste('Wavenumber (cm'^{-1},')',sep = "")))
plot_hyperspec
ggsave(filename=paste(outpath,"/Fig2_A.svg", sep=""),
       plot=plot_hyperspec, limitsize=T,width = 8,height = 6)


#### Fig2.B 差谱600~1800
spc_melt_summary <- fread(paste(filepath,"/","Fig2_B .txt", sep=""),header = T,sep="\t")
spc_melt_summary$group_strain <- factor(spc_melt_summary$group_strain,levels = c("HGE-Pg","HGE-Pg-SnF2","HGE-Ss","Control"))
spc_melt_summary <- filter(spc_melt_summary,wavenumber < 1800,wavenumber > 600)

color_cluster <- c("#A50026","#4575B4","#42B540FF","#BDBDBDFF")
plot_pvalue <- ggplot(data = spc_melt_summary,aes(x = wavenumber, y = value, group = group_strain))+
  geom_ribbon(aes(ymin = value - sd, ymax = value + sd),alpha = 0.3,) +
  geom_line(aes(color = group_strain),size = 0.8) + 
  scale_fill_manual(values = color_cluster,aesthetics = c("fill","color"))+
  geom_hline(yintercept = 0.004,linetype = "dashed", colour = "gray",size = 0.8)+
  geom_hline(yintercept = 0,linetype = "dashed", colour = "gray",size = 0.8)+
  geom_hline(yintercept = 0.008,linetype = "dashed", colour = "gray",size = 0.8)+
  geom_hline(yintercept = 0.012,linetype = "dashed", colour = "gray",size = 0.8)+
  geom_hline(yintercept = 0.02,linetype = "dashed", colour = "white",size = 0.8)+
  # geom_vline(aes(xintercept=754.092), colour="grey", linetype="dashed",size = 1)+
  # geom_vline(aes(xintercept=783), colour="grey", linetype="dashed",size = 1)+
  # geom_vline(aes(xintercept=853), colour="grey", linetype="dashed",size = 1)+
  # geom_vline(aes(xintercept=938), colour="grey", linetype="dashed",size = 1)+
  # geom_vline(aes(xintercept=1003), colour="grey", linetype="dashed",size = 1)+
  # geom_vline(aes(xintercept=1030), colour="grey", linetype="dashed",size = 1)+
  # geom_vline(aes(xintercept=1091), colour="grey", linetype="dashed",size = 1)+
  # geom_vline(aes(xintercept=1126), colour="grey", linetype="dashed",size = 1)+
  # geom_vline(aes(xintercept=1165), colour="grey", linetype="dashed",size = 1)+
  # geom_vline(aes(xintercept=1208), colour="grey", linetype="dashed",size = 1)+
  # geom_vline(aes(xintercept=1225), colour="grey", linetype="dashed",size = 1)+
  # geom_vline(aes(xintercept=1250), colour="grey", linetype="dashed",size = 1)+
  # geom_vline(aes(xintercept=1307), colour="grey", linetype="dashed",size = 1)+
  # geom_vline(aes(xintercept=1338), colour="grey", linetype="dashed",size = 1)+
  # geom_vline(aes(xintercept=1369), colour="grey", linetype="dashed",size = 1)+
  # geom_vline(aes(xintercept=1398), colour="grey", linetype="dashed",size = 1)+
  # geom_vline(aes(xintercept=1432), colour="grey", linetype="dashed",size = 1)+
  # geom_vline(aes(xintercept=1446), colour="grey", linetype="dashed",size = 1)+
  # geom_vline(aes(xintercept=1503), colour="grey", linetype="dashed",size = 1)+
  # geom_vline(aes(xintercept=1562), colour="grey", linetype="dashed",size = 1)+
  # geom_vline(aes(xintercept=1586.8), colour="grey", linetype="dashed",size = 1)+
  # geom_vline(aes(xintercept=1627.4), colour="grey", linetype="dashed",size = 1)+
  # geom_vline(aes(xintercept=1657), colour="grey", linetype="dashed",size = 1)+
  scale_x_continuous(breaks = seq(600,3200,200))+ default_theme()+
  scale_y_continuous(breaks = seq(0,0.1,0.004))+ 
  ylab("Normalized Intensity (a.u.)") +
  xlab(expression(paste('Wavenumber (cm'^{-1},')',sep = ""))) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(size = 15,angle = 45,hjust = 1,vjust = 1),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank())
plot_pvalue
ggsave(filename=paste(outpath,"/Fig2_B.svg", sep=""),
       plot=plot_pvalue, limitsize=T,width = 8,height = 6)

###Fig2.C heatmap
RICH_data <- fread(paste(filepath,"/","Fig2_C.txt", sep=""),header = T,sep="\t")
RICH_data$RICH_group <- factor(RICH_data$RICH_group,levels = c("Control","HGE-Ss","HGE-Pg-SnF2","HGE-Pg"))
levels(factor(RICH_data$wavenumber))

RICH_data <- RICH_data %>% dplyr::group_by(RICH_group) %>% 
  dplyr::mutate(A = rep(c(1:length(filename))[order(value[wavenumber == "1564"])],times = 325)) #

RICH_data <- RICH_data[order(RICH_data$A,RICH_data$RICH_group),]
RICH_data$filename <- factor(RICH_data$filename,levels = RICH_data$filename[!duplicated(RICH_data$filename)])


color_pigment_props <- c("white","#f0f7f2","#e4e7ed","lightgrey",
                         "#74ADD1","#1e592b","#B2DF8A","#FEE090",
                         "#FDAE61","#F46D43","#dc3e32","#A50026")
myPal <- colorRampPalette(color_pigment_props)(n=1000)
plot_tile <- ggplot(RICH_data, aes(x=factor(wavenumber), y=filename, fill=value)) + 
  geom_tile() + theme_bw() + coord_equal(ratio = 1) +
  scale_x_discrete(breaks = seq(600,1700,100)) +
  scale_fill_gradientn(colours =myPal) +
  #geom_text(aes(label=Freq), color=text_color) +# printing values
  theme(panel.grid = element_blank(),
        panel.border = element_blank(),
        legend.title = element_blank(),
        legend.key.width = unit(10,"pt"),
        legend.key.height = unit(80, "pt"),
        #legend.position = "nane",
        axis.ticks.y = element_blank(),
        axis.title = element_blank(),
        axis.text.x = element_text(size = 10,angle = 45,hjust = 0.8,vjust =1), 
        axis.text.y = element_blank(),        
        axis.text = element_text()
  )
plot_tile

ggsave(filename=paste(outpath,"/Fig2_C.png", sep=""),dpi = 300,
       plot=plot_tile, limitsize=T,width = 8,height = 8)

#### fig3—B
spc_melt_summary <- fread(paste(filepath,"/","Fig2_B.txt", sep=""),header = T,sep="\t")
spc_melt_summary$group_strain <- factor(spc_melt_summary$group_strain,levels = c("HGE-Pg","HGE-Pg-SnF2","HGE-Ss","Control"))
spc_melt_summary <- filter(spc_melt_summary,wavenumber < 1800,wavenumber > 600,group_strain == "HGE-Pg")

colnames(spc_melt_summary)
colnames(p_value_Data_i) <- c("p_value","wavenumber","group_strain","N","value","sd","se","ci","diff") 

all_data_table <- left_join(spc_melt_summary,p_value_Data_i[,1:3],by = c("wavenumber","group_strain"))
all_data_table$value <- all_data_table$value + 0.002

myPal <- c("#08306B","#313695","#4575B4","#74ADD1","white","#FEE090","#FEE090","#FEE090","#FEE090",
           "#FDAE61","#FDAE61","#FDAE61","#FDAE61","#FDAE61",
           "#F46D43","#F46D43","#F46D43","#F46D43","#F46D43","#F46D43",
           "#A50026","#A50026","#A50026","#A50026","#A50026","#A50026")
#myPal <- colorRampPalette(color_pigment_props)(n=900)
plot_pvalue <- ggplot(data = all_data_table,aes(x = wavenumber, y = value, group = group_strain))+
  #geom_point(aes(size = abs(p_value)),alpha=0.7) + scale_size(range=c(1,3)) + 
  geom_ribbon(aes(ymin = value - sd, ymax = value + sd),alpha = 0.6) +
  geom_tile(aes(height = 0.012,y=0.005, fill=value),alpha = 0.5) + 
  geom_line(color = "black",size = 0.8) + 
  geom_hline(yintercept = 0.002,linetype = "dashed", colour = "gray",size = 0.8)+
  #scale_fill_manual(values = color_cluster,aesthetics = c("fill","color"))+
  scale_fill_gradientn(colours =myPal) +
  scale_x_continuous(breaks = seq(600,3200,200))+ default_theme()+
  labs(fill = "Weight") + 
  ylab("Normalized Intensity (a.u.)") +
  xlab(expression(paste('Wavenumber (cm'^{-1},')',sep = ""))) +
  theme(
    #legend.position = "none",
    legend.key.width = unit(10,"pt"),
    legend.key.height = unit(60, "pt"),
    axis.text.x = element_text(size = 15,angle = 45,hjust = 1,vjust = 1),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank())
plot_pvalue
ggsave(filename=paste(outpath,"/Fig3_B.svg", sep=""),
       plot=plot_pvalue, limitsize=T,width = 8,height = 6)
