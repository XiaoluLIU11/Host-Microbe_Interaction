func_spc2hs <- function(spc,meta)
{
  spc_ori <- spc %>% filter(filename %in% meta$filename)
  meta_ori <- meta %>% filter(filename %in% spc$filename)
  wavenumber <- as.numeric(colnames(spc_ori[,-1]))
  spc2hs <- new("hyperSpec",data = meta_ori,spc = spc_ori[,-1],wavelength = wavenumber)
  
  return(spc2hs)
}
read_3data <- function(folderpath){
  all_spc <- fread(paste(folderpath,"/alldata_spc.txt",sep = ""), header = TRUE, sep = "\t")
  all_meta <- fread(paste(folderpath,"/alldata_meta.txt",sep = ""), header = TRUE, sep = "\t")
  data_info <- fread(paste(folderpath,"/data_info.txt",sep = ""), header = TRUE, sep = "\t")
  return(list(all_spc,all_meta,data_info))
}
data_clean <- function(spc,meta)
{
  spc_reor_filt_qc <- cal_QC(spc)
  spc_qc_meta <- left_join(spc_reor_filt_qc,meta,by = "filename")
  spc_meta_good <- spc_qc_meta %>% filter(Noise_mean < 0.1, Noise_sd < 0.1,Signal_max > 0.1)
  hs_ori <- func_spc2hs(spc,spc_meta_good)
  return(hs_ori)
}
pre_allpipeline <- function(hs_ori)
{
  #hs_baseline_sg <- pre_sg(hs_ori[,-1])
  hs_baseline_sg <- spc_waveumber(hs_ori,500,1750)
  hs_baseline <- pre_baseline(hs_baseline_sg,1)
  hs_baseline_sg_norm <- pre_norm(hs_baseline,method = "sum")
  #hs_baseline_sg_norm <- cbind(filename = hs_ori$filename,hs_baseline_sg_norm)
  return(hs_baseline_sg_norm)
}
plot_Meanspc_bygroup <- function(hyperspec,a = 0.0005)
{
  hyperspec_melt <- spc_melt(hyperspec,c("group_strain"),500,3200)
  hyperspec_melt_summary <- Rmisc::summarySE(hyperspec_melt, measurevar = "value", groupvars = c("group_strain","wavenumber"))
  n <- length(levels(factor(hyperspec_melt_summary$group_strain))) * a
  for(i in levels(factor(hyperspec_melt_summary$group_strain)))
  {
    inds <- which(hyperspec_melt_summary$group_strain == i )
    hyperspec_melt_summary[inds,]$value <- hyperspec_melt_summary[inds,]$value + n
    n <- n - a
  }
  
  plot_hyperspec <- ggplot(data = hyperspec_melt_summary,aes(x = wavenumber, y = value, group = factor(group_strain))) +
    geom_ribbon(aes(ymin = value - sd, ymax = value + sd, fill = factor(group_strain)),alpha = 0.3) +
    geom_line(aes(color = factor(group_strain)),size = 0.8) + theme_bw() + #facet_grid(. ~ label) +
    labs(y = "Normalized Intensity (a.u.)") + xlab(expression(paste('Wavenumber (cm'^{-1},')',sep = ""))) +
    scale_x_continuous(breaks = seq(500,3200,200)) +
    theme(
      panel.grid = element_blank(),
      panel.grid.minor = element_blank(),
      legend.title = element_blank(),
      legend.text = element_markdown(size = 15),
      #legend.position = "none",
      legend.background = element_blank(),
      text = element_text(color = "black"),
      axis.title.y = element_text(size = 20, angle = 90),
      axis.text.x = element_text(size = 15,angle = 45,hjust = 1,vjust = 1),
      axis.text.y = element_blank(),
      strip.background = element_rect(fill = "white"),
      strip.text = element_text(size = 20),
      axis.ticks = element_line(size = 1),
      axis.ticks.y = element_blank(),
      axis.ticks.length = unit(0.4,"lines"),
      axis.title = element_text(size = 20))
  plot_hyperspec
  
  return(plot_hyperspec)
}

org_model <- function(all_data_spc)
{
  train_data <- as.data.table(all_data_spc)
  factor_alldata <- select(train_data,group_Strain) #记得更改这里
  factor_alldata <- factor_alldata[!duplicated(factor_alldata),]
  training_set <- c()
  setkey(train_data,group_Strain) #和这里
  for (i in 1:nrow(factor_alldata))
  {
    training_set_i <- sample(train_data[factor_alldata[i,]],(0.7*nrow(train_data[factor_alldata[i,]]))%/%1)
    training_set <- rbind(training_set,training_set_i)
  }
  test_set <- train_data
  test_set <- filter(test_set, Number %in% setdiff(test_set$Number,training_set$Number))
  return(list(training_set,test_set))
}
tsne_plot <- function (tSNE_data) 
{
  plot <- ggplot(tSNE_data, aes(tSNE1, tSNE2, color = label,shape = factor(rep))) + 
    geom_point() + theme_bw() + stat_ellipse(level = 0.8,linetype = 2) + 
    theme(panel.grid = element_blank(), 
          panel.grid.minor = element_blank(), 
          legend.title = element_blank(), 
          legend.text = element_markdown(size = 15),
          legend.background = element_blank(), 
          text = element_text(color = "black"), 
          axis.text.x = element_text(size = 15,angle = 0), 
          axis.text.y = element_text(size = 15),
          axis.ticks.length = unit(0.4, "lines"),
          axis.title = element_text(size = 15))
  return(plot)
}


RICH_pipeline <- function(spc,meta)
{
  spc_mean_df_meta <- left_join(spc, meta, by = "Group")
  spc_mean_df_meta <- dplyr::select(spc_mean_df_meta,c(names(spc_mean_df_meta)[grep(pattern = "group", names(spc_mean_df_meta))],names(spc)[-1]))
  spc_mean_df_meta$group_D <- gsub("h","",spc_mean_df_meta$group_D)
  spc_mean_df_meta$group_D <- as.numeric(as.character(spc_mean_df_meta$group_D))
  
  mark1 <- 0.9
  mark2 <- 0.6
  
  plot_list <- mark1_list <- mark2_list <- list()
  
  for(j in levels(factor(spc_mean_df_meta$group_c)))
  {
    spc_mean_df_i_j <- filter(spc_mean_df_meta,group_c == j)
    i_j <- j
    cor_list <- cal_cor(spc_mean_df_i_j, "group_D", 0.9, 0.6)
    plot_list[[i_j]] <- cor_list$plot_cor +
      ggtitle(i_j) +
      theme(
        axis.title.x = element_text(size = 10),
        axis.title.y = element_text(size = 10),
        axis.text.x = element_text(size = 7),
        axis.text.y = element_text(size = 7),
        plot.title = element_text(hjust = 0.5, size = 10, face = "bold"),
      )
  }
  return(plot_list)
} 

org_peak_all <- function(peak_list)
{
  peak_all <- c()
  for(i in 1:length(peak_list))
  {
    peak_all <- c(peak_all,peak_list[[i]])
  }
  peak_all <- peak_all[!duplicated(peak_all)]
  peak_all <- sort(peak_all)
  inds <- c()
  for(i in 1:(length(peak_all)-1))
  {
    if(peak_all[i]+3 > peak_all[i+1])
    {
      inds[i] <- FALSE
    }else{
      inds[i] <- TRUE
    }
  }
  peak_all <- peak_all[c(inds,TRUE)]
  return(peak_all)
}
plot_findpeak <- function(data_table_mean,peak_num,stack_height = 0.001)
{
  hs_spc <- data_table_mean
  spc_mean_df_melt <- data.table(melt(hs_spc,id.vars = "Group",variable.name = "wavenumber",value.name = "intensity"))
  setkey(spc_mean_df_melt,Group,wavenumber)
  Group_levels <- levels(factor(spc_mean_df_melt$Group))
  
  n <- length(Group_levels) * stack_height
  for (i in 1:length(Group_levels)) {
    spc_mean_df_melt[which(spc_mean_df_melt$Group == Group_levels[i]), ]$intensity <-
      spc_mean_df_melt[which(spc_mean_df_melt$Group == Group_levels[i]), ]$intensity + n
    n <- n - stack_height
  }
  spc_mean_df_melt$wavenumber <- as.numeric(as.character(spc_mean_df_melt$wavenumber))
  spc_mean_df_melt$intensity_sel <- spc_mean_df_melt$intensity
  
  ogic <- spc_mean_df_melt$wavenumber %in% peak_num
  spc_mean_df_melt$intensity_sel[ogic == FALSE] <- NA
  
  
  plot_sel <- ggplot(spc_mean_df_melt, aes(x = wavenumber, y = intensity)) +
    geom_line(aes(group = Group)) +  
    geom_point(aes(x = wavenumber, y = intensity_sel),
               colour = "grey90", fill = "red", shape = 21, alpha = 0.8, size = 2.5) +
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
  return(plot_sel)
}
org_findpeak <- function(data_table_mean,peak_num){
  hs_spc <- data_table_mean
  spc_mean_df_melt <- data.table(melt(hs_spc,id.vars = "Group",variable.name = "wavenumber",value.name = "intensity"))
  setkey(spc_mean_df_melt,Group,wavenumber)
  Group_levels <- levels(factor(spc_mean_df_melt$Group))
  
  spc_mean_df_melt_peak <- spc_mean_df_melt %>% filter(wavenumber %in% peak_num)
  spc_mean_df_peak <- dcast(spc_mean_df_melt_peak,Group ~ wavenumber)
  return(spc_mean_df_peak)
}
plot_peaktrend_conc <- function(spc,meta,outpath)
{
  wavenume <- as.numeric(as.character(colnames(spc[,-1])))
  melt_spc <- data.table(melt(spc,id.vars = "Group",variable.name = "wavenumber",value.name = "intensity"))
  melt_spc_meta <- left_join(melt_spc,meta,by = "Group")
  melt_spc_meta <- data.frame(melt_spc_meta)
  for(j in 1:length(wavenume))
  {
    melt_spc_j <- dplyr::filter(melt_spc_meta,wavenumber == as.character(wavenume[j]))
    melt_spc_j$group_D <- factor(melt_spc_j$group_D,levels = c("0h","6h","12h","24h","48h","72h"))
    plot_list <- ggplot(melt_spc_j, aes(x = group_D, y = intensity,group = factor(group_c))) + 
      theme_bw() + geom_line(aes(group = factor(group_c), color = factor(group_c))) + 
      xlab("group_D") + 
      ylab("Metabolic activity (CD-ratio)") +
      ggtitle(wavenume[j]) +
      theme(
        axis.title.x = element_text(size = 10),
        axis.title.y = element_text(size = 10),
        axis.text.x = element_text(size = 7),
        axis.text.y = element_text(size = 7),
        plot.title = element_text(hjust = 0.5, size = 10, face = "bold"),
      )
    ggsave(filename = paste(outpath,"/",j,"_peaktrend_by_conc.png", sep = ""),plot = plot_list,
           width = 8, height = 6)
  }
}
plot_peaktrend_time <- function(spc,meta,outpath)
{
  wavenume <- as.numeric(as.character(colnames(spc[,-1])))
  melt_spc <- data.table(melt(spc,id.vars = "Group",variable.name = "wavenumber",value.name = "intensity"))
  melt_spc_meta <- left_join(melt_spc,meta,by = "Group")
  melt_spc_meta <- data.frame(melt_spc_meta)
  for(j in 1:length(wavenume))
  {
    melt_spc_j <- dplyr::filter(melt_spc_meta,wavenumber == as.character(wavenume[j]))
    melt_spc_j$group_D <- factor(melt_spc_j$group_D,levels = c("0h","6h","12h","24h","48h","72h"))
    plot_list <- ggplot(melt_spc_j, aes(x = factor(group_c), y = intensity,group = group_D)) + 
      theme_bw() + geom_line(aes(group = group_D, color = group_D)) + 
      xlab("group_D") + 
      ylab("Metabolic activity (CD-ratio)") +
      ggtitle(wavenume[j]) +
      theme(
        axis.title.x = element_text(size = 10),
        axis.title.y = element_text(size = 10),
        axis.text.x = element_text(size = 7),
        axis.text.y = element_text(size = 7),
        plot.title = element_text(hjust = 0.5, size = 10, face = "bold"),
      )
    ggsave(filename = paste(outpath,"/",j,"_peaktrend_by_time.png", sep = ""),plot = plot_list,
           width = 8, height = 6)
  }
}
plot_peaktrend_all <- function(spc,meta_filt,outpath)
{
  peak_list_all <- cal_findpeak(spc, step = 0.00001)
  peak_all_i <- org_peak_all(peak_list_all)
  spc_peak <- org_findpeak(spc,peak_all_i)
  
  trend_Conc_outpath <- func_initialization(outpath, output_folder = "trend_conc_output")
  plot_FP_list <- plot_peaktrend_conc(spc = spc_peak,meta_filt,trend_Conc_outpath)
  trend_Time_outpath <- func_initialization(outpath, output_folder = "trend_time_output")
  plot_FP_list <- plot_peaktrend_time(spc = spc_peak,meta_filt,trend_Time_outpath)
}

tsne_plot <- function (tSNE_data) 
{
  plot <- ggplot(tSNE_data, aes(tSNE1, tSNE2, color = label)) + 
    geom_point() + theme_bw() + stat_ellipse(level = 0.8) + 
    theme(panel.grid = element_blank(), 
          panel.grid.minor = element_blank(), 
          legend.title = element_blank(), 
          legend.text = element_markdown(size = 15),
          legend.background = element_blank(), 
          text = element_text(color = "black"), 
          axis.text.x = element_text(size = 15,angle = 0), 
          axis.text.y = element_text(size = 15),
          axis.ticks.length = unit(0.4, "lines"),
          axis.title = element_text(size = 15))
  return(plot)
}
cal_RICH_tt <- function(all_spc,Pvalue = 0.0001)
{
  all_spc_control <- filter(all_spc,RICH_group == "Control")
  all_spc_test <- filter(all_spc,RICH_group == "HGE-Pg")
  
  melt_summary <- Rmisc::summarySE(all_spc, measurevar = "value",groupvars = c("RICH_group", "wavenumber"))
  melt_summary_a <- melt_summary %>%
    group_by(wavenumber)  %>%
    mutate(diff = value - value[RICH_group == "Control"]) #%>%
    #filter(RICH_group != "control")
  
  wavenumber <- levels(factor(all_spc$wavenumber))
  p_value <- select_peak <- NULL
  for (i in wavenumber)
  {
    i_control <- filter(all_spc_control,wavenumber == i)
    i_test <- filter(all_spc_test,wavenumber == i)
    p <- t.test(i_control$value,i_test$value)$p.value
    if (p < Pvalue) {
      select_peak_i <- i
      select_peak <- c(select_peak,select_peak_i)
    }
    p_value <- c(p_value,p)
  }
  p_value <- as.data.frame(cbind(p_value,wavenumber))
  p_value$wavenumber <- as.numeric(as.character(p_value$wavenumber))
  p_value$p_value <- log10(as.numeric(as.character(p_value$p_value)))
  melt_summary_a$wavenumber <- as.numeric(as.character(melt_summary_a$wavenumber))
  p_value <- left_join(p_value,melt_summary_a,by="wavenumber")
  return(list(select_peak,p_value))
}

org_model <- function(all_data_spc)
{
  train_data <- as.data.table(all_data_spc)
  factor_alldata <- select(train_data,RICH_group) #记得更改这里
  factor_alldata <- factor_alldata[!duplicated(factor_alldata),]
  training_set <- c()
  setkey(train_data,RICH_group) #和这里
  for (i in 1:nrow(factor_alldata))
  {
    training_set_i <- sample(train_data[factor_alldata[i,]],(0.7*nrow(train_data[factor_alldata[i,]]))%/%1)
    training_set <- rbind(training_set,training_set_i)
  }
  test_set <- train_data
  test_set <- filter(test_set, filename %in% setdiff(test_set$filename,training_set$filename))
  return(list(training_set,test_set))
}

func_diff <- function (tb_data, meta) 
{
  tb_data <- left_join(meta, tb_data, by = "filename")
  hyperspec_melt <- spc_melt(tb_data, "group_strain", 500, 3200)
  hyperspec_melt$value <- as.numeric(as.character(hyperspec_melt$value))

  hyperspec_melt_control <- hyperspec_melt %>% filter(group_strain == "HGE")
  hyperspec_melt_summary <- Rmisc::summarySE(hyperspec_melt, measurevar = "value", groupvars = c("group_strain", "wavenumber"))
  hyperspec_mean <- Rmisc::summarySE(hyperspec_melt_control, measurevar = "value", groupvars = c("wavenumber"))
  hyperspec_melt_summary$value <- hyperspec_melt_summary$value - hyperspec_mean$value
  return(hyperspec_melt_summary)
}
