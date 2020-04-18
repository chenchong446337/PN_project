## The script are used to analyze the data from miniscope experiemnt analyzed by matlab
## updated on 01-28-2020
## data are analyzed with inscopix data processing images
## output: 1. The trace crossing the border; 2. cells proportion increased or decrease; 3. Heat maps
## 4. Number of cells active, inhibited
## updated: 03182020, add more data
## updated: 04012020, change the way to do z-score
## updated: 04142020, do k-mean analysis firstly and back cat the cells activation

## import needed library-----
library("ggplot2")
library("reshape2")
library("plyr")
library(openxlsx)
library(cowplot)
library(grid)
library(DescTools)
library(viridis)
library(tidyverse)
library(corrplot)
library(scales)
library(ggcorrplot)


## function to extract trace crossing the border
c_miniscope_matlab <- function(ID_trace, t_stim) {
  ## import and format the data
  dat_trace <- read.xlsx(ID_trace, colNames = F, rowNames = F)
  ## analyze the trace
  dat_trace<- dat_trace[dat_trace[,1]==1, ]
  dat_trace<- as.data.frame(t(dat_trace[,-1]))
  rownames(dat_trace) <- NULL
  colnames(dat_trace)<- NULL
  ## do z score of the whole trace
  dat_trace <- apply(dat_trace, 2, scale)
  
  ## divide ca2+ trace within -5 and 9s before and after stimuls, take -5 to -3 as baseline
  dat_stim <- vector(mode = "list", length = length(t_stim))
  stim_time<- seq(-5, 8.5, by=0.5)
  ## number of rows to be binned
  n <- 10 # 0.05*10=0.5
  
  ## extract cell activity when they cross the border
  for (i in seq_along(t_stim)){
    t1_p <- t_stim[i]
    dat_stim1 <- dat_trace[(t1_p-100):(t1_p+180-1),]
    dat_stim1 <- aggregate(dat_stim1,list(rep(1:(nrow(dat_stim1)%/%n+1),each=n,len=nrow(dat_stim1))),mean)[-1]
    dat_stim1_base <- dat_stim1[1:5,] ## baseline as -5 to -3
    dat_stim1_base_mean <- colMeans(dat_stim1_base, na.rm = T)
    dat_stim1_nor1 <- sweep(dat_stim1, 2, dat_stim1_base_mean, FUN = "-")
    rownames(dat_stim1_nor1) <- NULL
    colnames(dat_stim1_nor1)<- NULL
    dat_stim[[i]] <- dat_stim1_nor1
  }
  
  ## average the trace by cell number
  if (length(t_stim)>1){
    dat_cell_trace_average <- aaply(laply(dat_stim, as.matrix), c(2, 3), function(x) mean(x, na.rm=T))
    
  } else {
    dat_cell_trace_average <- data.matrix(dat_stim[[1]])
    
  }
  
  return(dat_cell_trace_average)
  
}

## Extract trace from each mice------
# for m3
path_trace_m3 <- as.list(list.files(path ="~cchen2/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/matlab_analysis/m3/trace_days/", pattern = "*.xlsx", full.names = T ))
t_stim_m3_d1 <- c(214, 698,1204,1757,2964) 
t_stim_m3_d2 <- c(274, 831,1364,1788,2621,3273)
t_stim_m3_d3 <- c(1885, 2516) 
t_stim_m3_d4 <- c(1813,3587)
t_stim_m3_d5 <- c(758) 
t_stim_m3_d6 <- c(132,507) 
t_stim_m3_d7 <- c(346) ## use first two 346, 1206, 1756, 2288

t_stim_m3 <- list(t_stim_m3_d1, t_stim_m3_d2, t_stim_m3_d3, t_stim_m3_d4, t_stim_m3_d5, t_stim_m3_d6, t_stim_m3_d7)

dat_trace_m3 <- mapply(c_miniscope_matlab, path_trace_m3, t_stim_m3, SIMPLIFY = F)

# for m7
path_trace_m7 <- as.list(list.files(path ="~cchen2/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/matlab_analysis/m7/trace_days/", pattern = "*.xlsx", full.names = T ))
t_stim_m7_d1 <- c( 605, 964, 1603, 2153, 2363, 2751, 3152) 
t_stim_m7_d2 <- c(798, 177,3300)
t_stim_m7_d3 <- c(360, 3352, 4044) 
t_stim_m7_d4 <- c(3236)
t_stim_m7_d5 <- c(317) 
t_stim_m7_d6 <- c(205)
t_stim_m7_d7 <- c(236) ## use first two 236, 1914,2777,3334,4091

t_stim_m7 <- list(t_stim_m7_d1, t_stim_m7_d2, t_stim_m7_d3, t_stim_m7_d4, t_stim_m7_d5, t_stim_m7_d6, t_stim_m7_d7)

dat_trace_m7 <- mapply(c_miniscope_matlab, path_trace_m7, t_stim_m7,SIMPLIFY = F)

# for m17
path_trace_m17 <- as.list(list.files(path ="~cchen2/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/matlab_analysis/m17/trace_days/", pattern = "*.xlsx", full.names = T ))
t_stim_m17_d1 <- c(815,1202,1656,2295,2580,3178) 
t_stim_m17_d2 <- c(1103,1678,2291,2853,3839)
t_stim_m17_d3 <- c(437,850,1210,1728,2927) 
t_stim_m17_d4 <- c(682,1745)
t_stim_m17_d5 <- c(1114) 
t_stim_m17_d6 <- c(617) 
t_stim_m17_d7 <- c(157) # 157, 2536, 2961

t_stim_m17 <- list(t_stim_m17_d1, t_stim_m17_d2, t_stim_m17_d3, t_stim_m17_d4, t_stim_m17_d5, t_stim_m17_d6, t_stim_m17_d7)

dat_trace_m17 <- mapply(c_miniscope_matlab, path_trace_m17, t_stim_m17, SIMPLIFY = F)

# for m18
path_trace_m18 <- as.list(list.files(path ="~cchen2/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/matlab_analysis/m18/trace_days/", pattern = "*.xlsx", full.names = T ))
t_stim_m18_d1 <- c(490, 950,1291,1714,2205,2610,3088,3521) 
t_stim_m18_d2 <- c(236, 898,1528,2419)
t_stim_m18_d3 <- c(784,1271,2012,3088) 
t_stim_m18_d4 <- c(493,1580,3874)
t_stim_m18_d5 <- c(195,3676) 
t_stim_m18_d6 <- c(467) 
t_stim_m18_d7 <- c(124) ## use the first two of 124,2238,3000

t_stim_m18 <- list(t_stim_m18_d1, t_stim_m18_d2, t_stim_m18_d3, t_stim_m18_d4, t_stim_m18_d5, t_stim_m18_d6, t_stim_m18_d7)

dat_trace_m18 <- mapply(c_miniscope_matlab, path_trace_m18, t_stim_m18, SIMPLIFY = F)

## combine data and do k-means analysis-----
dat_cell_trace <- vector(mode = "list", 7)
for (i in 1:7) {
  dat_cell_trace[[i]]<- cbind(dat_trace_m3[[i]], dat_trace_m7[[i]], dat_trace_m17[[i]], dat_trace_m18[[i]] )
  
}

dat_cell_trace_re <- vector(mode = "list", 7)
for (i in 1:length(dat_cell_trace)) {
  dat_cell_trace_d <- dat_cell_trace[[i]]
  cell_id <- sprintf("Cell%s",seq(1:ncol(dat_cell_trace_d)))
  colnames(dat_cell_trace_d)<- cell_id
  
  dat_cell_trace_d_cluster <- unname(kmeans(t(dat_cell_trace_d), centers = 3)[[1]])
  dat_cell_trace_d<- as.data.frame(dat_cell_trace_d)
  stim_time <- seq(-5, 8.5, by=0.5)
  dat_cell_trace_d$Time <- stim_time
  dat_cell_trace_d_re <- melt(dat_cell_trace_d, id.vars ='Time')
  dat_cell_trace_d_re$Group <- rep(dat_cell_trace_d_cluster, each=length(stim_time))

  ## sort data by the value in each group
  dat_cell_d_sort <- as.numeric(names(sort(tapply(dat_cell_trace_d_re$value, dat_cell_trace_d_re$Group, mean), decreasing = T)))
  
  dat_cell_trace_d_re$Group[dat_cell_trace_d_re$Group == dat_cell_d_sort[1]] ="Excited"
  dat_cell_trace_d_re$Group[dat_cell_trace_d_re$Group == dat_cell_d_sort[2]] ="Neutral"
  dat_cell_trace_d_re$Group[dat_cell_trace_d_re$Group == dat_cell_d_sort[3]] ="Inhibited"
  
  dat_cell_trace_d_re$Group <- factor(dat_cell_trace_d_re$Group, levels = c("Excited", "Neutral", "Inhibited"))
  rep_time <- c(ncol(dat_trace_m3[[i]]), ncol(dat_trace_m7[[i]]), ncol(dat_trace_m17[[i]]), ncol(dat_trace_m18[[i]]))
  dat_cell_trace_d_re$ID <- rep(c("m3", "m7", "m17", "m18"), rep_time*length(stim_time))
  
  dat_cell_trace_re[[i]] <- dat_cell_trace_d_re
}

## heat plot of d3, d6 and d7----
## set the range of z score during all days
score_range <- range(mapply(function(x) x$value, dat_cell_trace_re, SIMPLIFY = T))

for (i in c(3, 6, 7)) {
  dat_trace <- dat_cell_trace_re[[i]]
  p_heat <- ggplot(dat_trace, aes(Time, variable,fill= value))+ 
    geom_tile(height=2)+
    facet_grid(rows = vars(Group), scales = "free_y")+
    #scale_fill_gradient2(limits= score_range,low = "navy", high = "red4", mid = "white")+
    scale_fill_gradientn(limits= score_range, colours = c("navy", "white", "red4"), values = rescale(c(score_range[1], 0, score_range[2])))+
    labs(x="Time relative to crossing (s)", y="Number of cells")+
    geom_vline(xintercept = 0, col= 'red', linetype=2)+
    theme(axis.line.x = element_line(),
          axis.line.y = element_line(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank(),
          axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
    theme(legend.position = "none")
  assign(paste0("p_heat_d", i), p_heat)
}

## combine the heat plot
p_heat_com <- plot_grid(p_heat_d3, p_heat_d6, p_heat_d7, nrow = 1)

## plot trace by group-----
dat_cell_trace_pre_sta <- ddply(rbind(dat_cell_trace_re[[2]],dat_cell_trace_re[[3]]), .(Time, Group), summarise,n=length(value),mean=mean(value),sd=sd(value),se=sd(value)/sqrt(length(value)))
dat_cell_trace_con_sta <- ddply(rbind(dat_cell_trace_re[[5]], dat_cell_trace_re[[6]]), .(Time, Group), summarise,n=length(value),mean=mean(value),sd=sd(value),se=sd(value)/sqrt(length(value)))
dat_cell_trace_test_sta <- ddply(dat_cell_trace_re[[7]], .(Time, Group), summarise,n=length(value),mean=mean(value),sd=sd(value),se=sd(value)/sqrt(length(value)))

dat_cell_trace_sta <- rbind(dat_cell_trace_pre_sta, dat_cell_trace_con_sta, dat_cell_trace_test_sta)

dat_cell_trace_sta$Day <- rep(c("Pre", "Cond.", "Test"),each= length(stim_time)*3)
dat_cell_trace_sta$Day <- factor(dat_cell_trace_sta$Day, levels = c("Pre", "Cond.", "Test"))
dat_cell_trace_sta$Group <- factor(dat_cell_trace_sta$Group, levels = c("Neutral", "Excited", "Inhibited"))
range_trace_plot <- range(dat_cell_trace_sta$mean)

## group by day
p_trace <- ggplot(dat_cell_trace_sta, aes(Time, mean, colour=Group))+
  geom_line()+
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), alpha=0.3, colour="gray")+
  facet_grid(cols = vars(Day))+
  scale_colour_manual(values=c("#999999", "darkred", "navy"))+
  labs(x="Time relative to crossing (s)", y="Z score")+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  geom_vline(xintercept = 0, col= 'red', linetype=2)+
  scale_y_continuous(limits = c(range_trace_plot[1]- 0.2, range_trace_plot[2]+0.2))+
  theme(legend.title = element_blank())

## for EI ratio analysis-----
dat_cell_area <- c()
cell_area_day <- c("Rm","Pre", "Pre", "Rm","Cond.","Cond.", "Test")
for (i in c(2,3,5,6,7)){
  rep_time <- c(ncol(dat_trace_m3[[i]]), ncol(dat_trace_m7[[i]]), ncol(dat_trace_m17[[i]]), ncol(dat_trace_m18[[i]]))
  dat_trace <- dat_cell_trace_re[[i]]
  dat_trace$ID <- rep(c("m3", "m7", "m17", "m18"), rep_time*length(stim_time))
  dat_trace_sta <- ddply(dat_trace, .(ID, Time, Group), summarise,n=length(value),mean=mean(value),sd=sd(value),se=sd(value)/sqrt(length(value)))
  Time_area <- which(stim_time>=0)
  dat_anti_area <- as.data.frame(tapply(dat_trace_sta$mean, INDEX = list(dat_trace_sta$ID, dat_trace_sta$Group), 
                                        function (x) AUC(stim_time[Time_area], x[Time_area])))
  dat_anti_area$ID <- rownames(dat_anti_area)
  rownames(dat_anti_area)<-NULL
  dat_anti_area$Day <- cell_area_day[i]
  dat_anti_area$ratio <- abs(dat_anti_area$Excited/dat_anti_area$Inhibited)
  dat_cell_area <- rbind(dat_cell_area, dat_anti_area)
  }

dat_cell_area$Day <- factor(dat_cell_area$Day, levels = c("Pre", "Cond.","Test"))
dat_cell_area_com <- ddply(dat_cell_area, .(ID, Day), summarise, value=mean(ratio, na.rm = T))
dat_cell_area_sta <- ddply(dat_cell_area_com,.(Day), summarise,n=length(value),mean=mean(value),sd=sd(value),se=sd(value)/sqrt(length(value)))


p_EI_ratio <- ggplot(dat_cell_area_sta, aes(Day, mean, colour=Day, fill=Day))+
  geom_bar(stat="identity", color="black", position=position_dodge())+
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.1, position=position_dodge(.9), colour="black") +
  scale_fill_manual(values=c( "mediumseagreen","indianred4", "dodgerblue4"))+
  labs(x="", y="E/I ratio")+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  scale_y_continuous(limits = c(0, 3), expand = c(0,0))+
  theme(legend.title = element_blank())+
  theme(legend.position = "none")

setwd("~cchen2/Documents/neuroscience/Pn\ project/Figure/PDF/")
cairo_pdf("p_EI_ratio.pdf", width = 50/25.6, height = 65/25.6, family = "Arial")
p_EI_ratio
dev.off()
## portion of cell catlog------
dat_cell_cat <- mapply(function (x) subset(x, x$Time==0), dat_cell_trace_re, SIMPLIFY = F )

dat_cell_cat_table <- mapply(function(x) prop.table(table(x$ID, x$Group), 1), dat_cell_cat, SIMPLIFY = F)
dat_cell_cat_com <- as.data.frame(do.call(rbind, dat_cell_cat_table))
dat_cell_cat_com$ID <- rep(c(rownames(dat_cell_cat_com)[1:4]), 7)
dat_cell_cat_com$Day <- rep(c("Rm", "Pre", "Pre", "Rm", "Cond.", "Cond.", "Test"), each=4)

dat_cell_cat_com <- subset(dat_cell_cat_com, dat_cell_cat_com$Day!="Rm")

dat_cell_cat_com_re <- melt(dat_cell_cat_com, id.vars = c("ID", "Day"))
dat_cell_cat_com_sta <- ddply(dat_cell_cat_com_re, .(ID, Day, variable), summarise, value=mean(value))
dat_cell_cat_com_sta1 <- ddply(dat_cell_cat_com_sta, .(Day, variable), summarise,n=length(value),mean=mean(value),sd=sd(value),se=sd(value)/sqrt(length(value)))

dat_cell_cat_com_sta1$Day <- factor(dat_cell_cat_com_sta1$Day, levels = c("Pre", "Cond.", "Test"))
dat_cell_cat_com_sta1$variable <- factor(dat_cell_cat_com_sta1$variable, levels = c("Excited", "Neutral", "Inhibited"))

p_cat<- ggplot(dat_cell_cat_com_sta1, aes(Day, mean, group=variable, fill=variable))+
  geom_bar(stat="identity", color="black", position=position_dodge())+
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.2, position=position_dodge(.9)) +
  scale_fill_manual(values=c( "darkred","#999999", "navy"))+
  labs(x="", y="% of all neurons")+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  scale_y_continuous(limits = c(0,1),expand = c(0,0))+
  theme(legend.title = element_blank())

setwd("~cchen2/Documents/neuroscience/Pn\ project/Figure/PDF/")
cairo_pdf("p_cat.pdf", width = 95/25.6, height = 65/25.6, family = "Arial")
p_cat
dev.off()
## within and cross cells distance-----
## function to trim cell position file with the manually confirmed cells
cc_trim_position <- function(path_trace, path_position){
  dat_trace <- read.xlsx(path_trace, colNames = F, rowNames = F)
  ## choose the valid cell position
  dat_position <- read.csv(path_position)
  cell_valid <- dat_trace[,1]
  dat_position <- dat_position[cell_valid==1,]
  return(dat_position)
}

# for m3
path_trace_m3 <- as.list(list.files(path ="~cchen2/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/matlab_analysis/m3/trace_days/", pattern = "*.xlsx", full.names = T ))
path_position_m3 <- as.list(list.files(path ="~cchen2/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/matlab_analysis/m3/cell_positon/",pattern = "*.csv", full.names = T))
dat_position_m3 <- mapply(cc_trim_position, path_trace_m3, path_position_m3, SIMPLIFY = F)

# for m7
path_trace_m7 <- as.list(list.files(path ="~cchen2/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/matlab_analysis/m7/trace_days/", pattern = "*.xlsx", full.names = T ))
path_position_m7 <- as.list(list.files(path ="~cchen2/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/matlab_analysis/m7/cell_position/",pattern = "*.csv", full.names = T))
dat_position_m7 <- mapply(cc_trim_position, path_trace_m7, path_position_m7, SIMPLIFY = F)

# for m17
path_trace_m17 <- as.list(list.files(path ="~cchen2/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/matlab_analysis/m17/trace_days/", pattern = "*.xlsx", full.names = T ))
path_position_m17 <- as.list(list.files(path ="~cchen2/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/matlab_analysis/m17/cell_position/",pattern = "*.csv", full.names = T))
dat_position_m17 <- mapply(cc_trim_position, path_trace_m17, path_position_m17, SIMPLIFY = F)

# for m18
path_trace_m18 <- as.list(list.files(path ="~cchen2/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/matlab_analysis/m18/trace_days/", pattern = "*.xlsx", full.names = T ))
path_position_m18 <- as.list(list.files(path ="~cchen2/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/matlab_analysis/m18/cell_position/",pattern = "*.csv", full.names = T))
dat_position_m18 <- mapply(cc_trim_position, path_trace_m18, path_position_m18, SIMPLIFY = F)

dat_position_com <- vector(mode = "list", 7)
for (i in 1:7){
  dat_position <- rbind(dat_position_m3[[i]], dat_position_m7[[i]], dat_position_m17[[i]], dat_position_m18[[i]])
  rownames(dat_position) <- NULL
  dat_position$ID<- dat_cell_cat[[i]]$ID
  dat_position$Group <- dat_cell_cat[[i]]$Group
  dat_position_com[[i]] <- dat_position
}

## calculate the cross cells distance
cc_dist_fun <- function(cell_position){
  ## cross cell distance
  cross_distance<- by(cell_position[,1:2], cell_position[,"ID"], function(x) as.matrix(dist(x, upper = T)), simplify = F)
  for (i in 1: length(cross_distance)){
    y <- cross_distance[[i]]
    y[y==0] <- NA
    cross_distance[[i]] = y
  }
  cross_distance_mean <- unlist(mapply(function(x) colMeans(x, na.rm = T), cross_distance, SIMPLIFY = T))
  names(cross_distance_mean) <- sub('.*\\.', '', names(cross_distance_mean))
  cross_distance_mean <- cross_distance_mean[order(as.numeric(names(cross_distance_mean)))]
  
  ## within group distance
  within_distance <- by(cell_position[,1:2], list(cell_position[,"ID"], cell_position$Group), function(x) as.matrix(dist(x)), simplify = F)
  within_distance<- Filter(Negate(is.null), within_distance)
  for (i in 1: length(within_distance)){
    x <- within_distance[[i]]
    x[x==0]<- NA
    within_distance[[i]] = x
  }
  within_distance_mean <- unlist(mapply(function(x) colMeans(x, na.rm = T), within_distance, SIMPLIFY = T))
  within_distance_mean<- within_distance_mean[order(as.numeric(names(within_distance_mean)))]
  
  ## number 2.51 from Corder science paper
  cell_position$cross_dis <- cross_distance_mean*2.51
  cell_position$within_dis <- within_distance_mean* 2.51
  return(cell_position)
}

dat_cell_distance_list <- mapply(cc_dist_fun, dat_position_com, SIMPLIFY = F)
dat_cell_distance <- do.call(rbind, dat_cell_distance_list)
rep_time <- mapply(nrow, dat_cell_distance_list)
dat_cell_distance$Day <- rep(c("Rm", "Pre", "Pre", "Rm", "Cond.", "Cond.","Test"), rep_time)
dat_cell_distance <- subset(dat_cell_distance, dat_cell_distance$Day!="Rm")
dat_cell_distance$Group <- factor(dat_cell_distance$Group, levels =c("Excited", "Neutral", "Inhibited") )
dat_cell_distance$Day <- factor(dat_cell_distance$Day, levels = c("Pre", "Cond.", "Test"))

ggplot(dat_cell_distance, aes(within_dis, cross_dis, colour=Group))+
  geom_point()+
  geom_abline(slope = 1,intercept = 0, colour="gray", linetype=2)+
  facet_grid(cols = vars(Day))+
  scale_colour_manual(values=c( "darkred", "#999999","navy"))+
  labs(x="Within-group distance (µm)", y="Across-group distance (µm)  ")+
  theme(axis.line = element_line(colour = "black"),
        axis.text.x = element_text(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  scale_x_continuous(limits = c(0, 350))+
  scale_y_continuous(limits = c(0, 350),expand = c(0,0))+
  theme(legend.title = element_blank())

dat_cell_distance_re <- melt(dat_cell_distance[,-c(1,2)], id.vars = c('Day', "ID", "Group"))
dat_cell_distance_sta <- ddply(dat_cell_distance_re, .(ID,Day,Group,variable), summarise,value=mean(value, na.rm = T))
dat_cell_distance_sta1 <- ddply(dat_cell_distance_sta, .(Day, Group,variable), summarise,n=length(value),mean=mean(value, na.rm = T),sd=sd(value, na.rm = T),se=sd(value, na.rm = T)/sqrt(length(value)))

ggplot(dat_cell_distance_sta1, aes(Day, mean, group=variable, fill=variable))+
  geom_bar(stat="identity", color="black", position=position_dodge())+
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.2, position=position_dodge(.9)) +
  facet_grid(rows = vars(Group))+
  scale_fill_manual(values=c( "black","gray"))+
  labs(x="", y="Between cell distance (µm)")+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  scale_y_continuous( limits = c(0, 150),expand = c(0,0))+
  theme(legend.title = element_blank())

## cells correlation analysis-----
cc_cor_fun <- function(cell_trace){
  dat_cell_trace<- subset(cell_trace, cell_trace$Group!="Neutral")
  ID_cell <- unique(dat_cell_trace$ID)
  dat_cell_cor <- vector(mode = "list", length(ID_cell))
  
  for (i in seq_along(ID_cell)){
    dat_cell_trace_cor <- subset(dat_cell_trace, dat_cell_trace$ID==ID_cell[i])
    dat_cell_trace_cor <- dcast(dat_cell_trace_cor[,1:3], Time~variable)
    res_cor <- cor(dat_cell_trace_cor[,-1])
    res_cor[res_cor==1]<- NA
    dat_cell_cor[[i]] <- res_cor
  }
  return(dat_cell_cor)
}

dat_cell_cor_list <- mapply(cc_cor_fun, dat_cell_trace_re, SIMPLIFY = F)

## cor matrix heatmap for m3
p_cor_d3<- ggcorrplot(dat_cell_cor_list[[3]][[1]],hc.order = TRUE, outline.col = "white", ggtheme = ggplot2::theme_void, 
                     show.legend = F, colors = c("navy", "white", "red4")) 

p_cor_d6<- ggcorrplot(dat_cell_cor_list[[6]][[1]],hc.order = TRUE, outline.col = "white", ggtheme = ggplot2::theme_void, 
                      show.legend = F, colors = c("navy", "white", "red4")) 

p_cor_d7<- ggcorrplot(dat_cell_cor_list[[7]][[1]],hc.order = TRUE, outline.col = "white", ggtheme = ggplot2::theme_void, 
                      show.legend = F, colors = c("navy", "white", "red4")) 

## compare the correlation between cells
dat_cell_cor <- c()
dat_cell_cor_sta <- c()
for (i in c(2,3,5,6,7)){
  cor_value <- unlist(dat_cell_cor_list[[i]])
  cor_value <- abs(cor_value[!is.na(cor_value)])
  dat_cor <- data.frame(Day=cell_area_day[i], value=cor_value)
  dat_cell_cor <- rbind(dat_cell_cor, dat_cor)
  
  ## for mean of cor
  cor_value_mean <- mapply(function(x) mean(colMeans(abs(x), na.rm=T)), dat_cell_cor_list[[i]])
  cor_value_max <- mapply(function(x) mean(apply(x, 2, function(x) max(abs(x), na.rm = T))), dat_cell_cor_list[[i]])

  mouse_ID <- c(1: length(dat_cell_cor_list[[i]]))
  dat_cor_mean <- data.frame(Day = cell_area_day[i], ID= mouse_ID, value=cor_value_mean, value_max=cor_value_max)
  dat_cell_cor_sta <- rbind(dat_cell_cor_sta, dat_cor_mean)
  }

## plot the trace to show correlation, cell 18, cell11, and -cell12 m3 d6
trace_cell_plot <- subset(dat_cell_trace_re[[6]], dat_cell_trace_re[[6]]$variable=="Cell18"|dat_cell_trace_re[[6]]$variable=="Cell11"|dat_cell_trace_re[[6]]$variable=="Cell12")

p_trace_em <- ggplot(trace_cell_plot, aes(Time, value, color=variable))+
  geom_line()+
  scale_color_manual(values = c("red",  "blue", "black"))+
  theme_void()+
  annotate(x=c(5,7), y=c(-1,-1),"path")+
  annotate(x=c(7,7), y=c(-1,0),"path")+
  theme(legend.position = 'none')

setwd("~cchen2/Documents/neuroscience/Pn\ project/Figure/PDF/")
cairo_pdf("p_trace_em.pdf", width = 65/25.6, height = 65/25.6, family = "Arial")
p_trace_em
dev.off()  

## cummulative plot of matrix
p_cor_cum<- ggplot(dat_cell_cor, aes(value, group=Day, colour=Day))+
  stat_ecdf(geom = "step")+
  scale_color_manual(values=c( "mediumseagreen","indianred4", "dodgerblue4"))+
  labs(x="", y="Cummulative fraction of neurons")+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  scale_y_continuous(expand = c(0,0))+
  theme(legend.title = element_blank())+
  theme(legend.position = c(0.2, 0.8))

setwd("~cchen2/Documents/neuroscience/Pn\ project/Figure/PDF/")
cairo_pdf("p_cor_cum.pdf", width = 70/25.6, height = 65/25.6, family = "Arial")
p_cor_cum
dev.off()

## bar plot
dat_cell_cor_sta <- subset(dat_cell_cor_sta, dat_cell_cor_sta$value>0.5)
dat_cell_cor_sta_re <- melt(dat_cell_cor_sta, id.vars = c("ID", "Day"))
dat_cell_cor_sta1 <- ddply(dat_cell_cor_sta_re, .(ID, Day, variable), summarise, value=mean(value))
dat_cell_cor_sta2 <- ddply(dat_cell_cor_sta1, .(Day, variable), summarise,n=length(value),mean=mean(value),sd=sd(value),se=sd(value)/sqrt(length(value)))


p_cor_mean<- ggplot(subset(dat_cell_cor_sta2, dat_cell_cor_sta2$variable=="value"), aes(Day, mean, fill=Day))+
  geom_bar(stat="identity")+
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.1, position=position_dodge(.9), colour="black") +
  scale_fill_manual(values=c( "mediumseagreen","indianred4", "dodgerblue4"))+
  labs(x="", y="Correlation coefficient (r)")+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  scale_y_continuous( limits = c(0, 1),expand = c(0,0))+
  theme(legend.title = element_blank())+
  theme(legend.position = "none")

p_cor_best <- ggplot(subset(dat_cell_cor_sta2, dat_cell_cor_sta2$variable=="value_max"), aes(Day, mean, fill=Day))+
  geom_bar(stat="identity")+
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.1, position=position_dodge(.9), colour="black") +
  scale_fill_manual(values=c( "mediumseagreen","indianred4", "dodgerblue4"))+
  labs(x="", y="Best-match correlation (r) ")+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  scale_y_continuous( limits = c(0, 1),expand = c(0,0))+
  theme(legend.title = element_blank())+
  theme(legend.position = "none")

p_cor_bar <- plot_grid(p_cor_mean, p_cor_best, nrow = 1)

setwd("~cchen2/Documents/neuroscience/Pn\ project/Figure/PDF/")
cairo_pdf("p_cor_bar.pdf", width = 82/25.6, height = 65/25.6, family = "Arial")
p_cor_bar
dev.off()


## cell covariance analysis-------
cc_cov_fun <- function(cell_trace){
  dat_cell_trace<- subset(cell_trace, cell_trace$Group!="Neutral")
  ID_cell <- unique(dat_cell_trace$ID)
  dat_cell_cov <- vector(mode = "list", length(ID_cell))
  
  for (i in seq_along(ID_cell)){
    dat_cell_trace_cov <- subset(dat_cell_trace, dat_cell_trace$ID==ID_cell[i])
    dat_cell_trace_cov <- dcast(dat_cell_trace_cov[,1:3], Time~variable)
    res_cov <- cov(dat_cell_trace_cov[,-1])
    res_cov[res_cov==1]<- NA
    dat_cell_cov[[i]] <- res_cov
  }
  return(dat_cell_cov)
}

dat_cell_cov_list <- mapply(cc_cov_fun, dat_cell_trace_re, SIMPLIFY = F)

## cov matrix heatmap for m3
cov_range <- range(c(dat_cell_cov_list[[3]][[1]], dat_cell_cov_list[[6]][[1]], dat_cell_cov_list[[7]][[1]]))

p_cov_d3<- ggcorrplot(dat_cell_cov_list[[3]][[1]],hc.order = TRUE, outline.col = "white", ggtheme = ggplot2::theme_void, 
                      show.legend = F)+
  scale_fill_gradient2(limit = c(cov_range[1], cov_range[2]), low = "navy", high =  "red4", mid = "white")

p_cov_d6<- ggcorrplot(dat_cell_cov_list[[6]][[1]],hc.order = TRUE, outline.col = "white", ggtheme = ggplot2::theme_void, 
                      show.legend = F)+
  scale_fill_gradient2(limit = c(cov_range[1], cov_range[2]), low = "navy", high =  "red4", mid = "white")

p_cov_d7<- ggcorrplot(dat_cell_cov_list[[7]][[1]],hc.order = TRUE, outline.col = "white", ggtheme = ggplot2::theme_void, 
                      show.legend = F)+
  scale_fill_gradient2(limit = c(cov_range[1], cov_range[2]), low = "navy", high =  "red4", mid = "white")

p_cov_com <- plot_grid(p_cov_d3, p_cov_d6, p_cov_d7, nrow = 1)

setwd("~cchen2/Documents/neuroscience/Pn\ project/Figure/PDF/")
cairo_pdf("p_cov_com.pdf", width = 170/25.6, height = 60/25.6, family = "Arial")
p_cov_com
dev.off()


## compare the covrelation between cells
dat_cell_cov <- c()
dat_cell_cov_sta <- c()
for (i in c(2,3,5,6,7)){
  cov_value <- unlist(dat_cell_cov_list[[i]])
  cov_value <- abs(cov_value[!is.na(cov_value)])
  dat_cov <- data.frame(Day=cell_area_day[i], value=cov_value)
  dat_cell_cov <- rbind(dat_cell_cov, dat_cov)
  
  ## for mean of cov
  cov_value_mean <- mapply(function(x) mean(colMeans(abs(x), na.rm=T)), dat_cell_cov_list[[i]])
  cov_value_max <- mapply(function(x) mean(apply(x, 2, function(x) max(abs(x), na.rm = T))), dat_cell_cov_list[[i]])
  
  mouse_ID <- c(1: length(dat_cell_cov_list[[i]]))
  dat_cov_mean <- data.frame(Day = cell_area_day[i], ID= mouse_ID, value=cov_value_mean, value_max=cov_value_max)
  dat_cell_cov_sta <- rbind(dat_cell_cov_sta, dat_cov_mean)
}


## cummulative plot of matrix
p_cov_cum<- ggplot(dat_cell_cov, aes(value, group=Day, colour=Day))+
  stat_ecdf(geom = "step")+
  scale_color_manual(values=c( "mediumseagreen","indianred4", "dodgerblue4"))+
  labs(x="", y="Cummulative fraction of neurons")+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  scale_y_continuous(expand = c(0,0))+
  theme(legend.title = element_blank())+
  theme(legend.position = c(0.8, 0.4))

setwd("~cchen2/Documents/neuroscience/Pn\ project/Figure/PDF/")
cairo_pdf("p_cov_cum.pdf", width = 70/25.6, height = 65/25.6, family = "Arial")
p_cov_cum
dev.off()

## bar plot
dat_cell_cov_sta_re <- melt(dat_cell_cov_sta, id.vars = c("ID", "Day"))
dat_cell_cov_sta1 <- ddply(dat_cell_cov_sta_re, .(ID, Day, variable), summarise, value=mean(value))
dat_cell_cov_sta2 <- ddply(dat_cell_cov_sta1, .(Day, variable), summarise,n=length(value),mean=mean(value),sd=sd(value),se=sd(value)/sqrt(length(value)))


p_cov_mean<- ggplot(subset(dat_cell_cov_sta2, dat_cell_cov_sta2$variable=="value"), aes(Day, mean, fill=Day))+
  geom_bar(stat="identity")+
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.1, position=position_dodge(.9), colour="black") +
  scale_fill_manual(values=c( "mediumseagreen","indianred4", "dodgerblue4"))+
  labs(x="", y="Max covariance (z-score)")+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  scale_y_continuous( limits = c(0,0.8),expand = c(0,0))+
  theme(legend.title = element_blank())+
  theme(legend.position = "none")

p_cov_best <- ggplot(subset(dat_cell_cov_sta2, dat_cell_cov_sta2$variable=="value_max"), aes(Day, mean, fill=Day))+
  geom_bar(stat="identity")+
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.1, position=position_dodge(.9), colour="black") +
  scale_fill_manual(values=c( "mediumseagreen","indianred4", "dodgerblue4"))+
  labs(x="", y="Max covariance (z-score)")+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  scale_y_continuous( limits = c(0, 1.5),expand = c(0,0))+
  theme(legend.title = element_blank())+
  theme(legend.position = "none")

p_cov_bar <- plot_grid(p_cov_mean, p_cov_best, nrow = 1)

setwd("~cchen2/Documents/neuroscience/Pn\ project/Figure/PDF/")
cairo_pdf("p_cov_bar.pdf", width = 82/25.6, height = 65/25.6, family = "Arial")
p_cov_bar
dev.off()
