## The script are used to analyze the data from miniscope experiemnt analyzed by matlab
## updated on 01-28-2020
## data are analyzed with inscopix data processing images
## output: 1. The trace crossing the border; 2. cells proportion increased or decrease; 3. Heat maps
## 4. Number of cells active, inhibited
## updated: 03182020, add more data
## updated: 04012020, change the way to do z-score
## updated: 04142020, do k-mean analysis firstly and back cat the cells activation
## updated: since the cells activate and inhibited randomly, stop doing average

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
library(magrittr)

cc_miniscope_matlab <- function(ID_trace, t_stim) {
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
  names(dat_stim) <- sprintf("Trail%s",seq(1:length(t_stim)))
  # dat_trace <- do.call(cbind, dat_stim)
  # ## average the trace by cell number
  # if (length(t_stim)>1){
  #   dat_cell_trace_average <- aaply(laply(dat_stim, as.matrix), c(2, 3), function(x) mean(x, na.rm=T))
  #   
  # } else {
  #   dat_cell_trace_average <- data.matrix(dat_stim[[1]])
  #   
  # }
  
  return(dat_stim)
  
}


## compare crossing during one day to proof that the random cell activation ot depressed------
path_trace_m3 <- "~cchen2/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/matlab_analysis/m3/trace_days/m3_trace_d2.xlsx"
t_stim_m3 <- c(274, 831,1364)
dat_trace_m3 <- cc_miniscope_matlab(path_trace_m3, t_stim_m3)

path_trace_m7 <- "~cchen2/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/matlab_analysis/m7/trace_days/m7_trace_d2.xlsx"
t_stim_m7 <- c(798, 177,3300)
dat_trace_m7 <- cc_miniscope_matlab(path_trace_m7, t_stim_m7)


path_trace_m17 <- "~cchen2/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/matlab_analysis/m17/trace_days/m17_trace_d2.xlsx"
t_stim_m17 <- c(1103,1678,2291)
dat_trace_m17 <- cc_miniscope_matlab(path_trace_m17, t_stim_m17)


path_trace_m18 <- "~cchen2/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/matlab_analysis/m18/trace_days/m18_trace_d2.xlsx"
t_stim_m18 <- c(236, 898,1528)
dat_trace_m18 <- cc_miniscope_matlab(path_trace_m18, t_stim_m18)

## for heat map plot
dat_cell_trace <- vector(mode = "list", 3)
for (i in 1:3) {
  dat_cell_trace[[i]]<- cbind(dat_trace_m3[[i]], dat_trace_m7[[i]], dat_trace_m17[[i]], dat_trace_m18[[i]] )
  
}

dat_cell_trace_re <- vector(mode = "list", 3)
for (i in 1:length(dat_cell_trace)) {
  dat_cell_trace_d <- dat_cell_trace[[i]]
  cell_id <- sprintf("Cell%s",seq(1:ncol(dat_cell_trace_d)))
  colnames(dat_cell_trace_d)<- cell_id
  
  dat_cell_trace_d_cluster <- unname(kmeans(t(dat_cell_trace_d), centers = 3)[[1]])
  dat_cell_trace_d<- as.data.frame(dat_cell_trace_d)
  stim_time<- seq(-2, 4.5, by=0.5)
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

## align by cells
score_range <- range(mapply(function(x) x$value, dat_cell_trace_re, SIMPLIFY = T))

## plot by aligh each cells on cross1
dat_trace <- dat_cell_trace_re[[1]]
dat_trace_sta <- ddply(dat_trace, .(variable, Group), summarise,mean=mean(value), sum=sum(value))
dat_trace_sta <- dat_trace_sta[order(dat_trace_sta[,'mean']),]
dat_trace$variable <- factor(dat_trace$variable, levels = dat_trace_sta$variable)

for (i in c(1:3)) {
  dat_trace <- dat_cell_trace_re[[i]]
  dat_trace$variable <- factor(dat_trace$variable, levels = dat_trace_sta$variable)
  p_heat <- ggplot(dat_trace, aes(Time, variable,fill= value))+ 
    geom_tile(height=2)+
    # facet_grid(rows = vars(Group), scales = "free_y")+
    #scale_fill_gradient2(limits= score_range,low = "navy", high = "red4", mid = "white")+
    scale_fill_gradientn(limits= score_range, colours = c("navy", "white", "red3"), values = rescale(c(score_range[1], 0, score_range[2])))+
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
  assign(paste0("p_heat_cross", i), p_heat)
}

p_heat_cross <- plot_grid(p_heat_cross1, p_heat_cross2,p_heat_cross3, nrow = 1)

setwd("~cchen2/Documents/neuroscience/Pn\ project/Figure/PDF/")
cairo_pdf("p_anti_cross.pdf", width = 120/25.6, height = 65/25.6, family = "Arial")
p_heat_cross
dev.off()



## compare pinprick to show the random activation of cells-----
path_m3_pin <- "~cchen2/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/matlab_analysis/m3/pin_prick/m3_trace_pinprick.xlsx"
t_pin_m3 <- c(337,665,916)

dat_pin_m3 <-cc_miniscope_matlab(path_m3_pin, t_pin_m3)
## for m7
path_m7_pin <- "~cchen2/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/matlab_analysis/m7/pin_prick/m7_trace_pinprick.xlsx"
t_pin_m7 <- c(493, 1323, 1867) 
dat_pin_m7 <- cc_miniscope_matlab (path_m7_pin, t_pin_m7)

dat_cell_trace <- vector(mode = "list", 3)
for (i in 1:3) {
  dat_cell_trace[[i]]<- cbind(dat_pin_m3[[i]], dat_pin_m7[[i]] )
  
}

dat_cell_trace_re <- vector(mode = "list", 3)
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
  rep_time <- c(ncol(dat_pin_m3[[i]]), ncol(dat_pin_m7[[i]]))
  dat_cell_trace_d_re$ID <- rep(c("m3", "m7"), rep_time*length(stim_time))
  
  dat_cell_trace_re[[i]] <- dat_cell_trace_d_re
}

## align by cells
score_range <- range(mapply(function(x) x$value, dat_cell_trace_re, SIMPLIFY = T))

## plot by aligh each cells on cross1
dat_trace <- dat_cell_trace_re[[1]]
dat_trace_sta <- ddply(dat_trace, .(variable, Group), summarise,mean=mean(value), sum=sum(value))
dat_trace_sta <- dat_trace_sta[order(dat_trace_sta[,'mean']),]
dat_trace$variable <- factor(dat_trace$variable, levels = dat_trace_sta$variable)

for (i in c(1:3)) {
  dat_trace <- dat_cell_trace_re[[i]]
  dat_trace$variable <- factor(dat_trace$variable, levels = dat_trace_sta$variable)
  p_heat <- ggplot(dat_trace, aes(Time, variable,fill= value))+ 
    geom_tile(height=2)+
    # facet_grid(rows = vars(Group), scales = "free_y")+
    #scale_fill_gradient2(limits= score_range,low = "navy", high = "red4", mid = "white")+
    scale_fill_gradientn(limits= score_range, colours = c("navy", "white", "red3"), values = rescale(c(score_range[1], 0, score_range[2])))+
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
  assign(paste0("p_heat_pin", i), p_heat)
}

p_heat_pin <- plot_grid(p_heat_pin1, p_heat_pin2,p_heat_pin3, nrow = 1)


setwd("~cchen2/Documents/neuroscience/Pn\ project/Figure/PDF/")
cairo_pdf("p_anti_pin.pdf", width = 120/25.6, height = 65/25.6, family = "Arial")
p_heat_pin
dev.off()


dat_cell_trace_stim1_sta <- ddply(dat_cell_trace_re[[1]], .(Time, Group), summarise,n=length(value),mean=mean(value),sd=sd(value),se=sd(value)/sqrt(length(value)))
dat_cell_trace_stim2_sta <- ddply(dat_cell_trace_re[[2]], .(Time, Group), summarise,n=length(value),mean=mean(value),sd=sd(value),se=sd(value)/sqrt(length(value)))
dat_cell_trace_stim3_sta <- ddply(dat_cell_trace_re[[3]], .(Time, Group), summarise,n=length(value),mean=mean(value),sd=sd(value),se=sd(value)/sqrt(length(value)))

dat_cell_trace_sta <- rbind(dat_cell_trace_stim1_sta,  dat_cell_trace_stim2_sta, dat_cell_trace_stim3_sta)

dat_cell_trace_sta$Day <- rep(c("stim1",  "stim2", "stim3"),each= length(stim_time)*3)
dat_cell_trace_sta$Day <- factor(dat_cell_trace_sta$Day, levels = c("stim1",  "stim2", "stim3"))
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
  #scale_y_continuous(limits = c(range_trace_plot[1]- 0.2, range_trace_plot[2]+0.2))+
  theme(legend.title = element_blank())

dat_cell_area <- c()
cell_area_day <- c("stim1",  "stim2", "stim3")
for (i in c(1:3)){
  rep_time <- c(ncol(dat_pin_m3[[i]]), ncol(dat_pin_m7[[i]]))
  dat_trace <- dat_cell_trace_re[[i]]
  dat_trace$ID <- rep(c("m3", "m7"), rep_time*length(stim_time))
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

dat_cell_area$Day <- factor(dat_cell_area$Day, levels = c("stim1","stim2","stim3"))
dat_cell_area_com <- ddply(dat_cell_area, .(ID, Day), summarise, value=mean(ratio, na.rm = T))
dat_cell_area_sta <- ddply(dat_cell_area_com,.(Day), summarise,n=length(value),mean=mean(value),sd=sd(value),se=sd(value)/sqrt(length(value)))


p_EI_ratio <- ggplot(dat_cell_area_sta, aes(Day, mean, colour=Day, fill=Day))+
  geom_bar(stat="identity", color="black", position=position_dodge())+
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.1, position=position_dodge(.9), colour="black") +
  # scale_fill_manual(values=c( "mediumseagreen","indianred4", "dodgerblue4"))+
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


## for heat pain to proof random activation-----
path_m3_har <- "~cchen2/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/matlab_analysis/Hargreave/m3_trace_har.xlsx"
t_har_m3 <- c(1171, 2159, 3579)

dat_har_m3 <-c_miniscope_matlab(path_m3_har, t_har_m3)

## for m7
path_m7_har <- "~cchen2/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/matlab_analysis/Hargreave/m7_trace_har.xlsx"
t_har_m7 <- c(756, 1917, 2947)

dat_har_m7 <-c_miniscope_matlab(path_m7_har, t_har_m7)


dat_cell_trace <- vector(mode = "list", 3)
for (i in 1:3) {
  dat_cell_trace[[i]]<- cbind(dat_har_m3[[i]], dat_har_m7[[i]] )
  
}

dat_cell_trace_re <- vector(mode = "list", 3)
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
  rep_time <- c(ncol(dat_har_m3[[i]]), ncol(dat_har_m7[[i]]))
  dat_cell_trace_d_re$ID <- rep(c("m3", "m7"), rep_time*length(stim_time))
  
  dat_cell_trace_re[[i]] <- dat_cell_trace_d_re
}

## align by cells
score_range <- range(mapply(function(x) x$value, dat_cell_trace_re, SIMPLIFY = T))

## plot by aligh each cells on cross1
dat_trace <- dat_cell_trace_re[[1]]
dat_trace_sta <- ddply(dat_trace, .(variable, Group), summarise,mean=mean(value), sum=sum(value))
dat_trace_sta <- dat_trace_sta[order(dat_trace_sta[,'mean']),]
dat_trace$variable <- factor(dat_trace$variable, levels = dat_trace_sta$variable)

for (i in c(1:3)) {
  dat_trace <- dat_cell_trace_re[[i]]
  dat_trace$variable <- factor(dat_trace$variable, levels = dat_trace_sta$variable)
  p_heat <- ggplot(dat_trace, aes(Time, variable,fill= value))+ 
    geom_tile(height=2)+
    # facet_grid(rows = vars(Group), scales = "free_y")+
    #scale_fill_gradient2(limits= score_range,low = "navy", high = "red4", mid = "white")+
    scale_fill_gradientn(limits= score_range, colours = c("navy", "white", "red3"), values = rescale(c(score_range[1], 0, score_range[2])))+
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
  assign(paste0("p_heat_har", i), p_heat)
}

p_heat_har <- plot_grid(p_heat_har1, p_heat_har2,p_heat_har3, nrow = 1)

setwd("~cchen2/Documents/neuroscience/Pn\ project/Figure/PDF/")
cairo_pdf("p_anti_har.pdf", width = 120/25.6, height = 65/25.6, family = "Arial")
p_heat_har
dev.off()

## calculate the spontaneous firing frequency of the cells-------
#(take the day3 at this moment, may need to get data when they in homecage)

cc_firing_fun <- function(path_trace, path_peak){
  dat_trace <- read.xlsx(path_trace, colNames = F, rowNames = F)
  cell_valid <- dat_trace$X1
  
  # total time
  n_frame <- ncol(dat_trace)-1
  # 30 is the acquiring frequency
  t_period <- n_frame/30 # s
  dat_peak <- read.csv(path_peak)
  
  ## choose the valide cell number
  dat_peak <- dat_peak[, cell_valid==1]
  
  cell_rate <- apply(dat_peak, 2, function(x) length(which(!is.na(x)))/t_period)
  cell_rate_mean <- mean(cell_rate)
  cell_rate_max <- max(cell_rate)
  return(c(cell_rate_mean, cell_rate_max))
}

## for m3
path_trace_m3 <- "~cchen2/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/matlab_analysis/m3/trace_days/m3_trace_d3.xlsx"
path_peak_m3 <- "~cchen2/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/matlab_analysis/m3/peaks_days/m3_peaks_d3.csv"

dat_rate_m3 <- cc_firing_fun(path_trace_m3, path_peak_m3)

## for m7
path_trace_m7 <- "~cchen2/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/matlab_analysis/m7/trace_days/m7_trace_d3.xlsx"
path_peak_m7 <- "~cchen2/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/matlab_analysis/m7/peaks_days/m7_peaks_d3.csv"

dat_rate_m7 <- cc_firing_fun(path_trace_m7, path_peak_m7)

## for m17
path_trace_m17 <- "~cchen2/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/matlab_analysis/m17/trace_days/m17_trace_d3.xlsx"
path_peak_m17 <- "~cchen2/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/matlab_analysis/m17/peaks_days/m17_peaks_d3.csv"

dat_rate_m17 <- cc_firing_fun(path_trace_m17, path_peak_m17)

## for m18
path_trace_m18 <- "~cchen2/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/matlab_analysis/m18/trace_days/m18_trace_d3.xlsx"
path_peak_m18 <- "~cchen2/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/matlab_analysis/m18/peaks_days/m18_peaks_d3.csv"

dat_rate_m18 <- cc_firing_fun(path_trace_m18, path_peak_m18)


## combine data and plot
dat_rate <- as.data.frame(rbind(dat_rate_m3, dat_rate_m7, dat_rate_m17, dat_rate_m18))
colnames(dat_rate) <- c("rate_mean", "rate_max")
dat_rate$ID <- c("m3","m7","m17", "m18")
dat_rate$Group <- "Ctrl"

p_rate_mean <- ggplot(dat_rate, aes(Group, rate_mean))+
  geom_boxplot()+
  geom_jitter(width = 0.2, shape=1)+
  labs(x="", y="Mean firing rate (Hz)")+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  theme(legend.title = element_blank())


p_rate_max <- ggplot(dat_rate, aes(Group, rate_max))+
  geom_boxplot()+
  geom_jitter(width = 0.2, shape=1)+
  labs(x="", y="Max firing rate (Hz)")+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  theme(legend.title = element_blank())

p_rate <- plot_grid(p_rate_mean, p_rate_max, nrow = 1)

setwd("~cchen2/Documents/neuroscience/Pn\ project/Figure/PDF/")
cairo_pdf("p_firing_rate.pdf", width = 60/25.6, height = 60/25.6, family = "Arial")
p_rate
dev.off()


## randomly select time point to show spontaneous activation and deactivation-----
t_stim_m3 <- seq(300, 4000, by=500) 

t_stim_m7 <- seq(300, 4000, by=500) 

t_stim_m17 <- seq(300, 4000, by=500) 

t_stim_m18 <- seq(300, 4000, by=500) 

path_trace_list <- list(path_trace_m3, path_trace_m7, path_trace_m17, path_trace_m18)
t_stim_list <- list(m3= t_stim_m3, m7 = t_stim_m7, m17=t_stim_m17, m18=t_stim_m18)


dat_cell_trace <- mapply(cc_miniscope_matlab, path_trace_list, t_stim_list, SIMPLIFY = F)
dat_cell_trace_trail <- vector(mode = "list", length = length(t_stim_m3))
for (i in c(1: length(t_stim_m3))){
  dat_cell_trail<- do.call(cbind,sapply(dat_cell_trace,'[[',i))
  cell_id <- sprintf("Cell%s",seq(1:ncol(dat_cell_trail)))
  colnames(dat_cell_trail)<- cell_id
  
  dat_cell_trail_cluster <- unname(kmeans(t(dat_cell_trail), centers = 3)[[1]])
  dat_cell_trail<- as.data.frame(dat_cell_trail)
  stim_time <- seq(-5, 8.5, by=0.5)
  dat_cell_trail$Time <- stim_time
  dat_cell_trail_re <- melt(dat_cell_trail, id.vars ='Time')
  dat_cell_trail_re$Group <- rep(dat_cell_trail_cluster, each=length(stim_time))
  
  ## sort data by the value in each group
  dat_cell_d_sort <- as.numeric(names(sort(tapply(dat_cell_trail_re$value, dat_cell_trail_re$Group, mean), decreasing = T)))
  
  dat_cell_trail_re$Group[dat_cell_trail_re$Group == dat_cell_d_sort[1]] ="Excited"
  dat_cell_trail_re$Group[dat_cell_trail_re$Group == dat_cell_d_sort[2]] ="Neutral"
  dat_cell_trail_re$Group[dat_cell_trail_re$Group == dat_cell_d_sort[3]] ="Inhibited"
  
  dat_cell_trail_re$Group <- factor(dat_cell_trail_re$Group, levels = c("Excited", "Neutral", "Inhibited"))
  rep_time <- c(ncol(dat_cell_trace[[1]][[1]]), ncol(dat_cell_trace[[2]][[1]]), ncol(dat_cell_trace[[3]][[1]]), ncol(dat_cell_trace[[4]][[1]]))
  dat_cell_trail_re$ID <- rep(c("m3", "m7", "m17", "m18"), rep_time*length(stim_time))
  
  dat_cell_trace_trail[[i]] <- dat_cell_trail_re
}

## mean by trails
sta_fun <- function(dat_trace_re){
  dat_trace_re_sta <- ddply(dat_trace_re, .(Time, Group), summarise,n=length(value),value1=mean(value),sd=sd(value),se=sd(value)/sqrt(length(value)))
}
dat_cell_trace_sta <- do.call(rbind, mapply(sta_fun, dat_cell_trace_trail, SIMPLIFY = F))
dat_cell_trace_sta1 <- ddply(dat_cell_trace_sta, .(Time, Group),summarise,n=length(value1),mean=mean(value1),sd=sd(value1),se=sd(value1)/sqrt(length(value1)))
dat_cell_trace_sta1$Group <- factor(dat_cell_trace_sta1$Group, levels = c("Neutral", "Excited", "Inhibited"))
range_trace_plot <- range(dat_cell_trace_sta1$mean)

p_trace <- ggplot(dat_cell_trace_sta1, aes(Time, mean, colour=Group))+
  geom_line()+
  geom_ribbon(aes(ymin=mean-se, ymax=mean+se, fill=Group), alpha=0.1, linetype=0)+
  scale_colour_manual(values=c("#999999", "darkred", "navy"))+
  scale_fill_manual(values=c("#999999", "darkred", "navy"))+
  labs(x="Time relative to crossing (s)", y="Z score")+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  geom_vline(xintercept = 0, col= 'red', linetype=2)+
  # scale_y_continuous(limits = c(range_trace_plot[1]- 0.2, range_trace_plot[2]+0.3))+
  theme(legend.title = element_blank(), legend.position = c(0.2, 0.75))

## calculate the proportion of cell cat
cat_fun <- function(dat_trace_re){
  dat_cell_trace_ca <- subset(dat_trace_re, dat_trace_re$Time==0)
  dat_cat <- prop.table(table(dat_cell_trace_ca$ID, dat_cell_trace_ca$Group), 1)
  dat_cat <- as.data.frame(dat_cat)
}

dat_cat <- do.call(rbind, mapply(cat_fun, dat_cell_trace_trail, SIMPLIFY = F))
dat_cat_sta <- ddply(dat_cat, .(Var1, Var2), summarise, value=mean(Freq))
dat_cat_sta1 <- ddply(dat_cat_sta, .(Var2), summarise, n=length(value),mean=mean(value),sd=sd(value),se=sd(value)/sqrt(length(value)))
dat_cat_sta1$Group <- "Ctrl"
dat_cat_sta1$Var2 <- factor(dat_cat_sta1$Var2, levels = c("Inhibited", "Neutral" ,"Excited"))

p_cat_bar <- ggplot(dat_cat_sta1, aes(x=Group, y=mean, fill=Var2)) +
  geom_bar(stat="identity")+
  scale_fill_manual(values=c( "navy", "#999999","darkred"))+
  labs(x="", y="Fraction of all cells")+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  theme(legend.title = element_blank())

score_range <- range(dat_cell_trace_d_re$value)

dat_trace_sta <- ddply(dat_cell_trace_d_re, .(variable, Group), summarise,mean=mean(value), sum=sum(value))
dat_trace_sta <- dat_trace_sta[order(dat_trace_sta[,'mean']),]
dat_cell_trace_d_re$variable <- factor(dat_cell_trace_d_re$variable, levels = dat_trace_sta$variable)
p_heat <- ggplot(dat_cell_trace_d_re, aes(Time, variable,fill= value))+ 
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
          axis.title=element_text(family = "Arial",size = 12, face ="plain"))


setwd("~cchen2/Documents/neuroscience/Pn\ project/Figure/PDF/")
cairo_pdf("p_heat_ctrl.pdf", width = 60/25.6, height = 65/25.6, family = "Arial")
p_heat
dev.off()
## plot trace by group
dat_cell_trace_sta <- ddply(dat_cell_trace_d_re, .(Time, Group), summarise,n=length(value),mean=mean(value),sd=sd(value),se=sd(value)/sqrt(length(value)))
  
dat_cell_trace_sta$Group <- factor(dat_cell_trace_sta$Group, levels = c("Neutral", "Excited", "Inhibited"))
range_trace_plot <- range(dat_cell_trace_sta$mean)

## calculate the proportion of cell types
dat_cell_trace_ca <- subset(dat_cell_trace_d_re, dat_cell_trace_d_re$Time==0)
dat_cat <- as.data.frame(prop.table(table(dat_cell_trace_ca$Group)))
colnames(dat_cat)[1] <- "Cat"
dat_cat$Cat <- factor(dat_cat$Cat, levels = c("Inhibited", "Neutral" ,"Excited"))

dat_cat$Group <- "Ctrl"
dat_cat$label_y <- cumsum(dat_cat$Freq)
dat_cat$Freq <- round(dat_cat$Freq, 2)
p_cat_bar <- ggplot(dat_cat, aes(x=Group, y=Freq, fill=Cat)) +
  geom_bar(stat="identity")+
  geom_text(aes(y=label_y, label=Freq), vjust=1.6, color="white", size=3.5)+
  scale_fill_manual(values=c( "navy", "#999999","darkred"))+
  labs(x="", y="Fraction of all cells")+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  theme(legend.title = element_blank())

setwd("~cchen2/Documents/neuroscience/Pn\ project/Figure/PDF/")
cairo_pdf("p_bar_cat.pdf", width = 50/25.6, height = 65/25.6, family = "Arial")
p_cat_bar
dev.off()

## group by day
p_trace <- ggplot(dat_cell_trace_sta, aes(Time, mean, colour=Group))+
    geom_line()+
    geom_ribbon(aes(ymin=mean-se, ymax=mean+se, fill=Group), alpha=0.1, linetype=0)+
    scale_colour_manual(values=c("#999999", "darkred", "navy"))+
    scale_fill_manual(values=c("#999999", "darkred", "navy"))+
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
    theme(legend.title = element_blank(), legend.position = c(0.2, 0.75))

setwd("~cchen2/Documents/neuroscience/Pn\ project/Figure/PDF/")
cairo_pdf("p_trace_ctrl.pdf", width = 70/25.6, height = 65/25.6, family = "Arial")
p_trace
dev.off()
  
## calculate the E-I ratio-----
  dat_trace_sta <- ddply(dat_cell_trace_d_re, .(ID, Time, Group), summarise,n=length(value),mean=mean(value),sd=sd(value),se=sd(value)/sqrt(length(value)))
  Time_area <- which(stim_time>=0)
  dat_anti_area <- as.data.frame(tapply(dat_trace_sta$mean, INDEX = list(dat_trace_sta$ID, dat_trace_sta$Group), 
                                        function (x) AUC(stim_time[Time_area], x[Time_area])))
  dat_anti_area$ID <- rownames(dat_anti_area)
  rownames(dat_anti_area)<-NULL
  dat_anti_area$ratio <- abs(dat_anti_area$Excited/dat_anti_area$Inhibited)
  dat_anti_area$Group <- "Ctrl"
  
  ctrl_cor <- cor(dat_anti_area[,c(1,3)])
  p_EI_cor <- ggplot(dat_anti_area, aes(Excited, Inhibited))+
    geom_point()+
    geom_smooth(method = "lm", colour="red")+
    labs(x="AUC of excitation (z-score^2)", y="AUC of inhibiton (z-score^2)")+
    theme(axis.line.x = element_line(),
          axis.line.y = element_line(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank(),
          axis.title=element_text(family = "Arial",size = 12, face ="plain"))
  
  setwd("~cchen2/Documents/neuroscience/Pn\ project/Figure/PDF/")
  cairo_pdf("p_EI_cor.pdf", width = 65/25.6, height = 65/25.6, family = "Arial")
  p_EI_cor
  dev.off()
  
  p_EI_ratio <- ggplot(dat_anti_area, aes(Group, ratio))+
    geom_boxplot()+
    geom_jitter(width = 0.2, shape=1)+
    labs(x="", y="E/I ratio")+
    theme(axis.line.x = element_line(),
          axis.line.y = element_line(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank(),
          axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
    scale_y_continuous(limits = c(1.0, 1.4), expand = c(0,0))+
    theme(legend.title = element_blank())+
    theme(legend.position = "none")

  setwd("~cchen2/Documents/neuroscience/Pn\ project/Figure/PDF/")
  cairo_pdf("p_bar_EI.pdf", width = 30/25.6, height = 65/25.6, family = "Arial")
  p_EI_ratio
  dev.off()  

## calculate the firing frequency-----
cc_pin_peak <- function(path_ID, t_stim, path_peak){
    cell_valid <- read.xlsx(path_ID, colNames = F, rowNames = F)[,1]
    dat_peak <- read.csv(path_peak, header = T) %>%
      .[, which(cell_valid==1)] 
    dat_peak_vector <- c(dat_peak[!is.na(dat_peak)])
    
    ctrl_firing_rate <- 
      c(mapply(function (x) (x-100):(x-60), as.list(t_stim))) %>%
      intersect(., dat_peak_vector) %>%
      length(.)/ncol(dat_peak)/2
    
    pin_firing_rate <- 
      c(mapply(function (x) x:(x+40), as.list(t_stim))) %>%
      intersect(., dat_peak_vector) %>%
      length(.)/ncol(dat_peak)/2
    return(c(ctrl_firing_rate, pin_firing_rate))
  }
## caculate the firing freq
  ## for m3
  path_trace_m3 <- as.list(list.files(path ="~cchen2/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/matlab_analysis/m3/trace_days/", pattern = "*.xlsx", full.names = T ))
  t_stim_m3_d1 <- c(214, 698,1204,1757,2964) 
  t_stim_m3_d2 <- c(274, 831,1364,1788,2621,3273)
  t_stim_m3_d3 <- c(1885, 2516) 
  t_stim_m3_d4 <- c(1813,3587)
  t_stim_m3_d5 <- c(758) 
  t_stim_m3_d6 <- c(132,507) 
  t_stim_m3_d7 <- c(346) ## use first two 346, 1206, 1756, 2288, 2970, 3262, 3942
  
  t_stim_m3 <- list(t_stim_m3_d1, t_stim_m3_d2, t_stim_m3_d3, t_stim_m3_d4, t_stim_m3_d5, t_stim_m3_d6, t_stim_m3_d7)
  
  path_peak_m3 <- as.list(list.files(path ="~cchen2/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/matlab_analysis/m3/peaks_days/", pattern = "*.csv", full.names = T ))
  
  dat_freq_m3 <- mapply(cc_pin_peak, path_trace_m3, t_stim_m3, path_peak_m3, SIMPLIFY = F) %>%
    do.call(rbind, .)
  
  ## for m7
  path_trace_m7 <- as.list(list.files(path ="~cchen2/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/matlab_analysis/m7/trace_days/", pattern = "*.xlsx", full.names = T ))
  t_stim_m7_d1 <- c( 605, 964, 1603, 2153, 2363, 2751, 3152) 
  t_stim_m7_d2 <- c(798, 177,3300)
  t_stim_m7_d3 <- c(360, 3352, 4044) 
  t_stim_m7_d4 <- c(3236)
  t_stim_m7_d5 <- c(317) 
  t_stim_m7_d6 <- c(205)
  t_stim_m7_d7 <- c(236) ## use first two 236, 1914,2777,3334,4091
  
  t_stim_m7 <- list(t_stim_m7_d1, t_stim_m7_d2, t_stim_m7_d3, t_stim_m7_d4, t_stim_m7_d5, t_stim_m7_d6, t_stim_m7_d7)
  path_peak_m7 <- as.list(list.files(path ="~cchen2/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/matlab_analysis/m7/peaks_days/", pattern = "*.csv", full.names = T ))
  
  dat_freq_m7 <- mapply(cc_pin_peak, path_trace_m7, t_stim_m7, path_peak_m7, SIMPLIFY = F) %>%
    do.call(rbind, .)
  
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
  path_peak_m17 <- as.list(list.files(path ="~cchen2/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/matlab_analysis/m17/peaks_days/", pattern = "*.csv", full.names = T ))
  
  dat_freq_m17 <- mapply(cc_pin_peak, path_trace_m17, t_stim_m17, path_peak_m17, SIMPLIFY = F) %>%
    do.call(rbind, .)
  
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
  path_peak_m18 <- as.list(list.files(path ="~cchen2/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/matlab_analysis/m18/peaks_days/", pattern = "*.csv", full.names = T ))
  
  dat_freq_m18 <- mapply(cc_pin_peak, path_trace_m18, t_stim_m18, path_peak_m18, SIMPLIFY = F) %>%
    do.call(rbind, .)