## only compare the pre and cond condition


## import needed library-----
library("reshape2")
library(plyr)
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
library(zoo)
library(pracma)
library(gtools)

c_miniscope_matlab <- function(ID_trace, t_stim) {
  ## import and format the data
  dat_trace <- read.xlsx(ID_trace, colNames = F, rowNames = F) %>% 
    filter(., X1 ==1) %>% 
    select(., -c(1,2)) %>% 
    t() %>% 
    as.data.frame()
  
  rownames(dat_trace) <- NULL
  colnames(dat_trace)<- NULL
  ## do z score of the whole trace
  dat_trace <- apply(dat_trace, 2, scale)
  
  ## divide ca2+ trace within -2 and 6s before and after stimuls, take -2 to 0 as baseline
  dat_stim <- vector(mode = "list", length = length(t_stim))
  stim_time<- seq(-2, 6.5, by=0.5)
  ## number of rows to be binned
  n <- 10 # 0.05*10=0.5
  
  ## extract cell activity when they cross the border
  for (i in seq_along(t_stim)){
    t1_p <- t_stim[i]
    dat_stim1 <- dat_trace[(t1_p-40):(t1_p+140-1),] 
    dat_stim1 <- aggregate(dat_stim1,list(rep(1:(nrow(dat_stim1)%/%n+1),each=n,len=nrow(dat_stim1))),mean)[-1]
    
    dat_stim[[i]] <- dat_stim1[1:5,] %>%  ## baseline as -2 to 0
      colMeans(., na.rm = T) %>% 
      sweep(dat_stim1, 2, ., FUN = "-")
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
path_trace_m3 <- as.list(list.files(path ="~cchen/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/matlab_analysis/m3/trace_days/", pattern = "*.xlsx", full.names = T ))[c(3,6)]
t_stim_m3 <- list(t_stim_m3_d3 = c(1885),t_stim_m3_d6 = c(132))
dat_trace_m3 <- mapply(c_miniscope_matlab, path_trace_m3, t_stim_m3, SIMPLIFY = F)

# for m7
path_trace_m7 <- as.list(list.files(path ="~cchen/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/matlab_analysis/m7/trace_days/", pattern = "*.xlsx", full.names = T ))[c(3,6)]
t_stim_m7 <- list(t_stim_m7_d3 = c(360), t_stim_m7_d6 = c(205))
dat_trace_m7 <- mapply(c_miniscope_matlab, path_trace_m7, t_stim_m7,SIMPLIFY = F)

# for m17
path_trace_m17 <- as.list(list.files(path ="~cchen/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/matlab_analysis/m17/trace_days/", pattern = "*.xlsx", full.names = T ))[c(3,6)]
t_stim_m17 <- list(t_stim_m17_d3 = c(437), t_stim_m17_d6 = c(617) )
dat_trace_m17 <- mapply(c_miniscope_matlab, path_trace_m17, t_stim_m17, SIMPLIFY = F)

# for m18
path_trace_m18 <- as.list(list.files(path ="~cchen/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/matlab_analysis/m18/trace_days/", pattern = "*.xlsx", full.names = T ))[c(3,6)]
t_stim_m18 <- list(t_stim_m18_d3 = c(784), t_stim_m18_d6 = c(467))
dat_trace_m18 <- mapply(c_miniscope_matlab, path_trace_m18, t_stim_m18, SIMPLIFY = F)

# for m855
path_trace_m855 <- as.list(list.files(path ="~cchen/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/matlab_analysis/m855/trace_days", pattern = "*.xlsx", full.names = T ))[c(3,6)]
t_stim_m855 <- list(t_stim_m855_d3 = c(392) *2, t_stim_m855_d6 = c(64) *2)
dat_trace_m855 <- mapply(c_miniscope_matlab, path_trace_m855, t_stim_m855, SIMPLIFY = F)

# for m857
path_trace_m857 <- as.list(list.files(path ="~cchen/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/matlab_analysis/m857/trace_days", pattern = "*.xlsx", full.names = T ))[c(1,4)]
t_stim_m857 <- list(t_stim_m857_d3 = c(847) *2, t_stim_m857_d6 = c(723) *2)
dat_trace_m857 <- mapply(c_miniscope_matlab, path_trace_m857, t_stim_m857, SIMPLIFY = F)


## combine data and do k-means analysis-----

dat_cell_trace <- vector(mode = "list", 2)
for (i in 1:2) {
  dat_cell_trace[[i]] <- cbind(dat_trace_m3[[i]], dat_trace_m7[[i]], dat_trace_m17[[i]], dat_trace_m18[[i]], dat_trace_m855[[i]], dat_trace_m857[[i]])
  
}

## overview of the data
mouse_ID <- c("m3", "m7", "m17", "m18", "m855", "m857")


dat_cell_trace_re <- vector(mode = "list", 2)
for (i in 1:length(dat_cell_trace)) {
  dat_cell_trace_d <- dat_cell_trace[[i]]
  cell_id <- sprintf("Cell%s",seq(1:ncol(dat_cell_trace_d)))
  colnames(dat_cell_trace_d)<- cell_id
  
  dat_cell_trace_d_cluster <- unname(kmeans(t(dat_cell_trace_d), centers = 3, nstart = 50)[[1]])
  stim_time<- seq(-2, 6.5, by=0.5)
  
  dat_cell_trace_d_re<- as.data.frame(dat_cell_trace_d) %>% 
    mutate(., Time= stim_time) %>% 
    melt(., id.vars='Time') %>% 
    mutate(., Group= rep(dat_cell_trace_d_cluster, each=length(stim_time)) )
  
  ## sort data by the value in each group
  dat_cell_d_sort <- as.numeric(names(sort(tapply(dat_cell_trace_d_re$value, dat_cell_trace_d_re$Group, mean), decreasing = T)))
  
  dat_cell_trace_d_re$Group[dat_cell_trace_d_re$Group == dat_cell_d_sort[1]] ="Excited"
  dat_cell_trace_d_re$Group[dat_cell_trace_d_re$Group == dat_cell_d_sort[2]] ="Neutral"
  dat_cell_trace_d_re$Group[dat_cell_trace_d_re$Group == dat_cell_d_sort[3]] ="Inhibited"
  
  rep_time <- c(ncol(dat_trace_m3[[i]]), ncol(dat_trace_m7[[i]]), ncol(dat_trace_m17[[i]]), ncol(dat_trace_m18[[i]]), ncol(dat_trace_m855[[i]]), ncol(dat_trace_m857[[i]]))
  dat_cell_trace_d_re <- mutate(dat_cell_trace_d_re, Group= factor(Group, levels = c("Excited", "Neutral", "Inhibited"))) %>% 
    mutate(., ID = rep(mouse_ID, rep_time*length(stim_time)) )
  
  dat_cell_trace_re[[i]] <- dat_cell_trace_d_re
}




heat_m_ID <- "m855"
score_rang1 <- dat_cell_trace_re[[1]] %>% 
  filter(ID == heat_m_ID)

score_range <- dat_cell_trace_re[[2]] %>% 
  filter(ID ==heat_m_ID) %>% 
  rbind(., score_rang1) %>% 
  .$value %>% 
  range()

group_day <- c("Pre", "Cond")

for (i in c(1: 2)) {
  dat_trace <- dat_cell_trace_re[[i]] %>% 
    filter(ID == heat_m_ID) 
  
  dat_trace_sta <- dat_trace %>% 
    ddply(., .(variable, Group), summarise,mean=mean(value), sum=sum(value))
  dat_trace_sta <- dat_trace_sta[order(dat_trace_sta[,'mean']),]
  dat_trace$variable <- factor(dat_trace$variable, levels = dat_trace_sta$variable)
  p_heat <- ggplot(dat_trace, aes(Time, variable,fill= value))+ 
    geom_tile(height=2)+
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
  assign(str_c("p_heat_", group_day[i]), p_heat)
}

## combine the heat plot
p_heat_com <- plot_grid(p_heat_Pre, p_heat_Cond, nrow = 1)

setwd("~cchen/Documents/neuroscience/Pn\ project/Figure/PDF/")
cairo_pdf("p_anti_heat.pdf", width = 80/25.6, height = 55/25.6, family = "Arial")
p_heat_com
dev.off()


dat_cell_trace_pre_sta_m855 <-  ddply(dat_cell_trace_re[[1]], .(ID,Time, Group), summarise,mean_value=mean(value), sd=sd(value),se=sd(value)/sqrt(length(value))) %>%
  filter(ID == "m855") %>% 
  mutate(Day ="Pre")

dat_cell_trace_test_sta_m855 <- ddply(dat_cell_trace_re[[2]], .(ID,Time, Group), summarise,mean_value=mean(value), sd=sd(value),se=sd(value)/sqrt(length(value))) %>%
  filter(ID == "m855") %>% 
  mutate(Day ="Cond")


p_trace_m855 <- rbind(dat_cell_trace_pre_sta_m855, dat_cell_trace_test_sta_m855) %>% 
  mutate(Day = factor(Day, levels = c("Pre",  "Cond")), Group = factor(Group,levels = c( "Excited","Neutral", "Inhibited") )) %>% 
  ggplot(., aes(x = Time, y = mean_value, group = Group,colour= Group))+
  geom_line()+
  geom_ribbon(aes(ymin=mean_value-se, ymax=mean_value+se, fill=Group), alpha=0.1, linetype=0)+
  facet_grid(cols = vars(Day))+
  scale_colour_manual(values=c( "darkred","#999999", "navy"))+
  scale_fill_manual(values=c( "darkred", "#999999","navy"))+
  labs(x="Time relative to crossing (s)", y="Z score")+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  geom_vline(xintercept = 0, col= 'red', linetype=2)+
  theme(legend.title = element_blank(), legend.position = 'none')

setwd("~cchen/Documents/neuroscience/Pn\ project/Figure/PDF/")
cairo_pdf("p_trace.pdf", width = 90/25.6, height = 62/25.6, family = "Arial")
p_trace_m855
dev.off()

## for the sum
dat_cell_trace_sum <- rbind(dat_cell_trace_pre_sta_m855, dat_cell_trace_test_sta_m855) %>%
  as_tibble() %>% 
  subset(., .$Group!="Neutral") %>%
  ddply(., .(Time, Day), summarise,n=length(mean_value),mean=mean(mean_value),sd=sd(mean_value),se=sd(mean_value)/sqrt(length(mean_value))) %>%
  mutate(Day = factor(Day, levels = c("Pre",  "Cond")))

p_trace_sum <- ggplot(dat_cell_trace_sum, aes(Time, mean, colour=Day, group = Day))+
  geom_line()+
  geom_ribbon(data=subset(dat_cell_trace_sum,Time>=0 ),aes(x=Time,ymax=mean),ymin=0,alpha=0.3) +
  scale_colour_manual(values=c("deepskyblue4", "indianred"))+
  scale_fill_manual(values=c("deepskyblue4", "indianred"))+
  labs(x="Time relative to crossing (s)", y="z score")+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  theme(legend.title = element_blank(), legend.position = c(0.2, 0.8))+
  geom_hline(yintercept = 0, linetype="dashed", color = "gray")

setwd("~cchen/Documents/neuroscience/Pn\ project/Figure/PDF/")
cairo_pdf("p_trace_sum.pdf", width = 75/25.6, height = 60/25.6, family = "Arial")
p_trace_sum
dev.off()


## raster plot the peak ------------
cc_peak_fun <- function(path_peak, frame_crossing){
  ID <- str_extract(path_peak, regex("m\\d+"))
  frame_crossing_range <- c((frame_crossing-40):(frame_crossing+140))
  dat_raster <- read.csv(path_peak) %>% 
    as_tibble() %>% 
    gather() %>% 
    mutate(value = ifelse(value %in% frame_crossing_range, value, NA)) %>% 
    drop_na() %>% 
    mutate(value = value - frame_crossing) %>% 
    mutate(key = str_c(key, ID))
  return(dat_raster)
}

## for m3
path_peak_m3 <- as.list(list.files(path ="~cchen/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/matlab_analysis/m3/peaks_days", pattern = "*.csv", full.names = T ))[c(3, 6)]
dat_rate_m3 <- mapply(cc_peak_fun,path_peak_m3, t_stim_m3, SIMPLIFY = F) 

## for m7
path_peak_m7 <- as.list(list.files(path ="~cchen/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/matlab_analysis/m7/peaks_days/", pattern = "*.csv", full.names = T ))[c(3, 6)]
dat_rate_m7 <- mapply(cc_peak_fun,path_peak_m7, t_stim_m7, SIMPLIFY = F) 

## for m17
path_peak_m17 <- as.list(list.files(path ="~cchen/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/matlab_analysis/m17/peaks_days/", pattern = "*.csv", full.names = T ))[c(3, 6)]
dat_rate_m17 <- mapply(cc_peak_fun,path_peak_m17, t_stim_m17, SIMPLIFY = F) 

## for m18
path_peak_m18 <- as.list(list.files(path ="~cchen/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/matlab_analysis/m18/peaks_days/", pattern = "*.csv", full.names = T ))[c(3, 6)]
dat_rate_m18 <- mapply(cc_peak_fun,path_peak_m18, t_stim_m18, SIMPLIFY = F) 

## for m855
path_peak_m855 <- as.list(list.files(path ="~cchen/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/matlab_analysis/m855/peaks_days/", pattern = "*.csv", full.names = T ))[c(3, 6)]
dat_rate_m855 <- mapply(cc_peak_fun,path_peak_m855, t_stim_m855, SIMPLIFY = F) 

## for m855
path_peak_m857 <- as.list(list.files(path ="~cchen/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/matlab_analysis/m857/peaks_days/", pattern = "*.csv", full.names = T ))[c(1, 4)]
dat_rate_m857 <- mapply(cc_peak_fun,path_peak_m857, t_stim_m857, SIMPLIFY = F) 

dat_raster_pre <- rbind(dat_rate_m3[[1]],dat_rate_m7[[1]],dat_rate_m17[[1]],dat_rate_m18[[1]],dat_rate_m855[[1]],dat_rate_m857[[1]] ) %>%
  mutate(Group ="Pre")

dat_raster_post <- rbind(dat_rate_m3[[2]],dat_rate_m7[[2]],dat_rate_m17[[2]],dat_rate_m18[[2]],dat_rate_m855[[2]],dat_rate_m857[[2]] ) %>% 
  mutate(Group ="Cond")

p_raster <- rbind(dat_raster_pre, dat_raster_post) %>%  
  mutate(Group = factor(Group, levels = c("Pre", "Cond"))) %>% 
  mutate(value = value/20) %>% 
  ggplot(., aes(value, key, fill= Group, color= Group))+
  geom_tile(size = 1)+
  facet_grid(rows = vars(Group))+
  labs(x="Time relative to crossing (s)", y = "Neuron ID")+
  scale_fill_manual(values=c("deepskyblue4", "indianred"))+
  scale_colour_manual(values=c("deepskyblue4", "indianred"))+
  theme(axis.line.x = element_line(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.line.y  = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  theme(legend.title = element_blank(), legend.position = "none")+
  geom_vline(xintercept = 0, color="red")

setwd("~cchen/Documents/neuroscience/Pn\ project/Figure/PDF/")
cairo_pdf("p_raster.pdf", width = 90/25.6, height = 60/25.6, family = "Arial")
p_raster
dev.off()

## calculate the spikes before and after crossing------
cc_firing_fun <- function(path_trace, path_peak, t_stim) {
  cell_valid <- read.xlsx(path_trace, colNames = F, rowNames = F) %>%
    select (., X1)

  ctrl_rang <- c((t_stim-42): (t_stim-2))
  test_rang <- c(t_stim: (t_stim +140))
  
  
  dat_peak_before <- read.csv(path_peak) %>% 
    select(which(cell_valid==1)) %>%
    apply(., 2, function(x) intersect(x, ctrl_rang), simplify = F) %>% 
    sapply(., length)/2 %>% 
    round(., digits = 3)
  
  
  
  dat_peak_after <- read.csv(path_peak) %>% 
    select(which(cell_valid==1)) %>% 
    apply(., 2, function(x) intersect(x, test_rang), simplify = F) %>% 
    sapply(., length)/7 %>% 
    round(., digits = 3)
  
  diff_rate <- unname(dat_peak_after - dat_peak_before) %>% 
    .[.!=0]
  
  diff_rate1 <- ifelse(length(diff_rate)>0, mean(diff_rate), 0)
  
  
  return(diff_rate1)
}

## for m3
dat_rate_m3 <- mapply(cc_firing_fun, path_trace_m3, path_peak_m3, t_stim_m3)

## for m7
dat_rate_m7 <- mapply(cc_firing_fun, path_trace_m7, path_peak_m7, t_stim_m7)

## for m17
dat_rate_m17 <- mapply(cc_firing_fun, path_trace_m17, path_peak_m17, t_stim_m17)

## for m18
dat_rate_m18 <- mapply(cc_firing_fun, path_trace_m18, path_peak_m18, t_stim_m18)

## for m855
dat_rate_m855 <- mapply(cc_firing_fun, path_trace_m855, path_peak_m855, t_stim_m855)

## for m857

dat_rate_m857 <- mapply(cc_firing_fun, path_trace_m857, path_peak_m857, t_stim_m857)

dat_firing <- tibble(value = c(dat_rate_m3, dat_rate_m7, dat_rate_m17, dat_rate_m18, dat_rate_m855, dat_rate_m857)) %>% 
  mutate(ID = rep(mouse_ID, each = 2)) %>% 
  mutate(Day = rep(c("Pre", "Cond"), length(mouse_ID))) %>% 
  mutate(Day = factor(Day, levels = c("Pre", "Cond")))

p_firing <- dat_firing %>% 
  ggplot(., aes(Day, value, colour=Day))+
  geom_boxplot(outlier.shape = NA)+
  geom_line(aes(group=ID), colour="gray90")+
  geom_jitter(width = 0.2, shape=1, size=1)+
  scale_colour_manual(values=c("deepskyblue4", "indianred"))+
  labs(x="", y="Diff. of firing frequency (Hz)")+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  theme(legend.position = 'none')

t_firng <- dat_firing %>% 
  wilcox.test(value~Day,.)

setwd("~cchen/Documents/neuroscience/Pn\ project/Figure/PDF/")
cairo_pdf("p_firing.pdf", width = 40/25.6, height = 60/25.6, family = "Arial")
p_firing
dev.off() 


## plot the ca2+ transient events----
cc_duration_fun <- function(path_trace, path_peak, t_stim) {
  
  mouse_ID <- str_extract(path_trace, regex("m\\d+"))
  dat_trace_cell <- read.xlsx(path_trace, colNames = F, rowNames = F) 
  
  duration_rang <- c((t_stim-40): (t_stim+140))
  
  
  dat_peak_before <- read.csv(path_peak) %>% 
    select(which(dat_trace_cell$X1==1)) %>%
    apply(., 2, function(x) intersect(x, duration_rang), simplify = F) 
  
  dat_peak<- c()
  
  dat_peak_trace <- c()
  for (i in c(1: length(dat_peak_before))) {
    if (isempty(dat_peak_before[[i]])) next
    peak_name <- as.numeric(str_extract(names(dat_peak_before)[i], "(\\d)+"))
    peak_time <- dat_peak_before[[i]] 
    cell_trace <- dat_trace_cell %>%
      as_tibble() %>% 
      select(-1) %>% 
      slice(peak_name) %>% 
      select(all_of(peak_time)) %>% 
      t()
    
    dat_peak<- rbind(dat_peak, cell_trace)
    
    #  cell_trace1 <- dat_trace_cell %>%
    #    as_tibble() %>%
    #    select(-1) %>%
    #    slice(peak_name) %>%
    #    select(all_of(duration_rang)) %>%
    #    unlist(., use.names=FALSE)
    # 
    # 
    # plot(cell_trace1, type="l", main = mouse_ID)
    # points((peak_time - duration_rang[1]), cell_trace, col="red")

    ## get the trace for average
    
    for (j in seq_along(peak_time)){
      if (cell_trace[j]< 0.015|peak_time[j]< 100) next
      peak_duration <- c((peak_time[j]-100): (peak_time[j] + 300))
      dat_cell_peak_trace <- dat_trace_cell %>%
        as_tibble() %>% 
        select(-1) %>% 
        slice(peak_name) %>% 
        select(all_of(peak_duration)) %>% 
        unlist(., use.names=FALSE)
      dat_peak_trace <- cbind(dat_peak_trace, dat_cell_peak_trace)
      
      plot(dat_cell_peak_trace, type="l", main = mouse_ID)
      points(100, dat_cell_peak_trace[100], col = "red")
      
      
    }
    
  }
  return(list(dat_peak, dat_peak_trace))
}


path_peak_m3 <- as.list(list.files(path ="~cchen/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/matlab_analysis/m3/m3_07212022/peak_path/", pattern = "*.csv", full.names = T ))
path_trace_m3 <- as.list(list.files(path ="~cchen/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/matlab_analysis/m3/m3_07212022/trace_path", pattern = "*.xlsx", full.names = T ))

dat_amp_cond_m3 <- cc_duration_fun(path_trace_m3[[2]],path_peak_m3[[2]], t_stim_m3[[2]] )
dat_amp_pre_m3 <- cc_duration_fun(path_trace_m3[[1]],path_peak_m3[[1]], t_stim_m3[[1]] )


dat_amp_cond_m7 <- cc_duration_fun(path_trace_m7[[2]],path_peak_m7[[2]], t_stim_m7[[2]] )
dat_amp_pre_m7 <- cc_duration_fun(path_trace_m7[[1]],path_peak_m7[[1]], t_stim_m7[[1]] )

dat_amp_cond_m18 <- cc_duration_fun(path_trace_m18[[2]],path_peak_m18[[2]], t_stim_m18[[2]] )
dat_amp_pre_m18 <- cc_duration_fun(path_trace_m18[[1]],path_peak_m18[[1]], t_stim_m18[[1]] )

dat_amp_cond_m17 <- cc_duration_fun(path_trace_m17[[2]],path_peak_m17[[2]], t_stim_m17[[2]] )
dat_amp_pre_m17 <- cc_duration_fun(path_trace_m17[[1]],path_peak_m17[[1]], t_stim_m17[[1]] )

dat_amp_cond_m855 <- cc_duration_fun(path_trace_m855[[2]],path_peak_m855[[2]], t_stim_m855[[2]] )
dat_amp_pre_m855 <- cc_duration_fun(path_trace_m855[[1]],path_peak_m855[[1]], t_stim_m855[[1]] )

dat_amp_cond_m857 <- cc_duration_fun(path_trace_m857[[2]],path_peak_m857[[2]], t_stim_m857[[2]] )
dat_amp_pre_m857 <- cc_duration_fun(path_trace_m857[[1]],path_peak_m857[[1]], t_stim_m857[[1]] )

pre_amp <- c(dat_amp_pre_m3[[1]],dat_amp_pre_m7[[1]],dat_amp_pre_m17[[1]],dat_amp_pre_m18[[1]], dat_amp_pre_m855[[1]], dat_amp_pre_m857[[1]] )
cond_amp <- c(dat_amp_cond_m3[[1]], dat_amp_cond_m7[[1]],dat_amp_cond_m17[[1]],dat_amp_cond_m18[[1]], dat_amp_cond_m855[[1]], dat_amp_cond_m857[[1]] )

p_amp <- tibble(value = c(pre_amp, cond_amp)) %>% 
  mutate(Group = c(rep("Pre", length(pre_amp)), rep("Cond", length(cond_amp)))) %>% 
  filter(value > 0.015) %>% 
  ggplot(., aes(Group, value, color = Group))+
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(width = 0.2, shape=1)

## plot the trace
dat_trace_pre <- dat_amp_pre_m7[[2]] %>% 
  as_tibble() %>% 
  mutate(Time = 1:401) %>% 
  pivot_longer(-Time) 

p_trace_pre <- dat_trace_pre %>% 
  ggplot(., aes(Time, value, color= name))+
  geom_line()+
  facet_grid(rows = vars(name))+
  theme_void()


