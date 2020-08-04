## Project note------
## The script are used to analyze the data from miniscope experiemnt analyzed by matlab
## updated on 01-28-2020
## data are analyzed with inscopix data processing images
## output: 1. The trace crossing the border; 2. cells proportion increased or decrease; 3. Heat maps
## 4. Number of cells active, inhibited
## updated: 03182020, add more data
## updated: 04012020, change the way to do z-score
## updated: 04142020, do k-mean analysis firstly and back cat the cells activation
## updated: 05062020, only analyze the first across 
## updated:05112020, add code for ctrl experiment (mouse in homecage)
## updated: 05222020, calculate the spiking frequency before and after crossing

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

## function for analysis-----
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
  #stim_time<- seq(-5, 8.5, by=0.5)
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
path_trace_m3 <- as.list(list.files(path ="~cchen2/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/matlab_analysis/m3/trace_days/", pattern = "*.xlsx", full.names = T ))
t_stim_m3_d1 <- c(214) 
t_stim_m3_d2 <- c(274)
t_stim_m3_d3 <- c(1885) 
t_stim_m3_d4 <- c(1813)
t_stim_m3_d5 <- c(758) 
t_stim_m3_d6 <- c(132) 
t_stim_m3_d7 <- c(346) ## use first two 346, 1206, 1756, 2288, 2970, 3262, 3942

t_stim_m3 <- list(t_stim_m3_d1, t_stim_m3_d2, t_stim_m3_d3, t_stim_m3_d4, t_stim_m3_d5, t_stim_m3_d6, t_stim_m3_d7)

dat_trace_m3 <- mapply(c_miniscope_matlab, path_trace_m3, t_stim_m3, SIMPLIFY = F)

# for m7
path_trace_m7 <- as.list(list.files(path ="~cchen2/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/matlab_analysis/m7/trace_days/", pattern = "*.xlsx", full.names = T ))
t_stim_m7_d1 <- c(605) 
t_stim_m7_d2 <- c(798)
t_stim_m7_d3 <- c(360) 
t_stim_m7_d4 <- c(3236)
t_stim_m7_d5 <- c(317) 
t_stim_m7_d6 <- c(205)
t_stim_m7_d7 <- c(236) ## use first two 236, 1914,2777,3334,4091

t_stim_m7 <- list(t_stim_m7_d1, t_stim_m7_d2, t_stim_m7_d3, t_stim_m7_d4, t_stim_m7_d5, t_stim_m7_d6, t_stim_m7_d7)

dat_trace_m7 <- mapply(c_miniscope_matlab, path_trace_m7, t_stim_m7,SIMPLIFY = F)

# for m17
path_trace_m17 <- as.list(list.files(path ="~cchen2/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/matlab_analysis/m17/trace_days/", pattern = "*.xlsx", full.names = T ))
t_stim_m17_d1 <- c(815) 
t_stim_m17_d2 <- c(1103)
t_stim_m17_d3 <- c(437) 
t_stim_m17_d4 <- c(682)
t_stim_m17_d5 <- c(1114) 
t_stim_m17_d6 <- c(617) 
t_stim_m17_d7 <- c(157) # 157, 2536, 2961

t_stim_m17 <- list(t_stim_m17_d1, t_stim_m17_d2, t_stim_m17_d3, t_stim_m17_d4, t_stim_m17_d5, t_stim_m17_d6, t_stim_m17_d7)

dat_trace_m17 <- mapply(c_miniscope_matlab, path_trace_m17, t_stim_m17, SIMPLIFY = F)

# for m18
path_trace_m18 <- as.list(list.files(path ="~cchen2/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/matlab_analysis/m18/trace_days/", pattern = "*.xlsx", full.names = T ))
t_stim_m18_d1 <- c(490) 
t_stim_m18_d2 <- c(236)
t_stim_m18_d3 <- c(784) 
t_stim_m18_d4 <- c(493)
t_stim_m18_d5 <- c(195) 
t_stim_m18_d6 <- c(467) 
t_stim_m18_d7 <- c(124) ## use the first two of 124,2238,3000

t_stim_m18 <- list(t_stim_m18_d1, t_stim_m18_d2, t_stim_m18_d3, t_stim_m18_d4, t_stim_m18_d5, t_stim_m18_d6, t_stim_m18_d7)

dat_trace_m18 <- mapply(c_miniscope_matlab, path_trace_m18, t_stim_m18, SIMPLIFY = F)

# for m855
path_trace_m855 <- as.list(list.files(path ="~cchen2/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/matlab_analysis/m855/trace_days", pattern = "*.xlsx", full.names = T ))
t_stim_m855_d1 <- c(223)*2
t_stim_m855_d2 <- c(513)*2
t_stim_m855_d3 <- c(392) *2
t_stim_m855_d4 <- c(550)*2
t_stim_m855_d5 <- c(67)*2 
t_stim_m855_d6 <- c(64) *2
t_stim_m855_d7 <- c(41)*2

t_stim_m855 <- list(t_stim_m855_d1, t_stim_m855_d2, t_stim_m855_d3, t_stim_m855_d4, t_stim_m855_d5, t_stim_m855_d6, t_stim_m855_d7)

dat_trace_m855 <- mapply(c_miniscope_matlab, path_trace_m855, t_stim_m855, SIMPLIFY = F)

## combine data and do k-means analysis-----

dat_cell_trace <- vector(mode = "list", 7)
for (i in 1:7) {
  dat_cell_trace[[i]] <- cbind(dat_trace_m3[[i]], dat_trace_m7[[i]], dat_trace_m17[[i]], dat_trace_m18[[i]], dat_trace_m855[[i]])

}


dat_cell_trace_re <- vector(mode = "list", 7)

# k-menas clustering
mouse_ID <- c("m3", "m7", "m17", "m18", "m855")
for (i in 1:length(dat_cell_trace)) {
  dat_cell_trace_d <- dat_cell_trace[[i]]
  cell_id <- sprintf("Cell%s",seq(1:ncol(dat_cell_trace_d)))
  colnames(dat_cell_trace_d)<- cell_id
  
  dat_cell_trace_d_cluster <- unname(kmeans(t(dat_cell_trace_d), centers = 3, nstart = 20)[[1]])
  stim_time<- seq(-2, 6.5, by=0.5)
  #stim_time<- seq(-5, 8.5, by=0.5)
  
  dat_cell_trace_d_re<- as.data.frame(dat_cell_trace_d) %>% 
    mutate(., Time= stim_time) %>% 
    melt(., id.vars='Time') %>% 
    mutate(., Group= rep(dat_cell_trace_d_cluster, each=length(stim_time)) )
    
  ## sort data by the value in each group
  dat_cell_d_sort <- as.numeric(names(sort(tapply(dat_cell_trace_d_re$value, dat_cell_trace_d_re$Group, mean), decreasing = T)))
  
  dat_cell_trace_d_re$Group[dat_cell_trace_d_re$Group == dat_cell_d_sort[1]] ="Excited"
  dat_cell_trace_d_re$Group[dat_cell_trace_d_re$Group == dat_cell_d_sort[2]] ="Neutral"
  dat_cell_trace_d_re$Group[dat_cell_trace_d_re$Group == dat_cell_d_sort[3]] ="Inhibited"
  
  rep_time <- c(ncol(dat_trace_m3[[i]]), ncol(dat_trace_m7[[i]]), ncol(dat_trace_m17[[i]]), ncol(dat_trace_m18[[i]]), ncol(dat_trace_m855[[i]]) )
  dat_cell_trace_d_re <- mutate(dat_cell_trace_d_re, Group= factor(Group, levels = c("Excited", "Neutral", "Inhibited"))) %>% 
    mutate(., ID = rep(mouse_ID, rep_time*length(stim_time)) )

  dat_cell_trace_re[[i]] <- dat_cell_trace_d_re
}


## heat plot of d3, d6 and d7----
## set the range of z score during all days
score_range <- range(mapply(function(x) x$value, dat_cell_trace_re, SIMPLIFY = T))

for (i in c(3, 6, 7)) {
  dat_trace <- dat_cell_trace_re[[i]]
  dat_trace_sta <- ddply(dat_trace, .(variable, Group), summarise,mean=mean(value), sum=sum(value))
  dat_trace_sta <- dat_trace_sta[order(dat_trace_sta[,'mean']),]
  dat_trace$variable <- factor(dat_trace$variable, levels = dat_trace_sta$variable)
  p_heat <- ggplot(dat_trace, aes(Time, variable,fill= value))+ 
    geom_tile(height=2)+
    #facet_grid(rows = vars(Group), scales = "free_y")+
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

setwd("~cchen2/Documents/neuroscience/Pn\ project/Figure/PDF/")
cairo_pdf("p_anti_heat.pdf", width = 170/25.6, height = 65/25.6, family = "Arial")
p_heat_com
dev.off()

## plot trace by group-----
dat_cell_trace_pre_sta <- ddply(dat_cell_trace_re[[2]], .(ID,Time, Group), summarise,value=mean(value)) %>%
  rbind(.,ddply(dat_cell_trace_re[[3]], .(ID,Time, Group), summarise,value=mean(value)) ) %>%
  ddply(., .(ID, Time, Group),summarise,value=mean(value) )

dat_cell_trace_con_sta <- ddply(dat_cell_trace_re[[5]], .(ID,Time, Group), summarise,value=mean(value)) %>%
  rbind(.,ddply(dat_cell_trace_re[[6]], .(ID,Time, Group), summarise,value=mean(value)) ) %>%
  ddply(., .(ID, Time, Group),summarise,value=mean(value) )

dat_cell_trace_test_sta <- ddply(dat_cell_trace_re[[7]], .(ID,Time, Group), summarise,value=mean(value))

dat_cell_trace_sta <- ddply(dat_cell_trace_pre_sta, .(Time, Group), summarise,n=length(value),mean=mean(value),sd=sd(value),se=sd(value)/sqrt(length(value))) %>%
  rbind(., ddply(dat_cell_trace_con_sta, .(Time, Group), summarise,n=length(value),mean=mean(value),sd=sd(value),se=sd(value)/sqrt(length(value)))) %>%
  rbind(., ddply(dat_cell_trace_test_sta, .(Time, Group), summarise,n=length(value),mean=mean(value),sd=sd(value),se=sd(value)/sqrt(length(value)))) %>%
  mutate(., Day = rep(c("Pre", "Cond.", "Test"),each= length(stim_time)*3))

dat_cell_trace_sta$Day <- factor(dat_cell_trace_sta$Day, levels = c("Pre", "Cond.", "Test"))
dat_cell_trace_sta$Group <- factor(dat_cell_trace_sta$Group, levels = c("Neutral", "Excited", "Inhibited"))
range_trace_plot <- range(dat_cell_trace_sta$mean)

## group by day
p_trace <- ggplot(dat_cell_trace_sta, aes(Time, mean, colour=Group))+
  geom_line()+
  geom_ribbon(aes(ymin=mean-se, ymax=mean+se, fill=Group), alpha=0.1, linetype=0)+
  facet_grid(cols = vars(Day))+
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
  theme(legend.title = element_blank(), legend.position = 'none')

setwd("~cchen2/Documents/neuroscience/Pn\ project/Figure/PDF/")
cairo_pdf("p_trace.pdf", width = 165/25.6, height = 65/25.6, family = "Arial")
p_trace
dev.off()

## plot E and I points to show correlation-----
## points to show the correlation
dat_cell_cor <- NULL
cell_area_day <- c("Rm","Pre", "Pre", "Rm","Cond.","Cond.", "Test")

for (i in c(3, 6, 7)) {
  dat_trace <- subset(dat_cell_trace_re[[i]], dat_cell_trace_re[[i]]$ID=="m3") %>%
    subset(., .$Group != "Neutral") %>%
    ddply(., .(Time, Group), summarise, value=mean(value)) %>%
    dcast(., Time~Group) %>%
    mutate(., Day = cell_area_day[i])
  dat_cell_cor <- rbind(dat_cell_cor, dat_trace)
  
}
dat_cell_cor$Day <- factor(dat_cell_cor$Day, levels = c("Pre", "Cond.", "Test"))

p_trace_cor <- ggplot(dat_cell_cor, aes(Excited, Inhibited, colour=Day))+
  geom_point()+
  geom_smooth(method = "lm", lwd=0.8)+
  scale_colour_manual(values=c("seagreen", "indianred", "deepskyblue4"))+  
  labs(x="Excited (z-score)", y="Inhibited (z-score)")+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  theme(legend.title = element_blank(), legend.position = c(0.8, 0.8))

setwd("~cchen2/Documents/neuroscience/Pn\ project/Figure/PDF/")
cairo_pdf("p_trace_cor.pdf", width = 60/25.6, height = 65/25.6, family = "Arial")
p_trace_cor
dev.off()
## correlation for each mice
dat_cell_trace_cor <- NULL
cell_area_day <- c("Rm","Pre", "Pre", "Rm","Cond.","Cond.", "Test")
for (i in c(2,3,5,6,7)){
  dat_trace_cor <- dat_cell_trace_re[[i]] %>%
    ddply(., .(ID, Time, Group), summarise, value=mean(value,na.rm = T)) %>%
    subset(., .$Group!="Neutral") %>%
    dcast(., Time + ID ~ Group) %>%
    ddply(., .(ID), summarise, "corr" = cor(Excited, Inhibited, method = "spearman")) %>%
    mutate(., Day = cell_area_day[i])
  dat_cell_trace_cor <- rbind(dat_cell_trace_cor, dat_trace_cor)
  }

dat_cor_value <- ddply(dat_cell_trace_cor, .(ID, Day), summarise, "corr" = mean(corr)) %>%
  mutate(Day = factor(Day, levels = c("Pre", "Cond.", "Test")))
  
p_cor <-ggplot(dat_cor_value, aes(Day, corr, fill=Day))+
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(width = 0.2, shape=1)+
  scale_fill_manual(values=c("seagreen", "indianred", "deepskyblue4"))+
  labs(x="", y="Correlation coefficient (r)")+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  theme(legend.title = element_blank(), legend.position = "none")

setwd("~cchen2/Documents/neuroscience/Pn\ project/Figure/PDF/")
cairo_pdf("p_cor.pdf", width = 40/25.6, height = 65/25.6, family = "Arial")
p_cor
dev.off() 

  
## calculate the sum of active and inhibited trace----
dat_cell_trace_sum <- rbind(dat_cell_trace_pre_sta, dat_cell_trace_con_sta, dat_cell_trace_test_sta) %>%
  as_tibble() %>% 
  mutate(Day = rep(c("Pre", "Cond.", "Test"), c(nrow(dat_cell_trace_pre_sta),nrow(dat_cell_trace_con_sta),nrow(dat_cell_trace_test_sta)))) %>%
  subset(., .$Group!="Neutral") %>%
  ddply(., .(ID,Time, Day), summarise, value=sum(value)) %>%
  ddply(., .(Time, Day), summarise,n=length(value),mean=mean(value),sd=sd(value),se=sd(value)/sqrt(length(value))) %>%
  mutate(Day = factor(Day, levels = c("Pre", "Cond.", "Test")))

p_trace_sum <- ggplot(dat_cell_trace_sum, aes(Time, mean, colour=Day))+
  geom_line()+
  geom_ribbon(aes(ymin=mean-se, ymax=mean+se, fill=Day), alpha=0.1, linetype=0)+
  scale_colour_manual(values=c("seagreen", "indianred", "deepskyblue4"))+
  scale_fill_manual(values=c("seagreen", "indianred", "deepskyblue4"))+
  labs(x="Time relative to crossing (s)", y="AUC (z score)")+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  geom_vline(xintercept = 0, col= 'red', linetype=2)+
  theme(legend.title = element_blank(), legend.position = 'top')

setwd("~cchen2/Documents/neuroscience/Pn\ project/Figure/PDF/")
cairo_pdf("p_trace_sum.pdf", width = 60/25.6, height = 70/25.6, family = "Arial")
p_trace_sum
dev.off()

## for EI change analysis-----
dat_cell_area <- c()
cell_area_day <- c("Rm","Pre", "Pre", "Rm","Cond.","Cond.", "Test")
for (i in c(2,3,5,6,7)){
  rep_time <- c(ncol(dat_trace_m3[[i]]), ncol(dat_trace_m7[[i]]), ncol(dat_trace_m17[[i]]), ncol(dat_trace_m18[[i]]), ncol(dat_trace_m855[[i]]))
  dat_trace <- dat_cell_trace_re[[i]]
  dat_trace$ID <- rep(mouse_ID, rep_time*length(stim_time))
  dat_trace_sta <- ddply(dat_trace, .(ID, Time, Group), summarise,n=length(value),mean=mean(value),sd=sd(value),se=sd(value)/sqrt(length(value)))
  Time_area <- which(stim_time>=0)
  dat_anti_area <- as.data.frame(tapply(dat_trace_sta$mean, INDEX = list(dat_trace_sta$ID, dat_trace_sta$Group), 
                                        function (x) AUC(stim_time[Time_area], x[Time_area])))
  dat_anti_area$ID <- rownames(dat_anti_area)
  rownames(dat_anti_area)<-NULL
  dat_anti_area$Day <- cell_area_day[i]
  dat_anti_area$ratio <- abs(dat_anti_area$Excited/dat_anti_area$Inhibited)
  dat_anti_area$sum <- dat_anti_area$Excited + dat_anti_area$Inhibited
  dat_cell_area <- rbind(dat_cell_area, dat_anti_area)
}

dat_cell_area$Day <- factor(dat_cell_area$Day, levels = c("Pre", "Cond.","Test"))
dat_cell_area_sta <-  ddply(dat_cell_area, .(ID, Day), summarise, value=mean(ratio, na.rm = T))

p_EI_ratio <- ggplot(dat_cell_area_sta, aes(Day, value, fill=Day))+
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(width = 0.2, shape=1)+
  scale_fill_manual(values=c("seagreen", "indianred", "deepskyblue4"))+
  labs(x="", y="E/I ratio")+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  scale_y_continuous(limits = c(0, 4), expand = c(0,0))+
  theme(legend.title = element_blank(), legend.position = "none")+
  annotate(x=c(1,1,2,2), y=c(3.4,3.5,3.5,3.4),"path")+
  annotate("text",x=1.5,y=3.5, label="***", size=5)


## statistic test
t_EI<-aov(value~Day, data=dat_cell_area_sta)
summary(t_EI)
pairwise.t.test(dat_cell_area_sta$value, dat_cell_area_sta$Day, paired = T)


setwd("~cchen2/Documents/neuroscience/Pn\ project/Figure/PDF/")
cairo_pdf("p_EI_ratio.pdf", width = 40/25.6, height = 65/25.6, family = "Arial")
p_EI_ratio
dev.off() 
## plot the sum of E-I----
dat_cell_sum_sta <-  ddply(dat_cell_area, .(ID, Day), summarise, "sum"= mean(sum, na.rm = T))


p_EI_sum<- ggplot(dat_cell_sum_sta, aes(Day, sum, fill=Day))+
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(width = 0.2, shape=1)+
  scale_fill_manual(values=c("seagreen", "indianred", "deepskyblue4"))+
  labs(x="", y="AUG (z-score)")+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  #scale_y_continuous(limits = c(0, 10), expand = c(0,0))+
  theme(legend.title = element_blank(), legend.position = "none")+
  annotate(x=c(1,1,2,2), y=c(9,9.1,9.1,9),"path")+
  annotate("text",x=1.5,y=9.1, label="*", size=5)
  

## statistic test
t_sum_EI<-aov(sum~Day, data=dat_cell_sum_sta)
summary(t_EI)
pairwise.t.test(dat_cell_sum_sta$sum, dat_cell_sum_sta$Day, paired = T)

## show the change of E, I and sum
p_EIS_change <- dat_cell_area %>%
  select("Excited", "Inhibited","ID","Day", "sum") %>%
  melt(., id.vars= c('Day', "ID")) %>%
  ddply(., .(ID, Day, variable), summarise, value=mean(value, na.rm = T)) %>%
  ggplot(., aes(variable, value, fill=Day))+
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(position=position_jitterdodge(), shape=1, size=1)+
  scale_fill_manual(values=c("seagreen", "indianred", "deepskyblue4"))+
  labs(x="", y="AUC (z-score)")+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  theme(legend.title = element_blank(), legend.position = "top")

setwd("~cchen2/Documents/neuroscience/Pn\ project/Figure/PDF/")
cairo_pdf("p_EIS_change.pdf", width = 70/25.6, height = 72/25.6, family = "Arial")
p_EIS_change
dev.off()
  
## portion of cell catlog------

dat_cell_cat <- mapply(function (x) subset(x, x$Time==0), dat_cell_trace_re, SIMPLIFY = F ) 

dat_cell_cat1 <- mapply(function(x) prop.table(table(x$ID, x$Group), 1), dat_cell_cat, SIMPLIFY = F) %>%
  do.call(rbind, .) %>%
  as.data.frame() %>%
  mutate(ID = rep(rownames(.)[1:4], 7), Day=rep(c("Rm", "Pre", "Pre", "Rm", "Cond.", "Cond.", "Test"), each=4)) %>%
  subset(., .$Day!="Rm") %>%
  melt(., id.vars = c("ID", "Day")) %>%
  ddply(., .(ID, Day, variable), summarise, "Prop"=mean(value)) %>%
  mutate(Day = factor(Day, levels = c("Pre", "Cond.", "Test"))) %>%
  mutate(variable=factor(variable, levels = c("Excited", "Neutral", "Inhibited")))

p_cat<- ggplot(dat_cell_cat1, aes(variable, Prop, fill=Day))+
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(position=position_jitterdodge(), shape=1, size=1)+
  scale_fill_manual(values=c("seagreen", "indianred", "deepskyblue4"))+
  labs(x="", y="% of all neurons")+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  scale_y_continuous(limits = c(0,1),expand = c(0,0))+
  theme(legend.title = element_blank(), legend.position = "top")

setwd("~cchen2/Documents/neuroscience/Pn\ project/Figure/PDF/")
cairo_pdf("p_cat.pdf", width = 70/25.6, height = 72/25.6, family = "Arial")
p_cat
dev.off()

## calculate the spikes before and after crossing------
cc_firing_fun <- function(path_trace, path_peak, t_stim) {
  cell_valid <- read.xlsx(path_trace, colNames = F, rowNames = F) %>% 
    select (., X1)
  
  dat_peak <- read.csv(path_peak) %>% 
    select(which(cell_valid==1)) %>% 
    unlist() %>% 
    na.omit()
  
  ctrl_rang <- range(t_stim-40, t_stim)
  test_rang <- range(t_stim, t_stim +80)
  
  ctrl_freq <- length(dat_peak[dat_peak>= ctrl_rang[1] & dat_peak < ctrl_rang[2]])/2 # 2s before crossing
  test_freq <- length(dat_peak[dat_peak>= test_rang[1] & dat_peak < test_rang[2]])/4 # 4s after crossing
  return(c(ctrl_freq, test_freq))
}

## for m3

path_peak_m3 <- as.list(list.files(path ="~cchen2/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/matlab_analysis/m3/peaks_days/", pattern = "*.csv", full.names = T ))

dat_rate_m3 <- t(mapply(cc_firing_fun, path_trace_m3, path_peak_m3, t_stim_m3))

## for m7
path_peak_m7 <- as.list(list.files(path ="~cchen2/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/matlab_analysis/m7/peaks_days/", pattern = "*.csv", full.names = T ))

dat_rate_m7 <- t(mapply(cc_firing_fun, path_trace_m7, path_peak_m7, t_stim_m7))

## for m17
path_peak_m17 <- as.list(list.files(path ="~cchen2/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/matlab_analysis/m17/peaks_days/", pattern = "*.csv", full.names = T ))
dat_rate_m17 <- t(mapply(cc_firing_fun, path_trace_m17, path_peak_m17, t_stim_m17))

## for m18
path_peak_m18 <- as.list(list.files(path ="~cchen2/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/matlab_analysis/m18/peaks_days/", pattern = "*.csv", full.names = T ))

dat_rate_m18 <- t(mapply(cc_firing_fun, path_trace_m18, path_peak_m18, t_stim_m18))

## combine data and plot
dat_firing <- rbind(dat_rate_m3, dat_rate_m7, dat_rate_m17, dat_rate_m18) %>% 
  as_tibble() %>% 
  rename(Freq_before = V1, Freq_after = V2) %>% 
  mutate(Day = rep(str_c("D", 1:7), 4)) %>% 
  mutate(ID = rep(mouse_ID, each=7)) %>% 
  mutate(Group = rep(cell_area_day, 4)) %>% 
  filter(Group != "Rm") %>% 
  gather(variable, value, -ID, -Day, -Group) %>% 
  ddply(., .(ID, Group, variable, Day), summarise, value=mean(value)) %>% 
  ddply(., .(ID, Group, variable), summarise, value=mean(value)) %>% 
  mutate(Group= factor(Group, levels = c("Pre", "Cond.", "Test")), 
         variable=factor(variable, levels = c("Freq_before", "Freq_after")))

p_firing <- ggplot(dat_firing, aes(Group, value, fill=variable))+
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(position=position_jitterdodge(), shape=1, size=1)+
  scale_fill_manual(values=c("seagreen", "indianred"))+
  labs(x="", y="Firing rate (Hz)")+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  # scale_y_continuous(limits = c(0,3),expand = c(0,0))+
  theme(legend.title = element_blank(), legend.position = "top")

## within and cross cells distance-----
## function to trim cell position file with the manually confirmed cells
cc_trim_position <- function(path_trace, path_position){
  cell_valid <- read.xlsx(path_trace, colNames = F, rowNames = F) %>%
    .[,1]
  dat_position <- read.csv(path_position) %>%
    .[cell_valid==1,]
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
  dat_position_com[[i]] <- rbind(dat_position_m3[[i]], dat_position_m7[[i]], dat_position_m17[[i]], dat_position_m18[[i]]) %>%
    mutate(ID = dat_cell_cat[[i]]$ID , Group = dat_cell_cat[[i]]$Group )
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
rep_time <- mapply(nrow, dat_cell_distance_list)

dat_cell_distance <- do.call(rbind, dat_cell_distance_list) %>%
  mutate(Day = rep(c("Rm", "Pre", "Pre", "Rm", "Cond.", "Cond.","Test"), rep_time)) %>%
  subset(., .$Day!="Rm") %>%
  mutate(Group = factor(Group, levels =c("Excited", "Neutral", "Inhibited")), Day = factor(Day, levels = c("Pre", "Cond.", "Test")) )
  

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


dat_cell_distance_sta <- melt(dat_cell_distance[,-c(1,2)], id.vars = c('Day', "ID", "Group")) %>%
  ddply(., .(ID,Day,Group,variable), summarise,value=mean(value, na.rm = T))
  

ggplot(dat_cell_distance_sta, aes(Group, value, fill=variable))+
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(position=position_jitterdodge(), shape=1, size=1)+
  facet_grid(cols = vars(Day))+
  scale_fill_manual(values=c( "black","gray"))+
  labs(x="", y="Between cell distance (µm)")+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  theme(legend.title = element_blank())


## cell covariance analysis-------
cc_cov_fun <- function(cell_trace){
  dat_cell_trace<- subset(cell_trace, cell_trace$Group!="Neutral")
  ID_cell <- unique(dat_cell_trace$ID)
  dat_cell_cov <- vector(mode = "list", length(ID_cell))
  
  for (i in seq_along(ID_cell)){
    dat_cell_trace_cov <- subset(dat_cell_trace, dat_cell_trace$ID==ID_cell[i]) %>%
      select("Time", "variable","value") %>%
      dcast(., Time~variable)
    if (ncol(dat_cell_trace_cov)<3){
      res_cov <- NA
    } else {
      res_cov <- cov(dat_cell_trace_cov[,-1])
      res_cov[res_cov==1]<- NA
    }
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

p_cov_d7<- ggcorrplot(dat_cell_cov_list[[7]][[1]],hc.order = TRUE, outline.col = "white", ggtheme = ggplot2::theme_void)+
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
  c_trim <- function(x){
    x[upper.tri(x)]<- NA
    c(x[!is.na(x)])
  }
  cov_value <- abs(unlist(mapply(c_trim, dat_cell_cov_list[[i]] )))
  dat_cov <- data.frame(Day=cell_area_day[i], value=cov_value)
  dat_cell_cov <- rbind(dat_cell_cov, dat_cov)
  
  ## for mean of cov
  cov_value_mean <- dat_cell_cov_list[[i]] %>%
    mapply(function(x) ifelse(is.null(ncol(x)), NA,mean(colMeans(abs(x), na.rm=T))), .)

  cov_value_max <- dat_cell_cov_list[[i]] %>%
    mapply(function(x) ifelse(is.null(ncol(x)),NA,mean(apply(x, 2, function(x) max(abs(x), na.rm = T)))), .)
  
  mouse_ID <- c(1: length(dat_cell_cov_list[[i]]))
  dat_cov_mean <- data.frame(Day = cell_area_day[i], ID= mouse_ID, value=cov_value_mean, value_max=cov_value_max)
  dat_cell_cov_sta <- rbind(dat_cell_cov_sta, dat_cov_mean)
}


## cummulative plot of matrix
dat_cell_cov$Day <- factor(dat_cell_cov$Day, levels = c("Pre", "Cond.", "Test"))
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

dat_cell_cov_sta1 <- melt(dat_cell_cov_sta, id.vars = c("ID", "Day")) %>%
  ddply(., .(ID, Day, variable), summarise, value=mean(value, na.rm = T)) %>%
  mutate(Day = factor(Day, levels = c("Pre", "Cond.", "Test")))
  

p_cov_mean<- subset(dat_cell_cov_sta1, dat_cell_cov_sta1$variable=="value") %>%
  ggplot(., aes(Day, value, fill=Day))+
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(position=position_jitterdodge(), shape=1, size=1)+
  scale_fill_manual(values=c("seagreen", "indianred", "deepskyblue4"))+
  labs(x="", y="Mean covariance (z-score)")+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  theme(legend.title = element_blank(), legend.position = "none")

p_cov_best <- subset(dat_cell_cov_sta1, dat_cell_cov_sta1$variable=="value_max") %>%
  ggplot(., aes(Day, value, fill=Day))+
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(position=position_jitterdodge(), shape=1, size=1)+
  scale_fill_manual(values=c("seagreen", "indianred", "deepskyblue4"))+
  labs(x="", y="Max covariance (z-score)")+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  theme(legend.title = element_blank(), legend.position = "none")

p_cov_bar <- plot_grid(p_cov_mean, p_cov_best, nrow = 1)

setwd("~cchen2/Documents/neuroscience/Pn\ project/Figure/PDF/")
cairo_pdf("p_cov_bar.pdf", width = 82/25.6, height = 65/25.6, family = "Arial")
p_cov_bar
dev.off()

## cell correlation analysis-----
cc_cor_fun <- function(cell_trace){
  dat_cell_trace<- subset(cell_trace, cell_trace$Group!="Neutral")
  ID_cell <- unique(dat_cell_trace$ID)
  dat_cell_cor <- vector(mode = "list", length(ID_cell))
  
  for (i in seq_along(ID_cell)){
    dat_cell_trace_cor <- subset(dat_cell_trace, dat_cell_trace$ID==ID_cell[i]) %>%
      select("Time", "variable","value") %>%
      dcast(., Time~variable)
    if (ncol(dat_cell_trace_cor)<3){
      res_cor <- NA
    } else {
      res_cor <- cor(dat_cell_trace_cor[,-1])
      res_cor[res_cor==1]<- NA
    }
    dat_cell_cor[[i]] <- res_cor
  }
  return(dat_cell_cor)
}

dat_cell_cor_list <- mapply(cc_cor_fun, dat_cell_trace_re, SIMPLIFY = F)


## cor matrix heatmap for m3
cor_range <- range(c(dat_cell_cor_list[[3]][[1]], dat_cell_cor_list[[6]][[1]], dat_cell_cor_list[[7]][[1]]))

p_cor_d3<- ggcorrplot(dat_cell_cor_list[[3]][[1]],hc.order = TRUE, outline.col = "white")+
  theme_void()+
  scale_fill_gradient2(limit = c(cor_range[1], cor_range[2]), low = "navy", high =  "red4", mid = "white")


p_cor_d6<- ggcorrplot(dat_cell_cor_list[[6]][[1]],hc.order = TRUE, outline.col = "white")+
  theme_void()+
  scale_fill_gradient2(limit = c(cor_range[1], cor_range[2]), low = "navy", high =  "red4", mid = "white")

p_cor_d7<- ggcorrplot(dat_cell_cor_list[[7]][[1]],hc.order = TRUE, outline.col = "white")+
  theme_void()+
  scale_fill_gradient2(limit = c(cor_range[1], cor_range[2]), low = "navy", high =  "red4", mid = "white")

p_cor_com <- plot_grid(p_cor_d3, p_cor_d6, p_cor_d7, nrow = 1)

setwd("~cchen2/Documents/neuroscience/Pn\ project/Figure/PDF/")
cairo_pdf("p_cor_com.pdf", width = 200/25.6, height = 65/25.6, family = "Arial")
p_cor_com
dev.off()


## compare the correlation between cells
dat_cell_cor <- c()
dat_cell_cor_sta <- c()
for (i in c(2,3,5,6,7)){
  c_trim <- function(x){
    x[upper.tri(x)]<- NA
    c(x[!is.na(x)])
  }
  cor_value <- abs(unlist(mapply(c_trim, dat_cell_cor_list[[i]] )))
  dat_cor <- data.frame(Day=cell_area_day[i], value=cor_value)
  dat_cell_cor <- rbind(dat_cell_cor, dat_cor)
  
  ## for mean of cor
  cor_value_mean <- dat_cell_cor_list[[i]] %>%
    mapply(function(x) ifelse(is.null(ncol(x)), NA,mean(colMeans(abs(x), na.rm=T))), .)
  
  cor_value_max <- dat_cell_cor_list[[i]] %>%
    mapply(function(x) ifelse(is.null(ncol(x)),NA,mean(apply(x, 2, function(x) max(abs(x), na.rm = T)))), .)
  
  mouse_ID <- c(1: length(dat_cell_cor_list[[i]]))
  dat_cor_mean <- data.frame(Day = cell_area_day[i], ID= mouse_ID, value=cor_value_mean, value_max=cor_value_max)
  dat_cell_cor_sta <- rbind(dat_cell_cor_sta, dat_cor_mean)
}


## cummulative plot of matrix
dat_cell_cor$Day <- factor(dat_cell_cor$Day, levels = c("Pre", "Cond.", "Test"))
p_cor_cum<- ggplot(dat_cell_cor, aes(value, group=Day, colour=Day))+
  stat_ecdf(geom = "step")+
  scale_colour_manual(values=c("seagreen", "indianred", "deepskyblue4"))+
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
cairo_pdf("p_cor_cum.pdf", width = 80/25.6, height = 65/25.6, family = "Arial")
p_cor_cum
dev.off()

## bar plot

dat_cell_cor_sta1 <- melt(dat_cell_cor_sta, id.vars = c("ID", "Day")) %>%
  ddply(., .(ID, Day, variable), summarise, value=mean(value, na.rm = T)) %>%
  mutate(Day = factor(Day, levels = c("Pre", "Cond.", "Test")))


p_cor_mean<- subset(dat_cell_cor_sta1, dat_cell_cor_sta1$variable=="value") %>%
  ggplot(., aes(Day, value, fill=Day))+
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(position=position_jitterdodge(), shape=1, size=1)+
  scale_fill_manual(values=c("seagreen", "indianred", "deepskyblue4"))+
  labs(x="", y="Mean corariance (z-score)")+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  theme(legend.title = element_blank(), legend.position = "none")

p_cor_best <- subset(dat_cell_cor_sta1, dat_cell_cor_sta1$variable=="value_max") %>%
  ggplot(., aes(Day, value, fill=Day))+
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(position=position_jitterdodge(), shape=1, size=1)+
  scale_fill_manual(values=c("seagreen", "indianred", "deepskyblue4"))+
  labs(x="", y="Max corariance (z-score)")+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  theme(legend.title = element_blank(), legend.position = "none")

p_cor_bar <- plot_grid(p_cor_mean, p_cor_best, nrow = 1)

setwd("~cchen2/Documents/neuroscience/Pn\ project/Figure/PDF/")
cairo_pdf("p_cor_bar.pdf", width = 90/25.6, height = 65/25.6, family = "Arial")
p_cor_bar
dev.off()


## for mice in homecage----
## Proof the homoestatsis activation of ACC-pn cells
## better to use the ca transit from homecage
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
  #stim_time<- seq(-5, 8.5, by=0.5)
  stim_time<- seq(-2, 6.5, by=0.5)
  ## number of rows to be binned
  n <- 10 # 0.05*10=0.5
  
  ## extract cell activity when they cross the border
  for (i in seq_along(t_stim)){
    t1_p <- t_stim[i]
    dat_stim1 <- dat_trace[(t1_p-40):(t1_p+140-1),]
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

n_trial <- 50
path_trace_m3 <- "~cchen2/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/matlab_analysis/m3/trace_days/m3_trace_d2.xlsx"
t_stim_m3 <- read.xlsx(path_trace_m3, colNames = F, rowNames = F) %>%
  ncol() %>%
  c(100:(.-200)) %>%
  sample(., n_trial)
dat_trace_m3 <- cc_miniscope_matlab(path_trace_m3, t_stim_m3)

path_trace_m7 <- "~cchen2/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/matlab_analysis/m7/trace_days/m7_trace_d3.xlsx"
t_stim_m7 <- read.xlsx(path_trace_m7, colNames = F, rowNames = F) %>%
  ncol() %>%
  c(100:(.-200)) %>%
  sample(., n_trial)
dat_trace_m7 <- cc_miniscope_matlab(path_trace_m7, t_stim_m7)


path_trace_m17 <- "~cchen2/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/matlab_analysis/m17/trace_days/m17_trace_d3.xlsx"
t_stim_m17 <- read.xlsx(path_trace_m17, colNames = F, rowNames = F) %>%
  ncol() %>%
  c(100:(.-200)) %>%
  sample(., n_trial)
dat_trace_m17 <- cc_miniscope_matlab(path_trace_m17, t_stim_m17)


path_trace_m18 <- "~cchen2/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/matlab_analysis/m18/trace_days/m18_trace_d3.xlsx"
t_stim_m18 <- read.xlsx(path_trace_m18, colNames = F, rowNames = F) %>%
  ncol() %>%
  c(100:(.-200)) %>%
  sample(., n_trial)
dat_trace_m18 <- cc_miniscope_matlab(path_trace_m18, t_stim_m18)

## for heat map plot
dat_cell_trace <- vector(mode = "list", n_trial)
for (i in 1:n_trial) {
  dat_cell_trace[[i]]<- cbind(dat_trace_m3[[i]], dat_trace_m7[[i]], dat_trace_m17[[i]], dat_trace_m18[[i]] )
  
}

mouse_ID <- c("m3", "m7", "m17", "m18")
dat_cell_trace_re <- vector(mode = "list", n_trial)
for (i in 1:length(dat_cell_trace)) {
  dat_cell_trace_d <- dat_cell_trace[[i]]
  cell_id <- sprintf("Cell%s",seq(1:ncol(dat_cell_trace_d)))
  colnames(dat_cell_trace_d)<- cell_id
  
  dat_cell_trace_d_cluster <- unname(kmeans(t(dat_cell_trace_d), centers = 3)[[1]])
  dat_cell_trace_d<- as.data.frame(dat_cell_trace_d)
  #stim_time<- seq(-5, 8.5, by=0.5)
  stim_time<- seq(-2, 6.5, by=0.5)
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
  dat_cell_trace_d_re$ID <- rep(mouse_ID, rep_time*length(stim_time))
  
  dat_cell_trace_re[[i]] <- dat_cell_trace_d_re
}

## heatmap plot
dat_trace <- dat_cell_trace_re[[1]]
dat_trace_sta <- ddply(dat_trace, .(variable, Group), summarise,mean=mean(value), sum=sum(value))
dat_trace_sta <- dat_trace_sta[order(dat_trace_sta[,'mean']),]
dat_trace$variable <- factor(dat_trace$variable, levels = dat_trace_sta$variable)

score_range <- range(dat_cell_trace_re[[1]]$value)
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
  theme(legend.position = 'none')

setwd("~cchen2/Documents/neuroscience/Pn\ project/Figure/PDF/")
cairo_pdf("p_heat_ctrl.pdf", width = 45/25.6, height = 62/25.6, family = "Arial")
p_heat
dev.off()

## average of trace by ID and group
num_trial_sta<- mapply(function(x) ddply(x, .(ID, Time, Group), summarise,mean=mean(value)), dat_cell_trace_re, SIMPLIFY = F) %>%
  mapply(nrow, .)

dat_trace_trial_sta<- mapply(function(x) ddply(x, .(ID, Time, Group), summarise,value1=mean(value)), dat_cell_trace_re, SIMPLIFY = F) %>%
  do.call(rbind, .) %>%
  mutate(., Trial=rep(sprintf("Trial%s",seq(1:n_trial)), num_trial_sta)) %>%
  subset(., .$Group!="Neutral") %>%
  ddply(., .(ID, Time, Trial), summarise,value2=sum(value1)) %>%
  ddply(., .(ID, Time), summarise,value2=mean(value2)) 

dat_trace_sta <- mapply(function(x) ddply(x, .(ID, Time, Group), summarise,value1=mean(value)), dat_cell_trace_re, SIMPLIFY = F) %>%
  do.call(rbind, .) %>%
  mutate(., Trial=rep(sprintf("Trial%s",seq(1:n_trial)), num_trial_sta)) %>%
  ddply(., .(ID,Time, Group), summarise,value2=mean(value1)) %>%
  ddply(., .(Time, Group), summarise,n=length(value2),mean=mean(value2),sd=sd(value2),se=sd(value2)/sqrt(length(value2)))
dat_trace_sta$Group <-  factor(dat_trace_sta$Group, levels = c("Neutral", "Excited", "Inhibited"))

p_trace <- ggplot(dat_trace_sta, aes(Time, mean, colour=Group))+
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
  theme(legend.position = 'none')
  #theme(legend.title = element_blank())
  

p_trace_sum <- ddply(dat_trace_trial_sta, .(Time), summarise,n=length(value2),mean=mean(value2),sd=sd(value2),se=sd(value2)/sqrt(length(value2))) %>%
  mutate(., xfill = ifelse(Time >= 0, Time, NA)) %>%
  ggplot(., aes(Time, mean))+
  geom_line()+
  geom_ribbon(aes(ymin=mean-se, ymax=mean+se), alpha=0.3, linetype=0)+
  geom_area(aes(x=xfill),alpha=0.2, fill="red"  )+
  labs(x="Time relative to crossing (s)", y="Z score")+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  geom_vline(xintercept = 0, col= 'red', linetype=2)


Time_area <- which(stim_time>0)

dat_anti_area <- tapply(dat_trace_trial_sta$value2, INDEX =dat_trace_trial_sta$ID, 
                                      function (x) AUC(stim_time[Time_area], x[Time_area])) %>%
  as.data.frame() %>%
  mutate(., Group="Ctrl")

colnames(dat_anti_area)[1] <- c( "Value")

p_EI_sum<- ggplot(dat_anti_area, aes(Group, Value))+
  geom_boxplot()+
  geom_jitter(width = 0.2, shape=1)+
  labs(x="", y="AUC (z-score")+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  scale_y_continuous(limits = c(2,4), expand = c(0,0))+
  theme(legend.title = element_blank())+
  theme(legend.position = "none")

p_ctrl_combine <- plot_grid(p_trace, p_trace_sum, p_EI_sum, nrow = 1, rel_widths = c(1,1,0.4))

setwd("~cchen2/Documents/neuroscience/Pn\ project/Figure/PDF/")
cairo_pdf("p_ctrl_com.pdf", width = 170/25.6, height = 62/25.6, family = "Arial")
p_ctrl_combine
dev.off()

## for pin prick-------
path_m3_pin <- "~cchen2/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/matlab_analysis/m3/pin_prick/m3_trace_pinprick.xlsx"
t_pin_m3 <- c(337,665,916,1303,1541,1767,2032,2277,2496,2824)

dat_pin_m3 <-cc_miniscope_matlab(path_m3_pin, t_pin_m3)

## for m7
path_m7_pin <- "~cchen2/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/matlab_analysis/m7/pin_prick/m7_trace_pinprick.xlsx"
t_pin_m7 <- c(493, 1323, 1867, 2415, 3035, 3556, 4079, 4519, 5020, 5517) 

dat_pin_m7 <-cc_miniscope_matlab(path_m7_pin, t_pin_m7)

dat_cell_trace <- vector(mode = "list", length = length(t_pin_m3))
for (i in 1:length(t_pin_m3)) {
  dat_cell_trace[[i]] <- cbind(dat_pin_m3[[i]], dat_pin_m7[[i]] )
  
}


dat_cell_trace_re <- vector(mode = "list", length = length(t_pin_m3))

# k-means clustering
mouse_ID_pin <- c("m3", "m7")
for (i in 1:length(dat_cell_trace)) {
  dat_cell_trace_d <- dat_cell_trace[[i]]
  cell_id <- sprintf("Cell%s",seq(1:ncol(dat_cell_trace_d)))
  colnames(dat_cell_trace_d)<- cell_id
  
  dat_cell_trace_d_cluster <- unname(kmeans(t(dat_cell_trace_d), centers = 3, nstart = 20)[[1]])
  stim_time<- seq(-2, 6.5, by=0.5)
  #stim_time<- seq(-5, 8.5, by=0.5)
  
  dat_cell_trace_d_re<- as.data.frame(dat_cell_trace_d) %>% 
    mutate(., Time= stim_time) %>% 
    melt(., id.vars='Time') %>% 
    mutate(., Group= rep(dat_cell_trace_d_cluster, each=length(stim_time)) )
  
  ## sort data by the value in each group
  dat_cell_d_sort <- as.numeric(names(sort(tapply(dat_cell_trace_d_re$value, dat_cell_trace_d_re$Group, mean), decreasing = T)))
  
  dat_cell_trace_d_re$Group[dat_cell_trace_d_re$Group == dat_cell_d_sort[1]] ="Excited"
  dat_cell_trace_d_re$Group[dat_cell_trace_d_re$Group == dat_cell_d_sort[2]] ="Neutral"
  dat_cell_trace_d_re$Group[dat_cell_trace_d_re$Group == dat_cell_d_sort[3]] ="Inhibited"
  
  rep_time <- c(ncol(dat_pin_m3[[i]]), ncol(dat_pin_m7[[i]]))
  dat_cell_trace_d_re <- mutate(dat_cell_trace_d_re, Group= factor(Group, levels = c("Excited", "Neutral", "Inhibited"))) %>% 
    mutate(., ID = rep(mouse_ID_pin, rep_time*length(stim_time)) )
  
  dat_cell_trace_re[[i]] <- dat_cell_trace_d_re
}

## plot trace by group
dat_cell_pin_sta <- mapply(function (x) ddply(x, .(ID,Time, Group), summarise,value=mean(value)), dat_cell_trace_re, SIMPLIFY = F) %>% 
  do.call(rbind, .) %>% 
  ddply(., .(ID,Time, Group), summarise,value=mean(value)) %>% 
  mutate(Group = factor(Group, levels = c("Neutral", "Excited", "Inhibited")))




## group by day
p_trace <- ddply(dat_cell_pin_sta, .(Time, Group), summarise,n=length(value),mean=mean(value),sd=sd(value),se=sd(value)/sqrt(length(value))) %>% 
  ggplot(., aes(Time, mean, colour=Group))+
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
  theme(legend.title = element_blank(), legend.position = 'none')



## calculate the sum of active and inhibited trace
p_trace_sum <- filter(dat_cell_pin_sta, Group!="Neutral") %>%
  ddply(., .(ID,Time), summarise, value=sum(value)) %>%
  ddply(., .(Time), summarise,n=length(value),mean=mean(value),sd=sd(value),se=sd(value)/sqrt(length(value))) %>% 
  ggplot(., aes(Time, mean))+
  geom_line()+
  geom_ribbon(aes(ymin=mean-se, ymax=mean+se), alpha=0.1, linetype=0)+
  labs(x="Time relative to crossing (s)", y="AUC (z score)")+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  geom_vline(xintercept = 0, col= 'red', linetype=2)+
  theme(legend.title = element_blank(), legend.position = 'top')


## for heat pain----
path_m3_har <- "~cchen2/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/matlab_analysis/Hargreave/m3_trace_har.xlsx"
t_har_m3 <- c(1171, 2159, 3579, 4708, 5693, 6910)

dat_har_m3 <-cc_miniscope_matlab(path_m3_har, t_har_m3)

## for m7
path_m7_har <- "~cchen2/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/matlab_analysis/Hargreave/m7_trace_har.xlsx"
t_har_m7 <- c(756, 1917, 2947, 3541, 4402, 5450, 6937)

dat_har_m7 <-cc_miniscope_matlab(path_m7_har, t_har_m7)

dat_cell_trace <- vector(mode = "list", length = length(t_har_m3))
for (i in 1:length(t_har_m3)) {
  dat_cell_trace[[i]] <- cbind(dat_har_m3[[i]], dat_har_m7[[i]] )
  
}


dat_cell_trace_re <- vector(mode = "list", length = length(t_har_m3))

# k-menas clustering
mouse_ID_har <- c("m3", "m7")
for (i in 1:length(dat_cell_trace)) {
  dat_cell_trace_d <- dat_cell_trace[[i]]
  cell_id <- sprintf("Cell%s",seq(1:ncol(dat_cell_trace_d)))
  colnames(dat_cell_trace_d)<- cell_id
  
  dat_cell_trace_d_cluster <- unname(kmeans(t(dat_cell_trace_d), centers = 3, nstart = 20)[[1]])
  stim_time<- seq(-2, 6.5, by=0.5)
  #stim_time<- seq(-5, 8.5, by=0.5)
  
  dat_cell_trace_d_re<- as.data.frame(dat_cell_trace_d) %>% 
    mutate(., Time= stim_time) %>% 
    melt(., id.vars='Time') %>% 
    mutate(., Group= rep(dat_cell_trace_d_cluster, each=length(stim_time)) )
  
  ## sort data by the value in each group
  dat_cell_d_sort <- as.numeric(names(sort(tapply(dat_cell_trace_d_re$value, dat_cell_trace_d_re$Group, mean), decreasing = T)))
  
  dat_cell_trace_d_re$Group[dat_cell_trace_d_re$Group == dat_cell_d_sort[1]] ="Excited"
  dat_cell_trace_d_re$Group[dat_cell_trace_d_re$Group == dat_cell_d_sort[2]] ="Neutral"
  dat_cell_trace_d_re$Group[dat_cell_trace_d_re$Group == dat_cell_d_sort[3]] ="Inhibited"
  
  rep_time <- c(ncol(dat_har_m3[[i]]), ncol(dat_har_m7[[i]]))
  dat_cell_trace_d_re <- mutate(dat_cell_trace_d_re, Group= factor(Group, levels = c("Excited", "Neutral", "Inhibited"))) %>% 
    mutate(., ID = rep(mouse_ID_har, rep_time*length(stim_time)) )
  
  dat_cell_trace_re[[i]] <- dat_cell_trace_d_re
}

## plot trace by group
dat_cell_har_sta <- mapply(function (x) ddply(x, .(ID,Time, Group), summarise,value=mean(value)), dat_cell_trace_re, SIMPLIFY = F) %>% 
  do.call(rbind, .) %>% 
  ddply(., .(ID,Time, Group), summarise,value=mean(value)) %>% 
  mutate(Group = factor(Group, levels = c("Neutral", "Excited", "Inhibited")))




## group by day
p_trace <- ddply(dat_cell_har_sta, .(Time, Group), summarise,n=length(value),mean=mean(value),sd=sd(value),se=sd(value)/sqrt(length(value))) %>% 
  ggplot(., aes(Time, mean, colour=Group))+
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
  theme(legend.title = element_blank(), legend.position = 'none')



## calculate the sum of active and inhibited trace

p_trace_sum <- filter(dat_cell_har_sta, Group!="Neutral") %>%
  ddply(., .(ID,Time), summarise, value=sum(value)) %>%
  ddply(., .(Time), summarise,n=length(value),mean=mean(value),sd=sd(value),se=sd(value)/sqrt(length(value))) %>% 
  ggplot(., aes(Time, mean))+
  geom_line()+
  geom_ribbon(aes(ymin=mean-se, ymax=mean+se), alpha=0.1, linetype=0)+
  labs(x="Time relative to crossing (s)", y="AUC (z score)")+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  geom_vline(xintercept = 0, col= 'red', linetype=2)+
  theme(legend.title = element_blank(), legend.position = 'top')
