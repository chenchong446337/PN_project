## Project note------
## The script are used to analyze the data from miniscope experiemnt analyzed by matlab
## updated: 01172021, only compare the pre-test and test

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
  # dat_trace[dat_trace < 0] <- 0
  
  
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
path_trace_m3 <- as.list(list.files(path ="~cchen/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/matlab_analysis/m3/trace_days/", pattern = "*.xlsx", full.names = T ))[c(3,7)]
t_stim_m3 <- list(t_stim_m3_d3 = c(1885), t_stim_m3_d7 = c(346))
dat_trace_m3 <- mapply(c_miniscope_matlab, path_trace_m3, t_stim_m3, SIMPLIFY = F)

# for m7
path_trace_m7 <- as.list(list.files(path ="~cchen/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/matlab_analysis/m7/trace_days/", pattern = "*.xlsx", full.names = T ))[c(3,7)]
t_stim_m7 <- list(t_stim_m7_d3 = c(360), t_stim_m7_d7 = c(236))
dat_trace_m7 <- mapply(c_miniscope_matlab, path_trace_m7, t_stim_m7,SIMPLIFY = F)

# for m17
path_trace_m17 <- as.list(list.files(path ="~cchen/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/matlab_analysis/m17/trace_days/", pattern = "*.xlsx", full.names = T ))[c(3,7)]
t_stim_m17 <- list(t_stim_m17_d3 = c(437), t_stim_m17_d7 = c(157) )
dat_trace_m17 <- mapply(c_miniscope_matlab, path_trace_m17, t_stim_m17, SIMPLIFY = F)

# for m18
path_trace_m18 <- as.list(list.files(path ="~cchen/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/matlab_analysis/m18/trace_days/", pattern = "*.xlsx", full.names = T ))[c(3,7)]
t_stim_m18 <- list(t_stim_m18_d3 = c(784), t_stim_m18_d7 = c(124))
dat_trace_m18 <- mapply(c_miniscope_matlab, path_trace_m18, t_stim_m18, SIMPLIFY = F)

# for m855
path_trace_m855 <- as.list(list.files(path ="~cchen/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/matlab_analysis/m855/trace_days", pattern = "*.xlsx", full.names = T ))[c(3,7)]
t_stim_m855 <- list(t_stim_m855_d3 = c(392) *2, t_stim_m855_d7 = c(41)*2)
dat_trace_m855 <- mapply(c_miniscope_matlab, path_trace_m855, t_stim_m855, SIMPLIFY = F)

# for m857
path_trace_m857 <- as.list(list.files(path ="~cchen/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/matlab_analysis/m857/trace_days", pattern = "*.xlsx", full.names = T ))[c(1, 5)]
t_stim_m857 <- list(t_stim_m857_d3 = c(847) *2, t_stim_m857_d7 = c(167)*2)
dat_trace_m857 <- mapply(c_miniscope_matlab, path_trace_m857, t_stim_m857, SIMPLIFY = F)




## combine data and do k-means analysis-----

dat_cell_trace <- vector(mode = "list", 2)
for (i in 1:2) {
  dat_cell_trace[[i]] <- cbind(dat_trace_m3[[i]], dat_trace_m7[[i]], dat_trace_m17[[i]], dat_trace_m18[[i]], dat_trace_m855[[i]], dat_trace_m857[[i]])
  
}

## overview of the data
mouse_ID <- c("m3", "m7", "m17", "m18", "m855", "m857")
n_cells <- c(ncol(dat_trace_m3[[1]]), ncol(dat_trace_m7[[1]]), ncol(dat_trace_m17[[1]]), ncol(dat_trace_m18[[1]]), ncol(dat_trace_m855[[1]]), ncol(dat_trace_m857[[1]]))
mean(n_cells)
se=sd(n_cells)/sqrt(length(n_cells))

p_number <- tibble(num = n_cells, Group = "Num") %>% 
  ggplot(., aes(x = Group, y = num))+
  geom_boxplot()+
  geom_jitter(width = 0.2, shape=1)+
  labs(x="", y="No. of cells")+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  scale_y_continuous(limits = c(0, 100))+
  theme(legend.position = "none")

setwd("~cchen/Documents/neuroscience/Pn\ project/Figure/PDF/")
cairo_pdf("p_number.pdf", width = 30/25.6, height = 60/25.6, family = "Arial")
p_number
dev.off()



dat_cell_trace_re <- vector(mode = "list", 2)
for (i in 1:length(dat_cell_trace)) {
  dat_cell_trace_d <- dat_cell_trace[[i]]
  cell_id <- sprintf("Cell%s",seq(1:ncol(dat_cell_trace_d)))
  colnames(dat_cell_trace_d)<- cell_id
  
  dat_cell_trace_d_cluster <- unname(kmeans(t(dat_cell_trace_d), centers = 3)[[1]])
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

## set the range of z score during all days
# score_range <- rbind(dat_cell_trace_re[[1]], dat_cell_trace_re[[2]]) %>% 
#   .$value %>% 
#   range()

  
heat_m_ID <- "m857"
score_rang1 <- dat_cell_trace_re[[1]] %>% 
  filter(ID == heat_m_ID)

score_range <- dat_cell_trace_re[[2]] %>% 
  filter(ID ==heat_m_ID) %>% 
  rbind(., score_rang1) %>% 
  .$value %>% 
  range()
  
group_day <- c("Pre", "Test")

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
    labs(x="Time relative to crossing (s)", y="Number of cells", title  = heat_m_ID)+
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
p_heat_com <- plot_grid(p_heat_Pre, p_heat_Test, nrow = 1)

setwd("~cchen/Documents/neuroscience/Pn\ project/Figure/PDF/")
cairo_pdf("p_anti_heat.pdf", width = 80/25.6, height = 55/25.6, family = "Arial")
p_heat_com
dev.off()

## plot trace by group-----
dat_cell_trace_pre_sta <-  ddply(dat_cell_trace_re[[1]], .(ID,Time, Group), summarise,value=mean(value)) 
dat_cell_trace_test_sta <- ddply(dat_cell_trace_re[[2]], .(ID,Time, Group), summarise,value=mean(value))

dat_cell_trace_sta <- ddply(dat_cell_trace_pre_sta, .(Time, Group), summarise,n=length(value),mean=mean(value),sd=sd(value),se=sd(value)/sqrt(length(value))) %>%
  rbind(., ddply(dat_cell_trace_test_sta, .(Time, Group), summarise,n=length(value),mean=mean(value),sd=sd(value),se=sd(value)/sqrt(length(value)))) %>%
  mutate(., Day = rep(c("Pre",  "Test"),each= length(stim_time)*3)) %>% 
  mutate(Day = factor(Day, levels = c("Pre",  "Test")), Group = factor(Group,levels = c("Neutral", "Excited", "Inhibited") ))

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

setwd("~cchen/Documents/neuroscience/Pn\ project/Figure/PDF/")
cairo_pdf("p_trace.pdf", width = 90/25.6, height = 62/25.6, family = "Arial")
p_trace
dev.off()


## calculate the sum of active and inhibited trace----
dat_cell_trace_sum <- rbind(dat_cell_trace_pre_sta, dat_cell_trace_test_sta) %>%
  as_tibble() %>% 
  mutate(Day = rep(c("Pre",  "Test"), c(nrow(dat_cell_trace_pre_sta),nrow(dat_cell_trace_test_sta)))) %>%
  #subset(., .$Group!="Neutral") %>%
  ddply(., .(ID,Time, Day), summarise, value=sum(value)) %>%
  ddply(., .(Time, Day), summarise,n=length(value),mean=mean(value),sd=sd(value),se=sd(value)/sqrt(length(value))) %>%
  mutate(Day = factor(Day, levels = c("Pre",  "Test")))

p_trace_sum <- ggplot(dat_cell_trace_sum, aes(Time, mean, colour=Day))+
  geom_line()+
  geom_ribbon(aes(ymin=mean-se, ymax=mean+se, fill=Day), alpha=0.2, linetype=0)+
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
  geom_vline(xintercept = 0, col= 'red', linetype=2)+
  theme(legend.title = element_blank(), legend.position = c(0.2, 0.8))

## only plot m18
dat_cell_trace_pre_sta_m18 <-  ddply(dat_cell_trace_re[[1]], .(ID,Time, Group), summarise,mean_value=mean(value), sd=sd(value),se=sd(value)/sqrt(length(value))) %>%
  filter(ID == "m855") %>% 
  mutate(Day ="Pre")

dat_cell_trace_test_sta_m18 <- ddply(dat_cell_trace_re[[2]], .(ID,Time, Group), summarise,mean_value=mean(value), sd=sd(value),se=sd(value)/sqrt(length(value))) %>%
  filter(ID == "m855") %>% 
  mutate(Day ="Test")


p_trace_m18 <- rbind(dat_cell_trace_pre_sta_m18, dat_cell_trace_test_sta_m18) %>% 
  mutate(Day = factor(Day, levels = c("Pre",  "Test")), Group = factor(Group,levels = c( "Excited","Neutral", "Inhibited") )) %>% 
  ggplot(., aes(x = Time, y = mean_value, group = Group,colour= Day))+
  geom_line()+
  facet_grid(cols = vars(Day))+
  geom_ribbon(aes(ymin=mean_value-se, ymax=mean_value+se, fill=Day), alpha=0.1, linetype=0)+
  geom_line()+
  geom_ribbon(aes(ymin=mean_value-se, ymax=mean_value+se, fill=Day), alpha=0.2, linetype=0)+
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
  geom_vline(xintercept = 0, col= 'red', linetype=2)+
  theme(legend.title = element_blank(), legend.position = c(0.2, 0.8))

## for the sum
dat_cell_trace_sum <- rbind(dat_cell_trace_pre_sta, dat_cell_trace_test_sta) %>%
  as_tibble() %>% 
  mutate(Day = rep(c("Pre",  "Test"), c(nrow(dat_cell_trace_pre_sta),nrow(dat_cell_trace_test_sta)))) %>%
  #subset(., .$Group!="Neutral") %>%
  filter(ID =="m855") %>% 
  ddply(., .(Time, Day), summarise,n=length(value),mean=mean(value),sd=sd(value),se=sd(value)/sqrt(length(value))) %>%
  mutate(Day = factor(Day, levels = c("Pre",  "Test")))

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
cairo_pdf("p_trace_m18.pdf", width = 90/25.6, height = 60/25.6, family = "Arial")
p_trace_m18
dev.off()

setwd("~cchen/Documents/neuroscience/Pn\ project/Figure/PDF/")
cairo_pdf("p_trace_sum.pdf", width = 75/25.6, height = 60/25.6, family = "Arial")
p_trace_sum
dev.off()


## for EI change analysis-----
dat_cell_area <- c()
group_day <- c("Pre",  "Test")
for (i in c(1,2)){
  rep_time <- c(ncol(dat_trace_m3[[i]]), ncol(dat_trace_m7[[i]]), ncol(dat_trace_m17[[i]]), ncol(dat_trace_m18[[i]]), ncol(dat_trace_m855[[i]]), ncol(dat_trace_m857[[i]]))
  dat_trace <- dat_cell_trace_re[[i]]
  dat_trace$ID <- rep(mouse_ID, rep_time*length(stim_time))
  dat_trace_sta <- ddply(dat_trace, .(ID, Time, Group), summarise,n=length(value),mean=mean(value),sd=sd(value),se=sd(value)/sqrt(length(value)))
  Time_area <- which(stim_time>=0)
  dat_anti_area <- as.data.frame(tapply(dat_trace_sta$mean, INDEX = list(dat_trace_sta$ID, dat_trace_sta$Group), 
                                        function (x) AUC(stim_time[Time_area], x[Time_area])))
  dat_anti_area$ID <- rownames(dat_anti_area)
  rownames(dat_anti_area)<-NULL
  dat_anti_area$Day <- group_day[i]
  dat_anti_area$ratio <- abs(dat_anti_area$Excited/dat_anti_area$Inhibited)
  dat_anti_area$sum <- dat_anti_area$Excited + dat_anti_area$Inhibited
  dat_cell_area <- rbind(dat_cell_area, dat_anti_area)
}

dat_cell_area$Day <- factor(dat_cell_area$Day, levels = c("Pre","Test"))

p_EI_ratio <- dat_cell_area %>% 
  select(ID,Day, ratio) %>% 
  ggplot(., aes(Day, ratio, color=Day))+
  geom_boxplot(outlier.shape = NA)+
  geom_line(aes(group = ID), color="gray90")+
  geom_jitter(width = 0.2, shape=1)+
  scale_colour_manual(values=c("deepskyblue4", "indianred"))+
  labs(x="", y="E/I ratio")+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  scale_y_continuous(limits = c(0, 4), expand = c(0,0))+
  theme(legend.title = element_blank(), legend.position = "none")
#annotate(x=c(1,1,2,2), y=c(3.4,3.5,3.5,3.4),"path")+
#annotate("text",x=1.5,y=3.5, label="***", size=5)


## statistic test
t_EI<-dat_cell_area %>% 
  select(Day, ID, ratio) %>% 
  wilcox.test(ratio~Day, ., paired = TRUE)


setwd("~cchen/Documents/neuroscience/Pn\ project/Figure/PDF/")
cairo_pdf("p_EI_ratio.pdf", width = 40/25.6, height = 65/25.6, family = "Arial")
p_EI_ratio
dev.off() 
## plot the sum of E-I----

dat_EI_sum_sta <- dat_cell_area %>% 
  select(Day, ID, sum) %>% 
  ddply(., .(Day), summarise,n=length(sum),mean=mean(sum),sd=sd(sum),se=sd(sum)/sqrt(length(sum)))

p_sum <- dat_cell_area %>% 
  select(Day, ID, sum) %>% 
  wilcox.test(sum~Day, .)

p_EI_sum<- dat_cell_area %>% 
  select(Day, ID, sum) %>% 
  ggplot(., aes(Day, sum, colour=Day))+
  geom_boxplot(outlier.shape = NA)+
  geom_line(aes(group=ID), colour="gray90")+
  geom_jitter(width = 0.2, shape=1)+
  scale_colour_manual(values=c("deepskyblue4", "indianred" ))+
  labs(x="", y="AUC (z-score)")+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  theme(legend.title = element_blank(), legend.position = "none")
#annotate(x=c(1,1,2,2), y=c(9,9.1,9.1,9),"path")+
#annotate("text",x=1.5,y=9.1, label="*", size=5)


## statistic test
t_sum_EI<-dat_cell_area %>% 
  select(Day, ID, sum) %>% 
  t.test(sum~Day, ., paired=T)

p_trace_sum_com <- plot_grid(p_EI_ratio, p_EI_sum, nrow = 1) %>% 
  plot_grid(p_trace_sum, ., nrow = 2)

setwd("~cchen/Documents/neuroscience/Pn\ project/Figure/PDF/")
cairo_pdf("p_trace_sum_com.pdf", width = 75/25.6, height = 110/25.6, family = "Arial")
p_trace_sum_com
dev.off()

## show the change of E, I and sum
p_EIS_change <- dat_cell_area %>%
  select("Excited", "Inhibited","ID","Day", "sum") %>%
  melt(., id.vars= c('Day', "ID")) %>%
  ddply(., .(ID, Day, variable), summarise, value=mean(value, na.rm = T)) %>%
  ggplot(., aes(interaction(Day, variable), value, fill=Day))+
  geom_boxplot(outlier.shape = NA)+
  geom_line(aes(group=interaction(ID, variable)), colour="gray90")+
  geom_point(aes(colour = Day,shape = Day),position=position_jitterdodge(jitter.width = .2), size =2)+
  scale_fill_manual(values = c("white", "white"))+
  scale_shape_manual(values=c(20, 18))+
  scale_colour_manual(values=c("deepskyblue4", "indianred"))+
  labs(x="", y="AUC (z-score)")+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  theme(legend.title = element_blank(), legend.position = "none")

p_EIS_change <- dat_cell_area %>%
  select("Excited", "Inhibited","ID","Day", "sum") %>%
  melt(., id.vars= c('Day', "ID")) %>% 
  aov(value~Day,.)

summary(p_EIS_change)

setwd("~cchen/Documents/neuroscience/Pn\ project/Figure/PDF/")
cairo_pdf("p_EIS_change.pdf", width = 80/25.6, height = 65/25.6, family = "Arial")
p_EIS_change
dev.off()

## portion of cell catlog------

dat_cell_cat <- mapply(function (x) subset(x, x$Time==0), dat_cell_trace_re, SIMPLIFY = F ) 



p_cat<- mapply(function(x) prop.table(table(x$ID, x$Group), 1), dat_cell_cat, SIMPLIFY = F) %>%
  do.call(rbind, .) %>%
  as.data.frame() %>%
  mutate(ID = rep(rownames(.)[1:length(mouse_ID)], 2), Day=rep(c("Pre", "Test"), each=length(mouse_ID))) %>%
  melt(., id.vars = c("ID", "Day")) %>%
  ddply(., .(ID, Day, variable), summarise, "Prop"=mean(value)) %>%
  mutate(Day = factor(Day, levels = c("Pre", "Test"))) %>%
  mutate(variable=factor(variable, levels = c("Excited", "Neutral", "Inhibited"))) %>% 
  ggplot(., aes(interaction(Day,variable), Prop, fill=Day))+
  geom_boxplot(outlier.shape = NA)+
  geom_line(aes(group=interaction(ID, variable)), colour="gray90")+
  geom_point(aes(colour = Day,shape = Day),position=position_jitterdodge(jitter.width = .2), size =2)+
  scale_fill_manual(values = c("white", "white"))+
  scale_shape_manual(values=c(20, 18))+
  scale_colour_manual(values=c( "deepskyblue4", "indianred"))+
  labs(x="", y="% of all neurons")+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  scale_y_continuous(limits = c(0,1),expand = c(0,0))+
  theme(legend.position = c(0.15, 0.8), legend.title = element_blank())

t_excited<- mapply(function(x) prop.table(table(x$ID, x$Group), 1), dat_cell_cat, SIMPLIFY = F) %>%
  do.call(rbind, .) %>%
  as.data.frame() %>%
  mutate(ID = rep(rownames(.)[1:length(mouse_ID)], 2), Day=rep(c("Pre", "Test"), each=length(mouse_ID))) %>%
  melt(., id.vars = c("ID", "Day")) %>%
  ddply(., .(ID, Day, variable), summarise, "Prop"=mean(value)) %>% 
  ddply(., .(Day, variable), summarise, mean=mean(Prop),sd=sd(Prop),se=sd(Prop)/sqrt(length(Prop)))
  #aov(Prop~Day + variable, .)

summary(t_excited)

t_inhibited<- mapply(function(x) prop.table(table(x$ID, x$Group), 1), dat_cell_cat, SIMPLIFY = F) %>%
  do.call(rbind, .) %>%
  as.data.frame() %>%
  mutate(ID = rep(rownames(.)[1:length(mouse_ID)], 2), Day=rep(c("Pre", "Test"), each=length(mouse_ID))) %>%
  melt(., id.vars = c("ID", "Day")) %>%
  ddply(., .(ID, Day, variable), summarise, "Prop"=mean(value)) %>% 
  filter(variable =="Inhibited") %>% 
  wilcox.test(Prop~Day,paired = T, .)



setwd("~cchen/Documents/neuroscience/Pn\ project/Figure/PDF/")
cairo_pdf("p_cat.pdf", width = 80/25.6, height = 70/25.6, family = "Arial")
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
  
  ctrl_rang <- range(t_stim-42, t_stim-2)
  test_rang <- range(t_stim, t_stim + 40)
  
  ctrl_freq <- length(dat_peak[dat_peak>= ctrl_rang[1] & dat_peak < ctrl_rang[2]])/2 # 2s before crossing
  test_freq <- length(dat_peak[dat_peak>= test_rang[1] & dat_peak < test_rang[2]])/2 # 4s after crossing
  diff_freq <- test_freq - ctrl_freq
  
  return(c(ctrl_freq, test_freq, diff_freq))
}

## for m3
path_peak_m3 <- as.list(list.files(path ="~cchen/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/matlab_analysis/m3/peaks_days", pattern = "*.csv", full.names = T ))[c(3, 7)]
dat_rate_m3 <- t(mapply(cc_firing_fun, path_trace_m3, path_peak_m3, t_stim_m3))

## for m7
path_peak_m7 <- as.list(list.files(path ="~cchen/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/matlab_analysis/m7/peaks_days/", pattern = "*.csv", full.names = T ))[c(3, 7)]
dat_rate_m7 <- t(mapply(cc_firing_fun, path_trace_m7, path_peak_m7, t_stim_m7))

## for m17
path_peak_m17 <- as.list(list.files(path ="~cchen/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/matlab_analysis/m17/peaks_days/", pattern = "*.csv", full.names = T ))[c(3, 7)]
dat_rate_m17 <- t(mapply(cc_firing_fun, path_trace_m17, path_peak_m17, t_stim_m17))

## for m18
path_peak_m18 <- as.list(list.files(path ="~cchen/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/matlab_analysis/m18/peaks_days/", pattern = "*.csv", full.names = T ))[c(3, 7)]
dat_rate_m18 <- t(mapply(cc_firing_fun, path_trace_m18, path_peak_m18, t_stim_m18))

## for m855
path_peak_m855 <- as.list(list.files(path ="~cchen/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/matlab_analysis/m855/peaks_days/", pattern = "*.csv", full.names = T ))[c(3, 7)]
dat_rate_m855 <- t(mapply(cc_firing_fun, path_trace_m855, path_peak_m855, t_stim_m855))

## for m857
path_peak_m857 <- as.list(list.files(path ="~cchen/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/matlab_analysis/m857/peaks_days/", pattern = "*.csv", full.names = T ))
dat_rate_m857 <- t(mapply(cc_firing_fun, path_trace_m857, path_peak_m857, t_stim_m857))


## combine data and plot
mouse_ID <- c("m3", "m7", "m17", "m18", "m855", "m857")
group_day <- c("Pre", "Test")
dat_firing <- rbind(dat_rate_m3, dat_rate_m7, dat_rate_m17, dat_rate_m18, dat_rate_m855,dat_rate_m857 ) %>% 
  as_tibble() %>% 
  rename(Freq_before = V1, Freq_after = V2, diff = V3) %>% 
  mutate(ID = rep(mouse_ID, each=2)) %>% 
  mutate(Group = rep(group_day, length(mouse_ID))) %>% 
  gather(variable, value, -ID,  -Group) %>% 
  ddply(., .(ID, Group, variable), summarise, value=mean(value)) %>% 
  mutate(Group= factor(Group, levels = c("Pre",  "Test")))

p_firing <- dat_firing %>% 
  filter(variable=="diff") %>% 
  ggplot(., aes(Group, value, group=Group))+
  geom_boxplot(outlier.shape = NA)+
  geom_line(aes(group=ID), colour="gray90")+
  geom_jitter(aes(colour = Group, shape = Group),width = 0.2,  size=2)+
  scale_shape_manual(values=c(20, 18))+
  scale_color_manual(values=c("deepskyblue4", "indianred"))+
  labs(x="", y="Diff. of firing frequency (Hz)")+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  scale_y_continuous(limits = c(-4,4))+
  theme(legend.position = 'none')

p_test <- dat_firing %>% 
  filter(variable=="diff") %>% 
  wilcox.test(value~Group, ., paired = TRUE)


setwd("~cchen/Documents/neuroscience/Pn\ project/Figure/PDF/")
cairo_pdf("p_firing.pdf", width = 40/25.6, height = 65/25.6, family = "Arial")
p_firing
dev.off() 


## Global ID only plot day3 and day7 (08092020)-----

cc_globalID_fun <- function(path_ID, path_trace, t_stim){
  Global_ID <- read.csv(path_ID, header = F) %>% 
    na_if(., 0) %>% 
    drop_na()
  
  
  ## read.xlsx with no header and remove the first col
  cc_read.xlsx<- function(x){
    dat <- read.xlsx(x, colNames = F) %>% 
      select(., -1)
    dat
  }
  
  ## do z score of the whole trace
  dat_trace <- mapply(cc_read.xlsx, path_trace,SIMPLIFY = F) %>% 
    mapply(function(x) apply(x, 2, scale), ., SIMPLIFY = F)
  
  stim_time<- seq(-2, 6.5, by=0.5)
  ## number of rows to be binned
  n <- 10 # 0.05*10=0.5
  
  ## create a 2 days list and put the activation of each cell in it
  cc_trace_pick <- function(dat_trace, cell_pick){
    dat_trace_pick <- t(dat_trace[cell_pick,])
    dat_trace_pick
  }
  
  dat_cell_trace_day <- mapply(cc_trace_pick, dat_trace, as.list(Global_ID))
  
  cc_trace_extract <- function(dat_cell_trace, t_stim_day){
    dat_stim <- vector(mode = "list", length = length(t_stim_day))
    for (i in seq_along(t_stim_day)){
      t1_p <- t_stim_day[i]
      dat_stim1 <- dat_cell_trace[(t1_p-40):(t1_p+140-1),]
      dat_stim1 <- aggregate(dat_stim1,list(rep(1:(nrow(dat_stim1)%/%n+1),each=n,len=nrow(dat_stim1))),mean)[-1]
      dat_stim[[i]] <- dat_stim1[1:5,] %>%  ## baseline as -2 to 0
        colMeans(., na.rm = T) %>% 
        sweep(dat_stim1, 2, ., FUN = "-")
    }
    ## average the trace by cell number
    if (length(t_stim_day)>1){
      dat_cell_trace_average <- aaply(laply(dat_stim, as.matrix), c(2, 3), function(x) mean(x, na.rm=T))
      
    } else {
      dat_cell_trace_average <- data.matrix(dat_stim[[1]])
      
    }
    return(dat_cell_trace_average)
  }
  dat_cell_trace_day_extract <- mapply(cc_trace_extract, dat_cell_trace_day, t_stim, SIMPLIFY = F)
  
  return(dat_cell_trace_day_extract)
}

## for m3
global_ID_m3 <- "~cchen/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/matlab_analysis/m3/Global_id/m3_anti_pin_global_ID_37.csv"
dat_global_trace_m3 <- cc_globalID_fun(global_ID_m3, path_trace_m3, t_stim_m3)

## for m7
global_ID_m7 <- "~cchen/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/matlab_analysis/m7/Global_id/m7_anti_pin_global_ID_37.csv"
dat_global_trace_m7 <- cc_globalID_fun(global_ID_m7, path_trace_m7, t_stim_m7)

## for m17
global_ID_m17 <- "~cchen/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/matlab_analysis/m17/Global_ID/m17_global_ID_37.csv"
dat_global_trace_m17 <- cc_globalID_fun(global_ID_m17, path_trace_m17, t_stim_m17)

## for m18
global_ID_m18 <- "~cchen/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/matlab_analysis/m18/Global_ID/m18_global_ID_37.csv"
dat_global_trace_m18 <- cc_globalID_fun(global_ID_m18, path_trace_m18, t_stim_m18)

## for m855
global_ID_m855 <- "~cchen/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/matlab_analysis/m855/Global_ID/m855_global_ID_37.csv"
dat_global_trace_m855 <- cc_globalID_fun(global_ID_m855, path_trace_m855, t_stim_m855)

## for m857
global_ID_m857 <- "~cchen/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/matlab_analysis/m857/Global_ID/m857_global_ID_37.csv"
dat_global_trace_m857 <- cc_globalID_fun(global_ID_m857, path_trace_m857, t_stim_m857)

## combin data
dat_global_trace_d3 <- cbind(dat_global_trace_m3[[1]], dat_global_trace_m7[[1]], dat_global_trace_m17[[1]], dat_global_trace_m18[[1]], dat_global_trace_m855[[1]], dat_global_trace_m857[[1]])
colnames(dat_global_trace_d3) <- str_c("Cell", 1:ncol(dat_global_trace_d3))

dat_global_trace_d7 <- cbind(dat_global_trace_m3[[2]], dat_global_trace_m7[[2]], dat_global_trace_m17[[2]], dat_global_trace_m18[[2]], dat_global_trace_m855[[2]], dat_global_trace_m857[[2]])
colnames(dat_global_trace_d7) <- str_c("Cell", 1:ncol(dat_global_trace_d7))

dat_cell_trace_global <- vector(mode = "list", 2)
dat_cell_trace_global[[1]] <- dat_global_trace_d3
dat_cell_trace_global[[2]] <- dat_global_trace_d7


dat_cell_trace_global_re <- vector(mode = "list", 2)

# k-menas clustering
for (i in 1:length(dat_cell_trace_global)) {
  dat_cell_trace_d <- dat_cell_trace_global[[i]]
  
  
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
  
  dat_cell_trace_d_re <- mutate(dat_cell_trace_d_re, Group= factor(Group, levels = c("Excited", "Neutral", "Inhibited")))
  dat_cell_trace_global_re[[i]] <- dat_cell_trace_d_re
}

score_range <- do.call(rbind, dat_cell_trace_global_re) %>% 
  .$value %>% 
  range()

## align by cells
dat_trace_sta <- dat_cell_trace_global_re[[1]] %>% 
  ddply(., .(variable, Group), summarise,mean=mean(value), sum=sum(value)) %>% 
  arrange(., mean)

for (i in c(1, 2)) {
  dat_trace <- dat_cell_trace_global_re[[i]]
  dat_trace$variable <- factor(dat_trace$variable, levels = dat_trace_sta$variable)
  p_heat <- ggplot(dat_trace, aes(Time, variable,fill= value))+ 
    geom_tile(height=2)+
    #facet_grid(rows = vars(Group), scales = "free_y")+
    #scale_fill_gradient2(limits= score_range,low = "navy", high = "red4", mid = "white")+
    scale_fill_gradientn(limits= score_range, colours = c("navy", "white", "red4"), values = rescale(c(score_range[1], 0, score_range[2])))+
    labs(x="Time relative to crossing (s)", y="Cell ID")+
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
p_heat_com_global <- plot_grid(p_heat_d1, p_heat_d2, nrow = 1)


## Global ID only plot from d5 to d7 (06192022)-----

cc_globalID_fun <- function(path_ID, path_trace, t_stim){
  Global_ID <- read.csv(path_ID, header = F) %>% 
    as_tibble() %>% 
    select(5:7) %>% 
    na_if(., 0) %>% 
    drop_na()
  
  
  ## read.xlsx with no header and remove the first col
  cc_read.xlsx<- function(x){
    dat <- read.xlsx(x, colNames = F) %>% 
      select(., -1)
    dat
  }
  
  ## do z score of the whole trace
  dat_trace <- mapply(cc_read.xlsx, path_trace,SIMPLIFY = F) %>% 
    mapply(function(x) apply(x, 2, scale), ., SIMPLIFY = F)
  
  stim_time<- seq(-2, 6.5, by=0.5)
  ## number of rows to be binned
  n <- 10 # 0.05*10=0.5
  
  ## create a 2 days list and put the activation of each cell in it
  cc_trace_pick <- function(dat_trace, cell_pick){
    dat_trace_pick <- t(dat_trace[cell_pick,])
    dat_trace_pick
  }
  
  dat_cell_trace_day <- mapply(cc_trace_pick, dat_trace, as.list(Global_ID))
  
  cc_trace_extract <- function(dat_cell_trace, t_stim_day){
    dat_stim <- vector(mode = "list", length = length(t_stim_day))
    for (i in seq_along(t_stim_day)){
      t1_p <- t_stim_day[i]
      dat_stim1 <- dat_cell_trace[(t1_p-40):(t1_p+140-1),]
      dat_stim1 <- aggregate(dat_stim1,list(rep(1:(nrow(dat_stim1)%/%n+1),each=n,len=nrow(dat_stim1))),mean)[-1]
      dat_stim[[i]] <- dat_stim1[1:5,] %>%  ## baseline as -2 to 0
        colMeans(., na.rm = T) %>% 
        sweep(dat_stim1, 2, ., FUN = "-")
    }
    ## average the trace by cell number
    if (length(t_stim_day)>1){
      dat_cell_trace_average <- aaply(laply(dat_stim, as.matrix), c(2, 3), function(x) mean(x, na.rm=T))
      
    } else {
      dat_cell_trace_average <- data.matrix(dat_stim[[1]])
      
    }
    return(dat_cell_trace_average)
  }
  dat_cell_trace_day_extract <- mapply(cc_trace_extract, dat_cell_trace_day, t_stim, SIMPLIFY = F)
  
  return(dat_cell_trace_day_extract)
}

## for m3
global_ID_m3 <- "~cchen/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/matlab_analysis/m3/Global_id/m3_anti_pin_global_ID.csv"
path_trace_m3 <- as.list(list.files(path ="~cchen/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/matlab_analysis/m3/trace_days/", pattern = "*.xlsx", full.names = T ))[c(5:7)]
t_stim_m3 <- list(t_stim_m3_d5 = c(758),t_stim_m3_d6 = c(132), t_stim_m3_d7 = c(346))
dat_global_trace_m3 <- cc_globalID_fun(global_ID_m3, path_trace_m3, t_stim_m3)

## for m7
global_ID_m7 <- "~cchen/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/matlab_analysis/m7/Global_id/m7_anti_pin_global_ID.csv"
path_trace_m7 <- as.list(list.files(path ="~cchen/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/matlab_analysis/m7/trace_days/", pattern = "*.xlsx", full.names = T ))[c(5:7)]
t_stim_m7 <- list(t_stim_m7_d5 = c(317), t_stim_m7_d6 = c(205), t_stim_m7_d7 = c(236))
dat_global_trace_m7 <- cc_globalID_fun(global_ID_m7, path_trace_m7, t_stim_m7)

## for m17
global_ID_m17 <- "~cchen/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/matlab_analysis/m17/Global_ID/m17_global_ID.csv"
path_trace_m17 <- as.list(list.files(path ="~cchen/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/matlab_analysis/m17/trace_days/", pattern = "*.xlsx", full.names = T ))[c(5:7)]
t_stim_m17 <- list(t_stim_m17_d5 = c(1114), t_stim_m17_d6 = c(617), t_stim_m17_d7 = c(157) )
dat_global_trace_m17 <- cc_globalID_fun(global_ID_m17, path_trace_m17, t_stim_m17)

## for m18
global_ID_m18 <- "~cchen/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/matlab_analysis/m18/Global_ID/m18_global_ID.csv"
path_trace_m18 <- as.list(list.files(path ="~cchen/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/matlab_analysis/m18/trace_days/", pattern = "*.xlsx", full.names = T ))[c(5:7)]
t_stim_m18 <- list(t_stim_m18_d5 = c(195), t_stim_m18_d6 = c(467),  t_stim_m18_d7 = c(124))
dat_global_trace_m18 <- cc_globalID_fun(global_ID_m18, path_trace_m18, t_stim_m18)

## for m855
global_ID_m855 <- "~cchen/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/matlab_analysis/m855/Global_ID/m855_global_ID.csv"
path_trace_m855 <- as.list(list.files(path ="~cchen/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/matlab_analysis/m855/trace_days", pattern = "*.xlsx", full.names = T ))[c(5:7)]
t_stim_m855 <- list(t_stim_m855_d5 = c(67) *2, t_stim_m855_d6 = c(64) *2, t_stim_m855_d7 = c(41)*2)
dat_global_trace_m855 <- cc_globalID_fun(global_ID_m855, path_trace_m855, t_stim_m855)

## for m857
global_ID_m857 <- "~cchen/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/matlab_analysis/m857/Global_ID/m857_global_ID.csv"
path_trace_m857 <- as.list(list.files(path ="~cchen/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/matlab_analysis/m857/trace_days", pattern = "*.xlsx", full.names = T ))[c(3:5)]
t_stim_m857 <- list(t_stim_m857_d5 = c(1900) *2, t_stim_m857_d6 = c(723) *2, t_stim_m857_d7 = c(167)*2)
dat_global_trace_m857 <- cc_globalID_fun(global_ID_m857, path_trace_m857, t_stim_m857)

## combin data
dat_global_trace_d5 <- cbind(dat_global_trace_m3[[1]], dat_global_trace_m7[[1]], dat_global_trace_m17[[1]], dat_global_trace_m18[[1]], dat_global_trace_m857[[1]], dat_global_trace_m857[[1]])
colnames(dat_global_trace_d5) <- str_c("Cell", 1:ncol(dat_global_trace_d5))

dat_global_trace_d6 <- cbind(dat_global_trace_m3[[2]], dat_global_trace_m7[[2]], dat_global_trace_m17[[2]], dat_global_trace_m18[[2]], dat_global_trace_m855[[2]], dat_global_trace_m855[[2]])
colnames(dat_global_trace_d6) <- str_c("Cell", 1:ncol(dat_global_trace_d6))

dat_global_trace_d7 <- cbind(dat_global_trace_m3[[3]], dat_global_trace_m7[[3]], dat_global_trace_m17[[3]], dat_global_trace_m18[[3]], dat_global_trace_m855[[3]], dat_global_trace_m855[[3]])
colnames(dat_global_trace_d7) <- str_c("Cell", 1:ncol(dat_global_trace_d7))

dat_cell_trace_global <- vector(mode = "list", 3)
dat_cell_trace_global[[1]] <- dat_global_trace_d5
dat_cell_trace_global[[2]] <- dat_global_trace_d6
dat_cell_trace_global[[3]] <- dat_global_trace_d7


dat_cell_trace_global_re <- vector(mode = "list", 3)

# k-menas clustering
for (i in 1:length(dat_cell_trace_global)) {
  dat_cell_trace_d <- dat_cell_trace_global[[i]]
  
  
  dat_cell_trace_d_cluster <- unname(kmeans(t(dat_cell_trace_d), centers = 3)[[1]])
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
  
  dat_cell_trace_d_re <- mutate(dat_cell_trace_d_re, Group= factor(Group, levels = c("Excited", "Neutral", "Inhibited")))
  dat_cell_trace_global_re[[i]] <- dat_cell_trace_d_re
}

score_range <- do.call(rbind, dat_cell_trace_global_re) %>% 
  .$value %>% 
  range()

## align by cells
dat_trace_sta <- dat_cell_trace_global_re[[1]] %>% 
  ddply(., .(variable, Group), summarise,mean=mean(value), sum=sum(value)) %>% 
  arrange(., mean)

for (i in c(1:3)) {
  dat_trace <- dat_cell_trace_global_re[[i]]
  dat_trace$variable <- factor(dat_trace$variable, levels = dat_trace_sta$variable)
  p_heat <- ggplot(dat_trace, aes(Time, variable,fill= value))+ 
    geom_tile(height=2)+
    #facet_grid(rows = vars(Group), scales = "free_y")+
    #scale_fill_gradient2(limits= score_range,low = "navy", high = "red4", mid = "white")+
    scale_fill_gradientn(limits= score_range, colours = c("navy", "white", "red4"), values = rescale(c(score_range[1], 0, score_range[2])))+
    labs(x="Time relative to crossing (s)", y="Cell ID")+
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
p_heat_com_global <- plot_grid(p_heat_d1, p_heat_d2, p_heat_d3, nrow = 1)

setwd("~cchen/Documents/neuroscience/Pn\ project/Figure/PDF/")
cairo_pdf("p_anti_heat_global.pdf", width = 100/25.6, height = 55/25.6, family = "Arial")
p_heat_com_global
dev.off()

library(ggvenn)
## calculate the proportion of cells active during 
cell_id_active_d5 <- dat_cell_trace_global_re[[1]] %>% 
  filter(Group != "Neutral") %>% 
  distinct(., variable)

cell_id_active_d6 <- dat_cell_trace_global_re[[2]] %>% 
  filter(Group != "Neutral") %>% 
  distinct(., variable)


cell_id_active_d7 <- dat_cell_trace_global_re[[3]] %>% 
  filter(Group != "Neutral") %>% 
  distinct(., variable)

cell_ID_overlay <- list(D5 = cell_id_active_d5$variable, D6= cell_id_active_d6$variable, D7=cell_id_active_d7$variable)

p_overlap <- ggvenn(cell_ID_overlay, fill_color = c("turquoise4", "darkslategray4", "indianred"))

setwd("~cchen/Documents/neuroscience/Pn\ project/Figure/PDF/")
cairo_pdf("p_overlap.pdf", width = 100/25.6, height = 100/25.6, family = "Arial")
p_overlap
dev.off()


n_overlay <- intersect(intersect(distinct(cell_id_active_d5, variable), distinct(cell_id_active_d6, variable)), distinct(cell_id_active_d7, variable))



dat_cell_trace_global_sta <- ddply(dat_cell_trace_global_re[[1]], .(Time, Group), summarise,n=length(value),mean=mean(value),sd=sd(value),se=sd(value)/sqrt(length(value))) %>%
  rbind(., ddply(dat_cell_trace_global_re[[2]], .(Time, Group), summarise,n=length(value),mean=mean(value),sd=sd(value),se=sd(value)/sqrt(length(value)))) %>%
  mutate(., Day = rep(c("Pre",  "Test"),each= length(stim_time)*3)) %>% 
  mutate(Day = factor(Day, levels = c("Pre",  "Test")), Group = factor(Group,levels = c("Neutral", "Excited", "Inhibited") ))

range_trace_plot <- range(dat_cell_trace_global_sta$mean)

## group by day
p_trace_global <- ggplot(dat_cell_trace_global_sta, aes(Time, mean, colour=Group))+
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

## plot the summy
dat_cell_trace_sum <- rbind(dat_cell_trace_pre_sta, dat_cell_trace_d4_sta, dat_cell_trace_d6_sta) %>%
  as_tibble() %>% 
  mutate(Day = rep(c("Pre",  "D4", "D6"), c(nrow(dat_cell_trace_pre_sta),nrow(dat_cell_trace_d4_sta), nrow(dat_cell_trace_d6_sta)))) %>%
  subset(., .$Group!="Neutral") %>%
  ddply(., .(ID,Time, Day), summarise, value=sum(value)) %>%
  ddply(., .(Time, Day), summarise,n=length(value),mean=mean(value),sd=sd(value),se=sd(value)/sqrt(length(value))) %>%
  mutate(Day = factor(Day, levels = c("Pre",  "D4", "D6")))

p_trace_sum <- ggplot(dat_cell_trace_sum, aes(Time, mean, colour=Day))+
  geom_line()+
  geom_ribbon(aes(ymin=mean-se, ymax=mean+se, fill=Day), alpha=0.2, linetype=0)+
  #scale_colour_manual(values=c("deepskyblue4", "indianred"))+
  #scale_fill_manual(values=c("deepskyblue4", "indianred"))+
  labs(x="Time relative to crossing (s)", y="z score")+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  geom_vline(xintercept = 0, col= 'red', linetype=2)+
  theme(legend.title = element_blank(), legend.position = c(0.2, 0.8))

setwd("~cchen/Documents/neuroscience/Pn\ project/Figure/PDF/")
cairo_pdf("p_trace_sum_compare.pdf", width = 90/25.6, height = 60/25.6, family = "Arial")
p_trace_sum
dev.off()

setwd("~cchen/Documents/neuroscience/Pn\ project/Figure/PDF/")
cairo_pdf("p_trace_global.pdf", width = 100/25.6, height = 65/25.6, family = "Arial")
p_trace_global
dev.off()

## calculate the spikes before and after crossing------
cc_firing_fun <- function(path_trace, path_peak, t_stim) {
  cell_valid <- read.xlsx(path_trace, colNames = F, rowNames = F) %>% 
    select (., X1)
  
  dat_peak <- read.csv(path_peak) %>% 
    select(which(cell_valid==1)) %>% 
    unlist() %>% 
    na.omit()
  
  ctrl_rang <- range(t_stim-42, t_stim-2)
  test_rang <- range(t_stim, t_stim +40)
  
  ctrl_freq <- length(dat_peak[dat_peak>= ctrl_rang[1] & dat_peak < ctrl_rang[2]])/2 # 2s before crossing
  test_freq <- length(dat_peak[dat_peak>= test_rang[1] & dat_peak < test_rang[2]])/2 # 4s after crossing
  diff_freq <- test_freq - ctrl_freq
  
  return(c(ctrl_freq, test_freq, diff_freq))
}

## for m3
path_peak_m3 <- as.list(list.files(path ="~cchen/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/matlab_analysis/m3/peaks_days", pattern = "*.csv", full.names = T ))[c(3, 7)]
dat_rate_m3 <- t(mapply(cc_firing_fun, path_trace_m3, path_peak_m3, t_stim_m3))

## for m7
path_peak_m7 <- as.list(list.files(path ="~cchen/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/matlab_analysis/m7/peaks_days/", pattern = "*.csv", full.names = T ))[c(3, 7)]
dat_rate_m7 <- t(mapply(cc_firing_fun, path_trace_m7, path_peak_m7, t_stim_m7))

## for m17
path_peak_m17 <- as.list(list.files(path ="~cchen/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/matlab_analysis/m17/peaks_days/", pattern = "*.csv", full.names = T ))[c(3, 7)]
dat_rate_m17 <- t(mapply(cc_firing_fun, path_trace_m17, path_peak_m17, t_stim_m17))

## for m18
path_peak_m18 <- as.list(list.files(path ="~cchen/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/matlab_analysis/m18/peaks_days/", pattern = "*.csv", full.names = T ))[c(3, 7)]
dat_rate_m18 <- t(mapply(cc_firing_fun, path_trace_m18, path_peak_m18, t_stim_m18))

## for m855
path_peak_m855 <- as.list(list.files(path ="~cchen/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/matlab_analysis/m855/peaks_days/", pattern = "*.csv", full.names = T ))[c(3, 7)]
dat_rate_m855 <- t(mapply(cc_firing_fun, path_trace_m855, path_peak_m855, t_stim_m855))

## for m855
path_peak_m857 <- as.list(list.files(path ="~cchen/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/matlab_analysis/m857/peaks_days/", pattern = "*.csv", full.names = T ))
dat_rate_m857 <- t(mapply(cc_firing_fun, path_trace_m857, path_peak_m857, t_stim_m857))


## combine data and plot
dat_firing <- rbind(dat_rate_m3, dat_rate_m7, dat_rate_m17, dat_rate_m18, dat_rate_m855, dat_rate_m857) %>% 
  as_tibble() %>% 
  rename(Freq_before = V1, Freq_after = V2, diff = V3) %>% 
  mutate(ID = rep(mouse_ID, each=2)) %>% 
  mutate(Group = rep(group_day, length(mouse_ID))) %>% 
  gather(variable, value, -ID,  -Group) %>% 
  ddply(., .(ID, Group, variable), summarise, value=mean(value)) %>% 
  mutate(Group= factor(Group, levels = c("Pre",  "Test")))

p_firing <- dat_firing %>% 
  filter(variable=="diff") %>% 
  ggplot(., aes(Group, value, colour=Group))+
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
  scale_y_continuous(limits = c(-4,4))+
  theme(legend.position = 'none')

t_firng <- dat_firing %>% 
  filter(variable=="diff") %>% 
  wilcox.test(value~Group,.)

setwd("~cchen/Documents/neuroscience/Pn\ project/Figure/PDF/")
cairo_pdf("p_firing.pdf", width = 40/25.6, height = 65/25.6, family = "Arial")
p_firing
dev.off() 

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
cov_range <- range(c(dat_cell_cov_list[[1]][[1]], dat_cell_cov_list[[2]][[1]]))

p_cov_d3<- ggcorrplot(dat_cell_cov_list[[1]][[1]], hc.order = TRUE, outline.color = "gray")+
  scale_fill_gradient2(limit = c(cov_range[1], cov_range[2]), low = "navy", high =  "red4", mid = "white")+
  theme_void()+
  theme(legend.position = "none")

p_cov_d7<- ggcorrplot(dat_cell_cov_list[[2]][[1]], hc.order = TRUE, outline.color = "gray")+
  scale_fill_gradient2(limit = c(cov_range[1], cov_range[2]), low = "navy", high =  "red4", mid = "white")+
  theme_void()+
  theme(legend.position = "none")
p_cov_com <- plot_grid(p_cov_d3, p_cov_d7, nrow = 1)

setwd("~cchen/Documents/neuroscience/Pn\ project/Figure/PDF/")
cairo_pdf("p_cov_com.pdf", width = 170/25.6, height = 60/25.6, family = "Arial")
p_cov_com
dev.off()

## compare the covrelation between cells
dat_cell_cov <- c()
dat_cell_cov_sta <- c()
for (i in c(1,2)){
  c_trim <- function(x){
    x[upper.tri(x)]<- NA
    c(x[!is.na(x)])
  }
  cov_value <- abs(unlist(mapply(c_trim, dat_cell_cov_list[[i]] )))
  dat_cov <- data.frame(Day=group_day[i], value=cov_value)
  dat_cell_cov <- rbind(dat_cell_cov, dat_cov)
  
  ## for mean of cov
  cov_value_mean <- dat_cell_cov_list[[i]] %>%
    mapply(function(x) ifelse(is.null(ncol(x)), NA,mean(colMeans(abs(x), na.rm=T))), .)
  
  cov_value_max <- dat_cell_cov_list[[i]] %>%
    mapply(function(x) ifelse(is.null(ncol(x)),NA,mean(apply(x, 2, function(x) max(abs(x), na.rm = T)))), .)
  
  mouse_ID <- c(1: length(dat_cell_cov_list[[i]]))
  dat_cov_mean <- data.frame(Day = group_day[i], ID= mouse_ID, value=cov_value_mean, value_max=cov_value_max)
  dat_cell_cov_sta <- rbind(dat_cell_cov_sta, dat_cov_mean)
}


## cumulative plot of matrix
dat_cell_cov$Day <- factor(dat_cell_cov$Day, levels = c("Pre",  "Test"))
p_cov_cum <- ggplot(dat_cell_cov, aes(value, group=Day, colour=Day))+
  stat_ecdf(geom = "step")+
  scale_color_manual(values=c( "deepskyblue4", "indianred"))+
  labs(x="Covariance (z-score)", y="Cummulative probability")+
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

setwd("~cchen/Documents/neuroscience/Pn\ project/Figure/PDF/")
cairo_pdf("p_cov_cum.pdf", width = 70/25.6, height = 56/25.6, family = "Arial")
p_cov_cum
dev.off()

## bar plot

dat_cell_cov_sta1 <- melt(dat_cell_cov_sta, id.vars = c("ID", "Day")) %>%
  ddply(., .(ID, Day, variable), summarise, value=mean(value, na.rm = T)) %>%
  mutate(Day = factor(Day, levels = c("Pre",  "Test")))

dat_cell_cov_sta2 <- dat_cell_cov_sta1 %>% 
  ddply(., .(Day, variable), summarise, n=length(value),mean=mean(value),sd=sd(value),se=sd(value)/sqrt(length(value)))

p_cov_mean <- subset(dat_cell_cov_sta1, dat_cell_cov_sta1$variable=="value") %>%
  ggplot(., aes(Day, value, colour=Day))+
  geom_boxplot(outlier.shape = NA)+
  geom_line(aes(group=ID), colour="gray90")+
  geom_jitter(width = 0.2, shape=1)+
  scale_colour_manual(values=c("deepskyblue4", "indianred"))+
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
  ggplot(., aes(Day, value, colour=Day))+
  geom_boxplot(outlier.shape = NA)+
  geom_line(aes(group=ID), colour="gray90")+
  geom_jitter(position=position_jitterdodge(), shape=1, size=1)+
  scale_colour_manual(values=c("deepskyblue4", "indianred"))+
  labs(x="", y="Max covariance (z-score)")+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  scale_y_continuous(limits = c(0.3, 1))+
  theme(legend.title = element_blank(), legend.position = "none")

p_cov_bar <- plot_grid(p_cov_mean, p_cov_best, nrow = 1)

t_cov_mean <- dat_cell_cov_sta1 %>% 
  filter(variable == "value") %>% 
  t.test(value~Day,., paired = T)

t_vov_max <- dat_cell_cov_sta1 %>% 
  filter(variable == "value_max") %>% 
  t.test(value~Day,., paired = T)

setwd("~cchen/Documents/neuroscience/Pn\ project/Figure/PDF/")
cairo_pdf("p_cov_best.pdf", width = 35/25.6, height = 60/25.6, family = "Arial")
p_cov_best
dev.off()

setwd("~cchen/Documents/neuroscience/Pn\ project/Figure/PDF/")
cairo_pdf("p_cov_mean.pdf", width = 35/25.6, height = 60/25.6, family = "Arial")
p_cov_mean
dev.off()

## analyze the data from IT neurons to prove the specifivity of PT neurons during PAC-----
c_miniscope_inscopix <- function(ID_trace, t_stim) {
  ## import and format the data
  dat_trace <- read.csv(ID_trace, skip = 1) %>% 
    as_tibble() %>% 
    select(!starts_with("r")) %>% 
    rename(Time = names(.)[1])%>%
    mutate_at(vars(-("Time")),scale)
  
  ## check the time of recording frequency
  inter_time <- dat_trace$Time[2] - dat_trace$Time[1]
  
  if (inter_time < 1){
    dat_trace <- dat_trace %>% 
      filter(row_number() %% 2 == 1)
    
  } else{
    dat_trace <- dat_trace
  }
  
  ## analyze the trace when mouse crossing border
  time_crossing <- which.min(abs(dat_trace$Time - t_stim))
  
  dat_trace_trunct <- dat_trace[(time_crossing-10):(time_crossing+ 32),] %>% 
    as_tibble() %>% 
    mutate_at(vars(-("Time")), function(x) x- mean(x[1:10])) %>% 
    select(-Time)
  # mutate(mean = rowMeans(select(.,starts_with("a")))) %>% 
  # mutate(Time = seq(-2, 6.5, 0.2))
  
  #ggplot(dat_trace_trunct, aes(Time, mean)) + geom_line()
  
  return(dat_trace_trunct)
  
}


## Extract trace from each mice
# for m18
path_trace_m18 <- as.list(list.files(path ="~cchen/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/nvue/IT_PT/m18", pattern = ".*Green.*\\.csv", full.names = T ))[c(3,7)]
t_stim_m18 <- list(t_stim_m18_d3 = c(1843/10), t_stim_m18_d7 = c(423/10))
dat_trace_m18 <- mapply(c_miniscope_inscopix, path_trace_m18, t_stim_m18, SIMPLIFY = F)

## for m19
path_trace_m19 <- as.list(list.files(path ="~cchen/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/nvue/IT_PT/m19", pattern = ".*Green.*\\.csv", full.names = T ))[c(3,7)]
t_stim_m19 <- list(t_stim_m19_d3 = c(210/10), t_stim_m19_d7 = c(212/10))
dat_trace_m19 <- mapply(c_miniscope_inscopix, path_trace_m19, t_stim_m19, SIMPLIFY = F)

## for m20
path_trace_m20 <- as.list(list.files(path ="~cchen/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/nvue/IT_PT/m20", pattern = ".*Green.*\\.csv", full.names = T ))[c(3,7)]
t_stim_m20 <- list(t_stim_m20_d3 = c(218/10), t_stim_m20_d7 = c(249/10))
dat_trace_m20 <- mapply(c_miniscope_inscopix, path_trace_m20, t_stim_m20, SIMPLIFY = F)

## combine data and do k-means analysis

dat_cell_trace <- vector(mode = "list", 2)
for (i in 1:2) {
  dat_cell_trace[[i]] <- cbind(dat_trace_m18[[i]], dat_trace_m19[[i]], dat_trace_m20[[i]])
  
}


dat_cell_trace_re <- vector(mode = "list", 2)

mouse_ID <- c("m18", "m19", "m20")
for (i in 1:length(dat_cell_trace)) {
  dat_cell_trace_d <- dat_cell_trace[[i]]
  cell_id <- sprintf("Cell%s",seq(1:ncol(dat_cell_trace_d)))
  colnames(dat_cell_trace_d)<- cell_id
  
  dat_cell_trace_d_cluster <- unname(kmeans(t(dat_cell_trace_d), centers = 3, nstart = 50)[[1]])
  stim_time<- seq(-2, 6.5, by=0.2)
  
  dat_cell_trace_d_re<- dat_cell_trace_d %>% 
    mutate(., Time= stim_time) %>% 
    pivot_longer(-(Time), names_to = "Cell") %>% 
    slice(mixedorder(Cell)) %>% 
    mutate(., Group= rep(dat_cell_trace_d_cluster,each = length(stim_time)) )
  
  ## sort data by the value in each group
  dat_cell_d_sort <- as.numeric(names(sort(tapply(dat_cell_trace_d_re$value, dat_cell_trace_d_re$Group, mean), decreasing = T)))
  
  dat_cell_trace_d_re$Group[dat_cell_trace_d_re$Group == dat_cell_d_sort[1]] ="Excited"
  dat_cell_trace_d_re$Group[dat_cell_trace_d_re$Group == dat_cell_d_sort[2]] ="Neutral"
  dat_cell_trace_d_re$Group[dat_cell_trace_d_re$Group == dat_cell_d_sort[3]] ="Inhibited"
  
  rep_time <- c(ncol(dat_trace_m18[[i]]), ncol(dat_trace_m19[[i]]), ncol(dat_trace_m20[[i]]))
  dat_cell_trace_d_re <- mutate(dat_cell_trace_d_re, Group= factor(Group, levels = c("Excited", "Neutral", "Inhibited"))) %>% 
    mutate(., ID = rep(mouse_ID, rep_time*length(stim_time)) )
  
  dat_cell_trace_re[[i]] <- dat_cell_trace_d_re
}

## set the range of z score during all days
score_range <- rbind(dat_cell_trace_re[[1]], dat_cell_trace_re[[2]]) %>% 
  .$value %>% 
  range()
group_day <- c("Pre", "Test")

for (i in c(1: 2)) {
  dat_trace <- dat_cell_trace_re[[i]]
  dat_trace_sta <- ddply(dat_trace, .(Cell, Group), summarise,mean=mean(value), sum=sum(value))
  dat_trace_sta <- dat_trace_sta[order(dat_trace_sta[,'mean']),]
  dat_trace$Cell <- factor(dat_trace$Cell, levels = dat_trace_sta$Cell)
  p_heat <- ggplot(dat_trace, aes(Time, Cell,fill= value))+ 
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
p_heat_com <- plot_grid(p_heat_Pre, p_heat_Test, nrow = 1)

setwd("~cchen/Documents/neuroscience/Pn\ project/Figure/PDF/")
cairo_pdf("p_anti_heat_IT.pdf", width = 90/25.6, height = 55/25.6, family = "Arial")
p_heat_com
dev.off()

dat_cell_trace_pre_sta <-  ddply(dat_cell_trace_re[[1]], .(ID,Time, Group), summarise,value=mean(value)) 
dat_cell_trace_test_sta <- ddply(dat_cell_trace_re[[2]], .(ID,Time, Group), summarise,value=mean(value))

dat_cell_trace_sta <- ddply(dat_cell_trace_pre_sta, .(Time, Group), summarise,n=length(value),mean=mean(value),sd=sd(value),se=sd(value)/sqrt(length(value))) %>%
  rbind(., ddply(dat_cell_trace_test_sta, .(Time, Group), summarise,n=length(value),mean=mean(value),sd=sd(value),se=sd(value)/sqrt(length(value)))) %>%
  mutate(., Day = rep(c("Pre",  "Test"),each= length(stim_time)*3)) %>% 
  mutate(Day = factor(Day, levels = c("Pre",  "Test")), Group = factor(Group,levels = c("Neutral", "Excited", "Inhibited") ))

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

dat_cell_trace_sum <- rbind(dat_cell_trace_pre_sta, dat_cell_trace_test_sta) %>%
  as_tibble() %>% 
  mutate(Day = rep(c("Pre",  "Test"), c(nrow(dat_cell_trace_pre_sta),nrow(dat_cell_trace_test_sta)))) %>%
  filter(Group!="Neutral") %>%
  ddply(., .(ID,Time, Day), summarise, value=sum(value)) %>%
  ddply(., .(Time, Day), summarise,n=length(value),mean=mean(value),sd=sd(value),se=sd(value)/sqrt(length(value))) %>%
  mutate(Day = factor(Day, levels = c("Pre",  "Test")))

# p_trace <- dat_cell_trace_sum %>% 
#   filter(ID =="m20") %>% 
#   ggplot(., aes(Time, value, colour = Day))+
#   geom_line()

p_trace_sum <- ggplot(dat_cell_trace_sum, aes(Time, mean, colour=Day))+
  geom_line()+
  geom_ribbon(aes(ymin=mean-se, ymax=mean+se, fill=Day), alpha=0.2, linetype=0)+
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
  geom_vline(xintercept = 0, col= 'red', linetype=2)+
  theme(legend.title = element_blank(), legend.position = c(0.2, 0.8))

setwd("~cchen/Documents/neuroscience/Pn\ project/Figure/PDF/")
cairo_pdf("p_trace_sum_it.pdf", width = 80/25.6, height = 60/25.6, family = "Arial")
p_trace_sum
dev.off()

dat_cell_area <- c()
group_day <- c("Pre",  "Test")
for (i in c(1,2)){
  rep_time <- c(ncol(dat_trace_m18[[i]]), ncol(dat_trace_m19[[i]]), ncol(dat_trace_m20[[i]]))
  dat_trace <- dat_cell_trace_re[[i]]
  dat_trace$ID <- rep(mouse_ID, rep_time*length(stim_time))
  dat_trace_sta <- ddply(dat_trace, .(ID, Time, Group), summarise,n=length(value),mean=mean(value),sd=sd(value),se=sd(value)/sqrt(length(value)))
  Time_area <- which(stim_time>=0)
  dat_anti_area <- as.data.frame(tapply(dat_trace_sta$mean, INDEX = list(dat_trace_sta$ID, dat_trace_sta$Group), 
                                        function (x) AUC(stim_time[Time_area], x[Time_area])))
  dat_anti_area$ID <- rownames(dat_anti_area)
  rownames(dat_anti_area)<-NULL
  dat_anti_area$Day <- group_day[i]
  dat_anti_area$ratio <- abs(dat_anti_area$Excited/dat_anti_area$Inhibited)
  dat_anti_area$sum <- dat_anti_area$Excited + dat_anti_area$Inhibited
  dat_cell_area <- rbind(dat_cell_area, dat_anti_area)
}

dat_cell_area$Day <- factor(dat_cell_area$Day, levels = c("Pre","Test"))

dat_EI_sum_sta <- dat_cell_area %>% 
  select(Day, ID, sum) %>% 
  ddply(., .(Day), summarise,n=length(sum),mean=mean(sum),sd=sd(sum),se=sd(sum)/sqrt(length(sum)))

p_sum <- dat_cell_area %>% 
  select(Day, ID, sum) %>% 
  wilcox.test(sum~Day, .)

p_EI_sum<- dat_cell_area %>% 
  select(Day, ID, sum) %>% 
  ggplot(., aes(Day, sum, colour=Day))+
  geom_boxplot(outlier.shape = NA)+
  geom_line(aes(group=ID), colour="gray90")+
  geom_jitter(width = 0.2, shape=1)+
  scale_colour_manual(values=c("deepskyblue4", "indianred" ))+
  labs(x="", y="AUC (z-score)")+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  theme(legend.title = element_blank(), legend.position = "none")

setwd("~cchen/Documents/neuroscience/Pn\ project/Figure/PDF/")
cairo_pdf("p_trace_sum_com_it.pdf", width = 40/25.6, height = 60/25.6, family = "Arial")
p_EI_sum
dev.off()


## cell covariance analysis
cc_cov_fun <- function(cell_trace){
  dat_cell_trace<- subset(cell_trace, cell_trace$Group!="Neutral")
  ID_cell <- unique(dat_cell_trace$ID)
  dat_cell_cov <- vector(mode = "list", length(ID_cell))
  
  for (i in seq_along(ID_cell)){
    dat_cell_trace_cov <- subset(dat_cell_trace, dat_cell_trace$ID==ID_cell[i]) %>%
      select("Time", "Cell","value") %>%
      dcast(., Time~Cell)
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
cov_range <- range(c(dat_cell_cov_list[[1]][[1]], dat_cell_cov_list[[2]][[1]]))

p_cov_d3<- ggcorrplot(dat_cell_cov_list[[1]][[1]], hc.order = TRUE, outline.color = "gray")+
  scale_fill_gradient2(limit = c(cov_range[1], cov_range[2]), low = "navy", high =  "red4", mid = "white")+
  theme_void()+
  theme(legend.position = "none")

p_cov_d7<- ggcorrplot(dat_cell_cov_list[[2]][[1]], hc.order = TRUE, outline.color = "gray")+
  scale_fill_gradient2(limit = c(cov_range[1], cov_range[2]), low = "navy", high =  "red4", mid = "white")+
  theme_void()+
  theme(legend.position = "none")
p_cov_com <- plot_grid(p_cov_d3, p_cov_d7, nrow = 1)


dat_cell_cov <- c()
dat_cell_cov_sta <- c()
for (i in c(1,2)){
  c_trim <- function(x){
    x[upper.tri(x)]<- NA
    c(x[!is.na(x)])
  }
  cov_value <- abs(unlist(mapply(c_trim, dat_cell_cov_list[[i]] )))
  dat_cov <- data.frame(Day=group_day[i], value=cov_value)
  dat_cell_cov <- rbind(dat_cell_cov, dat_cov)
  
  ## for mean of cov
  cov_value_mean <- dat_cell_cov_list[[i]] %>%
    mapply(function(x) ifelse(is.null(ncol(x)), NA,mean(colMeans(abs(x), na.rm=T))), .)
  
  cov_value_max <- dat_cell_cov_list[[i]] %>%
    mapply(function(x) ifelse(is.null(ncol(x)),NA,mean(apply(x, 2, function(x) max(abs(x), na.rm = T)))), .)
  
  mouse_ID <- c(1: length(dat_cell_cov_list[[i]]))
  dat_cov_mean <- data.frame(Day = group_day[i], ID= mouse_ID, value=cov_value_mean, value_max=cov_value_max)
  dat_cell_cov_sta <- rbind(dat_cell_cov_sta, dat_cov_mean)
}


## cumulative plot of matrix
dat_cell_cov$Day <- factor(dat_cell_cov$Day, levels = c("Pre",  "Test"))
p_cov_cum <- ggplot(dat_cell_cov, aes(value, group=Day, colour=Day))+
  stat_ecdf(geom = "step")+
  scale_color_manual(values=c( "deepskyblue4", "indianred"))+
  labs(x="Covariance (z-score)", y="Cummulative probability")+
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

dat_cell_cov_sta1 <- melt(dat_cell_cov_sta, id.vars = c("ID", "Day")) %>%
  ddply(., .(ID, Day, variable), summarise, value=mean(value, na.rm = T)) %>%
  mutate(Day = factor(Day, levels = c("Pre",  "Test")))

dat_cell_cov_sta2 <- dat_cell_cov_sta1 %>% 
  ddply(., .(Day, variable), summarise, n=length(value),mean=mean(value),sd=sd(value),se=sd(value)/sqrt(length(value)))

p_cov_mean <- subset(dat_cell_cov_sta1, dat_cell_cov_sta1$variable=="value") %>%
  ggplot(., aes(Day, value, colour=Day))+
  geom_boxplot(outlier.shape = NA)+
  geom_line(aes(group=ID), colour="gray90")+
  geom_jitter(width = 0.2, shape=1)+
  scale_colour_manual(values=c("deepskyblue4", "indianred"))+
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
  ggplot(., aes(Day, value, colour=Day))+
  geom_boxplot(outlier.shape = NA)+
  geom_line(aes(group=ID), colour="gray90")+
  geom_jitter(position=position_jitterdodge(), shape=1, size=1)+
  scale_colour_manual(values=c("deepskyblue4", "indianred"))+
  labs(x="", y="Max covariance (z-score)")+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  scale_y_continuous(limits = c(0.3, 1))+
  theme(legend.title = element_blank(), legend.position = "none")

p_cov_bar <- plot_grid(p_cov_mean, p_cov_best, nrow = 1)

## correaltion analysis single trace between mean activity-----
cc_cov_fun <- function(cell_trace){
  dat_cell_trace<- subset(cell_trace, cell_trace$Group!="Neutral")
  # dat_cell_trace<- cell_trace
  ID_cell <- unique(dat_cell_trace$ID)
  dat_cell_cov <- vector(mode = "list", length(ID_cell))
  
  for (i in seq_along(ID_cell)){
    dat_cell_trace_cov <- subset(dat_cell_trace, dat_cell_trace$ID==ID_cell[i]) %>%
      select("Time", "variable","value") %>%
      dcast(., Time~variable) %>% 
      select(-Time)
    res_cov <- rep(0, ncol(dat_cell_trace_cov))
    if (ncol(dat_cell_trace_cov)<2){
      res_cov <- NA
    } else {
      for (j in 1:ncol(dat_cell_trace_cov)){
        x <- dat_cell_trace_cov[,j]
        y <- rowMeans(dat_cell_trace_cov[,-j])
        res_cov[j] <- cor(x, y, method = "spearman")
        
      }
      
    }
    dat_cell_cov[[i]] <- res_cov
  }
  return(dat_cell_cov)
}


dat_cell_cov_list <- mapply(cc_cov_fun, dat_cell_trace_re, SIMPLIFY = F)

cor_sync <- unlist( dat_cell_cov_list[[1]])

cor_after <- unlist(dat_cell_cov_list[[2]])

cor_value <- c(cor_sync, cor_after)
cor_group = c(rep("Pre", length(cor_sync)), rep("Post", length(cor_after)))
dat_cor <- tibble(Group = cor_group, cor = cor_value) %>% 
  mutate(Group = factor(Group, levels = c("Pre", "Post"))) %>% 
  ggplot(., aes(Group, cor, group = Group))+
  geom_violin()+
  #geom_boxplot(outlier.shape = NA)+
  geom_jitter(aes(colour= Group,shape = Group),width = 0.2, size = 2)+
  scale_color_manual(values=c("darkcyan", "indianred"))+
  scale_shape_manual(values=c(20, 18))+
  labs(x="", y="Neuronal synchrony")+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  theme(legend.title = element_blank(), legend.position = "none")+
  geom_hline(yintercept = 0, linetype="dashed", color= "gray")

# t_cor <- wilcox.test(cor~Group, dat_cor)

setwd("~cchen/Documents/neuroscience/Pn\ project/Figure/PDF/")
cairo_pdf("dat_cor.pdf", width = 45/25.6, height = 60/25.6, family = "Arial")
dat_cor
dev.off()
## analyze the data from IT neurons to prove the specifivity of IT neurons during PAC-----
c_miniscope_inscopix <- function(ID_trace, t_stim) {
  ## import and format the data
  dat_trace <- read.csv(ID_trace, skip = 1) %>% 
    as_tibble() %>% 
    select(!starts_with("r")) %>% 
    rename(Time = names(.)[1])%>%
    mutate_at(vars(-("Time")),scale)
  
  ## check the time of recording frequency
  inter_time <- dat_trace$Time[2] - dat_trace$Time[1]
  
  if (inter_time < 1){
    dat_trace <- dat_trace %>% 
      filter(row_number() %% 2 == 1)
    
  } else{
    dat_trace <- dat_trace
  }
  
  ## analyze the trace when mouse crossing border
  time_crossing <- which.min(abs(dat_trace$Time - t_stim))
  
  dat_trace_trunct <- dat_trace[(time_crossing-10):(time_crossing+ 32),] %>% 
    as_tibble() %>% 
    mutate_at(vars(-("Time")), function(x) x- mean(x[1:10])) %>% 
    select(-Time)
  # mutate(mean = rowMeans(select(.,starts_with("a")))) %>% 
  # mutate(Time = seq(-2, 6.5, 0.2))
  
  #ggplot(dat_trace_trunct, aes(Time, mean)) + geom_line()
  
  return(dat_trace_trunct)
  
}


## Extract trace from each mice
# for m18
path_trace_m18 <- as.list(list.files(path ="~cchen/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/nvue/IT_PT/m18", pattern = ".*Green.*\\.csv", full.names = T ))[c(3,7)]
t_stim_m18 <- list(t_stim_m18_d3 = c(1843/10), t_stim_m18_d7 = c(423/10))
dat_trace_m18 <- mapply(c_miniscope_inscopix, path_trace_m18, t_stim_m18, SIMPLIFY = F)

## for m19
path_trace_m19 <- as.list(list.files(path ="~cchen/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/nvue/IT_PT/m19", pattern = ".*Green.*\\.csv", full.names = T ))[c(3,7)]
t_stim_m19 <- list(t_stim_m19_d3 = c(210/10), t_stim_m19_d7 = c(212/10))
dat_trace_m19 <- mapply(c_miniscope_inscopix, path_trace_m19, t_stim_m19, SIMPLIFY = F)

## for m20
path_trace_m20 <- as.list(list.files(path ="~cchen/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/nvue/IT_PT/m20", pattern = ".*Green.*\\.csv", full.names = T ))[c(3,7)]
t_stim_m20 <- list(t_stim_m20_d3 = c(218/10), t_stim_m20_d7 = c(249/10))
dat_trace_m20 <- mapply(c_miniscope_inscopix, path_trace_m20, t_stim_m20, SIMPLIFY = F)

## combine data and do k-means analysis

dat_cell_trace <- vector(mode = "list", 2)
for (i in 1:2) {
  dat_cell_trace[[i]] <- cbind(dat_trace_m18[[i]], dat_trace_m19[[i]], dat_trace_m20[[i]])
  
}


dat_cell_trace_re <- vector(mode = "list", 2)

mouse_ID <- c("m18", "m19", "m20")
for (i in 1:length(dat_cell_trace)) {
  dat_cell_trace_d <- dat_cell_trace[[i]]
  cell_id <- sprintf("Cell%s",seq(1:ncol(dat_cell_trace_d)))
  colnames(dat_cell_trace_d)<- cell_id
  
  dat_cell_trace_d_cluster <- unname(kmeans(t(dat_cell_trace_d), centers = 3, nstart = 50)[[1]])
  stim_time<- seq(-2, 6.5, by=0.2)
  
  dat_cell_trace_d_re<- dat_cell_trace_d %>% 
    mutate(., Time= stim_time) %>% 
    pivot_longer(-(Time), names_to = "Cell") %>% 
    slice(mixedorder(Cell)) %>% 
    mutate(., Group= rep(dat_cell_trace_d_cluster,each = length(stim_time)) )
  
  ## sort data by the value in each group
  dat_cell_d_sort <- as.numeric(names(sort(tapply(dat_cell_trace_d_re$value, dat_cell_trace_d_re$Group, mean), decreasing = T)))
  
  dat_cell_trace_d_re$Group[dat_cell_trace_d_re$Group == dat_cell_d_sort[1]] ="Excited"
  dat_cell_trace_d_re$Group[dat_cell_trace_d_re$Group == dat_cell_d_sort[2]] ="Neutral"
  dat_cell_trace_d_re$Group[dat_cell_trace_d_re$Group == dat_cell_d_sort[3]] ="Inhibited"
  
  rep_time <- c(ncol(dat_trace_m18[[i]]), ncol(dat_trace_m19[[i]]), ncol(dat_trace_m20[[i]]))
  dat_cell_trace_d_re <- mutate(dat_cell_trace_d_re, Group= factor(Group, levels = c("Excited", "Neutral", "Inhibited"))) %>% 
    mutate(., ID = rep(mouse_ID, rep_time*length(stim_time)) )
  
  dat_cell_trace_re[[i]] <- dat_cell_trace_d_re
}

## set the range of z score during all days

heat_m_ID <- "m19"
score_rang1 <- dat_cell_trace_re[[1]] %>% 
  filter(ID == heat_m_ID)

score_range <- dat_cell_trace_re[[2]] %>% 
  filter(ID ==heat_m_ID) %>% 
  rbind(., score_rang1) %>% 
  .$value %>% 
  range()
group_day <- c("Pre", "Test")

for (i in c(1: 2)) {
  dat_trace <- dat_cell_trace_re[[i]] %>% 
    filter(ID == heat_m_ID)
  dat_trace_sta <- ddply(dat_trace, .(Cell, Group), summarise,mean=mean(value), sum=sum(value))
  dat_trace_sta <- dat_trace_sta[order(dat_trace_sta[,'mean']),]
  dat_trace$Cell <- factor(dat_trace$Cell, levels = dat_trace_sta$Cell)
  p_heat <- ggplot(dat_trace, aes(Time, Cell,fill= value))+ 
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
p_heat_com <- plot_grid(p_heat_Pre, p_heat_Test, nrow = 1)

setwd("~cchen/Documents/neuroscience/Pn\ project/Figure/PDF/")
cairo_pdf("p_anti_heat_IT.pdf", width = 90/25.6, height = 55/25.6, family = "Arial")
p_heat_com
dev.off()

dat_cell_trace_pre_sta <-  ddply(dat_cell_trace_re[[1]], .(ID,Time, Group), summarise,value=mean(value)) 
dat_cell_trace_test_sta <- ddply(dat_cell_trace_re[[2]], .(ID,Time, Group), summarise,value=mean(value))

dat_cell_trace_sta <- ddply(dat_cell_trace_pre_sta, .(Time, Group), summarise,n=length(value),mean=mean(value),sd=sd(value),se=sd(value)/sqrt(length(value))) %>%
  rbind(., ddply(dat_cell_trace_test_sta, .(Time, Group), summarise,n=length(value),mean=mean(value),sd=sd(value),se=sd(value)/sqrt(length(value)))) %>%
  mutate(., Day = rep(c("Pre",  "Test"),each= length(stim_time)*3)) %>% 
  mutate(Day = factor(Day, levels = c("Pre",  "Test")), Group = factor(Group,levels = c("Neutral", "Excited", "Inhibited") ))

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

dat_cell_trace_sum <- rbind(dat_cell_trace_pre_sta, dat_cell_trace_test_sta) %>%
  as_tibble() %>% 
  mutate(Day = rep(c("Pre",  "Test"), c(nrow(dat_cell_trace_pre_sta),nrow(dat_cell_trace_test_sta)))) %>%
  filter(Group!="Neutral") %>%
  ddply(., .(ID,Time, Day), summarise, value=sum(value)) %>%
  ddply(., .(Time, Day), summarise,n=length(value),mean=mean(value),sd=sd(value),se=sd(value)/sqrt(length(value))) %>%
  mutate(Day = factor(Day, levels = c("Pre",  "Test")))

# p_trace <- dat_cell_trace_sum %>% 
#   filter(ID =="m20") %>% 
#   ggplot(., aes(Time, value, colour = Day))+
#   geom_line()

p_trace_sum <- ggplot(dat_cell_trace_sum, aes(Time, mean, colour=Day))+
  geom_line()+
  geom_ribbon(aes(ymin=mean-se, ymax=mean+se, fill=Day), alpha=0.2, linetype=0)+
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
  geom_vline(xintercept = 0, col= 'red', linetype=2)+
  theme(legend.title = element_blank(), legend.position = c(0.2, 0.8))

setwd("~cchen/Documents/neuroscience/Pn\ project/Figure/PDF/")
cairo_pdf("p_trace_sum_it.pdf", width = 80/25.6, height = 60/25.6, family = "Arial")
p_trace_sum
dev.off()

dat_cell_area <- c()
group_day <- c("Pre",  "Test")
for (i in c(1,2)){
  rep_time <- c(ncol(dat_trace_m18[[i]]), ncol(dat_trace_m19[[i]]))
  dat_trace <- dat_cell_trace_re[[i]]
  dat_trace$ID <- rep(mouse_ID, rep_time*length(stim_time))
  dat_trace_sta <- ddply(dat_trace, .(ID, Time, Group), summarise,n=length(value),mean=mean(value),sd=sd(value),se=sd(value)/sqrt(length(value)))
  Time_area <- which(stim_time>=0)
  dat_anti_area <- as.data.frame(tapply(dat_trace_sta$mean, INDEX = list(dat_trace_sta$ID, dat_trace_sta$Group), 
                                        function (x) AUC(stim_time[Time_area], x[Time_area])))
  dat_anti_area$ID <- rownames(dat_anti_area)
  rownames(dat_anti_area)<-NULL
  dat_anti_area$Day <- group_day[i]
  dat_anti_area$ratio <- abs(dat_anti_area$Excited/dat_anti_area$Inhibited)
  dat_anti_area$sum <- dat_anti_area$Excited + dat_anti_area$Inhibited
  dat_cell_area <- rbind(dat_cell_area, dat_anti_area)
}

dat_cell_area$Day <- factor(dat_cell_area$Day, levels = c("Pre","Test"))

dat_EI_sum_sta <- dat_cell_area %>% 
  select(Day, ID, sum) %>% 
  ddply(., .(Day), summarise,n=length(sum),mean=mean(sum),sd=sd(sum),se=sd(sum)/sqrt(length(sum)))

p_sum <- dat_cell_area %>% 
  select(Day, ID, sum) %>% 
  wilcox.test(sum~Day, .)

p_EI_sum<- dat_cell_area %>% 
  select(Day, ID, sum) %>% 
  ggplot(., aes(Day, sum, colour=Day))+
  geom_boxplot(outlier.shape = NA)+
  geom_line(aes(group=ID), colour="gray90")+
  geom_jitter(width = 0.2, shape=1)+
  scale_colour_manual(values=c("deepskyblue4", "indianred" ))+
  labs(x="", y="AUC (z-score)")+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  theme(legend.title = element_blank(), legend.position = "none")

setwd("~cchen/Documents/neuroscience/Pn\ project/Figure/PDF/")
cairo_pdf("p_trace_sum_com_pt.pdf", width = 40/25.6, height = 60/25.6, family = "Arial")
p_EI_sum
dev.off()


## cell covariance analysis
cc_cov_fun <- function(cell_trace){
  dat_cell_trace<- subset(cell_trace, cell_trace$Group!="Neutral")
  ID_cell <- unique(dat_cell_trace$ID)
  dat_cell_cov <- vector(mode = "list", length(ID_cell))
  
  for (i in seq_along(ID_cell)){
    dat_cell_trace_cov <- subset(dat_cell_trace, dat_cell_trace$ID==ID_cell[i]) %>%
      select("Time", "Cell","value") %>%
      dcast(., Time~Cell)
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
cov_range <- range(c(dat_cell_cov_list[[1]][[1]], dat_cell_cov_list[[2]][[1]]))

p_cov_d3<- ggcorrplot(dat_cell_cov_list[[1]][[1]], hc.order = TRUE, outline.color = "gray")+
  scale_fill_gradient2(limit = c(cov_range[1], cov_range[2]), low = "navy", high =  "red4", mid = "white")+
  theme_void()+
  theme(legend.position = "none")

p_cov_d7<- ggcorrplot(dat_cell_cov_list[[2]][[1]], hc.order = TRUE, outline.color = "gray")+
  scale_fill_gradient2(limit = c(cov_range[1], cov_range[2]), low = "navy", high =  "red4", mid = "white")+
  theme_void()+
  theme(legend.position = "none")
p_cov_com <- plot_grid(p_cov_d3, p_cov_d7, nrow = 1)


dat_cell_cov <- c()
dat_cell_cov_sta <- c()
for (i in c(1,2)){
  c_trim <- function(x){
    x[upper.tri(x)]<- NA
    c(x[!is.na(x)])
  }
  cov_value <- abs(unlist(mapply(c_trim, dat_cell_cov_list[[i]] )))
  dat_cov <- data.frame(Day=group_day[i], value=cov_value)
  dat_cell_cov <- rbind(dat_cell_cov, dat_cov)
  
  ## for mean of cov
  cov_value_mean <- dat_cell_cov_list[[i]] %>%
    mapply(function(x) ifelse(is.null(ncol(x)), NA,mean(colMeans(abs(x), na.rm=T))), .)
  
  cov_value_max <- dat_cell_cov_list[[i]] %>%
    mapply(function(x) ifelse(is.null(ncol(x)),NA,mean(apply(x, 2, function(x) max(abs(x), na.rm = T)))), .)
  
  mouse_ID <- c(1: length(dat_cell_cov_list[[i]]))
  dat_cov_mean <- data.frame(Day = group_day[i], ID= mouse_ID, value=cov_value_mean, value_max=cov_value_max)
  dat_cell_cov_sta <- rbind(dat_cell_cov_sta, dat_cov_mean)
}


## cumulative plot of matrix
dat_cell_cov$Day <- factor(dat_cell_cov$Day, levels = c("Pre",  "Test"))
p_cov_cum <- ggplot(dat_cell_cov, aes(value, group=Day, colour=Day))+
  stat_ecdf(geom = "step")+
  scale_color_manual(values=c( "deepskyblue4", "indianred"))+
  labs(x="Covariance (z-score)", y="Cummulative probability")+
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

dat_cell_cov_sta1 <- melt(dat_cell_cov_sta, id.vars = c("ID", "Day")) %>%
  ddply(., .(ID, Day, variable), summarise, value=mean(value, na.rm = T)) %>%
  mutate(Day = factor(Day, levels = c("Pre",  "Test")))

dat_cell_cov_sta2 <- dat_cell_cov_sta1 %>% 
  ddply(., .(Day, variable), summarise, n=length(value),mean=mean(value),sd=sd(value),se=sd(value)/sqrt(length(value)))

p_cov_mean <- subset(dat_cell_cov_sta1, dat_cell_cov_sta1$variable=="value") %>%
  ggplot(., aes(Day, value, colour=Day))+
  geom_boxplot(outlier.shape = NA)+
  geom_line(aes(group=ID), colour="gray90")+
  geom_jitter(width = 0.2, shape=1)+
  scale_colour_manual(values=c("deepskyblue4", "indianred"))+
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
  ggplot(., aes(Day, value, colour=Day))+
  geom_boxplot(outlier.shape = NA)+
  geom_line(aes(group=ID), colour="gray90")+
  geom_jitter(position=position_jitterdodge(), shape=1, size=1)+
  scale_colour_manual(values=c("deepskyblue4", "indianred"))+
  labs(x="", y="Max covariance (z-score)")+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  scale_y_continuous(limits = c(0.3, 1))+
  theme(legend.title = element_blank(), legend.position = "none")

p_cov_bar <- plot_grid(p_cov_mean, p_cov_best, nrow = 1)

## plot the neurons trace of m3 mice----
p_trace_m3 <- read.xlsx("~cchen/Documents/neuroscience/Pn project/Data_analysis/miniscope/matlab_analysis/m3/trace_days//m3_trace_d3.xlsx", colNames = F, rowNames = F) %>% 
  filter(., X1 ==1) %>% 
  select(., -c(1,2)) %>% 
  t() %>% 
  apply(., 2, scale) %>% 
  as_tibble() %>% 
  select(sample(1: ncol(.), 10)) %>% 
  slice(1:2400) %>% 
  mutate(Time = seq(from = 0.05, to = 120, by = 0.05)) %>% 
  pivot_longer(-Time, ) %>% 
  ggplot(., aes(Time, value, color = name))+
  geom_line()+
  facet_grid(rows = vars(name))+
  theme_void()+
  theme(legend.position = "none")+
  theme(strip.text.y = element_blank())
  # annotate(x=c(110,110,120,120), y=c(0),"path")+
  # annotate(x=c(120), y=c(0,0,5,5),"path")








p_trace_m18_d7 <- read.xlsx("~cchen/Documents/neuroscience/Pn project/Data_analysis/miniscope/matlab_analysis/m18/trace_days//m18_trace_d7.xlsx", colNames = F, rowNames = F) %>% 
  filter(., X1 ==1) %>% 
  select(., -c(1,2)) %>% 
  t() %>% 
  apply(., 2, scale) %>% 
  as_tibble() %>% 
  select(1:20) %>% 
  slice(24:1224) %>% 
  mutate(Time = seq(from = 0, to = 60, by = 0.05)) %>% 
  pivot_longer(-Time, ) %>% 
  ggplot(., aes(Time, value, color = name))+
  geom_line()+
  facet_grid(rows = vars(name))+
  theme_void()+
  theme(legend.position = "none")+
  geom_vline(xintercept = 5, linetype="dotted", 
             color = "red", size=.5)

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
