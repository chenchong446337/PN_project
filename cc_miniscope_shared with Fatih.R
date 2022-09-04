## share with Fatih

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

# for m855 (for m855 and m857, the behavior videos were recorded at 10 Hz, while ca2+ imaging at 20Hz. That's why below I used frame*2)
path_trace_m855 <- as.list(list.files(path ="~cchen/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/matlab_analysis/m855/trace_days", pattern = "*.xlsx", full.names = T ))[c(3,7)]
t_stim_m855 <- list(t_stim_m855_d3 = c(392) *2, t_stim_m855_d7 = c(41)*2)
dat_trace_m855 <- mapply(c_miniscope_matlab, path_trace_m855, t_stim_m855, SIMPLIFY = F)

# for m857
path_trace_m857 <- as.list(list.files(path ="~cchen/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/matlab_analysis/m857/trace_days", pattern = "*.xlsx", full.names = T ))[c(1, 4)]
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


heat_m_ID <- "m18"
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

setwd("~cchen2/Documents/neuroscience/Pn\ project/Figure/PDF/")
cairo_pdf("p_trace.pdf", width = 90/25.6, height = 62/25.6, family = "Arial")
p_trace
dev.off()


## calculate the sum of active and inhibited trace----
dat_cell_trace_sum <- rbind(dat_cell_trace_pre_sta, dat_cell_trace_test_sta) %>%
  as_tibble() %>% 
  mutate(Day = rep(c("Pre",  "Test"), c(nrow(dat_cell_trace_pre_sta),nrow(dat_cell_trace_test_sta)))) %>%
  subset(., .$Group!="Neutral") %>%
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

setwd("~cchen2/Documents/neuroscience/Pn\ project/Figure/PDF/")
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



setwd("~cchen2/Documents/neuroscience/Pn\ project/Figure/PDF/")
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
  
  ctrl_rang <- range(t_stim-22, t_stim-2)
  test_rang <- range(t_stim, t_stim +40)
  
  ctrl_freq <- length(dat_peak[dat_peak>= ctrl_rang[1] & dat_peak < ctrl_rang[2]])/1 # 2s before crossing
  test_freq <- length(dat_peak[dat_peak>= test_rang[1] & dat_peak < test_rang[2]])/2 # 4s after crossing
  diff_freq <- test_freq - ctrl_freq
  
  return(c(ctrl_freq, test_freq, diff_freq))
}



