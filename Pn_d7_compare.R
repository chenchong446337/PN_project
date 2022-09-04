## This script is used to analyze the Ca2+ transit change during anticipation violation
## created on 04202020
## function to extract trace crossing the border
## updated on 06182020, 1. timecourse, and add more analysis


c_miniscope_matlab_d7 <- function(ID_trace, t_stim) {
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
  
  return(dat_stim)
  
}
## for m3
path_trace_m3 <- "~cchen2/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/matlab_analysis/m3/trace_days/m3_trace_d7.xlsx"
t_stim_m3_d7 <- c(346, 1206)

dat_stim_m3_d7 <- c_miniscope_matlab_d7(path_trace_m3, t_stim_m3_d7)


## for m7
path_trace_m7 <- "~cchen2/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/matlab_analysis/m7/trace_days/m7_trace_d7.xlsx"
t_stim_m7_d7 <- c(236, 1914)

dat_stim_m7_d7 <- c_miniscope_matlab_d7(path_trace_m7, t_stim_m7_d7)

## for m17
path_trace_m17 <- "~cchen2/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/matlab_analysis/m17/trace_days/m17_trace_d7.xlsx"
t_stim_m17_d7 <- c(157, 2113)

dat_stim_m17_d7 <- c_miniscope_matlab_d7(path_trace_m17, t_stim_m17_d7)


## for m18
path_trace_m18 <- "~cchen2/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/matlab_analysis/m18/trace_days/m18_trace_d7.xlsx"
t_stim_m18_d7 <- c(124,2238)

dat_stim_m18_d7 <- c_miniscope_matlab_d7(path_trace_m18, t_stim_m18_d7)

## for m855
path_trace_m855 <- "~cchen2/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/matlab_analysis/m855/trace_days/m855_trace_d7.xlsx"
t_stim_m855_d7 <- c(82,1572)

dat_stim_m855_d7 <- c_miniscope_matlab_d7(path_trace_m855, t_stim_m855_d7)

## for m857
path_trace_m857 <- "~cchen2/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/matlab_analysis/m857/trace_days/m857_trace_d7.xlsx"
t_stim_m857_d7 <- c(334,1892)

dat_stim_m857_d7 <- c_miniscope_matlab_d7(path_trace_m857, t_stim_m857_d7)


## combine data and do k-means analysis-----
dat_cell_trace <- vector(mode = "list", 2)
for (i in 1:2) {
  dat_cell_trace[[i]]<- cbind(dat_stim_m3_d7[[i]], dat_stim_m7_d7[[i]], dat_stim_m17_d7[[i]], dat_stim_m18_d7[[i]],dat_stim_m855_d7[[i]], dat_stim_m857_d7[[i]])
  
}

dat_cell_trace_re <- vector(mode = "list", 2)
for (i in 1:length(dat_cell_trace)) {
  dat_cell_trace_d <- dat_cell_trace[[i]]
  cell_id <- sprintf("Cell%s",seq(1:ncol(dat_cell_trace_d)))
  colnames(dat_cell_trace_d)<- cell_id
  
  dat_cell_trace_d_cluster <- unname(kmeans(t(dat_cell_trace_d), centers = 3, nstart = 50)[[1]])
  dat_cell_trace_d<- as.data.frame(dat_cell_trace_d)
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
  rep_time <- c(ncol(dat_stim_m3_d7[[i]]), ncol(dat_stim_m7_d7[[i]]), ncol(dat_stim_m17_d7[[i]]), ncol(dat_stim_m18_d7[[i]]),ncol(dat_stim_m855_d7[[i]]), ncol(dat_stim_m857_d7[[i]]))
  dat_cell_trace_d_re$ID <- rep(c("m3", "m7", "m17", "m18", "m855", "857"), rep_time*length(stim_time))
  
  dat_cell_trace_re[[i]] <- dat_cell_trace_d_re
}

score_range <- range(mapply(function(x) x$value, dat_cell_trace_re, SIMPLIFY = T))
## align by the cell ID

dat_trace_sta <- dat_cell_trace_re[[1]] %>%
  ddply(.,.(variable, Group),summarise,mean=mean(value), sum=sum(value)) %>%
  .[order(.[,'mean']),]

for (i in c(1, 2)) {
  dat_trace <- dat_cell_trace_re[[i]]
  dat_trace$variable <- factor(dat_trace$variable, levels = dat_trace_sta$variable)
  p_heat <- ggplot(dat_trace, aes(Time, variable,fill= value))+ 
    geom_tile(height=2)+
    # facet_grid(rows = vars(Group), scales = "free_y")+
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

p_heat_com <- plot_grid(p_heat_d1, p_heat_d2, nrow = 1)

setwd("~cchen2/Documents/neuroscience/Pn\ project/Figure/PDF/")
cairo_pdf("p_anti_heat_d7_2.pdf", width = 100/25.6, height = 60/25.6, family = "Arial")
p_heat_com
dev.off()

dat_cell_trace_pre_sta <- ddply(dat_cell_trace_re[[1]], .(Time, Group), summarise,n=length(value),mean=mean(value),sd=sd(value),se=sd(value)/sqrt(length(value)))
dat_cell_trace_test_sta <- ddply(dat_cell_trace_re[[2]], .(Time, Group), summarise,n=length(value),mean=mean(value),sd=sd(value),se=sd(value)/sqrt(length(value)))

dat_cell_trace_sta <- rbind(dat_cell_trace_pre_sta,  dat_cell_trace_test_sta)

dat_cell_trace_sta$Day <- rep(c("Pre",  "Test"),each= length(stim_time)*3)
dat_cell_trace_sta$Day <- factor(dat_cell_trace_sta$Day, levels = c("Pre",  "Test"))
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
  theme(legend.title = element_blank())

## plot sum of E-I
mouse_ID <- c("m3", "m7", "m17", "m18", "m855", "m857")
dat_cell_trace_pre1 <- ddply(dat_cell_trace_re[[1]], .(ID,Time, Group), summarise,value=mean(value))
dat_cell_trace_test1 <- ddply(dat_cell_trace_re[[2]], .(ID,Time, Group), summarise,value=mean(value))

dat_cell_trace_sum <-  rbind(dat_cell_trace_pre1, dat_cell_trace_test1) %>% 
  mutate( Day = rep(c("Pre",  "Test"), c(nrow(dat_cell_trace_pre1), nrow(dat_cell_trace_test1)))) %>%
  subset(., .$Group!="Neutral") %>%
  ddply(., .(ID,Time, Day), summarise, value=sum(value)) %>%
  ddply(., .(Time, Day), summarise,n=length(value),mean=mean(value),sd=sd(value),se=sd(value)/sqrt(length(value))) %>%
  mutate(Day = factor(Day, levels = c("Pre", "Test")))

p_trace_sum <- ggplot(dat_cell_trace_sum, aes(Time, mean, colour=Day))+
  geom_line()+
  geom_ribbon(aes(ymin=mean-se, ymax=mean+se, fill=Day), alpha=0.1, linetype=0)+
  scale_colour_manual(values=c("seagreen", "indianred"))+
  scale_fill_manual(values=c("seagreen", "indianred"))+
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
cairo_pdf("p_trace_d7_2.pdf", width = 130/25.6, height = 65/25.6, family = "Arial")
p_trace
dev.off()

setwd("~cchen2/Documents/neuroscience/Pn\ project/Figure/PDF/")
cairo_pdf("p_trace_d7_sum.pdf", width = 80/25.6, height = 65/25.6, family = "Arial")
p_trace_sum
dev.off()


dat_cell_area <- c()
cell_area_day <- c("Pre",  "Test")
for (i in c(1:2)){
  rep_time <- c(ncol(dat_stim_m3_d7[[i]]), ncol(dat_stim_m7_d7[[i]]), ncol(dat_stim_m17_d7[[i]]), ncol(dat_stim_m18_d7[[i]]), ncol(dat_stim_m855_d7[[i]]),ncol(dat_stim_m857_d7[[i]]))
  dat_trace <- dat_cell_trace_re[[i]]
  dat_trace$ID <- rep(c("m3", "m7", "m17", "m18", "m855", "m857"), rep_time*length(stim_time))
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

dat_cell_area <- dat_cell_area %>% 
  replace(is.na(.), 0) %>% 
  mutate(sum = Excited + Inhibited) %>% 
  mutate(Day = factor(Day, levels = c("Pre", "Test")))
dat_cell_area_sta <- dat_cell_area %>% 
  ddply(.,.(Day), summarise,n=length(ratio),mean=mean(ratio),sd=sd(ratio),se=sd(ratio)/sqrt(length(ratio)))


p_EI_ratio <- dat_cell_area %>% 
  ggplot(., aes(Day, ratio, colour=Day))+
  geom_boxplot(outlier.shape = NA)+
  geom_line(aes(group=ID), colour="gray90")+
  geom_jitter(width = 0.2, shape=1)+
  labs(x="", y="E/I ratio")+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  theme(legend.title = element_blank())+
  theme(legend.position = "none")

p_EI_ratio <- dat_cell_area %>% 
  ggplot(., aes(Day, sum, colour=Day))+
  geom_boxplot(outlier.shape = NA)+
  geom_line(aes(group=ID), colour="gray90")+
  geom_jitter(width = 0.2, shape=1)+
  labs(x="", y="E/I ratio")+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  theme(legend.title = element_blank())+
  theme(legend.position = "none")


setwd("~cchen2/Documents/neuroscience/Pn\ project/Figure/PDF/")
cairo_pdf("p_EI_d7_comp_2.pdf", width = 42/25.6, height = 65/25.6, family = "Arial")
p_EI_ratio
dev.off()

## for only m3, m7, stim1 to stim5-----
c_miniscope_matlab_d7 <- function(ID_trace, t_stim) {
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
  
  return(dat_stim)
  
}
## for m3
path_trace_m3 <- "~cchen2/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/matlab_analysis/m3/trace_days/m3_trace_d7.xlsx"
t_stim_m3_d7 <- c(346, 1206, 1756, 2288)

dat_stim_m3_d7 <- c_miniscope_matlab_d7(path_trace_m3, t_stim_m3_d7)


## for m7
path_trace_m7 <- "~cchen2/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/matlab_analysis/m7/trace_days/m7_trace_d7.xlsx"
t_stim_m7_d7 <- c(236, 1914,2777,3334)

dat_stim_m7_d7 <- c_miniscope_matlab_d7(path_trace_m7, t_stim_m7_d7)



dat_cell_trace <- vector(mode = "list", 4)
for (i in 1:4) {
  dat_cell_trace[[i]]<- cbind(dat_stim_m3_d7[[i]], dat_stim_m7_d7[[i]] )
  
}

dat_cell_trace_re <- vector(mode = "list", 4)
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
  rep_time <- c(ncol(dat_stim_m3_d7[[i]]), ncol(dat_stim_m7_d7[[i]]))
  dat_cell_trace_d_re$ID <- rep(c("m3", "m7"), rep_time*length(stim_time))
  
  dat_cell_trace_re[[i]] <- dat_cell_trace_d_re
}

score_range <- range(mapply(function(x) x$value, dat_cell_trace_re, SIMPLIFY = T))

for (i in c(1:4)) {
  dat_trace <- dat_cell_trace_re[[i]]
  dat_trace_sta <- ddply(dat_trace, .(variable, Group), summarise,mean=mean(value), sum=sum(value))
  dat_trace_sta <- dat_trace_sta[order(dat_trace_sta[,'mean']),]
  dat_trace$variable <- factor(dat_trace$variable, levels = dat_trace_sta$variable)
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
  assign(paste0("p_d7_stim", i), p_heat)
}

p_heat_com <- plot_grid(p_d7_stim1, p_d7_stim2, p_d7_stim3, p_d7_stim4, nrow = 1)


dat_cell_trace_stim1_sta <- ddply(dat_cell_trace_re[[1]], .(Time, Group), summarise,n=length(value),mean=mean(value),sd=sd(value),se=sd(value)/sqrt(length(value)))
dat_cell_trace_stim2_sta <- ddply(dat_cell_trace_re[[2]], .(Time, Group), summarise,n=length(value),mean=mean(value),sd=sd(value),se=sd(value)/sqrt(length(value)))
dat_cell_trace_stim3_sta <- ddply(dat_cell_trace_re[[3]], .(Time, Group), summarise,n=length(value),mean=mean(value),sd=sd(value),se=sd(value)/sqrt(length(value)))
dat_cell_trace_stim4_sta <- ddply(dat_cell_trace_re[[4]], .(Time, Group), summarise,n=length(value),mean=mean(value),sd=sd(value),se=sd(value)/sqrt(length(value)))

dat_cell_trace_sta <- rbind(dat_cell_trace_stim1_sta,  dat_cell_trace_stim2_sta, dat_cell_trace_stim3_sta, dat_cell_trace_stim4_sta)

dat_cell_trace_sta$Day <- rep(c("stim1",  "stim2", "stim3", "stim4"),each= length(stim_time)*3)
dat_cell_trace_sta$Day <- factor(dat_cell_trace_sta$Day, levels = c("stim1",  "stim2", "stim3", "stim4"))
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


dat_cell_area <- c()
cell_area_day <- c("stim1",  "stim2", "stim3", "stim4")
for (i in c(1:4)){
  rep_time <- c(ncol(dat_stim_m3_d7[[i]]), ncol(dat_stim_m7_d7[[i]]))
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

dat_cell_area$Day <- factor(dat_cell_area$Day, levels = c("stim1","stim2","stim3","stim4"))
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


## compare corssing back(08132020)--------
c_miniscope_matlab_d7 <- function(ID_trace, t_stim) {
  dat_trace <- read.xlsx(ID_trace, colNames = F, rowNames = F) %>% 
    filter(., X1 ==1) %>% 
    select(., -c(1,2)) %>% 
    t() %>% 
    as.data.frame()
  
  rownames(dat_trace) <- NULL
  colnames(dat_trace)<- NULL
  ## do z score of the whole trace
  dat_trace <- apply(dat_trace, 2, scale)
  
  ## divide ca2+ trace within -5 and 9s before and after stimuls, take -5 to -3 as baseline
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
  
  return(dat_stim)
  
}
## for m3
path_trace_m3 <- "~cchen/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/matlab_analysis/m3/trace_days/m3_trace_d7.xlsx"
t_stim_m3_d7 <- c(346, 1003)

dat_stim_m3_d7 <- c_miniscope_matlab_d7(path_trace_m3, t_stim_m3_d7)


## for m7
path_trace_m7 <- "~cchen/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/matlab_analysis/m7/trace_days/m7_trace_d7.xlsx"
t_stim_m7_d7 <- c(236, 1632)

dat_stim_m7_d7 <- c_miniscope_matlab_d7(path_trace_m7, t_stim_m7_d7)

## for m17
path_trace_m17 <- "~cchen/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/matlab_analysis/m17/trace_days/m17_trace_d7.xlsx"
t_stim_m17_d7 <- c(157, 2113)

dat_stim_m17_d7 <- c_miniscope_matlab_d7(path_trace_m17, t_stim_m17_d7)


## for m18
path_trace_m18 <- "~cchen/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/matlab_analysis/m18/trace_days/m18_trace_d7.xlsx"
t_stim_m18_d7 <- c(124,1923)

dat_stim_m18_d7 <- c_miniscope_matlab_d7(path_trace_m18, t_stim_m18_d7)

## for m855
path_trace_m855 <- "~cchen/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/matlab_analysis/m855/trace_days/m855_trace_d7.xlsx"
t_stim_m855_d7 <- c(84,1572)

dat_stim_m855_d7 <- c_miniscope_matlab_d7(path_trace_m855, t_stim_m855_d7)

## for m857
path_trace_m857 <- "~cchen/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/matlab_analysis/m857/trace_days/m857_trace_d7.xlsx"
t_stim_m857_d7 <- c(334,1892)

dat_stim_m857_d7 <- c_miniscope_matlab_d7(path_trace_m857, t_stim_m857_d7)


## combine data and do k-means analysis
dat_cell_trace <- vector(mode = "list", 2)
for (i in 1:2) {
  dat_cell_trace[[i]]<- cbind(dat_stim_m3_d7[[i]], dat_stim_m7_d7[[i]], dat_stim_m17_d7[[i]], dat_stim_m18_d7[[i]], dat_stim_m855_d7[[i]],dat_stim_m857_d7[[i]] )
  
}

dat_cell_trace_re <- vector(mode = "list", 2)
mouse_ID <- c("m3", "m7", "m17", "m18", "m855", "m857")

for (i in 1:length(dat_cell_trace)) {
  dat_cell_trace_d <- dat_cell_trace[[i]]
  cell_id <- sprintf("Cell%s",seq(1:ncol(dat_cell_trace_d)))
  colnames(dat_cell_trace_d)<- cell_id
  
  dat_cell_trace_d_cluster <- unname(kmeans(t(dat_cell_trace_d), centers = 3)[[1]])
  dat_cell_trace_d<- as.data.frame(dat_cell_trace_d)
  stim_time <- seq(-2, 6.5, by=0.5)
  dat_cell_trace_d$Time <- stim_time
  dat_cell_trace_d_re <- melt(dat_cell_trace_d, id.vars ='Time')
  dat_cell_trace_d_re$Group <- rep(dat_cell_trace_d_cluster, each=length(stim_time))
  
  ## sort data by the value in each group
  dat_cell_d_sort <- as.numeric(names(sort(tapply(dat_cell_trace_d_re$value, dat_cell_trace_d_re$Group, mean), decreasing = T)))
  
  dat_cell_trace_d_re$Group[dat_cell_trace_d_re$Group == dat_cell_d_sort[1]] ="Excited"
  dat_cell_trace_d_re$Group[dat_cell_trace_d_re$Group == dat_cell_d_sort[2]] ="Neutral"
  dat_cell_trace_d_re$Group[dat_cell_trace_d_re$Group == dat_cell_d_sort[3]] ="Inhibited"
  
  dat_cell_trace_d_re$Group <- factor(dat_cell_trace_d_re$Group, levels = c("Excited", "Neutral", "Inhibited"))
  rep_time <- c(ncol(dat_stim_m3_d7[[i]]), ncol(dat_stim_m7_d7[[i]]), ncol(dat_stim_m17_d7[[i]]), ncol(dat_stim_m18_d7[[i]]), ncol(dat_stim_m855_d7[[i]]), ncol(dat_stim_m857_d7[[i]]))
  dat_cell_trace_d_re$ID <- rep(mouse_ID, rep_time*length(stim_time))
  
  dat_cell_trace_re[[i]] <- dat_cell_trace_d_re
}

## align by cells
score_range <- range(mapply(function(x) x$value, dat_cell_trace_re, SIMPLIFY = T))
dat_trace_sta <- dat_cell_trace_re[[1]] %>%
  ddply(.,.(variable, Group),summarise,mean=mean(value), sum=sum(value)) %>%
  arrange(mean)
for (i in c(1, 2)) {
  dat_trace <- dat_cell_trace_re[[i]]
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

p_heat_com_d7 <- plot_grid(p_heat_d1, p_heat_d2, nrow = 1)

setwd("~cchen2/Documents/neuroscience/Pn\ project/Figure/PDF/")
cairo_pdf("p_anti_heat_d7.pdf", width = 90/25.6, height = 55/25.6, family = "Arial")
p_heat_com_d7
dev.off()

dat_cell_trace_pre_sta <-  ddply(dat_cell_trace_re[[1]], .(ID,Time, Group), summarise,value=mean(value)) 
dat_cell_trace_test_sta <- ddply(dat_cell_trace_re[[2]], .(ID,Time, Group), summarise,value=mean(value))

dat_cell_trace_sta <- ddply(dat_cell_trace_pre_sta, .(Time, Group), summarise,n=length(value),mean=mean(value),sd=sd(value),se=sd(value)/sqrt(length(value))) %>%
  rbind(., ddply(dat_cell_trace_test_sta, .(Time, Group), summarise,n=length(value),mean=mean(value),sd=sd(value),se=sd(value)/sqrt(length(value)))) %>%
  mutate(., Day = rep(c("Pre",  "Test"),each= length(stim_time)*3)) %>% 
  mutate(Day = factor(Day, levels = c("Pre",  "Test")), Group = factor(Group,levels = c("Neutral", "Excited", "Inhibited") ))

range_trace_plot <- range(dat_cell_trace_sta$mean)

## group by day
p_trace <-ggplot(dat_cell_trace_sta, aes(Time, mean, colour=Group))+
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
cairo_pdf("p_trace_d7.pdf", width = 90/25.6, height = 62/25.6, family = "Arial")
p_trace
dev.off()

## calculate the sum of active and inhibited trace
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
  scale_colour_manual(values=c("indianred", "darkslategray"))+
  scale_fill_manual(values=c("indianred", "darkslategray"))+
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

setwd("~cchen2/Documents/neuroscience/Pn\ project/Figure/PDF/")
cairo_pdf("p_trace_sum_d7.pdf", width = 80/25.6, height = 60/25.6, family = "Arial")
p_trace_sum
dev.off()


dat_cell_area <- c()
cell_area_day <- c("Pre",  "Test")
for (i in c(1:2)){
  rep_time <- c(ncol(dat_stim_m3_d7[[i]]), ncol(dat_stim_m7_d7[[i]]), ncol(dat_stim_m17_d7[[i]]), ncol(dat_stim_m18_d7[[i]]), ncol(dat_stim_m855_d7[[i]]), ncol(dat_stim_m857_d7[[i]]))
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

dat_cell_area$Day <- factor(dat_cell_area$Day, levels = c("Pre","Test"))
dat_cell_area_com <- ddply(dat_cell_area, .(ID, Day), summarise, value=mean(ratio, na.rm = T))
dat_cell_area_sta <- ddply(dat_cell_area_com,.(Day), summarise,n=length(value),mean=mean(value),sd=sd(value),se=sd(value)/sqrt(length(value)))


dat_cell_area$Day <- factor(dat_cell_area$Day, levels = c("Pre","Test"))

p_EI_ratio <- dat_cell_area %>% 
  select(ID,Day, ratio) %>% 
  ggplot(., aes(Day, ratio, color=Day))+
  geom_boxplot(outlier.shape = NA)+
  geom_line(aes(group = ID), color="gray90")+
  geom_jitter(width = 0.2, shape=1)+
  scale_colour_manual(values=c("indianred", "darkslategray"))+
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
annotate(x=c(1,1,2,2), y=c(3.5,3.6,3.6,3.5),"path")+
annotate("text",x=1.5,y=3.6, label="*", size=5)


## statistic test
t_EI<-dat_cell_area %>% 
  select(Day, ID, ratio) %>% 
  t.test(ratio~Day, ., paired = TRUE)
setwd("~cchen2/Documents/neuroscience/Pn\ project/Figure/PDF/")
cairo_pdf("p_EI_d7_comp.pdf", width = 40/25.6, height = 65/25.6, family = "Arial")
p_EI_ratio
dev.off()


## plot the sum of E-I

p_EI_sum<- dat_cell_area %>% 
  select(Day, ID, sum) %>% 
  ggplot(., aes(Day, sum, group=Day))+
  geom_boxplot(outlier.shape = NA)+
  geom_line(aes(group=ID), colour="gray90")+
  geom_jitter(aes(colour = Day, shape = Day),width = 0.2, size = 2)+
  scale_colour_manual(values=c("indianred", "darkslategray"))+
  scale_shape_manual(values=c(20, 18))+
  labs(x="", y="AUG (z-score)")+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  theme(legend.title = element_blank(), legend.position = "none")+
annotate(x=c(1,1,2,2), y=c(12,12.5,12.5,12),"path")+
annotate("text",x=1.5,y=12.5, label="**", size=5)

## statistic test
t_sum<-dat_cell_area %>% 
  select(Day, ID, sum) %>% 
  wilcox.test(sum~Day,.)
  t.test(sum~Day, ., paired = TRUE)

setwd("~cchen2/Documents/neuroscience/Pn\ project/Figure/PDF/")
cairo_pdf("p_sum_d7_comp.pdf", width = 40/25.6, height = 65/25.6, family = "Arial")
p_EI_sum
dev.off()

p_EIS_change <- dat_cell_area %>%
  select("Excited", "Inhibited","ID","Day", "sum") %>%
  melt(., id.vars= c('Day', "ID")) %>%
  ddply(., .(ID, Day, variable), summarise, value=mean(value, na.rm = T)) %>%
  ggplot(., aes(interaction(Day, variable), value, colour=Day))+
  geom_boxplot(outlier.shape = NA)+
  geom_line(aes(group=interaction(ID, variable)), colour="gray90") +
  geom_jitter(width = 0.2, shape=1, size=1)+
  scale_colour_manual(values=c("indianred", "darkslategray" ))+
  labs(x="", y="AUC (z-score)")+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  theme(legend.title = element_blank(), legend.position = "none")

setwd("~cchen2/Documents/neuroscience/Pn\ project/Figure/PDF/")
cairo_pdf("p_EIS_change_d7.pdf", width = 95/25.6, height = 60/25.6, family = "Arial")
p_EIS_change
dev.off()

dat_cell_cat <- mapply(function (x) subset(x, x$Time==0), dat_cell_trace_re, SIMPLIFY = F ) 
p_cat<- mapply(function(x) prop.table(table(x$ID, x$Group), 1), dat_cell_cat, SIMPLIFY = F) %>%
  do.call(rbind, .) %>%
  as.data.frame() %>%
  mutate(ID = rep(rownames(.)[1:length(mouse_ID)], 2), Day=rep(c("Pre", "Test"), each=length(mouse_ID))) %>%
  melt(., id.vars = c("ID", "Day")) %>%
  ddply(., .(ID, Day, variable), summarise, "Prop"=mean(value)) %>%
  mutate(Day = factor(Day, levels = c("Pre", "Test"))) %>%
  mutate(variable=factor(variable, levels = c("Excited", "Neutral", "Inhibited"))) %>% 
  ggplot(., aes(interaction(Day,variable), Prop, colour=Day))+
  geom_boxplot(outlier.shape = NA)+
  geom_line(aes(group=interaction(ID, variable)), colour="gray90")+
  geom_jitter(position=position_jitterdodge(), shape=1, size=1)+
  scale_colour_manual(values=c( "indianred", "darkslategray"))+
  labs(x="", y="% of all neurons")+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  scale_y_continuous(limits = c(0,1),expand = c(0,0))+
  theme(legend.title = element_blank(), legend.position = "none")

setwd("~cchen2/Documents/neuroscience/Pn\ project/Figure/PDF/")
cairo_pdf("p_cat_d7.pdf", width = 90/25.6, height = 60/25.6, family = "Arial")
p_cat
dev.off()
## plot by aligning the cells to see if they have opposite rule
dat_cell_trace_pre_sta <-  ddply(dat_cell_trace_re[[1]], .(ID,Time, Group), summarise,value=mean(value)) 
dat_cell_trace_test_sta <- dat_cell_trace_re[[2]] %>% 
  mutate(Group = dat_cell_trace_re[[1]]$Group ) %>% 
  ddply(., .(ID,Time, Group), summarise,value=mean(value))

dat_cell_trace_sta <- ddply(dat_cell_trace_pre_sta, .(Time, Group), summarise,n=length(value),mean=mean(value),sd=sd(value),se=sd(value)/sqrt(length(value))) %>%
  rbind(., ddply(dat_cell_trace_test_sta, .(Time, Group), summarise,n=length(value),mean=mean(value),sd=sd(value),se=sd(value)/sqrt(length(value)))) %>%
  mutate(., Day = rep(c("Pre",  "Test"),each= length(stim_time)*3)) %>% 
  mutate(Day = factor(Day, levels = c("Pre",  "Test")), Group = factor(Group,levels = c("Neutral", "Excited", "Inhibited") ))

range_trace_plot <- range(dat_cell_trace_sta$mean)

p_trace <-ggplot(dat_cell_trace_sta, aes(Time, mean, colour=Day))+
  geom_line()+
  geom_ribbon(aes(ymin=mean-se, ymax=mean+se, fill=Day), alpha=0.1, linetype=0)+
  facet_grid(cols = vars(Group))+
  scale_colour_manual(values=c("indianred", "darkslategray"))+
  scale_fill_manual(values=c("indianred", "darkslategray"))+
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
cairo_pdf("p_trace_d7_trace.pdf", width = 120/25.6, height = 65/25.6, family = "Arial")
p_trace
dev.off()

## calculate the spikes before and after crossing------
cc_firing_fun <- function(path_trace, path_peak, t_stim) {
  cell_valid <- read.xlsx(path_trace, colNames = F, rowNames = F) %>% 
    select (., X1)
  
  dat_peak <- read.csv(path_peak) %>% 
    select(which(cell_valid==1)) %>% 
    unlist() %>% 
    na.omit()
  
  ctrl_freq <- rep(0, 2)
  test_freq <- rep(0, 2)
  
  for (i in seq_along(t_stim)){
    ctrl_rang <- range(t_stim[i]-22, t_stim[i]-2)
    test_rang <- range(t_stim[i], t_stim[i] +40)
    
    ctrl_freq[i] <- length(dat_peak[dat_peak>= ctrl_rang[1] & dat_peak < ctrl_rang[2]])/1 # 2s before crossing
    test_freq[i] <- length(dat_peak[dat_peak>= test_rang[1] & dat_peak < test_rang[2]])/2 # 4s after crossing

  }
  
  dat_rate <- data.frame(Before= ctrl_freq, After = test_freq)
  return(dat_rate)
}

## for m3
path_peak_m3 <- list.files(path ="~cchen2/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/matlab_analysis/m3/peaks_days", pattern = "*.csv", full.names = T )[7]
dat_rate_m3 <- cc_firing_fun(path_trace_m3, path_peak_m3, t_stim_m3_d7)
## for m7
path_peak_m7 <- list.files(path ="~cchen2/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/matlab_analysis/m7/peaks_days/", pattern = "*.csv", full.names = T )[7]
dat_rate_m7 <- cc_firing_fun(path_trace_m7, path_peak_m7, t_stim_m7_d7)

## for m17
path_peak_m17 <- list.files(path ="~cchen2/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/matlab_analysis/m17/peaks_days/", pattern = "*.csv", full.names = T )[7]
dat_rate_m17 <- cc_firing_fun(path_trace_m17, path_peak_m17, t_stim_m17_d7)
## for m18
path_peak_m18 <- list.files(path ="~cchen2/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/matlab_analysis/m18/peaks_days/", pattern = "*.csv", full.names = T )[7]
dat_rate_m18 <- cc_firing_fun(path_trace_m18, path_peak_m18, t_stim_m18_d7)
## for m855
path_peak_m855 <- list.files(path ="~cchen2/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/matlab_analysis/m855/peaks_days/", pattern = "*.csv", full.names = T )[7]
dat_rate_m855 <- cc_firing_fun(path_trace_m855, path_peak_m855, t_stim_m855_d7)

## combine data and plot
dat_firing <- rbind(dat_rate_m3, dat_rate_m7, dat_rate_m17, dat_rate_m18, dat_rate_m855) %>% 
  as_tibble() %>% 
  mutate(diff = After - Before) %>% 
  mutate(ID = rep(mouse_ID, each=2)) %>% 
  mutate(Group = rep(c("Cross", "Cross_back"), length(mouse_ID))) %>% 
  gather(variable, value, -ID,  -Group) %>% 
  ddply(., .(ID, Group, variable), summarise, value=mean(value)) %>% 
  mutate(Group= factor(Group, levels = c("Cross", "Cross_back")))

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

setwd("~cchen2/Documents/neuroscience/Pn\ project/Figure/PDF/")
cairo_pdf("p_firing.pdf", width = 40/25.6, height = 65/25.6, family = "Arial")
p_firing
dev.off() 

## covariance analysis
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

p_cov_d3<- ggcorrplot(dat_cell_cov_list[[1]][[1]], hc.order = TRUE, outline.col = "white")+
  scale_fill_gradient2(limit = c(cov_range[1], cov_range[2]), low = "navy", high =  "red4", mid = "white")+
  theme_void()+
  theme(legend.position = "none")

p_cov_d7<- ggcorrplot(dat_cell_cov_list[[2]][[1]], hc.order = TRUE, outline.col = "white")+
  scale_fill_gradient2(limit = c(cov_range[1], cov_range[2]), low = "navy", high =  "red4", mid = "white")+
  theme_void()+
  theme(legend.position = "none")
p_cov_com <- plot_grid(p_cov_d3, p_cov_d7, nrow = 1)

setwd("~cchen2/Documents/neuroscience/Pn\ project/Figure/PDF/")
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
dat_cell_cov$Day <- factor(dat_cell_cov$Day, levels = c("Pre",  "Test"))
p_cov_cum<- ggplot(dat_cell_cov, aes(value, group=Day, colour=Day))+
  stat_ecdf(geom = "step")+
  scale_color_manual(values=c( "deepskyblue4", "indianred"))+
  labs(x="", y="Cummulative probability")+
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
  mutate(Day = factor(Day, levels = c("Pre",  "Test")))


p_cov_mean<- subset(dat_cell_cov_sta1, dat_cell_cov_sta1$variable=="value") %>%
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
  theme(legend.title = element_blank(), legend.position = "none")

p_cov_bar <- plot_grid(p_cov_mean, p_cov_best, nrow = 1)

setwd("~cchen2/Documents/neuroscience/Pn\ project/Figure/PDF/")
cairo_pdf("p_cov_bar.pdf", width = 82/25.6, height = 65/25.6, family = "Arial")
p_cov_bar
dev.off()


## compare the first to last crossing-----
c_miniscope_matlab_d7 <- function(ID_trace, t_stim) {
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
  
  return(dat_stim)
  
}
## for m3
path_trace_m3 <- "~cchen2/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/matlab_analysis/m3/trace_days/m3_trace_d7.xlsx"
t_stim_m3_d7 <- c(346, 1756)

dat_stim_m3_d7 <- c_miniscope_matlab_d7(path_trace_m3, t_stim_m3_d7)


## for m7
path_trace_m7 <- "~cchen2/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/matlab_analysis/m7/trace_days/m7_trace_d7.xlsx"
t_stim_m7_d7 <- c(236, 2777)

dat_stim_m7_d7 <- c_miniscope_matlab_d7(path_trace_m7, t_stim_m7_d7)


## for m17
path_trace_m17 <- "~cchen2/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/matlab_analysis/m17/trace_days/m17_trace_d7.xlsx"
t_stim_m17_d7 <- c(157, 2536)

dat_stim_m17_d7 <- c_miniscope_matlab_d7(path_trace_m17, t_stim_m17_d7)

## for m18
path_trace_m18 <- "~cchen2/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/matlab_analysis/m18/trace_days/m18_trace_d7.xlsx"
t_stim_m18_d7 <- c(124, 2238)

dat_stim_m18_d7 <- c_miniscope_matlab_d7(path_trace_m18, t_stim_m18_d7)

## for m855
path_trace_m855 <- "~cchen2/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/matlab_analysis/m855/trace_days/m855_trace_d7.xlsx"
t_stim_m855_d7 <- c(82, 3062)

dat_stim_m855_d7 <- c_miniscope_matlab_d7(path_trace_m855, t_stim_m855_d7)

## for m857
path_trace_m857 <- "~cchen2/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/matlab_analysis/m857/trace_days/m857_trace_d7.xlsx"
t_stim_m857_d7 <- c(167*2, 1153*2)

dat_stim_m857_d7 <- c_miniscope_matlab_d7(path_trace_m857, t_stim_m857_d7)


dat_cell_trace <- vector(mode = "list", 2)
for (i in 1:2) {
  dat_cell_trace[[i]]<- cbind(dat_stim_m3_d7[[i]], dat_stim_m7_d7[[i]], dat_stim_m17_d7[[i]], dat_stim_m18_d7[[i]],dat_stim_m855_d7[[i]], dat_stim_m857_d7[[i]])
  
}

dat_cell_trace_re <- vector(mode = "list", 2)
for (i in 1:length(dat_cell_trace)) {
  dat_cell_trace_d <- dat_cell_trace[[i]]
  cell_id <- sprintf("Cell%s",seq(1:ncol(dat_cell_trace_d)))
  colnames(dat_cell_trace_d)<- cell_id
  
  dat_cell_trace_d_cluster <- unname(kmeans(t(dat_cell_trace_d), centers = 3)[[1]])
  dat_cell_trace_d<- as.data.frame(dat_cell_trace_d)
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
  rep_time <- c(ncol(dat_stim_m3_d7[[i]]), ncol(dat_stim_m7_d7[[i]]), ncol(dat_stim_m17_d7[[i]]), ncol(dat_stim_m18_d7[[i]]),ncol(dat_stim_m855_d7[[i]]), ncol(dat_stim_m857_d7[[i]]))
  dat_cell_trace_d_re$ID <- rep(c("m3", "m7", "m17", "m18", "m855", "857"), rep_time*length(stim_time))
  
  dat_cell_trace_re[[i]] <- dat_cell_trace_d_re
}

score_range <- range(mapply(function(x) x$value, dat_cell_trace_re, SIMPLIFY = T))

for (i in c(1:2)) {
  dat_trace <- dat_cell_trace_re[[i]]
  dat_trace_sta <- ddply(dat_trace, .(variable, Group), summarise,mean=mean(value), sum=sum(value))
  dat_trace_sta <- dat_trace_sta[order(dat_trace_sta[,'mean']),]
  dat_trace$variable <- factor(dat_trace$variable, levels = dat_trace_sta$variable)
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
  assign(paste0("p_d7_stim", i), p_heat)
}

p_heat_com <- plot_grid(p_d7_stim1, p_d7_stim2, nrow = 1)


dat_cell_trace_stim1_sta <- ddply(dat_cell_trace_re[[1]], .(ID, Time, Group), summarise,value=mean(value))
dat_cell_trace_stim2_sta <- ddply(dat_cell_trace_re[[2]], .(ID, Time, Group), summarise,value=mean(value))

dat_cell_trace_sta <- ddply(dat_cell_trace_stim1_sta, .(Time, Group), summarise,n=length(value),mean=mean(value),sd=sd(value),se=sd(value)/sqrt(length(value))) %>%
  rbind(., ddply(dat_cell_trace_stim2_sta, .(Time, Group), summarise,n=length(value),mean=mean(value),sd=sd(value),se=sd(value)/sqrt(length(value)))) %>%
  mutate(., Day = rep(c("stim1",  "stim2"),each= length(stim_time)*3)) %>% 
  mutate(Day = factor(Day, levels = c("stim1",  "stim2")), Group = factor(Group,levels = c("Neutral", "Excited", "Inhibited") ))

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
dat_cell_trace_sum <- rbind(dat_cell_trace_stim1_sta, dat_cell_trace_stim2_sta) %>%
  as_tibble() %>% 
  mutate(Day = rep(c("stim1",  "stim2"), c(nrow(dat_cell_trace_stim1_sta),nrow(dat_cell_trace_stim2_sta)))) %>%
  subset(., .$Group!="Neutral") %>%
  ddply(., .(ID,Time, Day), summarise, value=sum(value)) %>%
  ddply(., .(Time, Day), summarise,n=length(value),mean=mean(value),sd=sd(value),se=sd(value)/sqrt(length(value))) %>%
  mutate(Day = factor(Day, levels = c("stim1",  "stim2")))

p_trace_sum <- ggplot(dat_cell_trace_sum, aes(Time, mean, colour=Day))+
  geom_line()+
  geom_ribbon(aes(ymin=mean-se, ymax=mean+se, fill=Day), alpha=0.2, linetype=0)+
  scale_colour_manual(values=c( "indianred", "royalblue4"))+
  scale_fill_manual(values=c("indianred", "royalblue4"))+
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

setwd("~cchen2/Documents/neuroscience/Pn\ project/Figure/PDF/")
cairo_pdf("p_trace_sum_stim.pdf", width = 80/25.6, height = 60/25.6, family = "Arial")
p_trace_sum
dev.off()

mouse_ID <- c("m3", "m7", "m17", "m18", "m855", "m857")

dat_cell_area <- c()
group_day <- c("stim1",  "stim2")
for (i in c(1,2)){
  rep_time <- c(ncol(dat_stim_m3_d7[[i]]), ncol(dat_stim_m7_d7[[i]]), ncol(dat_stim_m17_d7[[i]]), ncol(dat_stim_m18_d7[[i]]), ncol(dat_stim_m855_d7[[i]]), ncol(dat_stim_m857_d7[[i]]))
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
  dat_cell_area <- rbind(dat_cell_area, dat_anti_area) %>% 
    replace_na(list(Excited = 0, Inhibited= 0, sum = 0)) 
}

dat_cell_area$Day <- factor(dat_cell_area$Day, levels = c("stim1","stim2"))

dat_EI_sum_sta <- dat_cell_area %>% 
  select(Day, ID, sum) %>% 
  ddply(., .(Day), summarise,n=length(sum),mean=mean(sum),sd=sd(sum),se=sd(sum)/sqrt(length(sum)))

p_sum <- dat_cell_area %>% 
  select(Day, ID, sum) %>% 
  wilcox.test(sum~Day, .)

p_EI_sum<- dat_cell_area %>% 
  select(Day, ID, sum) %>% 
  ggplot(., aes(Day, sum, group=Day))+
  geom_boxplot(outlier.shape = NA)+
  geom_line(aes(group=ID), colour="gray90")+
  geom_jitter(aes(colour = Day, shape = Day),width = 0.2, size = 2)+
  scale_colour_manual(values=c("indianred", "royalblue4" ))+
  scale_shape_manual(values=c(20, 18))+
  labs(x="", y="AUG (z-score)")+
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

setwd("~cchen2/Documents/neuroscience/Pn\ project/Figure/PDF/")
cairo_pdf("p_EI_sum_stim.pdf", width = 40/25.6, height = 65/25.6, family = "Arial")
p_EI_sum
dev.off()

