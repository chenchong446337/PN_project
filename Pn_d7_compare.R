## This script is used to analyze the Ca2+ transit change during anticipation violation
## created on 04202020
## function to extract trace crossing the border


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
t_stim_m3_d7 <- c(346, 1206)

dat_stim_m3_d7 <- c_miniscope_matlab_d7(path_trace_m3, t_stim_m3_d7)


## for m7
path_trace_m7 <- "~cchen2/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/matlab_analysis/m7/trace_days/m7_trace_d7.xlsx"
t_stim_m7_d7 <- c(236, 1914)

dat_stim_m7_d7 <- c_miniscope_matlab_d7(path_trace_m7, t_stim_m7_d7)

## for m17
path_trace_m17 <- "~cchen2/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/matlab_analysis/m17/trace_days/m17_trace_d7.xlsx"
t_stim_m17_d7 <- c(157, 2536)

dat_stim_m17_d7 <- c_miniscope_matlab_d7(path_trace_m17, t_stim_m17_d7)


## for m18
path_trace_m18 <- "~cchen2/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/matlab_analysis/m18/trace_days/m18_trace_d7.xlsx"
t_stim_m18_d7 <- c(124,2238)

dat_stim_m18_d7 <- c_miniscope_matlab_d7(path_trace_m18, t_stim_m18_d7)



## combine data and do k-means analysis-----
dat_cell_trace <- vector(mode = "list", 2)
for (i in 1:2) {
  dat_cell_trace[[i]]<- cbind(dat_stim_m3_d7[[i]], dat_stim_m7_d7[[i]], dat_stim_m17_d7[[i]], dat_stim_m18_d7[[i]] )
  
}

dat_cell_trace_re <- vector(mode = "list", 2)
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
  rep_time <- c(ncol(dat_stim_m3_d7[[i]]), ncol(dat_stim_m7_d7[[i]]), ncol(dat_stim_m17_d7[[i]]), ncol(dat_stim_m18_d7[[i]]))
  dat_cell_trace_d_re$ID <- rep(c("m3", "m7", "m17", "m18"), rep_time*length(stim_time))
  
  dat_cell_trace_re[[i]] <- dat_cell_trace_d_re
}

dat_cell_trace_re[[2]]$Group<-dat_cell_trace_re[[1]]$Group
score_range <- range(mapply(function(x) x$value, dat_cell_trace_re, SIMPLIFY = T))

for (i in c(1, 2)) {
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

setwd("~cchen2/Documents/neuroscience/Pn\ project/Figure/PDF/")
cairo_pdf("p_trace_d7_2.pdf", width = 130/25.6, height = 65/25.6, family = "Arial")
p_trace
dev.off()

dat_cell_area <- c()
cell_area_day <- c("Pre",  "Test")
for (i in c(1:2)){
  rep_time <- c(ncol(dat_stim_m3_d7[[i]]), ncol(dat_stim_m7_d7[[i]]), ncol(dat_stim_m17_d7[[i]]), ncol(dat_stim_m18_d7[[i]]))
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

dat_cell_area$Day <- factor(dat_cell_area$Day, levels = c("Pre","Test"))
dat_cell_area_com <- ddply(dat_cell_area, .(ID, Day), summarise, value=mean(ratio, na.rm = T))
dat_cell_area_sta <- ddply(dat_cell_area_com,.(Day), summarise,n=length(value),mean=mean(value),sd=sd(value),se=sd(value)/sqrt(length(value)))


p_EI_ratio <- ggplot(dat_cell_area_sta, aes(Day, mean, colour=Day, fill=Day))+
  geom_bar(stat="identity", color="black", position=position_dodge())+
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.1, position=position_dodge(.9), colour="black") +
  geom_jitter(data = dat_cell_area_com, aes(Day, value),width = 0.1, colour="black", shape=1 )+
  scale_fill_manual(values=c( "indianred4", "tomato3"))+
  labs(x="", y="E/I ratio")+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  scale_y_continuous(limits = c(0, 2), expand = c(0,0))+
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


## compare corssing back--------
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
t_stim_m3_d7 <- c(346, 1003)

dat_stim_m3_d7 <- c_miniscope_matlab_d7(path_trace_m3, t_stim_m3_d7)


## for m7
path_trace_m7 <- "~cchen2/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/matlab_analysis/m7/trace_days/m7_trace_d7.xlsx"
t_stim_m7_d7 <- c(236, 1632)

dat_stim_m7_d7 <- c_miniscope_matlab_d7(path_trace_m7, t_stim_m7_d7)

## for m17
path_trace_m17 <- "~cchen2/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/matlab_analysis/m17/trace_days/m17_trace_d7.xlsx"
t_stim_m17_d7 <- c(157, 2113)

dat_stim_m17_d7 <- c_miniscope_matlab_d7(path_trace_m17, t_stim_m17_d7)


## for m18
path_trace_m18 <- "~cchen2/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/matlab_analysis/m18/trace_days/m18_trace_d7.xlsx"
t_stim_m18_d7 <- c(124,1923)

dat_stim_m18_d7 <- c_miniscope_matlab_d7(path_trace_m18, t_stim_m18_d7)



## combine data and do k-means analysis-----
dat_cell_trace <- vector(mode = "list", 2)
for (i in 1:2) {
  dat_cell_trace[[i]]<- cbind(dat_stim_m3_d7[[i]], dat_stim_m7_d7[[i]], dat_stim_m17_d7[[i]], dat_stim_m18_d7[[i]] )
  
}

dat_cell_trace_re <- vector(mode = "list", 2)
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
  rep_time <- c(ncol(dat_stim_m3_d7[[i]]), ncol(dat_stim_m7_d7[[i]]), ncol(dat_stim_m17_d7[[i]]), ncol(dat_stim_m18_d7[[i]]))
  dat_cell_trace_d_re$ID <- rep(c("m3", "m7", "m17", "m18"), rep_time*length(stim_time))
  
  dat_cell_trace_re[[i]] <- dat_cell_trace_d_re
}

## align by cells
dat_cell_trace_re[[2]]$Group<-dat_cell_trace_re[[1]]$Group
score_range <- range(mapply(function(x) x$value, dat_cell_trace_re, SIMPLIFY = T))

for (i in c(1, 2)) {
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
  assign(paste0("p_heat_d", i), p_heat)
}

p_heat_com_d7 <- plot_grid(p_heat_d1, p_heat_d2, nrow = 1)
setwd("~cchen2/Documents/neuroscience/Pn\ project/Figure/PDF/")
cairo_pdf("p_anti_heat_d7.pdf", width = 100/25.6, height = 60/25.6, family = "Arial")
p_heat_com_d7
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

setwd("~cchen2/Documents/neuroscience/Pn\ project/Figure/PDF/")
cairo_pdf("p_trace_d7.pdf", width = 130/25.6, height = 65/25.6, family = "Arial")
p_trace
dev.off()

dat_cell_area <- c()
cell_area_day <- c("Pre",  "Test")
for (i in c(1:2)){
  rep_time <- c(ncol(dat_stim_m3_d7[[i]]), ncol(dat_stim_m7_d7[[i]]), ncol(dat_stim_m17_d7[[i]]), ncol(dat_stim_m18_d7[[i]]))
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

dat_cell_area$Day <- factor(dat_cell_area$Day, levels = c("Pre","Test"))
dat_cell_area_com <- ddply(dat_cell_area, .(ID, Day), summarise, value=mean(ratio, na.rm = T))
dat_cell_area_sta <- ddply(dat_cell_area_com,.(Day), summarise,n=length(value),mean=mean(value),sd=sd(value),se=sd(value)/sqrt(length(value)))


p_EI_ratio <- ggplot(dat_cell_area_sta, aes(Day, mean, colour=Day, fill=Day))+
  geom_bar(stat="identity", color="black", position=position_dodge())+
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.1, position=position_dodge(.9), colour="black") +
  geom_jitter(data = dat_cell_area_com, aes(Day, value),width = 0.1, colour="black", shape=1 )+
  # scale_fill_manual(values=c( "mediumseagreen","indianred4", "dodgerblue4"))+
  labs(x="", y="E/I ratio")+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  scale_y_continuous(limits = c(0, 2), expand = c(0,0))+
  theme(legend.title = element_blank())+
  theme(legend.position = "none")

setwd("~cchen2/Documents/neuroscience/Pn\ project/Figure/PDF/")
cairo_pdf("p_EI_d7_comp.pdf", width = 42/25.6, height = 65/25.6, family = "Arial")
p_EI_ratio
dev.off()
