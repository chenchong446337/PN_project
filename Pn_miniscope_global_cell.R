## This script is used to analyze cell activity change before and after learning anticipation of pain relief
## average the trace of identified cells from Day1-3, and Day5-6, since not all cells are identified by Matlab
## created on 02202020
## last updated 03062020


## library
library("ggplot2")
library("reshape2")
library("plyr")
library(openxlsx)
library(cowplot)
library(grid)
library(DescTools)
library(zoo)

## make a function for global ID analysis-----
cc_globalID_fun <- function(path_ID, path_trace, t_stim){
  Global_ID <- read.csv(path_ID, header = F)
  ## split into two groups (day1-3) and Day5-6, choose the one both have at leat one no-0 element
  
  Global_ID$sum1_3 <- rowSums(Global_ID[,2:3])
  Global_ID$sum5_6 <- rowSums(Global_ID[, 5:6])
  
  dat_global_ID <- subset(Global_ID, sum1_3>0 & sum5_6>0)
  dat_global_ID<- dat_global_ID[,c(1:7)]
  dat_global_ID[dat_global_ID==0]<- NA
  
  ## read.xlsx with no header and remove the first col
  cc_read.xlsx<- function(x){
    dat <- read.xlsx(x, colNames = F)
    dat<- dat[,-1]
    dat
  }
  dat_trace <- mapply(cc_read.xlsx, path_trace,SIMPLIFY = F)
  stim_time<- seq(-5, 8.5, by=0.5)
  ## number of rows to be binned
  n <- 10 # 0.05*10=0.5
  
  ## activity of each cell during each day
  dat_cell_trace_day <- vector(mode = "list", length = nrow(dat_global_ID))
  
  for (i in c(1:nrow(dat_global_ID))){
    cell_trace_day <- matrix(0, nrow = length(stim_time), ncol = 7)
    for (j in c(1:7)){
      cell_id <- dat_global_ID[i,j]
      if (is.na(cell_id)) {
        dat_stim_mean <- rep(NA, length(stim_time))
      } else{
        cell_trace <- dat_trace[[j]][cell_id,]
        cell_stim <- t_stim[[j]]
        dat_stim <- matrix(0, nrow = length(stim_time), ncol = length(cell_stim))
        
        for(z in seq_along(cell_stim)){
          t1_p <- cell_stim[z]
          dat_stim1 <- unname(unlist(c(cell_trace[(t1_p-100):(t1_p+180-1)])))
          dat_stim1<- rollapply(dat_stim1, n, mean, by = n, align = "left", partial = TRUE)
          dat_stim1_base <- unname(unlist(c(cell_trace[(t1_p-100):(t1_p-60-1)])))
          dat_stim1_base_mean <- mean(dat_stim1_base)
          dat_stim1_base_sd <- sd(dat_stim1_base)
          dat_stim1_nor <- (dat_stim1 - dat_stim1_base_mean)/dat_stim1_base_sd
          dat_stim[,z] <- dat_stim1_nor
        }
        
        dat_stim_mean <- rowMeans(dat_stim)
        
      }
      cell_trace_day[,j]<- dat_stim_mean
    }
    dat_cell_trace_day[[i]]<- cell_trace_day
  }
  ## calculate the mean activity of each cells for day2-3 and day5-6
  cell_active_learning <- mapply(function(x) rowMeans(x[,2:3], na.rm=TRUE), dat_cell_trace_day)
  cell_active_anti <- mapply(function(x) rowMeans(x[,5:6], na.rm=TRUE), dat_cell_trace_day)
  cell_active_anti_pain <- mapply(function(x) x[,7], dat_cell_trace_day)
  ## statistical test to see which cells activation increased
  t_compare <- which(stim_time>=0)
  t_cells_increase<- sapply(1:ncol(cell_active_learning),function(x) wilcox.test(cell_active_anti[t_compare,x],cell_active_learning[t_compare,x], paired = T, alternative = "greater")$p.value)
  t_cells_decrease<- sapply(1:ncol(cell_active_learning),function(x) wilcox.test(cell_active_anti[t_compare,x],cell_active_learning[t_compare,x], paired = T, alternative = "less")$p.value)
  
  ## compare the cells number of activated or decreased
  ID_cell_increase <- which(t_cells_increase<0.001)
  ID_cell_decrease <- which(t_cells_decrease< 0.001)
  ID_cell_unchange <- setdiff(c(1:ncol(cell_active_anti)), c(ID_cell_increase, ID_cell_decrease))
  n_cell_increase <- length(which(t_cells_increase<0.001))
  n_cell_decrease <- length(which(t_cells_decrease<0.001))
  
  cell_increase_learning_mean <- rowMeans(cell_active_learning[,ID_cell_increase])
  cell_increase_anti_mean <- rowMeans(cell_active_anti[, ID_cell_increase])
  cell_increase_anti_pain_mean <- rowMeans(cell_active_anti_pain[,ID_cell_increase], na.rm = T)
  cell_increase_act <- c(cell_increase_learning_mean, cell_increase_anti_mean, cell_increase_anti_pain_mean)
  cell_increase_time <- rep(stim_time, 3)
  cell_increase_group <- rep(c("Ctrl", "Anti","Pain"), each=length(cell_increase_anti_mean))
  dat_cell_increase <- data_frame(Group=cell_increase_group, Time= cell_increase_time, Value = cell_increase_act)
  cell_increase_acu <- mapply(function(x) AUC(stim_time[t_compare],x[t_compare]) , list(cell_increase_learning_mean, cell_increase_anti_mean, cell_increase_anti_pain_mean))

  ## compare the cells which show increase activity during learning, anti and pain
  cell_decrease_learning_mean <- rowMeans(cell_active_learning[,ID_cell_decrease])
  cell_decrease_anti_mean <- rowMeans(cell_active_anti[, ID_cell_decrease])
  cell_decrease_anti_pain_mean <- rowMeans(cell_active_anti_pain[,ID_cell_decrease], na.rm = T)
  cell_decrease_act <- c(cell_decrease_learning_mean, cell_decrease_anti_mean, cell_decrease_anti_pain_mean)
  cell_decrease_time <- rep(stim_time, 3)
  cell_decrease_group <- rep(c("Ctrl", "Anti","Pain"), each=length(cell_decrease_anti_mean))
  dat_cell_decrease <- data_frame(Group=cell_decrease_group, Time= cell_decrease_time, Value = cell_decrease_act)
  cell_decrease_acu <- mapply(function(x) AUC(stim_time[t_compare],x[t_compare]) , list(cell_decrease_learning_mean, cell_decrease_anti_mean, cell_decrease_anti_pain_mean))
  
  ## compare the cells which show un-changed activity during learning, anti and pain
  cell_unchange_learning_mean <- rowMeans(cell_active_learning[,ID_cell_unchange])
  cell_unchange_anti_mean <- rowMeans(cell_active_anti[, ID_cell_unchange])
  cell_unchange_anti_pain_mean <- rowMeans(cell_active_anti_pain[,ID_cell_unchange], na.rm = T)
  cell_unchange_act <- c(cell_unchange_learning_mean, cell_unchange_anti_mean, cell_unchange_anti_pain_mean)
  cell_unchange_time <- rep(stim_time, 3)
  cell_unchange_group <- rep(c("Ctrl", "Anti","Pain"), each=length(cell_unchange_anti_mean))
  dat_cell_unchange <- data_frame(Group=cell_unchange_group, Time= cell_unchange_time, Value = cell_unchange_act)
  cell_unchange_acu <- mapply(function(x) AUC(stim_time[t_compare],x[t_compare]) , list(cell_unchange_learning_mean, cell_unchange_anti_mean, cell_unchange_anti_pain_mean))
  
  return(list(dat_cell_trace_day, dat_cell_increase, dat_cell_decrease, dat_cell_unchange))
}


## 1. Test the cell activity across sessions------
## for mice m3
Global_ID_m3 <- "~cchen2/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/matlab_analysis/m3/Global_id/m3_global_ID.csv"
path_trace_m3 <- list.files(path ="~cchen2/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/matlab_analysis/m3/trace_days/", pattern = "*.xlsx" )
setwd("~cchen2/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/matlab_analysis/m3/trace_days")
t_stim_m3_d1 <- c(214, 698,1204,1757,2964) 
t_stim_m3_d2 <- c(274, 831,1364,1788,2621,3273)
t_stim_m3_d3 <- c(1885, 2516) 
t_stim_m3_d4 <- c(1813,3587)
t_stim_m3_d5 <- c(758) 
t_stim_m3_d6 <- c(132,507) 
t_stim_m3_d7 <- c(346, 1206) ## use first two 346, 1206, 1756, 2288
t_stim_m3 <- list(t_stim_m3_d1, t_stim_m3_d2, t_stim_m3_d3, t_stim_m3_d4, t_stim_m3_d5, t_stim_m3_d6, t_stim_m3_d7)

dat_global_ID_m3 <- cc_globalID_fun(Global_ID_m3, path_trace_m3, t_stim_m3)

## for mice m7
Global_ID_m7 <- "~cchen2/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/matlab_analysis/m7/Global_ID/m7_global_ID.csv"
path_trace_m7 <- list.files(path ="~cchen2/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/matlab_analysis/m7/trace_days/", pattern = "*.xlsx" )
setwd("~cchen2/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/matlab_analysis/m7/trace_days")
t_stim_m7_d1 <- c( 605, 964, 1603, 2153, 2363, 2751, 3152) 
t_stim_m7_d2 <- c(798, 177,3300)
t_stim_m7_d3 <- c(360, 3352, 4044) 
t_stim_m7_d4 <- c(3236)
t_stim_m7_d5 <- c(317) 
t_stim_m7_d6 <- c(205)
t_stim_m7_d7 <- c(236,1914, 2777, 3334) ## use first two 236, 1914,2777,3334,4091
t_stim_m7 <- list(t_stim_m7_d1, t_stim_m7_d2, t_stim_m7_d3, t_stim_m7_d4, t_stim_m7_d5, t_stim_m7_d6, t_stim_m7_d7)


dat_global_ID_m7 <- cc_globalID_fun(Global_ID_m7, path_trace_m7, t_stim_m7)


## combine the data and make a plot
## for cells increase activation-----
dat_increase_com <- rbind(as.data.frame(dat_global_ID_m3[[2]]), as.data.frame(dat_global_ID_m7[[2]]))
dat_increase_com$m_ID <- rep(c("m3", "m7"), each=nrow(as.data.frame(dat_global_ID_m3[[2]])))
dat_increase_com_sta <- ddply(dat_increase_com, .(Time, Group), summarise,n=length(Value),mean=mean(Value),sd=sd(Value),se=sd(Value)/sqrt(length(Value)))

p_cell_increase <- ggplot(dat_increase_com_sta, aes(Time, mean, group=Group, colour=Group))+
  geom_line()+
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), alpha=0.3, colour="gray")+
  labs(x="Time relative to crossing (s)", y="Z score")+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  geom_vline(xintercept = 0, col= 'red', linetype=2)+
  geom_hline(yintercept = 0, col="grey", linetype=2)+
  theme(legend.title = element_blank())


## for cells with decreased activation---
dat_decrease_com <- rbind(as.data.frame(dat_global_ID_m3[[3]]), as.data.frame(dat_global_ID_m7[[3]]))
dat_decrease_com$m_ID <- rep(c("m3", "m7"), each=nrow(as.data.frame(dat_global_ID_m3[[2]])))
dat_decrease_com_sta <- ddply(dat_decrease_com, .(Time, Group), summarise,n=length(Value),mean=mean(Value),sd=sd(Value),se=sd(Value)/sqrt(length(Value)))


p_cell_decrease <- ggplot(dat_decrease_com_sta, aes(Time, mean, group=Group, colour=Group))+
  geom_line()+
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), alpha=0.3, colour="gray")+
  labs(x="Time relative to crossing (s)", y="Z score")+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  geom_vline(xintercept = 0, col= 'red', linetype=2)+
  geom_hline(yintercept = 0, col="grey", linetype=2)+
  theme(legend.title = element_blank())

## for cells with unchanged activation----
dat_stable_com <- rbind(as.data.frame(dat_global_ID_m3[[4]]), as.data.frame(dat_global_ID_m7[[4]]))
dat_stable_com$m_ID <- rep(c("m3", "m7"), each=nrow(as.data.frame(dat_global_ID_m3[[2]])))
dat_stable_com_sta <- ddply(dat_stable_com, .(Time, Group), summarise,n=length(Value),mean=mean(Value),sd=sd(Value),se=sd(Value)/sqrt(length(Value)))

p_cell_stable <- ggplot(dat_stable_com_sta, aes(Time, mean, group=Group, colour=Group))+
  geom_line()+
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), alpha=0.3, colour="gray")+
  labs(x="Time relative to crossing (s)", y="Z score")+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  geom_vline(xintercept = 0, col= 'red', linetype=2)+
  geom_hline(yintercept = 0, col="grey", linetype=2)+
  theme(legend.title = element_blank())

## combine data of all kinds of cells and plot with same y axis
dat_com <- rbind(dat_increase_com_sta, dat_decrease_com_sta, dat_stable_com_sta)
dat_com$cell_type <- rep(c("Increased", "Decreased", "Stable"), each=nrow(dat_decrease_com_sta))
dat_com$cell_type<- factor(dat_com$cell_type, c("Increased", "Decreased", "Stable"))
p_com <- ggplot(dat_com, aes(Time, mean, group=Group, colour=Group))+
  geom_line()+
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), alpha=0.3, colour="gray")+
  labs(x="Time relative to crossing (s)", y="Z score")+
  facet_wrap(.~cell_type)+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  geom_vline(xintercept = 0, col= 'red', linetype=2)+
  geom_hline(yintercept = 0, col="grey", linetype=2)+
  theme(legend.title = element_blank())

## save the figure
setwd("~cchen2/Documents/neuroscience/Pn\ project/Figure/PDF/")
cairo_pdf("p_com.pdf", width = 260/25.6, height = 80/25.6, family = "Arial")
p_com
dev.off()



## correlation test between activity increased and decreased cells during diff situation-----
dat_cell_cor_m3 <- cbind(as.data.frame(dat_global_ID_m3[c(2,3)]))
dat_cell_cor_m3 <- dat_cell_cor_m3[,c(1:3, 6)]
colnames(dat_cell_cor_m3)[3:4] <- c("Value_in", "Value_de")
dat_cell_cor_m3$Group<- factor(dat_cell_cor_m3$Group, c("Ctrl", "Anti", "Pain"))

p_cor<- ggplot(dat_cell_cor_m3, aes(Value_de, Value_in,))+
  geom_point(shape=1)+
  geom_smooth(method='lm',colour="black")+
  labs(x="Z score of activation decreased cells", y="Z score of activation increased cells")+
  facet_wrap(.~Group, scales = "free")+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  theme(legend.position = "none")
  

dat_cell_cor_p<- by(dat_cell_cor_m3, dat_cell_cor$Group, FUN = function(X) cor.test(X$Value_de, X$Value_in, method = "spearman")$p.value)
dat_cell_cor_r <- by(dat_cell_cor_m3, dat_cell_cor$Group, FUN = function(X) cor(X$Value_de, X$Value_in, method = "spearman"))

## save the figure
setwd("~cchen2/Documents/neuroscience/Pn\ project/Figure/PDF/")
cairo_pdf("p_cor.pdf", width = 230/25.6, height = 80/25.6, family = "Arial")
p_cor
dev.off()

## plot the trace of increased activity cells and decreased, compared the time of peaks
dat_cell_com_comb <- rbind(dat_increase_com, dat_decrease_com)
dat_cell_com_comb$Activation <- rep(c("Increase", "Decrease"), each=nrow(dat_increase_com))

ggplot(dat_cell_com_comb, aes(Time,Value, colour=Activation ))+
  geom_line()+
  facet_wrap(.~m_ID+Group, scales = "free")

## 2. Compare cell activation during anticipation of pain relief and pain-----
cc_globalID_fun2 <- function(path_ID, path_trace, t_stim){
  Global_ID <- read.csv(path_ID, header = F)
  ## split into two groups (day1-3) and Day5-6, choose the one both have at leat one no-0 element
  
  Global_ID$sum5_6 <- rowSums(Global_ID[, 5:6])
  
  dat_global_ID <- subset(Global_ID, sum5_6>0 & Global_ID[,8]>0)
  dat_global_ID<- dat_global_ID[,c(5, 6, 8)]
  dat_global_ID[dat_global_ID==0]<- NA
  
  ## read.xlsx with no header and remove the first col
  cc_read.xlsx<- function(x){
    dat <- read.xlsx(x, colNames = F)
    dat<- dat[,-1]
    dat
  }
  dat_trace <- mapply(cc_read.xlsx, path_trace,SIMPLIFY = F)
  stim_time<- seq(-5, 8.5, by=0.5)
  ## number of rows to be binned
  n <- 10 # 0.05*10=0.5
  
  ## activity of each cell during anti (day5-6) and pinprick
  dat_cell_trace_day <- vector(mode = "list", length = nrow(dat_global_ID))
  
  for (i in c(1:nrow(dat_global_ID))){
    cell_trace_day <- matrix(0, nrow = length(stim_time), ncol = length(dat_trace))
    for (j in c(1:length(dat_trace))){
      cell_id <- dat_global_ID[i,j]
      if (is.na(cell_id)) {
        dat_stim_mean <- rep(NA, length(stim_time))
      } else{
        cell_trace <- dat_trace[[j]][cell_id,]
        cell_stim <- t_stim[[j]]
        dat_stim <- matrix(0, nrow = length(stim_time), ncol = length(cell_stim))
        
        for(z in seq_along(cell_stim)){
          t1_p <- cell_stim[z]
          dat_stim1 <- as.numeric(cell_trace[(t1_p-100):(t1_p+180-1)])
          dat_stim1<- rollapply(dat_stim1, n, mean, by = n, align = "left", partial = TRUE)
          dat_stim1_base <- as.numeric(cell_trace[(t1_p-100):(t1_p-60-1)])
          dat_stim1_base_mean <- mean(dat_stim1_base)
          dat_stim1_base_sd <- sd(dat_stim1_base)
          dat_stim1_nor <- (dat_stim1 - dat_stim1_base_mean)/dat_stim1_base_sd
          dat_stim[,z] <- dat_stim1_nor
        }
        
        dat_stim_mean <- rowMeans(dat_stim)
        
      }
      cell_trace_day[,j]<- dat_stim_mean
    }
    dat_cell_trace_day[[i]]<- cell_trace_day
  }
  ## get the ID of cells which has increased and decreased activation during anti pain relief
  cell_active_anti <- mapply(function(x) rowMeans(x[,1:2], na.rm=TRUE), dat_cell_trace_day)
  cell_active_pin <- mapply(function(x) identity(x[,3]), dat_cell_trace_day)
  
  ## statistical test to see which cells activation increased or decreased
  row_trim1 <- which(stim_time <= -3)
  row_trim2 <- which(stim_time>0 & stim_time<=5)
  cell_active_anti_p_increase <- unname(apply(cell_active_anti, 2, function(x) wilcox.test(x[row_trim2], x[row_trim1],alternative = "greater")$p.value))
  cell_active_anti_p_decrease <- unname(apply(cell_active_anti, 2, function(x) wilcox.test(x[row_trim2], x[row_trim1],alternative = "less")$p.value))
  
  cell_active_pin_p_increase <- unname(apply(cell_active_pin, 2, function(x) wilcox.test(x[row_trim2], x[row_trim1],alternative = "greater")$p.value))
  cell_active_pin_p_decrease <- unname(apply(cell_active_pin, 2, function(x) wilcox.test(x[row_trim2], x[row_trim1],alternative = "less")$p.value))
  
  ## compare the cells number of activated or decreased
  ID_cell_increase_anti <- which(cell_active_anti_p_increase<0.001)
  ID_cell_decrease_anti <- which(cell_active_anti_p_decrease< 0.001)
  ID_cell_unchange_anti <- setdiff(c(1:ncol(cell_active_anti)), c(ID_cell_increase_anti, ID_cell_decrease_anti))
  
  ID_cell_increase_pin <- which(cell_active_pin_p_increase<0.001)
  ID_cell_decrease_pin <- which(cell_active_pin_p_decrease< 0.001)
  ID_cell_unchange_pin <- setdiff(c(1:ncol(cell_active_pin)), c(ID_cell_increase_pin, ID_cell_decrease_pin))
  
  cell_anti_increase_anti_mean <- rowMeans(cell_active_anti[, ID_cell_increase_anti])
  cell_anti_increase_pin_mean <- rowMeans(cell_active_pin[,ID_cell_increase_anti], na.rm = T)
  cell_anti_increase_act <- c(cell_anti_increase_anti_mean, cell_anti_increase_pin_mean)
  
  cell_pin_increase_anti_mean <- rowMeans(cell_active_anti[, ID_cell_increase_pin])
  cell_pin_increase_pin_mean <- rowMeans(cell_active_pin[,ID_cell_increase_pin], na.rm = T)
  cell_pin_increase_act <- c(cell_pin_increase_anti_mean, cell_pin_increase_pin_mean)
  
  
  cell_increase_time <- rep(stim_time, 2)
  cell_increase_group <- rep(c("Anti","Pin"), each=length(cell_anti_increase_anti_mean))
  
  ## compare the cells which show decrease activity during learning, anti and pain
  cell_anti_decrease_anti_mean <- rowMeans(cell_active_anti[, ID_cell_decrease_anti])
  cell_anti_decrease_pin_mean <- rowMeans(cell_active_pin[,ID_cell_decrease_anti], na.rm = T)
  cell_anti_decrease_act <- c(cell_anti_decrease_anti_mean, cell_anti_decrease_pin_mean)
  
  cell_pin_decrease_anti_mean <- rowMeans(cell_active_anti[, ID_cell_decrease_pin])
  cell_pin_decrease_pin_mean <- rowMeans(cell_active_pin[,ID_cell_decrease_pin], na.rm = T)
  cell_pin_decrease_act <- c(cell_pin_decrease_anti_mean, cell_pin_decrease_pin_mean)
  
  ## compare the cells which show un-changed activity during learning, anti and pain
  cell_anti_unchange_anti_mean <- rowMeans(cell_active_anti[, ID_cell_unchange_anti])
  cell_anti_unchange_pin_mean <- rowMeans(cell_active_pin[,ID_cell_unchange_anti], na.rm = T)
  cell_anti_unchange_act <- c(cell_anti_unchange_anti_mean, cell_anti_unchange_pin_mean)
  
  cell_pin_unchange_anti_mean <- rowMeans(cell_active_anti[, ID_cell_unchange_pin])
  cell_pin_unchange_pin_mean <- rowMeans(cell_active_pin[,ID_cell_unchange_pin], na.rm = T)
  cell_pin_unchange_act <- c(cell_pin_unchange_anti_mean, cell_pin_unchange_pin_mean)
  
  
  dat_anti_pin <- data.frame(Time = cell_increase_time, Group = cell_increase_group, anti_in = cell_anti_increase_act,
                             anti_de = cell_anti_decrease_act, anti_un = cell_anti_unchange_act)
  
  dat_pin_anti <- data.frame(Time = cell_increase_time, Group = cell_increase_group, pin_in = cell_pin_increase_act,
                             pin_de = cell_pin_decrease_act, pin_un = cell_pin_unchange_act)
  
  return(list(dat_anti_pin, dat_pin_anti))
}

## for m3
Global_ID_m3 <- "~cchen2/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/matlab_analysis/m3/Global_id/m3_anti_pin_global_ID.csv"
path_trace_m3 <- list.files(path ="~cchen2/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/matlab_analysis/m3/trace_days/", pattern = "*.xlsx" )[5:6]
path_trace_m3_pin <- "~cchen2/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/matlab_analysis/m3/pin_prick/m3_trace_pinprick.xlsx"
path_trace_m3_comp <- c(path_trace_m3, path_trace_m3_pin)
setwd("~cchen2/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/matlab_analysis/m3/trace_days")

t_stim_m3_d5 <- c(758) 
t_stim_m3_d6 <- c(132,507) 
t_pin_m3 <- c(337,665,916,1303,1541,1767,2032,2277,2496,2824)

t_stim_m3 <- list( t_stim_m3_d5, t_stim_m3_d6, t_pin_m3)

dat_anti_pin_m3 <- cc_globalID_fun2(Global_ID_m3, path_trace_m3_comp, t_stim_m3)

## For mice m7
Global_ID_m7 <- "~cchen2/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/matlab_analysis/m7/Global_ID/m7_anti_pin_global_ID.csv"
path_trace_m7 <- list.files(path ="~cchen2/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/matlab_analysis/m7/trace_days/", pattern = "*.xlsx" )[5:6]
path_trace_m7_pin <- "~cchen2/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/matlab_analysis/m7/pin_prick/m7_trace_pinprick.xlsx"
path_trace_m7_comp <- c(path_trace_m7, path_trace_m7_pin)
setwd("~cchen2/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/matlab_analysis/m7/trace_days")

t_stim_m7_d5 <- c(317) 
t_stim_m7_d6 <- c(205)
t_pin_m7 <- c(493, 1323, 1867, 2415, 3035, 3556, 4079, 4519, 5020, 5517) 

t_stim_m7 <- list(t_stim_m7_d5, t_stim_m7_d6, t_pin_m7)

dat_anti_pin_m7 <- cc_globalID_fun2(Global_ID_m7, path_trace_m7_comp, t_stim_m7)



## combine all data and plot the trace to compare the cell activation during anti and pinprick
dat_pin_anti <- rbind(dat_anti_pin_m3[[1]], dat_anti_pin_m7[[1]])
dat_pin_anti$m_ID <- rep(c("m3", "m7"), each=nrow(dat_anti_pin_m3[[1]]))

dat_pin_anti_re <- melt(dat_pin_anti, id.vars = c("Group", "Time", "m_ID"))
dat_pin_anti_sta<- ddply(dat_pin_anti_re, .(Time, Group, variable), summarise,n=length(value),mean=mean(value),sd=sd(value),se=sd(value)/sqrt(length(value)))

p_anti_in <- ggplot(subset(dat_pin_anti_sta, dat_pin_anti_sta$variable=="anti_in"), aes(Time, mean, colour=Group))+
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), alpha=0.3, colour="gray")+
  geom_line()+
  labs(x="Time relative to crossing (s)", y="Z score")+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  geom_vline(xintercept = 0, col= 'red', linetype=2)+
  geom_hline(yintercept = 0, col="grey", linetype=2)+
  theme(legend.title = element_blank())

p_anti_de <- ggplot(subset(dat_pin_anti_sta, dat_pin_anti_sta$variable=="anti_de"), aes(Time, mean, colour=Group))+
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), alpha=0.3, colour="gray")+
  geom_line()+
  labs(x="Time relative to crossing (s)", y="Z score")+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  geom_vline(xintercept = 0, col= 'red', linetype=2)+
  geom_hline(yintercept = 0, col="grey", linetype=2)+
  theme(legend.title = element_blank())

p_anti_un <- ggplot(subset(dat_pin_anti_sta, dat_pin_anti_sta$variable=="anti_un"), aes(Time, mean, colour=Group))+
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), alpha=0.3, colour="gray")+
  geom_line()+
  labs(x="Time relative to crossing (s)", y="Z score")+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  geom_vline(xintercept = 0, col= 'red', linetype=2)+
  geom_hline(yintercept = 0, col="grey", linetype=2)+
  theme(legend.title = element_blank())

p_anti_pin_com<- plot_grid(p_anti_in, p_anti_de, p_anti_un, nrow = 1)
## for pin activated cells
dat_pin_anti <- rbind(dat_anti_pin_m3[[2]], dat_anti_pin_m7[[2]])
dat_pin_anti$m_ID <- rep(c("m3", "m7"), each=nrow(dat_anti_pin_m3[[2]]))

dat_pin_anti_re <- melt(dat_pin_anti, id.vars = c("Group", "Time", "m_ID"))
dat_pin_anti_sta<- ddply(dat_pin_anti_re, .(Time, Group, variable), summarise,n=length(value),mean=mean(value),sd=sd(value),se=sd(value)/sqrt(length(value)))

p_pin_in <- ggplot(subset(dat_pin_anti_sta, dat_pin_anti_sta$variable=="pin_in"), aes(Time, mean, colour=Group))+
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), alpha=0.3, colour="gray")+
  geom_line()+
  labs(x="Time relative to crossing (s)", y="Z score")+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  geom_vline(xintercept = 0, col= 'red', linetype=2)+
  geom_hline(yintercept = 0, col="grey", linetype=2)+
  theme(legend.title = element_blank())

p_pin_de <- ggplot(subset(dat_pin_anti_sta, dat_pin_anti_sta$variable=="pin_de"), aes(Time, mean, colour=Group))+
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), alpha=0.3, colour="gray")+
  geom_line()+
  labs(x="Time relative to crossing (s)", y="Z score")+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  geom_vline(xintercept = 0, col= 'red', linetype=2)+
  geom_hline(yintercept = 0, col="grey", linetype=2)+
  theme(legend.title = element_blank())

p_pin_un <- ggplot(subset(dat_pin_anti_sta, dat_pin_anti_sta$variable=="pin_un"), aes(Time, mean, colour=Group))+
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), alpha=0.3, colour="gray")+
  geom_line()+
  labs(x="Time relative to crossing (s)", y="Z score")+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  geom_vline(xintercept = 0, col= 'red', linetype=2)+
  geom_hline(yintercept = 0, col="grey", linetype=2)+
  theme(legend.title = element_blank())

p_pin_anti_com<- plot_grid(p_pin_in, p_pin_de, p_pin_un, nrow = 1)

p_sub_anti_pin <- plot_grid(p_anti_pin_com, p_pin_anti_com, ncol = 1)

setwd("~cchen2/Documents/neuroscience/Pn\ project/Figure/PDF/")
cairo_pdf("p_sub_anti_pin.pdf", width = 300/25.6, height = 150/25.6, family = "Arial")
p_sub_anti_pin
dev.off()


## golbal ID alignment, only for the fist stimulation-------
cc_globalID_fun <- function(path_ID, path_trace, t_stim){
  dat_global_ID <- read.csv(path_ID, header = F) %>% 
    mutate(., sum1_3 = rowSums(.[, 2:3]), sum5_6 = rowSums(.[, 5:6])) %>% 
    filter(., sum1_3 > 0 & sum5_6 >0 & V7 >0) %>% 
    select(., c(1:7)) %>% 
    na_if(., 0)
  
  
  ## read.xlsx with no header and remove the first col, and do the z-score
  cc_read.xlsx<- function(x){
    dat <- read.xlsx(x, colNames = F) %>% 
      select(., -c(1,2)) %>% 
      t() %>% 
      as.data.frame() %>% 
      apply(., 2, scale)
    return(dat)
  }
  dat_trace <- mapply(cc_read.xlsx, path_trace,SIMPLIFY = F)
  stim_time<- seq(-2, 6.5, by=0.5)
  ## number of rows to be binned
  n <- 10 # 0.05*10=0.5
  
  ## activity of each cell during each day
  dat_cell_trace_day <- vector(mode = "list", length = nrow(dat_global_ID))
  
  for (i in c(1:nrow(dat_global_ID))){
    cell_trace_day <- matrix(0, nrow = length(stim_time), ncol = 7)
    
    for (j in c(1:7)){
      cell_id <- dat_global_ID[i,j]
      if (is.na(cell_id)) {
        dat_stim_mean <- rep(NA, length(stim_time))
      } else{
        cell_trace <- dat_trace[[j]][,cell_id]
        t1_p <- t_stim[[j]]
        
        dat_stim1 <- cell_trace[(t1_p-40):(t1_p+140-1)] %>% 
          matrix(., n) %>% 
          colMeans() 
        # dat_stim1 <- aggregate(dat_stim1,list(rep(1:(nrow(dat_stim1)%/%n+1),each=n,len=nrow(dat_stim1))),mean)[-1]
        
        dat_stim_mean <- dat_stim1 - mean(dat_stim1[1:5]) ## baseline as -2 to 0
        
      }
      cell_trace_day[,j]<- dat_stim_mean
    }
    dat_cell_trace_day[[i]]<- cell_trace_day
  }
  return(dat_cell_trace_day)
}

## analyze for each mcie------
Global_ID_m3 <- "~cchen2/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/matlab_analysis/m3/Global_id/m3_global_ID.csv"
path_trace_m3 <- list.files(path ="~cchen2/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/matlab_analysis/m3/trace_days/", pattern = "*.xlsx", full.names = T )
t_stim_m3_d1 <- c(214) 
t_stim_m3_d2 <- c(274)
t_stim_m3_d3 <- c(1885) 
t_stim_m3_d4 <- c(1813)
t_stim_m3_d5 <- c(758) 
t_stim_m3_d6 <- c(132) 
t_stim_m3_d7 <- c(346) ## use first two 346, 1206, 1756, 2288, 2970, 3262, 3942
t_stim_m3 <- list(t_stim_m3_d1, t_stim_m3_d2, t_stim_m3_d3, t_stim_m3_d4, t_stim_m3_d5, t_stim_m3_d6, t_stim_m3_d7)

dat_m3_global <- cc_globalID_fun(Global_ID_m3, path_trace_m3, t_stim_m3)

## for m7
Global_ID_m7 <- "~cchen2/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/matlab_analysis/m7/Global_ID/m7_global_ID.csv"
path_trace_m7 <- list.files(path ="~cchen2/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/matlab_analysis/m7/trace_days/", pattern = "*.xlsx", full.names = T )
t_stim_m7_d1 <- c(605) 
t_stim_m7_d2 <- c(798)
t_stim_m7_d3 <- c(360) 
t_stim_m7_d4 <- c(3236)
t_stim_m7_d5 <- c(317) 
t_stim_m7_d6 <- c(205)
t_stim_m7_d7 <- c(236) ## use first two 236, 1914,2777,3334,4091
t_stim_m7 <- list(t_stim_m7_d1, t_stim_m7_d2, t_stim_m7_d3, t_stim_m7_d4, t_stim_m7_d5, t_stim_m7_d6, t_stim_m7_d7)

dat_m7_global <- cc_globalID_fun(Global_ID_m7, path_trace_m7, t_stim_m7)

## for m17
Global_ID_m17 <- "~cchen2/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/matlab_analysis/m17/Global_ID/m17_global_ID.csv"
path_trace_m17 <- list.files(path ="~cchen2/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/matlab_analysis/m17/trace_days/", pattern = "*.xlsx", full.names = T )
t_stim_m17_d1 <- c(815) 
t_stim_m17_d2 <- c(1103)
t_stim_m17_d3 <- c(437) 
t_stim_m17_d4 <- c(682)
t_stim_m17_d5 <- c(1114) 
t_stim_m17_d6 <- c(617) 
t_stim_m17_d7 <- c(157) # 157, 2536, 2961
t_stim_m17 <- list(t_stim_m17_d1, t_stim_m17_d2, t_stim_m17_d3, t_stim_m17_d4, t_stim_m17_d5, t_stim_m17_d6, t_stim_m17_d7)

dat_m17_global <- cc_globalID_fun(Global_ID_m17, path_trace_m17, t_stim_m17)

## for m18
Global_ID_m18 <- "~cchen2/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/matlab_analysis/m18/Global_ID/m18_global_ID.csv"
path_trace_m18 <- list.files(path ="~cchen2/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/matlab_analysis/m18/trace_days/", pattern = "*.xlsx", full.names = T )
t_stim_m18_d1 <- c(490) 
t_stim_m18_d2 <- c(236)
t_stim_m18_d3 <- c(784) 
t_stim_m18_d4 <- c(493)
t_stim_m18_d5 <- c(195) 
t_stim_m18_d6 <- c(467) 
t_stim_m18_d7 <- c(124) ## use the first two of 124,2238,3000
t_stim_m18 <- list(t_stim_m18_d1, t_stim_m18_d2, t_stim_m18_d3, t_stim_m18_d4, t_stim_m18_d5, t_stim_m18_d6, t_stim_m18_d7)

dat_m18_global <- cc_globalID_fun(Global_ID_m18, path_trace_m18, t_stim_m18)

stim_time<- seq(-2, 6.5, by=0.5)
mouse_ID <- c("m3", "m7", "m17", "m18")
con_day <- c("Rm", "Pre","Pre","Rm","Cond.","Cond.", "Test")
dat_global_cell_trace <- NULL
for (i in c(3, 6,7)) {
  trace_com <- cbind(mapply(function (x) x[,i], dat_m3_global), mapply(function (x) x[,i], dat_m7_global), mapply(function (x) x[,i], dat_m17_global), mapply(function (x) x[,i], dat_m18_global)) %>% 
    as.data.frame() 
  rep_col <- colnames(trace_com)[colSums(is.na(trace_com)) > 0]
  
  if (length(rep_col) >0) {
    trace_com_rep <- cbind(mapply(function (x) x[,i-1], dat_m3_global), mapply(function (x) x[,i-1], dat_m7_global), mapply(function (x) x[,i-1], dat_m17_global), mapply(function (x) x[,i-1], dat_m18_global)) %>% 
      as.data.frame()
    trace_com[, rep_col] <- trace_com_rep[, rep_col]
  } else {
    trace_com <- trace_com
  }
  
  trace_com <- set_colnames(trace_com, sprintf("Cell%s",seq(1:ncol(trace_com)))) %>% 
    mutate(., Time= stim_time) %>% 
    melt(., id.vars="Time") %>% 
    mutate(., ID = rep(mouse_ID, c(length(dat_m3_global), length(dat_m7_global), length(dat_m17_global), length(dat_m18_global))* length(stim_time))) %>% 
    mutate(., Day = con_day[i])
  dat_global_cell_trace <- rbind(dat_global_cell_trace, trace_com)
}

## align by the activatio of cells
dat_trace <- filter(dat_global_cell_trace, Day=="Pre")
dat_trace_sta <- ddply(dat_trace, .(variable), summarise,value=mean(value), sum=sum(value)) %>% 
  arrange(., value)

dat_trace$variable <- factor(dat_trace$variable, levels = dat_trace_sta$variable)

dat_global_cell_trace <- mutate(dat_global_cell_trace,variable= factor(variable, levels = dat_trace_sta$variable),
                                Day = factor(Day, levels = c("Pre", "Cond.", "Test")))
score_range <- range(dat_global_cell_trace$value)

p_heat <- ggplot(dat_global_cell_trace, aes(Time, variable,fill= value))+ 
  geom_tile(height=2)+
  facet_grid(cols = vars(Day), scales = "free_y")+
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
