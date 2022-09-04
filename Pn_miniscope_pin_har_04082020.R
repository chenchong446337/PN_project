## update log
## 05052020 calculate firing freq of ctrl and pin
## 05102020 for expectation of pain comparing


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

## function for data with reduced frame

## for pinprick, ctrl and expectation of pain-----
path_m3_pin <- "~cchen2/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/matlab_analysis/m3/pin_prick/m3_trace_pinprick.xlsx"
t_pin_m3 <- c(337,665,916,1303,1541,1767,2032,2277,2496,2824)
# m3 ctrl
path_m3_pin_ctrl <- "~cchen2/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/matlab_analysis/m3/pin_prick/m3_trace_pin_ctrl.xlsx"
t_pin_ctrl_m3 <- c(134, 528, 820, 1084, 1261, 1493, 1729, 1964, 2142,2367)*2

path_pin_m3 <- list(path_m3_pin_ctrl, path_m3_pin)
stim_pin_m3 <- list(t_pin_ctrl_m3,t_pin_m3 )

dat_pin_m3 <-mapply(c_miniscope_matlab, path_pin_m3, stim_pin_m3, SIMPLIFY = F)
## for m7
path_m7_pin <- "~cchen2/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/matlab_analysis/m7/pin_prick/m7_trace_pinprick.xlsx"
t_pin_m7 <- c(493, 1323, 1867, 2415, 3035, 3556, 4079, 4519, 5020, 5517) 
## m7 ctrl
path_m7_pin_ctrl <- "~cchen2/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/matlab_analysis/m7/pin_prick/m7_trace_pin_ctrl.xlsx"
t_pin_ctrl_m7 <- c(195,594, 901, 1163, 1536, 1818, 2137, 2450, 2624, 3008 )*2

path_pin_m7 <- list(path_m7_pin_ctrl, path_m7_pin)
stim_pin_m7 <- list(t_pin_ctrl_m7, t_pin_m7)

dat_pin_m7 <-mapply(c_miniscope_matlab, path_pin_m7, stim_pin_m7, SIMPLIFY = F)


## for m16
path_m16_pin <- "~cchen2/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/matlab_analysis/m16/pin_prick/m16_trace_pin.xlsx"
t_pin_m16 <- c(115, 414, 641,965, 1272, 1597, 1879, 2130, 2438,2679)*2

## m16 ctrl
path_m16_pin_ctrl <- "~cchen2/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/matlab_analysis/m16/pin_prick/m16_trace_pin_ctrl.xlsx"
t_pin_ctrl_m16 <- c(192, 498, 779, 1068, 1359, 1624, 1918, 2209, 2592, 2858)*2 

path_pin_m16 <- list(path_m16_pin_ctrl, path_m16_pin)
stim_pin_m16 <- list(t_pin_ctrl_m16, t_pin_m16)

dat_pin_m16 <-mapply(c_miniscope_matlab, path_pin_m16, stim_pin_m16, SIMPLIFY = F)
## for m17
path_m17_pin <- "~cchen2/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/matlab_analysis/m17/pin_prick/m17_trace_pin.xlsx"
t_pin_m17 <- c(262, 497, 803, 1057, 1404, 1676, 2020, 2248, 2543, 2805 )*2

## m17 ctrl
path_m17_pin_ctrl <- "~cchen2/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/matlab_analysis/m17/pin_prick/m17_trace_pin_ctrl.xlsx"
t_pin_ctrl_m17 <- c(216, 605, 982, 1198, 1604, 1829, 2098, 2348, 2520, 2743)*2

path_pin_m17 <- list(path_m17_pin_ctrl, path_m17_pin )
stim_pin_m17 <- list(t_pin_ctrl_m17, t_pin_m17)

dat_pin_m17 <-mapply(c_miniscope_matlab, path_pin_m17, stim_pin_m17, SIMPLIFY = F)
## for m18
path_m18_pin <- "~cchen2/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/matlab_analysis/m18/pin_prick/m18_trace_pin.xlsx"
t_pin_m18 <- c(139, 681, 1012, 1365, 1663, 1900, 2110, 2324, 2566, 2769)*2
## m18 ctrl
path_m18_pin_ctrl <- "~cchen2/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/matlab_analysis/m18/pin_prick/m18_trace_pin_ctrl.xlsx"
t_pin_ctrl_m18 <- c(225, 465, 883, 1385, 1561, 1735, 1897, 2085, 2370, 2691 ) *2

path_pin_m18 <- list(path_m18_pin_ctrl, path_m18_pin)
stim_pin_m18 <- list(t_pin_ctrl_m18, t_pin_m18)

dat_pin_m18 <-mapply(c_miniscope_matlab, path_pin_m18, stim_pin_m18, SIMPLIFY = F)

## for the exp pain
path_exp_pain <- as.list(list.files("~cchen2/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/pinprick_pain_exp", pattern = str_c("^m.*.xlsx"), full.names = T ))
dat_stim_type <- read.csv("~cchen2/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/pinprick_pain_exp/random_pinprick.csv", row.names = 1)

t_exp_pain <- read.xlsx("~cchen2/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/pinprick_pain_exp/pin_event.xlsx") %>% 
  .[dat_stim_type =="Pain"] %>% 
  matrix(., nrow = 10) %>%
  as.data.frame() %>% 
  '*' (2) %>% 
  as.list()

t_exp_ctrl <- read.xlsx("~cchen2/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/pinprick_pain_exp/pin_event.xlsx") %>% 
  .[dat_stim_type!="Pain"] %>% 
  matrix(., nrow = 10) %>%
  as.data.frame() %>% 
  '*' (2) %>% 
  as.list()

dat_exp_ctrl1 <- mapply(c_miniscope_matlab, path_exp_pain, t_exp_ctrl, SIMPLIFY = F)
dat_exp_ctrl <- do.call(cbind,dat_exp_ctrl1)
  
dat_exp_pain <- mapply(c_miniscope_matlab, path_exp_pain, t_exp_pain, SIMPLIFY = F) %>% 
  do.call(cbind,.)

## for hargveast experiment
path_har <- as.list(list.files("~cchen2/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/har_06112020/", pattern = str_c("^m.*.xlsx"), full.names = T ))

t_stim_har <- read.xlsx("~cchen2/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/har_06112020/2020_06_01_har_events.xlsx") 
t_stim_har <- split(t_stim_har$Frame,t_stim_har$Mice)

dat_har <- mapply(c_miniscope_matlab, path_har, t_stim_har)
dat_cell_har <- do.call(cbind,dat_har)

## combine data and do k-means analysis-----

dat_cell_pin <- vector(mode = "list", 5)

for (i in c(1,2)) {
  dat_cell_pin[[i]] <- cbind(dat_pin_m3[[i]], dat_pin_m7[[i]], dat_pin_m16[[i]],dat_pin_m17[[i]], dat_pin_m18[[i]] )
  
}

dat_cell_pin[[4]] <- dat_exp_ctrl
dat_cell_pin[[5]] <- dat_exp_pain
dat_cell_pin[[3]] <- dat_cell_har

dat_cell_pin_re <- vector(mode = "list", 5)

# k-means clustering
mouse_ID <- c("m3", "m7", "m16","m17", "m18")
for (i in 1:length(dat_cell_pin)) {
  dat_cell_pin_d <- dat_cell_pin[[i]]
  cell_id <- sprintf("Cell%s",seq(1:ncol(dat_cell_pin_d)))
  colnames(dat_cell_pin_d)<- cell_id
  
  dat_cell_pin_d_cluster <- unname(kmeans(t(dat_cell_pin_d), centers = 3, nstart = 20)[[1]])
  stim_time<- seq(-2, 6.5, by=0.5)
  #stim_time<- seq(-5, 8.5, by=0.5)
  
  dat_cell_pin_d_re<- as.data.frame(dat_cell_pin_d) %>% 
    mutate(., Time= stim_time) %>% 
    melt(., id.vars='Time') %>% 
    mutate(., Group= rep(dat_cell_pin_d_cluster, each=length(stim_time)) )
  
  ## sort data by the value in each group
  dat_cell_d_sort <- as.numeric(names(sort(tapply(dat_cell_pin_d_re$value, dat_cell_pin_d_re$Group, mean), decreasing = T)))
  
  dat_cell_pin_d_re$Group[dat_cell_pin_d_re$Group == dat_cell_d_sort[1]] ="Excited"
  dat_cell_pin_d_re$Group[dat_cell_pin_d_re$Group == dat_cell_d_sort[2]] ="Neutral"
  dat_cell_pin_d_re$Group[dat_cell_pin_d_re$Group == dat_cell_d_sort[3]] ="Inhibited"
  
  if (i <3){
    rep_time <- c(ncol(dat_pin_m3[[i]]), ncol(dat_pin_m7[[i]]), ncol(dat_pin_m16[[i]]),ncol(dat_pin_m17[[i]]), ncol(dat_pin_m18[[i]]))
  } else if (i>3) {
    rep_time <- mapply(ncol, dat_exp_ctrl1)
  } else {
    rep_time <- mapply(ncol, dat_har)
  }
  dat_cell_pin_d_re <- mutate(dat_cell_pin_d_re, Group= factor(Group, levels = c("Excited", "Neutral", "Inhibited"))) %>% 
    mutate(., ID = rep(mouse_ID, rep_time*length(stim_time)) )
  
  dat_cell_pin_re[[i]] <- dat_cell_pin_d_re
}

score_range <- range(mapply(function(x) x$value, dat_cell_pin_re, SIMPLIFY = T))
pin_group <- c("Touch", "Pin","Har", "Shuffle_touch", "Shuffle_pin")

for (i in c(1:5)) {
  dat_pin <- dat_cell_pin_re[[i]]
  dat_pin_sta <- ddply(dat_pin, .(variable, Group), summarise,mean=mean(value), sum=sum(value))
  dat_pin_sta <- dat_pin_sta[order(dat_pin_sta[,'mean']),]
  dat_pin$variable <- factor(dat_pin$variable, levels = dat_pin_sta$variable)
  p_heat <- ggplot(dat_pin, aes(Time, variable,fill= value))+ 
    geom_raster() +
    #facet_grid(rows = vars(Group), scales = "free_y")+
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
  assign(paste0("p_heat_", pin_group[i]), p_heat)
}

## combine the heat plot
p_heat_com <- plot_grid(p_heat_Touch, p_heat_Pin, p_heat_Har, p_heat_Shuffle_touch,p_heat_Shuffle_pin, nrow = 1)

setwd("~cchen2/Documents/neuroscience/Pn\ project/Figure/PDF/")
cairo_pdf("p_pin_heat.pdf", width = 180/25.6, height = 60/25.6, family = "Arial")
p_heat_com
dev.off()

## plot trace by group

dat_cell_pin_sta <- mapply(function(x) ddply(x, .(ID,Time, Group), summarise, value=mean(value)), dat_cell_pin_re, SIMPLIFY = F) %>% 
  mapply(function (x) ddply(x, .(Time, Group), summarise, n=length(value),mean=mean(value),sd=sd(value),se=sd(value)/sqrt(length(value))),., SIMPLIFY = F) %>% 
  do.call(rbind,.) %>% 
  mutate(., Day = rep(c("Touch", "Pin","Har", "Shuffle_touch", "Shuffle_pin"),each= length(stim_time)*3)) %>% 
  mutate(Day=factor(Day, levels = c("Touch", "Pin","Har", "Shuffle_touch", "Shuffle_pin"))) %>% 
  mutate(Group = factor(Group, levels = c("Neutral", "Excited", "Inhibited")))

range_pin_plot <- range(dat_cell_pin_sta$mean)

## group by day
p_pin <- ggplot(dat_cell_pin_sta, aes(Time, mean, colour=Group))+
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
cairo_pdf("p_pin.pdf", width = 160/25.6, height = 60/25.6, family = "Arial")
p_pin
dev.off()

## plot E and I points to show correlation
## points to show the correlation
dat_cell_cor <- NULL

for (i in c(1:5)) {
  dat_pin <- dat_cell_pin_re[[i]] %>% 
    filter(ID=="m18") %>% 
    filter(Group!="Neutral") %>% 
    ddply(., .(Time, Group), summarise, value=mean(value)) %>%
    dcast(., Time~Group) %>%
    mutate(., Day = pin_group[i])
  dat_cell_cor <- rbind(dat_cell_cor, dat_pin)
  
}

dat_cell_cor <- dat_cell_cor %>% 
  mutate(Day= factor(Day, levels = c("Touch", "Pin","Har", "Shuffle_touch", "Shuffle_pin")))

p_trace_cor <- ggplot(dat_cell_cor, aes(Excited, Inhibited, colour=Day))+
  geom_point()+
  geom_smooth(method = "lm", lwd=0.8)+
  scale_colour_manual(values=c("seagreen", "indianred", "deepskyblue4", "deeppink3", "lightseagreen"))+  
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
cairo_pdf("p_trace_cor.pdf", width = 82/25.6, height = 60/25.6, family = "Arial")
p_trace_cor
dev.off()

## correlation for each mice
dat_cell_pin_cor <- NULL
for (i in c(1:5)){
  dat_pin_cor <- dat_cell_pin_re[[i]] %>%
    ddply(., .(ID, Time, Group), summarise, value=mean(value)) %>%
    subset(., .$Group!="Neutral") %>%
    dcast(., Time + ID ~ Group) %>%
    ddply(., .(ID), summarise, "corr" = cor(Excited, Inhibited, method = "spearman")) %>%
    mutate(., Day = pin_group[i])
  dat_cell_pin_cor <- rbind(dat_cell_pin_cor, dat_pin_cor)
}

dat_cor_value <- ddply(dat_cell_pin_cor, .(ID, Day), summarise, "corr" = mean(corr)) %>%
  mutate(Day = factor(Day, levels = c("Touch", "Pin","Har", "Shuffle_touch", "Shuffle_pin")))

p_cor <-ggplot(dat_cor_value, aes(Day, corr, fill=Day))+
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(width = 0.2, shape=1)+
  scale_fill_manual(values=c("seagreen", "indianred","deepskyblue4", "deeppink3", "lightseagreen"))+
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
cairo_pdf("p_cor.pdf", width = 40/25.6, height = 60/25.6, family = "Arial")
p_cor
dev.off()

## calculate the sum of active and inhibited trace
dat_cell_pin_sum <-  mapply(function(x) ddply(x, .(ID,Time, Group), summarise, value=mean(value)), dat_cell_pin_re, SIMPLIFY = F) %>%
  do.call(rbind,.) %>% 
  subset(., .$Group!="Neutral") %>%
  mutate( Day = rep(c("Touch", "Pin","Har", "Shuffle_touch", "Shuffle_pin"), each= length(stim_time)*2*length(mouse_ID))) %>%
  ddply(., .(ID,Time, Day), summarise, value=sum(value)) %>%
  ddply(., .(Time, Day), summarise,n=length(value),mean=mean(value),sd=sd(value),se=sd(value)/sqrt(length(value))) %>%
  mutate(Day = factor(Day, levels = c("Touch", "Pin","Har", "Shuffle_touch", "Shuffle_pin")))

p_pin_sum <- ggplot(dat_cell_pin_sum, aes(Time, mean, colour=Day))+
  geom_line()+
  geom_ribbon(aes(ymin=mean-se, ymax=mean+se, fill=Day), alpha=0.1, linetype=0)+
  scale_colour_manual(values=c("seagreen", "indianred", "deepskyblue4","deeppink3","lightseagreen"))+
  scale_fill_manual(values=c("seagreen", "indianred", "deepskyblue4","deeppink3","lightseagreen"))+
  labs(x="Time relative to crossing (s)", y="AUC (z score)")+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  geom_vline(xintercept = 0, col= 'red', linetype=2)+
  theme(legend.title = element_blank(), legend.position = "none")

setwd("~cchen2/Documents/neuroscience/Pn\ project/Figure/PDF/")
cairo_pdf("p_pin_sum.pdf", width = 90/25.6, height = 60/25.6, family = "Arial")
p_pin_sum
dev.off()
## for EI change analysis
dat_cell_area <- c()
for (i in c(1:5)){
  dat_pin_sta <- ddply(dat_cell_pin_re[[i]], .(ID, Time, Group), summarise,n=length(value),mean=mean(value),sd=sd(value),se=sd(value)/sqrt(length(value)))
  Time_area <- which(stim_time>=0)
  dat_anti_area <- as.data.frame(tapply(dat_pin_sta$mean, INDEX = list(dat_pin_sta$ID, dat_pin_sta$Group), 
                                        function (x) AUC(stim_time[Time_area], x[Time_area]))) %>% 
    mutate(ID = rownames(.)) %>% 
    mutate(Day=pin_group[i]) %>% 
    mutate(ratio = abs(Excited/Inhibited), sum = Excited + Inhibited)
  
  rownames(dat_anti_area)<-NULL
  
  dat_cell_area <- rbind(dat_cell_area, dat_anti_area)
}

dat_cell_area_sta <-  ddply(dat_cell_area, .(ID, Day), summarise, value=mean(ratio, na.rm = T)) %>% 
  mutate(Day = factor(Day, levels = c("Touch", "Pin","Har", "Shuffle_touch", "Shuffle_pin")))

p_EI_ratio <- ggplot(dat_cell_area_sta, aes(Day, value, fill=Day))+
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(width = 0.2, shape=1)+
  scale_fill_manual(values=c("seagreen", "indianred",  "deepskyblue4", "deeppink3", "lightseagreen"))+
  labs(x="", y="E/I ratio")+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  scale_y_continuous(limits = c(0, 3), expand = c(0,0))+
  theme(legend.title = element_blank(), legend.position = "none")

setwd("~cchen2/Documents/neuroscience/Pn\ project/Figure/PDF/")
cairo_pdf("p_EI_ratio.pdf", width = 40/25.6, height = 60/25.6, family = "Arial")
p_EI_ratio
dev.off()

## plot the sum of E-I
dat_cell_sum_sta <-  ddply(dat_cell_area, .(ID, Day), summarise, "sum"= mean(sum, na.rm = T)) %>% 
  mutate(Day = factor(Day, levels = c("Touch", "Pin","Har", "Shuffle_touch", "Shuffle_pin")))


p_EI_sum <- ggplot(dat_cell_sum_sta, aes(Day, sum, fill=Day))+
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(width = 0.2, shape=1)+
  scale_fill_manual(values=c("seagreen", "indianred", "deepskyblue4", "deeppink3", "lightseagreen"))+
  labs(x="", y="AUG (z-score)")+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        axis.text.x = element_text(angle = 40, vjust=1, hjust=1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  #scale_y_continuous(limits = c(0, 10), expand = c(0,0))+
  theme(legend.title = element_blank(), legend.position = "none")

test_pin_sum <- pairwise.wilcox.test(dat_cell_sum_sta$sum, dat_cell_sum_sta$Day)
t_ei_sum <-aov(sum~Day, data=dat_cell_sum_sta)
summary(t_ei_sum)
TukeyHSD(t_ei_sum, which = "Day")

setwd("~cchen2/Documents/neuroscience/Pn\ project/Figure/PDF/")
cairo_pdf("p_pin_EI_sum.pdf", width = 60/25.6, height = 75/25.6, family = "Arial")
p_EI_sum
dev.off()

## show the change of E, I and sum
p_EIS_change <- dat_cell_area %>%
  select("Excited", "Inhibited","ID","Day", "sum") %>%
  melt(., id.vars= c('Day', "ID")) %>%
  ddply(., .(ID, Day, variable), summarise, value=mean(value, na.rm = T)) %>%
  mutate(Day = factor(Day, levels = c("Touch", "Pin","Har", "Shuffle_touch", "Shuffle_pin"))) %>% 
  ggplot(., aes(variable, value, fill=Day))+
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(position=position_jitterdodge(), shape=1, size=1)+
  scale_fill_manual(values=c("seagreen", "indianred", "deepskyblue4", "deeppink3", "lightseagreen"))+
  labs(x="", y="AUC (z-score)")+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  theme(legend.position = "none")
  #theme(legend.title = element_blank(), legend.position = "top")

setwd("~cchen2/Documents/neuroscience/Pn\ project/Figure/PDF/")
cairo_pdf("p_EIS_change.pdf", width = 82/25.6, height = 60/25.6, family = "Arial")
p_EIS_change
dev.off()

## calculate the firing frequency before and after pinprick-----
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

## for m3
path_m3_pin <- "~cchen2/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/matlab_analysis/m3/pin_prick/m3_trace_pinprick.xlsx"
t_pin_m3 <- c(337,665,916,1303,1541,1767,2032,2277,2496,2824)
path_peak_m3 <- "~cchen2/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/matlab_analysis/m3/pin_prick/m3_pin_peaks.csv"

path_m7_pin <- "~cchen2/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/matlab_analysis/m7/pin_prick/m7_trace_pinprick.xlsx"
t_pin_m7 <- c(493, 1323, 1867, 2415, 3035, 3556, 4079, 4519, 5020, 5517) 
path_peak_m7 <- "~cchen2/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/matlab_analysis/m7/pin_prick/m7_pin_peaks.csv"

## combine data to analyze
path_ID_list <- list(path_m3_pin, path_m7_pin)
t_stim_list <- list(t_pin_m3, t_pin_m7)
path_peak_list <- list(path_peak_m3, path_peak_m7)

dat_firng <- mapply(cc_pin_peak, path_ID_list, t_stim_list, path_peak_list, SIMPLIFY = F) %>%
  do.call(rbind, .)

## for heat pain------
path_m3_har <- "~cchen2/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/matlab_analysis/Hargreave/m3_trace_har.xlsx"
t_har_m3 <- c(1171, 2159, 3579, 4708, 5693, 6910)

dat_har_m3 <-c_miniscope_matlab(path_m3_har, t_har_m3)

## for m7
path_m7_har <- "~cchen2/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/matlab_analysis/Hargreave/m7_trace_har.xlsx"
t_har_m7 <- c(756, 1917, 2947, 3541, 4402, 5450, 6937)

dat_har_m7 <-c_miniscope_matlab(path_m7_har, t_har_m7)

## combine har data
dat_cell_trace_har <- cbind(dat_har_m3, dat_har_m7)
cell_id <- sprintf("Cell%s",seq(1:ncol(dat_cell_trace_har)))
colnames(dat_cell_trace_har)<- cell_id

dat_cell_trace_har_cluster <- unname(kmeans(t(dat_cell_trace_har), centers = 3)[[1]])
dat_cell_trace_har<- as.data.frame(dat_cell_trace_har)

stim_time<- seq(-2, 6.5, by=0.5)
dat_cell_trace_har$Time <- stim_time

dat_cell_trace_har_re <- melt(dat_cell_trace_har, id.vars ='Time')
dat_cell_trace_har_re$Group <- rep(dat_cell_trace_har_cluster, each=length(stim_time))

## sort data by the value in each group
dat_cell_har_sort <- as.numeric(names(sort(tapply(dat_cell_trace_har_re$value, dat_cell_trace_har_re$Group, sum), decreasing = T)))

dat_cell_trace_har_re$Group[dat_cell_trace_har_re$Group == dat_cell_har_sort[1]] ="Excited"
dat_cell_trace_har_re$Group[dat_cell_trace_har_re$Group == dat_cell_har_sort[2]] ="Neutral"
dat_cell_trace_har_re$Group[dat_cell_trace_har_re$Group == dat_cell_har_sort[3]] ="Inhibited"

dat_cell_trace_har_re$Group <- factor(dat_cell_trace_har_re$Group, levels = c("Excited", "Neutral", "Inhibited"))


p_heat_har <- ggplot(dat_cell_trace_har_re, aes(Time, variable,fill= value))+ 
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


## plot the har trace by groups
dat_cell_trace_har_re_sta <- ddply(dat_cell_trace_har_re, .(Time, Group), summarise,n=length(value),mean=mean(value),sd=sd(value),se=sd(value)/sqrt(length(value)))
dat_cell_trace_har_re_sta$Group <- factor(dat_cell_trace_har_re_sta$Group, levels = c("Neutral", "Excited", "Inhibited"))
p_trace_har <- ggplot(dat_cell_trace_har_re_sta, aes(Time, mean, colour=Group))+
  geom_line()+
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), alpha=0.3, colour="gray")+
  scale_colour_manual(values=c("#999999", "darkred", "navy"))+
  labs(x="Time relative to crossing (s)", y="Z score")+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  #scale_y_continuous(limits = c(range_trace_plot[1]- 0.2, range_trace_plot[2]+0.2))+
  geom_vline(xintercept = 0, col= 'red', linetype=2)+
  theme(legend.position = 'none')

## for har trace
Time_area <- which(stim_time>=0)

## area under curve for the whole population
dat_har_area <- tapply(dat_cell_trace_har_re_sta$mean, INDEX = dat_cell_trace_har_re_sta$Group, 
                       function (x) AUC(stim_time[Time_area], x[Time_area]), simplify = T)
#dat_har_area$ratio <- abs(dat_anti_area$Excited/ dat_anti_area$Inhibited)
har_EI_ratio <- abs(unname(dat_har_area[2]/dat_har_area[3]))
dat_har_EI_ratio <- data.frame(Group="har", ratio = har_EI_ratio)


setwd("~cchen2/Documents/neuroscience/Pn\ project/Figure/PDF/")
cairo_pdf("p_har_trace.pdf", width = 65/25.6, height = 60/25.6, family = "Arial")
p_trace_har
dev.off()

setwd("~cchen2/Documents/neuroscience/Pn\ project/Figure/PDF/")
cairo_pdf("p_har_heat.pdf", width = 85/25.6, height = 60/25.6, family = "Arial")
p_heat_har
dev.off()

## combine har and pin ratio for plot
dat_pin_har <- rbind(dat_pin_EI_ratio, dat_har_EI_ratio)
p_EI_ratio_pin <- ggplot(dat_pin_har, aes(Group, ratio, fill=Group))+
  geom_bar(stat="identity")+
  labs(x="", y="E/I ratio")+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  scale_y_continuous(expand = c(0,0), limits = c(0,4))+
  theme(legend.position = 'none')

setwd("~cchen2/Documents/neuroscience/Pn\ project/Figure/PDF/")
cairo_pdf("p_EI_pin.pdf", width = 32/25.6, height = 60/25.6, family = "Arial")
p_EI_ratio_pin
dev.off()


## compare pain and expectation of pain-------
path_exp_pain <- as.list(list.files("~cchen2/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/pinprick_pain_exp", pattern = str_c("^m.*.xlsx"), full.names = T ))
dat_stim_type <- read.csv("~cchen2/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/pinprick_pain_exp/random_pinprick.csv", row.names = 1)

t_exp_pain <- read.xlsx("~cchen2/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/pinprick_pain_exp/pin_event.xlsx") %>% 
  .[dat_stim_type=="Pain"] %>% 
  matrix(., nrow = 10) %>%
  as.data.frame() %>% 
  as.list()

t_exp_ctrl <- read.xlsx("~cchen2/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/pinprick_pain_exp/pin_event.xlsx") %>% 
  .[dat_stim_type!="Pain"] %>% 
  matrix(., nrow = 10) %>%
  as.data.frame() %>% 
  as.list()

dat_exp_pin <- mapply(c_miniscope_matlab, path_exp_pain, t_exp_pain, SIMPLIFY = F)
dat_exp_ctrl <- mapply(c_miniscope_matlab, path_exp_pain, t_exp_ctrl, SIMPLIFY = F)

## combine data and do k-means analysis
dat_cell_exp_pin <- list(do.call(cbind, dat_exp_pin), do.call(cbind, dat_exp_ctrl))

dat_cell_exp_pin_re <- vector(mode = "list", 2)

# k-means clustering
mouse_ID <- c("m3", "m7", "m16","m17", "m18")
for (i in 1:length(dat_cell_exp_pin)) {
  dat_cell_exp_pin_d <- dat_cell_exp_pin[[i]]
  cell_id <- sprintf("Cell%s",seq(1:ncol(dat_cell_exp_pin_d)))
  colnames(dat_cell_exp_pin_d)<- cell_id
  
  dat_cell_exp_pin_d_cluster <- unname(kmeans(t(dat_cell_exp_pin_d), centers = 3, nstart = 20)[[1]])
  stim_time<- seq(-2, 6.5, by=0.5)
  #stim_time<- seq(-5, 8.5, by=0.5)
  
  dat_cell_exp_pin_d_re<- as.data.frame(dat_cell_exp_pin_d) %>% 
    mutate(., Time= stim_time) %>% 
    melt(., id.vars='Time') %>% 
    mutate(., Group= rep(dat_cell_exp_pin_d_cluster, each=length(stim_time)) )
  
  ## sort data by the value in each group
  dat_cell_d_sort <- as.numeric(names(sort(tapply(dat_cell_exp_pin_d_re$value, dat_cell_exp_pin_d_re$Group, mean), decreasing = T)))
  
  dat_cell_exp_pin_d_re$Group[dat_cell_exp_pin_d_re$Group == dat_cell_d_sort[1]] ="Excited"
  dat_cell_exp_pin_d_re$Group[dat_cell_exp_pin_d_re$Group == dat_cell_d_sort[2]] ="Neutral"
  dat_cell_exp_pin_d_re$Group[dat_cell_exp_pin_d_re$Group == dat_cell_d_sort[3]] ="Inhibited"
  
  if (i ==1){
    rep_time <- mapply(ncol, dat_exp_pin)
  } else {
    rep_time <- mapply(ncol, dat_exp_ctrl)
  }
  dat_cell_exp_pin_d_re <- mutate(dat_cell_exp_pin_d_re, Group= factor(Group, levels = c("Excited", "Neutral", "Inhibited"))) %>% 
    mutate(., ID = rep(mouse_ID, rep_time*length(stim_time)) )
  
  dat_cell_exp_pin_re[[i]] <- dat_cell_exp_pin_d_re
}

score_range <- range(mapply(function(x) x$value, dat_cell_exp_pin_re, SIMPLIFY = T))
exp_pin_group <- c("pin", "exp_pin")
for (i in c(1,2)) {
  dat_exp_pin <- dat_cell_exp_pin_re[[i]]
  dat_exp_pin_sta <- ddply(dat_exp_pin, .(variable, Group), summarise,mean=mean(value), sum=sum(value))
  dat_exp_pin_sta <- dat_exp_pin_sta[order(dat_exp_pin_sta[,'mean']),]
  dat_exp_pin$variable <- factor(dat_exp_pin$variable, levels = dat_exp_pin_sta$variable)
  p_heat <- ggplot(dat_exp_pin, aes(Time, variable,fill= value))+ 
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
  assign(paste0("p_heat_", exp_pin_group[i]), p_heat)
}

## combine the heat plot
p_heat_com <- plot_grid(p_heat_pin, p_heat_exp_pin, nrow = 1)

## plot trace by group
dat_cell_exp_pin_test_sta <- ddply(dat_cell_exp_pin_re[[1]], .(ID,Time, Group), summarise,value=mean(value))

dat_cell_exp_pin_ctrl_test_sta <- ddply(dat_cell_exp_pin_re[[2]], .(ID,Time, Group), summarise,value=mean(value))

dat_cell_exp_pin_sta <- ddply(dat_cell_exp_pin_test_sta, .(Time, Group), summarise,n=length(value),mean=mean(value),sd=sd(value),se=sd(value)/sqrt(length(value))) %>%
  rbind(., ddply(dat_cell_exp_pin_ctrl_test_sta, .(Time, Group), summarise,n=length(value),mean=mean(value),sd=sd(value),se=sd(value)/sqrt(length(value)))) %>%
  mutate(., Day = rep(c("pin", "exp_pin"),each= length(stim_time)*3))

dat_cell_exp_pin_sta$Day <- factor(dat_cell_exp_pin_sta$Day, levels = c("pin", "exp_pin"))
dat_cell_exp_pin_sta$Group <- factor(dat_cell_exp_pin_sta$Group, levels = c("Neutral", "Excited", "Inhibited"))
range_exp_pin_plot <- range(dat_cell_exp_pin_sta$mean)

## group by day
p_exp_pin <- ggplot(dat_cell_exp_pin_sta, aes(Time, mean, colour=Group))+
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


## plot E and I points to show correlation-----
## points to show the correlation
dat_cell_cor <- NULL

for (i in c(1,2)) {
  dat_exp_pin <- subset(dat_cell_exp_pin_re[[i]], dat_cell_exp_pin_re[[i]]$ID=="m3") %>%
    subset(., .$Group != "Neutral") %>%
    ddply(., .(Time, Group), summarise, value=mean(value)) %>%
    dcast(., Time~Group) %>%
    mutate(., Day = exp_pin_group[i])
  dat_cell_cor <- rbind(dat_cell_cor, dat_exp_pin)
  
}
dat_cell_cor$Day <- factor(dat_cell_cor$Day, levels = c("pin", "exp_pin"))

p_trace_cor <- ggplot(dat_cell_cor, aes(Excited, Inhibited, colour=Day))+
  geom_point()+
  geom_smooth(method = "lm", lwd=0.8)+
  scale_colour_manual(values=c("seagreen", "indianred"))+  
  labs(x="Excited (z-score)", y="Inhibited (z-score)")+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  theme(legend.title = element_blank(), legend.position = c(0.8, 0.8))


## correlation for each mice
dat_cell_exp_pin_cor <- NULL
for (i in c(1,2)){
  dat_exp_pin_cor <- dat_cell_exp_pin_re[[i]] %>%
    ddply(., .(ID, Time, Group), summarise, value=mean(value)) %>%
    subset(., .$Group!="Neutral") %>%
    dcast(., Time + ID ~ Group) %>%
    ddply(., .(ID), summarise, "corr" = cor(Excited, Inhibited, method = "spearman")) %>%
    mutate(., Day = exp_pin_group[i])
  dat_cell_exp_pin_cor <- rbind(dat_cell_exp_pin_cor, dat_exp_pin_cor)
}

dat_cor_value <- ddply(dat_cell_exp_pin_cor, .(ID, Day), summarise, "corr" = mean(corr)) %>%
  mutate(Day = factor(Day, levels = c("pin", "exp_pin")))

p_cor <-ggplot(dat_cor_value, aes(Day, corr, fill=Day))+
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(width = 0.2, shape=1)+
  scale_fill_manual(values=c("seagreen", "indianred"))+
  labs(x="", y="Correlation coefficient (r)")+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  theme(legend.title = element_blank(), legend.position = "none")


## calculate the sum of active and inhibited trace
dat_cell_exp_pin_sum <- rbind(dat_cell_exp_pin_test_sta, dat_cell_exp_pin_ctrl_test_sta) %>%
  subset(., .$Group!="Neutral") %>%
  mutate( Day = rep(c("pin", "exp_pin"), each= length(stim_time)*2*length(mouse_ID))) %>%
  ddply(., .(ID,Time, Day), summarise, value=sum(value)) %>%
  ddply(., .(Time, Day), summarise,n=length(value),mean=mean(value),sd=sd(value),se=sd(value)/sqrt(length(value))) %>%
  mutate(Day = factor(Day, levels = c("pin", "exp_pin")))

p_exp_pin_sum <- ggplot(dat_cell_exp_pin_sum, aes(Time, mean, colour=Day))+
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

## for EI change analysis
dat_cell_area <- c()
for (i in c(1,2)){
  rep_time <- mapply(ncol, dat_exp_ctrl)
  dat_exp_pin <- dat_cell_exp_pin_re[[i]]
  dat_exp_pin$ID <- rep(mouse_ID, rep_time*length(stim_time))
  dat_exp_pin_sta <- ddply(dat_exp_pin, .(ID, Time, Group), summarise,n=length(value),mean=mean(value),sd=sd(value),se=sd(value)/sqrt(length(value)))
  Time_area <- which(stim_time>=0)
  dat_anti_area <- as.data.frame(tapply(dat_exp_pin_sta$mean, INDEX = list(dat_exp_pin_sta$ID, dat_exp_pin_sta$Group), 
                                        function (x) AUC(stim_time[Time_area], x[Time_area])))
  dat_anti_area$ID <- rownames(dat_anti_area)
  rownames(dat_anti_area)<-NULL
  dat_anti_area$Day <- exp_pin_group[i]
  dat_anti_area$ratio <- abs(dat_anti_area$Excited/dat_anti_area$Inhibited)
  dat_anti_area$sum <- dat_anti_area$Excited + dat_anti_area$Inhibited
  dat_cell_area <- rbind(dat_cell_area, dat_anti_area)
}

dat_cell_area$Day <- factor(dat_cell_area$Day, levels = c("pin", "exp_pin"))
dat_cell_area_sta <-  ddply(dat_cell_area, .(ID, Day), summarise, value=mean(ratio, na.rm = T))

p_EI_ratio <- ggplot(dat_cell_area_sta, aes(Day, value, fill=Day))+
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(width = 0.2, shape=1)+
  scale_fill_manual(values=c("seagreen", "indianred"))+
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

## plot the sum of E-I
dat_cell_sum_sta <-  ddply(dat_cell_area, .(ID, Day), summarise, "sum"=mean(sum, na.rm = T))


p_EI_sum<- ggplot(dat_cell_sum_sta, aes(Day, sum, fill=Day))+
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(width = 0.2, shape=1)+
  scale_fill_manual(values=c("seagreen", "indianred"))+
  labs(x="", y="AUG (z-score)")+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  #scale_y_continuous(limits = c(0, 10), expand = c(0,0))+
  theme(legend.title = element_blank(), legend.position = "none")

## show the change of E, I and sum
p_EIS_change <- dat_cell_area %>%
  select("Excited", "Inhibited","ID","Day", "sum") %>%
  melt(., id.vars= c('Day', "ID")) %>%
  ddply(., .(ID, Day, variable), summarise, value=mean(value, na.rm = T)) %>%
  ggplot(., aes(variable, value, fill=Day))+
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(position=position_jitterdodge(), shape=1, size=1)+
  scale_fill_manual(values=c("seagreen", "indianred"))+
  labs(x="", y="AUC (z-score)")+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  theme(legend.title = element_blank(), legend.position = "top")

## compare ctrl and exp-ctrl-----
dat_cell_exp_ctrl_re <- list(dat_cell_pin_re[[2]], dat_cell_exp_pin_re[[2]])

score_range <- range(mapply(function(x) x$value, dat_cell_exp_ctrl_re, SIMPLIFY = T))
exp_ctrl_group <- c("ctrl", "exp_ctrl")
for (i in c(1,2)) {
  dat_exp_ctrl <- dat_cell_exp_ctrl_re[[i]]
  dat_exp_ctrl_sta <- ddply(dat_exp_ctrl, .(variable, Group), summarise,mean=mean(value), sum=sum(value))
  dat_exp_ctrl_sta <- dat_exp_ctrl_sta[order(dat_exp_ctrl_sta[,'mean']),]
  dat_exp_ctrl$variable <- factor(dat_exp_ctrl$variable, levels = dat_exp_ctrl_sta$variable)
  p_heat <- ggplot(dat_exp_ctrl, aes(Time, variable,fill= value))+ 
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
  assign(paste0("p_heat_", exp_ctrl_group[i]), p_heat)
}

## combine the heat plot
p_heat_com <- plot_grid(p_heat_ctrl, p_heat_exp_ctrl, nrow = 1)

## plot trace by group
dat_cell_exp_ctrl_test_sta <- ddply(dat_cell_exp_ctrl_re[[1]], .(ID,Time, Group), summarise,value=mean(value))

dat_cell_exp_ctrl_ctrl_test_sta <- ddply(dat_cell_exp_ctrl_re[[2]], .(ID,Time, Group), summarise,value=mean(value))

dat_cell_exp_ctrl_sta <- ddply(dat_cell_exp_ctrl_test_sta, .(Time, Group), summarise,n=length(value),mean=mean(value),sd=sd(value),se=sd(value)/sqrt(length(value))) %>%
  rbind(., ddply(dat_cell_exp_ctrl_ctrl_test_sta, .(Time, Group), summarise,n=length(value),mean=mean(value),sd=sd(value),se=sd(value)/sqrt(length(value)))) %>%
  mutate(., Day = rep(c("ctrl", "exp_ctrl"),each= length(stim_time)*3))

dat_cell_exp_ctrl_sta$Day <- factor(dat_cell_exp_ctrl_sta$Day, levels = c("ctrl", "exp_ctrl"))
dat_cell_exp_ctrl_sta$Group <- factor(dat_cell_exp_ctrl_sta$Group, levels = c("Neutral", "Excited", "Inhibited"))
range_exp_ctrl_plot <- range(dat_cell_exp_ctrl_sta$mean)

## group by day
p_exp_ctrl <- ggplot(dat_cell_exp_ctrl_sta, aes(Time, mean, colour=Group))+
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


## plot E and I points to show correlation
## points to show the correlation
dat_cell_cor <- NULL

for (i in c(1,2)) {
  dat_exp_ctrl <- subset(dat_cell_exp_ctrl_re[[i]], dat_cell_exp_ctrl_re[[i]]$ID=="m3") %>%
    subset(., .$Group != "Neutral") %>%
    ddply(., .(Time, Group), summarise, value=mean(value)) %>%
    dcast(., Time~Group) %>%
    mutate(., Day = exp_ctrl_group[i])
  dat_cell_cor <- rbind(dat_cell_cor, dat_exp_ctrl)
  
}
dat_cell_cor$Day <- factor(dat_cell_cor$Day, levels = c("ctrl", "exp_ctrl"))

p_trace_cor <- ggplot(dat_cell_cor, aes(Excited, Inhibited, colour=Day))+
  geom_point()+
  geom_smooth(method = "lm", lwd=0.8)+
  scale_colour_manual(values=c("seagreen", "indianred"))+  
  labs(x="Excited (z-score)", y="Inhibited (z-score)")+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  theme(legend.title = element_blank(), legend.position = c(0.8, 0.8))


## correlation for each mice
dat_cell_exp_ctrl_cor <- NULL
for (i in c(1,2)){
  dat_exp_ctrl_cor <- dat_cell_exp_ctrl_re[[i]] %>%
    ddply(., .(ID, Time, Group), summarise, value=mean(value)) %>%
    subset(., .$Group!="Neutral") %>%
    dcast(., Time + ID ~ Group) %>%
    ddply(., .(ID), summarise, "corr" = cor(Excited, Inhibited, method = "spearman")) %>%
    mutate(., Day = exp_ctrl_group[i])
  dat_cell_exp_ctrl_cor <- rbind(dat_cell_exp_ctrl_cor, dat_exp_ctrl_cor)
}

dat_cor_value <- ddply(dat_cell_exp_ctrl_cor, .(ID, Day), summarise, "corr" = mean(corr)) %>%
  mutate(Day = factor(Day, levels = c("ctrl", "exp_ctrl")))

p_cor <-ggplot(dat_cor_value, aes(Day, corr, fill=Day))+
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(width = 0.2, shape=1)+
  scale_fill_manual(values=c("seagreen", "indianred"))+
  labs(x="", y="Correlation coefficient (r)")+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  theme(legend.title = element_blank(), legend.position = "none")


## calculate the sum of active and inhibited trace
dat_cell_exp_ctrl_sum <- rbind(dat_cell_exp_ctrl_test_sta, dat_cell_exp_ctrl_ctrl_test_sta) %>%
  subset(., .$Group!="Neutral") %>%
  mutate( Day = rep(c("ctrl", "exp_ctrl"), each= length(stim_time)*2*length(mouse_ID))) %>%
  ddply(., .(ID,Time, Day), summarise, value=sum(value)) %>%
  ddply(., .(Time, Day), summarise,n=length(value),mean=mean(value),sd=sd(value),se=sd(value)/sqrt(length(value))) %>%
  mutate(Day = factor(Day, levels = c("ctrl", "exp_ctrl")))

p_exp_ctrl_sum <- ggplot(dat_cell_exp_ctrl_sum, aes(Time, mean, colour=Day))+
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

## for EI change analysis
dat_cell_area <- c()
for (i in c(1,2)){
  if (i ==1) {
    rep_time <- c(ncol(dat_pin_m3[[2]]), ncol(dat_pin_m7[[2]]), ncol(dat_pin_m16[[2]]),ncol(dat_pin_m17[[2]]), ncol(dat_pin_m18[[2]]))
  } else{
    rep_time <- mapply(ncol, dat_exp_ctrl)
  }
  
  dat_exp_ctrl11 <- dat_cell_exp_ctrl_re[[i]]
  dat_exp_ctrl11$ID <- rep(mouse_ID, rep_time*length(stim_time))
  dat_exp_ctrl11_sta <- ddply(dat_exp_ctrl11, .(ID, Time, Group), summarise,n=length(value),mean=mean(value),sd=sd(value),se=sd(value)/sqrt(length(value)))
  Time_area <- which(stim_time>=0)
  dat_anti_area <- as.data.frame(tapply(dat_exp_ctrl_sta$mean, INDEX = list(dat_exp_ctrl_sta$ID, dat_exp_ctrl_sta$Group), 
                                        function (x) AUC(stim_time[Time_area], x[Time_area])))
  dat_anti_area$ID <- rownames(dat_anti_area)
  rownames(dat_anti_area)<-NULL
  dat_anti_area$Day <- exp_ctrl_group[i]
  dat_anti_area$ratio <- abs(dat_anti_area$Excited/dat_anti_area$Inhibited)
  dat_anti_area$sum <- dat_anti_area$Excited + dat_anti_area$Inhibited
  dat_cell_area <- rbind(dat_cell_area, dat_anti_area)
}

dat_cell_area$Day <- factor(dat_cell_area$Day, levels = c("ctrl", "exp_ctrl"))
dat_cell_area_sta <-  ddply(dat_cell_area, .(ID, Day), summarise, value=mean(ratio, na.rm = T))

p_EI_ratio <- ggplot(dat_cell_area_sta, aes(Day, value, fill=Day))+
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(width = 0.2, shape=1)+
  scale_fill_manual(values=c("seagreen", "indianred"))+
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

## plot the sum of E-I
dat_cell_sum_sta <-  ddply(dat_cell_area, .(ID, Day), summarise, "sum"=mean(sum, na.rm = T))


p_EI_sum<- ggplot(dat_cell_sum_sta, aes(Day, sum, fill=Day))+
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(width = 0.2, shape=1)+
  scale_fill_manual(values=c("seagreen", "indianred"))+
  labs(x="", y="AUG (z-score)")+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  #scale_y_continuous(limits = c(0, 10), expand = c(0,0))+
  theme(legend.title = element_blank(), legend.position = "none")

## show the change of E, I and sum
p_EIS_change <- dat_cell_area %>%
  select("Excited", "Inhibited","ID","Day", "sum") %>%
  melt(., id.vars= c('Day', "ID")) %>%
  ddply(., .(ID, Day, variable), summarise, value=mean(value, na.rm = T)) %>%
  ggplot(., aes(variable, value, fill=Day))+
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(position=position_jitterdodge(), shape=1, size=1)+
  scale_fill_manual(values=c("seagreen", "indianred"))+
  labs(x="", y="AUC (z-score)")+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  theme(legend.title = element_blank(), legend.position = "top")


## cell catalogs of the neurons during different behavior-----

dat_cell_cat <- mapply(function (x) subset(x, x$Time==0), dat_cell_pin_re, SIMPLIFY = F ) %>% 
  mapply(function(x) prop.table(table(x$ID, x$Group), 1), ., SIMPLIFY = F) %>%
  do.call(rbind, .) %>%
  as.data.frame() %>%
  mutate(ID = rep(rownames(.)[1:5], 5), Day=rep(c("Touch", "Pin","Har", "Shuffle_touch", "Shuffle_pin"), each=5)) %>%
  melt(., id.vars = c("ID", "Day")) %>%
  ddply(., .(ID, Day, variable), summarise, "Prop"=mean(value)) %>%
  mutate(Day = factor(Day, levels = c("Touch", "Pin","Har", "Shuffle_touch", "Shuffle_pin"))) %>%
  mutate(variable=factor(variable, levels = c("Excited", "Neutral", "Inhibited")))

p_cat<- ggplot(dat_cell_cat, aes(variable, Prop, fill=Day))+
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


## correlation analysis------
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

dat_cell_cor_list <- mapply(cc_cor_fun, dat_cell_pin_re, SIMPLIFY = F)


## cor matrix heatmap for m3
cor_range <- range(c(dat_cell_cor_list[[1]][[1]], dat_cell_cor_list[[2]][[1]], dat_cell_cor_list[[3]][[1]], dat_cell_cor_list[[4]][[1]],  dat_cell_cor_list[[5]][[1]] ))

p_cor_touch<- ggcorrplot(dat_cell_cor_list[[1]][[1]],hc.order = TRUE, outline.col = "white")+
  theme_void()+
  scale_fill_gradient2(limit = c(cor_range[1], cor_range[2]), low = "navy", high =  "red4", mid = "white")


p_cor_pin<- ggcorrplot(dat_cell_cor_list[[2]][[1]],hc.order = TRUE, outline.col = "white")+
  theme_void()+
  scale_fill_gradient2(limit = c(cor_range[1], cor_range[2]), low = "navy", high =  "red4", mid = "white")

p_cor_har<- ggcorrplot(dat_cell_cor_list[[3]][[1]],hc.order = TRUE, outline.col = "white")+
  theme_void()+
  scale_fill_gradient2(limit = c(cor_range[1], cor_range[2]), low = "navy", high =  "red4", mid = "white")

p_cor_shuffle_touch<- ggcorrplot(dat_cell_cor_list[[4]][[1]],hc.order = TRUE, outline.col = "white")+
  theme_void()+
  scale_fill_gradient2(limit = c(cor_range[1], cor_range[2]), low = "navy", high =  "red4", mid = "white")

p_cor_shuffle_pin<- ggcorrplot(dat_cell_cor_list[[5]][[1]],hc.order = TRUE, outline.col = "white")+
  theme_void()+
  scale_fill_gradient2(limit = c(cor_range[1], cor_range[2]), low = "navy", high =  "red4", mid = "white")

p_cor_com <- plot_grid(p_cor_touch, p_cor_pin, p_cor_har,p_cor_shuffle_touch,p_cor_shuffle_pin, nrow = 2)



dat_cell_cor <- c()
dat_cell_cor_sta <- c()
for (i in c(1:5)){
  c_trim <- function(x){
    x[upper.tri(x)]<- NA
    c(x[!is.na(x)])
  }
  cor_value <- abs(unlist(mapply(c_trim, dat_cell_cor_list[[i]] )))
  dat_cor <- data.frame(Day=pin_group[i], value=cor_value)
  dat_cell_cor <- rbind(dat_cell_cor, dat_cor)
  
  ## for mean of cor
  cor_value_mean <- dat_cell_cor_list[[i]] %>%
    mapply(function(x) ifelse(is.null(ncol(x)), NA,mean(colMeans(abs(x), na.rm=T))), .)
  
  cor_value_max <- dat_cell_cor_list[[i]] %>%
    mapply(function(x) ifelse(is.null(ncol(x)),NA,mean(apply(x, 2, function(x) max(abs(x), na.rm = T)))), .)
  
  mouse_ID <- c(1: length(dat_cell_cor_list[[i]]))
  dat_cor_mean <- data.frame(Day =pin_group[i] , ID= mouse_ID, value=cor_value_mean, value_max=cor_value_max)
  dat_cell_cor_sta <- rbind(dat_cell_cor_sta, dat_cor_mean)
}


dat_cell_cor$Day <- factor(dat_cell_cor$Day, levels = c("Touch", "Pin","Har", "Shuffle_touch", "Shuffle_pin"))
p_cor_cum<- ggplot(dat_cell_cor, aes(value, group=Day, colour=Day))+
  stat_ecdf(geom = "step")+
  #scale_colour_manual(values=c("seagreen", "indianred", "deepskyblue4"))+
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


dat_cell_cor_sta1 <- melt(dat_cell_cor_sta, id.vars = c("ID", "Day")) %>%
  ddply(., .(ID, Day, variable), summarise, value=mean(value, na.rm = T)) %>%
  mutate(Day = factor(Day, levels = c("Touch", "Pin","Har", "Shuffle_touch", "Shuffle_pin")))


p_cor_mean<- subset(dat_cell_cor_sta1, dat_cell_cor_sta1$variable=="value") %>%
  ggplot(., aes(Day, value, fill=Day))+
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(position=position_jitterdodge(), shape=1, size=1)+
  # scale_fill_manual(values=c("seagreen", "indianred", "deepskyblue4"))+
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
  # scale_fill_manual(values=c("seagreen", "indianred", "deepskyblue4"))+
  labs(x="", y="Max corariance (z-score)")+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  theme(legend.title = element_blank(), legend.position = "none")

