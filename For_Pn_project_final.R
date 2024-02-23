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

## for Fig1-------
# function to analyze the behavior
cc_anti_behavior <- function(path_ID, frame_rate, length_pix) {
  dat_anti <- read.csv(path_ID, header = F, stringsAsFactors = F) %>% 
    .[-c(1:3), ] %>% 
    apply(., 2, as.numeric) %>%   ## only select the position of head
    as.data.frame() %>% 
    mutate(V2 = V2- V8, V3 = V3- V9, V5 = V5-V8, V6 = V6 -V9) %>%  ##normalize from FIJI
    select(str_c("V", 1:6)) 
  
  colnames(dat_anti)<- c("Frame", "Head_x","Head_y", "Prob", "Tail_x", "Tail_y")
  
  ## remove the frame for the first few frame
  dat_anti <- dat_anti[seq(1, nrow(dat_anti), frame_rate/10), ] %>% 
    subset(., Prob >0.5 | Frame> 60)
  
  ## change the dim from pix to mm (the dim of one chammer is 165mm)
  # value 387 for anti experiment 20200309
  # value 447 for miniscope 20191209
  r_pix_mm <- 165 / length_pix
  dat_anti <- dat_anti %>% 
    mutate(Head_x = abs(Head_x * r_pix_mm), Head_y = abs(Head_y * r_pix_mm)) %>% 
    mutate(Tail_x = abs(Tail_x * r_pix_mm), Tail_y = abs(Tail_y * r_pix_mm))
  
  chamber_div <- 165.5 ## the lenght of two chambels and the gap between
  
  ## assign each frame in hot as 1 (48) and nor plate(30) as 2
  dat_anti$chamber <- ifelse(dat_anti$Head_y < chamber_div, "Hot","Nor")
  ratio_time <- unname(by(dat_anti[,3], dat_anti$chamber, length)[2]/nrow(dat_anti))
  
  p_anti <- ggplot(dat_anti, aes(Head_y, Head_x, colour=chamber))+
    geom_point(shape=1)+
    theme_void()+
    theme(legend.position = "none")
  
  ## calculate how many times the mouse cross the border, and latency
  digit <- ifelse(dat_anti$Head_y< chamber_div, 1,2)
  
  for (i in c(2: (length(digit)-20))) {
    digit[i]<- ifelse(digit[i]!= digit[i-1] & digit[i+20] == digit[i], digit[i], digit[i-1] )
  }
  
  dat_anti$digit<- digit
  cross_digit <- diff(digit)
  cross_hot_nor <- which(cross_digit==1)+1
  
  cross_nor_hot <- which(cross_digit==-1)+1
  
  total_cross <- length(c(cross_hot_nor, cross_nor_hot))
  # calculate the latency of fist moving 
  # frame_rate <- 10 ##10 Hz
  cross_latency_1st  <- ifelse(digit[1]==2,2, cross_hot_nor[1]/10)
  cross_latency_1st <- ifelse(cross_latency_1st <2, 2, cross_latency_1st)
  
  # latency of cross back from nor to hot
  cross_latency_2nd <- ifelse(length(cross_nor_hot)==0, nrow(dat_anti)/10 - cross_latency_1st, 
                              cross_nor_hot[1]/10 - cross_latency_1st)
  
  
  ## calculate the distance and speed of movement
  x_diff <- c(0, diff(dat_anti$Head_x))
  y_diff <- c(0, diff(dat_anti$Head_y))
  dis_frame <- sqrt(x_diff^2 + y_diff^2)
  move_speed <- dis_frame*10
  
  move_speed_t <- zoo(move_speed) %>% 
    rollapply(., width = 10, by = 10, FUN = mean, na.rm = TRUE, align = "left") %>% 
    as.vector()
  peak_threshold <- mean(move_speed_t) +  sd(move_speed_t)
  peak_speed <- findpeaks(move_speed_t, nups = 2,minpeakdistance = 2, minpeakheight = peak_threshold) %>% 
    .[,2]*10 %>% 
    sort()
  latency_accerlation <- (peak_speed[which(peak_speed> cross_hot_nor[1])][1] - cross_hot_nor[1])/10
  # total distance, normalization
  total_distance <- round(sum(dis_frame), digits = 2)
  
  # speed in different chambers
  speed_hot <- mean(move_speed[which(dat_anti$chamber == "Hot")])
  speed_nor <- mean(move_speed[which(dat_anti$chamber == "Nor")])
  speed_compare <- speed_nor/speed_hot
  # average speed (3s) before crossing the chambers
  speed_to_nor_1st <- ifelse(cross_hot_nor[1]< 20, mean(move_speed[0:cross_hot_nor[1]]) ,mean(move_speed[(cross_hot_nor[1]-20):cross_hot_nor[1]]))
  speed_to_hot_1st <- ifelse(length(cross_nor_hot)==0, NA, mean(move_speed[(cross_nor_hot[1]-20):cross_nor_hot[1]]))
  
  dat_anti_result <- tibble(ratio=ratio_time, latency1= cross_latency_1st, latency2 = cross_latency_2nd,
                            total_distance = total_distance, speed1 = speed_to_nor_1st, speed2= speed_to_hot_1st, total_cross= total_cross,
                            speed_compare=speed_compare, latency_accerlation= latency_accerlation)
  return(dat_anti_result)
}

# for compare the ctrl and cond. group--
dat_anti_ctrl <- vector(mode = "list", 7)
dat_anti_con <- vector(mode = "list", 7)
length_pix <- c(363,362,363,364,362,364,365)
day_group <- c("Rm","Rm", "Pre", "D4","D5","D6", "Post")
for(i in c(1:7)){
  path_ctrl <- str_c("~cchen/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/wt_ctrl_07042020/wt_ctrl/T",i)
  path_con <-  str_c("~cchen/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/wt_ctrl_07042020/wt_con/T",i)
  dat_anti_ctrl[[i]]<- list.files(path_ctrl,pattern = ".csv", full.names = T ) %>% 
    as.list() %>% 
    mapply(cc_anti_behavior, ., frame_rate = 10, length_pix=length_pix[i], SIMPLIFY = F) %>% 
    do.call(rbind,.) %>% 
    as_tibble() %>% 
    add_column(Day = day_group[i], Group="Ctrl", ID=str_c("m",1:9), .before = 'ratio')
  
  dat_anti_con[[i]]<- list.files(path_con,pattern = ".csv", full.names = T ) %>% 
    as.list() %>% 
    mapply(cc_anti_behavior, ., frame_rate = 10, length_pix=length_pix[i], SIMPLIFY = F) %>% 
    do.call(rbind,.) %>% 
    as_tibble() %>% 
    add_column(Day = day_group[i], Group="Cond.", ID=str_c("m",1:10), .before = 'ratio')
}

dat_anti_wt <- do.call(rbind, dat_anti_ctrl) %>% 
  rbind(., do.call(rbind, dat_anti_con)) %>% 
  filter(Day!="Rm") %>% 
  pivot_longer(-c(Day, Group, ID), names_to = "variable", values_to = "value" , values_drop_na = TRUE) %>% 
  ddply(., .(ID, variable, Day, Group), summarise, 'value'= mean(value)) %>% 
  mutate(Day=factor(Day, levels = c("Pre", "D4", "D5", "D6", "Post")), Group=factor(Group, levels = c("Ctrl", "Cond.")))

dat_anti_wt_sta <- ddply(dat_anti_wt,.(variable,Day, Group), summarise,n=length(value),mean=mean(value),sd=sd(value),se=sd(value)/sqrt(length(value))) %>% 
  filter(Day=="Post")


dat_ratio <- dat_anti_wt %>% 
  filter(variable=="ratio" ) %>% 
  mutate(value=value*100)

p_ratio <- dat_anti_wt %>% 
  filter(variable=="ratio" ) %>% 
  mutate(value=value*100) %>% 
  ddply(.,.(Day, Group), summarise, mean= mean(value), se=sd(value)/sqrt(length(value))) %>% 
  ggplot(data = ., aes(x = Day, y = mean, color = Group))+
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=0.1) +
  geom_point(aes(colour = Group, shape = Group),size=2)+
  geom_line(aes(group = Group))+
  geom_jitter(data= dat_ratio, aes(x= Day, y = value, color= Group, shape = Group), width = 0.2, alpha = 0.4)+
  scale_shape_manual(values=c(20, 18))+
  scale_color_manual(values=c("darkcyan", "indianred"))+
  labs(x="", y="Time spend in chamber 2 (%)")+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  #scale_y_continuous(limits = c(0, 100), expand = c(0, 0))+
  theme(legend.position = "none")

dat_latency <- dat_anti_wt %>% 
  filter(variable=="latency2" )

p_latency2 <- dat_anti_wt %>% 
  filter(variable=="latency2" ) %>% 
  ddply(.,.(Day, Group), summarise, mean= mean(value), se=sd(value)/sqrt(length(value))) %>% 
  ggplot(data = ., aes(x = Day, y = mean, color = Group))+
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=0.1) +
  geom_point(aes(colour = Group, shape = Group),size=2)+
  geom_line(aes(group = Group))+
  geom_jitter(data= dat_latency, aes(x= Day, y = value, color= Group, shape = Group), width = 0.2, alpha = 0.4)+
  scale_shape_manual(values=c(20, 18))+
  scale_color_manual(values=c("darkcyan", "indianred"))+
  labs(x="", y="Latency of crossing back (s)")+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  scale_y_continuous(limits = c(0, 220), expand = c(0, 0))+
  theme(legend.position = "none")

p_wt_pac <- plot_grid(p_ratio, p_latency2, nrow = 1)

setwd("~cchen/Documents/neuroscience/Pn\ project/Figure/PDF/")
cairo_pdf("p_wt_pac.pdf", width = 180/25.6, height = 60/25.6, family = "Arial")
p_wt_pac
dev.off()

p_ratio <- dat_anti_wt %>% 
  filter(variable=="ratio" ) %>% 
  mutate(value=value*100) %>% 
  ggplot(., aes(interaction(Day, Group), value, fill= Day))+
  geom_boxplot( outlier.shape = NA)+
  geom_line(aes(group=interaction(ID, Group)), colour="gray90")+
  geom_point(aes(colour = Group,shape = Group),position=position_jitterdodge(jitter.width = .2), size = 2)+
  scale_fill_manual(values = c("white", "white"))+
  scale_shape_manual(values=c(20, 18))+
  scale_color_manual(values=c("darkcyan", "indianred"))+
  labs(x="", y="Time spend in 2nd plate (%)")+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  scale_y_continuous(limits = c(0, 100), expand = c(0, 0))+
  theme(legend.position = "none")

p_latency2 <- dat_anti_wt %>% 
  filter(variable=="latency2") %>% 
  ggplot(., aes(interaction(Day, Group), value, fill= Day))+
  geom_boxplot( outlier.shape = NA)+
  geom_line(aes(group=interaction(ID, Group)), colour="gray90")+
  geom_point(aes(colour = Group,shape = Group),position=position_jitterdodge(jitter.width = .2))+
  scale_fill_manual(values = c("white", "white"))+
  scale_shape_manual(values=c(20, 18))+
  scale_color_manual(values=c("darkcyan", "indianred"))+
  labs(x="", y="Latency of crossing back (s)")+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  theme(legend.position = "none")

p_pain_compare1 <- plot_grid(p_latency2, p_ratio, nrow = 1)
setwd("~cchen/Documents/neuroscience/Pn\ project/Figure/PDF/")
cairo_pdf("p_pain_compare1.pdf", width = 95/25.6, height = 60/25.6, family = "Arial")
p_pain_compare1
dev.off()

## for pain behavir, like licking, and wearing
# mannualy analyzed by chong (07062020)
dat_wt_manual <- read.xlsx("~cchen2/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/wt_ctrl_07042020/wt_ctrl_07042020_manualy.xlsx") %>% 
  mutate('Rearing'= Total_rearing/Total_frame*600, 'First_rearing'=First_rearing*0.1, 
         'Guarding'= (Guarding - Frame_1st_crossing)*0.1, 'Acceleration'= (Acceleration -Frame_1st_crossing )*0.1) %>% 
  pivot_longer(-c(Group, ID), names_to = "variable", values_to = 'value') %>% 
  mutate(Group=factor(Group, levels=c("Ctrl", "Cond.")))

dat_wt_manual_sta <- dat_wt_manual %>% 
  ddply(., .( variable, Group), summarise,n=length(value),mean=mean(value),sd=sd(value),se=sd(value)/sqrt(length(value)))


p_licking <- dat_wt_manual %>% 
  filter(variable=="First_licking") %>% 
  ggplot(., aes(Group, value, group = Group))+
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(aes(colour= Group,shape = Group),width = 0.2, size = 2)+
  scale_color_manual(values=c("darkcyan", "indianred"))+
  scale_shape_manual(values=c(20, 18))+
  labs(x="", y="Latency to 1st licking (s)")+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  scale_y_continuous(limits = c(0, 50), expand = c(0,0))+
  theme(legend.position = "none")+
  annotate(x=c(1,1,2,2), y=c(40,42,42,40),"path")+
  annotate("text",x=1.5,y=42, label="**", size=5)

p_licking_test <- dat_wt_manual %>% 
  filter(variable=="First_licking") %>% 
  wilcox.test(value~Group, .)

p_1st_rearing <- dat_wt_manual %>% 
  filter(variable=="First_rearing") %>% 
  ggplot(., aes(Group, value, group=Group))+
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(aes(colour= Group,shape = Group),width = 0.2, size = 2)+
  scale_color_manual(values=c("darkcyan", "indianred"))+
  scale_shape_manual(values=c(20, 18))+
  labs(x="", y="Latency to first rearing (s)")+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  theme(legend.position = "none")+
  annotate(x=c(1,1,2,2), y=c(60,62,62,60),"path")+
  annotate("text",x=1.5,y=62, label="*", size=5)

p_1st_rearing_test <- dat_wt_manual %>% 
  filter(variable=="First_rearing") %>% 
  wilcox.test(value~Group,.)

p_rearing <- dat_wt_manual %>% 
  filter(variable=="Rearing") %>% 
  ggplot(., aes(Group, value, group=Group))+
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(aes(colour= Group,shape = Group),width = 0.2, size = 2
              )+
  scale_color_manual(values=c("darkcyan", "indianred"))+
  scale_shape_manual(values=c(20, 18))+
  labs(x="", y="# rearing per min ")+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  scale_y_continuous(limits = c(1, 10))+
  theme(legend.position = "none")+
  annotate(x=c(1,1,2,2), y=c(9,9.1,9.1,9),"path")+
  annotate("text",x=1.5,y=9.1, label="*", size=5)

p_jump <- dat_wt_manual %>% 
  filter(variable=="Jump") %>% 
  ggplot(., aes(Group, value/10, group = Group))+
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(aes(colour= Group,shape = Group),width = 0.2, size = 2)+
  scale_color_manual(values=c("darkcyan", "indianred"))+
  scale_shape_manual(values=c(20, 18))+
  labs(x="", y="Latency of 1st jump (s)")+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  scale_y_continuous(limits = c(30, 200))+
  theme(legend.position = "none")+
  annotate(x=c(1,1,2,2), y=c(185,190,190,185),"path")+
  annotate("text",x=1.5,y=190, label="*", size=5)

p_jump_test<- dat_wt_manual %>% 
  filter(variable=="Jump") %>% 
  wilcox.test(value~Group,.)
p_pain_compare2 <- plot_grid(p_licking, p_1st_rearing, p_jump, nrow = 1)
setwd("~cchen2/Documents/neuroscience/Pn\ project/Figure/PDF/")
cairo_pdf("p_pain_compare2.pdf", width = 90/25.6, height = 60/25.6, family = "Arial")
p_pain_compare2
dev.off()



## manually sorted pain behavior after naloxone injection
dat_nalox_pain <- read.xlsx("~cchen2/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/naloxone_08042020/wt_nalox_08042020_manualy.xlsx") %>% 
  mutate(First_licking = (First_licking - Frame_1st_crossing)/10, First_rearing = (First_rearing-Frame_1st_crossing)/10, 
         Guarding = (Guarding - Frame_1st_crossing)/10, Acceleration = (Acceleration- Frame_1st_crossing)/10, Jump = (Jump- Frame_1st_crossing)/10) %>% 
  pivot_longer(-c(Group, ID), names_to = "variable", values_to = "value" , values_drop_na = TRUE) %>% 
  mutate(Group = factor(Group, levels = c("Ctrl", "Nalox.")))

dat_nalox_pain_sta <- dat_nalox_pain %>% 
  ddply(.,.(variable, Group), summarise,n=length(value),mean=mean(value),sd=sd(value),se=sd(value)/sqrt(length(value)))

p_licking <- dat_nalox_pain %>% 
  filter(variable=="First_licking") %>% 
  ggplot(., aes(Group, value, group=Group))+
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(aes(colour= Group,shape = Group), width = 0.2)+
  scale_shape_manual(values=c(18, 17))+
  scale_color_manual(values=c("indianred", "mediumpurple4"))+
  labs(x="", y="Latency to 1st licking (s)")+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  scale_y_continuous(limits = c(0, 70), expand = c(0,0))+
  theme(legend.position = "none")+
  annotate(x=c(1,1,2,2), y=c(62,63,63,62),"path")+
  annotate("text",x=1.5,y=63, label="**", size=5)

p_licking_test <- dat_nalox_pain %>% 
  filter(variable=="First_licking") %>% 
  wilcox.test(value~Group, .)

p_1st_rearing <- dat_nalox_pain %>% 
  filter(variable=="First_rearing") %>% 
  ggplot(., aes(Group, value, group=Group))+
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(aes(colour= Group,shape = Group), width = 0.2)+
  scale_shape_manual(values=c(18, 17))+
  scale_color_manual(values=c("indianred", "mediumpurple4"))+
  labs(x="", y="Latency to first rearing (s)")+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  theme(legend.position = "none")

p_1st_rearing_test <- dat_nalox_pain %>% 
  filter(variable=="First_rearing") %>% 
  wilcox.test(value~Group,.)


p_jump <- dat_nalox_pain %>% 
  filter(variable=="Jump1") %>% 
  ggplot(., aes(Group, value/10, group=Group))+
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(aes(colour= Group,shape = Group), width = 0.2)+
  scale_shape_manual(values=c(18, 17))+
  scale_color_manual(values=c("indianred", "mediumpurple4"))+
  labs(x="", y="Latency to 1st jump (s) ")+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  theme(legend.position = "none")+
  annotate(x=c(1,1,2,2), y=c(190,195,195,190),"path")+
  annotate("text",x=1.5,y=195, label="**", size=5)

p_jump_test<- dat_nalox_pain %>% 
  filter(variable=="Jump1") %>% 
  wilcox.test(value~Group,.)

p_pain_compare2 <- plot_grid(p_licking, p_1st_rearing, p_jump, nrow = 1)

setwd("~cchen2/Documents/neuroscience/Pn\ project/Figure/PDF/")
cairo_pdf("p_con_nalox.pdf", width = 100/25.6, height = 60/25.6, family = "Arial")
p_pain_compare2
dev.off()



## For Figure s2----
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

## Extract trace from each mice(from D3-D6)------
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

# for m857
path_trace_m857 <- as.list(list.files(path ="~cchen2/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/matlab_analysis/m857/trace_days", pattern = "*.xlsx", full.names = T ))
t_stim_m857_d1 <- c(223)*2
t_stim_m857_d2 <- c(513)*2
t_stim_m857_d3 <- c(847) *2
t_stim_m857_d4 <- c(252)*2
t_stim_m857_d5 <- c(1900)*2 
t_stim_m857_d6 <- c(723) *2
t_stim_m857_d7 <- c(167)*2

t_stim_m857 <- list(t_stim_m857_d1, t_stim_m857_d2, t_stim_m857_d3, t_stim_m857_d4, t_stim_m857_d5, t_stim_m857_d6, t_stim_m857_d7)

dat_trace_m857 <- mapply(c_miniscope_matlab, path_trace_m855, t_stim_m855, SIMPLIFY = F)

## combine data and do k-means analysis-----

dat_cell_trace <- vector(mode = "list", 7)
for (i in 1:7) {
  dat_cell_trace[[i]] <- cbind(dat_trace_m3[[i]], dat_trace_m7[[i]], dat_trace_m17[[i]], dat_trace_m18[[i]], dat_trace_m855[[i]], dat_trace_m857[[i]])
  
}


dat_cell_trace_re <- vector(mode = "list", 7)

# k-menas clustering
mouse_ID <- c("m3", "m7", "m17", "m18", "m855", "m857")
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
  
  rep_time <- c(ncol(dat_trace_m3[[i]]), ncol(dat_trace_m7[[i]]), ncol(dat_trace_m17[[i]]), ncol(dat_trace_m18[[i]]), ncol(dat_trace_m855[[i]]), ncol(dat_trace_m855[[i]]) )
  dat_cell_trace_d_re <- mutate(dat_cell_trace_d_re, Group= factor(Group, levels = c("Excited", "Neutral", "Inhibited"))) %>% 
    mutate(., ID = rep(mouse_ID, rep_time*length(stim_time)) )
  
  dat_cell_trace_re[[i]] <- dat_cell_trace_d_re
}

dat_cell_trace_pre_sta <-  ddply(dat_cell_trace_re[[3]], .(ID,Time, Group), summarise,value=mean(value)) 
dat_cell_trace_d4_sta <- ddply(dat_cell_trace_re[[4]], .(ID,Time, Group), summarise,value=mean(value))
dat_cell_trace_d5_sta <- ddply(dat_cell_trace_re[[5]], .(ID,Time, Group), summarise,value=mean(value))
dat_cell_trace_d6_sta <- ddply(dat_cell_trace_re[[6]], .(ID,Time, Group), summarise,value=mean(value))

score_range <- rbind(dat_cell_trace_re[[3]], dat_cell_trace_re[[4]], dat_cell_trace_re[[5]], dat_cell_trace_re[[6]]) %>% 
  .$value %>% 
  range()
group_day <- c("D3",  "D4", "D5", "D6")

for (i in c(3: 6)) {
  dat_trace <- dat_cell_trace_re[[i]]
  dat_trace_sta <- ddply(dat_trace, .(variable, Group), summarise,mean=mean(value), sum=sum(value))
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
  assign(str_c("p_heat_d", i), p_heat)
}

## combine the heat plot
p_heat_com <- plot_grid(p_heat_d3, p_heat_d4,p_heat_d5, p_heat_d6, nrow = 1)

setwd("~cchen2/Documents/neuroscience/Pn\ project/Figure/PDF/")
cairo_pdf("p_heat_com.pdf", width = 120/25.6, height = 55/25.6, family = "Arial")
p_heat_com
dev.off()

dat_cell_trace_sum <- rbind(dat_cell_trace_pre_sta, dat_cell_trace_d4_sta, dat_cell_trace_d5_sta, dat_cell_trace_d6_sta) %>%
  as_tibble() %>% 
  mutate(Day = rep(c("D3",  "D4", "D5", "D6"), c(nrow(dat_cell_trace_pre_sta),nrow(dat_cell_trace_d4_sta), nrow(dat_cell_trace_d5_sta), nrow(dat_cell_trace_d6_sta)))) %>%
  subset(., .$Group!="Neutral") %>%
  ddply(., .(ID,Time, Day), summarise, value=sum(value)) %>%
  ddply(., .(Time, Day), summarise,n=length(value),mean=mean(value),sd=sd(value),se=sd(value)/sqrt(length(value))) %>%
  mutate(Day = factor(Day, levels = c("D3",  "D4", "D5", "D6")))

p_trace_sum <- ggplot(dat_cell_trace_sum, aes(Time, mean, colour=Day))+
  geom_line()+
  geom_ribbon(aes(ymin=mean-se, ymax=mean+se, fill=Day), alpha=0.2, linetype=0)+
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
cairo_pdf("p_trace_sum_compare.pdf", width = 90/25.6, height = 60/25.6, family = "Arial")
p_trace_sum
dev.off()

## calculate the AUC
dat_cell_area <- c()
group_day <- c("D3",  "D4", "D5", "D6")
for (i in c(3:6)){
  rep_time <- c(ncol(dat_trace_m3[[i]]), ncol(dat_trace_m7[[i]]), ncol(dat_trace_m17[[i]]), ncol(dat_trace_m18[[i]]), ncol(dat_trace_m855[[i]]), ncol(dat_trace_m857[[i]]))
  dat_trace <- dat_cell_trace_re[[i]]
  dat_trace$ID <- rep(mouse_ID, rep_time*length(stim_time))
  dat_trace_sta <- ddply(dat_trace, .(ID, Time, Group), summarise,n=length(value),mean=mean(value),sd=sd(value),se=sd(value)/sqrt(length(value)))
  Time_area <- which(stim_time>=0)
  dat_anti_area <- as.data.frame(tapply(dat_trace_sta$mean, INDEX = list(dat_trace_sta$ID, dat_trace_sta$Group), 
                                        function (x) AUC(stim_time[Time_area], x[Time_area])))
  dat_anti_area$ID <- rownames(dat_anti_area)
  rownames(dat_anti_area)<-NULL
  dat_anti_area$Day <- str_c("D", i)
  dat_anti_area$ratio <- abs(dat_anti_area$Excited/dat_anti_area$Inhibited)
  dat_anti_area$sum <- dat_anti_area$Excited + dat_anti_area$Inhibited
  dat_cell_area <- rbind(dat_cell_area, dat_anti_area)
}

dat_cell_area$Day <- factor(dat_cell_area$Day, levels = c("D3",  "D4", "D5", "D6"))


dat_EI_sum_sta <- dat_cell_area %>% 
  select(Day, ID, sum) %>% 
  ddply(., .(Day), summarise,n=length(sum),mean=mean(sum, na.rm = T),sd=sd(sum,na.rm = T),se=sd(sum, na.rm = T)/sqrt(length(sum)))

p_EI_sum<- ggplot(dat_EI_sum_sta, aes(x = Day, y = mean, colour = Day))+
  geom_bar(stat="identity", fill = "white")+
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=0.1) +
  scale_colour_manual(values=c("deepskyblue4", "indianred" ))+
  labs(x="", y="AUG (z-score)")+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  theme(legend.title = element_blank(), legend.position = "none")

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
  labs(x="", y="AUG (z-score)")+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  theme(legend.title = element_blank(), legend.position = "none")

## for fig s1------
dat_anti_ctrl <- vector(mode = "list", 7)
dat_anti_con <- vector(mode = "list", 7)
length_pix <- c(363,362,363,364,362,364,365)
day_group <- c("Rm","Rm", "Pre", "D4","D5","D6", "Test")
for(i in c(1:7)){
  path_ctrl <- str_c("~cchen2/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/wt_ctrl_07042020/wt_ctrl/T",i)
  path_con <-  str_c("~cchen2/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/wt_ctrl_07042020/wt_con/T",i)
  dat_anti_ctrl[[i]]<- list.files(path_ctrl,pattern = ".csv", full.names = T ) %>% 
    as.list() %>% 
    mapply(cc_anti_behavior, ., frame_rate = 10, length_pix=length_pix[i], SIMPLIFY = F) %>% 
    do.call(rbind,.) %>% 
    add_column(Day = day_group[i], Group="Ctrl", ID=str_c("m",1:9), .before = 'ratio')
  
  dat_anti_con[[i]]<- list.files(path_con,pattern = ".csv", full.names = T ) %>% 
    as.list() %>% 
    mapply(cc_anti_behavior, ., frame_rate = 10, length_pix=length_pix[i], SIMPLIFY = F) %>% 
    do.call(rbind,.) %>% 
    add_column(Day = day_group[i], Group="Cond.", ID=str_c("m",1:10), .before = 'ratio')
}

dat_anti_wt <- do.call(rbind, dat_anti_ctrl) %>% 
  rbind(., do.call(rbind, dat_anti_con)) %>% 
  filter(Day!="Rm") %>% 
  pivot_longer(-c(Day, Group, ID), names_to = "variable", values_to = "value" , values_drop_na = TRUE) %>% 
  ddply(., .(ID, variable, Day, Group), summarise, 'value'= mean(value)) %>% 
  mutate(Day=factor(Day, levels = c("Pre", "D4", "D5", "D6", "Test")), Group=factor(Group, levels = c("Ctrl", "Cond."))) 

p_latency1 <- dat_anti_wt %>% 
  filter(variable=="latency1") %>% 
  ddply(., .(Group, Day) , summarise,n=length(value),mean=mean(value),sd=sd(value),se=sd(value)/sqrt(length(value))) %>% 
  ggplot(., aes(Day, mean, group=Group, colour = Group))+
  geom_line()+
  geom_point(aes(shape = Group, colour = Group), size = 2)+ 
  geom_errorbar(aes(ymin=mean-se,ymax=mean+se, colour=Group),width=.2)+
  scale_shape_manual(values=c(20, 18))+
  scale_color_manual(values=c("darkcyan", "indianred"))+
  labs(x="", y="Latency (s)")+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  scale_y_continuous(limits = c(0, 150), expand = c(0, 0))+
  theme(legend.position = "none")

aov_test_latency1 <- dat_anti_wt %>% 
  filter(variable=="latency1") %>% 
  aov(value ~ Group + Day + Group:Day, .)
summary(aov_test_latency1 )

TukeyHSD(aov_test_latency1, which = "Group:Day")


p_latency2 <- dat_anti_wt %>% 
  filter(variable=="latency2") %>% 
  ddply(., .(Group, Day) , summarise,n=length(value),mean=mean(value),sd=sd(value),se=sd(value)/sqrt(length(value))) %>% 
  ggplot(., aes(Day, mean, group=Group, colour = Group))+
  geom_line()+
  geom_point(aes(shape = Group, colour = Group), size = 2)+ 
  geom_errorbar(aes(ymin=mean-se,ymax=mean+se, colour=Group),width=.2)+
  scale_shape_manual(values=c(20, 18))+
  scale_color_manual(values=c("darkcyan", "indianred"))+
  labs(x="", y="Latency (s)")+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  scale_y_continuous(limits = c(0, 150))+
  theme(legend.position = "none")

aov_test_latency2 <- dat_anti_wt %>% 
  filter(variable=="latency2") %>% 
  aov(value ~ Group + Day + Group:Day, .)
summary(aov_test_latency2 )

TukeyHSD(aov_test_latency2, which = "Group:Day")

p_ratio <- dat_anti_wt %>% 
  filter(variable=="ratio") %>% 
  ddply(., .(Group, Day) , summarise,n=length(value),mean=mean(value),sd=sd(value),se=sd(value)/sqrt(length(value))) %>% 
  ggplot(., aes(Day, mean, group=Group, colour = Group))+
  geom_line()+
  geom_point(aes(shape = Group, colour = Group), size = 2)+ 
  geom_errorbar(aes(ymin=mean-se,ymax=mean+se, colour=Group),width=.2)+
  scale_shape_manual(values=c(20, 18))+
  scale_color_manual(values=c("darkcyan", "indianred"))+
  labs(x="", y="Time spend in 2nd plate (%)")+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  theme(legend.position = "none")

aov_test_ratio <- dat_anti_wt %>% 
  filter(variable=="ratio") %>% 
  aov(value ~ Group + Day + Group:Day, .)
summary(aov_test_ratio )

TukeyHSD(aov_test_ratio, which = "Group:Day")

