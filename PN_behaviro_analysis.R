## For anti behavior tracking data analysis
## behavior videos are analyzed with deeplabcut
## created on 03192020
## last updated 03242020

cc_anti_behavior <- function(path_ID, frame_rate, length_pix) {
  dat_anti <- read.csv(path_ID, header = F, stringsAsFactors = F) %>% 
    .[-c(1:3), ] %>% 
    apply(., 2, as.numeric) %>%   ## only select the position of head
    as.data.frame()


  ## normalize the data, zero point value from fiji
  dat_anti$V2<- dat_anti$V2 - dat_anti$V8
  dat_anti$V3<- dat_anti$V3 - dat_anti$V9
  dat_anti$V5 <- dat_anti$V5 - dat_anti$V8
  dat_anti$V6 <- dat_anti$V6 - dat_anti$V9
  dat_anti <- dat_anti[,1:6]
  
  #colnames(dat_anti)<- c("Frame", "Head_x","Head_y","likehood", 'Tail_x',"Tail_y", "likehood_tail")
  colnames(dat_anti)<- c("Frame", "Head_x","Head_y", "Prob", "Tail_x", "Tail_y")
  
  ## remove the frame for the first few frame
  dat_anti <- dat_anti[seq(1, nrow(dat_anti), frame_rate/10), ] %>% 
    subset(., Prob >0.5 | Frame> 60)
  
  ## reduce the frame rate to 5hz to mimimize the varaition 
  ##n<- 2
  ## dat_anti<- aggregate(dat_anti,list(rep(1:(nrow(dat_anti)%/%n+1),each=n,len=nrow(dat_anti))),mean)[-1]
  
  
  ## change the dim from pix to mm (the dim of one chammer is 165mm)
  # value 387 for anti experiment 20200309
  # value 447 for miniscope 20191209
  r_pix_mm <- 165 / length_pix
  dat_anti$Head_x <- abs(dat_anti$Head_x * r_pix_mm)
  dat_anti$Head_y <- abs(dat_anti$Head_y * r_pix_mm)
  
  
  dat_anti$Tail_x <- abs(dat_anti$Tail_x * r_pix_mm)
  dat_anti$Tail_y <- abs(dat_anti$Tail_y * r_pix_mm)
  
  chamber_div <- 165.5 ## the lenght of two chambels and the gap between
  
  ## assign each frame in hot as 1 (48) and nor plate(30) as 2
  dat_anti$chamber <- ifelse(dat_anti$Head_y < chamber_div, "Hot","Nor")
  ratio_time <- unname(by(dat_anti[,3], dat_anti$chamber, length)[2]/nrow(dat_anti))
  
  p_anti <- ggplot(dat_anti, aes(Head_y, Head_x))+
    geom_point(shape=1, colour="brown4")+
    theme_void()+
    theme(legend.position = "none")
  
  ## calculate how many times the mouse cross the border, and latency
  digit <- ifelse(dat_anti$Head_y< chamber_div, 1,2)
  
  for (i in c(2: (length(digit)-15))) {
    digit[i]<- ifelse(digit[i]!= digit[i-1] & digit[i+15] == digit[i], digit[i], digit[i-1] )
  }
  
  dat_anti$digit<- digit
  cross_digit <- diff(digit)
  cross_hot_nor <- which(cross_digit==1)+1
  
  cross_nor_hot <- which(cross_digit==-1)+1
  
  total_cross <- length(c(cross_hot_nor, cross_nor_hot))
  # calculate the latency of fist moving 
  # frame_rate <- 10 ##10 Hz
  cross_latency_1st  <- cross_hot_nor[1]/10
  
  # latency of cross back from nor to hot
  cross_latency_2nd <- ifelse(length(cross_nor_hot)==0, nrow(dat_anti)/10 - cross_latency_1st, 
                              (cross_nor_hot[1] - cross_hot_nor[1])/10)
  
  
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

## for opto anti behavior------
## for mice m20
setwd("~cchen2/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/anti_opto_behavior/m20")

path_trace_m20 <- list.files( pattern = "*.csv" )
dat_behavior_m20 <- mapply(cc_anti_behavior, path_trace_m20, frame_rate = 10,length_pix=387, SIMPLIFY = F)

## for mice m19
setwd("~cchen2/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/anti_opto_behavior/m19")

path_trace_m19 <- list.files( pattern = "*.csv" )
dat_behavior_m19 <- mapply(cc_anti_behavior, path_trace_m19, frame_rate = 10, length_pix=387,SIMPLIFY = F)

## for mice m22
setwd("~cchen2/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/anti_opto_behavior/m22")

path_trace_m22 <- list.files( pattern = "*.csv" )
dat_behavior_m22 <- mapply(cc_anti_behavior, path_trace_m22, frame_rate = 10, length_pix=387,SIMPLIFY = F)

## for mice m23
setwd("~cchen2/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/anti_opto_behavior/m23")

path_trace_m23 <- list.files( pattern = "*.csv" )
dat_behavior_m23 <- mapply(cc_anti_behavior, path_trace_m23, frame_rate = 10, length_pix=387,SIMPLIFY = F)

## for mice m63
setwd("~cchen2/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/anti_opto_behavior/m63")

path_trace_m63 <- list.files( pattern = "*.csv" )
dat_behavior_m63 <- mapply(cc_anti_behavior, path_trace_m63, frame_rate = 10, length_pix=387,SIMPLIFY = F)

## for mice m851
setwd("~cchen2/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/anti_opto_behavior/m851")

path_trace_m851 <- list.files( pattern = "*.csv" )
dat_behavior_m851 <- mapply(cc_anti_behavior, path_trace_m851, frame_rate = 10, length_pix=387,SIMPLIFY = F)


## for miniscope implanted behavior----
# for m16
setwd("~cchen2/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/matlab_analysis/behavior_tracking/m16")

path_trace_m16 <- list.files( pattern = "*.csv" )
dat_behavior_m16 <- mapply(cc_anti_behavior, path_trace_m16, frame_rate = 30, length_pix=387,SIMPLIFY = F)

# for m17
setwd("~cchen2/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/matlab_analysis/behavior_tracking/m17")

path_trace_m17 <- list.files( pattern = "*.csv" )
dat_behavior_m17 <- mapply(cc_anti_behavior, path_trace_m17, frame_rate = 30, length_pix=387,SIMPLIFY = F)

# for m18
setwd("~cchen2/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/matlab_analysis/behavior_tracking/m18")

path_trace_m18 <- list.files( pattern = "*.csv" )
dat_behavior_m18 <- mapply(cc_anti_behavior, path_trace_m18, frame_rate = 30, length_pix=387,SIMPLIFY = F)

# for m3
setwd("~cchen2/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/matlab_analysis/behavior_tracking/m3")

path_trace_m3 <- list.files( pattern = "*.csv" )
dat_behavior_m3 <- mapply(cc_anti_behavior, path_trace_m3, frame_rate = 30, length_pix=478,SIMPLIFY = F)

# for m7
setwd("~cchen2/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/matlab_analysis/behavior_tracking/m7")

path_trace_m7 <- list.files( pattern = "*.csv" )
dat_behavior_m7 <- mapply(cc_anti_behavior, path_trace_m7, frame_rate = 30, length_pix=478,SIMPLIFY = F)

# for m27
setwd("~cchen2/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/matlab_analysis/behavior_tracking/m27")

path_trace_m27 <- list.files( pattern = "*.csv" )
dat_behavior_m27 <- mapply(cc_anti_behavior, path_trace_m27, frame_rate = 30, length_pix=393,SIMPLIFY = F)


## for HM4Di behavior------
# for m10
setwd("~cchen2/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/anti_hm4di/m10")

path_trace_m10 <- list.files( pattern = "*.csv" )
dat_behavior_m10 <- mapply(cc_anti_behavior, path_trace_m10, frame_rate = 10, length_pix=455, SIMPLIFY = F)

## for m11
setwd("~cchen2/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/anti_hm4di/m11")

path_trace_m11 <- list.files( pattern = "*.csv" )
dat_behavior_m11 <- mapply(cc_anti_behavior, path_trace_m11, frame_rate = 10, length_pix=455, SIMPLIFY = F)

## for m12
setwd("~cchen2/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/anti_hm4di/m12")

path_trace_m12 <- list.files( pattern = "*.csv" )
dat_behavior_m12 <- mapply(cc_anti_behavior, path_trace_m12, frame_rate = 10, length_pix=455, SIMPLIFY = F)

# for m13
setwd("~cchen2/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/anti_hm4di/m13")

path_trace_m13 <- list.files( pattern = "*.csv" )
dat_behavior_m13 <- mapply(cc_anti_behavior, path_trace_m13, frame_rate = 10, length_pix=455, SIMPLIFY = F)


## combine all anti data and plot----
mouse_anti_list <- list(dat_behavior_m16, dat_behavior_m17, dat_behavior_m18, dat_behavior_m3, dat_behavior_m7, dat_behavior_m27)

## for 1st latency, day2-3 as Ctrl, day5-6 as anti and day7 in pain
latency_1st_ctrl_d2 <- sapply(mouse_anti_list, function(x) x[[2]][[2]])
latency_1st_ctrl_d3 <- sapply(mouse_anti_list, function(x) x[[3]][[2]])
latency_1st_ctrl <- rowMeans(cbind(latency_1st_ctrl_d2, latency_1st_ctrl_d3))

latency_1st_ctrl_d5 <- sapply(mouse_anti_list, function(x) x[[5]][[2]])
latency_1st_ctrl_d6 <- sapply(mouse_anti_list, function(x) x[[6]][[2]])
latency_1st_anti <- rowMeans(cbind(latency_1st_ctrl_d5, latency_1st_ctrl_d6))

latency_1st_anti_pain <- sapply(mouse_anti_list, function(x) x[[7]][[2]])

anti_latency <- c(latency_1st_ctrl, latency_1st_anti, latency_1st_anti_pain)
anti_group <- rep(c("Ctrl", "Anti", "Anti_pain"), each=length(latency_1st_ctrl))
dat_anti_1st_latency <- data.frame(Group = anti_group, Latency = anti_latency)
dat_anti_1st_latency$Group <- factor(dat_anti_1st_latency$Group, levels = c("Ctrl", "Anti", "Anti_pain"))
# dat_anti_1st_latency_sta <- ddply(dat_anti_1st_latency, .(Group), summarise,n=length(Latency),mean=mean(Latency),sd=sd(Latency),se=sd(Latency)/sqrt(length(Latency)))

p_anti_1st_latency <- ggplot(dat_anti_1st_latency, aes(Group, Latency, colour=Group))+
  geom_boxplot(outlier.shape =NA )+
  geom_jitter(width = 0.25, shape=1)+
  labs(x="", y="Latency of 1st border crossing (s)")+
  theme(axis.line = element_line(colour = "black"),
        axis.text.x = element_text(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  theme(legend.position = 'none')

## for 2nd latencyv(crossing back), day2-3 as Ctrl, day5-6 as anti and day7 in pain
latency_2nd_ctrl_d2 <- sapply(mouse_anti_list, function(x) x[[2]][[3]])
latency_2nd_ctrl_d3 <- sapply(mouse_anti_list, function(x) x[[3]][[3]])
latency_2nd_ctrl <- rowMeans(cbind(latency_2nd_ctrl_d2, latency_2nd_ctrl_d3))

latency_2nd_ctrl_d5 <- sapply(mouse_anti_list, function(x) x[[5]][[3]])
latency_2nd_ctrl_d6 <- sapply(mouse_anti_list, function(x) x[[6]][[3]])
latency_2nd_anti <- rowMeans(cbind(latency_2nd_ctrl_d5, latency_2nd_ctrl_d6))

latency_2nd_anti_pain <- sapply(mouse_anti_list, function(x) x[[7]][[3]])

anti_2nd_latency <- c(latency_2nd_ctrl, latency_2nd_anti, latency_2nd_anti_pain)
anti_group <- rep(c("Ctrl", "Anti", "Anti_pain"), each=length(latency_2nd_ctrl))
dat_anti_2nd_latency <- data.frame(Group = anti_group, Latency = anti_2nd_latency)
dat_anti_2nd_latency$Group <- factor(dat_anti_2nd_latency$Group, levels = c("Ctrl", "Anti", "Anti_pain"))
# dat_anti_2nd_latency_sta <- ddply(dat_anti_2nd_latency, .(Group), summarise,n=length(Latency),mean=mean(Latency),sd=sd(Latency),se=sd(Latency)/sqrt(length(Latency)))

p_anti_2nd_latency <- ggplot(dat_anti_2nd_latency, aes(Group, Latency, colour=Group))+
  geom_boxplot(outlier.shape =NA )+
  geom_jitter(width = 0.25, shape=1)+
  labs(x="", y="Latency of crossing back (s)")+
  theme(axis.line = element_line(colour = "black"),
        axis.text.x = element_text(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  theme(legend.position = 'none')

## for 1st speed, day2-3 as Ctrl, day5-6 as anti and day7 in pain
speed_1st_ctrl_d2 <- sapply(mouse_anti_list, function(x) x[[2]][[5]])
speed_1st_ctrl_d3 <- sapply(mouse_anti_list, function(x) x[[3]][[5]])
speed_1st_ctrl <- rowMeans(cbind(speed_1st_ctrl_d2, speed_1st_ctrl_d3))

speed_1st_ctrl_d5 <- sapply(mouse_anti_list, function(x) x[[5]][[5]])
speed_1st_ctrl_d6 <- sapply(mouse_anti_list, function(x) x[[6]][[5]])
speed_1st_anti <- rowMeans(cbind(speed_1st_ctrl_d5, speed_1st_ctrl_d6))

speed_1st_anti_pain <- sapply(mouse_anti_list, function(x) x[[7]][[5]])

anti_speed <- c(speed_1st_ctrl, speed_1st_anti, speed_1st_anti_pain)
anti_group <- rep(c("Ctrl", "Anti", "Anti_pain"), each=length(speed_1st_ctrl))
dat_anti_1st_speed <- data.frame(Group = anti_group, speed = anti_speed)
dat_anti_1st_speed$Group <- factor(dat_anti_1st_speed$Group, levels = c("Ctrl", "Anti", "Anti_pain"))
# dat_anti_1st_speed_sta <- ddply(dat_anti_1st_speed, .(Group), summarise,n=length(speed),mean=mean(speed),sd=sd(speed),se=sd(speed)/sqrt(length(speed)))

p_anti_1st_speed <- ggplot(dat_anti_1st_speed, aes(Group, speed, colour=Group))+
  geom_boxplot(outlier.shape =NA )+
  geom_jitter(width = 0.25, shape=1)+
  labs(x="", y="Speed of crossing (mm/s)")+
  theme(axis.line = element_line(colour = "black"),
        axis.text.x = element_text(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  theme(legend.position = 'none')

## for 2nd speed (crossing back), day2-3 as Ctrl, day5-6 as anti and day7 in pain
speed_2nd_ctrl_d2 <- sapply(mouse_anti_list, function(x) x[[2]][[6]])
speed_2nd_ctrl_d3 <- sapply(mouse_anti_list, function(x) x[[3]][[6]])
speed_2nd_ctrl <- rowMeans(cbind(speed_2nd_ctrl_d2, speed_2nd_ctrl_d3), na.rm = T)

speed_2nd_ctrl_d5 <- sapply(mouse_anti_list, function(x) x[[5]][[6]])
speed_2nd_ctrl_d6 <- sapply(mouse_anti_list, function(x) x[[6]][[6]])[-4]
speed_2nd_anti <- rowMeans(cbind(speed_2nd_ctrl_d5, speed_2nd_ctrl_d6), na.rm = T)
speed_2nd_anti <- speed_2nd_anti[!is.na(speed_2nd_anti)]
speed_2nd_anti_pain <- sapply(mouse_anti_list, function(x) x[[7]][[6]])

anti_2nd_speed <- c(speed_2nd_ctrl, speed_2nd_anti, speed_2nd_anti_pain)
anti_group_2nd <- rep(c("Ctrl", "Anti", "Anti_pain"), c(length(speed_2nd_ctrl), length(speed_2nd_anti), length(speed_2nd_anti_pain)))
dat_anti_2nd_speed <- data.frame(Group = anti_group_2nd, speed = anti_2nd_speed)
dat_anti_2nd_speed$Group <- factor(dat_anti_2nd_speed$Group, levels = c("Ctrl", "Anti", "Anti_pain"))
# dat_anti_2nd_speed_sta <- ddply(dat_anti_2nd_speed, .(Group), summarise,n=length(speed),mean=mean(speed),sd=sd(speed),se=sd(speed)/sqrt(length(speed)))

p_anti_2nd_speed <- ggplot(dat_anti_2nd_speed, aes(Group, speed, colour=Group))+
  geom_boxplot(outlier.shape =NA )+
  geom_jitter(width = 0.25, shape=1)+
  labs(x="", y="Speed of crossing (mm/s)")+
  theme(axis.line = element_line(colour = "black"),
        axis.text.x = element_text(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  theme(legend.position = 'none')


## for total number of crossing, day2-3 as Ctrl, day5-6 as anti and day7 in pain
num_cross_ctrl_d2 <- sapply(mouse_anti_list, function(x) x[[2]][[7]])
num_cross_ctrl_d3 <- sapply(mouse_anti_list, function(x) x[[3]][[7]])
num_cross_ctrl <- rowMeans(cbind(num_cross_ctrl_d2, num_cross_ctrl_d3))

num_cross_ctrl_d5 <- sapply(mouse_anti_list, function(x) x[[5]][[7]])
num_cross_ctrl_d6 <- sapply(mouse_anti_list, function(x) x[[6]][[7]])
num_cross_anti <- rowMeans(cbind(num_cross_ctrl_d5, num_cross_ctrl_d6))

num_cross_anti_pain <- sapply(mouse_anti_list, function(x) x[[7]][[7]])

anti_num_cross <- c(num_cross_ctrl, num_cross_anti, num_cross_anti_pain)
anti_group <- rep(c("Ctrl", "Anti", "Anti_pain"), each=length(num_cross_ctrl))
dat_anti_num_cross <- data.frame(Group = anti_group, num = anti_num_cross)
dat_anti_num_cross$Group <- factor(dat_anti_num_cross$Group, levels = c("Ctrl", "Anti", "Anti_pain"))
# dat_anti_1st_speed_sta <- ddply(dat_anti_1st_speed, .(Group), summarise,n=length(speed),mean=mean(speed),sd=sd(speed),se=sd(speed)/sqrt(length(speed)))

p_anti_num_cross <- ggplot(dat_anti_num_cross, aes(Group, num, colour=Group))+
  geom_boxplot(outlier.shape =NA )+
  geom_jitter(width = 0.25, shape=1)+
  labs(x="", y="Total crossing #")+
  theme(axis.line = element_line(colour = "black"),
        axis.text.x = element_text(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  theme(legend.position = 'none')

## for ratio of time in Nor chammber, day2-3 as Ctrl, day5-6 as anti and day7 in pain
ratio_time_ctrl_d2 <- sapply(mouse_anti_list, function(x) x[[2]][[1]])
ratio_time_ctrl_d3 <- sapply(mouse_anti_list, function(x) x[[3]][[1]])
ratio_time_ctrl <- rowMeans(cbind(ratio_time_ctrl_d2, ratio_time_ctrl_d3))

ratio_time_ctrl_d5 <- sapply(mouse_anti_list, function(x) x[[5]][[1]])
ratio_time_ctrl_d6 <- sapply(mouse_anti_list, function(x) x[[6]][[1]])
ratio_time_anti <- rowMeans(cbind(ratio_time_ctrl_d5, ratio_time_ctrl_d6))

ratio_time_anti_pain <- sapply(mouse_anti_list, function(x) x[[7]][[1]])

anti_ratio_time <- c(ratio_time_ctrl, ratio_time_anti, ratio_time_anti_pain)
anti_group <- rep(c("Ctrl", "Anti", "Anti_pain"), each=length(ratio_time_ctrl))
dat_anti_ratio_time <- data.frame(Group = anti_group, ratio = anti_ratio_time)
dat_anti_ratio_time$Group <- factor(dat_anti_ratio_time$Group, levels = c("Ctrl", "Anti", "Anti_pain"))
# dat_anti_1st_speed_sta <- ddply(dat_anti_1st_speed, .(Group), summarise,n=length(speed),mean=mean(speed),sd=sd(speed),se=sd(speed)/sqrt(length(speed)))

p_anti_ratio_time <- ggplot(dat_anti_ratio_time, aes(Group, ratio, colour=Group))+
  geom_boxplot(outlier.shape =NA )+
  geom_jitter(width = 0.25, shape=1)+
  labs(x="", y="Time spend in 2nd plate (%)")+
  theme(axis.line = element_line(colour = "black"),
        axis.text.x = element_text(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  theme(legend.position = 'none')

## combine all figures
p_anti_com <- plot_grid(p_anti_1st_latency, p_anti_num_cross, p_anti_2nd_latency, 
                        p_anti_1st_speed, p_anti_ratio_time, p_anti_2nd_speed, nrow = 2)

setwd("~cchen2/Documents/neuroscience/Pn\ project/Figure/PDF/")
cairo_pdf("p_anti_com.pdf", width = 140/25.6, height = 120/25.6, family = "Arial")
p_anti_com
dev.off()



## compare ctrl and opto, combine all data and plot-----
mouse_opto <- list(dat_behavior_m19, dat_behavior_m20, dat_behavior_m23, dat_behavior_m63, dat_behavior_m851)
mouse_wt <- list(dat_behavior_m16, dat_behavior_m17, dat_behavior_m18, dat_behavior_m22, dat_behavior_m3, dat_behavior_m7, dat_behavior_m27)

## for time spend in nor plate
prop_wt <- sapply(mouse_wt, function(x) x[[7]][[1]])
prop_opto <- sapply(mouse_opto, function(x) x[[7]][[1]])
prop_comb <- c(prop_wt, prop_opto)
prop_group <- rep(c("Ctrl", "Opto"), c(length(prop_wt), length(prop_opto)))

dat_anti_prop <- data.frame(Group = prop_group, Value = prop_comb)
dat_anti_prop$Group <- factor(dat_anti_prop$Group, levels = c("Ctrl", "Opto"))

p_anti_prop <- ggplot(dat_anti_prop, aes(Group, Value, colour=Group))+
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(width = 0.2, shape=1)+
  labs(x="", y="Time spend in 2nd plate (%)")+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  theme(legend.position = 'none')


## for 1st latency to nor chamber
latency_1st_wt <- sapply(mouse_wt, function(x) x[[7]][[2]])
latency_1st_opto <- sapply(mouse_opto, function(x) x[[7]][[2]])
latency_1st_comb <- c(latency_1st_wt, latency_1st_opto)
latency_1st_group <- rep(c("Ctrl", "Opto"), c(length(latency_1st_wt), length(latency_1st_opto)))

dat_anti_latency_1st <- data.frame(Group = latency_1st_group, Latency = latency_1st_comb)
dat_anti_latency_1st$Group <- factor(dat_anti_latency_1st$Group, levels = c("Ctrl", "Opto"))

p_anti_latency_1st <- ggplot(dat_anti_latency_1st, aes(Group, Latency, colour=Group))+
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(width = 0.2, shape=1)+
  labs(x="", y="Latency of 1st border crossing")+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  theme(legend.position = 'none')

## For the first latency to nor
latency_2nd_wt <- sapply(mouse_wt, function(x) x[[7]][[3]])
latency_2nd_opto <- sapply(mouse_opto, function(x) x[[7]][[3]])
latency_2nd_comb <- c(latency_2nd_wt, latency_2nd_opto)
latency_2nd_group <- rep(c("Ctrl", "Opto"), c(length(latency_2nd_wt), length(latency_2nd_opto)))

dat_anti_latency_2nd <- data.frame(Group = latency_2nd_group, Latency = latency_2nd_comb)
dat_anti_latency_2nd$Group <- factor(dat_anti_latency_2nd$Group, levels = c("Ctrl", "Opto"))

p_anti_latency_2nd <- ggplot(dat_anti_latency_2nd, aes(Group, Latency, colour=Group))+
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(width = 0.2, shape=1)+
  labs(x="", y="Latency of crossing back")+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  theme(legend.position = 'none')

## for 1st speed to nor chamber
speed_1st_wt <- sapply(mouse_wt, function(x) x[[7]][[5]])
speed_1st_opto <- sapply(mouse_opto, function(x) x[[7]][[5]])
speed_1st_comb <- c(speed_1st_wt, speed_1st_opto)
speed_1st_group <- rep(c("Ctrl", "Opto"), c(length(speed_1st_wt), length(speed_1st_opto)))

dat_anti_speed_1st <- data.frame(Group = speed_1st_group, speed = speed_1st_comb)
dat_anti_speed_1st$Group <- factor(dat_anti_speed_1st$Group, levels = c("Ctrl", "Opto"))

p_anti_speed_1st <- ggplot(dat_anti_speed_1st, aes(Group, speed, colour=Group))+
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(width = 0.2, shape=1)+
  labs(x="", y="speed of 1st border crossing")+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  theme(legend.position = 'none')


## Speed of crossing back
speed_2nd_wt <- sapply(mouse_wt, function(x) x[[7]][[6]])
speed_2nd_opto <- sapply(mouse_opto, function(x) x[[7]][[6]])
speed_2nd_comb <- c(speed_2nd_wt, speed_2nd_opto)
speed_2nd_group <- rep(c("Ctrl", "Opto"), c(length(speed_2nd_wt), length(speed_2nd_opto)))

dat_anti_speed_2nd <- data.frame(Group = speed_2nd_group, speed = speed_2nd_comb)
dat_anti_speed_2nd$Group <- factor(dat_anti_speed_2nd$Group, levels = c("Ctrl", "Opto"))

p_anti_speed_2nd <- ggplot(dat_anti_speed_2nd, aes(Group, speed, colour=Group))+
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(width = 0.2, shape=1)+
  labs(x="", y="speed of crossing back")+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  theme(legend.position = 'none')


## total crossing time
n_cross_wt <- sapply(mouse_wt, function(x) x[[7]][[7]])
n_cross_opto <- sapply(mouse_opto, function(x) x[[7]][[7]])
n_cross_comb <- c(n_cross_wt, n_cross_opto)
n_cross_group <- rep(c("Ctrl", "Opto"), c(length(n_cross_wt), length(n_cross_opto)))

dat_anti_n_cross <- data.frame(Group = n_cross_group, Value = n_cross_comb)
dat_anti_n_cross$Group <- factor(dat_anti_n_cross$Group, levels = c("Ctrl", "Opto"))

p_anti_n_cross <- ggplot(dat_anti_n_cross, aes(Group, Value, colour=Group))+
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(width = 0.2, shape=1)+
  labs(x="", y="Time spend in 2nd plate (%)")+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  theme(legend.position = 'none')

p_opto_com <- plot_grid(p_anti_latency_1st, p_anti_prop, p_anti_latency_2nd,p_anti_speed_2nd, nrow = 2)

setwd("~cchen2/Documents/neuroscience/Pn\ project/Figure/PDF/")
cairo_pdf("p_opto_com.pdf", width = 80/25.6, height = 120/25.6, family = "Arial")
p_opto_com
dev.off()

## for HM4Di behavior------
# for m10
setwd("~cchen2/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/anti_hm4di/m10")

path_trace_m10 <- list.files( pattern = "*.csv" )
dat_behavior_m10 <- mapply(cc_anti_behavior, path_trace_m10, frame_rate = 10, length_pix=455, SIMPLIFY = F)

## for m11
setwd("~cchen2/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/anti_hm4di/m11")

path_trace_m11 <- list.files( pattern = "*.csv" )
dat_behavior_m11 <- mapply(cc_anti_behavior, path_trace_m11, frame_rate = 10, length_pix=455, SIMPLIFY = F)

## for m12
setwd("~cchen2/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/anti_hm4di/m12")

path_trace_m12 <- list.files( pattern = "*.csv" )
dat_behavior_m12 <- mapply(cc_anti_behavior, path_trace_m12, frame_rate = 10, length_pix=455, SIMPLIFY = F)

# for m13
setwd("~cchen2/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/anti_hm4di/m13")

path_trace_m13 <- list.files( pattern = "*.csv" )
dat_behavior_m13 <- mapply(cc_anti_behavior, path_trace_m13, frame_rate = 10, length_pix=455, SIMPLIFY = F)

## combine and compare
mouse_HM4di <- list(dat_behavior_m10, dat_behavior_m11)
mouse_ctrl <- list(dat_behavior_m12, dat_behavior_m13)

## for time spend in nor plate
prop_ctrl <- sapply(mouse_ctrl, function(x) x[[7]][[1]])
prop_HM4di <- sapply(mouse_HM4di, function(x) x[[7]][[1]])
prop_comb <- c(prop_ctrl, prop_HM4di)
prop_group <- rep(c("Ctrl", "HM4di"), c(length(prop_ctrl), length(prop_HM4di)))

dat_anti_prop <- data.frame(Group = prop_group, Value = prop_comb)
dat_anti_prop$Group <- factor(dat_anti_prop$Group, levels = c("Ctrl", "HM4di"))

p_anti_prop <- ggplot(dat_anti_prop, aes(Group, Value, colour=Group))+
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(width = 0.2, shape=1)+
  labs(x="", y="Time spend in 2nd plate (%)")+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  theme(legend.position = 'none')


## for WT ctrl to prove the placebo analgesia-----
dat_anti_ctrl <- vector(mode = "list", 7)
dat_anti_con <- vector(mode = "list", 7)
length_pix <- c(363,362,363,364,362,364,365)
day_group <- c("Rm","Pre", "Pre", "Rm","Cond.","Cond.", "Test")
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
    add_column(Day = day_group[i], Group="Con", ID=str_c("m",1:10), .before = 'ratio')
}

dat_anti_wt <- do.call(rbind, dat_anti_ctrl) %>% 
  rbind(., do.call(rbind, dat_anti_con)) %>% 
  pivot_longer(-c(Day, Group, ID), names_to = "variable", values_to = "value" , values_drop_na = TRUE) %>% 
  ddply(., .(ID, variable, Day, Group), summarise, 'value'= mean(value)) %>% 
  filter(Day!="Rm") %>% 
  mutate(Day=factor(Day, levels = c("Pre", "Cond.", "Test")), Group=factor(Group, levels = c("Ctrl", "Con")))

dat_anti_wt_sta <- ddply(dat_anti_wt,.(variable,Day, Group), summarise, 'value'= mean(value))

p_ratio <- dat_anti_wt %>% 
  filter(variable=="ratio" & Day!="Cond.") %>% 
  mutate(value=value*100) %>% 
  ggplot(., aes(Day, value, fill=Group))+
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(position=position_jitterdodge(), shape=1)+
  scale_fill_manual(values=c("seagreen", "indianred"))+
  labs(x="", y="Time spend in 2nd plate (%)")+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  scale_y_continuous(limits = c(0, 100), expand = c(0, 0))+
  theme(legend.position = "none")+
  annotate(x=c(1.8,1.8,2.2,2.2), y=c(85,87,87,85)+5,"path")+
  annotate("text",x=2,y=92, label="***", size=5)

p_ratio_test <- dat_anti_wt %>% 
  filter(variable=="ratio" & Day=="Test") %>% 
  wilcox.test(value~Group, .)

p_latency1 <- dat_anti_wt %>% 
  filter(variable=="latency1"& Day!="Cond.") %>% 
  ggplot(., aes(Day, value, fill=Group))+
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(position=position_jitterdodge(), shape=1)+
  scale_fill_manual(values=c("seagreen", "indianred"))+
  labs(x="", y="Latency (s)")+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
 scale_y_continuous(limits = c(0, 150), expand = c(0, 0))+
  theme(legend.position = "none")+
  annotate(x=c(1.8,1.8,2.2,2.2), y=c(100,102,102,100),"path")+
  annotate("text",x=2,y=102, label="*", size=5)

p_latency1_test <- dat_anti_wt %>% 
  filter(variable=="latency1"& Day=="Test") %>% 
  wilcox.test(value~Group, .)

p_latency2 <- dat_anti_wt %>% 
  filter(variable=="latency2"& Day!="Cond.") %>% 
  ggplot(., aes(Day, value, fill=Group))+
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(position=position_jitterdodge(), shape=1)+
  scale_fill_manual(values=c("seagreen", "indianred"))+
  labs(x="", y="Latency (s)")+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  scale_y_continuous(limits = c(5, 80))+
  theme(legend.position = "none")+
  annotate(x=c(1.8,1.8,2.2,2.2), y=c(70,72,72,70),"path")+
  annotate("text",x=2,y=72, label="**", size=5)

p_latency2_test <- dat_anti_wt %>% 
  filter(variable=="latency2"& Day=="Test") %>% 
  wilcox.test(value~Group,.)

p_speed_compare <- dat_anti_wt %>% 
  filter(variable=="speed_compare"& Day!="Cond.") %>% 
  ggplot(., aes(Day, value, fill=Group))+
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(position=position_jitterdodge(), shape=1)+
  scale_fill_manual(values=c("seagreen", "indianred"))+
  labs(x="", y="Ratio of moving speed")+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  theme(legend.position = "none")+
  annotate(x=c(1.8,1.8,2.2,2.2), y=c(2,2.02,2.02,2),"path")+
  annotate("text",x=2,y=2.02, label="***", size=5)

p_speed_compare_test <- dat_anti_wt %>% 
  filter(variable=="speed_compare"& Day=="Test") %>% 
  wilcox.test(value~Group,.)

p_total_cross <- dat_anti_wt %>% 
  filter(variable=="total_cross"& Day!="Cond.") %>% 
  ggplot(., aes(Day, value, fill=Group))+
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(position=position_jitterdodge(), shape=1)+
  scale_fill_manual(values=c("seagreen", "indianred"))+
  labs(x="", y="Ratio of moving speed")+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  theme(legend.position = "none")+
  annotate(x=c(1.8,1.8,2.2,2.2), y=c(2,2.02,2.02,2),"path")+
  annotate("text",x=2,y=2.02, label="***", size=5)

p_total_cross_test <- dat_anti_wt %>% 
  filter(variable=="total_cross"& Day=="Test") %>% 
  wilcox.test(value~Group,.)

p_pain_compare1 <- plot_grid(p_ratio, p_latency1, p_latency2, p_speed_compare, ncol = 2)


## for pain behavir, like licking, and wearing-----
# mannualy analyzed by chong (07062020)
dat_wt_manual <- read.xlsx("~cchen2/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/wt_ctrl_07042020/wt_ctrl_07042020_manualy.xlsx") %>% 
  mutate('Rearing'= Total_rearing/Total_frame*600, 'First_rearing'=First_rearing*0.1, 
         'Guarding'= (Guarding - Frame_1st_crossing)*0.1, 'Acceleration'= (Acceleration -Frame_1st_crossing )*0.1) %>% 
  pivot_longer(-c(Group, ID), names_to = "variable", values_to = 'value') %>% 
  mutate(Group=factor(Group, levels=c("Ctrl", "Con")))

p_licking <- dat_wt_manual %>% 
  filter(variable=="First_licking") %>% 
  ggplot(., aes(Group, value, fill=Group))+
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(position=position_jitterdodge(), shape=1)+
  scale_fill_manual(values=c("seagreen", "indianred"))+
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
  ggplot(., aes(Group, value, fill=Group))+
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(position=position_jitterdodge(), shape=1)+
  scale_fill_manual(values=c("seagreen", "indianred"))+
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
  ggplot(., aes(Group, value, fill=Group))+
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(position=position_jitterdodge(), shape=1)+
  scale_fill_manual(values=c("seagreen", "indianred"))+
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

p_rearing_test<- dat_wt_manual %>% 
  filter(variable=="Rearing") %>% 
  wilcox.test(value~Group,.)

p_Guarding <- dat_wt_manual %>% 
  filter(variable=="Guarding") %>% 
  ggplot(., aes(Group, value, fill=Group))+
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(position=position_jitterdodge(), shape=1)+
  scale_fill_manual(values=c("seagreen", "indianred"))+
  labs(x="", y="Latency to 1st guarding (s) ")+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  scale_y_continuous(limits = c(2, 25))+
  theme(legend.position = "none")+
  annotate(x=c(1,1,2,2), y=c(23,24,24,23),"path")+
  annotate("text",x=1.5,y=24, label="*", size=5)

p_Guarding_test<- dat_wt_manual %>% 
  filter(variable=="Guarding") %>% 
  wilcox.test(value~Group,.)

p_Acceleration <- dat_wt_manual %>% 
  filter(variable=="Acceleration") %>% 
  ggplot(., aes(Group, value, fill=Group))+
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(position=position_jitterdodge(), shape=1)+
  scale_fill_manual(values=c("seagreen", "indianred"))+
  labs(x="", y="Latency to 1st acceleration (s) ")+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  theme(legend.position = "none")+
  annotate(x=c(1,1,2,2), y=c(75,77,77,75),"path")+
  annotate("text",x=1.5,y=77, label="***", size=5)

p_Acceleration_test<- dat_wt_manual %>% 
  filter(variable=="Acceleration") %>% 
  wilcox.test(value~Group,.)
p_pain_compare2 <- plot_grid(p_licking, p_1st_rearing, p_rearing, p_Acceleration, nrow = 1)
p_con <- plot_grid(p_pain_compare1, p_pain_compare2, nrow = 2, rel_heights = c(2, 1))

setwd("~cchen2/Documents/neuroscience/Pn\ project/Figure/PDF/")
cairo_pdf("p_con.pdf", width = 120/25.6, height = 180/25.6, family = "Arial")
p_con
dev.off()
## heat map plot to show the path----

cc_behavior_reshape <- function(path_ID, frame_rate, length_pix) {
  dat_anti <- read.csv(path_ID, header = F, stringsAsFactors = F)
  
  dat_anti <- dat_anti[-c(1:3), ]
  ## only slect the position of head
  dat_anti <- as.data.frame(apply(dat_anti, 2, as.numeric))
  
  ## normalize the data, zero point value from fiji
  dat_anti$V2<- dat_anti$V2 - dat_anti$V8
  dat_anti$V3<- dat_anti$V3 - dat_anti$V9
  dat_anti$V5 <- dat_anti$V5 - dat_anti$V8
  dat_anti$V6 <- dat_anti$V6 - dat_anti$V9
  dat_anti <- dat_anti[,1:6]
  
  colnames(dat_anti)<- c("Frame", "Head_x","Head_y", "Prob", "Tail_x", "Tail_y")
  
  ## remove the frame for the first few frame
  dat_anti <- dat_anti[seq(1, nrow(dat_anti), frame_rate/10), ]
  dat_anti <- subset(dat_anti, Prob >0.5 | Frame> 60)
  
  r_pix_mm <- 165 / length_pix
  dat_anti$Head_x <- abs(dat_anti$Head_x * r_pix_mm)
  dat_anti$Head_y <- abs(dat_anti$Head_y * r_pix_mm)
  
  dat_anti$Tail_x <- abs(dat_anti$Tail_x * r_pix_mm)
  dat_anti$Tail_y <- abs(dat_anti$Tail_y * r_pix_mm)
  
  chamber_div <- 165.5 ## the lenght of two chambels and the gap between
  
  ## assign each frame in hot as 1 (48) and nor plate(30) as 2
  dat_anti$chamber <- ifelse(dat_anti$Head_y < chamber_div, "Hot","Nor")
  return(dat_anti)
 
}

## for m1 day1
path_ctrl <- list(str_c("~cchen2/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/wt_ctrl_07042020/wt_ctrl/T",1),
                  str_c("~cchen2/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/wt_ctrl_07042020/wt_ctrl/T",5),
                  str_c("~cchen2/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/wt_ctrl_07042020/wt_ctrl/T",7))
cc_plot <-function(x){
  p_heat <- ggplot(x, aes(x=Head_y, y=Head_x) ) +
    geom_hex() +
    scale_fill_viridis()+
    theme_void()
  
  return(p_heat)
}
  
dat_heat_ctrl <- path_ctrl %>% 
  mapply(function(x) list.files(x, pattern = ".csv", full.names = T),., SIMPLIFY = F) %>% 
  mapply(function(x) x[1], ., SIMPLIFY = F) %>% 
  mapply(cc_behavior_reshape, ., frame_rate = 10, length_pix=364, SIMPLIFY = F) %>% 
  mapply(cc_plot,., SIMPLIFY = F) 

p_heat_ctrl <- plot_grid(dat_heat_ctrl[[1]], dat_heat_ctrl[[2]], dat_heat_ctrl[[3]], ncol = 1)
  
setwd("~cchen2/Documents/neuroscience/Pn\ project/Figure/PDF/")
cairo_pdf("p_heat_ctrl.pdf", width = 90/25.6, height = 120/25.6, family = "Arial")
p_heat_ctrl
dev.off()

path_con <- list(str_c("~cchen2/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/wt_ctrl_07042020/wt_con/T",1),
                  str_c("~cchen2/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/wt_ctrl_07042020/wt_con/T",5),
                  str_c("~cchen2/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/wt_ctrl_07042020/wt_con/T",7))


dat_heat_con <- path_con %>% 
  mapply(function(x) list.files(x, pattern = ".csv", full.names = T),., SIMPLIFY = F) %>% 
  mapply(function(x) x[1], ., SIMPLIFY = F) %>% 
  mapply(cc_behavior_reshape, ., frame_rate = 10, length_pix=364, SIMPLIFY = F) %>% 
  mapply(cc_plot,., SIMPLIFY = F) 

p_heat_con <- plot_grid(dat_heat_con[[1]], dat_heat_con[[2]], dat_heat_con[[3]], ncol = 1)

setwd("~cchen2/Documents/neuroscience/Pn\ project/Figure/PDF/")
cairo_pdf("p_heat_con.pdf", width = 90/25.6, height = 120/25.6, family = "Arial")
p_heat_con
dev.off()

path_con <-  str_c("~cchen2/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/wt_ctrl_07042020/wt_con/T",i)

