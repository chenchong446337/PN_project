

cc_dcb_tpp <- function(path_ID, frame_rate){
  
  dat_anti <- read.csv(path_ID, header = F, stringsAsFactors = F) %>% 
    .[-c(1:3), ] %>% 
    apply(., 2, as.numeric) %>%   ## only select the position of head
    as.data.frame() %>% 
    mutate(body_x = (V2 + V5)/2, body_y = (V3+ V6)/2, refer_y = mean(V12)) %>% 
    select(V1,body_x,body_y, refer_y)
  

  
  dat_anti <- dat_anti[seq(1, nrow(dat_anti), frame_rate/10), ] %>% 
    mutate(chamber = ifelse(body_y > refer_y, "left", "right"))
  
  ## for the right side
  ratio_time <- unname(by(dat_anti[,3], dat_anti$chamber, length)[1]/nrow(dat_anti))
  
  p_anti <- ggplot(dat_anti, aes(-body_y, body_x,colour=chamber))+
    geom_point(shape=1, size=1)+
    theme_void()+
    theme(legend.position = "none")
  
  return(ratio_time)
}

## for data from 30_30 TPP
dat_TPP_ctrl_30 <- list.files("~cchen2/Documents/neuroscience/Pn\ project/Data_analysis/Dor_dcb/TPP_ctrl/30_30",pattern = ".csv", full.names = T ) %>% 
  as.list() %>% 
  mapply(cc_dcb_tpp, ., frame_rate = 10)

dat_TPP_caspas3_30 <- list.files("~cchen2/Documents/neuroscience/Pn\ project/Data_analysis/Dor_dcb/TPP_caspas3/30_30",pattern = ".csv", full.names = T ) %>% 
  as.list() %>% 
  mapply(cc_dcb_tpp, ., frame_rate = 10)

## for 30_20 TPP----
dat_TPP_ctrl_20 <- list.files("~cchen2/Documents/neuroscience/Pn\ project/Data_analysis/Dor_dcb/TPP_ctrl/30_20",pattern = ".csv", full.names = T ) %>% 
  as.list() %>% 
  mapply(cc_dcb_tpp, ., frame_rate = 10)
  
dat_TPP_caspas3_20 <- list.files("~cchen2/Documents/neuroscience/Pn\ project/Data_analysis/Dor_dcb/TPP_caspas3/30_20",pattern = ".csv", full.names = T ) %>% 
  as.list() %>% 
  mapply(cc_dcb_tpp, ., frame_rate = 10)



## for data from 30_40 TPP
dat_TPP_ctrl_40 <- list.files("~cchen2/Documents/neuroscience/Pn\ project/Data_analysis/Dor_dcb/TPP_ctrl/30_40",pattern = ".csv", full.names = T ) %>% 
  as.list() %>% 
  mapply(cc_dcb_tpp, ., frame_rate = 10)

dat_TPP_caspas3_40 <- list.files("~cchen2/Documents/neuroscience/Pn\ project/Data_analysis/Dor_dcb/TPP_caspas3/30_40",pattern = ".csv", full.names = T ) %>% 
  as.list() %>% 
  mapply(cc_dcb_tpp, ., frame_rate = 10)

dat_TPP_ctrl <- c(dat_TPP_ctrl_20, dat_TPP_ctrl_30, dat_TPP_ctrl_40)
dat_TPP_caspas <- c(dat_TPP_caspas3_20, dat_TPP_caspas3_30, dat_TPP_caspas3_40)

dat_group <- c(rep("Ctrl", length(dat_TPP_ctrl)), rep("Caspas3", length(dat_TPP_caspas)))
dat_temp <- c(rep(c(20, 30, 40), each = length(dat_TPP_ctrl_20)), rep(c(20, 30, 40), each = length(dat_TPP_caspas3_20)))

p_ratio <- tibble(Group = dat_group, Ratio = c(dat_TPP_ctrl, dat_TPP_caspas), Temp = dat_temp ) %>% 
  mutate(Group = factor(Group, levels = c("Ctrl", "Caspas3"))) %>% 
  mutate(Temp = factor(Temp, levels = c(20, 30, 40))) %>% 
  mutate(Ratio = Ratio *100) %>% 
  ggplot(., aes(interaction(Group, Temp),Ratio, fill=Group))+
  geom_boxplot(outlier.shape = NA)+
  geom_point(aes(colour = Group,shape = Group),position=position_jitterdodge(jitter.width = .2))+
  scale_fill_manual(values = c("white", "white" ))+
  scale_shape_manual(values=c(16,17))+
  scale_color_manual(values=c("darkcyan", "indianred"))+
  labs(x="Test temperature", y="Time at test temperature (%)")+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  theme(legend.title = element_blank(), legend.position = c(0.2, 0.8))

setwd("~cchen2/Documents/neuroscience/Pn\ project/Figure/PDF/DCN")
cairo_pdf("p_ratio.pdf", width = 90/25.6, height = 60/25.6, family = "Arial")
p_ratio
dev.off()

## for data from 40_50 TPP
dat_TPP_ctrl_50 <- list.files("~cchen2/Documents/neuroscience/Pn\ project/Data_analysis/Dor_dcb/TPP_ctrl/40_50",pattern = ".csv", full.names = T ) %>% 
  as.list() %>% 
  mapply(cc_dcb_tpp, ., frame_rate = 10)

dat_TPP_caspas3_50 <- list.files("~cchen2/Documents/neuroscience/Pn\ project/Data_analysis/Dor_dcb/TPP_caspas3/40_50",pattern = ".csv", full.names = T ) %>% 
  as.list() %>% 
  mapply(cc_dcb_tpp, ., frame_rate = 10)

dat_tpp_50 <- c(dat_TPP_ctrl_50, dat_TPP_caspas3_50)
dat_group_50 <- c(rep("Ctrl", length(dat_TPP_ctrl_50)), rep("Caspas3", length(dat_TPP_caspas3_50)))

p_ratio_50 <- tibble(Group = dat_group_50, Ratio = dat_tpp_50) %>% 
  mutate(Group = factor(Group, levels = c("Ctrl", "Caspas3"))) %>% 
  mutate(Ratio = Ratio *100) %>% 
  ggplot(., aes(Group, Ratio,  fill = Group))+
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(aes(colour= Group,shape = Group),width = 0.2)+
  scale_fill_manual(values = c("white", "white" ))+
  scale_shape_manual(values=c(16,17))+
  scale_color_manual(values=c("darkcyan", "indianred"))+
  labs(x="Test temperature", y="Time at test temperature (%)")+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  theme(legend.title = element_blank(), legend.position = "none")

setwd("~cchen2/Documents/neuroscience/Pn\ project/Figure/PDF/DCN")
cairo_pdf("p_ratio_50.pdf", width = 40/25.6, height = 60/25.6, family = "Arial")
p_ratio_50
dev.off()


## for the thermal track data analysis-----

cc_dcb_track <- function(path_ID){
  frame_rate=24
  ID_mouse <- str_match(path_ID, "m\\d{3}")
  dat_anti <- read.csv(path_ID, header = F, stringsAsFactors = F) %>% 
    .[-c(1:3), ] %>% 
    apply(., 2, as.numeric) %>%   
    as.data.frame() %>% 
    slice(-(1: frame_rate)) %>% # drop the first second video
    #mutate(body_x = (V2 + V5)/2, body_y = (V3+ V6)/2) %>% 
    rename(body_x = V2, body_y = V3) %>% 
    select(body_x,body_y, V8, V11) %>% 
    mutate(body_x = ifelse(body_x < mean(V11), mean(V11), ifelse(body_x > mean(V8), mean(V8), body_x))) %>% 
    round(., 1) %>% 
    slice_head(n = frame_rate*60*30) # only analysis for 30 min
  
  tem_zone <- seq(5, 48, length.out = 25)
  track_zone <- round(seq(mean(dat_anti$V11), mean(dat_anti$V8), length.out = 26),1)
 
   dat_anti1 <- dat_anti[seq(1, nrow(dat_anti), frame_rate/10), ] %>% 
    mutate(Time = seq(1,10*60*30 )) %>% 
    group_by(zone=cut(body_x, breaks= track_zone )) %>% 
    mutate(zone = as.numeric(zone)) %>% 
    mutate(Zone_group = ifelse(zone < 7, "Cold", ifelse(zone>21, "Hot", "Nor")))
  
  ## calculate the time spend in hot and cold zone  
  dat_anti_duration <- dat_anti1 %>% 
    group_by(Time_group = cut(Time, breaks = seq(0, 18000, 300))) %>% 
    mutate(Time_group = as.numeric(Time_group)*0.5) %>% 
    ddply(., .(Time_group,zone,Zone_group), summarise, duration = length(body_x)*0.1) %>% 
    mutate(ID = ID_mouse )
  
  return(dat_anti_duration)
}

## for the EGFP group

dat_track_EGFP <- list.files("~cchen2/Documents/neuroscience/Pn\ project/Data_analysis/Dor_dcb/Thermal_Track/EGFP", pattern = ".csv", full.names = T) %>% 
  as.list() %>% 
  mapply(cc_dcb_track, ., SIMPLIFY = F) %>% 
  do.call(rbind,.)

## for the caspas3 group

dat_track_Caspas3 <- list.files("~cchen2/Documents/neuroscience/Pn\ project/Data_analysis/Dor_dcb/Thermal_Track/Caspas3", pattern = ".csv", full.names = T) %>% 
  as.list() %>% 
  mapply(cc_dcb_track, ., SIMPLIFY = F) %>% 
  do.call(rbind,.)

dat_track <- rbind(dat_track_EGFP, dat_track_Caspas3) %>% 
  mutate(Group = c(rep("EGFP", nrow(dat_track_EGFP)), rep("Caspas3", nrow(dat_track_Caspas3))))

## combine 
p_zone <- dat_track %>% 
  ddply(., .(ID,zone, Group), summarise, stay = sum(duration)) %>% 
  ddply(., .(zone, Group), summarise,n=length(stay),mean=mean(stay),sd=sd(stay),se=sd(stay)/sqrt(length(stay))) %>% 
  mutate(Group=factor(Group, levels = c("EGFP", "Caspas3"))) %>% 
  ggplot(., aes(zone, mean, colour = Group))+
  geom_errorbar(aes(ymin=mean-se-se, ymax=mean+se), width=0.1)+
  geom_point()+
  geom_line()

   

p_heat_cold <- dat_track %>% 
  filter(Zone_group == "Cold" | Zone_group =="Hot") %>% 
  ddply(., .(Group, ID), summarise, latency = sum(duration)) %>% 
  mutate(Group=factor(Group, levels = c("EGFP", "Caspas3"))) %>% 
  ggplot(., aes(Group, latency, colour= Group))+
  geom_boxplot()+
  geom_jitter()

p_heat_cold_cum <- dat_track %>% 
  filter(Zone_group == "Cold" | Zone_group =="Hot") %>% 
  ddply(., .(ID, Time_group,Group), summarise, latency = sum(duration)) %>% 
  group_by(ID) %>% 
  mutate(cum = cumsum(latency)) %>% 
  mutate(Group=factor(Group, levels = c("EGFP", "Caspas3"))) %>%
  ddply(., .( Time_group,Group), summarise, mean = mean(cum)) %>% 
  ggplot(., aes(Time_group, mean, colour= Group))+
  geom_line()



## for the hot-cold test----
p_cold <- read.xlsx("~cchen2/Documents/neuroscience/Pn\ project/Data_analysis/Dor_dcb/DOR_DCB_12152020.xlsx", sheet = "Heat_cold") %>% 
  filter(Temp=="Cold") %>% 
  ggplot(., aes(Group, Withdraw, fill=Group))+
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(aes(colour= Group,shape = Group),width = 0.2)+
  scale_fill_manual(values = c("white", "white" ))+
  scale_shape_manual(values=c(16,17))+
  scale_color_manual(values=c("darkcyan", "indianred"))+
  labs(x="", y="Withdraw latency (s)")+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  theme(legend.title = element_blank(), legend.position = "none")


p_hot <- read.xlsx("~cchen2/Documents/neuroscience/Pn\ project/Data_analysis/Dor_dcb/DOR_DCB_12152020.xlsx", sheet = "Heat_cold") %>% 
  filter(Temp=="Hot") %>% 
  ggplot(., aes(Group, Withdraw, fill=Group))+
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(aes(colour= Group,shape = Group),width = 0.2)+
  scale_fill_manual(values = c("white", "white" ))+
  scale_shape_manual(values=c(16,17))+
  scale_color_manual(values=c("darkcyan", "indianred"))+
  labs(x="", y="Withdraw latency (s)")+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  theme(legend.title = element_blank(), legend.position = "none")

p_cold_hot <- plot_grid(p_hot, p_cold, ncol = 2)

setwd("~cchen2/Documents/neuroscience/Pn\ project/Figure/PDF/DCN")
cairo_pdf("p_cold_hot.pdf", width = 90/25.6, height = 60/25.6, family = "Arial")
p_cold_hot
dev.off()

## for hotplate test------
p_hotplate <- read.xlsx("~cchen2/Documents/neuroscience/Pn\ project/Data_analysis/Dor_dcb/DOR_DCB_12152020.xlsx", sheet = "Hotplate") %>% 
  mutate(Group = factor(Group, levels = c("EGFP", "Caspas3"))) %>% 
  ggplot(., aes(Group, Latency_jump, fill=Group))+
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(aes(colour= Group,shape = Group),width = 0.2)+
  scale_fill_manual(values = c("white", "white" ))+
  scale_shape_manual(values=c(16,17))+
  scale_color_manual(values=c("skyblue4", "tomato3"))+
  labs(x="", y="Latency to 1st jumping (s)")+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  theme(legend.title = element_blank(), legend.position = "none")

setwd("~cchen2/Documents/neuroscience/Pn\ project/Figure/PDF/DCN")
cairo_pdf("p_hotplate.pdf", width = 42/25.6, height = 60/25.6, family = "Arial")
p_hotplate
dev.off()

## for hotplate test------
p_vonfrey <- read.xlsx("~cchen2/Documents/neuroscience/Pn\ project/Data_analysis/Dor_dcb/DOR_DCB_12152020.xlsx", sheet = "Vonfrey") %>% 
  mutate(Group = factor(Group, levels = c("EGFP", "Caspas3"))) %>% 
  ggplot(., aes(Group, Threshold, fill=Group))+
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(aes(colour= Group,shape = Group),width = 0.2)+
  scale_fill_manual(values = c("white", "white" ))+
  scale_shape_manual(values=c(16,17))+
  scale_color_manual(values=c("skyblue4", "tomato3"))+
  labs(x="", y="pain threshold (g)")+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  theme(legend.title = element_blank(), legend.position = "none")

setwd("~cchen2/Documents/neuroscience/Pn\ project/Figure/PDF/DCN")
cairo_pdf("p_vonfrey.pdf", width = 42/25.6, height = 60/25.6, family = "Arial")
p_vonfrey
dev.off()

p_von_mutiple <- read.xlsx("~cchen2/Documents/neuroscience/Pn\ project/Data_analysis/Dor_dcb/DOR_DCB_12152020.xlsx", sheet = "VonFrey2") %>% 
  mutate(Group = factor(Group, levels = c("EGFP", "Caspas3"))) %>% 
  ddply(., .(Group, Num_von), summarise,mean=mean(Withdraw),sd=sd(Withdraw),se=sd(Withdraw)/sqrt(length(Withdraw))) %>% 
  ggplot(., aes(Num_von, mean, color=Group))+
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=0.1) +
  geom_point(aes(colour = Group, shape = Group),size=2)+
  geom_line()+
  scale_shape_manual(values=c(16,17))+
  scale_color_manual(values=c("skyblue4", "tomato3"))+
  labs(x="Force (g)", y="Withdraw frequency")+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  scale_x_continuous(trans='log10')+
  theme(legend.title = element_blank(), legend.position = c(0.2, 0.8))

setwd("~cchen2/Documents/neuroscience/Pn\ project/Figure/PDF/DCN")
cairo_pdf("p_von_multiple.pdf", width = 90/25.6, height = 60/25.6, family = "Arial")
p_von_mutiple
dev.off()
## for splash test----
p_splash <- read.xlsx("~cchen2/Documents/neuroscience/Pn\ project/Data_analysis/Dor_dcb/DOR_DCB_12152020.xlsx", sheet = "Splash_test") %>% 
  mutate(Group = factor(Group, levels = c("EGFP", "Caspas3"))) %>% 
  ggplot(., aes(Group, Duration, fill=Group))+
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(aes(colour= Group,shape = Group), width = 0.25)+
  scale_fill_manual(values = c("white", "white" ))+
  scale_shape_manual(values=c(16,17))+
  scale_color_manual(values=c("skyblue4", "tomato3"))+
  labs(x="", y="Time of grooming (s)")+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  theme(legend.title = element_blank(), legend.position = "none")

setwd("~cchen2/Documents/neuroscience/Pn\ project/Figure/PDF/DCN")
cairo_pdf("p_splash.pdf", width = 40/25.6, height = 60/25.6, family = "Arial")
p_splash
dev.off()
## Raster plot the splashing behavior data
p_EGFP<- read.xlsx("~cchen2/Documents/neuroscience/Pn\ project/Data_analysis/Dor_dcb/DOR_DCB_12152020.xlsx", sheet = "Spalsh_test2") %>% 
  mutate(Type = factor(Type, levels = c(1,2))) %>% 
  mutate(Group = factor(Group, levels = c("EGFP", "Caspas3"))) %>% 
  filter(Group == "EGFP") %>% 
  ggplot(., aes(Bin, ID, fill=Type, colour=Type))+
  geom_tile()+
  scale_colour_manual(values=c( "white", "skyblue4" ))+
  scale_fill_manual(values=c("white", "skyblue4"  ))+
  labs(x = "", y = "Mouse ID")+
  theme_classic()+
  theme(legend.position = "none")

p_caspas3<- read.xlsx("~cchen2/Documents/neuroscience/Pn\ project/Data_analysis/Dor_dcb/DOR_DCB_12152020.xlsx", sheet = "Spalsh_test2") %>% 
  mutate(Type = factor(Type, levels = c(1,2))) %>% 
  mutate(Group = factor(Group, levels = c("EGFP", "Caspas3"))) %>% 
  filter(Group == "Caspas3") %>% 
  ggplot(., aes(Bin, ID, fill=Type, colour=Type))+
  geom_tile()+
  scale_colour_manual(values=c( "white", "tomato3" ))+
  scale_fill_manual(values=c("white", "tomato3"  ))+
  labs(x = "Time (s)", y = "Mouse ID")+
  theme_classic()+
  theme(legend.position = "none")

p_splash_raster <- plot_grid(p_EGFP, p_caspas3, nrow = 2)

setwd("~cchen2/Documents/neuroscience/Pn\ project/Figure/PDF/")
cairo_pdf("p_splash_raster.pdf", width = 95/25.6, height = 65/25.6, family = "Arial")
p_splash_raster
dev.off()
## for the tail suspension test----
p_suspension<- read.xlsx("~cchen2/Documents/neuroscience/Pn\ project/Data_analysis/Dor_dcb/DOR_DCB_12152020.xlsx", sheet = "Tail_suspension") %>% 
  mutate(Group = factor(Group, levels = c("EGFP", "Caspas3"))) %>% 
  ggplot(., aes(Group, Duration, fill=Group))+
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(aes(colour= Group,shape = Group), width = 0.25)+
  scale_fill_manual(values = c("white", "white" ))+
  scale_shape_manual(values=c(16,17))+
  scale_color_manual(values=c("skyblue4", "tomato3"))+
  labs(x="", y="pain threshold (g)")+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  theme(legend.title = element_blank(), legend.position = "none")

setwd("~cchen2/Documents/neuroscience/Pn\ project/Figure/PDF/DCN")
cairo_pdf("p_suspension.pdf", width = 40/25.6, height = 60/25.6, family = "Arial")
p_suspension
dev.off()

## for the forced swimming test----
p_FST_freq<- read.xlsx("~cchen2/Documents/neuroscience/Pn\ project/Data_analysis/Dor_dcb/DOR_DCB_12152020.xlsx", sheet = "DOR_CRE_FST") %>% 
  mutate(Group = factor(Group, levels = c("EGFP", "Caspas3"))) %>% 
  ggplot(., aes(Group, Freq_immobility, fill=Group))+
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(aes(colour= Group,shape = Group), width = 0.25)+
  scale_fill_manual(values = c("white", "white" ))+
  scale_shape_manual(values=c(16,17))+
  scale_color_manual(values=c("skyblue4", "tomato3"))+
  labs(x="", y="Immobility frequency (#/min)")+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  theme(legend.title = element_blank(), legend.position = "none")

p_FST_duration<- read.xlsx("~cchen2/Documents/neuroscience/Pn\ project/Data_analysis/Dor_dcb/DOR_DCB_12152020.xlsx", sheet = "DOR_CRE_FST") %>% 
  mutate(Group = factor(Group, levels = c("EGFP", "Caspas3"))) %>% 
  ggplot(., aes(Group, Duration_immobility, fill=Group))+
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(aes(colour= Group,shape = Group), width = 0.25)+
  scale_fill_manual(values = c("white", "white" ))+
  scale_shape_manual(values=c(16,17))+
  scale_color_manual(values=c("skyblue4", "tomato3"))+
  labs(x="", y="Immobility duration (s)")+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  theme(legend.title = element_blank(), legend.position = "none")

p_fst <- plot_grid(p_FST_freq, p_FST_duration, ncol = 2)
setwd("~cchen2/Documents/neuroscience/Pn\ project/Figure/PDF/DCN")
cairo_pdf("p_fst.pdf", width = 80/25.6, height = 60/25.6, family = "Arial")
p_fst
dev.off()

## DOR-cre the open field experiment-------
dat_open <- read.xlsx("~cchen2/Documents/neuroscience/Pn\ project/Data_analysis/Dor_dcb/DCB_caspas_01212021.xlsx", sheet = 1) %>% 
  mutate(Group = factor(Group, levels = c("EGFP", "Caspas")))

p_oepn_distance <- dat_open %>% 
  select(ID, Group, Bin, Distance) %>% 
  ddply(., .(Group, Bin), summarise,n=length(Distance),mean=mean(Distance),sd=sd(Distance),se=sd(Distance)/sqrt(length(Distance))) %>% 
  ggplot(., aes(Bin, mean, colour = Group,group = Group)) +
  geom_point(aes(shape = Group), size=2)+
  geom_line() +
  geom_errorbar(aes(ymin=mean-se,ymax=mean+se), width=.2)+
  scale_shape_manual(values=c(16,17))+
  scale_color_manual(values=c("skyblue4", "tomato3"))+
  labs(x="Time (min)", y="Travel distance (cm)")+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  theme(legend.title = element_blank(), legend.position = c)

p_oepn_velocity <- dat_open %>% 
  select(ID, Group, Bin, Velocity) %>% 
  ddply(., .(Group, Bin), summarise,n=length(Velocity),mean=mean(Velocity),sd=sd(Velocity),se=sd(Velocity)/sqrt(length(Velocity))) %>% 
  ggplot(., aes(Bin, mean, color= Group, group = Group)) +
  geom_point(aes(shape = Group), size=2)+
  geom_line() +
  geom_errorbar(aes(ymin=mean-se,ymax=mean+se), width=.2)+
  scale_shape_manual(values=c(16,17))+
  scale_color_manual(values=c("skyblue4", "tomato3"))+
  labs(x="Time (min)", y="Velocity (cm/s)")+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  theme(legend.title = element_blank(), legend.position = c(0.8, 0.8))

setwd("~cchen2/Documents/neuroscience/Pn\ project/Figure/PDF/DCN")
cairo_pdf("p_open_velocity.pdf", width = 90/25.6, height = 60/25.6, family = "Arial")
p_oepn_velocity
dev.off()

p_oepn_center <- dat_open %>% 
  select(ID, Group, Bin, Center_freq) %>% 
  ddply(., .(Group, Bin), summarise,n=length(Center_freq),mean=mean(Center_freq),sd=sd(Center_freq),se=sd(Center_freq)/sqrt(length(Center_freq))) %>% 
  ggplot(., aes(Bin, mean, color= Group, group = Group)) +
  geom_point(size=2)+
  geom_line() +
  geom_errorbar(aes(ymin=mean-se,ymax=mean+se), width=.2)

p_oepn_center_dur <- dat_open %>% 
  select(ID, Group, Bin, Center_duration) %>% 
  ddply(., .(Group, Bin), summarise,n=length(Center_duration),mean=mean(Center_duration),sd=sd(Center_duration),se=sd(Center_duration)/sqrt(length(Center_duration))) %>% 
  ggplot(., aes(Bin, mean, color= Group, group = Group)) +
  geom_point(size=2)+
  geom_line() +
  geom_errorbar(aes(ymin=mean-se,ymax=mean+se), width=.2)


p_oepn_center_velocity <- dat_open %>% 
  select(ID, Group, Bin, Center_velocity) %>% 
  ddply(., .(Group, Bin), summarise,n=length(Center_velocity),mean=mean(Center_velocity),sd=sd(Center_velocity),se=sd(Center_velocity)/sqrt(length(Center_velocity))) %>% 
  ggplot(., aes(Bin, mean, color= Group, group = Group)) +
  geom_point(size=2)+
  geom_line() +
  geom_errorbar(aes(ymin=mean-se,ymax=mean+se), width=.2)


p_oepn_center_dur <- dat_open %>% 
  select(ID, Group, Bin, Center_freq) %>% 
  ddply(., .(ID, Group), summarise, duration = sum(Center_freq)) %>% 
  ggplot(., aes(Group, duration, color= Group)) +
  geom_boxplot()+
  geom_jitter()


## Time in the center zone
p_openfield_center <- read.xlsx("~cchen2/Documents/neuroscience/Pn\ project/Data_analysis/Dor_dcb/DOR_DCB_12152020.xlsx", sheet = "DOR_CRE_Openfield") %>% 
  mutate(Group = factor(Group, levels = c("EGFP", "Caspas3"))) %>% 
  ggplot(., aes(Group, Center_time, fill=Group))+
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(aes(colour= Group,shape = Group), width = 0.25)+
  scale_fill_manual(values = c("white", "white" ))+
  scale_shape_manual(values=c(16,17))+
  scale_color_manual(values=c("skyblue4", "tomato3"))+
  labs(x="", y="Center duration (s)")+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  theme(legend.title = element_blank(), legend.position = "none")

setwd("~cchen2/Documents/neuroscience/Pn\ project/Figure/PDF/DCN")
cairo_pdf("p_openfield_center.pdf", width = 40/25.6, height = 60/25.6, family = "Arial")
p_openfield_center
dev.off()

## for EPM test-----
p_EPM_open <- read.xlsx("~cchen2/Documents/neuroscience/Pn\ project/Data_analysis/Dor_dcb/DOR_DCB_12152020.xlsx", sheet = "DOR_CRE_EPM") %>% 
  mutate(Group = factor(Group, levels = c("EGFP", "Caspas3"))) %>% 
  ggplot(., aes(Group, Open_arm, fill=Group))+
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(aes(colour= Group,shape = Group), width = 0.25)+
  scale_fill_manual(values = c("white", "white" ))+
  scale_shape_manual(values=c(16,17))+
  scale_color_manual(values=c("skyblue4", "tomato3"))+
  labs(x="", y="Time in open arm (s)")+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  theme(legend.title = element_blank(), legend.position = "none")

setwd("~cchen2/Documents/neuroscience/Pn\ project/Figure/PDF/DCN")
cairo_pdf("p_epm_open.pdf", width = 40/25.6, height = 60/25.6, family = "Arial")
p_EPM_open
dev.off()

## Rotarod experiment of mice with Caspas3 injection-----
p_rotarod <- read.xlsx("~cchen2/Documents/neuroscience/Pn\ project/Data_analysis/Dor_dcb/DOR_DCB_12152020.xlsx", sheet = "Rotarod") %>% 
  mutate(Group = factor(Group, levels = c("EGFP", "Caspas3"))) %>% 
  mutate(Latency = rowMeans(select(.,starts_with("L")), na.rm = T)) %>% 
  ggplot(., aes(Group, Latency, fill = Group))+
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(aes(colour= Group,shape = Group),width = 0.2)+
  scale_fill_manual(values = c("white", "white" ))+
  scale_shape_manual(values=c(16,17))+
  scale_color_manual(values=c("darkcyan", "indianred"))+
  labs(x="", y="Latency to fall (s)")+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  theme(legend.title = element_blank(), legend.position = "none")

setwd("~cchen2/Documents/neuroscience/Pn\ project/Figure/PDF/DCN")
cairo_pdf("p_latency_fall.pdf", width = 40/25.6, height = 60/25.6, family = "Arial")
p_rotarod
dev.off()


## openfiled of Dor flox mice with cre or egfp injection-----
p_open_center <- read.xlsx("~cchen2/Documents/neuroscience/Pn\ project/Data_analysis/Dor_dcb/DOR_DCB_12152020.xlsx", sheet = "Dor_flox_openfield") %>% 
  mutate(Group = factor(Group, levels = c("eGFP", "Cre"))) %>% 
  ggplot(., aes(Group, Center, fill = Group))+
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(aes(colour= Group,shape = Group),width = 0.2)+
  scale_fill_manual(values = c("white", "white" ))+
  scale_shape_manual(values=c(16,17))+
  scale_color_manual(values=c("seagreen4", "deepskyblue4"))+
  labs(x="", y="Center time (s)")+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  theme(legend.title = element_blank(), legend.position = "none")

p_open_distance <- read.xlsx("~cchen2/Documents/neuroscience/Pn\ project/Data_analysis/Dor_dcb/DOR_DCB_12152020.xlsx", sheet = "Dor_flox_openfield") %>% 
  mutate(Group = factor(Group, levels = c("eGFP", "Cre"))) %>%
  mutate(Distance = Distance/100) %>% 
  ggplot(., aes(Group, Distance, fill = Group))+
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(aes(colour= Group,shape = Group),width = 0.2)+
  scale_fill_manual(values = c("white", "white" ))+
  scale_shape_manual(values=c(16,17))+
  scale_color_manual(values=c("seagreen4", "deepskyblue4"))+
  labs(x="", y="Distance moved (m)")+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  theme(legend.title = element_blank(), legend.position = "none")
t_open_center <- read.xlsx("~cchen2/Documents/neuroscience/Pn\ project/Data_analysis/Dor_dcb/DOR_DCB_12152020.xlsx", sheet = "Dor_flox_openfield") %>% 
  wilcox.test(Center~Group,.)

setwd("~cchen2/Documents/neuroscience/Pn\ project/Figure/PDF/DCN")
cairo_pdf("p_open_center.pdf", width = 40/25.6, height = 60/25.6, family = "Arial")
p_open_center
dev.off()

setwd("~cchen2/Documents/neuroscience/Pn\ project/Figure/PDF/DCN")
cairo_pdf("p_open_distance.pdf", width = 40/25.6, height = 60/25.6, family = "Arial")
p_open_distance
dev.off()

## Splash test of Dor flox mice with cre or egfp injection-----
p_splash <- read.xlsx("~cchen2/Documents/neuroscience/Pn\ project/Data_analysis/Dor_dcb/DOR_DCB_Flox_mice.xlsx", sheet = "Splash test") %>% 
  mutate(Group = factor(Group, levels = c("eGFP", "Cre"))) %>% 
  ggplot(., aes(Group, Grooming_time, fill = Group))+
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(aes(colour= Group,shape = Group),width = 0.2)+
  scale_fill_manual(values = c("white", "white" ))+
  scale_shape_manual(values=c(16,17))+
  scale_color_manual(values=c("seagreen4", "deepskyblue4"))+
  labs(x="", y="Grooming time (s)")+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  theme(legend.title = element_blank(), legend.position = "none")

t_splash_dor_flox <- read.xlsx("~cchen2/Documents/neuroscience/Pn\ project/Data_analysis/Dor_dcb/DOR_DCB_Flox_mice.xlsx", sheet = "Splash test") %>% 
  mutate(Group = factor(Group, levels = c("eGFP", "Cre"))) %>% 
  wilcox.test(Grooming_time~Group,.)

setwd("~cchen2/Documents/neuroscience/Pn\ project/Figure/PDF/DCN")
cairo_pdf("p_splash.pdf", width = 40/25.6, height = 60/25.6, family = "Arial")
p_splash
dev.off()

## plot to show less time spend in center
p_cre <- read.xlsx("~cchen2/Documents/neuroscience/Pn\ project/Data_analysis/Dor_dcb/Dor_flox_openfield/Heatmap_example.xlsx", sheet = 1) %>% 
  ggplot(., aes(X_center, Y_center))+
  # stat_density_2d(aes(fill = ..density..), geom = "raster", contour = FALSE) +
  geom_path(colour = "deepskyblue4")+
  theme_void()
  

p_EGFP<- read.xlsx("~cchen2/Documents/neuroscience/Pn\ project/Data_analysis/Dor_dcb/Dor_flox_openfield/Heatmap_example.xlsx", sheet = 2) %>% 
  ggplot(., aes(X_center, Y_center))+
  # stat_density_2d(aes(fill = ..density..), geom = "raster", contour = FALSE) +
  geom_path(colour = "seagreen4")+
  theme_void()

p_open_example <- plot_grid(p_EGFP, p_cre, ncol = 2)


setwd("~cchen2/Documents/neuroscience/Pn\ project/Figure/PDF/DCN")
cairo_pdf("p_dorflox_open.pdf", width = 80/25.6, height = 40/25.6, family = "Arial")
p_open_example
dev.off()

## DOR-flox mice, openfield velocity------

p_oepn_velocity <- read.xlsx("~cchen2/Documents/neuroscience/Pn\ project/Data_analysis/Dor_dcb/DOR_DCB_Flox_mice.xlsx", sheet = "Openfield") %>% 
  select(ID, Group, Bin, Velocity) %>% 
  ddply(., .(Group, Bin), summarise,n=length(Velocity),mean=mean(Velocity),sd=sd(Velocity),se=sd(Velocity)/sqrt(length(Velocity))) %>% 
  ggplot(., aes(Bin, mean, color= Group, group = Group)) +
  geom_point(aes(shape = Group), size=2)+
  geom_line() +
  geom_errorbar(aes(ymin=mean-se,ymax=mean+se), width=.2)+
  scale_shape_manual(values=c(16,17))+
  scale_color_manual(values=c("seagreen4", "deepskyblue4"))+
  labs(x="Time (min)", y="Velocity (cm/s)")+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  theme(legend.title = element_blank(), legend.position = c(0.8, 0.8))

setwd("~cchen2/Documents/neuroscience/Pn\ project/Figure/PDF/DCN")
cairo_pdf("p_dor_open_velocity.pdf", width = 90/25.6, height = 60/25.6, family = "Arial")
p_oepn_velocity
dev.off()


## for DORflox mice with CRS ------
# WT with CRS to prove depression
p_open_center_wt <- read.xlsx("~cchen2/Documents/neuroscience/Pn\ project/Data_analysis/Dor_dcb/DOR_flox_with_CRS.xlsx", sheet = "WT_openfield") %>% 
  mutate(Group = factor(Group, levels = c("Ctrl", "CRS"))) %>% 
  ggplot(., aes(Group, Center, fill = Group))+
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(aes(colour= Group,shape = Group),width = 0.2)+
  scale_fill_manual(values = c("white", "white" ))+
  scale_shape_manual(values=c(16,17))+
  scale_color_manual(values=c("seagreen4", "deepskyblue4"))+
  labs(x="", y="Center time (s)")+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  theme(legend.title = element_blank(), legend.position = "none")

## DORflox with CRS openfield
p_open_center <- read.xlsx("~cchen2/Documents/neuroscience/Pn\ project/Data_analysis/Dor_dcb/DOR_flox_with_CRS.xlsx", sheet = "Openfield") %>% 
  mutate(Group = factor(Group, levels = c("eGFP", "Cre"))) %>% 
  ggplot(., aes(Group, Center, fill = Group))+
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(aes(colour= Group,shape = Group),width = 0.2)+
  scale_fill_manual(values = c("white", "white" ))+
  scale_shape_manual(values=c(16,17))+
  scale_color_manual(values=c("seagreen4", "deepskyblue4"))+
  labs(x="", y="Center time (s)")+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  theme(legend.title = element_blank(), legend.position = "none")

## DORflox with CRS SPT-----
p_SPT <- read.xlsx("~cchen2/Documents/neuroscience/Pn\ project/Data_analysis/Dor_dcb/DOR_flox_with_CRS.xlsx", sheet = "SPT") %>% 
  mutate(Group = factor(Group, levels = c("eGFP", "Cre"))) %>% 
  ggplot(., aes(Group, score, fill = Group))+
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(aes(colour= Group,shape = Group),width = 0.2)+
  scale_fill_manual(values = c("white", "white" ))+
  scale_shape_manual(values=c(16,17))+
  scale_color_manual(values=c("seagreen4", "deepskyblue4"))+
  labs(x="", y="Sucrose prefernce (%)")+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  theme(legend.title = element_blank(), legend.position = "none")

## DORflox with CRS tst-----
p_TST <- read.xlsx("~cchen2/Documents/neuroscience/Pn\ project/Data_analysis/Dor_dcb/DOR_flox_with_CRS.xlsx", sheet = "TST") %>% 
  mutate(Group = factor(Group, levels = c("eGFP", "Cre"))) %>% 
  ggplot(., aes(Group, Immobility, fill = Group))+
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(aes(colour= Group,shape = Group),width = 0.2)+
  scale_fill_manual(values = c("white", "white" ))+
  scale_shape_manual(values=c(16,17))+
  scale_color_manual(values=c("seagreen4", "deepskyblue4"))+
  labs(x="", y="Immobility (s)")+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  theme(legend.title = element_blank(), legend.position = "none")

## DORflox with CRS FST-----
p_FST <- read.xlsx("~cchen2/Documents/neuroscience/Pn\ project/Data_analysis/Dor_dcb/DOR_flox_with_CRS.xlsx", sheet = "FST") %>% 
  mutate(Group = factor(Group, levels = c("eGFP", "Cre"))) %>% 
  ggplot(., aes(Group, Immobility, fill = Group))+
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(aes(colour= Group,shape = Group),width = 0.2)+
  scale_fill_manual(values = c("white", "white" ))+
  scale_shape_manual(values=c(16,17))+
  scale_color_manual(values=c("seagreen4", "deepskyblue4"))+
  labs(x="", y="Immobility (s)")+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  theme(legend.title = element_blank(), legend.position = "none")

p_crs <- plot_grid(p_open_center_wt, p_open_center, p_SPT, p_TST, p_FST, ncol = 3)


setwd("~cchen2/Documents/neuroscience/Pn\ project/Figure/PDF/DCN")
cairo_pdf("p_CRS.pdf", width = 135/25.6, height = 120/25.6, family = "Arial")
p_crs
dev.off()

## DOR-flox mice, rotarod test(mice with CRS)----
p_rotarod <- read.xlsx("~cchen2/Documents/neuroscience/Pn\ project/Data_analysis/Dor_dcb/DOR_DCB_Flox_mice.xlsx", sheet = "Rotarod_crs") %>% 
  mutate(Group = factor(Group, levels = c("eGFP", "Cre"))) %>% 
  ddply(., .(Group, Mouse), summarise,n=length(Latency),mean=mean(Latency)) %>% 
  ggplot(., aes(Group, mean, fill = Group))+
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(aes(colour= Group,shape = Group),width = 0.2)+
  scale_fill_manual(values = c("white", "white" ))+
  scale_shape_manual(values=c(16,17))+
  scale_color_manual(values=c("seagreen4", "deepskyblue4"))+
  labs(x="", y="Latency to fall (s)")+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  theme(legend.title = element_blank(), legend.position = "none")

## DOR-cre mice, rotarod test(mice with DREAD injection)----
p_rotarod <- read.xlsx("~cchen2/Documents/neuroscience/Pn\ project/Data_analysis/Dor_dcb/DOR_DCB_Flox_mice.xlsx", sheet = "Rotarod_dread") %>% 
  mutate(Group = factor(Group, levels = c("hM4Di", "hM3Dq"))) %>% 
  ddply(., .(Group, Mouse), summarise,n=length(Latency),mean=mean(Latency)) %>% 
  ggplot(., aes(Group, mean, fill = Group))+
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(aes(colour= Group,shape = Group),width = 0.2)+
  scale_fill_manual(values = c("white", "white" ))+
  scale_shape_manual(values=c(16,17))+
  scale_color_manual(values=c("seagreen4", "deepskyblue4"))+
  labs(x="", y="Latency to fall (s)")+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  scale_y_continuous(limits = c(0, 250))+
  theme(legend.title = element_blank(), legend.position = "none")


setwd("~cchen2/Documents/neuroscience/Pn\ project/Figure/PDF/DCN")
cairo_pdf("p_rotarod_dread.pdf", width = 40/25.6, height = 60/25.6, family = "Arial")
p_rotarod
dev.off()

## DOR-flox mice, SPT (mice without CRS)----
p_SPT <- read.xlsx("~cchen2/Documents/neuroscience/Pn\ project/Data_analysis/Dor_dcb/DOR_DCB_Flox_mice.xlsx", sheet = "SPT") %>% 
  mutate(Group = factor(Group, levels = c("eGFP", "Cre"))) %>% 
  filter(Volume_suc >0 | Volume_wat >0) %>% 
  ggplot(., aes(Group, Score, fill = Group))+
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(aes(colour= Group,shape = Group),width = 0.2)+
  scale_fill_manual(values = c("white", "white" ))+
  scale_shape_manual(values=c(16,17))+
  scale_color_manual(values=c("seagreen4", "deepskyblue4"))+
  labs(x="", y="Sucrose prefernce (%)")+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  theme(legend.title = element_blank(), legend.position = "none")

## Dorcre mcice with DREAD injection and local perfusion-----
dat_dor_dread_center <- read.xlsx("~cchen2/Documents/neuroscience/Pn\ project/Data_analysis/Dor_dcb/DORcre_DREAD_perfusion.xlsx", sheet = "Dread_perfusion_OF") %>% 
  as_tibble() %>% 
  filter(Quality =="P") %>% 
  mutate(Solution = factor(Solution, level = c("Saline", "CNO"))) %>% 
  ggplot(., aes(Solution, Time_center, colour = Solution))+
  geom_boxplot()+
  geom_jitter(width = 0.2,color = "darkgray")+
  scale_color_manual(values = c("#0099B4FF","#AD002AFF"))+
  labs(x="", y="Center time (s)")+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  theme(legend.title = element_blank(), legend.position = "none")

dat_dor_dread_distance <- read.xlsx("~cchen2/Documents/neuroscience/Pn\ project/Data_analysis/Dor_dcb/DORcre_DREAD_perfusion.xlsx", sheet = "Dread_perfusion_OF") %>% 
  as_tibble() %>% 
  filter(Quality =="P") %>% 
  mutate(Solution = factor(Solution, level = c("Saline", "CNO"))) %>% 
  ggplot(., aes(Solution, Total_distance, colour = Solution))+
  geom_boxplot()+
  geom_jitter(color = "gray")+
  scale_color_manual(values = c("#0099B4FF","#AD002AFF"))+
  labs(x="", y="Center time (s)")+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  theme(legend.title = element_blank(), legend.position = "none")






