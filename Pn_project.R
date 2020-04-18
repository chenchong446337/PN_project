## the script for pN project in Stanford

library(openxlsx)
library(ggplot2)
library(cowplot)
library(extrafont)
library("plyr")
library("dplyr")
library(readABF)


## For CNO data------
dat_cno <- read.xlsx("~cchen2/Documents/neuroscience/Pn\ project/Data_analysis/Pn_project.xlsx", sheet = 2,colNames = T)
dat_cno$Group <- factor(dat_cno$Group,c("Saline", "CNO") )
dat_cno_sta <- ddply(dat_cno,.(Time,Group),summarise,n=length(Threshold),mean=mean(Threshold),sd=sd(Threshold),se=sd(Threshold)/sqrt(length(Threshold)))



p_cno <- ggplot(dat_cno_sta, aes(Time, mean, Group=Group, colour=Group))+
   geom_line()+
  geom_point(size=2)+
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=0.4) +
  geom_line()+
  labs(x="Time (min)", y="Withdrawal threshold (g)")+
  theme(axis.line = element_line(colour = "black"),
        axis.text.x = element_text(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(size = 12, face ="plain"))+
  theme(legend.position=c(0.05, 0.9), legend.title = element_blank())

dat_cno_treat <- subset(dat_cno, dat_cno$Time=="60"|dat_cno$Time=="90")

dat_cno_tret_sta <- ddply(dat_cno_treat,.(ID, Group),summarise,n=length(Threshold),mean=mean(Threshold),sd=sd(Threshold),se=sd(Threshold)/sqrt(length(Threshold)))

p_cno_60 <- ggplot(dat_cno_tret_sta, aes(Group, mean, colour=Group))+
  geom_boxplot(outlier.shape = "")+
  geom_jitter(width = 0.25, shape=1, size=2)+
  labs(x="", y="Withdrawal Threshold (g)")+
  #scale_color_manual(values = c("black", "red"))+
  theme(axis.line = element_line(colour = "black"),
        axis.text.x = element_text(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(size = 12, face ="plain"))+
  # scale_y_continuous(limits = c(0, 1))+
  theme( legend.position="none",legend.title = element_blank())

## compare the baseline

dat_cno_base <- subset(dat_cno, dat_cno$Time=="-30"|dat_cno$Time=="0")
dat_cno_base_sta <- ddply(dat_cno_base,.(ID, Group),summarise,n=length(Threshold),mean=mean(Threshold),sd=sd(Threshold),se=sd(Threshold)/sqrt(length(Threshold)))

p_cno_base <- ggplot(dat_cno_base_sta, aes(Group, mean, colour=Group))+
  geom_boxplot(outlier.shape = "")+
  geom_jitter(width = 0.25, shape=1, size=2)+
  labs(x="", y="Withdrawal Threshold (g)")+
  # scale_color_manual(values = c("black", "red"))+
  theme(axis.line = element_line(colour = "black"),
        axis.text.x = element_text(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(size = 12, face ="plain"))+
  # scale_y_continuous(limits = c(0, 1))+
  theme( legend.position="none",legend.title = element_blank())

p_cno_box <- plot_grid(p_cno_base, p_cno_60, ncol = 2)
p_cno_com <- plot_grid(p_cno, p_cno_box, ncol = 1)
t.test(mean~Group, data = dat_cno_tret_sta)
## plot for each mice

p_cno2 <- ggplot(dat_cno, aes(Time, Threshold, Group=ID, colour=Group))+
  geom_line()+
  geom_point(size=2)+
  geom_line()+
  labs(x="Time (min)", y="Withdrawal threshold (g)")+
  theme(axis.line = element_line(colour = "black"),
        axis.text.x = element_text(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(size = 12, face ="plain"))+
  theme(legend.position=c(0.05, 0.9), legend.title = element_blank())



cairo_pdf("p_cno.pdf", width = 160/25.6, height = 180/25.6, family = "Arial")
p_cno_com
dev.off()

cairo_pdf("p_cno1.pdf", width = 140/25.6, height = 90/25.6, family = "Arial")
p_cno2
dev.off()


## Overlay between DOr and traped cells in Pontine-----
dat_ish <- read.xlsx("C:/Users/cchen20/Box/Chong Chen's Files/BLA project/data/Pn proejct/Pn_project.xlsx", sheet = 3,colNames = T)

dat_ish$overlap <- dat_ish$Num_over_dor_trap/dat_ish$Num_trap*100
dat_ish$overlap1 <- dat_ish$Num_over_trap_vgat/dat_ish$Num_trap *100

p_over_dor<- ggplot(dat_ish, aes(Area, overlap, colour=Area)) +
  geom_boxplot()+
  geom_jitter(width = 0.25, shape=1, size=2)+
  labs(x="", y="Num. of Dor+ neuron/Num. of TRAPed neuron")+
  theme(axis.line = element_line(colour = "black"),
        axis.text.x = element_text(angle=45, vjust = 0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  theme(legend.position="none")

p_over_vgat <- ggplot(dat_ish, aes(Area, overlap1, colour=Area)) +
  geom_boxplot(outlier.shape = "")+
  geom_jitter(width = 0.25, shape=1, size=2)+
  labs(x="", y="Num. of Dor+ neuron/Num. of TRAPed neuron")+
  theme(axis.line = element_line(colour = "black"),
        axis.text.x = element_text(angle=45, vjust = 0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  theme(legend.position="none")

cairo_pdf("p_over_dor.pdf", width = 60/25.6, height = 70/25.6, family = "Arial")
p_over_dor
dev.off()

cairo_pdf("p_over_vgat.pdf", width = 60/25.6, height = 70/25.6, family = "Arial")
p_over_vgat
dev.off()

### analysis for offset anagesia------
dat_offset <- read.xlsx("~cchen2/Documents/neuroscience/Pn\ project/Data_analysis/Pn_project.xlsx", sheet = "offset analgesia",colNames = T)
dat_offset <- subset(dat_offset, dat_offset$Temperature!="52")
dat_offset$Temperature <- as.factor(dat_offset$Temperature)
dat_offset_sta <- ddply(dat_offset,.(ID, Temperature),summarise,n=length(Latency),mean=mean(Latency),sd=sd(Latency),se=sd(Latency)/sqrt(length(Latency)))

p_offset <- ggplot(dat_offset)+
  geom_boxplot(aes(x = Temperature, y = Latency, group = Temperature, colour=Temperature))+
  geom_point(aes(x = Temperature, y = Latency, colour=Temperature)) +
  geom_line(aes(x = Temperature, y = Latency, group = ID), colour="gray") +
  labs(x="", y="Latency (s)")+
  theme(axis.line = element_line(colour = "black"),
        axis.text.x = element_text(angle=45, vjust = 0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  scale_x_discrete(labels=c("48_before", "48_after"))+
  theme(legend.position="none", legend.title = element_blank())

cairo_pdf("p_offset.pdf", width = 75/25.6, height = 90/25.6, family = "Arial")
p_offset
dev.off()

### Plot for tail flick data------------
dat_flick <- read.xlsx("~cchen2/Documents/neuroscience/Pn\ project/Data_analysis/Pn_project.xlsx", sheet = "Tail_flick",colNames = T)
dat_flick$Group <-factor(dat_flick$Group, c("Saline", "CNO") )
dat_flick_sta <- ddply(dat_flick,.(Group, Temperature),summarise,n=length(Latency),mean=mean(Latency),sd=sd(Latency),se=sd(Latency)/sqrt(length(Latency)))
dat_flick_48 <- subset(dat_flick, dat_flick$Temperature=="48")
dat_flick_52 <- subset(dat_flick, dat_flick$Temperature=="52")

p_flick_48 <- ggplot(dat_flick_48)+
  geom_boxplot(aes(x = Group, y = Latency, colour=Group))+
  geom_point(aes(x = Group, y = Latency, colour=Group)) +
  geom_line(aes(x = Group, y = Latency, group = ID), colour="gray")+
  labs(x="", y="Latency (s)")+
  theme(axis.line = element_line(colour = "black"),
        axis.text.x = element_text(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  theme(legend.position="none", legend.title = element_blank())


p_flick_52 <- ggplot(dat_flick_52)+
  geom_boxplot(aes(x = Group, y = Latency, colour=Group))+
  geom_point(aes(x = Group, y = Latency, colour=Group)) +
  geom_line(aes(x = Group, y = Latency, group = ID), colour="gray")+
  labs(x="", y="Latency (s)")+
  theme(axis.line = element_line(colour = "black"),
        axis.text.x = element_text(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  theme(legend.position="none", legend.title = element_blank())

p_temp <- plot_grid(p_flick_48, p_flick_52)

cairo_pdf("p_temp.pdf", width = 110/25.6, height = 90/25.6, family = "Arial")
p_temp
dev.off()

## ephys for pn cells-----
dat_firing <- read.xlsx("~cchen2/Documents/neuroscience/Pn\ project/Data_analysis/Pn_project.xlsx", sheet = "Ephys_firing_pn",colNames = T)
dat_firing <- subset(dat_firing, dat_firing$Group!="Washout")
dat_firing$Group <- factor(dat_firing$Group, c("Ctrl", "Detra", "DAMGO"))

p_firing<- ggplot(data = dat_firing, aes(Current, Firing, colour=Group)) +
  geom_point()+
  geom_line()+
  labs(x="Current (pA)", y="Firing freq. (Hz)")+
  theme(axis.line = element_line(colour = "black"),
        axis.text.x = element_text(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  scale_y_continuous(limits = c(0, 160))+
  theme(legend.position=c(0.1, 0.8), legend.title = element_blank())


dat_rmp <-  read.xlsx("C:/Users/cchen20/Box/Chong Chen's Files/BLA project/data/Pn proejct/Pn_project.xlsx", sheet = "Ephys_RMP",colNames = T)
dat_rmp$Group <- factor(dat_rmp$Group, c("Ctrl", "Detra", "DAMGO"))

p_rmp <- ggplot(dat_rmp)+
  geom_boxplot(aes(x = Group, y = Current, colour=Group))+
  geom_point(aes(x = Group, y = Current, colour=Group)) +
  geom_line(aes(x = Group, y = Current, group = ID), colour="gray")+
  labs(x="", y="Current (pA)")+
  theme(axis.line = element_line(colour = "black"),
        axis.text.x = element_text(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  theme(legend.position="none", legend.title = element_blank())

p_ephys <- plot_grid(p_rmp, p_firing, ncol = 1)

cairo_pdf("p_ephys.pdf", width = 90/25.6, height = 160/25.6, family = "Arial")
p_ephys
dev.off()
## plot the ephys trace of current injection-----
work_dic <- "C:/Users/cchen20/Box/Chong Chen's Files/shibin/data analysis/"

trace_ctrl <- readABF("C:/Users/cchen20/Box/Chong Chen's Files/BLA project/data/Pn proejct/032019/19320009.abf")
trace_ctrl_n50 <- trace_ctrl[[4]][[1]][,1]
trace_ctrl_100 <- trace_ctrl[[4]][[4]][,1]

plot(trace_ctrl_n50, type = 'l', xlab = "", ylab = "", ylim = c(-90,10 ), axes = F, col="#F8766D")
lines(trace_ctrl_100, col="#F8766D")
segments(13000, -50, 13000, -30)

## for trace with deltra II
trace_deltra <- readABF("C:/Users/cchen20/Box/Chong Chen's Files/BLA project/data/Pn proejct/032019/19320012.abf")
trace_deltra_n50 <- trace_deltra[[4]][[1]][,1]
trace_deltra_100 <- trace_deltra[[4]][[4]][,1]

plot(trace_deltra_n50, type = 'l', xlab = "", ylab = "", ylim = c(-90,10 ), axes = F, col="#00BA38")
lines(trace_deltra_100, col="#00BA38")
segments(13000, -50, 13000, -30)

## for trace with DAMGO
trace_damgo <- readABF("C:/Users/cchen20/Box/Chong Chen's Files/BLA project/data/Pn proejct/032019/19320015.abf")
trace_damgo_n50 <- trace_damgo[[4]][[1]][,1]
trace_damgo_100 <- trace_damgo[[4]][[4]][,1]

plot(trace_damgo_n50, type = 'l', xlab = "", ylab = "", ylim = c(-90,10 ), axes = F, col="#619CFF")
lines(trace_damgo_100, col="#619CFF")
segments(13000, -50, 13000, -30)

## for the two hot plates data-----
dat_hot <- read.xlsx("~cchen2/Documents/neuroscience/Pn\ project/Data_analysis/Pn_project.xlsx", sheet = "two_plate",colNames = T)
dat_hot$Group <- factor(dat_hot$Group,c("Saline", "CNO") )
dat_hot$Ratio <- dat_hot$Ratio*100
dat_hot <- subset(dat_hot, dat_hot$T2!=30)
dat_hot$T2 <- as.factor(dat_hot$T2)
dat_hot_sta <- ddply(dat_hot,.(Group,T2),summarise,n=length(Ratio),mean=mean(Ratio),sd=sd(Ratio),se=sd(Ratio)/sqrt(length(Ratio)))
dat_hot_15 <- subset(dat_hot,dat_hot$T2=="15")




p_hot <- ggplot(dat_hot,aes(y = Ratio, x = T2, colour = Group)) + 
  geom_boxplot()+
  geom_point(position=position_jitterdodge(jitter.width=0.1), aes(group=Group), shape=1)+
  labs(x="Test temperature", y="Time at test temperature (%)")+
  theme(axis.line = element_line(colour = "black"),
        axis.text.x = element_text(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  scale_x_discrete(labels=c("15 °C", "40 °C"))+
  theme(legend.title = element_blank())

## plot for opto_CFA experiment-----
dat_opt_von <- read.xlsx("~cchen2/Documents/neuroscience/Pn\ project/Data_analysis/Pn_project.xlsx", sheet = "Von_CFA",colNames = T)
dat_opt_von$Time<-factor(dat_opt_von$Time,c("D0", "D1", "D3", "D5", "D7", "D9"))
dat_opt_von$Group<- factor(dat_opt_von$Group, c("Off","On"))


dat_opt_von_sta <- ddply(dat_opt_von, .(Time, Group), summarise,n=length(Threshold),mean=mean(Threshold),sd=sd(Threshold),se=sd(Threshold)/sqrt(length(Threshold)))


p_opt_von<- ggplot(dat_opt_von_sta, aes(x=Time, y=mean, group=Group)) + 
  geom_point(aes(colour=Group), size=2)+
  geom_line(aes(colour=Group)) +
  geom_errorbar(aes(ymin=mean-se,ymax=mean+se, colour=Group), width=.2) +
  xlab("Days after CFA injection")+ylab("Mechanical nociceptive threshold (g)") +
  theme_bw() +
  theme(axis.line = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title.y=element_text(family = "Arial",size = 12, face ="plain"))+
  theme(legend.position=c(0.8,0.85), legend.title = element_blank())
  
## for the hargvest test (thermol pain)
dat_opt_har <- read.xlsx("~cchen2/Documents/neuroscience/Pn\ project/Data_analysis/Pn_project.xlsx", sheet = "Har_CFA",colNames = T)
dat_opt_har$Time<-factor(dat_opt_har$Time,c("D0", "D1", "D3", "D5","D7", "D9"))
dat_opt_har$Group<- factor(dat_opt_har$Group, c("Off","On"))


dat_opt_har_sta <- ddply(dat_opt_har, .(Time, Group), summarise,n=length(Latency),mean=mean(Latency),sd=sd(Latency),se=sd(Latency)/sqrt(length(Latency)))


p_opt_har<- ggplot(dat_opt_har_sta, aes(x=Time, y=mean, group=Group)) + 
  geom_point(aes(colour=Group), size=2)+
  geom_line(aes(colour=Group)) +
  geom_errorbar(aes(ymin=mean-se,ymax=mean+se, colour=Group), width=.2) +
  xlab("Days after CFA injection")+ylab("Paw withdraw latency (s)") +
  theme_bw() +
  theme(axis.line = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title.y=element_text(family = "Arial",size = 12, face ="plain"))+
  theme(legend.position=c(0.8,0.85), legend.title = element_blank())

## combine the plot together
p_opt <- plot_grid(p_opt_von, p_opt_har, ncol = 1)

setwd("~cchen2/Documents/neuroscience/Pn\ project/Figure/PDF/")
cairo_pdf("pn_opto_behavior.pdf", width = 120/25.6, height = 160/25.6, family = "Arial")
p_opt
dev.off()



## plot for the firing with Delta II of ACC neurons-----

dat_firing <- read.xlsx("~cchen2/Documents/neuroscience/Pn\ project/Data_analysis/Pn_project.xlsx", sheet = "Ephys_firing_ACC",colNames = T)
dat_firing_sta <- ddply(dat_firing, .( Group, Current), summarise,n=length(Firing),mean=mean(Firing),sd=sd(Firing),se=sd(Firing)/sqrt(length(Firing)))
p_firing<- ggplot(data = dat_firing_sta, aes(Current, mean, colour=Group)) +
  geom_point(aes(colour=Group), size=2)+
  geom_line(aes(colour=Group)) +
  geom_errorbar(aes(ymin=mean-se,ymax=mean+se, colour=Group), width=.2) +
  labs(x="Current (pA)", y="Firing freq. (Hz)")+
  theme(axis.line = element_line(colour = "black"),
        axis.text.x = element_text(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  scale_y_continuous(limits = c(0, 40))+
  theme(legend.position=c(0.1, 0.8), legend.title = element_blank())

## for pain anticipation experiment(WT mice, 07162019)-----
dat_anti_con <- read.xlsx("~cchen2/Documents/neuroscience/Pn\ project/Data_analysis/Pn_project.xlsx", sheet = "Pain_anti",colNames = T)
dat_anti_con$ratio <- dat_anti_con$t_50/dat_anti_con$T_30 *100
dat_anti_con_tmove <- ddply(dat_anti_con, .(Time), summarise, n=length(t_move1),mean=mean(t_move1),sd=sd(t_move1),se=sd(t_move1)/sqrt(length(t_move1)))
dat_anti_con_ratio_sta <- ddply(dat_anti_con, .(Time), summarise, n=length(ratio),mean=mean(ratio),sd=sd(ratio),se=sd(ratio)/sqrt(length(ratio)))
  
p_tmove <- ggplot(data = dat_anti_con_tmove, aes(Time, mean))+
  geom_point(size=2)+
  geom_line() +
  geom_errorbar(aes(ymin=mean-se,ymax=mean+se), width=.2) +
  labs(x="Tranining days", y="Latency (s)")+
  theme(axis.line = element_line(colour = "black"),
        axis.text.x = element_text(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  scale_y_continuous(limits = c(0, 40))+
  theme(legend.position=c(0.1, 0.8), legend.title = element_blank())


p_ratio <- ggplot(data = dat_anti_con_ratio_sta, aes(Time, mean))+
  geom_point(size=2)+
  geom_line() +
  geom_errorbar(aes(ymin=mean-se,ymax=mean+se), width=.2) +
  labs(x="Tranining days", y="Proportion of time in 50 (%) ")+
  theme(axis.line = element_line(colour = "black"),
        axis.text.x = element_text(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  scale_y_continuous(limits = c(0, 20))
  #theme(legend.position=c(0.1, 0.8), legend.title = element_blank())

## analyze and plot data on test day

## analyze the anti experiment with HM4di injection----
dat_anti_hm4di <- read.xlsx("~cchen2/Documents/neuroscience/Pn\ project/Data_analysis/Pn_project.xlsx", sheet = "p_anti_cno",colNames = T)
dat_anti_hm4di$Group <- factor(dat_anti_hm4di$Group, c("Sham", "hm4DI") )

p_lick <- ggplot(dat_anti_hm4di, aes(x= Group, y = T_lick, colour=Group))+
  geom_boxplot()+
  geom_point()+
  labs(x="", y="Latency of first licking (s) ")+
  theme(axis.line = element_line(colour = "black"),
        axis.text.x = element_text(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  theme(legend.title = element_blank())
  
p_corss1 <- ggplot(dat_anti_hm4di, aes(x= Group, y = T_1_2, colour=Group))+
  geom_boxplot()+
  geom_point()+
  labs(x="", y="Latency of crossing (s)  ")+
  theme(axis.line = element_line(colour = "black"),
        axis.text.x = element_text(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  theme(legend.title = element_blank())

p_corss2 <- ggplot(dat_anti_hm4di, aes(x= Group, y = T_2_1, colour=Group))+
  geom_boxplot()+
  geom_point()+
  labs(x="", y="Latency of crossing (s)")+
  theme(axis.line = element_line(colour = "black"),
        axis.text.x = element_text(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  theme(legend.title = element_blank())

p_ratio2 <- ggplot(dat_anti_hm4di, aes(x= Group, y = R_z2, colour=Group))+
  geom_boxplot()+
  geom_point()+
  labs(x="", y="Proportion of staying in zone 2 (%)  ")+
  theme(axis.line = element_line(colour = "black"),
        axis.text.x = element_text(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  theme(legend.title = element_blank())

p_anti <- plot_grid(p_corss1, p_lick, p_corss2, p_ratio2, ncol = 2)

setwd("~cchen2/Documents/neuroscience/Pn\ project/Figure/PDF/")
cairo_pdf("p_anti.pdf", width = 160/25.6, height = 170/25.6, family = "Arial")
p_anti
dev.off()

## The work directory for saving figures----
setwd("~cchen2/Documents/neuroscience/Pn\ project/Figure/PDF/")
