## 1. only compare the first time when they cross the border during D3-D7

## 2. compare the cell activation change for different crossing during one day



## function to extract trace crossing the border
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
path_trace_m3 <- as.list(list.files(path ="~cchen/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/matlab_analysis/m3/trace_days/", pattern = "*.xlsx", full.names = T ))[c(3:7)]

t_stim_m3_d3 <- c(1885) 
t_stim_m3_d4 <- c(3587)
t_stim_m3_d5 <- c(758) 
t_stim_m3_d6 <- c(507) 
t_stim_m3_d7 <- c(346) ## use first two 346, 1206, 1756, 2288, 2970, 3262, 3942

t_stim_m3 <- list( t_stim_m3_d3, t_stim_m3_d4, t_stim_m3_d5, t_stim_m3_d6, t_stim_m3_d7)

dat_trace_m3 <- mapply(c_miniscope_matlab, path_trace_m3, t_stim_m3, SIMPLIFY = F)

# for m7
path_trace_m7 <- as.list(list.files(path ="~cchen/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/matlab_analysis/m7/trace_days/", pattern = "*.xlsx", full.names = T ))[c(3:7)]

t_stim_m7_d3 <- c(360) 
t_stim_m7_d4 <- c(3236)
t_stim_m7_d5 <- c(317) 
t_stim_m7_d6 <- c(205)
t_stim_m7_d7 <- c(236) ## use first two 236, 1914,2777,3334,4091

t_stim_m7 <- list( t_stim_m7_d3, t_stim_m7_d4, t_stim_m7_d5, t_stim_m7_d6, t_stim_m7_d7)

dat_trace_m7 <- mapply(c_miniscope_matlab, path_trace_m7, t_stim_m7,SIMPLIFY = F)

# for m17
path_trace_m17 <- as.list(list.files(path ="~cchen/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/matlab_analysis/m17/trace_days/", pattern = "*.xlsx", full.names = T ))[c(3:7)]

t_stim_m17_d3 <- c(437) 
t_stim_m17_d4 <- c(1745)
t_stim_m17_d5 <- c(1114) 
t_stim_m17_d6 <- c(617) 
t_stim_m17_d7 <- c(157) # 157, 2536, 2961

t_stim_m17 <- list(t_stim_m17_d3, t_stim_m17_d4, t_stim_m17_d5, t_stim_m17_d6, t_stim_m17_d7)

dat_trace_m17 <- mapply(c_miniscope_matlab, path_trace_m17, t_stim_m17, SIMPLIFY = F)

# for m18
path_trace_m18 <- as.list(list.files(path ="~cchen/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/matlab_analysis/m18/trace_days/", pattern = "*.xlsx", full.names = T ))[c(3:7)]

t_stim_m18_d3 <- c(784) 
t_stim_m18_d4 <- c(1580)
t_stim_m18_d5 <- c(195) 
t_stim_m18_d6 <- c(467) 
t_stim_m18_d7 <- c(124) ## use the first two of 124,2238,3000

t_stim_m18 <- list( t_stim_m18_d3, t_stim_m18_d4, t_stim_m18_d5, t_stim_m18_d6, t_stim_m18_d7)

dat_trace_m18 <- mapply(c_miniscope_matlab, path_trace_m18, t_stim_m18, SIMPLIFY = F)

path_trace_m855 <- as.list(list.files(path ="~cchen/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/matlab_analysis/m855/trace_days", pattern = "*.xlsx", full.names = T ))[c(3:7)]
t_stim_m855_d3 <- c(392*2) 
t_stim_m855_d4 <- c(550*2)
t_stim_m855_d5 <- c(67*2) 
t_stim_m855_d6 <- c(64*2) 
t_stim_m855_d7 <- c(41*2) ## use the first two of 124,2238,3000

t_stim_m855 <- list( t_stim_m855_d3, t_stim_m855_d4, t_stim_m855_d5, t_stim_m855_d6, t_stim_m855_d7)
dat_trace_m855 <- mapply(c_miniscope_matlab, path_trace_m855, t_stim_m855, SIMPLIFY = F)


# for m857
path_trace_m857 <- as.list(list.files(path ="~cchen/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/matlab_analysis/m857/trace_days", pattern = "*.xlsx", full.names = T ))[c(1:5)]
t_stim_m857_d3 <- c(847*2) 
t_stim_m857_d4 <- c(252*2)
t_stim_m857_d5 <- c(1900*2) 
t_stim_m857_d6 <- c(723*2) 
t_stim_m857_d7 <- c(167*2) ## use the first two of 124,2238,3000

t_stim_m857 <- list( t_stim_m857_d3, t_stim_m857_d4, t_stim_m857_d5, t_stim_m857_d6, t_stim_m857_d7)
dat_trace_m857 <- mapply(c_miniscope_matlab, path_trace_m857, t_stim_m857, SIMPLIFY = F)


## combine data and do k-means analysis-----
dat_cell_trace <- vector(mode = "list", 5)
for (i in 1:5) {
  dat_cell_trace[[i]]<- cbind(dat_trace_m3[[i]], dat_trace_m7[[i]], dat_trace_m17[[i]], dat_trace_m18[[i]], dat_trace_m855[[i]], dat_trace_m857[[i]] )
  
}

dat_cell_trace_re <- vector(mode = "list", 5)
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
  rep_time <- c(ncol(dat_trace_m3[[i]]), ncol(dat_trace_m7[[i]]), ncol(dat_trace_m17[[i]]), ncol(dat_trace_m18[[i]]), ncol(dat_trace_m855[[i]]), ncol(dat_trace_m857[[i]]))
  dat_cell_trace_d_re$ID <- rep(c("m3", "m7", "m17", "m18", "m855", "m857"), rep_time*length(stim_time))
  
  dat_cell_trace_re[[i]] <- dat_cell_trace_d_re
}

## heat plot of d3, d6 and d7----
## set the range of z score during all days
score_range <- range(mapply(function(x) x$value, dat_cell_trace_re, SIMPLIFY = T))

for (i in c(1:5)) {
  dat_trace <- dat_cell_trace_re[[i]]
  dat_trace_sta <- ddply(dat_trace, .(variable, Group), summarise,mean=mean(value), sum=sum(value))
  dat_trace_sta <- dat_trace_sta[order(dat_trace_sta[,'mean']),]
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

## combine the heat plot
p_heat_com <- plot_grid(p_heat_d1, p_heat_d2,p_heat_d3, p_heat_d4, p_heat_d5, nrow = 1)



## plot trace by group-----
dat_cell_trace_sta <- c()
for (i in c(1:5)){
  dat_trace_sta <- ddply(dat_cell_trace_re[[i]], .(ID,Time, Group), summarise,value1=mean(value)) %>% 
    ddply(., .(Time, Group), summarise, n=length(value1),mean=mean(value1),sd=sd(value1),se=sd(value1)/sqrt(length(value1)))
  
  dat_cell_trace_sta <- rbind(dat_cell_trace_sta, dat_trace_sta)
}


dat_cell_trace_sta$Day <- rep(c( "D3","D4","D5","D6","D7"),each= length(stim_time)*3)
dat_cell_trace_sta$Day <- factor(dat_cell_trace_sta$Day, levels = c( "D3","D4","D5","D6","D7"))
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


## calculate the sum of active and inhibited trace----
rep_tiem <- mapply(nrow, dat_cell_trace_re)
dat_cell_trace_sum <- do.call(rbind, dat_cell_trace_re) %>% 
  as_tibble() %>% 
  mutate(Day = rep(c("D3","D4","D5","D6","D7"), times = rep_tiem)) %>%
  subset(., .$Group!="Neutral") %>%
  ddply(., .(ID,Time, Day), summarise, value=sum(value)) %>%
  ddply(., .(Time, Day), summarise,n=length(value),mean=mean(value),sd=sd(value),se=sd(value)/sqrt(length(value))) %>%
  mutate(Day = factor(Day, levels = c("D3","D4","D5","D6","D7")))

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

dat_cell_area <- c()
dat_prop_group <- c()
cell_area_day <- c("D1", "D2", "D3","D4","D5","D6","D7")

for (i in c(1:7)){
  rep_time <- c(ncol(dat_trace_m3[[i]]), ncol(dat_trace_m7[[i]]), ncol(dat_trace_m17[[i]]), ncol(dat_trace_m18[[i]]))
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

dat_cell_area$Day <- factor(dat_cell_area$Day, levels = c("D1", "D2", "D3","D4","D5","D6","D7"))
dat_cell_area_com <- ddply(dat_cell_area, .( Day), summarise, value=mean(ratio, na.rm = T))
dat_cell_area_sta <- ddply(dat_cell_area_com,.(Day), summarise,n=length(value, na.rm=T),mean=mean(value),sd=sd(value),se=sd(value)/sqrt(length(value)))



## cells correlation analysis-----
cc_cov_fun <- function(cell_trace){
  dat_cell_trace<- subset(cell_trace, cell_trace$Group!="Neutral")
  ID_cell <- unique(dat_cell_trace$ID)
  dat_cell_cov <- vector(mode = "list", length(ID_cell))
  
  for (i in seq_along(ID_cell)){
    dat_cell_trace_cov <- subset(dat_cell_trace, dat_cell_trace$ID==ID_cell[i])
    dat_cell_trace_cov <- dcast(dat_cell_trace_cov[,1:3], Time~variable)
    if (ncol(dat_cell_trace_cov) < 3){
      res_cov <- NA
    } else{
      res_cov <- cov(dat_cell_trace_cov[,-1])
      
    }
    dat_cell_cov[[i]] <- res_cov
  }
  return(dat_cell_cov)
}


dat_cell_cov_list <- mapply(cc_cov_fun, dat_cell_trace_re, SIMPLIFY = F)
cov_range <- c()
for (i in c(1:7)){
  cov_range1 <- range(dat_cell_cov_list[[i]][[1]])
  cov_range <- range(cov_range, cov_range1)
  p_cov<- ggcorrplot(dat_cell_cov_list[[i]][[1]],hc.order = TRUE, outline.col = "white", ggtheme = ggplot2::theme_void,show.legend = F)+
    scale_fill_gradient2(limit = c(cov_range[1], cov_range[2]), low = "navy", high =  "red4", mid = "white")
  assign(paste0("p_cov_d", i), p_cov)
  
}

p_cov_com <- plot_grid(p_cov_d1, p_cov_d2,p_cov_d3, p_cov_d4, p_cov_d5,p_cov_d6, p_cov_d7, nrow = 2)



dat_cell_cov <- c()
dat_cell_cov_sta <- c()
for (i in c(1:7)){
  cov_value <- unlist(dat_cell_cov_list[[i]])
  cov_value <- abs(cov_value[!is.na(cov_value)])
  dat_cov <- data.frame(Day=cell_area_day[i], value=cov_value)
  dat_cell_cov <- rbind(dat_cell_cov, dat_cov)
  
  ## for mean of cov
  cov_value_mean <- mapply(function(x) ifelse(ncol(x)<2, NULL ,mean(colMeans(abs(x), na.rm=T))), dat_cell_cov_list[[i]])
  cov_value_max <- mapply(function(x) ifelse(ncol(x)<2, NA,mean(apply(x, 2, function(x) max(abs(x), na.rm = T)))), dat_cell_cov_list[[i]])
  
  mouse_ID <- c(1: length(dat_cell_cov_list[[i]]))
  dat_cov_mean <- data.frame(Day = cell_area_day[i], ID= mouse_ID, value=cov_value_mean, value_max=cov_value_max)
  dat_cell_cov_sta <- rbind(dat_cell_cov_sta, dat_cov_mean)
}


##2. for ca2+ transit during one day------
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
path_trace_m3 <- "~cchen/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/matlab_analysis/m3/trace_days/m3_trace_d2.xlsx"
t_stim_m3_d2 <- c(274, 831,1364)

dat_stim_m3_d2 <- c_miniscope_matlab_d7(path_trace_m3, t_stim_m3_d2)


## for m7
path_trace_m7 <- "~cchen/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/matlab_analysis/m7/trace_days/m7_trace_d2.xlsx"
t_stim_m7_d2 <- c(798, 177,3300)

dat_stim_m7_d2 <- c_miniscope_matlab_d7(path_trace_m7, t_stim_m7_d2)

## for m17
path_trace_m17 <- "~cchen/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/matlab_analysis/m17/trace_days/m17_trace_d2.xlsx"
t_stim_m17_d2 <- c(1103,1678,2291)

dat_stim_m17_d2 <- c_miniscope_matlab_d7(path_trace_m17, t_stim_m17_d2)


## for m18
path_trace_m18 <- "~cchen/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/matlab_analysis/m18/trace_days/m18_trace_d2.xlsx"
t_stim_m18_d2 <- c(236, 898,1528)

dat_stim_m18_d2 <- c_miniscope_matlab_d7(path_trace_m18, t_stim_m18_d2)

dat_cell_trace <- vector(mode = "list", 3)
for (i in 1:3) {
  dat_cell_trace[[i]]<- cbind(dat_stim_m3_d2[[i]], dat_stim_m7_d2[[i]], dat_stim_m17_d2[[i]], dat_stim_m18_d2[[i]] )
  
}

dat_cell_trace_re <- vector(mode = "list", 3)
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
  rep_time <- c(ncol(dat_stim_m3_d2[[i]]), ncol(dat_stim_m7_d2[[i]]), ncol(dat_stim_m17_d2[[i]]), ncol(dat_stim_m18_d2[[i]]))
  dat_cell_trace_d_re$ID <- rep(c("m3", "m7", "m17", "m18"), rep_time*length(stim_time))
  
  dat_cell_trace_re[[i]] <- dat_cell_trace_d_re
}

dat_cell_trace_re[[2]]$Group<-dat_cell_trace_re[[1]]$Group
dat_cell_trace_re[[3]]$Group<-dat_cell_trace_re[[1]]$Group

score_range <- range(mapply(function(x) x$value, dat_cell_trace_re, SIMPLIFY = T))

for (i in c(1:3)) {
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

p_heat_com <- plot_grid(p_heat_d1, p_heat_d2,p_heat_d3, nrow = 1)


dat_cell_trace_stim1_sta <- ddply(dat_cell_trace_re[[1]], .(Time, Group), summarise,n=length(value),mean=mean(value),sd=sd(value),se=sd(value)/sqrt(length(value)))
dat_cell_trace_stim2_sta <- ddply(dat_cell_trace_re[[2]], .(Time, Group), summarise,n=length(value),mean=mean(value),sd=sd(value),se=sd(value)/sqrt(length(value)))
dat_cell_trace_stim3_sta <- ddply(dat_cell_trace_re[[3]], .(Time, Group), summarise,n=length(value),mean=mean(value),sd=sd(value),se=sd(value)/sqrt(length(value)))

dat_cell_trace_sta <- rbind(dat_cell_trace_stim1_sta,  dat_cell_trace_stim2_sta, dat_cell_trace_stim3_sta)

dat_cell_trace_sta$Day <- rep(c("stim1",  "stim2", "stim3"),each= length(stim_time)*3)
dat_cell_trace_sta$Day <- factor(dat_cell_trace_sta$Day, levels = c("stim1",  "stim2", "stim3"))
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
cell_area_day <- c("stim1",  "stim2", "stim3")
for (i in c(1:3)){
  rep_time <- c(ncol(dat_stim_m3_d2[[i]]), ncol(dat_stim_m7_d2[[i]]), ncol(dat_stim_m17_d2[[i]]), ncol(dat_stim_m18_d2[[i]]))
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

dat_cell_area$Day <- factor(dat_cell_area$Day, levels = c("stim1","stim2","stim3"))
dat_cell_area_com <- ddply(dat_cell_area, .(ID, Day), summarise, value=mean(ratio, na.rm = T))
dat_cell_area_sta <- ddply(dat_cell_area_com,.(Day), summarise,n=length(value),mean=mean(value),sd=sd(value),se=sd(value)/sqrt(length(value)))


p_EI_ratio <- ggplot(dat_cell_area_sta, aes(Day, mean, colour=Day, fill=Day))+
  geom_bar(stat="identity", color="black", position=position_dodge())+
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.1, position=position_dodge(.9), colour="black") +
  geom_jitter(data = dat_cell_area_com, aes(Day, value),width = 0.1, colour="black", shape=1 )+
  #scale_fill_manual(values=c( "indianred4", "tomato3"))+
  labs(x="", y="E/I ratio")+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  #scale_y_continuous(limits = c(0, 2), expand = c(0,0))+
  theme(legend.title = element_blank())+
  theme(legend.position = "none")

## compare crossing back during conditioning------
path_trace_m3 <- "~cchen/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/matlab_analysis/m3/trace_days/m3_trace_d6.xlsx"
t_stim_m3_d6 <- c(132, 372)

dat_stim_m3_d2 <- c_miniscope_matlab_d7(path_trace_m3, t_stim_m3_d6)


## for m18
path_trace_m18 <- "~cchen/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/matlab_analysis/m18/trace_days/m18_trace_d5.xlsx"
t_stim_m18_d5 <- c(195, 3545)

dat_stim_m18_d2 <- c_miniscope_matlab_d7(path_trace_m18, t_stim_m18_d5)


dat_cell_trace <- vector(mode = "list", 2)
for (i in 1:2) {
  dat_cell_trace[[i]]<- cbind(dat_stim_m3_d2[[i]],  dat_stim_m18_d2[[i]] )
  
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
  rep_time <- c(ncol(dat_stim_m3_d2[[i]]),  ncol(dat_stim_m18_d2[[i]]))
  dat_cell_trace_d_re$ID <- rep(c("m3",  "m18"), rep_time*length(stim_time))
  
  dat_cell_trace_re[[i]] <- dat_cell_trace_d_re
}

dat_cell_trace_re[[2]]$Group<-dat_cell_trace_re[[1]]$Group

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
  assign(paste0("p_heat_d", i), p_heat)
}

p_heat_com <- plot_grid(p_heat_d1, p_heat_d2,nrow = 1)

setwd("~cchen/Documents/neuroscience/Pn\ project/Figure/PDF/")
cairo_pdf("p_anti_crossingback.pdf", width = 150/25.6, height = 65/25.6, family = "Arial")
p_heat_com
dev.off()

dat_cell_trace_stim1_sta <- ddply(dat_cell_trace_re[[1]], .(Time, Group), summarise,n=length(value),mean=mean(value),sd=sd(value),se=sd(value)/sqrt(length(value)))
dat_cell_trace_stim2_sta <- ddply(dat_cell_trace_re[[2]], .(Time, Group), summarise,n=length(value),mean=mean(value),sd=sd(value),se=sd(value)/sqrt(length(value)))

dat_cell_trace_sta <- rbind(dat_cell_trace_stim1_sta,  dat_cell_trace_stim2_sta)

dat_cell_trace_sta$Day <- rep(c("stim1","stim2"),each= length(stim_time)*3)
dat_cell_trace_sta$Day <- factor(dat_cell_trace_sta$Day, levels = c("stim1",  "stim2"))
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
cell_area_day <- c("stim1",  "stim2")
for (i in c(1:2)){
  rep_time <- c(ncol(dat_stim_m3_d2[[i]]),  ncol(dat_stim_m18_d2[[i]]))
  dat_trace <- dat_cell_trace_re[[i]]
  dat_trace$ID <- rep(c("m3",  "m18"), rep_time*length(stim_time))
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

dat_cell_area$Day <- factor(dat_cell_area$Day, levels = c("stim1","stim2"))
dat_cell_area_com <- ddply(dat_cell_area, .(ID, Day), summarise, value=mean(ratio, na.rm = T))
dat_cell_area_sta <- ddply(dat_cell_area_com,.(Day), summarise,n=length(value),mean=mean(value),sd=sd(value),se=sd(value)/sqrt(length(value)))


p_EI_ratio <- ggplot(dat_cell_area_sta, aes(Day, mean, colour=Day, fill=Day))+
  geom_bar(stat="identity", color="black", position=position_dodge())+
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.1, position=position_dodge(.9), colour="black") +
  geom_jitter(data = dat_cell_area_com, aes(Day, value),width = 0.1, colour="black", shape=1 )+
  #scale_fill_manual(values=c( "indianred4", "tomato3"))+
  labs(x="", y="E/I ratio")+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  #scale_y_continuous(limits = c(0, 2), expand = c(0,0))+
  theme(legend.title = element_blank())+
  theme(legend.position = "none")


## compare crossing back during pre-conditioning------
path_trace_m3 <- "~cchen/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/matlab_analysis/m3/trace_days/m3_trace_d3.xlsx"
t_stim_m3_d3 <- c(1885, 2289)

dat_stim_m3_d2 <- c_miniscope_matlab_d7(path_trace_m3, t_stim_m3_d3)
## for m7
path_trace_m7 <- "~cchen/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/matlab_analysis/m7/trace_days/m7_trace_d3.xlsx"
t_stim_m7_d3 <- c(360, 499)

dat_stim_m7_d2 <- c_miniscope_matlab_d7(path_trace_m7, t_stim_m7_d3)

## for m17
path_trace_m17 <- "~cchen/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/matlab_analysis/m17/trace_days/m17_trace_d3.xlsx"
t_stim_m17_d3 <- c(437, 608)

dat_stim_m17_d2 <- c_miniscope_matlab_d7(path_trace_m17, t_stim_m17_d3)


## for m18
path_trace_m18 <- "~cchen/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/matlab_analysis/m18/trace_days/m18_trace_d3.xlsx"
t_stim_m18_d3 <- c(784, 342)

dat_stim_m18_d2 <- c_miniscope_matlab_d7(path_trace_m18, t_stim_m18_d3)


dat_cell_trace <- vector(mode = "list", 2)
for (i in 1:2) {
  dat_cell_trace[[i]]<- cbind(dat_stim_m3_d2[[i]],dat_stim_m7_d2[[i]],dat_stim_m17_d2[[i]],  dat_stim_m18_d2[[i]] )
  
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
  rep_time <- c(ncol(dat_stim_m3_d2[[i]]), ncol(dat_stim_m7_d2[[i]]), ncol(dat_stim_m17_d2[[i]]), ncol(dat_stim_m18_d2[[i]]))
  dat_cell_trace_d_re$ID <- rep(c("m3","m7", "m17", "m18"), rep_time*length(stim_time))
  
  dat_cell_trace_re[[i]] <- dat_cell_trace_d_re
}

dat_cell_trace_re[[2]]$Group<-dat_cell_trace_re[[1]]$Group

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
  assign(paste0("p_heat_d", i), p_heat)
}

p_heat_com <- plot_grid(p_heat_d1, p_heat_d2,nrow = 1)

setwd("~cchen/Documents/neuroscience/Pn\ project/Figure/PDF/")
cairo_pdf("p_anti_crossingback.pdf", width = 150/25.6, height = 65/25.6, family = "Arial")
p_heat_com
dev.off()

dat_cell_trace_stim1_sta <- ddply(dat_cell_trace_re[[1]], .(Time, Group), summarise,n=length(value),mean=mean(value),sd=sd(value),se=sd(value)/sqrt(length(value)))
dat_cell_trace_stim2_sta <- ddply(dat_cell_trace_re[[2]], .(Time, Group), summarise,n=length(value),mean=mean(value),sd=sd(value),se=sd(value)/sqrt(length(value)))

dat_cell_trace_sta <- rbind(dat_cell_trace_stim1_sta,  dat_cell_trace_stim2_sta)

dat_cell_trace_sta$Day <- rep(c("stim1","stim2"),each= length(stim_time)*3)
dat_cell_trace_sta$Day <- factor(dat_cell_trace_sta$Day, levels = c("stim1",  "stim2"))
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
cell_area_day <- c("stim1",  "stim2")
for (i in c(1:2)){
  rep_time <- c(ncol(dat_stim_m3_d2[[i]]),ncol(dat_stim_m7_d2[[i]]), ncol(dat_stim_m17_d2[[i]]),  ncol(dat_stim_m18_d2[[i]]))
  dat_trace <- dat_cell_trace_re[[i]]
  dat_trace$ID <- rep(c("m3","m7", "m17", "m18"), rep_time*length(stim_time))
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

dat_cell_area$Day <- factor(dat_cell_area$Day, levels = c("stim1","stim2"))
dat_cell_area_com <- ddply(dat_cell_area, .(ID, Day), summarise, value=mean(ratio, na.rm = T))
dat_cell_area_sta <- ddply(dat_cell_area_com,.(Day), summarise,n=length(value),mean=mean(value, na.rm = T),sd=sd(value),se=sd(value)/sqrt(length(value)))


p_EI_ratio <- ggplot(dat_cell_area_sta, aes(Day, mean, colour=Day, fill=Day))+
  geom_bar(stat="identity", color="black", position=position_dodge())+
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.1, position=position_dodge(.9), colour="black") +
  geom_jitter(data = dat_cell_area_com, aes(Day, value),width = 0.1, colour="black", shape=1 )+
  #scale_fill_manual(values=c( "indianred4", "tomato3"))+
  labs(x="", y="E/I ratio")+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  #scale_y_continuous(limits = c(0, 2), expand = c(0,0))+
  theme(legend.title = element_blank())+
  theme(legend.position = "none")


## compare for d4-----

c_miniscope_matlab_d <- function(ID_trace, t_stim) {
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
path_trace_m3 <- "~cchen/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/matlab_analysis/m3/trace_days/m3_trace_d4.xlsx"
t_stim_m3_d4 <- c(1813, 3314)

dat_stim_m3_d4 <- c_miniscope_matlab_d(path_trace_m3, t_stim_m3_d4)



## for m17
path_trace_m17 <- "~cchen/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/matlab_analysis/m17/trace_days/m17_trace_d4.xlsx"
t_stim_m17_d4 <- c(682, 1325)

dat_stim_m17_d4 <- c_miniscope_matlab_d(path_trace_m17, t_stim_m17_d4)


## for m18
path_trace_m18 <- "~cchen/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/matlab_analysis/m18/trace_days/m18_trace_d4.xlsx"
t_stim_m18_d4 <- c(493,1390)

dat_stim_m18_d4 <- c_miniscope_matlab_d(path_trace_m18, t_stim_m18_d4)



## combine data and do k-means analysis-----
dat_cell_trace <- vector(mode = "list", 2)
for (i in 1:2) {
  dat_cell_trace[[i]]<- cbind(dat_stim_m3_d4[[i]],  dat_stim_m17_d4[[i]], dat_stim_m18_d4[[i]] )
  
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
  rep_time <- c(ncol(dat_stim_m3_d4[[i]]), ncol(dat_stim_m17_d4[[i]]), ncol(dat_stim_m18_d4[[i]]))
  dat_cell_trace_d_re$ID <- rep(c("m3", "m17", "m18"), rep_time*length(stim_time))
  
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

p_heat_com_d4 <- plot_grid(p_heat_d1, p_heat_d2, nrow = 1)
setwd("~cchen/Documents/neuroscience/Pn\ project/Figure/PDF/")
cairo_pdf("p_anti_heat_d4.pdf", width = 100/25.6, height = 60/25.6, family = "Arial")
p_heat_com_d4
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

setwd("~cchen/Documents/neuroscience/Pn\ project/Figure/PDF/")
cairo_pdf("p_trace_d4.pdf", width = 130/25.6, height = 65/25.6, family = "Arial")
p_trace
dev.off()

dat_cell_area <- c()
cell_area_day <- c("Pre",  "Test")
for (i in c(1:2)){
  rep_time <- c(ncol(dat_stim_m3_d4[[i]]), ncol(dat_stim_m17_d4[[i]]), ncol(dat_stim_m18_d4[[i]]))
  dat_trace <- dat_cell_trace_re[[i]]
  dat_trace$ID <- rep(c("m3",  "m17", "m18"), rep_time*length(stim_time))
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

setwd("~cchen/Documents/neuroscience/Pn\ project/Figure/PDF/")
cairo_pdf("p_EI_d4_comp.pdf", width = 42/25.6, height = 65/25.6, family = "Arial")
p_EI_ratio
dev.off()

## compare the Ca2+ dynamic from day3 to Day6-----
dat_cell_trace_re <- vector(mode = "list", 7)
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
  rep_time <- c(ncol(dat_trace_m3[[i]]), ncol(dat_trace_m7[[i]]), ncol(dat_trace_m17[[i]]), ncol(dat_trace_m18[[i]]))
  dat_cell_trace_d_re$ID <- rep(c("m3", "m7", "m17", "m18"), rep_time*length(stim_time))
  
  dat_cell_trace_re[[i]] <- dat_cell_trace_d_re
}

## heat plot of d3, d6 and d7----
## set the range of z score during all days
score_range <- range(mapply(function(x) x$value, dat_cell_trace_re, SIMPLIFY = T))

for (i in c(3:6)) {
  dat_trace <- dat_cell_trace_re[[i]]
  dat_trace_sta <- ddply(dat_trace, .(variable, Group), summarise,mean=mean(value), sum=sum(value))
  dat_trace_sta <- dat_trace_sta[order(dat_trace_sta[,'mean']),]
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

## combine the heat plot
p_heat_com <- plot_grid(p_heat_d3, p_heat_d4, p_heat_d5,p_heat_d6, nrow = 1)

## plot trace by group-----
dat_cell_trace_sta <- c()
for (i in c(3:6)){
  dat_trace_sta <- ddply(dat_cell_trace_re[[i]], .(Time, Group), summarise,n=length(value),mean=mean(value),sd=sd(value),se=sd(value)/sqrt(length(value)))
  
  dat_cell_trace_sta <- rbind(dat_cell_trace_sta, dat_trace_sta)
}

dat_cell_trace_sta$Day <- rep(c("D3","D4","D5","D6"),each= length(stim_time)*3)
dat_cell_trace_sta$Day <- factor(dat_cell_trace_sta$Day, levels = c( "D3","D4","D5","D6"))
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
dat_prop_group <- c()
cell_area_day <- c( "D3","D4","D5","D6")

for (i in c(3:6)){
  rep_time <- c(ncol(dat_trace_m3[[i]]), ncol(dat_trace_m7[[i]]), ncol(dat_trace_m17[[i]]), ncol(dat_trace_m18[[i]]))
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

dat_cell_area$Day <- factor(dat_cell_area$Day, levels = c("D3","D4","D5","D6"))
dat_cell_area_com <- ddply(dat_cell_area, .( Day), summarise, value=mean(ratio, na.rm = T))
dat_cell_area_sta <- ddply(dat_cell_area_com,.(Day), summarise,n=length(value, na.rm=T),mean=mean(value),sd=sd(value),se=sd(value)/sqrt(length(value)))




