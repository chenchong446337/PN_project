## Global ID for all deteceted cells
cc_globalID_fun <- function(path_ID, path_trace, t_stim){
  Global_ID <- read.csv(path_ID, header = F)
  
  ## change 0 to NA and only keep the completed rows
  Global_ID[Global_ID==0] <- NA
  Global_ID <- Global_ID[complete.cases(Global_ID), ]
 
  
  ## read.xlsx with no header and remove the first col
  cc_read.xlsx<- function(x){
    dat <- read.xlsx(x, colNames = F)
    dat<- dat[,-1]
    dat
  }
  dat_trace_raw <- mapply(cc_read.xlsx, path_trace,SIMPLIFY = F)
  
  ## do z score of the whole trace
  dat_trace <- mapply(function(x) apply(x, 2, scale), dat_trace_raw, SIMPLIFY = F)
  stim_time<- seq(-5, 8.5, by=0.5)
  ## number of rows to be binned
  n <- 10 # 0.05*10=0.5
  
  ## create a 7 days list and put the activation of each cell in it
  cc_trace_pick <- function(dat_trace, cell_pick){
    dat_trace_pick <- t(dat_trace[cell_pick,])
    dat_trace_pick
  }
  
  dat_cell_trace_day <- mapply(cc_trace_pick, dat_trace, as.list(Global_ID))
  
  cc_trace_extract <- function(dat_cell_trace, t_stim_day){
    dat_stim <- vector(mode = "list", length = length(t_stim_day))
    for (i in seq_along(t_stim_day)){
      t1_p <- t_stim_day[i]
      dat_stim1 <- dat_cell_trace[(t1_p-100):(t1_p+180-1),]
      dat_stim1 <- aggregate(dat_stim1,list(rep(1:(nrow(dat_stim1)%/%n+1),each=n,len=nrow(dat_stim1))),mean)[-1]
      dat_stim1_base <- dat_stim1[1:5,] ## baseline as -5 to -3
      dat_stim1_base_mean <- colMeans(dat_stim1_base, na.rm = T)
      dat_stim1_nor1 <- sweep(dat_stim1, 2, dat_stim1_base_mean, FUN = "-")
      rownames(dat_stim1_nor1) <- NULL
      colnames(dat_stim1_nor1)<- NULL
      dat_stim[[i]] <- dat_stim1_nor1
    }
    ## average the trace by cell number
    if (length(t_stim_day)>1){
      dat_cell_trace_average <- aaply(laply(dat_stim, as.matrix), c(2, 3), function(x) mean(x, na.rm=T))
      
    } else {
      dat_cell_trace_average <- data.matrix(dat_stim[[1]])
      
    }
    return(dat_cell_trace_average)
  }
  dat_cell_trace_day_extract <- mapply(cc_trace_extract, dat_cell_trace_day, t_stim, SIMPLIFY = F)
  return(dat_cell_trace_day_extract)
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
t_stim_m3_d7 <- c(346) ## use first two 346, 1206, 1756, 2288
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
t_stim_m7_d7 <- c(236) ## use first two 236, 1914,2777,3334,4091
t_stim_m7 <- list(t_stim_m7_d1, t_stim_m7_d2, t_stim_m7_d3, t_stim_m7_d4, t_stim_m7_d5, t_stim_m7_d6, t_stim_m7_d7)


dat_global_ID_m7 <- cc_globalID_fun(Global_ID_m7, path_trace_m7, t_stim_m7)


## for mice m18
Global_ID_m18 <- "~cchen2/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/matlab_analysis/m18/Global_ID/m18_global_ID.csv"
path_trace_m18 <- list.files(path ="~cchen2/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/matlab_analysis/m18/trace_days/", pattern = "*.xlsx" )
setwd("~cchen2/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/matlab_analysis/m18/trace_days")
t_stim_m18_d1 <- c(490, 950,1291,1714,2205,2610,3088,3521) 
t_stim_m18_d2 <- c(236, 898,1528,2419)
t_stim_m18_d3 <- c(784,1271,2012,3088) 
t_stim_m18_d4 <- c(493,1580,3874)
t_stim_m18_d5 <- c(195,3676) 
t_stim_m18_d6 <- c(467) 
t_stim_m18_d7 <- c(124) ## use the first two of 124,2238,3000
t_stim_m18 <- list(t_stim_m18_d1, t_stim_m18_d2, t_stim_m18_d3, t_stim_m18_d4, t_stim_m18_d5, t_stim_m18_d6, t_stim_m18_d7)

dat_global_ID_m18 <- cc_globalID_fun(Global_ID_m18, path_trace_m18, t_stim_m18)


## combine all data for analysis----
for (i in c(1:7)){
  stim_time<- seq(-5, 8.5, by=0.5)
  dat_extract <- as.data.frame(cbind(dat_global_ID_m3[[i]],dat_global_ID_m7[[i]],dat_global_ID_m18[[i]] ))
  cell_id <- sprintf("Cell%s",seq(1:ncol(dat_extract)))
  colnames(dat_extract) <- cell_id
  dat_extract$Time <- stim_time
  assign(paste0("dat_global_day", i), dat_extract)
}


## k means clustering of data for day6, and plot based on the d6 cluster

dat_global_day6_cluster <- unname(kmeans(t(dat_global_day6[,-ncol(dat_global_day6)]), 3)[[1]])

dat_global_day6_re <- melt(dat_global_day6, id.vars ='Time')
dat_global_day6_re$Group <- rep(dat_global_day6_cluster, each=nrow(dat_global_day6))

## sort data by the value in each group
dat_global_day6_sort <- as.numeric(names(sort(tapply(dat_global_day6_re$value, dat_global_day6_re$Group, mean), decreasing = T)))

dat_global_day6_re$Group[dat_global_day6_re$Group == dat_global_day6_sort[1]] ="Excited"
dat_global_day6_re$Group[dat_global_day6_re$Group == dat_global_day6_sort[2]] ="Neutral"
dat_global_day6_re$Group[dat_global_day6_re$Group == dat_global_day6_sort[3]] ="Inhibited"
order_value <- c("Excited", "Neutral", "Inhibited")
dat_global_day6_re <- dat_global_day6_re[order(match(dat_global_day6_re$Group, order_value)),]

## use the same oder for the other days
dat_global_day6_re_sort <- subset(dat_global_day6_re, dat_global_day6_re$Time==0)
order_value_cell <- dat_global_day6_re_sort$variable

## for day1
dat_global_day1_re <- melt(dat_global_day1, id.vars ='Time')
dat_global_day1_re <- dat_global_day1_re[order(match(dat_global_day1_re$variable, order_value_cell)),]

## for day2
dat_global_day2_re <- melt(dat_global_day2, id.vars ='Time')
dat_global_day2_re <- dat_global_day2_re[order(match(dat_global_day2_re$variable, order_value_cell)),]

## for day3
dat_global_day3_re <- melt(dat_global_day3, id.vars ='Time')
dat_global_day3_re <- dat_global_day3_re[order(match(dat_global_day3_re$variable, order_value_cell)),]

## for day4
dat_global_day4_re <- melt(dat_global_day4, id.vars ='Time')
dat_global_day4_re <- dat_global_day4_re[order(match(dat_global_day4_re$variable, order_value_cell)),]

## for day5
dat_global_day5_re <- melt(dat_global_day5, id.vars ='Time')
dat_global_day5_re <- dat_global_day5_re[order(match(dat_global_day5_re$variable, order_value_cell)),]

## for day7
dat_global_day7_re <- melt(dat_global_day7, id.vars ='Time')
dat_global_day7_re <- dat_global_day7_re[order(match(dat_global_day7_re$variable, order_value_cell)),]

score_range <- range(c(dat_global_day1_re$value, dat_global_day2_re$value, dat_global_day3_re$value, dat_global_day4_re$value,
                       dat_global_day5_re$value, dat_global_day6_re$value, dat_global_day7_re$value))


p_heat_d1 <- ggplot(dat_global_day1_re, aes(Time, variable,fill= value))+ 
  geom_tile()+
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

p_heat_d2 <- ggplot(dat_global_day2_re, aes(Time, variable,fill= value))+ 
  geom_tile()+
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

p_heat_d3 <- ggplot(dat_global_day3_re, aes(Time, variable,fill= value))+ 
  geom_tile()+
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

p_heat_d4 <- ggplot(dat_global_day4_re, aes(Time, variable,fill= value))+ 
  geom_tile()+
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

p_heat_d5 <- ggplot(dat_global_day5_re, aes(Time, variable,fill= value))+ 
  geom_tile()+
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

p_heat_d6 <- ggplot(dat_global_day6_re, aes(Time, variable,fill= value))+ 
  geom_tile()+
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

p_heat_d7 <- ggplot(dat_global_day7_re, aes(Time, variable,fill= value))+ 
  geom_tile()+
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

p_heat_global <- plot_grid(p_heat_d1, p_heat_d2, p_heat_d3, p_heat_d4, p_heat_d5, p_heat_d6, p_heat_d7, nrow = 2)


## only plot day3 and day7 (08092020)-----

cc_globalID_fun <- function(path_ID, path_trace, t_stim){
  Global_ID <- read.csv(path_ID, header = F) %>% 
    na_if(., 0) %>% 
    drop_na()

  
  ## read.xlsx with no header and remove the first col
  cc_read.xlsx<- function(x){
    dat <- read.xlsx(x, colNames = F) %>% 
      select(., -1)
    dat
  }
  
  ## do z score of the whole trace
  dat_trace <- mapply(cc_read.xlsx, path_trace,SIMPLIFY = F) %>% 
    mapply(function(x) apply(x, 2, scale), ., SIMPLIFY = F)
    
  stim_time<- seq(-2, 6.5, by=0.5)
  ## number of rows to be binned
  n <- 10 # 0.05*10=0.5
  
  ## create a 2 days list and put the activation of each cell in it
  cc_trace_pick <- function(dat_trace, cell_pick){
    dat_trace_pick <- t(dat_trace[cell_pick,])
    dat_trace_pick
  }
  
  dat_cell_trace_day <- mapply(cc_trace_pick, dat_trace, as.list(Global_ID))
  
  cc_trace_extract <- function(dat_cell_trace, t_stim_day){
    dat_stim <- vector(mode = "list", length = length(t_stim_day))
    for (i in seq_along(t_stim_day)){
      t1_p <- t_stim_day[i]
      dat_stim1 <- dat_cell_trace[(t1_p-40):(t1_p+140-1),]
      dat_stim1 <- aggregate(dat_stim1,list(rep(1:(nrow(dat_stim1)%/%n+1),each=n,len=nrow(dat_stim1))),mean)[-1]
      dat_stim[[i]] <- dat_stim1[1:5,] %>%  ## baseline as -2 to 0
        colMeans(., na.rm = T) %>% 
        sweep(dat_stim1, 2, ., FUN = "-")
    }
    ## average the trace by cell number
    if (length(t_stim_day)>1){
      dat_cell_trace_average <- aaply(laply(dat_stim, as.matrix), c(2, 3), function(x) mean(x, na.rm=T))
      
    } else {
      dat_cell_trace_average <- data.matrix(dat_stim[[1]])
      
    }
    return(dat_cell_trace_average)
  }
  dat_cell_trace_day_extract <- mapply(cc_trace_extract, dat_cell_trace_day, t_stim, SIMPLIFY = F)
  
  return(dat_cell_trace_day_extract)
}

## for m3
Global_ID_m3 <- "~cchen2/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/matlab_analysis/m3/Global_id/m3_anti_pin_global_ID_37.csv"
path_trace_m3 <- list.files(path ="~cchen2/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/matlab_analysis/m3/trace_days/", pattern = "*.xlsx", full.names = T)[c(3,7)]
t_stim_m3 <- list(t_stim_m3_d3 = c(1885), t_stim_m3_d7 = c(346))
dat_global_ID_m3 <- cc_globalID_fun(Global_ID_m3, path_trace_m3, t_stim_m3)


## for m3
Global_ID_m7 <- "~cchen2/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/matlab_analysis/m7/Global_id/m7_anti_pin_global_ID_37.csv"
path_trace_m7 <- list.files(path ="~cchen2/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/matlab_analysis/m7/trace_days/", pattern = "*.xlsx", full.names = T)[c(3,7)]
t_stim_m7 <- list(t_stim_m7_d3 = c(360), t_stim_m7_d7 = c(236))
dat_global_ID_m7 <- cc_globalID_fun(Global_ID_m7, path_trace_m7, t_stim_m7)

## for m17
Global_ID_m17 <- "~cchen2/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/matlab_analysis/m17/Global_ID/m17_global_ID_37.csv"
path_trace_m17 <- list.files(path ="~cchen2/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/matlab_analysis/m17/trace_days/", pattern = "*.xlsx", full.names = T)[c(3,7)]
t_stim_m17 <- list(t_stim_m17_d3 = c(437), t_stim_m17_d7 = c(157) )
dat_global_ID_m17 <- cc_globalID_fun(Global_ID_m17, path_trace_m17, t_stim_m17)

## for m18
Global_ID_m18 <- "~cchen2/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/matlab_analysis/m18/Global_ID/m18_global_ID_37.csv"
path_trace_m18 <- list.files(path ="~cchen2/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/matlab_analysis/m18/trace_days/", pattern = "*.xlsx", full.names = T)[c(3,7)]
t_stim_m18 <- list(t_stim_m18_d3 = c(784), t_stim_m18_d7 = c(124))
dat_global_ID_m18 <- cc_globalID_fun(Global_ID_m18, path_trace_m18, t_stim_m18)

## for m855
Global_ID_m855 <- "~cchen2/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/matlab_analysis/m855/Global_ID/m855_global_ID_37.csv"
path_trace_m855 <- list.files(path ="~cchen2/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/matlab_analysis/m855/trace_days/", pattern = "*.xlsx", full.names = T)[c(3,7)]
t_stim_m855 <- list(t_stim_m855_d3 = c(392) *2, t_stim_m855_d7 = c(41)*2)
dat_global_ID_m855 <- cc_globalID_fun(Global_ID_m855, path_trace_m855, t_stim_m855)

## combine data and do k-means analysis-----

dat_cell_trace <- vector(mode = "list", 2)
for (i in 1:2) {
  dat_cell_trace[[i]] <- cbind(dat_global_ID_m3[[i]], dat_global_ID_m7[[i]], dat_global_ID_m17[[i]], dat_global_ID_m18[[i]], dat_global_ID_m855[[i]])
  
}


dat_cell_trace_re <- vector(mode = "list", 2)

mouse_ID <- c("m3", "m7", "m17", "m18", "m855")
for (i in 1:length(dat_cell_trace)) {
  dat_cell_trace_d <- dat_cell_trace[[i]]
  cell_id <- sprintf("Cell%s",seq(1:ncol(dat_cell_trace_d)))
  colnames(dat_cell_trace_d)<- cell_id
  
  dat_cell_trace_d_cluster <- unname(kmeans(t(dat_cell_trace_d), centers = 3, nstart = 20)[[1]])
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
  
  rep_time <- c(ncol(dat_global_ID_m3[[i]]), ncol(dat_global_ID_m7[[i]]), ncol(dat_global_ID_m17[[i]]), ncol(dat_global_ID_m18[[i]]), ncol(dat_global_ID_m855[[i]]) )
  dat_cell_trace_d_re <- mutate(dat_cell_trace_d_re, Group= factor(Group, levels = c("Excited", "Neutral", "Inhibited"))) %>% 
    mutate(., ID = rep(mouse_ID, rep_time*length(stim_time)) )
  
  dat_cell_trace_re[[i]] <- dat_cell_trace_d_re
}

## heat plot of d3 and d7----
## set the range of z score during all days
score_range <- rbind(dat_cell_trace_re[[1]], dat_cell_trace_re[[2]]) %>% 
  .$value %>% 
  range()
group_day <- c("Pre", "Test")

dat_trace_sta <- dat_cell_trace_re[[1]] %>% 
  ddply(., .(variable, Group), summarise,mean=mean(value), sum=sum(value)) %>% 
  arrange(mean)

for (i in c(1: 2)) {
  dat_trace <- dat_cell_trace_re[[i]]
  
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


