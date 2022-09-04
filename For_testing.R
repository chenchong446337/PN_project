## the script for testing different ideas

## plot the average anti trace-----
path_trace_m3 <- "~cchen2/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/matlab_analysis/m3/trace_days/m3_trace_d2.xlsx"
t_stim_m3 <- read.xlsx(path_trace_m3, colNames = F, rowNames = F) %>%
  ncol() %>%
  c(100:(.-200)) %>%
  sample(., n_trial)
dat_trace_m3 <- c_miniscope_matlab(path_trace_m3, t_stim_m3)

path_trace_m7 <- "~cchen2/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/matlab_analysis/m7/trace_days/m7_trace_d3.xlsx"
t_stim_m7 <- read.xlsx(path_trace_m7, colNames = F, rowNames = F) %>%
  ncol() %>%
  c(100:(.-200)) %>%
  sample(., n_trial)
dat_trace_m7 <- c_miniscope_matlab(path_trace_m7, t_stim_m7)


path_trace_m17 <- "~cchen2/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/matlab_analysis/m17/trace_days/m17_trace_d3.xlsx"
t_stim_m17 <- read.xlsx(path_trace_m17, colNames = F, rowNames = F) %>%
  ncol() %>%
  c(100:(.-200)) %>%
  sample(., n_trial)
dat_trace_m17 <- c_miniscope_matlab(path_trace_m17, t_stim_m17)


path_trace_m18 <- "~cchen2/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/matlab_analysis/m18/trace_days/m18_trace_d3.xlsx"
t_stim_m18 <- read.xlsx(path_trace_m18, colNames = F, rowNames = F) %>%
  ncol() %>%
  c(100:(.-200)) %>%
  sample(., n_trial)
dat_trace_m18 <- c_miniscope_matlab(path_trace_m18, t_stim_m18)


dat_cell_trace <- cbind(dat_trace_m3, dat_trace_m7, dat_trace_m17, dat_trace_m18)


# k-menas clustering
mouse_ID <- c("m3", "m7", "m17", "m18")

cell_id <- sprintf("Cell%s",seq(1:ncol(dat_cell_trace)))
colnames(dat_cell_trace)<- cell_id
  
dat_cell_trace_cluster <- unname(kmeans(t(dat_cell_trace), centers = 3, nstart = 20)[[1]])
stim_time<- seq(-2, 6.5, by=0.5)
  #stim_time<- seq(-5, 8.5, by=0.5)
  
dat_cell_trace_re<- as.data.frame(dat_cell_trace) %>% 
    mutate(., Time= stim_time) %>% 
    melt(., id.vars='Time') %>% 
    mutate(., Group= rep(dat_cell_trace_cluster, each=length(stim_time)) )
  
## sort data by the value in each group
dat_cell_sort <- as.numeric(names(sort(tapply(dat_cell_trace_re$value, dat_cell_trace_re$Group, mean), decreasing = T)))
  
dat_cell_trace_re$Group[dat_cell_trace_re$Group == dat_cell_sort[1]] ="Excited"
dat_cell_trace_re$Group[dat_cell_trace_re$Group == dat_cell_sort[2]] ="Neutral"
dat_cell_trace_re$Group[dat_cell_trace_re$Group == dat_cell_sort[3]] ="Inhibited"
  
rep_time <- c(ncol(dat_trace_m3), ncol(dat_trace_m7), ncol(dat_trace_m17), ncol(dat_trace_m18))
dat_cell_trace_re <- mutate(dat_cell_trace_re, Group= factor(Group, levels = c("Excited", "Neutral", "Inhibited"))) %>% 
    mutate(., ID = rep(mouse_ID, rep_time*length(stim_time)) )
  

## plot trace by group
dat_cell_trace_sta <- ddply(dat_cell_trace_re, .(ID,Time, Group), summarise,value=mean(value)) %>% 
  mutate(Group = factor(Group, levels = c("Neutral", "Excited", "Inhibited")))

p_trace <- ddply(dat_cell_trace_sta, .(Time, Group), summarise,n=length(value),mean=mean(value),sd=sd(value),se=sd(value)/sqrt(length(value))) %>% 
  ggplot(., aes(Time, mean, colour=Group))+
  geom_line()+
  geom_ribbon(aes(ymin=mean-se, ymax=mean+se, fill=Group), alpha=0.1, linetype=0)+
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


## make a function for global ID analysis-----
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


##calculate the mean activity of each cells for day2-3 and day5-6
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

## compare crossing back on d7
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
## Extract trace from each mice------
# for m3
path_trace_m3 <- as.list(list.files(path ="~cchen2/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/matlab_analysis/m3/trace_days/", pattern = "*.xlsx", full.names = T ))
path_trace_m3[[8]] <- path_trace_m3[[7]]
t_stim_m3_d1 <- c(214) 
t_stim_m3_d2 <- c(274)
t_stim_m3_d3 <- c(1885) 
t_stim_m3_d4 <- c(1813)
t_stim_m3_d5 <- c(758) 
t_stim_m3_d6 <- c(132) 
t_stim_m3_d7 <- c(346) ## use first two 346, 1206, 1756, 2288, 2970, 3262, 3942
t_stim_m3_d8 <- c(1003)
t_stim_m3 <- list(t_stim_m3_d1, t_stim_m3_d2, t_stim_m3_d3, t_stim_m3_d4, t_stim_m3_d5, t_stim_m3_d6, t_stim_m3_d7, t_stim_m3_d8)

dat_trace_m3 <- mapply(c_miniscope_matlab, path_trace_m3, t_stim_m3, SIMPLIFY = F)

# for m7
path_trace_m7 <- as.list(list.files(path ="~cchen2/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/matlab_analysis/m7/trace_days/", pattern = "*.xlsx", full.names = T ))
path_trace_m7[[8]] <- path_trace_m7[[7]]

t_stim_m7_d1 <- c(605) 
t_stim_m7_d2 <- c(798)
t_stim_m7_d3 <- c(360) 
t_stim_m7_d4 <- c(3236)
t_stim_m7_d5 <- c(317) 
t_stim_m7_d6 <- c(205)
t_stim_m7_d7 <- c(236) ## use first two 236, 1914,2777,3334,4091
t_stim_m7_d8 <- c(1632)

t_stim_m7 <- list(t_stim_m7_d1, t_stim_m7_d2, t_stim_m7_d3, t_stim_m7_d4, t_stim_m7_d5, t_stim_m7_d6, t_stim_m7_d7, t_stim_m7_d8)

dat_trace_m7 <- mapply(c_miniscope_matlab, path_trace_m7, t_stim_m7,SIMPLIFY = F)

# for m17
path_trace_m17 <- as.list(list.files(path ="~cchen2/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/matlab_analysis/m17/trace_days/", pattern = "*.xlsx", full.names = T ))
path_trace_m17[[8]] <- path_trace_m17[[7]]

t_stim_m17_d1 <- c(815) 
t_stim_m17_d2 <- c(1103)
t_stim_m17_d3 <- c(437) 
t_stim_m17_d4 <- c(682)
t_stim_m17_d5 <- c(1114) 
t_stim_m17_d6 <- c(617) 
t_stim_m17_d7 <- c(157) # 157, 2536, 2961
t_stim_m17_d8 <- c(2113)

t_stim_m17 <- list(t_stim_m17_d1, t_stim_m17_d2, t_stim_m17_d3, t_stim_m17_d4, t_stim_m17_d5, t_stim_m17_d6, t_stim_m17_d7, t_stim_m17_d8)

dat_trace_m17 <- mapply(c_miniscope_matlab, path_trace_m17, t_stim_m17, SIMPLIFY = F)

# for m18
path_trace_m18 <- as.list(list.files(path ="~cchen2/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/matlab_analysis/m18/trace_days/", pattern = "*.xlsx", full.names = T ))
path_trace_m18[[8]] <- path_trace_m18[[7]]

t_stim_m18_d1 <- c(490) 
t_stim_m18_d2 <- c(236)
t_stim_m18_d3 <- c(784) 
t_stim_m18_d4 <- c(493)
t_stim_m18_d5 <- c(195) 
t_stim_m18_d6 <- c(467) 
t_stim_m18_d7 <- c(124) ## use the first two of 124,2238,3000
t_stim_m18_d8 <- c(1923)
t_stim_m18 <- list(t_stim_m18_d1, t_stim_m18_d2, t_stim_m18_d3, t_stim_m18_d4, t_stim_m18_d5, t_stim_m18_d6, t_stim_m18_d7, t_stim_m18_d8)

dat_trace_m18 <- mapply(c_miniscope_matlab, path_trace_m18, t_stim_m18, SIMPLIFY = F)

## combine data and do k-means analysis-----

dat_cell_trace <- vector(mode = "list", 8)
for (i in 1:8) {
  dat_cell_trace[[i]] <- cbind(dat_trace_m3[[i]], dat_trace_m7[[i]], dat_trace_m17[[i]], dat_trace_m18[[i]] )
  
}


dat_cell_trace_re <- vector(mode = "list", 8)

# k-menas clustering
mouse_ID <- c("m3", "m7", "m17", "m18")
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
  
  rep_time <- c(ncol(dat_trace_m3[[i]]), ncol(dat_trace_m7[[i]]), ncol(dat_trace_m17[[i]]), ncol(dat_trace_m18[[i]]))
  dat_cell_trace_d_re <- mutate(dat_cell_trace_d_re, Group= factor(Group, levels = c("Excited", "Neutral", "Inhibited"))) %>% 
    mutate(., ID = rep(mouse_ID, rep_time*length(stim_time)) )
  
  dat_cell_trace_re[[i]] <- dat_cell_trace_d_re
}


## heat plot of d3, d6 and d7----
## set the range of z score during all days
score_range <- range(mapply(function(x) x$value, dat_cell_trace_re, SIMPLIFY = T))

for (i in c(3, 6, 7, 8)) {
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

## combine the heat plot
p_heat_com <- plot_grid(p_heat_d3, p_heat_d6, p_heat_d7, nrow = 1)

setwd("~cchen2/Documents/neuroscience/Pn\ project/Figure/PDF/")
cairo_pdf("p_anti_heat.pdf", width = 170/25.6, height = 65/25.6, family = "Arial")
p_heat_com
dev.off()

## plot trace by group-----
dat_cell_trace_pre_sta <- ddply(dat_cell_trace_re[[2]], .(ID,Time, Group), summarise,value=mean(value)) %>%
  rbind(.,ddply(dat_cell_trace_re[[3]], .(ID,Time, Group), summarise,value=mean(value)) ) %>%
  ddply(., .(ID, Time, Group),summarise,value=mean(value) )

dat_cell_trace_con_sta <- ddply(dat_cell_trace_re[[5]], .(ID,Time, Group), summarise,value=mean(value)) %>%
  rbind(.,ddply(dat_cell_trace_re[[6]], .(ID,Time, Group), summarise,value=mean(value)) ) %>%
  ddply(., .(ID, Time, Group),summarise,value=mean(value) )

dat_cell_trace_test_sta <- ddply(dat_cell_trace_re[[7]], .(ID,Time, Group), summarise,value=mean(value))

dat_cell_trace_pain_sta <- ddply(dat_cell_trace_re[[8]], .(ID,Time, Group), summarise,value=mean(value))

dat_cell_trace_sta <- ddply(dat_cell_trace_pre_sta, .(Time, Group), summarise,n=length(value),mean=mean(value),sd=sd(value),se=sd(value)/sqrt(length(value))) %>%
  rbind(., ddply(dat_cell_trace_con_sta, .(Time, Group), summarise,n=length(value),mean=mean(value),sd=sd(value),se=sd(value)/sqrt(length(value)))) %>%
  rbind(., ddply(dat_cell_trace_test_sta, .(Time, Group), summarise,n=length(value),mean=mean(value),sd=sd(value),se=sd(value)/sqrt(length(value)))) %>%
  rbind(., ddply(dat_cell_trace_pain_sta, .(Time, Group), summarise,n=length(value),mean=mean(value),sd=sd(value),se=sd(value)/sqrt(length(value)))) %>%
  mutate(., Day = rep(c("Pre", "Cond.", "Test", "Pain"),each= length(stim_time)*3))

dat_cell_trace_sta$Day <- factor(dat_cell_trace_sta$Day, levels = c("Pre", "Cond.", "Test", "Pain"))
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
  theme(legend.title = element_blank(), legend.position = 'none')

setwd("~cchen2/Documents/neuroscience/Pn\ project/Figure/PDF/")
cairo_pdf("p_trace.pdf", width = 165/25.6, height = 65/25.6, family = "Arial")
p_trace
dev.off()

## plot E and I points to show correlation-----
## points to show the correlation
dat_cell_cor <- NULL
cell_area_day <- c("Rm","Pre", "Pre", "Rm","Cond.","Cond.", "Test", "Pain")

for (i in c(3, 6, 7, 8)) {
  dat_trace <- subset(dat_cell_trace_re[[i]], dat_cell_trace_re[[i]]$ID=="m3") %>%
    subset(., .$Group != "Neutral") %>%
    ddply(., .(Time, Group), summarise, value=mean(value)) %>%
    dcast(., Time~Group) %>%
    mutate(., Day = cell_area_day[i])
  dat_cell_cor <- rbind(dat_cell_cor, dat_trace)
  
}
dat_cell_cor$Day <- factor(dat_cell_cor$Day, levels = c("Pre", "Cond.", "Test", "Pain"))

p_trace_cor <- ggplot(dat_cell_cor, aes(Excited, Inhibited, colour=Day))+
  geom_point()+
  geom_smooth(method = "lm", lwd=0.8)+
  scale_colour_manual(values=c("seagreen", "indianred", "deepskyblue4", "Red"))+  
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
cairo_pdf("p_trace_cor.pdf", width = 60/25.6, height = 65/25.6, family = "Arial")
p_trace_cor
dev.off()
## correlation for each mice
dat_cell_trace_cor <- NULL
cell_area_day <- c("Rm","Pre", "Pre", "Rm","Cond.","Cond.", "Test", "Pain")
for (i in c(2,3,5,6,7, 8)){
  dat_trace_cor <- dat_cell_trace_re[[i]] %>%
    ddply(., .(ID, Time, Group), summarise, value=mean(value)) %>%
    subset(., .$Group!="Neutral") %>%
    dcast(., Time + ID ~ Group) %>%
    ddply(., .(ID), summarise, "corr" = cor(Excited, Inhibited, method = "spearman")) %>%
    mutate(., Day = cell_area_day[i])
  dat_cell_trace_cor <- rbind(dat_cell_trace_cor, dat_trace_cor)
}

dat_cor_value <- ddply(dat_cell_trace_cor, .(ID, Day), summarise, "corr" = mean(corr)) %>%
  mutate(Day = factor(Day, levels = c("Pre", "Cond.", "Test", "Pain")))

p_cor <-ggplot(dat_cor_value, aes(Day, corr, fill=Day))+
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(width = 0.2, shape=1)+
  scale_fill_manual(values=c("seagreen", "indianred", "deepskyblue4", "Red"))+
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
cairo_pdf("p_cor.pdf", width = 40/25.6, height = 65/25.6, family = "Arial")
p_cor
dev.off() 


## calculate the sum of active and inhibited trace----
dat_cell_trace_sum <- rbind(dat_cell_trace_pre_sta, dat_cell_trace_con_sta, dat_cell_trace_test_sta,dat_cell_trace_pain_sta ) %>%
  subset(., .$Group!="Neutral") %>%
  mutate(Day = rep(c("Pre", "Cond.", "Test", "Pain"), each= length(stim_time)*2*length(mouse_ID))) %>%
  ddply(., .(ID,Time, Day), summarise, value=sum(value)) %>%
  ddply(., .(Time, Day), summarise,n=length(value),mean=mean(value),sd=sd(value),se=sd(value)/sqrt(length(value))) %>%
  mutate(Day = factor(Day, levels = c("Pre", "Cond.", "Test", "Pain")))

p_trace_sum <- ggplot(dat_cell_trace_sum, aes(Time, mean, colour=Day))+
  geom_line()+
  geom_ribbon(aes(ymin=mean-se, ymax=mean+se, fill=Day), alpha=0.1, linetype=0)+
  scale_colour_manual(values=c("seagreen", "indianred", "deepskyblue4", "red"))+
  scale_fill_manual(values=c("seagreen", "indianred", "deepskyblue4", "red"))+
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
cairo_pdf("p_trace_sum.pdf", width = 60/25.6, height = 70/25.6, family = "Arial")
p_trace_sum
dev.off()


## for EI change analysis-----
dat_cell_area <- c()
cell_area_day <- c("Rm","Pre", "Pre", "Rm","Cond.","Cond.", "Test", "Pain")
for (i in c(2,3,5,6,7, 8)){
  rep_time <- c(ncol(dat_trace_m3[[i]]), ncol(dat_trace_m7[[i]]), ncol(dat_trace_m17[[i]]), ncol(dat_trace_m18[[i]]))
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

dat_cell_area$Day <- factor(dat_cell_area$Day, levels = c("Pre", "Cond.","Test", "Pain"))
dat_cell_area_sta <-  ddply(dat_cell_area, .(ID, Day), summarise, value=mean(ratio, na.rm = T))

p_EI_ratio <- ggplot(dat_cell_area_sta, aes(Day, value, fill=Day))+
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(width = 0.2, shape=1)+
  scale_fill_manual(values=c("seagreen", "indianred", "deepskyblue4", "Red"))+
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

setwd("~cchen2/Documents/neuroscience/Pn\ project/Figure/PDF/")
cairo_pdf("p_EI_ratio.pdf", width = 40/25.6, height = 65/25.6, family = "Arial")
p_EI_ratio
dev.off() 
## plot the sum of E-I----
dat_cell_sum_sta <-  ddply(dat_cell_area, .(ID, Day), summarise, "sum"=mean(sum, na.rm = T))


p_EI_sum<- ggplot(dat_cell_sum_sta, aes(Day, sum, fill=Day))+
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(width = 0.2, shape=1)+
  scale_fill_manual(values=c("seagreen", "indianred", "deepskyblue4", "red"))+
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
  scale_fill_manual(values=c("seagreen", "indianred", "deepskyblue4", "red"))+
  labs(x="", y="AUC (z-score)")+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  theme(legend.title = element_blank(), legend.position = "top")

setwd("~cchen2/Documents/neuroscience/Pn\ project/Figure/PDF/")
cairo_pdf("p_EIS_change.pdf", width = 70/25.6, height = 72/25.6, family = "Arial")
p_EIS_change
dev.off()

## for extract name----
str_split(x, boundary("word"), simplify = T)
str_match(y, "^m[\\d]")

## firing frequency-----
cc_firing_fun <- function(path_trace, path_peak){
  dat_trace <- read.xlsx(path_trace, colNames = F, rowNames = F)
  cell_valid <- dat_trace$X1
  
  # total time
  n_frame <- ncol(dat_trace)-1
  # 30 is the acquiring frequency
  t_period <- n_frame/30 # s
  dat_peak <- read.csv(path_peak) %>% 
    unlist() %>% 
    na.omit()
  
  ctrl_rang <- range(t_stim_m3_d1-40, t_stim_m3_d1)
  test_rang <- range(t_stim_m3_d1, t_stim_m3_d1 +80)
  
  ctrl_freq <- length(dat_peak[dat_peak>= ctrl_rang[1] & dat_peak < ctrl_rang[2]])/2 # 2s before crosing
  test_freq <- length(dat_peak[dat_peak>= test_rang[1] & dat_peak < test_rang[2]])/4 # 4s after crossing
  
  
   
  ## choose the valide cell number
  dat_peak <- dat_peak[, cell_valid==1]
  
  cell_rate <- apply(dat_peak, 2, function(x) length(which(!is.na(x)))/t_period)
  cell_rate_mean <- mean(cell_rate)
  cell_rate_max <- max(cell_rate)
  return(c(cell_rate_mean, cell_rate_max))
}


## Firing frequency, test all days----
cell_area_day <-  str_c('D', 1:7)
dat_firing <- rbind(dat_rate_m3, dat_rate_m7, dat_rate_m17, dat_rate_m18) %>% 
  as_tibble() %>% 
  rename(Freq_before = V1, Freq_after = V2) %>% 
  mutate(Day = rep(str_c("D", 1:7), 4)) %>% 
  mutate(ID = rep(mouse_ID, each=7)) %>% 
  mutate(Group = rep(cell_area_day, 4)) %>% 
  # filter(Group != "Rm") %>% 
  gather(variable, value, -ID, -Day, -Group) %>% 
  ddply(., .(ID, Group, variable, Day), summarise, value=mean(value)) %>% 
  ddply(., .(ID, Group, variable), summarise, value=mean(value)) %>% 
  mutate( variable=factor(variable, levels = c("Freq_before", "Freq_after")))


## plot the whole trace to show if there is something interesting
ID_trace <- path_trace_m3[[7]]
ID_trace <- path_exp_nal[[4]]
dat_trace <- read.xlsx(ID_trace, colNames = F, rowNames = F) %>% 
  filter(., X1 ==1) %>% 
  select(., -c(1,2)) %>% 
  t() %>% 
  as.data.frame()
rownames(dat_trace) <- NULL
colnames(dat_trace)<- NULL

colnames(dat_trace) <- str_c("Cell", 1:ncol(dat_trace))

## filter the trace
n = 10
length_row <- nrow(dat_trace)
dat_trace_filt <- 
  aggregate(dat_trace,list(rep(1:(length_row%/%n+1),each=n,len=length_row)),mean)[-1]

dat_trace_re <- dat_trace_filt %>% 
  as_tibble() %>% 
  mutate(Time = seq(0, by=0.5, length.out = nrow(dat_trace_filt))) %>%
  gather(variable, value, -Time)


ggplot(dat_trace_re, aes(Time, variable,fill= value))+ 
  geom_raster() +
  #facet_grid(rows = vars(Group), scales = "free_y")+
  labs(x="Time relative to crossing (s)", y="Number of cells")+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))

## plot the diff betwen head and tail to check the behavior----
path_ctrl <- str_c("~cchen2/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/wt_ctrl_07042020/wt_ctrl/T",1)
path_ctrl_path <- list.files(path_ctrl,pattern = ".csv", full.names = T ) 
path_ID <- path_ctrl_path[1]
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

#colnames(dat_anti)<- c("Frame", "Head_x","Head_y","likehood", 'Tail_x',"Tail_y", "likehood_tail")
colnames(dat_anti)<- c("Frame", "Head_x","Head_y", "Prob", "Tail_x", "Tail_y")

## remove the frame for the first few frame
dat_anti <- dat_anti[seq(1, nrow(dat_anti), frame_rate/10), ]
dat_anti <- subset(dat_anti, Prob >0.5 | Frame> 60)

## reduce the frame rate to 5hz to mimimize the varaition 
##n<- 2
## dat_anti<- aggregate(dat_anti,list(rep(1:(nrow(dat_anti)%/%n+1),each=n,len=nrow(dat_anti))),mean)[-1]


## change the dim from pix to mm (the dim of one chammer is 165mm)
# value 387 for anti experiment 20200309
# value 447 for miniscope 20191209
length_pix <- 366
r_pix_mm <- 165 / length_pix
dat_anti$Head_x <- abs(dat_anti$Head_x * r_pix_mm)
dat_anti$Head_y <- abs(dat_anti$Head_y * r_pix_mm)


dat_anti$Tail_x <- abs(dat_anti$Tail_x * r_pix_mm)
dat_anti$Tail_y <- abs(dat_anti$Tail_y * r_pix_mm)

chamber_div <- 165 ## the lenght of two chambels and the gap between

## assign each frame in hot as 1 (48) and nor plate(30) as 2
dat_anti$chamber <- ifelse(dat_anti$Head_y < chamber_div, "Hot","Nor")
ratio_time <- unname(by(dat_anti[,3], dat_anti$chamber, length)[2]/nrow(dat_anti))

p_anti <- ggplot(dat_anti, aes(Head_y, Head_x, color=chamber))+
  geom_point(shape=1)+
  theme_void()+
  theme(legend.position = "none")

ggplot(dat_anti, aes(x=Head_y, y=Head_x) ) +
  geom_hex() +
  scale_fill_viridis()+
  theme_void()


ggplot(dat_anti, aes(x=Head_y, y=Head_x) ) +
  stat_density_2d(aes(fill = ..density..), geom = "raster", contour = FALSE) +
  scale_fill_distiller(palette=4, direction=-1) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  theme(
    legend.position='none'
  )


dat_anti <- dat_anti %>% 
  mutate(dis = abs(sqrt((Head_x-Tail_x)^2 + (Head_y - Tail_y)^2)))

p_dis <- ggplot(dat_anti, aes(Frame, dis))+
  geom_line()
  
dat_anti_optic %>% 
  filter(variable=="latency2") %>%  
  
p_ratio <- dat_anti_optic_filter %>% 
  filter(variable=="ratio" ) %>% 
  ggplot(., aes(Group, change, color=Group))+
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(position=position_jitterdodge(), shape=1)+
  scale_color_manual(values=c("darkcyan", "indianred"))+
  labs(x="", y="Time spend in 2nd plate (%)")

p_latency2 <- dat_anti_optic_filter %>% 
  filter(variable=="latency2" ) %>% 
  ggplot(., aes(Group, change, color=Group))+
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(position=position_jitterdodge(), shape=1)+
  scale_color_manual(values=c("darkcyan", "indianred"))+
  labs(x="", y="Time spend in 2nd plate (%)")






dat_inhibit_manual <- read.xlsx("~cchen2/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/anti_opto_behavior/08252020/wt_optic_08252020_manualy.xlsx") %>% 
  mutate('Rearing'= (First_rearing - Frame_1st_crossing)/10, 'Licking'= (First_licking -Frame_1st_crossing)/10 , 'Jumping'= (Jump -Frame_1st_crossing)/10) %>% 
  pivot_longer(-c(Group, ID), names_to = "variable", values_to = 'value') %>% 
  mutate(Group=factor(Group, levels=c("EYFP", "Nphr", "Chr2")))

dat_inhibit_manual_sta <- dat_inhibit_manual %>% 
  ddply(., .( variable, Group), summarise,n=length(value),mean=mean(value),sd=sd(value),se=sd(value)/sqrt(length(value)))

p_licking <- dat_inhibit_manual %>% 
  filter(variable=="Licking") %>% 
  ggplot(., aes(Group, value, color=Group))+
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(position=position_jitterdodge(), shape=1)+
  scale_color_manual(values=c("springgreen4", "orange2"))+
  labs(x="", y="Latency to 1st licking (s)")+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  scale_y_continuous(limits = c(0, 60), expand = c(0,0))+
  theme(legend.position = "none")+
  annotate(x=c(1,1,2,2), y=c(55,56,56,55),"path")+
  annotate("text",x=1.5,y=56, label="**", size=5)

p_licking_test <- dat_inhibit_manual %>% 
  filter(variable=="Licking") %>% 
  wilcox.test(value~Group, .)

p_1st_rearing <- dat_inhibit_manual %>% 
  filter(variable=="Rearing") %>% 
  ggplot(., aes(Group, value, color=Group))+
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(position=position_jitterdodge(), shape=1)+
  scale_color_manual(values=c("springgreen4", "orange2"))+
  labs(x="", y="Latency to first rearing (s)")+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  theme(legend.position = "none")

p_1st_rearing_test <- dat_inhibit_manual %>% 
  filter(variable=="Rearing") %>% 
  wilcox.test(value~Group,.)

p_Jump <- dat_inhibit_manual %>% 
  filter(variable=="Jumping") %>% 
  ggplot(., aes(Group, value, color=Group))+
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(position=position_jitterdodge(), shape=1)+
  scale_color_manual(values=c("springgreen4", "orange2"))+
  labs(x="", y="Latency to 1st jumping ")+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  theme(legend.position = "none")+
  annotate(x=c(1,1,2,2), y=c(185,187,187,185),"path")+
  annotate("text",x=1.5,y=187, label="*", size=5)

p_jumping_test<- dat_inhibit_manual %>% 
  filter(variable=="Jumping") %>% 
  wilcox.test(value~Group,.)


p_num_rearing <- dat_inhibit_manual %>% 
  filter(variable=="num_rearing") %>% 
  ggplot(., aes(Group, value, color=Group))+
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(position=position_jitterdodge(), shape=1)+
  scale_color_manual(values=c("springgreen4", "orange2"))+
  labs(x="", y="# rearing / min ")+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  theme(legend.position = "none")

## for pn_cpp_0104----
cc_anti_cpp <- function(path_ID, frame_rate){
  
  dat_anti <- read.csv(path_ID, header = F, stringsAsFactors = F) %>% 
    .[-c(1:3), ] %>% 
    apply(., 2, as.numeric) %>%   ## only select the position of head
    as.data.frame() %>% 
    select(str_c("V", c(1:4,8:9))) 
  
  colnames(dat_anti)<- c("Frame", "Head_x","Head_y", "Prob",  "reference_x", "reference_y")
  
  ## remove the frame for the first few frame
  ## the length and width of the CPP box (426.6,416 )
  limit_x <- mean(dat_anti$reference_x)
  limit_y <- mean(dat_anti$reference_y)
  
  dat_anti <- dat_anti[seq(1, nrow(dat_anti), frame_rate/10), ] %>% 
    subset(., Prob >0.5 | Frame> 60) %>% 
    mutate(Head_x = ifelse(Head_x < limit_x , limit_x, Head_x), Head_y = ifelse(Head_y < (limit_y - 416), (limit_y-416), Head_y)) %>% 
    mutate(Head_x = ifelse(Head_x > limit_x + 426.6 , limit_x +426.6, Head_x), Head_y = ifelse(Head_y > (limit_y + 430), (limit_y+430), Head_y)) %>% 
    mutate(reference_y = mean(reference_y)) %>% 
    mutate(chamber = ifelse(Head_y > (reference_y +12), "right", "left"))
  
  ## for the right side
  ratio_time <- unname(by(dat_anti[,3], dat_anti$chamber, length)[1]/nrow(dat_anti))
  
  p_anti <- ggplot(dat_anti, aes(Head_y, Head_x,colour=chamber))+
    geom_point(shape=1, size=1)+
    theme_void()+
    theme(legend.position = "none")
  
  return(ratio_time)
}

### plot anti behavior of WT and cond. mice for every day-----

dat_anti_ctrl <- vector(mode = "list", 7)
dat_anti_con <- vector(mode = "list", 7)
length_pix <- c(363,362,363,364,362,364,365)
day_group <- c("Rm","Rm", "Pre", "D4","D5","D6", "Post")
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
  mutate(Day=factor(Day, levels = c("Pre", "D4", "D5", "D6", "Post")), Group=factor(Group, levels = c("Ctrl", "Cond.")))

dat_anti_wt_sta <- ddply(dat_anti_wt,.(variable,Day, Group), summarise,n=length(value),mean=mean(value),sd=sd(value),se=sd(value)/sqrt(length(value))) 

p_ratio <- dat_anti_wt_sta %>% 
  filter(variable=="ratio" ) %>% 
  ggplot(., aes(Day, mean, group= Group, colour = Group))+
  geom_point( aes(shape = Group), size = 2)+ 
  geom_errorbar(aes(ymin=mean-se,ymax=mean+se), width = .2)+
  geom_line( )+
  scale_shape_manual(values=c(20,18))+
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
  
p_latency1 <- dat_anti_wt_sta %>% 
  filter(variable=="latency1" ) %>% 
  ggplot(., aes(Day, mean, group= Group, colour = Group))+
  geom_point( aes(shape = Group), size = 2)+ 
  geom_errorbar(aes(ymin=mean-se,ymax=mean+se), width = .2)+
  geom_line( )+
  scale_shape_manual(values=c(20,18))+
  scale_color_manual(values=c("darkcyan", "indianred"))+
  labs(x="", y="Latency of 1st crossing (s)")+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  theme(legend.position = "none")

p_latency2 <- dat_anti_wt_sta %>% 
  filter(variable=="latency2" ) %>% 
  ggplot(., aes(Day, mean, group= Group, colour = Group))+
  geom_point( aes(shape = Group), size = 2)+ 
  geom_errorbar(aes(ymin=mean-se,ymax=mean+se), width = .2)+
  geom_line( )+
  scale_shape_manual(values=c(20,18))+
  scale_color_manual(values=c("darkcyan", "indianred"))+
  labs(x="", y="Latency of 1st crossing back (s)")+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  theme(legend.position = "none")

p_ratio <- dat_anti_wt %>% 
  filter(variable=="ratio" ) %>% 
  mutate(value=value*100) %>% 
  ggplot(., aes(Day, value, colour=Group))+
  stat_boxplot(position =position_dodge(width = 0.8),geom ='errorbar', width = 0.2)+
  geom_boxplot(outlier.shape = NA, position =position_dodge(width = 0.8), width = 0.75)+
  geom_jitter(position=position_jitterdodge(), shape=1, alpha = 0.5)+
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
  theme(legend.title = element_blank())

p_ratio_test <- dat_anti_wt %>% 
  filter(variable=="ratio" & Day=="Test") %>% 
  wilcox.test(value~Group, .)

p_latency1 <- dat_anti_wt %>% 
  filter(variable=="latency1") %>% 
  ggplot(., aes(Day, value, colour=Group))+
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(position=position_jitterdodge(), shape=1)+
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
  theme(legend.position = "none")+
  annotate(x=c(1.8,1.8,2.2,2.2), y=c(100,102,102,100),"path")+
  annotate("text",x=2,y=102, label="*", size=5)

p_latency1_test <- dat_anti_wt %>% 
  filter(variable=="latency1"& Day=="Test") %>% 
  wilcox.test(value~Group, .)

p_latency2 <- dat_anti_wt %>% 
  filter(variable=="latency2") %>% 
  ggplot(., aes(Day, value, colour=Group))+
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(position=position_jitterdodge(), shape=1)+
  scale_color_manual(values=c("darkcyan", "indianred"))+
  labs(x="", y="Latency (s)")+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  scale_y_continuous(limits = c(0, 90))+
  theme(legend.position = "none")+
  annotate(x=c(1.8,1.8,2.2,2.2), y=c(70,72,72,70),"path")+
  annotate("text",x=2,y=72, label="**", size=5)



