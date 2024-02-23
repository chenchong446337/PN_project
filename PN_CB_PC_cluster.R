## CLustering analysis for PC data

back_crossing_frame <- tibble(ID = c("m617", "m625", "m676", "m684", "m687", "m685"), Frame = c(1974, 2588, 1612, 1294, 2210, 2698), 
                              Last_frame= c(3016, 3168, 3106, 3180, 3726, 3640))

c_miniscope_matlab_ft <- function(file_trace) {
  ## import and format the data
  dat_trace1 <- raveio::read_mat(file_trace)
  ID_mouse <- str_extract(file_trace, regex("m\\d+"))
  
  t_crossing_back <- back_crossing_frame %>% 
    dplyr::filter(ID == ID_mouse)
  t_crossing <- c(dat_trace1[[1]][5], t_crossing_back$Frame, t_crossing_back$Last_frame)
  
  cross_ID <- dat_trace1$global.map %>% 
    as_tibble() %>% 
    select(V1, V4, V5) %>% 
    mutate_all(na_if, 0) %>% 
    drop_na() %>% 
    select(V1, V5) %>% 
    pull(V5)
  
  
  dat_stim_trace <- vector(mode = "list", length = length(t_crossing))
  
  
  for (i in seq_along(t_crossing)) {
    dat_trace <- dat_trace1[[3]][[5]] %>% 
      as_tibble() %>% 
      select(all_of(cross_ID)) %>% 
      apply(., 2, scale)
    
    ## number of rows to be binned
    n <- 2 # 0.05*2=0.1
    t1_p <- t_crossing[i]
    if (t1_p < 40){
      column_means = colMeans(dat_trace[1:t1_p, ])
      
      # Replicate these mean values 10 times
      n_reapt <- (40 - t1_p)
      new_rows = matrix(rep(column_means, each = n_reapt), nrow = n_reapt, ncol = ncol(dat_trace))
      
      # Combine the new matrix with the original one
      dat_trace = rbind(new_rows, dat_trace)
      
    } else {
      dat_trace <- dat_trace
    }
    
    
    ## extract cell activity when they cross the border
    #dat_stim1 <- dat_trace[(t1_p-40):(t1_p+140-1),] 
    t1_p <- ifelse(t1_p < 40, 40, t1_p)
    dat_stim1 <- dat_trace[(t1_p-40):(t1_p+40-1),] 
    
    dat_stim <- aggregate(dat_stim1,list(rep(1:(nrow(dat_stim1)%/%n+1),each=n,len=nrow(dat_stim1))),mean)[-1]
    # apply(., 2, scale) %>% 
    # as_tibble() %>% 
    # replace(is.na(.), 0)
    
    
    colnames(dat_stim) <- str_c(ID_mouse,"Cell", 1: ncol(dat_stim))
    dat_stim_trace[[i]] <- dat_stim
  }
  return(dat_stim_trace)
}


mouse_file <- as.list(list.files("~cchen/Documents/neuroscience/Pn\ project/Data_analysis/For\ revision/PC_data", full.names = T))

dat_cell_trace <- mapply(c_miniscope_matlab_ft, mouse_file, SIMPLIFY = F)

## do the clustering analysis

dat_cell_pre_trace <- lapply(dat_cell_trace, function(x) x[[1]]) %>% 
  do.call(cbind,.)

dat_cell_cond_trace <- lapply(dat_cell_trace, function(x) x[[2]]) %>% 
  do.call(cbind,.)

dat_cb_combine <- lapply(dat_cell_trace, function(x) x[[3]]) %>% 
  do.call(cbind,.) %>% 
  rbind(dat_cell_pre_trace, dat_cell_cond_trace,.)

fviz_nbclust(t(dat_cb_combine), kmeans, method = "silhouette")

dat_cell_trace_cluster <- kmeans(t(dat_cb_combine), centers = 2, iter.max = 10)[[1]]


# Heatmap plot
comp_time <- seq(-2, 1.9, by=0.1)
group_day <- c("Crossing", "Crossing_back", "Last_crossing")
score_range <- range(dat_cell_trace)



for(i in c(1:3)){
  cell_order <- lapply(dat_cell_trace, function(x)  x[[1]]) %>% 
    do.call(cbind,.) %>% 
    mutate(Time = comp_time) %>% 
    pivot_longer(-Time) %>% 
    mutate(cluster = map_dbl(name, ~dat_cell_trace_cluster[.x])) %>% 
    mutate(cluster = factor(cluster, levels = c(1,2))) %>% 
    group_by(name, cluster) %>%
    summarise(mean = mean(value), .groups = "drop") %>%
    arrange(cluster, mean)
  
  dat_rect <- lapply(dat_cell_trace, function(x)  x[[1]]) %>% 
    do.call(cbind,.) %>% 
    mutate(Time = comp_time) %>% 
    pivot_longer(-Time) %>% 
    mutate(cluster = map_dbl(name, ~dat_cell_trace_cluster[.x])) %>% 
    mutate(cluster = factor(cluster, levels = c(1, 2))) %>% 
    ddply(., .(cluster), summarise, n = length(name)/40) %>% 
    add_column(xmin = 2, xmax = 2.1) %>% 
    mutate(ymax = cumsum(n) ) %>% 
    mutate(ymin = c(1, lag(ymax)[-1]))
  
  
  p_heat <- lapply(dat_cell_trace, function(x)  x[[i]]) %>% 
    do.call(cbind,.) %>% 
    mutate(Time = comp_time) %>% 
    pivot_longer(-Time) %>% 
    mutate(name = factor(name, levels = cell_order$name)) %>% 
    ggplot(., aes(Time, name,fill= value))+ 
    geom_tile(height=2)+
    scale_fill_gradientn(limits= score_range, colours = c("navy", "white", "red4"), values = rescale(c(score_range[1], 0, score_range[2])))+
    labs(x="", y="")+
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
    theme(legend.position = "none")+
    annotate("rect", xmin = 2, xmax = 2.1, ymin = dat_rect$ymin, ymax = dat_rect$ymax-1, alpha = .2)
  
  assign(str_c("p_heat_", group_day[i]), p_heat)
  
}

p_heat <- plot_grid(p_heat_Crossing, p_heat_Crossing_back, p_heat_Last_crossing, nrow  = 1)

setwd("~cchen/Documents/neuroscience/Pn\ project/Figure/PDF/")
cairo_pdf("p_heat_com_cluster.pdf", width = 120/25.6, height = 60/25.6, family = "Arial")
p_heat
dev.off()    


# heatmap plot
p_cell_firing_corssing <- lapply(dat_cell_trace, function(x)  x[[1]]) %>% 
  do.call(cbind,.) %>% 
  add_column(Group = "Crossing")

p_cell_firing_corssing_back <- lapply(dat_cell_trace, function(x)  x[[2]]) %>% 
  do.call(cbind,.) %>% 
  add_column(Group = "Crossing_back")

p_cell_firing <-  lapply(dat_cell_trace, function(x)  x[[3]]) %>% 
  do.call(cbind,.) %>% 
  add_column(Group = "Last_crossing") %>% 
  rbind(p_cell_firing_corssing, p_cell_firing_corssing_back,.) %>% 
  as_tibble() %>% 
  pivot_longer(-Group) %>% 
  mutate(cluster = map_dbl(name, ~dat_cell_trace_cluster[.x])) %>% 
  mutate(Group = factor(Group, levels = c("Crossing", "Crossing_back", "Last_crossing"))) %>% 
  ddply(.,.(name, Group, cluster), summarise, mean = mean(value)) %>% 
  ggplot(., aes(Group, mean, group = Group, color = Group))+
  geom_violin()+
  #geom_boxplot(outlier.shape = NA)+
  #geom_line(aes(group=name), colour="gray90")+
  geom_jitter(aes(colour = Group, shape = Group),width = 0.2,  size=2, alpha= 0.3)+
  #scale_color_manual(values=c("#8491B4FF", "#00A087FF", "#3C5488FF"))+
  facet_grid(cols = vars(cluster))+
  labs(x="", y="Mean amplitude of Ca2+ activity \nduring first crossing (s.d.)")+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  scale_y_continuous(limits = c(-1, 6))+
  theme(legend.position = 'none')


setwd("~cchen/Documents/neuroscience/Pn\ project/Figure/PDF/")
cairo_pdf("p_cell_firing.pdf", width = 80/25.6, height = 65/25.6, family = "Arial" )
p_cell_firing
dev.off()

## compare PC in two clusters during Pre, Cond, and Post-----
c_miniscope_matlab_ft <- function(file_trace) {
  ## import and format the data
  dat_trace1 <- raveio::read_mat(file_trace)
  ID <- str_extract(file_trace, regex("m\\d+"))
  num_compare <- c(1,4, 5)
  t_crossing <- c(1,4, 5)
  
  cross_ID <- dat_trace1$global.map %>% 
    as_tibble() %>% 
    select(V1, V4, V5) %>% 
    mutate_all(na_if, 0) %>% 
    drop_na()
  
  dat_stim_trace <- vector(mode = "list", length = length(num_compare))
  
  
  for (i in seq_along(num_compare)) {
    global_cell <- pull(cross_ID[,i])
    dat_trace <- dat_trace1$traces[[num_compare[i]]] %>% 
      as_tibble() %>% 
      select(all_of(global_cell)) %>% 
      apply(., 2, scale)
    
    ## number of rows to be binned
    n <- 2 # 0.05*2=0.1
    t1_p <- dat_trace1$crossing[t_crossing[i]]
    
    if (t1_p < 40){
      column_means = colMeans(dat_trace[1:t1_p, ])
      
      # Replicate these mean values 10 times
      n_reapt <- (40 - t1_p)
      new_rows = matrix(rep(column_means, each = n_reapt), nrow = n_reapt, ncol = ncol(dat_trace))
      
      # Combine the new matrix with the original one
      dat_trace = rbind(new_rows, dat_trace)
      
    } else {
      dat_trace <- dat_trace
    }
    
    
    ## extract cell activity when they cross the border
    #dat_stim1 <- dat_trace[(t1_p-40):(t1_p+140-1),] 
    t1_p <- ifelse(t1_p < 40, 40, t1_p)
    dat_stim1 <- dat_trace[(t1_p-39):(t1_p+40),] 
    
    dat_stim <- aggregate(dat_stim1,list(rep(1:(nrow(dat_stim1)%/%n+1),each=n,len=nrow(dat_stim1))),mean)[-1]
    # apply(., 2, scale) %>% 
    # as_tibble() %>% 
    # replace(is.na(.), 0)
    
    
    colnames(dat_stim) <- str_c(ID,"Cell", 1: ncol(dat_stim))
    dat_stim_trace[[i]] <- dat_stim
  }
  return(dat_stim_trace)
}



mouse_file <- as.list(list.files("~cchen/Documents/neuroscience/Pn\ project/Data_analysis/For\ revision/PC_data", full.names = T))

dat_cell_trace <- mapply(c_miniscope_matlab_ft, mouse_file, SIMPLIFY = F)

dat_cell_pre_trace <- lapply(dat_cell_trace, function(x) x[[1]]) %>% 
  do.call(cbind,.)

dat_cell_cond_trace <- lapply(dat_cell_trace, function(x) x[[2]]) %>% 
  do.call(cbind,.)

dat_cb_combine <- lapply(dat_cell_trace, function(x) x[[3]]) %>% 
  do.call(cbind,.) %>% 
  rbind(dat_cell_pre_trace, dat_cell_cond_trace,.)

comp_time <- seq(-2, 1.9, by=0.1)
group_day <- c("Pre", "Cond", "Post")
score_range <- range(dat_cell_trace)

for(i in c(1:3)){
  cell_order <- lapply(dat_cell_trace, function(x) x[[2]]) %>% 
    do.call(cbind,.) %>% 
    mutate(Time = comp_time) %>% 
    pivot_longer(-Time) %>% 
    mutate(cluster = map_dbl(name, ~dat_cell_trace_cluster[.x])) %>% 
    mutate(cluster = factor(cluster, levels = c(1,2))) %>% 
    group_by(name, cluster) %>%
    summarise(mean = mean(value), .groups = "drop") %>%
    arrange(cluster, mean)
  
  dat_rect <- dat_cell_cond_trace %>% 
    mutate(Time = comp_time) %>% 
    pivot_longer(-Time) %>% 
    mutate(cluster = map_dbl(name, ~dat_cell_trace_cluster[.x])) %>% 
    mutate(cluster = factor(cluster, levels = c(1,2))) %>% 
    ddply(., .(cluster), summarise, n = length(name)/40) %>% 
    add_column(xmin = 2, xmax = 2.1) %>% 
    mutate(ymax = cumsum(n) ) %>% 
    mutate(ymin = c(1, lag(ymax)[-1]))
  
  p_heat <- lapply(dat_cell_trace, function(x)  x[[i]]) %>% 
    do.call(cbind,.) %>% 
    mutate(Time = comp_time) %>% 
    pivot_longer(-Time) %>% 
    mutate(name = factor(name, levels = cell_order$name)) %>% 
    ggplot(., aes(Time, name,fill= value))+ 
    geom_tile(height=2)+
    scale_fill_gradientn(limits= score_range, colours = c("navy", "white", "red4"), values = rescale(c(score_range[1], 0, score_range[2])))+
    labs(x="", y="")+
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
    theme(legend.position = "none")+
    annotate("rect", xmin = 2, xmax = 2.1, ymin = dat_rect$ymin, ymax = dat_rect$ymax-1, alpha = .2)
  
  assign(str_c("p_heat_", group_day[i]), p_heat)
  
}

p_heat <- plot_grid(p_heat_Pre, p_heat_Cond, p_heat_Post, nrow  = 1)

setwd("~cchen/Documents/neuroscience/Pn\ project/Figure/PDF/")
cairo_pdf("p_heat_pc.pdf", width = 140/25.6, height = 60/25.6, family = "Arial")
p_heat
dev.off()


## boxplot 

p_combine_cluster <- lapply(dat_cell_trace, function(x) x[[3]]) %>% 
  do.call(cbind,.) %>% 
  rbind(dat_cell_pre_trace, dat_cell_cond_trace,.) %>% 
  mutate(Group = rep(c("Pre", "Cond", "Post"), each = 40)) %>% 
  mutate(Time = rep(comp_time, 3)) %>% 
  pivot_longer(-c(Time, Group)) %>% 
  mutate(cluster = map_dbl(name, ~dat_cell_trace_cluster[.x])) %>% 
  ddply(.,.(name, Group, cluster), summarise, mean = mean(value)) %>% 
  mutate(Group = factor(Group, levels = c("Pre", "Cond", "Post"))) %>% 
  ggplot(., aes(Group, mean, group = Group, color = Group))+
  geom_violin()+
  #geom_boxplot(outlier.shape = NA)+
  #geom_line(aes(group=name), colour="gray90")+
  geom_jitter(aes(colour = Group, shape = Group),width = 0.2,  size=2, alpha= 0.3)+
  scale_color_manual(values=c("#8491B4FF", "#00A087FF", "#3C5488FF"))+
  facet_grid(cols = vars(cluster))+
  labs(x="Time (s)", y="Mean amplitude of Ca2+ activity \nduring first crossing (s.d.)")+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  scale_y_continuous(limits = c(-1, 8))+
  theme(legend.position = 'none')

t_combine_cluster <- lapply(dat_cell_trace, function(x) x[[3]]) %>% 
  do.call(cbind,.) %>% 
  rbind(dat_cell_pre_trace, dat_cell_cond_trace,.) %>% 
  mutate(Group = rep(c("Pre", "Cond", "Post"), each = 40)) %>% 
  mutate(Time = rep(comp_time, 3)) %>% 
  pivot_longer(-c(Time, Group)) %>% 
  mutate(cluster = map_dbl(name, ~dat_cell_trace_cluster[.x])) %>% 
  ddply(.,.(name, Group, cluster), summarise, mean = mean(value)) %>% 
  mutate(Group = factor(Group, levels = c("Pre", "Cond", "Post"))) %>% 
  filter(cluster==2) %>% 
  filter(Group!= "Pre") %>% 
  wilcox.test(mean~Group, .,paired = T)
## plot for each mouse
p_combine_cluster_mouse <- lapply(dat_cell_trace, function(x) x[[3]]) %>% 
  do.call(cbind,.) %>% 
  rbind(dat_cell_pre_trace, dat_cell_cond_trace,.) %>% 
  mutate(Group = rep(c("Pre", "Cond", "Post"), each = 40)) %>% 
  mutate(Time = rep(comp_time, 3)) %>% 
  pivot_longer(-c(Time, Group)) %>% 
  mutate(cluster = map_dbl(name, ~dat_cell_trace_cluster[.x])) %>% 
  ddply(.,.(name, Group, cluster), summarise, mean = mean(value)) %>% 
  mutate(ID = str_extract(name, regex("m\\d+"))) %>% 
  ddply(., .(ID,Group, cluster ), summarise, mean_acti = mean(mean)) %>% 
  mutate(Group = factor(Group, levels = c("Pre", "Cond", "Post"))) %>% 
  ggplot(., aes(Group, mean_acti, group = Group, color = Group))+
  geom_boxplot(outlier.shape = NA)+
  geom_line(aes(group=ID), colour="gray90")+
  geom_jitter(aes(colour = Group, shape = Group),width = 0.2,  size=2, alpha= 0.3)+
  scale_color_manual(values=c("#8491B4FF", "#00A087FF", "#3C5488FF"))+
  facet_grid(cols = vars(cluster))+
  labs(x="", y="Mean amplitude of Ca2+ activity \nduring first crossing (s.d.)")+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  scale_y_continuous(limits = c(-1, 4))+
  theme(legend.position = 'none')


t_combine_cluster_mouse <- lapply(dat_cell_trace, function(x) x[[3]]) %>% 
  do.call(cbind,.) %>% 
  rbind(dat_cell_pre_trace, dat_cell_cond_trace,.) %>% 
  mutate(Group = rep(c("Pre", "Cond", "Post"), each = 40)) %>% 
  mutate(Time = rep(comp_time, 3)) %>% 
  pivot_longer(-c(Time, Group)) %>% 
  mutate(cluster = map_dbl(name, ~dat_cell_trace_cluster[.x])) %>% 
  ddply(.,.(name, Group, cluster), summarise, mean = mean(value)) %>% 
  mutate(ID = str_extract(name, regex("m\\d+"))) %>% 
  ddply(., .(ID,Group, cluster ), summarise, mean_acti = mean(mean)) %>% 
  mutate(Group = factor(Group, levels = c("Pre", "Cond", "Post"))) %>% 
  filter(cluster==1) %>% 
  filter(Group !="Cond") %>% 
  wilcox.test(mean_acti~Group,., paired = T)

summary(t_combine_cluster_mouse)
TukeyHSD(t_combine_cluster_mouse)
setwd("~cchen/Documents/neuroscience/Pn\ project/Figure/PDF/")
cairo_pdf("p_combine_cluster.pdf", width = 80/25.6, height = 65/25.6, family = "Arial" )
p_combine_cluster_mouse
dev.off()

t_combine_cluster <- lapply(dat_cell_trace, function(x) x[[3]]) %>% 
  do.call(cbind,.) %>% 
  rbind(dat_cell_pre_trace, dat_cell_cond_trace,.) %>% 
  mutate(Group = rep(c("Pre", "Cond", "Post"), each = 40)) %>% 
  mutate(Time = rep(comp_time, 3)) %>% 
  pivot_longer(-c(Time, Group)) %>% 
  mutate(cluster = map_dbl(name, ~dat_cell_trace_cluster[.x])) %>% 
  ddply(.,.(name, Group, cluster), summarise, mean = mean(value)) %>% 
  filter(cluster ==2) %>% 
  aov(mean~ Group,.)

summary(t_combine_cluster)
TukeyHSD(t_combine_cluster)  

## amplitude of each neurons in this cluster

c_miniscope_matlab_ft_time <- function(file_trace) {
  ## import and format the data
  dat_trace1 <- raveio::read_mat(file_trace)
  ID <- str_extract(file_trace, regex("m\\d+"))
  num_compare <- c(1,4, 5)
  t_crossing <- c(1,4, 5)
  
  cross_ID <- dat_trace1$global.map %>% 
    as_tibble() %>% 
    select(V1, V4, V5) %>% 
    mutate_all(na_if, 0) %>% 
    drop_na()
  
  dat_spike_time <- list.files("~cchen/Documents/neuroscience/Pn\ project/Data_analysis/For\ revision/Event_time/", pattern = ID, full.names = T)
  
  
  dat_stim_trace <- vector(mode = "list", length = length(num_compare))
  
  
  for (i in seq_along(num_compare)) {
    global_cell <- pull(cross_ID[,i])
    dat_trace <- dat_trace1$traces[[num_compare[i]]] %>% 
      as_tibble() %>% 
      select(all_of(global_cell)) %>% 
      apply(., 2, scale) 
    
    dat_event_time <- read.xlsx(dat_spike_time[i], colNames = F) %>% 
      as_tibble() %>% 
      select(all_of(global_cell)) %>% 
      as.matrix()
    
    # Create a new matrix to hold the result
    dat_event_peak <- matrix(0, nrow = nrow(dat_event_time), ncol = ncol(dat_event_time))
    
    # Identify positions where event_time_matrix is 1
    event_positions <- dat_event_time == 1
    
    # Replace "1"s in result_matrix with corresponding values from data_trace_matrix
    dat_event_peak[event_positions] <- dat_trace[event_positions]
    
    # Optionally, convert result_matrix back to a data frame
    dat_event_peak <- as.data.frame(dat_event_peak)
    
    # result_df now contains the raw values from data_trace where event_time had "1"s, and "0" otherwise
    
    
    ## number of rows to be binned
    t1_p <- dat_trace1$crossing[t_crossing[i]]
    
    if (t1_p < 40){
      column_means = colMeans(dat_event_peak[1:t1_p, ])
      
      # Replicate these mean values 10 times
      n_reapt <- (40 - t1_p)
      new_rows = matrix(rep(column_means, each = n_reapt), nrow = n_reapt, ncol = ncol(dat_event_peak))
      
      # Combine the new matrix with the original one
      dat_event_peak = rbind(new_rows, dat_event_peak)
      
    } else {
      dat_event_peak <- dat_event_peak
    }
    
    
    ## extract cell activity when they cross the border
    #dat_stim1 <- dat_trace[(t1_p-40):(t1_p+140-1),] 
    t1_p <- ifelse(t1_p < 40, 40, t1_p)
    dat_stim <- dat_event_peak[(t1_p-39):(t1_p+40),] 
    colnames(dat_stim) <- str_c(ID,"Cell", 1: ncol(dat_stim))
    dat_stim_trace[[i]] <- dat_stim
  }
  return(dat_stim_trace)
}


dat_cell_trace_event <- mapply(c_miniscope_matlab_ft_time, mouse_file, SIMPLIFY = F)

dat_cell_pre_trace_event <- lapply(dat_cell_trace_event, function(x) x[[1]]) %>% 
  do.call(cbind,.) %>% 
  add_column(Group = "Pre") %>% 
  pivot_longer(-Group) %>% 
  mutate(value = na_if(value, 0)) %>% 
  drop_na()

dat_cell_cond_trace_event <- lapply(dat_cell_trace_event, function(x) x[[2]]) %>% 
  do.call(cbind,.) %>% 
  add_column(Group = "Cond") %>% 
  pivot_longer(-Group) %>% 
  mutate(value = na_if(value, 0)) %>% 
  drop_na()

dat_cb_combine_event <- lapply(dat_cell_trace_event, function(x) x[[3]]) %>% 
  do.call(cbind,.) %>% 
  add_column(Group = "Post") %>% 
  pivot_longer(-Group) %>% 
  mutate(value = na_if(value, 0)) %>% 
  drop_na() %>% 
  rbind(dat_cell_pre_trace_event, dat_cell_cond_trace_event,.)

h_cell_peak <- dat_cb_combine_event %>% 
  mutate(cluster = map_dbl(name, ~dat_cell_trace_cluster[.x])) %>% 
  filter(cluster ==1) %>% 
  mutate(Group = factor(Group, levels = c("Pre", "Cond", "Post"))) %>% 
  ggplot(., aes(x = value, col = Group, fill = Group))+
  geom_density(adjust=1.5, alpha=.4)

## calculate the mean and sd of spikes in Pre group
dat_cell_pre_trace_event_sta <- lapply(dat_cell_trace_event, function(x) x[[1]]) %>% 
  do.call(cbind,.) %>% 
  add_column(Group = "Pre") %>% 
  pivot_longer(-Group) %>% 
  mutate(value = na_if(value, 0)) %>% 
  drop_na() %>% 
  ddply(., .(Group), summarise, mean = mean(value), sd = sd(value))


p_hist_pre <- lapply(dat_cell_trace_event, function(x) x[[1]]) %>% 
  do.call(cbind,.) %>% 
  add_column(Group = "Pre") %>% 
  pivot_longer(-Group) %>% 
  mutate(value = na_if(value, 0)) %>% 
  drop_na() %>% 
  ggplot(.,aes(x=value))+
  geom_histogram()+
  geom_vline(xintercept = dat_cell_pre_trace_event_sta$mean, col="red", linetype = 2 )+
  geom_vline(xintercept = (dat_cell_pre_trace_event_sta$mean + dat_cell_pre_trace_event_sta$sd), col="blue", linetype = 2 )+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  theme(legend.title = element_blank(), legend.position = "none")

setwd("~cchen/Documents/neuroscience/Pn\ project/Figure/PDF/")
cairo_pdf("p_hist_pre.pdf", width = 80/25.6, height = 65/25.6, family = "Arial" )
p_hist_pre
dev.off()

## the mean + 1sd is 3 
p_peak_amp <- dat_cb_combine_event %>% 
  mutate(cluster = map_dbl(name, ~dat_cell_trace_cluster[.x])) %>% 
  mutate(Group = factor(Group, levels = c("Pre", "Cond", "Post"))) %>% 
  mutate(cluster = factor(cluster, levels = c(2,1))) %>% 
  filter(value > 3) %>% 
  ggplot(., aes(Group, value, col = Group))+
  geom_violin()+
  geom_jitter(aes(colour = Group, shape = Group),width = 0.2,  size=2, alpha= 0.3)+
  scale_color_manual(values=c("#8491B4FF", "#00A087FF", "#3C5488FF"))+
  facet_grid(cols = vars(cluster))+
  labs(x="", y="Peak amplitude (s.d)")+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  scale_y_continuous(limits = c(3, 15))+
  theme(legend.title = element_blank(), legend.position = "none")


## for individual mice 
p_peak_amp_mouse <- dat_cb_combine_event %>% 
  mutate(cluster = map_dbl(name, ~dat_cell_trace_cluster[.x])) %>% 
  mutate(Group = factor(Group, levels = c("Pre", "Cond", "Post"))) %>% 
  filter(value > 3) %>% 
  mutate(ID = str_extract(name, regex("m\\d+"))) %>% 
  ddply(.,.(ID, Group, cluster), summarise, mean_amp = mean(value)) %>% 
  ggplot(., aes(Group, mean_amp, col = Group))+
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(aes(colour = Group, shape = Group),width = 0.2,  size=2, alpha= 0.3)+
  scale_color_manual(values=c("#8491B4FF", "#00A087FF", "#3C5488FF"))+
  facet_grid(cols = vars(cluster))+
  labs(x="", y="Peak amplitude (s.d)")+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  scale_y_continuous(limits = c(3, 15))+
  theme(legend.title = element_blank(), legend.position = "none")

## calculate the firing frequency of these events with larger amplitude 
dat_cell_trace_cluster_num <- tibble(name = names(dat_cell_trace_cluster), cluster = dat_cell_trace_cluster) %>% 
  as_tibble() %>% 
  slice(rep(1:n(), each = 3)) %>% 
  add_column(Group = rep(c("Pre", "Cond", "Post"), nrow(.)/3))

p_peak_amp_frequency <- dat_cb_combine_event %>% 
  mutate(cluster = map_dbl(name, ~dat_cell_trace_cluster[.x])) %>% 
  filter(value > 3) %>% 
  ddply(.,.(name, Group, cluster), summarise, n = length(value)) %>% 
  mutate(Freq = n/4) %>% 
  # merge(dat_cell_trace_cluster_num,., by = c("name", "cluster", "Group"), all.x = TRUE) %>% 
  # select(name, Group, cluster, Freq) %>% 
  # mutate(Freq = ifelse(is.na(Freq), 0, Freq)) %>% 
  mutate(Group = factor(Group, levels = c("Pre", "Cond", "Post"))) %>% 
  mutate(cluster = factor(cluster, levels = c(2,1))) %>% 
  ggplot(., aes(Group, Freq, col = Group))+
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(aes(colour = Group, shape = Group),width = 0.2,  size=2, alpha= 0.3)+
  scale_color_manual(values=c("#8491B4FF", "#00A087FF", "#3C5488FF"))+
  facet_grid(cols = vars(cluster))+
  labs(x="", y="Freq. of Ca2+ events (Hz)")+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  theme(legend.title = element_blank(), legend.position = "none")

t_peak_amp_frequency <- dat_cb_combine_event %>% 
  mutate(cluster = map_dbl(name, ~dat_cell_trace_cluster[.x])) %>% 
  mutate(Group = factor(Group, levels = c("Pre", "Cond", "Post"))) %>% 
  filter(value > 3) %>% 
  ddply(.,.(name, Group, cluster), summarise, n = length(value)) %>% 
  mutate(Freq = n/4) %>% 
  # merge(dat_cell_trace_cluster_num,., by = c("name", "cluster", "Group"), all.x = TRUE) %>% 
  # select(name, Group, cluster, Freq) %>% 
  # mutate(Freq = ifelse(is.na(Freq), 0, Freq)) %>% 
  mutate(Group = factor(Group, levels = c("Pre", "Cond", "Post"))) %>% 
  filter(cluster==2) %>% 
  filter(Group != "Post") %>% 
  wilcox.test(Freq~Group,.)


setwd("~cchen/Documents/neuroscience/Pn\ project/Figure/PDF/")
cairo_pdf("p_peak_amp.pdf", width = 80/25.6, height = 65/25.6, family = "Arial" )
p_peak_amp
dev.off()


dat_peak_amp_sta <- dat_cb_combine_event %>% 
  mutate(cluster = map_dbl(name, ~dat_cell_trace_cluster[.x])) %>% 
  mutate(Group = factor(Group, levels = c("Pre", "Cond", "Post"))) %>% 
  filter(value > 3) %>%
  ddply(.,.(Group, cluster), summarise, n = length(value)) 


t_peak_amp_sta <- dat_cb_combine_event %>% 
  mutate(cluster = map_dbl(name, ~dat_cell_trace_cluster[.x])) %>% 
  mutate(Group = factor(Group, levels = c("Pre", "Cond", "Post"))) %>% 
  filter(value > 3) %>%
  ddply(.,.(Group, name, cluster), summarise, n_value = length(value)) %>% 
  filter(cluster==1) %>% 
  aov(n_value~ Group, .)

summary(t_peak_amp_sta)
TukeyHSD(t_peak_amp_sta)

t_peak_amp <- dat_cb_combine_event %>% 
  mutate(cluster = map_dbl(name, ~dat_cell_trace_cluster[.x])) %>% 
  mutate(Group = factor(Group, levels = c("Pre", "Cond", "Post"))) %>% 
  filter(cluster ==1) %>% 
  filter(value > 3) %>% 
  aov(value~Group,.)

summary(t_peak_amp)
TukeyHSD(t_peak_amp)
