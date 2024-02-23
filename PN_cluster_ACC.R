back_crossing_frame <- tibble(ID = c("m3", "m7", "m17", "m18", "m855", "m857"), Frame = c(1003, 1632, 2113,1923,1572,1892), 
                              Last_frame= c(1756, 2777, 2536, 2238, 3062, 2306)) ## m3, last one 3942

c_miniscope_matlab_ft <- function(file_trace) {
  ## import and format the data
  dat_trace1 <- raveio::read_mat(file_trace)
  ID_mouse <- str_extract(file_trace, regex("m\\d+"))
  
  t_crossing_back <- back_crossing_frame %>% 
    filter(ID== ID_mouse )
  t_crossing <- c(dat_trace1[[2]][5], t_crossing_back$Frame, t_crossing_back$Last_frame)
  
  cross_ID <- dat_trace1$global_map %>% 
    as_tibble() %>% 
    select(V1, V4, V5) %>% 
    mutate_all(na_if, 0) %>% 
    drop_na() %>% 
    pull(3)
  
  dat_stim_trace <- vector(mode = "list", length = length(t_crossing))
  
  
  for (i in seq_along(t_crossing)) {
    dat_trace <- dat_trace1[[8]] %>% 
      as_tibble() %>% 
      select(all_of(cross_ID)) %>% 
      apply(., 2, scale)
    
    ## number of rows to be binned
    n <- 2 # 0.05*2=0.1
    t1_p <- t_crossing[i]
    ## extract cell activity when they cross the border
    #dat_stim1 <- dat_trace[(t1_p-40):(t1_p+140-1),] 
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


## analyze for all mouse 
mouse_file <- as.list(list.files("~cchen/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/FD-processed", full.names = T))

dat_cell_trace <- mapply(c_miniscope_matlab_ft, mouse_file, SIMPLIFY = F)


dat_cell_pre_trace <- lapply(dat_cell_trace, function(x) x[[1]]) %>% 
  do.call(cbind,.)

dat_cell_cond_trace <- lapply(dat_cell_trace, function(x) x[[2]]) %>% 
  do.call(cbind,.)

dat_cb_combine <- lapply(dat_cell_trace, function(x) x[[3]]) %>% 
  do.call(cbind,.) %>% 
  rbind(dat_cell_pre_trace, dat_cell_cond_trace,.)


fviz_nbclust(t(dat_cb_combine), kmeans, method = "silhouette")

dat_cell_trace_cluster <- kmeans(t(dat_cb_combine), centers = 3, iter.max = 10)[[1]]


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

p_trace <- dat_cb_combine %>% 
  add_column(Group = rep(c("Crossing", "Crossing_back", "Last_crossing"), each = 40)) %>% 
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


## for Pre Cond and Post
c_miniscope_matlab_ft <- function(file_trace) {
  ## import and format the data
  dat_trace1 <- raveio::read_mat(file_trace)
  ID <- str_extract(file_trace, regex("m\\d+"))
  dat_trace1[[2]][1] <- ifelse(ID=="m857", 1694, dat_trace1[[2]][1])
  num_compare <- c(4,7, 8)
  t_crossing <- c(1,4, 5)
  
  cross_ID <- dat_trace1$global_map %>% 
    as_tibble() %>% 
    select(V1, V4, V5) %>% 
    mutate_all(na_if, 0) %>% 
    drop_na()
  
  dat_stim_trace <- vector(mode = "list", length = length(num_compare))
  
  
  for (i in seq_along(num_compare)) {
    global_cell <- pull(cross_ID[,i])
    dat_trace <- dat_trace1[[num_compare[i]]] %>% 
      as_tibble() %>% 
      select(all_of(global_cell)) %>% 
      apply(., 2, scale)
    
    ## number of rows to be binned
    n <- 2 # 0.05*2=0.1
    t1_p <- dat_trace1$crossing[t_crossing[i]]
    ## extract cell activity when they cross the border
    #dat_stim1 <- dat_trace[(t1_p-40):(t1_p+140-1),] 
    dat_stim1 <- dat_trace[(t1_p-40):(t1_p+40-1),] 
    
    dat_stim <- aggregate(dat_stim1,list(rep(1:(nrow(dat_stim1)%/%n+1),each=n,len=nrow(dat_stim1))),mean)[-1]
    # apply(., 2, scale) %>% 
    # as_tibble() %>% 
    # replace(is.na(.), 0)
    
    
    colnames(dat_stim) <- str_c(ID,"Cell", 1: ncol(dat_stim))
    dat_stim_trace[[i]] <- dat_stim
  }
  return(dat_stim_trace)
}


## analyze for all mouse 
mouse_file <- as.list(list.files("~cchen/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/FD-processed", full.names = T))

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
  filter(cluster==3) %>% 
  filter(Group!= "Post") %>% 
  wilcox.test(mean~Group, .,paired = T)

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