c_miniscope_matlab_ft <- function(file_trace) {
  ## import and format the data
  dat_trace1 <- raveio::read_mat(file_trace)[[1]]
  names(dat_trace1) <- c("traces", "Spatial", "global.map", "crossing", "choices")
  ID <- str_extract(file_trace, regex("m\\d+"))
  num_compare <- c(1,4, 5)
  t_crossing <- c(1,4, 5)
  
  cell_bad <- which(dat_trace1$choices==0)
  
  cross_ID <- dat_trace1$global.map %>% 
    as_tibble() %>% 
    select(V1, V4, V5) %>% 
    mutate_all(na_if, 0) %>% 
    drop_na() %>% 
    #filter(!V5 %in% cell_bad)
  
  dat_stim_trace <- vector(mode = "list", length = length(num_compare))
  
  
  for (i in seq_along(num_compare)) {
    global_cell <- pull(cross_ID[,i])
    dat_trace <- dat_trace1$traces[[num_compare[i]]][[1]] %>% 
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



mouse_file <- as.list(list.files("~cchen/Documents/neuroscience/Pn\ project/Data_analysis/For\ revision/PC_data3", full.names = T))

dat_cell_trace <- mapply(c_miniscope_matlab_ft, mouse_file, SIMPLIFY = F)



## for the first 2s vs last 2s
dat_trace_combine <- vector(mode = 'list', length = 3)

for (i in seq_along(group_day)){
  dat_trace_combine[[i]] <- lapply(dat_cell_trace, function(x)  x[[i]]) %>% 
    do.call(cbind,.) %>% 
    mutate(Time = comp_time) %>% 
    mutate(Group = group_day[i]) %>% 
    pivot_longer(-c(Time, Group)) 
  
}


p_trace_2s <- dat_trace_combine %>% 
  do.call(rbind,.) %>% 
  mutate(Time_group = ifelse(Time >0, 'T2', 'T1')) %>% 
  mutate(ID = str_extract(name, "m\\d+")) %>% 
  mutate(Group = factor(Group, levels = c("Pre", "Cond", "Post"))) %>% 
  ddply(., .(Time_group,name, Group, ID), summarise, mean_acti = mean(value)) %>% 
  ddply(.,.(Time_group,Group, ID), summarise, mean_acti1 = mean(mean_acti), n= length(mean_acti)) %>% 
  pivot_wider(names_from = Time_group, values_from = mean_acti1) %>% 
  mutate(diff = T2-T1) %>% 
  ggplot(., aes(Group, diff, group = Group))+
  geom_boxplot(outlier.shape = NA)+
  geom_line(aes(group=ID), colour="gray90")+
  geom_jitter(aes(colour = Group, shape = Group),width = 0.2,  size=2, alpha= 0.3)+
  scale_color_manual(values=c("#8491B4FF", "#00A087FF", "#3C5488FF"))+
  labs(x="", y="Population activity (∆F/F)")+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  scale_y_continuous(limits = c(-0.2, 0.6))+
  theme(legend.title = element_blank(), legend.position = "none")


p_trace_s2_cell <- dat_trace_combine %>% 
  do.call(rbind,.) %>% 
  mutate(Group = factor(Group, levels = c("Pre", "Cond", "Post"))) %>% 
  filter(Time <= 0) %>% 
  ddply(., .(name, Group), summarise, mean_acti = mean(value)) %>% 
  ggplot(., aes(Group, mean_acti, group = Group))+
  geom_violin()+
  geom_jitter(aes(colour = Group, shape = Group),width = 0.2,  size=2, alpha= 0.3)+
  scale_color_manual(values=c("#8491B4FF", "#00A087FF", "#3C5488FF"))+
  labs(x="", y="Population activity (∆F/F)")+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  scale_y_continuous(limits = c(-1.5, 8))+
  theme(legend.title = element_blank(), legend.position = "none")

t_trace_s2_cell <- dat_trace_combine %>% 
  do.call(rbind,.) %>% 
  mutate(Group = factor(Group, levels = c("Pre", "Cond", "Post"))) %>% 
  filter(Time <= 0) %>% 
  ddply(., .(name, Group), summarise, mean_acti = mean(value)) %>% 
  filter(Group != "Post") %>% 
  wilcox.test(mean_acti~Group, ., paired = T)


## clustering based on their activity on D7
dat_cell_pre_trace <- lapply(dat_cell_trace, function(x) x[[1]]) %>% 
  do.call(cbind,.)

dat_cell_cond_trace <- lapply(dat_cell_trace, function(x) x[[2]]) %>% 
  do.call(cbind,.)

dat_cb_combine <- lapply(dat_cell_trace, function(x) x[[3]]) %>% 
  do.call(cbind,.) %>% 
  rbind(dat_cell_pre_trace, dat_cell_cond_trace,.)

cairo_pdf("p_cluster_opto.pdf", width = 80/25.6, height = 60/25.6, family = "Arial")
fviz_nbclust(t(dat_cb_combine), kmeans, method = "silhouette")
dev.off() 


fviz_nbclust(t(dat_cb_combine), kmeans, method = "silhouette")

dat_cell_trace_cluster <- kmeans(t(dat_cb_combine), centers = 2, iter.max = 10)[[1]]


## clustering based on the activity of neurons during d7-----
## for comparing crossing, back-crossing and last crossing------
back_crossing_frame <- tibble(ID = c("m617", "m625", "m676", "m684", "m687", "m685"), Frame = c(1974, 2588, 1612, 1294, 2210, 2698), 
                              Last_frame= c(3016, 3168, 3106, 3180, 3726, 3640))

c_miniscope_matlab_ft <- function(file_trace) {
  dat_trace1 <- raveio::read_mat(file_trace)[[1]]
  
  ## import and format the data
  names(dat_trace1) <- c("traces", "Spatial", "global.map", "crossing")
  ID_mouse <- str_extract(file_trace, regex("m\\d+"))
  num_compare <- c(1,4, 5)
  t_crossing <- c(1,4, 5)
  
  
  t_crossing_back <- back_crossing_frame %>% 
    dplyr::filter(ID == ID_mouse)
  t_crossing <- c(dat_trace1$crossing[5], t_crossing_back$Frame, t_crossing_back$Last_frame)
  
  cross_ID <- dat_trace1$global.map%>% 
    as_tibble() %>% 
    select(V1, V4, V5) %>% 
    mutate_all(na_if, 0) %>% 
    drop_na() %>% 
    select(V1, V5) %>% 
    pull(V5)

    
  dat_stim_trace <- vector(mode = "list", length = length(t_crossing))
  
  
  for (i in seq_along(t_crossing)) {
    dat_trace <- dat_trace1$traces[[5]][[1]] %>% 
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



mouse_file <- as.list(list.files("~cchen/Documents/neuroscience/Pn\ project/Data_analysis/For\ revision/PC_data3", full.names = T))

dat_cell_trace <- mapply(c_miniscope_matlab_ft, mouse_file, SIMPLIFY = F)


dat_cell_pre_trace <- lapply(dat_cell_trace, function(x) x[[1]]) %>% 
  do.call(cbind,.)

dat_cell_cond_trace <- lapply(dat_cell_trace, function(x) x[[2]]) %>% 
  do.call(cbind,.)

dat_cb_combine <- lapply(dat_cell_trace, function(x) x[[3]]) %>% 
  do.call(cbind,.) %>% 
  rbind(dat_cell_pre_trace, dat_cell_cond_trace,.)


fviz_nbclust(t(dat_cb_combine), kmeans, method = "silhouette")

dat_cell_trace_cluster <- kmeans(t(dat_cb_combine), centers = 2, iter.max = 10)[[1]]

## get the number of each neurons in the cluster 
dat_cell_trace_cluster_num <- tibble(name = names(dat_cell_trace_cluster), cluster = dat_cell_trace_cluster) %>% 
  as_tibble() %>% 
  mutate(cell_ID =as.numeric(str_extract(name, "(?<=Cell)\\d+"))) %>% 
  mutate(mouse_ID = stringr::str_extract(name, ".*(?=Cell)")) %>% 
  select(-name)

## get the global map information

dat_global_cell <- vector(mode="list", length = length(mouse_file))

for (i in 1: length(mouse_file)) {
  file_trace <- mouse_file[[i]]
  dat_trace1 <- raveio::read_mat(file_trace)[[1]]
  names(dat_trace1) <- c("traces", "Spatial", "global.map", "crossing")
  ID <- str_extract(file_trace, regex("m\\d+"))
  
  cross_ID <- dat_trace1$global.map %>% 
    as_tibble() %>% 
    select(V1, V4, V5) %>% 
    rename(Pre = V1, Cond= V4, Post = V5) %>% 
    mutate_all(na_if, 0) %>% 
    drop_na() %>% 
    mutate(mouse_ID = ID) %>% 
    mutate(cell_ID = 1:nrow(.))
  
  dat_global_cell[[i]] <- cross_ID
  
}

dat_global_cell_combine <- dat_global_cell %>% 
  do.call(rbind,.) %>% 
  full_join(., dat_cell_trace_cluster_num, by = join_by(mouse_ID, cell_ID)) 

write.csv(dat_global_cell_combine, "global_cell_cluster2.csv")



####
## d^2, cross-day aligned neurons-----
back_crossing_frame <- tibble(ID_mouse = c("m617", "m625", "m676", "m684", "m687", "m685"), Frame_d3= c(436, 938,1670, 808, 1562, 774),
                              Frame_d6 = c(3570, 2116, 2794, 4100, 2912,2002),
                              Frame_d7 = c(1974, 2588, 1612, 1294, 1562, 2698))
c_miniscope_matlab_ft <- function(file_trace) {
  ## import and format the data
  dat_trace1 <- raveio::read_mat(file_trace)
  ID <- str_extract(file_trace, regex("m\\d+"))
  num_compare <- c(1, 4, 5)
  t_crossing <- c(1, 4, 5)
  
  
  t_crossing_back <- back_crossing_frame %>% 
    dplyr::filter(ID_mouse == ID) %>% 
    dplyr::select(Frame_d3, Frame_d6,Frame_d7) %>% 
    unlist() %>% 
    unname()
  
  
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
    t1_p_back <- t_crossing_back[i]
    
    ## for corssing less than 2s
    if (t1_p < 40){
      column_means = colMeans(dat_trace[1:t1_p, ])
      
      # Replicate these mean values 10 times
      n_reapt <- (40 - t1_p)
      new_rows = matrix(rep(column_means, each = n_reapt), nrow = n_reapt, ncol = ncol(dat_trace))
      
      # Combine the new matrix with the original one
      dat_trace_add = rbind(new_rows, dat_trace)
      
    } else {
      dat_trace_add <- dat_trace
    }
    ## extract cell activity when they cross the border
    t1_p <- ifelse(t1_p < 40, 40, t1_p)
    dat_stim1 <- dat_trace_add[(t1_p-40):(t1_p+40-1),] 
    
    dat_stim_cross <- aggregate(dat_stim1,list(rep(1:(nrow(dat_stim1)%/%n+1),each=n,len=nrow(dat_stim1))),mean)[-1]
    colnames(dat_stim_cross) <- str_c(ID,"Cell", 1: ncol(dat_stim_cross))
    
    ## extract cell activity when cross back
    dat_stim2 <- dat_trace[(t1_p_back-40):(t1_p_back+40-1),] 
    
    dat_stim_back <- aggregate(dat_stim2,list(rep(1:(nrow(dat_stim2)%/%n+1),each=n,len=nrow(dat_stim2))),mean)[-1]
    colnames(dat_stim_back) <- str_c(ID,"Cell", 1: ncol(dat_stim_back))
    
    
    ## calculate the d' of each cells
    comp_time <- seq(-2, 1.9, by=0.1)
    dat_stim_cross_sta <- dat_stim_cross %>% 
      as_tibble() %>% 
      add_column(Time = comp_time) %>% 
      pivot_longer(-Time) %>% 
      ddply(.,.(name), summarise,mean_cross=max(value)[1],sd_cross=sd(value))
    
    dat_stim_back_sta <- dat_stim_back %>% 
      as_tibble() %>% 
      add_column(Time = comp_time) %>% 
      pivot_longer(-Time) %>% 
      ddply(.,.(name), summarise,mean_back=max(value)[1],sd_back=sd(value))
    
    dat_stim_combine <- full_join(dat_stim_cross_sta, dat_stim_back_sta, by = 'name') %>% 
      mutate(sd_pool = sqrt(sd_cross^2 + sd_back^2)/sqrt(2)) %>% 
      mutate(sd_pool = ifelse(sd_pool < 1e-06, 1e-06, sd_pool)) %>% 
      mutate(d = (mean_cross - mean_back)/sd_pool) %>% 
      mutate(d2 = d^2) %>% 
      dplyr::select(name, d2) %>% 
      mutate(ID = ID)
    
    dat_stim_trace[[i]] <- dat_stim_combine
    
    
  }
  return(dat_stim_trace)
}

mouse_file <- as.list(list.files("~cchen/Documents/neuroscience/Pn\ project/Data_analysis/For\ revision/PC_data", full.names = T))

dat_cell_trace <- mapply(c_miniscope_matlab_ft, mouse_file, SIMPLIFY = F)


dat_cell_d_pre <- lapply(dat_cell_trace, function(x)  x[[1]]) %>% 
  do.call(rbind,.) %>% 
  mutate(Group = "Pre")

dat_cell_d_cond <- lapply(dat_cell_trace, function(x)  x[[2]]) %>% 
  do.call(rbind,.) %>% 
  mutate(Group = "Cond")

dat_cell_d_combine <- lapply(dat_cell_trace, function(x)  x[[3]]) %>% 
  do.call(rbind,.) %>% 
  mutate(Group = "Post") %>% 
  rbind(dat_cell_d_pre, dat_cell_d_cond,.) %>% 
  ddply(., .(ID, Group), summarise,n=length(d2),mean=mean(d2),sd=sd(d2),se=sd(d2)/sqrt(length(d2))) %>% 
  mutate(Group = factor(Group, levels = c("Pre", "Cond","Post")))


dat_cell_d_combine_sta <- lapply(dat_cell_trace, function(x)  x[[3]]) %>% 
  do.call(rbind,.) %>% 
  mutate(Group = "Post") %>% 
  rbind(dat_cell_d_pre,dat_cell_d_cond, .) %>% 
  ddply(., .(ID, Group), summarise,n=length(d2),value=mean(d2),sd=sd(d2),se=sd(d2)/sqrt(length(d2))) %>% 
  mutate(Group = factor(Group, levels = c("Pre", "Cond","Post"))) %>% 
  ddply(., .(Group),summarise,n=length(value),mean=mean(value),sd=sd(value),se=sd(value)/sqrt(length(value)))


p_d_combine <- ggplot(dat_cell_d_combine, aes(Group, mean, group = Group))+
  geom_boxplot(outlier.shape = NA)+
  geom_line(aes(group=ID), colour="gray90")+
  geom_jitter(aes(colour = Group, shape = Group),width = 0.2,  size=2)+
  scale_color_manual(values=c("#8491B4FF", "#00A087FF", "#3C5488FF"))+
  labs(x="", y="(d')^2")+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  scale_y_continuous(limits = c(-1, 20))+
  theme(legend.title = element_blank(), legend.position = "none")


t_d_combine <- dat_cell_d_combine %>% 
  filter(Group != "Post") %>% 
  wilcox.test(mean~Group, ., paired = T)

## compare by cells
p_d_combine_cell <- lapply(dat_cell_trace, function(x)  x[[3]]) %>% 
  do.call(rbind,.) %>% 
  mutate(Group = "Post") %>% 
  rbind(dat_cell_d_pre, dat_cell_d_cond, .) %>% 
  mutate(Group = factor(Group, levels = c("Pre", "Cond","Post"))) %>% 
  ggplot(., aes(Group, d2, group = Group))+
  geom_violin( )+
  geom_jitter(aes(colour = Group, shape = Group),width = 0.2,  size=2, alpha = 0.3 )+
  scale_color_manual(values=c("#8491B4FF", "#00A087FF", "#3C5488FF"))+
  labs(x="", y="(d')^2" )+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  scale_y_continuous(limits = c(-1, 500))+
  theme(legend.title = element_blank(), legend.position = "none")

t_d_combine_cell <- lapply(dat_cell_trace, function(x)  x[[3]]) %>% 
  do.call(rbind,.) %>% 
  mutate(Group = "Post") %>% 
  rbind(dat_cell_d_pre,dat_cell_d_cond, .) %>% 
  mutate(Group = factor(Group, levels = c("Pre","Cond" ,"Post"))) %>% 
  filter(Group !="Cond") %>% 
  wilcox.test(d2~ Group,., paired = T)



## check microzone based on binary calcium spikes
c_miniscope_matlab_ft_time <- function(file_trace) {
  ## import and format the data
  dat_trace1 <- raveio::read_mat(file_trace)[[1]]
  names(dat_trace1) <- c("traces", "Spatial", "global.map", "crossing", "choices")
  ID <- str_extract(file_trace, regex("m\\d+"))
  num_compare <- c(1,4, 5)
  t_crossing <- c(1,4, 5)
  
  cell_bad <- which(dat_trace1$choices==0)
  
  cross_ID <- dat_trace1$global.map %>% 
    as_tibble() %>% 
    select(V1, V4, V5) %>% 
    mutate_all(na_if, 0) %>% 
    drop_na() %>% 
    filter(!V5 %in% cell_bad)
  
  dat_spike_time <- list.files("~cchen/Documents/neuroscience/Pn\ project/Data_analysis/For\ revision/PC_spike_time5/", pattern = ID, full.names = T)
  
  
  dat_stim_trace <- vector(mode = "list", length = length(num_compare))
  
  
  for (i in seq_along(num_compare)) {
    global_cell <- pull(cross_ID[,i])
    
    
    dat_event <- read.xlsx(dat_spike_time[i], colNames = F) %>% 
      as_tibble() %>% 
      select(all_of(global_cell))
    
    ## number of rows to be binned
    n <- 2 # 0.05*2=0.1
    t1_p <- dat_trace1$crossing[t_crossing[i]]
    
    if (t1_p < 40){
      
      ## for event dat
      n_reapt <- (40 - t1_p)
      means <- dat_event %>% 
        slice(1:t1_p) %>% 
        summarize(across(everything(), mean, na.rm = TRUE))
      
      # Create a new data frame with 10 rows of mean values
      new_df <- map_df(means, ~ rep(., n_reapt))
      
      # Combine the new data frame with the original data frame
      dat_event <- bind_rows(new_df, dat_event)
      
      
    } else {
      
      dat_event <- dat_event
    }
    
    
    
    t1_p <- ifelse(t1_p < 40, 40, t1_p)
    
    dat_stim <- dat_event[(t1_p-40):(t1_p+40-1),] 
    

    
    
    colnames(dat_stim) <- str_c(ID,"Cell", 1: ncol(dat_stim))
    dat_stim_trace[[i]] <- dat_stim
    
    
  }
  return(dat_stim_trace)
}

mouse_file <- as.list(list.files("~cchen/Documents/neuroscience/Pn\ project/Data_analysis/For\ revision/PC_data5", full.names = T))


dat_cell_trace_event <- mapply(c_miniscope_matlab_ft_time, mouse_file, SIMPLIFY = F)

dat_spike <- dat_cell_trace_event[[6]] %>% 
  do.call(rbind,.) %>% 
  as.matrix() %>% 
  cor() %>% 
  round(.,2)

corrplot(dat_spike, method = "square",  order = "hclust", tl.col = "black", tl.srt = 45, tl.pos = "n")



## for comparing crossing, back-crossing and last crossing------
back_crossing_frame <- tibble(ID_mouse = c("m617", "m625", "m676", "m684", "m687", "m685"), Frame = c(1974, 2588, 1612, 1294, 2210, 2698), 
                              Last_frame= c(3016, 3168, 3106, 3180, 3726, 3640))

c_miniscope_matlab_ft <- function(file_trace) {
  ## import and format the data
  
  dat_trace1 <- raveio::read_mat(file_trace)[[1]]
  names(dat_trace1) <- c("traces", "Spatial", "global.map", "crossing", "choices")
  ID <- str_extract(file_trace, regex("m\\d+"))
  
  cell_bad <- which(dat_trace1$choices==0)

  cross_ID <- dat_trace1$global.map %>%
    as_tibble() %>%
    select(V1, V4, V5) %>%
    mutate_all(na_if, 0) %>%
    drop_na() %>%
    filter(!V5 %in% cell_bad) %>%
    pull(V5)

# cross_ID <- dat_trace1$global.map %>%
#   as_tibble() %>%
#   select(V5) %>%
#   filter(!V5 ==0) %>%
#   filter(!V5 %in% cell_bad) %>%
#   pull(V5)

  dat_spike_time <- list.files("~cchen/Documents/neuroscience/Pn\ project/Data_analysis/For\ revision/PC_spike_time5/", pattern = ID, full.names = T)
  
  
  t_crossing_back <- back_crossing_frame %>% 
    dplyr::filter(ID_mouse == ID)
  t_crossing <- c(dat_trace1$crossing[5], t_crossing_back$Frame, t_crossing_back$Last_frame)
  
  
  
  dat_stim_trace <- vector(mode = "list", length = length(t_crossing))
  
  
  for (i in seq_along(t_crossing)) {
    dat_event <- read.xlsx(dat_spike_time[3], colNames = F) %>% 
      as_tibble() %>% 
      select(all_of(cross_ID)) %>% 
      as.matrix()
    
    ## number of rows to be binned
    n <- 2 # 0.05*2=0.1
    t1_p <- t_crossing[i]
    
    if (t1_p < 40){
      column_means = colMeans(dat_event[1:t1_p, ])
      
      # Replicate these mean values 10 times
      n_reapt <- (40 - t1_p)
      new_rows = matrix(rep(column_means, each = n_reapt), nrow = n_reapt, ncol = ncol(dat_event))
      
      # Combine the new matrix with the original one
      dat_event = rbind(new_rows, dat_event)
      
    } else {
      dat_event <- dat_event
    }
    
    
    
    ## extract cell activity when they cross the border
    #dat_stim1 <- dat_trace[(t1_p-40):(t1_p+140-1),] 
    t1_p <- ifelse(t1_p < 40, 40, t1_p)
    dat_stim1 <- dat_event[(t1_p-40):(t1_p+40-1),] 
    dat_stim <- aggregate(dat_stim1,list(rep(1:(nrow(dat_stim1)%/%n+1),each=n,len=nrow(dat_stim1))),mean)[-1]
    
    
    
    colnames(dat_stim) <- str_c(ID,"Cell", 1: ncol(dat_stim))
    dat_stim_trace[[i]] <- dat_stim
  }
  return(dat_stim_trace)
}


## analyze for all mouse 
mouse_file <- as.list(list.files("~cchen/Documents/neuroscience/Pn\ project/Data_analysis/For\ revision/PC_data5", full.names = T))

dat_cell_trace <- mapply(c_miniscope_matlab_ft, mouse_file, SIMPLIFY = F)

dat_spike <- dat_cell_trace[[5]][[1]] %>% 
  #do.call(rbind,.) %>% 
  as_tibble() %>% 
  select_if(~any(. != 0)) %>% 
  as.matrix() %>% 
  cor() %>% 
  round(.,2) 


corrplot(dat_spike, method = "square",  order = "hclust",tl.col = "black", tl.srt = 45, tl.pos = "n",is.corr = FALSE)


## for ca2_ spike-----
mouse_file <- as.list(list.files("~cchen/Documents/neuroscience/Pn\ project/Data_analysis/For\ revision/PC_data5", full.names = T))

file_trace <- mouse_file[[2]]
dat_trace1 <- raveio::read_mat(file_trace)[[1]]
names(dat_trace1) <- c("traces", "Spatial", "global.map", "crossing", "choices")
ID <- str_extract(file_trace, regex("m\\d+"))

cell_bad <- which(dat_trace1$choices==0)

cross_ID <- dat_trace1$global.map %>%
  as_tibble() %>%
  select(V1, V4, V5) %>%
  mutate_all(na_if, 0) %>%
  drop_na() %>%
  filter(!V5 %in% cell_bad) %>%
  pull(V1)

# cross_ID <- dat_trace1$global.map %>%
#   as_tibble() %>%
#   select(V1) %>%
#   filter(!V1 ==0) %>%
#   pull(V1)

dat_spike_time <- list.files("~cchen/Documents/neuroscience/Pn\ project/Data_analysis/For\ revision/PC_spike_time5/", pattern = ID, full.names = T)

dat_event <- read.xlsx(dat_spike_time[1], colNames = F) %>% 
  as_tibble() %>% 
  select(all_of(cross_ID)) %>% 
  slice((695*2):(850*2)) %>% 
  select_if(~any(. != 0)) %>% 
  as.matrix()


dat_trace_d1_res <- dat_event %>%
  cor() %>% 
  round(.,2)

#library(corrplot)
corrplot(dat_trace_d1_res, method = "square",  order = "hclust",tl.col = "black", tl.srt = 45, tl.pos = "n",is.corr = FALSE)

