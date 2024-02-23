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
  group_day <- c("Pre", "Cond", "Post")
  
  
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
    dat_stim <- dat_trace[(t1_p-39):(t1_p+40),] 
    
    colnames(dat_stim) <- str_c(ID,"Cell", 1: ncol(dat_stim))
    
    peak_positions <- dat_stim %>% 
      apply(., 2, function(x) findpeaks(x, nups = 1, minpeakdistance = 2, minpeakheight = 3,threshold = mean(x)+ 2*sd(x)))## 0.01 for nonscaled data
    
    
    
    dat_trace_wave <- c()
    ## extract ca2+ events waves of each event
    for (j in 1:length(peak_positions)) {
      
      if (is.null(peak_positions[[j]]) ) {
        next  # Skip to the next iteration of the loop
      }
      peak_time <- peak_positions[[j]] 
      
      cell_name <- colnames(dat_stim)[j]
      trace_tem <- dat_stim[,j]
      
      dat_trace_cell <- vector(mode = 'list', length = nrow(peak_time))
      for (d in 1: nrow(peak_time)) {
        dat_event1 <- trace_tem[peak_time[d,3]: peak_time[d,4]] %>% 
          as_tibble() %>% 
          add_column(Time = (peak_time[d,3]: peak_time[d,4]) - peak_time[d,2] ) %>% 
          add_column(name = cell_name) 
        
        dat_trace_cell[[d]] <- dat_event1
        
        
      }
      
      dat_trace_wave <- dat_trace_cell %>% 
        do.call(rbind,.) %>% 
        rbind(dat_trace_wave,.)
      
    }
    
    dat_stim_trace[[i]] <- dat_trace_wave %>% 
      add_column(Group = group_day[i])
  }
  return(dat_stim_trace)
}



mouse_file <- as.list(list.files("~cchen/Documents/neuroscience/Pn\ project/Data_analysis/For\ revision/PC_data", full.names = T))

dat_cell_trace <- mapply(c_miniscope_matlab_ft, mouse_file, SIMPLIFY = F)

dat_cell_pre_trace <- lapply(dat_cell_trace, function(x) x[[1]]) %>% 
  do.call(rbind,.)


dat_cell_cond_trace <- lapply(dat_cell_trace, function(x) x[[2]]) %>% 
  do.call(rbind,.) 

p_cell_acti <- lapply(dat_cell_trace, function(x) x[[3]]) %>% 
  do.call(rbind,.) %>% 
  rbind(dat_cell_pre_trace, dat_cell_cond_trace,.) %>%
  mutate(cluster = map_dbl(name, ~dat_cell_trace_cluster[.x])) %>% 
  mutate(cluster = factor(cluster, levels = c(1,2))) %>% 
  ddply(.,.(Time, Group, cluster), summarise, mean_acti = mean(value), se=sd(value)/sqrt(length(value))) %>% 
  drop_na() %>% 
  mutate(Group = factor(Group, levels = c("Pre", "Cond", "Post"))) %>% 
  filter(cluster==2) %>% 
  ggplot(., aes(Time, mean_acti))+
  geom_line(aes(color = Group, group = Group))+
  geom_ribbon(aes(ymin = mean_acti - se, ymax = mean_acti + se, fill= Group), alpha = 0.5)+
  theme_void()+
  theme(legend.title = element_blank(), legend.position = c(0.8, 0.8))+
  annotate("segment", x = 5, xend = 9, y = 0.5, yend = 0.5)+
  annotate("segment", x = 9, xend = 9, y = 0.5, yend = 1)
  
setwd("~cchen/Documents/neuroscience/Pn\ project/Figure/PDF/")
cairo_pdf("p_cell_acti.pdf", width = 80/25.6, height = 65/25.6, family = "Arial" )
p_cell_acti
dev.off()

## calculate peak amplitude
p_peak <- lapply(dat_cell_trace, function(x) x[[3]]) %>% 
  do.call(rbind,.) %>% 
  rbind(dat_cell_pre_trace, dat_cell_cond_trace,.) %>%
  mutate(cluster = map_dbl(name, ~dat_cell_trace_cluster[.x])) %>% 
  mutate(cluster = factor(cluster, levels = c(1,2))) %>% 
  drop_na() %>% 
  mutate(Group = factor(Group, levels = c("Pre", "Cond", "Post"))) %>% 
  filter(cluster==1) %>% 
  filter(Time ==0) %>% 
  ggplot(.,aes(Group, value, col= Group))+
  geom_boxplot()+
  geom_jitter()


## compare crossing and back
back_crossing_frame <- tibble(ID_mouse = c("m617", "m625", "m676", "m684", "m687", "m685"), Frame_d3= c(436, 938,1670, 808, 1562, 774),
                              Frame_d6 = c(3570, 2116, 2794, 4100, 2912,2002),
                              Frame_d7 = c(1974, 2588, 1612, 1294, 1562, 2698))


c_miniscope_matlab_ft <- function(file_trace) {
  ## import and format the data
  dat_trace1 <- raveio::read_mat(file_trace)
  ID <- str_extract(file_trace, regex("m\\d+"))
  num_compare <- c(1,4, 5)
 
  
  cross_ID <- dat_trace1$global.map %>% 
    as_tibble() %>% 
    select(V1, V4, V5) %>% 
    mutate_all(na_if, 0) %>% 
    drop_na()
  
  t_crossing <- back_crossing_frame %>% 
    filter(ID_mouse == ID) %>% 
    select(-ID_mouse) 
  
  dat_stim_trace <- vector(mode = "list", length = length(num_compare))
  group_day <- c("Pre1", "Cond1", "Post1")
  
  
  for (i in seq_along(num_compare)) {
    global_cell <- pull(cross_ID[,i])
    dat_trace <- dat_trace1$traces[[num_compare[i]]] %>% 
      as_tibble() %>% 
      select(all_of(global_cell)) %>% 
      apply(., 2, scale)
    
    ## number of rows to be binned
    n <- 2 # 0.05*2=0.1
    #t1_p <- dat_trace1$crossing[t_crossing[i]]
    t1_p <- t_crossing[,i] %>% 
      pull()
    
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
    dat_stim <- dat_trace[(t1_p-39):(t1_p+40),] 
    
    colnames(dat_stim) <- str_c(ID,"Cell", 1: ncol(dat_stim))
    
    peak_positions <- dat_stim %>% 
      apply(., 2, function(x) findpeaks(x, nups = 1, minpeakdistance = 2, minpeakheight = 3.5 ,threshold = mean(x)+ 2*sd(x)))## 0.01 for nonscaled data
    
    
    
    dat_trace_wave <- c()
    ## extract ca2+ events waves of each event
    for (j in 1:length(peak_positions)) {
      
      if (is.null(peak_positions[[j]]) ) {
        next  # Skip to the next iteration of the loop
      }
      peak_time <- peak_positions[[j]] 
      
      cell_name <- colnames(dat_stim)[j]
      trace_tem <- dat_stim[,j]
      
      dat_trace_cell <- vector(mode = 'list', length = nrow(peak_time))
      for (d in 1: nrow(peak_time)) {
        dat_event1 <- trace_tem[peak_time[d,3]: peak_time[d,4]] %>% 
          as_tibble() %>% 
          add_column(Time = (peak_time[d,3]: peak_time[d,4]) - peak_time[d,2] ) %>% 
          add_column(ID = str_c(ID, cell_name)) 
        
        dat_trace_cell[[d]] <- dat_event1
        
        
      }
      
      dat_trace_wave <- dat_trace_cell %>% 
        do.call(rbind,.) %>% 
        rbind(dat_trace_wave,.)
      
    }
    
    dat_stim_trace[[i]] <- dat_trace_wave %>% 
      add_column(Group = group_day[i])
  }
  return(dat_stim_trace)
}

mouse_file <- as.list(list.files("~cchen/Documents/neuroscience/Pn\ project/Data_analysis/For\ revision/PC_data", full.names = T))

dat_cell_trace1 <- mapply(c_miniscope_matlab_ft, mouse_file, SIMPLIFY = F)

dat_cell_pre_trace1 <- lapply(dat_cell_trace1, function(x) x[[1]]) %>% 
  do.call(rbind,.)


dat_cell_cond_trace1 <- lapply(dat_cell_trace1, function(x) x[[2]]) %>% 
  do.call(rbind,.) 

p_cell_acti <- lapply(dat_cell_trace1, function(x) x[[3]]) %>% 
  do.call(rbind,.) %>% 
  rbind(dat_cell_pre_trace1, dat_cell_cond_trace1,.) %>%
  ddply(.,.(Time, Group), summarise, mean_acti = mean(value), length = length(value),se=sd(value)/sqrt(length(value))) %>% 
  drop_na() %>% 
  filter(Time >= -5 & Time <= 10) %>% 
  mutate(Group = factor(Group, levels = c("Pre1", "Cond1", "Post1"))) %>% 
  ggplot(., aes(Time, mean_acti))+
  geom_line(aes(color = Group, group = Group))+
  geom_ribbon(aes(ymin = mean_acti - se, ymax = mean_acti + se, fill= Group), alpha = 0.5)+
  theme_void()+
  theme(legend.title = element_blank(), legend.position = c(0.8, 0.8))+
  annotate("segment", x = 5, xend = 9, y = 0.5, yend = 0.5)+
  annotate("segment", x = 9, xend = 9, y = 0.5, yend = 1)



dat_cell_pre_trace_combine <- lapply(dat_cell_trace1, function(x) x[[1]]) %>% 
  do.call(rbind,.) %>% 
  rbind(dat_cell_pre_trace,.) %>% 
  ddply(.,.(Time, Group), summarise, mean_acti = mean(value), se=sd(value)/sqrt(length(value))) %>% 
  drop_na() %>% 
  mutate(Group = factor(Group, levels = c("Pre", "Pre1"))) %>% 
  ggplot(., aes(Time, mean_acti))+
  geom_line(aes(color = Group, group = Group))+
  geom_ribbon(aes(ymin = mean_acti - se, ymax = mean_acti + se, fill= Group), alpha = 0.5)


dat_cell_cond_trace_combine <- lapply(dat_cell_trace1, function(x) x[[2]]) %>% 
  do.call(rbind,.) %>% 
  rbind(dat_cell_cond_trace,.) %>% 
  ddply(.,.(Time, Group), summarise, mean_acti = mean(value), se=sd(value)/sqrt(length(value))) %>% 
  drop_na() %>% 
  mutate(Group = factor(Group, levels = c("Cond", "Cond1"))) %>% 
  ggplot(., aes(Time, mean_acti))+
  geom_line(aes(color = Group, group = Group))+
  geom_ribbon(aes(ymin = mean_acti - se, ymax = mean_acti + se, fill= Group), alpha = 0.5)

dat_cell_post_trace <- lapply(dat_cell_trace, function(x) x[[3]]) %>% 
  do.call(rbind,.) 

dat_cell_cond_trace_combine <- lapply(dat_cell_trace1, function(x) x[[3]]) %>% 
  do.call(rbind,.) %>% 
  rbind(dat_cell_post_trace,.) %>% 
  ddply(.,.(Time, Group), summarise, mean_acti = mean(value), se=sd(value)/sqrt(length(value))) %>% 
  drop_na() %>% 
  mutate(Group = factor(Group, levels = c("Post", "Post1"))) %>% 
  ggplot(., aes(Time, mean_acti))+
  geom_line(aes(color = Group, group = Group))+
  geom_ribbon(aes(ymin = mean_acti - se, ymax = mean_acti + se, fill= Group), alpha = 0.5)

## compare the difference 
dat_cell_pre_trace_diff <- lapply(dat_cell_trace1, function(x) x[[1]]) %>% 
  do.call(rbind,.) %>% 
  rbind(dat_cell_pre_trace,.) %>% 
  ddply(.,.(Time, Group), summarise, mean_acti = mean(value)) %>% 
  pivot_wider(names_from = Group, values_from = mean_acti) %>% 
  drop_na() %>% 
  mutate(diff = Pre1 - Pre) %>% 
  mutate(Group = "Pre") %>% 
  select(Time, diff, Group)

dat_cell_cond_trace_diff <- lapply(dat_cell_trace1, function(x) x[[2]]) %>% 
  do.call(rbind,.) %>% 
  rbind(dat_cell_cond_trace,.) %>% 
  ddply(.,.(Time, Group), summarise, mean_acti = mean(value)) %>% 
  pivot_wider(names_from = Group, values_from = mean_acti) %>% 
  drop_na() %>% 
  mutate(diff = Cond1 - Cond) %>% 
  mutate(Group = "Cond") %>% 
  select(Time, diff, Group)

p_diff <- lapply(dat_cell_trace1, function(x) x[[3]]) %>% 
  do.call(rbind,.) %>% 
  rbind(dat_cell_post_trace,.) %>% 
  ddply(.,.(Time, Group), summarise, mean_acti = mean(value)) %>% 
  pivot_wider(names_from = Group, values_from = mean_acti) %>% 
  drop_na() %>% 
  mutate(diff = Post1 - Post) %>% 
  mutate(Group = "Post") %>% 
  select(Time, diff, Group) %>% 
  rbind(dat_cell_pre_trace_diff, dat_cell_cond_trace_diff,.) %>% 
  ggplot(., aes(Time, diff, color = Group))+
  geom_line()

