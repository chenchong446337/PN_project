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

comp_time <- seq(-2, 1.9, by=0.1)
group_day <- c("Pre", "Cond", "Post")
score_range <- range(dat_cell_trace)

for(i in c(1:3)){
  cell_order <- lapply(dat_cell_trace, function(x)  x[[i]]) %>% 
    do.call(cbind,.) %>% 
    mutate(Time = comp_time) %>% 
    pivot_longer(-Time) %>% 
    ddply(.,.(name), summarise, mean = mean(value)) %>% 
    arrange(mean)
  
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
    theme(legend.position = "none")
  assign(str_c("p_heat_", group_day[i]), p_heat)
  
}

p_heat <- plot_grid(p_heat_Pre, p_heat_Cond, p_heat_Post, nrow  = 1)

setwd("~cchen/Documents/neuroscience/Pn\ project/Figure/PDF/")
cairo_pdf("p_heat_pc.pdf", width = 140/25.6, height = 60/25.6, family = "Arial")
p_heat
dev.off()

## plot to show the amplitude of events during this period
dat_trace_combine <- vector(mode = 'list', length = 3)

for (i in seq_along(group_day)){
  dat_trace_combine[[i]] <- lapply(dat_cell_trace, function(x)  x[[i]]) %>% 
    do.call(cbind,.) %>% 
    mutate(Time = comp_time) %>% 
    mutate(Group = group_day[i]) %>% 
    pivot_longer(-c(Time, Group)) 
  
}

p_event_pre <- dat_trace_combine %>% 
  do.call(rbind,.) %>% 
  mutate(Group = factor(Group, levels = c("Pre", "Cond", "Post"))) %>% 
  ggplot(., aes(Time, value, group = name))+ 
  geom_line()+
  facet_grid(rows = vars(Group))+
  geom_hline(yintercept = 2.5, linetype=2, color='red')
  

## for individual mouse
dat_cell_pre_trace <- lapply(dat_cell_trace, function(x) x[[1]]) %>% 
  do.call(cbind,.) %>% 
  add_column(Group = "Pre") %>% 
  pivot_longer(-Group) %>% 
  mutate(ID = str_extract(name, "m\\d+")) 

dat_cell_cond_trace <- lapply(dat_cell_trace, function(x) x[[2]]) %>% 
  do.call(cbind,.) %>% 
  add_column(Group = "Cond") %>% 
  pivot_longer(-Group) %>% 
  mutate(ID = str_extract(name, "m\\d+")) 


p_cell_acti <- lapply(dat_cell_trace, function(x) x[[3]]) %>% 
  do.call(cbind,.) %>% 
  add_column(Group = "Post") %>% 
  pivot_longer(-Group) %>% 
  mutate(ID = str_extract(name, "m\\d+")) %>% 
  rbind(dat_cell_pre_trace, dat_cell_cond_trace,.) %>%
  mutate(Group = factor(Group, levels = c("Pre", "Cond", "Post"))) %>% 
  ddply(., .(name, Group, ID), summarise, mean_acti = mean(value)) %>% 
  ddply(.,.(Group, ID), summarise, mean_acti1 = mean(mean_acti), n= length(mean_acti)) %>% 
  ggplot(., aes(Group, mean_acti1, group = Group))+
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

t_cell_acti <- lapply(dat_cell_trace, function(x) x[[3]]) %>% 
  do.call(cbind,.) %>% 
  add_column(Group = "Post") %>% 
  pivot_longer(-Group) %>% 
  mutate(ID = str_extract(name, "m\\d+")) %>% 
  rbind(dat_cell_pre_trace, dat_cell_cond_trace,.) %>%
  mutate(Group = factor(Group, levels = c("Pre", "Cond", "Post"))) %>% 
  ddply(., .(name, Group, ID), summarise, mean_acti = mean(value)) %>% 
  ddply(.,.(Group, ID), summarise, mean_acti1 = mean(mean_acti)) %>% 
  aov(mean_acti1~ Group,.)

summary(t_cell_acti)
TukeyHSD(t_cell_acti,  ordered = TRUE)  
pairwise.wilcox.test(t_cell_acti$mean_acti1, t_cell_acti$Group, p.adjust.method = "bonferroni", paired = TRUE)

# for individual neurons
dat_cell_pre_trace <- lapply(dat_cell_trace, function(x) x[[1]]) %>% 
  do.call(cbind,.) %>% 
  add_column(Group = "Pre") %>% 
  pivot_longer(-Group)

dat_cell_cond_trace <- lapply(dat_cell_trace, function(x) x[[2]]) %>% 
  do.call(cbind,.) %>% 
  add_column(Group = "Cond") %>% 
  pivot_longer(-Group)

p_cell_acti <- lapply(dat_cell_trace, function(x) x[[3]]) %>% 
  do.call(cbind,.) %>% 
  add_column(Group = "Post") %>% 
  pivot_longer(-Group) %>% 
  rbind(dat_cell_pre_trace, dat_cell_cond_trace,.) %>%
  mutate(Group = factor(Group, levels = c("Pre", "Cond", "Post"))) %>% 
  #pivot_longer(-Group) %>% 
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


t_cell_acti <- lapply(dat_cell_trace, function(x) x[[3]]) %>% 
  do.call(cbind,.) %>% 
  add_column(Group = "Post") %>% 
  pivot_longer(-Group) %>% 
  rbind(dat_cell_pre_trace, dat_cell_cond_trace,.) %>%
  mutate(Group = factor(Group, levels = c("Pre", "Cond", "Post"))) %>% 
  #pivot_longer(-Group) %>% 
  ddply(., .(name, Group), summarise, mean_acti = mean(value)) %>% 
  ddply(.,.(Group), summarise, mean= mean(mean_acti))
  aov(mean_acti~ Group,.)

t_cell_acti <- lapply(dat_cell_trace, function(x) x[[3]]) %>% 
  do.call(cbind,.) %>% 
  add_column(Group = "Post") %>% 
  pivot_longer(-Group) %>% 
  rbind(dat_cell_pre_trace, dat_cell_cond_trace,.) %>%
  mutate(Group = factor(Group, levels = c("Pre", "Cond", "Post"))) %>% 
  #pivot_longer(-Group) %>% 
  ddply(., .(name, Group), summarise, mean_acti = mean(value)) %>% 
  filter(Group != "Cond") %>% 
  wilcox.test(mean_acti~Group, ., paired = T)



summary(t_cell_acti)
TukeyHSD(t_cell_acti)  


model <-  lapply(dat_cell_trace, function(x) x[[3]]) %>% 
  do.call(cbind,.) %>% 
  add_column(Group = "Post") %>% 
  pivot_longer(-Group) %>% 
  rbind(dat_cell_pre_trace, dat_cell_cond_trace,.) %>%
  mutate(Group = factor(Group, levels = c("Pre", "Cond", "Post"))) %>% 
  aov(formula = value ~ Group + Error(name/Group))
summary(model)

setwd("~cchen/Documents/neuroscience/Pn\ project/Figure/PDF/")
cairo_pdf("p_cell_acti_pc.pdf", width = 50/25.6, height = 65/25.6, family = "Arial")
p_cell_acti
dev.off()

## d^2, cross-day aligned neurons-----
back_crossing_frame <- tibble(ID_mouse = c("m617", "m625", "m676", "m684", "m687"), Frame_d3= c(436, 938,1670, 808, 1562),
                              Frame_d7 = c(1974, 2588, 1612, 1294, 1562))
c_miniscope_matlab_ft <- function(file_trace) {
  ## import and format the data
  dat_trace1 <- raveio::read_mat(file_trace)
  ID <- str_extract(file_trace, regex("m\\d+"))
  num_compare <- c(1, 5)
  t_crossing <- c(1, 5)
  t_crossing_back <- back_crossing_frame %>% 
    dplyr::filter(ID_mouse == ID) %>% 
    dplyr::select(Frame_d3, Frame_d7) %>% 
    unlist() %>% 
    unname()
  
  cross_ID <- dat_trace1$global.map %>% 
    as_tibble() %>% 
    select(V1, V5) %>% 
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
      ddply(.,.(name), summarise,mean_cross=mean(value),sd_cross=sd(value))
    
    dat_stim_back_sta <- dat_stim_back %>% 
      as_tibble() %>% 
      add_column(Time = comp_time) %>% 
      pivot_longer(-Time) %>% 
      ddply(.,.(name), summarise,mean_back=mean(value),sd_back=sd(value))
    
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

mouse_file <- as.list(list.files("~cchen/Documents/neuroscience/Pn\ project/Data_analysis/For\ revision/PC_data2", full.names = T))

dat_cell_trace <- mapply(c_miniscope_matlab_ft, mouse_file, SIMPLIFY = F)

dat_cell_d_pre <- lapply(dat_cell_trace, function(x)  x[[1]]) %>% 
  do.call(rbind,.) %>% 
  mutate(Group = "Pre")

dat_cell_d_combine <- lapply(dat_cell_trace, function(x)  x[[2]]) %>% 
  do.call(rbind,.) %>% 
  mutate(Group = "Post") %>% 
  rbind(dat_cell_d_pre, .) %>% 
  ddply(., .(ID, Group), summarise,n=length(d2),mean=mean(d2),sd=sd(d2),se=sd(d2)/sqrt(length(d2))) %>% 
  mutate(Group = factor(Group, levels = c("Pre", "Post")))


dat_cell_d_combine_sta <- lapply(dat_cell_trace, function(x)  x[[2]]) %>% 
  do.call(rbind,.) %>% 
  mutate(Group = "Post") %>% 
  rbind(dat_cell_d_pre, .) %>% 
  ddply(., .(ID, Group), summarise,n=length(d2),value=mean(d2),sd=sd(d2),se=sd(d2)/sqrt(length(d2))) %>% 
  mutate(Group = factor(Group, levels = c("Pre", "Post"))) %>% 
  ddply(., .(Group),summarise,n=length(value),mean=mean(value),sd=sd(value),se=sd(value)/sqrt(length(value)))


p_d_combine <- ggplot(dat_cell_d_combine, aes(Group, mean, group = Group))+
  geom_boxplot(outlier.shape = NA)+
  geom_line(aes(group=ID), colour="gray90")+
  geom_jitter(aes(colour = Group, shape = Group),width = 0.2,  size=2)+
  scale_color_manual(values=c("#8491B4FF", "#3C5488FF"))+
  scale_shape_manual(values=c(16,15))+
  labs(x="", y="(d')^2")+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  scale_y_continuous(limits = c(-1, 25))+
  theme(legend.title = element_blank(), legend.position = "none")

wilcox.test(mean~Group, dat_cell_d_combine, paired = T, alternative = "less")

## compare by cells
p_d_combine_cell <- lapply(dat_cell_trace, function(x)  x[[2]]) %>% 
  do.call(rbind,.) %>% 
  mutate(Group = "Post") %>% 
  rbind(dat_cell_d_pre, .) %>% 
  mutate(Group = factor(Group, levels = c("Pre", "Post"))) %>% 
  ggplot(., aes(Group, d2, group = Group))+
  geom_violin( )+
  geom_jitter(aes(colour = Group, shape = Group),width = 0.2,  size=2, alpha = 0.3 )+
  scale_color_manual(values=c("#8491B4FF", "#3C5488FF"))+
  scale_shape_manual(values=c(16,15))+
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

t_d_combine_cell <- lapply(dat_cell_trace, function(x)  x[[2]]) %>% 
  do.call(rbind,.) %>% 
  mutate(Group = "Post") %>% 
  rbind(dat_cell_d_pre, .) %>% 
  mutate(Group = factor(Group, levels = c("Pre", "Post"))) %>% 
  wilcox.test(d2~ Group,.)

setwd("~cchen/Documents/neuroscience/Pn\ project/Figure/PDF/")
cairo_pdf("p_d_combine_pc.pdf", width = 32/25.6, height = 65/25.6, family = "Arial")
p_d_combine
dev.off()

## for comparing crossing, back-crossing and last crossing------
back_crossing_frame <- tibble(ID = c("m617", "m625", "m676", "m684", "m687", "m685"), Frame = c(1974, 2588, 1612, 1294, 2210, 2698), 
                              Last_frame= c(3016, 3168, 3106, 3180, 3726, 3640))

c_miniscope_matlab_ft <- function(file_trace) {
  ## import and format the data
  dat_trace1 <- raveio::read_mat(file_trace)
  ID_mouse <- str_extract(file_trace, regex("m\\d+"))
  
  t_crossing_back <- back_crossing_frame %>% 
    dplyr::filter(ID == ID_mouse)
  t_crossing <- c(dat_trace1[[1]][5], t_crossing_back$Frame, t_crossing_back$Last_frame)
  
  
  
  dat_stim_trace <- vector(mode = "list", length = length(t_crossing))
  
  
  for (i in seq_along(t_crossing)) {
    dat_trace <- dat_trace1[[3]][[5]] %>% 
      as_tibble() %>% 
      #apply(., 2, function(x) runmed(x, k = 41)) %>% 
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


## analyze for all mouse 
mouse_file <- as.list(list.files("~cchen/Documents/neuroscience/Pn\ project/Data_analysis/For\ revision/PC_data2", full.names = T))

dat_cell_trace <- mapply(c_miniscope_matlab_ft, mouse_file, SIMPLIFY = F)
comp_time <- seq(-2, 1.9, by=0.1)
group_day <- c("Crossing", "Crossing_back", "Last_crossing")
score_range <- range(dat_cell_trace)


for(i in c(1:3)){
  cell_order <- lapply(dat_cell_trace, function(x)  x[[1]]) %>% 
    do.call(cbind,.) %>% 
    mutate(Time = comp_time) %>% 
    pivot_longer(-Time) %>% 
    ddply(.,.(name), summarise, mean = mean(value)) %>% 
    arrange(mean)
  
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
    theme(legend.position = "none")
  assign(str_c("p_heat_", group_day[i]), p_heat)
  
}

p_heat <- plot_grid(p_heat_Crossing, p_heat_Crossing_back, p_heat_Last_crossing, nrow  = 1)

setwd("~cchen/Documents/neuroscience/Pn\ project/Figure/PDF/")
cairo_pdf("p_heat_pc_crossing.pdf", width = 140/25.6, height = 60/25.6, family = "Arial")
p_heat
dev.off()

## calculate cell activity change mouse by mouse
cc_active_fun <- function (dat_trace){
  m_id <- dat_trace[[1]] %>% 
    colnames() %>%
    .[1] %>% 
    str_extract(., regex("m\\d+"))
  
  dat_trace_combine <- dat_trace %>% 
    do.call(rbind,.) %>% 
    as_tibble() %>% 
    add_column(Time = rep(comp_time, 3)) %>% 
    add_column(Group = rep(c("Crossing", "Crossing_back", "Last_crossing"), each = length(comp_time))) %>% 
    pivot_longer(-c(Time, Group)) 
  
  dat_trace_acti <- dat_trace_combine %>% 
    ddply(., .(name, Group), summarise, mean_acti = mean(value)) %>% 
    ddply(., .(Group), summarise, cell_acti = mean(mean_acti)) %>% 
    add_column(ID = m_id)
  
  # calculate percentage of cell numer increased activity
  ratio_pre_con <- dat_trace_combine %>% 
    ddply(., .(name, Group), summarise, mean_acti = mean(value)) %>% 
    pivot_wider(., names_from = Group, values_from = mean_acti) %>% 
    mutate(diff = Crossing_back - Crossing) %>% 
    dplyr::filter(diff > 0)
  
  ratio_pre_post <- dat_trace_combine %>% 
    ddply(., .(name, Group), summarise, mean_acti = mean(value)) %>% 
    pivot_wider(., names_from = Group, values_from = mean_acti) %>% 
    mutate(diff = Last_crossing - Crossing) %>% 
    dplyr::filter(diff > 0)
  
  ratio_con_post <- dat_trace_combine %>% 
    ddply(., .(name, Group), summarise, mean_acti = mean(value)) %>% 
    pivot_wider(., names_from = Group, values_from = mean_acti) %>% 
    mutate(diff = Last_crossing - Crossing_back) %>% 
    dplyr::filter(diff > 0)
  
  dat_ratio <- c(nrow(ratio_pre_con), nrow(ratio_pre_post), nrow(ratio_con_post))/ncol(dat_trace[[1]])
  
  return(list(dat_trace_acti, dat_ratio))
}


dat_acti_combine <- mapply(cc_active_fun, dat_cell_trace, SIMPLIFY = F)

p_dat_anti_cross <- lapply(dat_acti_combine, function(x) x[[1]]) %>% 
  do.call(rbind,.) %>% 
  mutate(Group = factor(Group, levels = c("Crossing", "Crossing_back", "Last_crossing"))) %>% 
  ggplot(., aes(Group, cell_acti, group = Group))+
  geom_boxplot(outlier.shape = NA)+
  geom_line(aes(group=ID), colour="gray90")+
  geom_jitter(aes(colour = Group, shape = Group),width = 0.2,  size=2)+
  #scale_color_manual(values=c("#8491B4FF", "#00A087FF", "#3C5488FF"))+
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

setwd("~cchen/Documents/neuroscience/Pn\ project/Figure/PDF/")
cairo_pdf("p_dat_anti_cross.pdf", width = 50/25.6, height = 65/25.6, family = "Arial")
p_dat_anti_cross
dev.off()

p_test <- lapply(dat_acti_combine, function(x) x[[1]]) %>% 
  do.call(rbind,.) %>% 
  mutate(Group = factor(Group, levels = c("Crossing", "Crossing_back", "Last_crossing"))) %>%
  aov(cell_acti~ Group,.)
summary(p_test)

TukeyHSD(p_test)

## compare for each cells
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
  ddply(.,.(Group, name), summarise, mean = mean(value)) %>% 
  mutate(Group = factor(Group, levels = c("Crossing", "Crossing_back", "Last_crossing"))) %>% 
  ggplot(., aes(Group, mean, group = Group))+
  geom_violin()+
  geom_jitter(aes(colour = Group, shape = Group),width = 0.2,  size=2, alpha = 0.3)+
  #scale_color_manual(values=c("#8491B4FF", "#00A087FF", "#3C5488FF"))+
  scale_y_continuous(limits = c(-1.5, 10))+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  scale_y_continuous(limits = c(-1.5, 10))+
  theme(legend.title = element_blank(), legend.position = "none")

setwd("~cchen/Documents/neuroscience/Pn\ project/Figure/PDF/")
cairo_pdf("p_cell_firing_pc_crossing.pdf", width = 50/25.6, height = 65/25.6, family = "Arial")
p_cell_firing
dev.off()

t_cell_firing <-  lapply(dat_cell_trace, function(x)  x[[3]]) %>% 
  do.call(cbind,.) %>% 
  add_column(Group = "Last_crossing") %>% 
  rbind(p_cell_firing_corssing, p_cell_firing_corssing_back,.) %>% 
  as_tibble() %>% 
  pivot_longer(-Group) %>% 
  ddply(.,.(Group, name), summarise, mean = mean(value)) %>%
  mutate(Group = factor(Group, levels = c("Crossing", "Crossing_back", "Last_crossing"))) %>% 
  aov(mean~Group,.)

summary(t_cell_firing)


## compare Ca2+ spikes ------
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
  
  ## for the spike frequency file
  dat_spike <- list.files("~cchen/Documents/neuroscience/Pn\ project/Data_analysis/For\ revision/PC_spike3/", pattern = ID, full.names = T)
  
  dat_stim_trace <- vector(mode = "list", length = length(num_compare))
  dat_event_day <- rep(0, 3)
  
  for (i in seq_along(num_compare)) {
    global_cell <- pull(cross_ID[,i])
 
    
    dat_event <- read.xlsx(dat_spike[i], colNames = F) %>% 
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

    dat_stim1 <- dat_event[(t1_p-40):(t1_p+40-1),] 
    
    dat_stim <- aggregate(dat_stim1,list(rep(1:(nrow(dat_stim1)%/%n+1),each=n,len=nrow(dat_stim1))),mean)[-1]
 
    
    
    colnames(dat_stim) <- str_c(ID,"Cell", 1: ncol(dat_stim))
    dat_stim_trace[[i]] <- dat_stim
    
    ## calculate the event rate
    dat_event1 <- dat_event[(t1_p-40):(t1_p+40-1),] %>% 
      apply(., 2, mean) %>% 
      mean()
    
    dat_event_day[i] <- dat_event1
    
  }
  return(list(dat_stim_trace, dat_event_day))
}



## analyze for all mouse 
mouse_file <- as.list(list.files("~cchen/Documents/neuroscience/Pn\ project/Data_analysis/For\ revision/PC_data", full.names = T))

dat_cell_trace <- mapply(c_miniscope_matlab_ft, mouse_file, SIMPLIFY = F)

## plot the firing rate

p_dat_firing <- lapply(dat_cell_trace, function(x)  x[[2]]) %>% 
  do.call(rbind,.) %>% 
  as_tibble() %>% 
  rename(Pre =V1, Cond= V2, Post = V3) %>% 
  mutate(ID = sort(c("m617", "m625", "m676", "m684", "m687", "m685"))) %>% 
  pivot_longer(-ID) %>%
  mutate(name = factor(name, levels = c("Pre", "Cond", "Post"))) %>% 
  ggplot(., aes(name, value, group = name))+
  geom_boxplot(outlier.shape = NA)+
  geom_line(aes(group=ID), colour="gray90")+
  geom_jitter(aes(colour = name, shape = name),width = 0.2,  size=2)+
  scale_color_manual(values=c("#8491B4FF", "#00A087FF", "#3C5488FF"))+
  labs(x="", y="Ca2+ event rate (Hz)")+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  scale_y_continuous(limits = c(0.4, 2))+
  theme(legend.title = element_blank(), legend.position = "none")

setwd("~cchen/Documents/neuroscience/Pn\ project/Figure/PDF/")
cairo_pdf("p_dat_firing.pdf", width = 40/25.6, height = 65/25.6, family = "Arial")
p_dat_firing
dev.off()

t_dat_firing <- lapply(dat_cell_trace, function(x)  x[[2]]) %>% 
  do.call(rbind,.) %>% 
  as_tibble() %>% 
  rename(Pre =V1, Cond= V2, Post = V3) %>% 
  mutate(ID = sort(c("m617", "m625", "m676", "m684", "m687", "m685"))) %>% 
  pivot_longer(-ID) %>%
  mutate(name = factor(name, levels = c("Pre", "Cond", "Post"))) %>% 
  aov(value~name,.)

summary(t_dat_firing)
## compare each cells
p_cell_firing_corssing <- lapply(dat_cell_trace, function(x)  x[[1]]) %>% 
  lapply(., function(x) x[[1]]) %>% 
  do.call(cbind,.) %>% 
  add_column(Group = "Pre") %>% 
  pivot_longer(-Group)

p_cell_firing_corssing_back <- lapply(dat_cell_trace, function(x)  x[[1]]) %>% 
  lapply(., function(x) x[[2]]) %>% 
  do.call(cbind,.) %>% 
  add_column(Group = "Cond") %>% 
  pivot_longer(-Group)


p_cell_firing <-  lapply(dat_cell_trace, function(x)  x[[1]]) %>% 
  lapply(., function(x) x[[3]]) %>% 
  do.call(cbind,.) %>% 
  add_column(Group = "Post") %>% 
  pivot_longer(-Group) %>% 
  rbind(p_cell_firing_corssing, p_cell_firing_corssing_back,.) %>% 
  as_tibble() %>% 
  ddply(.,.(Group, name), summarise, mean = mean(value)) %>% 
  mutate(Group = factor(Group, levels = c("Pre", "Cond", "Post"))) %>% 
  ggplot(., aes(Group, mean, group = Group))+
  geom_violin()+
  geom_jitter(aes(colour = Group, shape = Group),width = 0.2,  size=2, alpha = 0.3)+
  scale_color_manual(values=c("#8491B4FF", "#00A087FF", "#3C5488FF"))+
  
  #scale_color_manual(values=c("#8491B4FF", "#00A087FF", "#3C5488FF"))+
  labs(x="", y="Ca2+ event rate (Hz)")+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  scale_y_continuous(limits = c(-0.2, 4))+
  theme(legend.title = element_blank(), legend.position = "none")

setwd("~cchen/Documents/neuroscience/Pn\ project/Figure/PDF/")
cairo_pdf("p_cell_firing.pdf", width = 40/25.6, height = 65/25.6, family = "Arial" )
p_cell_firing
dev.off()

t_cell_firing <- lapply(dat_cell_trace, function(x)  x[[1]]) %>% 
  lapply(., function(x) x[[3]]) %>% 
  do.call(cbind,.) %>% 
  add_column(Group = "Post") %>% 
  pivot_longer(-Group) %>% 
  rbind(p_cell_firing_corssing, p_cell_firing_corssing_back,.) %>% 
  as_tibble() %>% 
  ddply(.,.(Group, name), summarise, mean = mean(value)) %>% 
  aov(mean~Group, .)

summary(t_cell_firing)
## Ca2+ spike, for comparing crossing, back-crossing and last crossing------
back_crossing_frame <- tibble(ID = c("m617", "m625", "m676"), Frame = c(1974, 2588, 1612), 
                              Last_frame= c(3016, 3168, 3106))

c_miniscope_matlab_ft <- function(file_trace) {
  ## import and format the data
  dat_trace1 <- raveio::read_mat(file_trace)
  ID_mouse <- str_extract(file_trace, regex("m\\d+"))
  
  t_crossing_back <- back_crossing_frame %>% 
    dplyr::filter(ID== ID_mouse )
  t_crossing <- c(dat_trace1[[1]][5], t_crossing_back$Frame, t_crossing_back$Last_frame)
  
  dat_spike <- list.files("~cchen/Documents/neuroscience/Pn\ project/Data_analysis/For\ revision/PC_spike3/", pattern = ID_mouse, full.names = T)
  
  
  dat_stim_trace <- vector(mode = "list", length = length(t_crossing))
  dat_event_day <- rep(0, 3)
  
  
  for (i in seq_along(t_crossing)) {
 
    
    dat_event <- read.xlsx(dat_spike[3], colNames = F) %>% 
      as_tibble() 
    
    ## number of rows to be binned
    n <- 2 # 0.05*2=0.1
    t1_p <- t_crossing[i]
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
 
    dat_stim1 <- dat_event[(t1_p-40):(t1_p+40-1),] 
    
    dat_stim <- aggregate(dat_stim1,list(rep(1:(nrow(dat_stim1)%/%n+1),each=n,len=nrow(dat_stim1))),mean)[-1]
    # apply(., 2, scale) %>% 
    # as_tibble() %>% 
    # replace(is.na(.), 0)
    
    
    colnames(dat_stim) <- str_c(ID_mouse,"Cell", 1: ncol(dat_stim))
    dat_stim_trace[[i]] <- dat_stim
    
    ## calculate the event rate
    dat_event1 <- dat_event[(t1_p-40):(t1_p+40-1),] %>% 
      apply(., 2, mean) %>% 
      mean()
    
    dat_event_day[i] <- dat_event1
  }
  return(list(dat_stim_trace, dat_event_day))
}


## analyze for all mouse 
mouse_file <- as.list(list.files("~cchen/Documents/neuroscience/Pn\ project/Data_analysis/For\ revision/PC_data2", full.names = T))

dat_cell_trace <- mapply(c_miniscope_matlab_ft, mouse_file, SIMPLIFY = F)


p_dat_firing_d7 <- lapply(dat_cell_trace, function(x)  x[[2]]) %>% 
  do.call(rbind,.) %>% 
  as_tibble() %>% 
  rename(Crossing =V1, Crossing_back= V2, Last_crossing = V3) %>% 
  mutate(ID = sort(c("m617", "m625", "m676", "m684", "m687"))) %>% 
  pivot_longer(-ID) %>%
  mutate(name = factor(name, levels = c("Crossing", "Crossing_back", "Last_crossing"))) %>% 
  ggplot(., aes(name, value, group = name))+
  geom_boxplot(outlier.shape = NA)+
  geom_line(aes(group=ID), colour="gray90")+
  geom_jitter(aes(colour = name, shape = name),width = 0.2,  size=2)+
  #scale_color_manual(values=c("#8491B4FF", "#00A087FF", "#3C5488FF"))+
  labs(x="", y="Ca2+ event rate (Hz)")+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  scale_y_continuous(limits = c(-0.1, 0.6))+
  theme(legend.title = element_blank(), legend.position = "none")




## plot the trace to show PC activity----
cell_order <- dat_cell_trace[[2]][[3]] %>% 
  apply(., 2, mean) 

cell_order_ID <- cell_order[order(cell_order)] %>% 
  names() %>% 
  tail(.,10) %>% 
  regmatches(., regexpr("cell(\\d+)", ., ignore.case = TRUE)) %>% 
  sub("cell", "", ., ignore.case = TRUE) %>% 
  as.numeric()



dat_trace1 <- as.list(list.files("~cchen/Documents/neuroscience/Pn\ project/Data_analysis/For\ revision/PC_data", full.names = T))[[5]] %>% 
  raveio::read_mat(.)
  




p_pc_trace <- dat_trace1[[3]][[4]] %>% 
  as_tibble() %>% 
  select(sample(ncol(.), min(10, ncol(.)))) %>% 
  apply(., 2, scale) %>% 
  as_tibble() %>% 
  rowid_to_column(., "Time") %>% 
  filter(Time < 1000) %>% 
  pivot_longer(-Time) %>% 
  ggplot(., aes(Time, value, color = name))+
  geom_line()+
  #scale_y_continuous(limits = plot_range)+
  facet_grid(rows = vars(name))+
  #theme_void()+
  theme(legend.position = "none")+
  theme(strip.text.y = element_blank()) 


dat_cell_trace_cluster_num <- tibble(name = names(dat_cell_trace_cluster), cluster = dat_cell_trace_cluster) %>% 
  as_tibble() %>% 
  mutate(cell_ID =as.numeric(str_extract(name, "(?<=Cell)\\d+"))) %>% 
  mutate(mouse_ID = stringr::str_extract(name, ".*(?=Cell)")) %>% 
  select(-name) 

dat_global_cell <- vector(mode="list", length = length(mouse_file))

for (i in 1: length(mouse_file)) {
  file_trace <- mouse_file[[i]]
  dat_trace1 <- raveio::read_mat(file_trace)
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
  full_join(., dat_cell_trace_cluster_num, by = join_by(mouse_ID, cell_ID)) %>% 
  filter(mouse_ID == "m684")

## plot the trace of m676 from Pre, Cond and post
dat_trace1 <- as.list(list.files("~cchen/Documents/neuroscience/Pn\ project/Data_analysis/For\ revision/PC_data", full.names = T))[[4]] %>% 
  raveio::read_mat(.)

day_compare <- c(1, 4, 5)
for (i in seq_along(day_compare) ) {
  
  t1_p <- dat_trace1$crossing[day_compare[i]]
  
  plot_range <- c(1:1000)
  
  plot_cell_ID <- dat_global_cell_combine %>% 
    arrange(cluster) %>% 
    select(i) %>% 
    pull() %>% 
    .[1:12]


  
  p_pc_trace <- dat_trace1[[3]][[day_compare[i]]] %>% 
    as_tibble() %>% 
    select(all_of(plot_cell_ID)) %>% 
    apply(., 2, scale) %>% 
    as_tibble() %>% 
    slice(plot_range) %>% 
    rowid_to_column(., "Time") %>% 
    pivot_longer(-Time) %>% 
    mutate(name = factor(name, levels = c(str_c("V",plot_cell_ID)))) %>% 
    ggplot(., aes(Time, value, color = name))+
    geom_line()+
    #scale_y_continuous(limits = plot_range)+
    facet_grid(rows = vars(name))+
    theme_void()+
    theme(legend.position = "none")+
    theme(strip.text.y = element_blank()) +
    annotate("rect", xmin = t1_p-40, xmax = t1_p + 40, ymin = -Inf, ymax = Inf, alpha = .2,fill = "blue")
  
  assign(str_c("p_trace", day_compare[i]), p_pc_trace)
  
  
}


p_trace <- plot_grid(p_trace1, p_trace4, p_trace5, nrow = 1)






setwd("~cchen/Documents/neuroscience/Pn\ project/Figure/PDF/")
cairo_pdf("p_PC_trace.pdf", width = 60/25.6, height = 65/25.6, family = "Arial")
p_trace
dev.off()

## for comparing crossing, back-crossing and last crossing------
c_miniscope_matlab_ft <- function(file_trace) {
  ## import and format the data
  dat_trace1 <- raveio::read_mat(file_trace)
  ID_mouse <- str_extract(file_trace, regex("m\\d+"))
  
  
  length_post <- dat_trace1[[3]][[5]] %>% 
    as_tibble() %>% 
    nrow()
  t_crossing_post <- sample(c(40:(length_post-40)), 100)
  t_crossing_back <- sample(c(40:(length_post-40)), 100)
  t_crossing_last <- sample(c(40:(length_post-40)), 100)
  
  t_crossing <- tibble(cross = t_crossing_post, cross_back = t_crossing_back, last_cross = t_crossing_last)
  group_day <- c("Crossing", "Crossing_back", "Last_crossing")
  
  
  
  
  dat_stim_trace <- vector(mode="list", length = 3)
  
  
  for (i in c(1:3)) {
    dat_trace <- dat_trace1[[3]][[5]] %>% 
      as_tibble() %>% 
      apply(., 2, scale)
    
    n <- 2 # 0.05*2=0.1
    
    mean_cell_acti <- vector(mode="list", length = 100)
    
    for (j in c(1:100)) {
      t1_p <- t_crossing[j,i] %>% 
        unlist()
      
      dat_stim1 <- dat_trace[(t1_p-40):(t1_p+40-1),] 
      
      dat_stim <- aggregate(dat_stim1,list(rep(1:(nrow(dat_stim1)%/%n+1),each=n,len=nrow(dat_stim1))),mean)[-1]
      # apply(., 2, scale) %>% 
      # as_tibble() %>% 
      # replace(is.na(.), 0)
      
      mean_cell_acti[[j]] <- dat_stim %>% 
        apply(., 2, mean) 
      
    }
    
    dat_stim_trace[[i]] <- mean_cell_acti %>%
      do.call(rbind, .) %>%
      apply(., 2, mean) %>%
      tibble::enframe(name = "variable", value = "mean") %>%
      pivot_longer(cols = mean, names_to = "metric", values_to = "value") %>% 
      mutate(Group = group_day[i]) %>% 
      select(-metric)
    
  }
  
  dat_stim_trace_cm <- dat_stim_trace %>% 
    do.call(rbind,.)
  return(dat_stim_trace_cm)
}


mouse_file <- as.list(list.files("~cchen/Documents/neuroscience/Pn\ project/Data_analysis/For\ revision/PC_data2", full.names = T))

dat_cell_trace <- mapply(c_miniscope_matlab_ft, mouse_file, SIMPLIFY = F)




p_dat_anti_cross <- dat_cell_trace%>% 
  do.call(rbind,.) %>% 
  as_tibble() %>% 
  mutate(Group = factor(Group, levels = c("Crossing", "Crossing_back", "Last_crossing"))) %>% 
  ggplot(., aes(Group, value, color= Group))+
  geom_violin()+
  geom_jitter(aes(colour = Group, shape = Group),width = 0.2,  size=2)+
  #scale_color_manual(values=c("#8491B4FF", "#00A087FF", "#3C5488FF"))+
  labs(x="", y="Population activity (∆F/F)")+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  scale_y_continuous(limits = c(-1.5, 10))+
  theme(legend.title = element_blank(), legend.position = "none")

setwd("~cchen/Documents/neuroscience/Pn\ project/Figure/PDF/")
cairo_pdf("p_dat_anti_cross.pdf", width = 50/25.6, height = 65/25.6, family = "Arial")
p_dat_anti_cross
dev.off()  

## random crossing of PC neurons, placebo state---------
c_miniscope_matlab_ft <- function(file_trace) {
  ## import and format the data
  dat_trace1 <- raveio::read_mat(file_trace)
  ID <- str_extract(file_trace, regex("m\\d+"))
  num_compare <- c(1,4, 5)
  
  length_pre <- dat_trace1$traces[[1]] %>% 
    nrow()
  t_crossing_pre <- sample(c(40:(length_pre-40)), 100)
  
  length_cond <- dat_trace1$traces[[4]] %>% 
    nrow()
  t_crossing_cond <- sample(c(40:(length_cond-40)), 100)
  
  length_post <- dat_trace1$traces[[5]] %>% 
    nrow()
  t_crossing_post <- sample(c(40:(length_post-40)), 100)
  
  t_crossing <- list(t_crossing_pre, t_crossing_cond, t_crossing_post)
  
  cross_ID <- dat_trace1$global.map %>% 
    as_tibble() %>% 
    select(V1, V4, V5) %>% 
    mutate_all(na_if, 0) %>% 
    drop_na()
  
  dat_stim_trace <- vector(mode = 'list', length = length(num_compare))
  
  
  for (i in seq_along(num_compare)) {
    global_cell <- pull(cross_ID[,i])
    dat_trace <- dat_trace1$traces[[num_compare[i]]] %>% 
      as_tibble() %>% 
      #select(all_of(global_cell)) %>% 
      #apply(., 2, function(x) runmed(x, k = 41)) %>% 
      apply(., 2, scale)
    
    ## number of rows to be binned
    n <- 2 # 0.05*2=0.1
    
    t_crossing_day <- t_crossing[[i]]
    t_crossing_mean <- c()
    
    for (j in seq_along(t_crossing_day)){
      t1_p <- t_crossing_day[j]
      
      dat_stim1 <- dat_trace[(t1_p-40):(t1_p+40-1),] 
      
      dat_stim <- aggregate(dat_stim1,list(rep(1:(nrow(dat_stim1)%/%n+1),each=n,len=nrow(dat_stim1))),mean)[-1] %>% 
        apply(., 2, mean) 
      
      t_crossing_mean <- cbind(t_crossing_mean,dat_stim)   
      
      
    }
    
    dat_stim_trace[[i]] <- t_crossing_mean %>% 
      as_tibble() %>% 
      apply(., 1, mean)
  }
  return(dat_stim_trace)
}

mouse_file <- as.list(list.files("~cchen/Documents/neuroscience/Pn\ project/Data_analysis/For\ revision/PC_data2", full.names = T))

dat_cell_trace <- mapply(c_miniscope_matlab_ft, mouse_file, SIMPLIFY = F)

dat_cell_pre_trace <- lapply(dat_cell_trace, function(x) x[[1]]) %>% 
  unlist() 

dat_cell_cond_trace <- lapply(dat_cell_trace, function(x) x[[2]]) %>% 
  unlist()

dat_cell_post_trace <- lapply(dat_cell_trace, function(x) x[[3]]) %>% 
  unlist()
p_cell_acti <- c(dat_cell_pre_trace, dat_cell_cond_trace, dat_cell_post_trace) %>% 
  as_tibble() %>% 
  add_column(Group = c(rep("Pre", length(dat_cell_pre_trace)), rep("Cond", length(dat_cell_cond_trace)), rep("Post", length(dat_cell_post_trace)))) %>% 
  mutate(Group = factor(Group, levels = c("Pre", "Cond", "Post"))) %>% 
  pivot_longer(-Group) %>% 
  ggplot(., aes(Group, value, group = Group))+
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

setwd("~cchen/Documents/neuroscience/Pn\ project/Figure/PDF/")
cairo_pdf("p_cell_anti_pc.pdf", width = 50/25.6, height = 65/25.6, family = "Arial")
p_cell_acti
dev.off()


## for random crossing control---------
# 1. for the d^2

c_miniscope_matlab_ft <- function(file_trace) {
  ## import and format the data
  dat_trace1 <- raveio::read_mat(file_trace)
  ID <- str_extract(file_trace, regex("m\\d+"))
  num_compare <- c(1, 5)
  length_pre <- dat_trace1$traces[[1]]  %>% 
    nrow()
  t_crossing_pre <- tibble(cross = sample(c(40:(length_pre-40)), 100), back = sample(c(40:(length_pre-40)), 100))
  
  length_post <- dat_trace1$traces[[5]]  %>% 
    nrow()
  t_crossing_post <- tibble(cross = sample(c(40:(length_post-40)), 100), back = sample(c(40:(length_post-40)), 100))
  
  t_crossing <- list(t_crossing_pre,  t_crossing_post)
  
  cross_ID <- dat_trace1$global.map %>% 
    as_tibble() %>% 
    dplyr::select(V1, V5) %>% 
    mutate_all(na_if, 0) %>% 
    drop_na()
  
  group <- c("Pre", "Post")
  
  dat_stim_trace <- rep(0, length(num_compare))
  dat_stim_trace_cell <- vector(mode = "list", length = length(num_compare))
  
  
  for (i in seq_along(num_compare)) {
    
    global_cell <- pull(cross_ID[,i])
    dat_trace <- dat_trace1$traces[[num_compare[i]]] %>% 
      as_tibble() %>% 
      #dplyr::select(all_of(global_cell)) %>% 
      #apply(., 2, function(x) runmed(x, k = 41)) %>% 
      apply(., 2, scale)
    
    dat_stim_combine_mean <- rep(0, 100)
    dat_stim_combine_cell <- NULL
    ## number of rows to be binned
    n <- 2 # 0.05*2=0.1
    
    for (j in c(1:100)) {
      t1_p <- unlist(t_crossing[[i]][j,1])
      t1_p_back <- unlist(t_crossing[[i]][j,2])
      
      dat_stim1 <- dat_trace[(t1_p-40):(t1_p+40-1),] 
      
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
        ddply(.,.(name), summarise,mean_cross=mean(value),sd_cross=sd(value))
      
      dat_stim_back_sta <- dat_stim_back %>% 
        as_tibble() %>% 
        add_column(Time = comp_time) %>% 
        pivot_longer(-Time) %>% 
        ddply(.,.(name), summarise,mean_back=mean(value),sd_back=sd(value))
      
      dat_stim_combine1 <- full_join(dat_stim_cross_sta, dat_stim_back_sta, by = 'name') %>% 
        mutate(sd_pool = sqrt(sd_cross^2 + sd_back^2)/sqrt(2)) %>% 
        mutate(sd_pool = ifelse(sd_pool < 1e-06, 1e-06, sd_pool)) %>% 
        mutate(d = (mean_cross - mean_back)/sd_pool) %>% 
        mutate(d2 = d^2) %>% 
        dplyr::select(name, d2) %>% 
        mutate(ID = ID) %>% 
        mutate(Group = j)
      
      dat_stim_combine_cell <- rbind(dat_stim_combine_cell, dat_stim_combine1)
      
      dat_stim_combine_mean[j] <- dat_stim_combine1$d2 %>% 
        mean()
      
    }
    
    
    dat_stim_trace[i] <- dat_stim_combine_mean %>% 
      mean()
    
    dat_stim_trace_cell[[i]] <- dat_stim_combine_cell %>% 
      ddply(.,.(name), summarise, value = mean(d2)) %>% 
      mutate(Group = group[i])
    
  }
  return(list(dat_stim_trace, dat_stim_trace_cell))
}

mouse_file <- as.list(list.files("~cchen/Documents/neuroscience/Pn\ project/Data_analysis/For\ revision/PC_data2", full.names = T))

dat_cell_trace <- mapply(c_miniscope_matlab_ft, mouse_file, SIMPLIFY = F)


dat_cell_d_pre <- lapply(dat_cell_trace, function(x)  x[[2]]) %>% 
  lapply(., function(x)  x[[2]]) %>% 
  do.call(rbind,.) %>% 
  mutate(Group = "Pre")

p_d_combine_cell <- lapply(dat_cell_trace, function(x)  x[[2]]) %>% 
  lapply(., function(x)  x[[2]]) %>% 
  do.call(rbind,.) %>% 
  mutate(Group = "Post") %>% 
  rbind(dat_cell_d_pre, .) %>% 
  mutate(Group = factor(Group, levels = c("Pre", "Post"))) %>% 
  ggplot(., aes(Group, value, group = Group))+
  geom_violin( )+
  geom_jitter(aes(colour = Group, shape = Group),width = 0.2,  size=2, alpha = 0.3 )+
  scale_color_manual(values=c("#8491B4FF", "#3C5488FF"))+
  scale_shape_manual(values=c(16,15))+
  labs(x="", y="(d')^2" )+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  scale_y_continuous(limits = c(-1, 150))+
  theme(legend.title = element_blank(), legend.position = "none")

t_d_combine_cell <- lapply(dat_cell_trace, function(x)  x[[2]]) %>% 
  lapply(., function(x)  x[[2]]) %>% 
  do.call(rbind,.) %>% 
  mutate(Group = "Post") %>% 
  rbind(dat_cell_d_pre, .) %>% 
  mutate(Group = factor(Group, levels = c("Pre", "Post"))) %>% 
  wilcox.test(value~Group, ., paired = T, alternative = "less")

setwd("~cchen/Documents/neuroscience/Pn\ project/Figure/PDF/")
cairo_pdf("p_d_combine_pc_random.pdf", width = 40/25.6, height = 65/25.6, family = "Arial")
p_d_combine_cell
dev.off()




## plot the cell number detected every day

## plot the label cells

num_m617 <- c(98, 95, 92, 41)
num_m625 <- c(74, 69, 79, 23)
num_m676 <- c(135, 159, 145, 32)
num_m684 <- c(111, 71, 77, 36)
num_m687 <- c(150, 140, 132, 47)

dat_num <- tibble(m617 = num_m617, m625=num_m625,m676 = num_m676, m684 = num_m684, m687 = num_m687) %>% 
  add_column(Group = c("Pre", "Cond", "Post", "Align")) %>%
  pivot_longer(-Group) %>% 
  mutate(Group = factor(Group, levels = c("Pre", "Cond", "Post", "Align"))) %>% 
  ggplot(., aes(Group, value, group = Group))+
  geom_boxplot(outlier.shape = NA)+
  #geom_line(aes(group=ID), colour="gray90")+
  geom_jitter(aes(colour = Group, shape = Group),width = 0.2,  size=2)+
  scale_color_manual(values=c("#8491B4FF", "#00A087FF", "#3C5488FF", "violetred3"))+
  labs(x="", y="Number of detected cells")+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  scale_y_continuous(limits = c(20, 200))+
  theme(legend.title = element_blank(), legend.position = "none")

setwd("~cchen/Documents/neuroscience/Pn\ project/Figure/PDF/")
cairo_pdf("p_number.pdf", width = 55/25.6, height = 65/25.6, family = "Arial")
dat_num
dev.off()

## compare the neuron activity and their time of crossing------
c_miniscope_matlab_ft <- function(file_trace) {
  ## import and format the data
  dat_trace1 <- raveio::read_mat(file_trace)
  ID <- str_extract(file_trace, regex("m\\d+"))
  num_compare <- c(1,3,4, 5)
  t_crossing <- c(1,3, 4, 5)
  
  cross_ID <- dat_trace1$global.map %>% 
    as_tibble() %>% 
    select(V1, V3, V4, V5) %>% 
    mutate_all(na_if, 0) %>% 
    drop_na()
  
  group_day <- c("Pre","D5", "Cond", "Post")
  
  dat_stim_trace <- vector(mode = "list", length = length(num_compare))
  dat_t_crossing <- tibble(Crossing = dat_trace1$crossing[t_crossing], Group = group_day, ID = ID)
  
  for (i in seq_along(num_compare)) {
    global_cell <- pull(cross_ID[,i])
    dat_trace <- dat_trace1$traces[[num_compare[i]]] %>% 
      as_tibble() %>% 
      #select(all_of(global_cell)) %>% 
      #apply(., 2, function(x) runmed(x, k = 41)) %>% 
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
  return(list(dat_stim_trace, dat_t_crossing))
}



mouse_file <- as.list(list.files("~cchen/Documents/neuroscience/Pn\ project/Data_analysis/For\ revision/PC_data2", full.names = T))

dat_cell_trace <- mapply(c_miniscope_matlab_ft, mouse_file, SIMPLIFY = F)

dat_cell_trace1 <- dat_cell_trace %>% 
  lapply(., function(x) x[[1]])
  
comp_time <- seq(-2, 1.9, by=0.1)
group_day <- c("Pre","D5", "Cond", "Post")
score_range <- range(dat_cell_trace1)

for(i in c(1:4)){
  cell_order <- lapply(dat_cell_trace1, function(x)  x[[i]]) %>% 
    do.call(cbind,.) %>% 
    mutate(Time = comp_time) %>% 
    pivot_longer(-Time) %>% 
    ddply(.,.(name), summarise, mean = mean(value)) %>% 
    arrange(mean)
  
  p_heat <- lapply(dat_cell_trace1, function(x)  x[[i]]) %>% 
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
    theme(legend.position = "none")
  assign(str_c("p_heat_", group_day[i]), p_heat)
  
}

p_heat <- plot_grid(p_heat_Pre, p_heat_Cond, p_heat_Post, nrow  = 1)

setwd("~cchen/Documents/neuroscience/Pn\ project/Figure/PDF/")
cairo_pdf("p_heat_pc.pdf", width = 140/25.6, height = 60/25.6, family = "Arial")
p_heat
dev.off()

## for individual mouse
dat_cell_pre_trace <- lapply(dat_cell_trace1, function(x) x[[1]]) %>% 
  do.call(cbind,.) %>% 
  add_column(Group = "Pre") %>% 
  pivot_longer(-Group) %>% 
  mutate(ID = str_extract(name, "m\\d+")) 

dat_cell_d5_trace <- lapply(dat_cell_trace1, function(x) x[[2]]) %>% 
  do.call(cbind,.) %>% 
  add_column(Group = "D5") %>% 
  pivot_longer(-Group) %>% 
  mutate(ID = str_extract(name, "m\\d+")) 

dat_cell_cond_trace <- lapply(dat_cell_trace1, function(x) x[[3]]) %>% 
  do.call(cbind,.) %>% 
  add_column(Group = "Cond") %>% 
  pivot_longer(-Group) %>% 
  mutate(ID = str_extract(name, "m\\d+")) 

dat_cell_ca_trace <- lapply(dat_cell_trace1, function(x) x[[4]]) %>% 
  do.call(cbind,.) %>% 
  add_column(Group = "Post") %>% 
  pivot_longer(-Group) %>% 
  mutate(ID = str_extract(name, "m\\d+")) %>% 
  rbind(dat_cell_pre_trace, dat_cell_d5_trace, dat_cell_cond_trace,.) %>% 
  mutate(Group = factor(Group, levels = c("Pre","D5", "Cond", "Post"))) %>% 
  ddply(., .(name, Group, ID), summarise, mean_acti = mean(value)) %>% 
  ddply(.,.(Group, ID), summarise, mean_acti1 = mean(mean_acti), n= length(mean_acti)) %>% 
  as_tibble()


dat_t_crossing <- lapply(dat_cell_trace, function(x) x[[2]]) %>% 
  do.call(rbind,.) %>% 
  mutate(Crossing = Crossing/10) %>% 
  mutate(Group = factor(Group, levels = c("Pre","D5", "Cond", "Post")))
  
  




p_Ca_cross_correlation <- left_join(dat_cell_ca_trace, dat_t_crossing) %>% 
  filter(Crossing < 200) %>% 
  ggplot(., aes(mean_acti1, Crossing))+
  geom_point(aes(colour = Group),  size=2)+
  geom_smooth(method = "lm", se=F)+
  labs(x="Mean activity of rACC-Pn neurons (s.d.)", y="Latency of 1st crossing (s)")+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  #scale_y_continuous(expand = c(0, 0), limits = c(0, 100)) +
  theme(legend.position = 'none')


cor.test(p_Ca_cross_correlation$mean_acti1, p_Ca_cross_correlation$Crossing)  

setwd("~cchen/Documents/neuroscience/Pn\ project/Figure/PDF/")
cairo_pdf("p_Ca_cross_correlation.pdf", width = 90/25.6, height = 60/25.6, family = "Arial")
p_Ca_cross_correlation
dev.off()


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
      ddply(.,.(name), summarise,mean_cross=mean(value),sd_cross=sd(value))
    
    dat_stim_back_sta <- dat_stim_back %>% 
      as_tibble() %>% 
      add_column(Time = comp_time) %>% 
      pivot_longer(-Time) %>% 
      ddply(.,.(name), summarise,mean_back=mean(value),sd_back=sd(value))
    
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
  filter(Group != "Cond") %>% 
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

setwd("~cchen/Documents/neuroscience/Pn\ project/Figure/PDF/")
cairo_pdf("p_d_combine.pdf", width = 45/25.6, height = 65/25.6, family = "Arial")
p_d_combine_cell
dev.off()


## do the clustering analysis crossing difference phases of PAC-----

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
  full_join(., dat_cell_trace_cluster_num, by = join_by(mouse_ID, cell_ID)) %>% 
  filter(cluster ==1)

write.csv(dat_global_cell_combine, "global_cell_cluster.csv")
  
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
  filter(value > 3) %>% 
  ggplot(., aes(Group, value, col = Group))+
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
## for d^2

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
      ddply(.,.(name), summarise,mean_cross=mean(value),sd_cross=sd(value))
    
    dat_stim_back_sta <- dat_stim_back %>% 
      as_tibble() %>% 
      add_column(Time = comp_time) %>% 
      pivot_longer(-Time) %>% 
      ddply(.,.(name), summarise,mean_back=mean(value),sd_back=sd(value))
    
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
  rbind(dat_cell_d_pre,dat_cell_d_cond,.) %>% 
  mutate(cluster = map_dbl(name, ~dat_cell_trace_cluster[.x])) %>% 
  mutate(Group = factor(Group, levels = c("Pre", "Cond","Post")))


dat_cell_d_combine_sta <- lapply(dat_cell_trace, function(x)  x[[3]]) %>% 
  do.call(rbind,.) %>% 
  mutate(Group = "Post") %>% 
  rbind(dat_cell_d_pre,dat_cell_d_cond,.) %>% 
  ddply(., .(ID, Group), summarise,n=length(d2),value=mean(d2),sd=sd(d2),se=sd(d2)/sqrt(length(d2))) %>% 
  mutate(Group = factor(Group, levels = c("Pre", "Cond","Post")))
ddply(., .(Group),summarise,n=length(value),mean=mean(value),sd=sd(value),se=sd(value)/sqrt(length(value)))



## compare by cells
p_d_combine_cell <- dat_cell_d_combine %>% 
  filter(cluster == 2) %>% 
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
  scale_y_continuous(limits = c(-1, 450))+
  theme(legend.title = element_blank(), legend.position = "none")


setwd("~cchen/Documents/neuroscience/Pn\ project/Figure/PDF/")
cairo_pdf("p_d_combine_cell.pdf", width = 50/25.6, height = 60/25.6, family = "Arial" )
p_d_combine_cell
dev.off()
t_d_combine_cell <- dat_cell_d_combine %>% 
  filter(cluster==1) %>% 
  filter(Group != "Post") %>% 
  wilcox.test(d2~ Group,., paired = T)


## for comparing crossing, back-crossing and last crossing------
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


## compare the firing frequency during crossing and crossing back
back_crossing_frame <- tibble(ID = c("m617", "m625", "m676", "m684", "m687", "m685"), Frame = c(1974, 2588, 1612, 1294, 2210, 2698), 
                              Last_frame= c(3016, 3168, 3106, 3180, 3726, 3640))

c_miniscope_matlab_ft <- function(file_trace) {
  ## import and format the data
  dat_trace1 <- raveio::read_mat(file_trace)
  ID <- str_extract(file_trace, regex("m\\d+"))
  num_compare <- c(1,4, 5)

  t_crossing <- c(dat_trace1[[1]][5], back_crossing_frame$Frame, back_crossing_frame$Last_frame)
  
  cross_ID <- dat_trace1$global.map %>% 
    as_tibble() %>% 
    select(V1, V4, V5) %>% 
    mutate_all(na_if, 0) %>% 
    drop_na() %>% 
    select(V1, V5) %>% 
    pull(V5)
  
  ## for the spike frequency file
  dat_spike <- list.files("~cchen/Documents/neuroscience/Pn\ project/Data_analysis/For\ revision/PC_spike3/", pattern = ID, full.names = T)
  
  dat_stim_trace <- vector(mode = "list", length = length(num_compare))
  dat_event_day <- rep(0, 3)
  
  for (i in seq_along(num_compare)) {

    
    dat_event <- read.xlsx(dat_spike[3], colNames = F) %>% 
      as_tibble() %>% 
      select(all_of(cross_ID))
    
    ## number of rows to be binned
    n <- 2 # 0.05*2=0.1
    t1_p <- t_crossing[i]
    
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
    
    dat_stim1 <- dat_event[(t1_p-40):(t1_p+40-1),] 
    
    dat_stim <- aggregate(dat_stim1,list(rep(1:(nrow(dat_stim1)%/%n+1),each=n,len=nrow(dat_stim1))),mean)[-1]
    
    
    
    colnames(dat_stim) <- str_c(ID,"Cell", 1: ncol(dat_stim))
    dat_stim_trace[[i]] <- dat_stim
    
    ## calculate the event rate
    dat_event1 <- dat_event[(t1_p-40):(t1_p+40-1),] %>% 
      apply(., 2, mean) %>% 
      mean()
    
    dat_event_day[i] <- dat_event1
    
  }
  return(dat_stim_trace)
}


mouse_file <- as.list(list.files("~cchen/Documents/neuroscience/Pn\ project/Data_analysis/For\ revision/PC_data", full.names = T))

dat_cell_trace <- mapply(c_miniscope_matlab_ft, mouse_file, SIMPLIFY = F)


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
  filter(cluster ==2) %>% 
  ggplot(., aes(Group, mean, group = Group, color = Group))+
  geom_violin()+
  #geom_boxplot(outlier.shape = NA)+
  #geom_line(aes(group=name), colour="gray90")+
  geom_jitter(aes(colour = Group, shape = Group),width = 0.2,  size=2, alpha= 0.3)+
  #scale_color_manual(values=c("#8491B4FF", "#00A087FF", "#3C5488FF"))+
  #facet_grid(rows = vars(cluster))+
  labs(x="", y="Ca2+ spike rate (Hz)")+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  scale_y_continuous(limits = c(-0.2, 4))+
  theme(legend.position = 'none')

t_spike_crossing <- lapply(dat_cell_trace, function(x)  x[[3]]) %>% 
  do.call(cbind,.) %>% 
  add_column(Group = "Last_crossing") %>% 
  rbind(p_cell_firing_corssing, p_cell_firing_corssing_back,.) %>% 
  as_tibble() %>% 
  pivot_longer(-Group) %>% 
  mutate(cluster = map_dbl(name, ~dat_cell_trace_cluster[.x])) %>% 
  mutate(Group = factor(Group, levels = c("Crossing", "Crossing_back", "Last_crossing"))) %>% 
  ddply(.,.(name, Group, cluster), summarise, mean = mean(value)) %>% 
  filter(cluster ==1) %>% 
  aov(mean~Group,.)

summary(t_spike_crossing)
TukeyHSD(t_spike_crossing)

setwd("~cchen/Documents/neuroscience/Pn\ project/Figure/PDF/")
cairo_pdf("p_cell_firing.pdf", width = 45/25.6, height = 60/25.6, family = "Arial" )
p_cell_firing
dev.off()

## for spike frequency analysis----
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
  
  ## for the spike frequency file
  dat_spike <- list.files("~cchen/Documents/neuroscience/Pn\ project/Data_analysis/For\ revision/PC_spike3/", pattern = ID, full.names = T)
  
  dat_stim_trace <- vector(mode = "list", length = length(num_compare))
  dat_event_day <- rep(0, 3)
  
  for (i in seq_along(num_compare)) {
    global_cell <- pull(cross_ID[,i])
    
    
    dat_event <- read.xlsx(dat_spike[i], colNames = F) %>% 
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
    
    dat_stim1 <- dat_event[(t1_p-40):(t1_p+40-1),] 
    
    dat_stim <- aggregate(dat_stim1,list(rep(1:(nrow(dat_stim1)%/%n+1),each=n,len=nrow(dat_stim1))),mean)[-1]
    
    
    
    colnames(dat_stim) <- str_c(ID,"Cell", 1: ncol(dat_stim))
    dat_stim_trace[[i]] <- dat_stim
    
    ## calculate the event rate
    dat_event1 <- dat_event[(t1_p-40):(t1_p+40-1),] %>% 
      apply(., 2, mean) %>% 
      mean()
    
    dat_event_day[i] <- dat_event1
    
  }
  return(list(dat_stim_trace, dat_event_day))
}


mouse_file <- as.list(list.files("~cchen/Documents/neuroscience/Pn\ project/Data_analysis/For\ revision/PC_data", full.names = T))

dat_cell_trace <- mapply(c_miniscope_matlab_ft, mouse_file, SIMPLIFY = F)

## plot the firing rate
p_cell_firng_pre <- lapply(dat_cell_trace, function(x)  x[[1]]) %>% 
  lapply(., function(x) x[[1]]) %>% 
  do.call(cbind,.) %>% 
  add_column(Group = "Pre") %>% 
  add_column(Time = comp_time) %>% 
  pivot_longer(-c(Time, Group))

p_cell_firing_cond<- lapply(dat_cell_trace, function(x)  x[[1]]) %>% 
  lapply(., function(x) x[[2]]) %>% 
  do.call(cbind,.) %>% 
  add_column(Group = "Cond") %>% 
  add_column(Time = comp_time) %>% 
  pivot_longer(-c(Time, Group))

dat_cell_spike_combine <- lapply(dat_cell_trace, function(x)  x[[1]]) %>% 
  lapply(., function(x) x[[3]]) %>% 
  do.call(cbind,.) %>% 
  add_column(Group = "Post") %>% 
  add_column(Time = comp_time) %>% 
  pivot_longer(-c(Time, Group)) %>% 
  rbind(p_cell_firng_pre, p_cell_firing_cond,.) %>% 
  as_tibble()


p_cell_firing <- dat_cell_spike_combine %>% 
  ddply(.,.(Group, name), summarise, mean = mean(value)) %>% 
  mutate(cluster = map_dbl(name, ~dat_cell_trace_cluster[.x])) %>% 
  mutate(Group = factor(Group, levels = c("Pre", "Cond", "Post"))) %>% 
  filter(cluster ==1) %>% 
  ggplot(., aes(Group, mean, group = Group))+
  geom_violin()+
  #geom_boxplot(outlier.shape = NA)+
  geom_jitter(aes(colour = Group, shape = Group),width = 0.2,  size=2, alpha= 0.3)+
  scale_color_manual(values=c("#8491B4FF", "#00A087FF", "#3C5488FF"))+
  labs(x="", y="Ca2+ spike rate (Hz)")+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  scale_y_continuous(limits = c(-0.2, 5))+
  theme(legend.title = element_blank(), legend.position = "none")

setwd("~cchen/Documents/neuroscience/Pn\ project/Figure/PDF/")
cairo_pdf("p_cell_firing.pdf", width = 40/25.6, height = 60/25.6, family = "Arial" )
p_cell_firing
dev.off()

## calculate the firing frequency of neurons in cluster 1 or 2 in Pre
dat_cell_firng_pre_sta <- lapply(dat_cell_trace, function(x)  x[[1]]) %>% 
  lapply(., function(x) x[[1]]) %>% 
  do.call(cbind,.) %>% 
  add_column(Group = "Pre") %>% 
  add_column(Time = comp_time) %>% 
  pivot_longer(-c(Time, Group)) %>% 
  mutate(cluster = map_dbl(name, ~dat_cell_trace_cluster[.x])) %>% 
  ddply(., .(name, cluster), summarise, mean_acti = mean(value)) %>% 
  ddply(.,.(cluster), summarise, mean= mean(mean_acti), sd= sd(mean_acti))
  
## heat map plot the firing frequency

dat_cell_spike_combine %>% 
  ddply(.,.(Group, name), summarise, mean = mean(value)) %>% 
  mutate(cluster = map_dbl(name, ~dat_cell_trace_cluster[.x])) %>% 
  mutate(Group = factor(Group, levels = c("Pre", "Cond", "Post"))) %>% 
  filter(cluster ==2) %>% 
  ddply(.,.(Group), summarise, mean_acti = mean(mean))

score_range <- range(dat_cell_spike_combine$value)

for(i in c(1:3)){
  cell_order <- dat_cell_spike_combine %>% 
    filter(Group == "Cond") %>% 
    mutate(cluster = map_dbl(name, ~dat_cell_trace_cluster[.x])) %>% 
    mutate(cluster = factor(cluster, levels = c(1,2))) %>% 
    group_by(name, cluster) %>%
    summarise(mean = mean(value), .groups = "drop") %>%
    arrange(cluster, mean)
  
  dat_rect <- dat_cell_spike_combine %>% 
    filter(Group == "Cond") %>% 
    mutate(cluster = map_dbl(name, ~dat_cell_trace_cluster[.x])) %>% 
    mutate(cluster = factor(cluster, levels = c(2,1))) %>% 
    ddply(., .(cluster), summarise, n = length(name)/40) %>% 
    add_column(xmin = 2, xmax = 2.1) %>% 
    mutate(ymax = cumsum(n) ) %>% 
    mutate(ymin = c(1, lag(ymax)[-1]))
  
  p_heat <- dat_cell_spike_combine %>% 
    filter(Group == group_day[[i]]) %>% 
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
cairo_pdf("p_cell_firing.pdf", width = 50/25.6, height = 60/25.6, family = "Arial" )
p_cell_firing
dev.off()

t_cell_firing <- dat_cell_spike_combine %>% 
  ddply(.,.(Group, name), summarise, mean = mean(value)) %>% 
  mutate(cluster = map_dbl(name, ~dat_cell_trace_cluster[.x])) %>% 
  filter(cluster==1) %>% 
  filter(Group != "Post") %>% 
  wilcox.test(mean~Group,., paired = T)

summary(t_cell_firing)
TukeyHSD(t_cell_firing)


## correlation analysis of PC neurons to identify clusters------
mouse_file <- as.list(list.files("~cchen/Documents/neuroscience/Pn\ project/Data_analysis/For\ revision/PC_data", full.names = T))
dat_trace1 <- raveio::read_mat(mouse_file[[6]])

cross_ID <- dat_trace1$global.map %>% 
  as_tibble() %>% 
  select(V1, V4, V5) %>% 
  mutate_all(na_if, 0) %>% 
  drop_na() %>% 
  pull(V5)

dat_trace_d1_res <- dat_trace1$traces[[5]] %>%
  as_tibble() %>% 
  select(all_of(cross_ID)) %>% 
  apply(., 2, scale) %>% 
  cor() %>% 
  round(.,2)

#library(corrplot)
corrplot(dat_trace_d1_res, method = "square",  order = "hclust", tl.col = "black", tl.srt = 45, tl.pos = "n")

# Perform hierarchical clustering
hc <- hclust(as.dist(1 - cor_matrix))  # Using 1 - correlation as distance

# Plot the dendrogram (optional)
plot(hc, labels = colnames(dat_trace_d1_res), main = "Dendrogram of Variables")


k <- 4 # Example: aiming for 5 clusters
clusters <- cutree(hc, k = k)

# Associate each variable with its cluster
variable_clusters <- data.frame(Variable = names(clusters), Cluster = clusters) %>% 
  filter(Cluster ==4) %>% 
  mutate(name = as.numeric(str_extract(Variable,"[0-9]+"))) %>% 
  pull(name)

# View variables by cluster
print(variable_clusters)

## For a mannualy sorted single videos-----

path_file <- "~cchen/Documents/neuroscience/Pn\ project/Data_analysis/For\ revision/m687_pc_single.mat"
dat_trace1 <- raveio::read_mat(path_file)[[2]]

cell_rm <- raveio::read_mat(path_file)[[1]] %>% 
  t() %>% 
  as_tibble() %>% 
  pull(V1)



dat_trace_d1_res <- dat_trace1 %>% 
  as_tibble() %>% 
  select(-all_of(cell_rm)) %>% 
  apply(., 2, scale) %>% 
  cor() %>% 
  round(.,2)

#library(corrplot)
corrplot(dat_trace_d1_res, method = "square",  order = "hclust", tl.col = "black", tl.srt = 45, tl.pos = "n")







