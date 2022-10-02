## crossed-day-aligned neurons, heatmap, venn plot------

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

## heatmap plot the activity, ordered by cell activity
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
cairo_pdf("p_heat.pdf", width = 140/25.6, height = 60/25.6, family = "Arial")
p_heat
dev.off()

## calculate cell activity change mouse by mouse (cell becomes more active)
cc_active_fun <- function (dat_trace){
  m_id <- dat_trace[[1]] %>% 
    colnames() %>%
    .[1] %>% 
    str_extract(., regex("m\\d+"))
  
  dat_trace_combine <- dat_trace %>% 
    do.call(rbind,.) %>% 
    as_tibble() %>% 
    add_column(Time = rep(comp_time, 3)) %>% 
    add_column(Group = rep(c("Pre", "Cond", "Post"), each = length(comp_time))) %>% 
    pivot_longer(-c(Time, Group)) 
  
  dat_trace_acti <- dat_trace_combine %>% 
    ddply(., .(name, Group), summarise, mean_acti = mean(value)) %>% 
    ddply(., .(Group), summarise, cell_acti = mean(mean_acti)) %>% 
    add_column(ID = m_id)
  
  # calculate percentage of cell numer increased activity
  ratio_pre_con <- dat_trace_combine %>% 
    ddply(., .(name, Group), summarise, mean_acti = mean(value)) %>% 
    pivot_wider(., names_from = Group, values_from = mean_acti) %>% 
    mutate(diff = Cond - Pre) %>% 
    dplyr::filter(diff > 0)
  
  ratio_pre_post <- dat_trace_combine %>% 
    ddply(., .(name, Group), summarise, mean_acti = mean(value)) %>% 
    pivot_wider(., names_from = Group, values_from = mean_acti) %>% 
    mutate(diff = Post - Pre) %>% 
    dplyr::filter(diff > 0)
  
  ratio_con_post <- dat_trace_combine %>% 
    ddply(., .(name, Group), summarise, mean_acti = mean(value)) %>% 
    pivot_wider(., names_from = Group, values_from = mean_acti) %>% 
    mutate(diff = Post - Cond) %>% 
    dplyr::filter(diff > 0)
  
  
  dat_ratio <- c(nrow(ratio_pre_con), nrow(ratio_pre_post), nrow(ratio_con_post))/ncol(dat_trace[[1]])
  
  name_con_pre <- ratio_pre_con$name
  name_post_pre <- ratio_pre_post$name
  name_post_con <- ratio_con_post$name
  return(list(dat_trace_acti, dat_ratio, name_con_pre, name_post_pre, name_post_con))
}


dat_acti_combine <- mapply(cc_active_fun, dat_cell_trace, SIMPLIFY = F)

p_dat_anti <- lapply(dat_acti_combine, function(x) x[[1]]) %>% 
  do.call(rbind,.) %>% 
  mutate(Group = factor(Group, levels = c("Pre", "Cond", "Post"))) %>% 
  ggplot(., aes(Group, cell_acti, group = Group))+
  geom_boxplot(outlier.shape = NA)+
  geom_line(aes(group=ID), colour="gray90")+
  geom_jitter(aes(colour = Group, shape = Group),width = 0.2,  size=2)+
  scale_color_manual(values=c("#8491B4FF", "#00A087FF", "#3C5488FF"))+
  labs(x="", y="Population activity (∆F/F)")+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  scale_y_continuous(limits = c(-0.3, 0.6))+
  theme(legend.title = element_blank(), legend.position = "none")

p_test <- lapply(dat_acti_combine, function(x) x[[1]]) %>% 
  do.call(rbind,.) %>% 
  mutate(Group = factor(Group, levels = c("Pre", "Cond", "Post"))) %>% 
  aov(cell_acti~ Group,.)

summary(p_test)



dat_cell_acti <- lapply(dat_acti_combine, function(x) x[[1]]) %>% 
  do.call(rbind,.) %>% 
  mutate(Group = factor(Group, levels = c("Pre", "Cond", "Post")))

pairwise.t.test(dat_cell_acti$cell_acti, dat_cell_acti$Group, paired = T, alternative = "greater")

setwd("~cchen/Documents/neuroscience/Pn\ project/Figure/PDF/")
cairo_pdf("p_dat_anti.pdf", width = 50/25.6, height = 65/25.6, family = "Arial")
p_dat_anti
dev.off()

## compare cells show increased activity
p_cell_portion <- lapply(dat_acti_combine, function(x) x[[2]]) %>% 
  do.call(rbind,.) %>% 
  as_tibble() %>% 
  rename("Cond>Pre" = V1, "Post>Pre" = V2, "Post>Cond" = V3) %>% 
  mutate(ID = sort(c("m3", "m7", "m17", "m18", "m855", "m857"))) %>% 
  pivot_longer(-ID) %>% 
  mutate(value = value*100) %>% 
  mutate(name = factor(name, levels = c("Cond>Pre",  "Post>Cond", "Post>Pre"))) %>% 
  ggplot(., aes(name, value, color = name))+
  geom_boxplot(outlier.shape = NA)+
  #geom_line(aes(group=ID), colour="gray90")+
  geom_jitter(aes(colour = name, shape = name),width = 0.2,  size=2)+
  scale_color_manual(values=c("#8491B4FF", "#00A087FF", "#3C5488FF"))+
  labs(x="", y="Proportion of neurons (%)")+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  scale_y_continuous(limits = c(20, 100))+
  geom_hline(yintercept = 50, linetype = 2, col = "red")+
  theme(legend.title = element_blank(), legend.position = "none")

setwd("~cchen/Documents/neuroscience/Pn\ project/Figure/PDF/")
cairo_pdf("p_cell_portion.pdf", width = 60/25.6, height = 65/25.6, family = "Arial")
p_cell_portion
dev.off()

cell_portion_sta <- lapply(dat_acti_combine, function(x) x[[2]]) %>% 
  do.call(rbind,.) %>% 
  as_tibble() %>% 
  rename("Cond>Pre" = V1, "Post>Pre" = V2, "Post>Cond" = V3) %>% 
  mutate(ID = sort(c("m3", "m7", "m17", "m18", "m855", "m857"))) %>% 
  pivot_longer(-ID) %>% 
  mutate(value = value*100) %>% 
  mutate(name = factor(name, levels = c("Cond>Pre",  "Post>Cond", "Post>Pre"))) %>% 
  ddply(., .(name), summarise,n=length(value),mean=mean(value),sd=sd(value),se=sd(value)/sqrt(length(value)))

## venn plot of the cells
dat_name_con_pre <- lapply(dat_acti_combine, function(x) x[[3]]) %>% 
  unlist()

dat_name_post_pre <- lapply(dat_acti_combine, function(x) x[[4]]) %>% 
  unlist()

dat_name_post_con <- lapply(dat_acti_combine, function(x) x[[5]]) %>% 
  unlist()
dat_overlap <- list(Pre = dat_name_con_pre, Cond= dat_name_post_pre, Post = dat_name_post_con)

ggvenn::ggvenn(dat_overlap,  c("Pre", "Cond", "Post"), fill_color = c("#8491B4FF", "#00A087FF", "blue"), stroke_color = "black")

## compare crossing, back-crossing and last crossing---------
back_crossing_frame <- tibble(ID = c("m3", "m7", "m17", "m18", "m855", "m857"), Frame = c(1003, 1632, 2113,1923,1572,1892), 
                              Last_frame= c(1756, 2777, 2536, 2238, 3062, 2306))

c_miniscope_matlab_ft <- function(file_trace) {
  ## import and format the data
  dat_trace1 <- raveio::read_mat(file_trace)
  ID_mouse <- str_extract(file_trace, regex("m\\d+"))
  
  t_crossing_back <- back_crossing_frame %>% 
    filter(ID== ID_mouse )
  t_crossing <- c(dat_trace1[[2]][5], t_crossing_back$Frame, t_crossing_back$Last_frame)
  
  
  
  dat_stim_trace <- vector(mode = "list", length = length(t_crossing))
  
  
  for (i in seq_along(t_crossing)) {
    dat_trace <- dat_trace1[[8]] %>% 
      as_tibble() %>% 
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
comp_time <- seq(-2, 1.9, by=0.1)
group_day <- c("Crossing", "Crossing_back", "Last_crossing")
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

p_heat <- plot_grid(p_heat_Crossing, p_heat_Crossing_back, p_heat_Last_crossing, nrow  = 1)


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
  scale_y_continuous(limits = c(-0.5, 0.6))+
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

dat_cell_acti <- lapply(dat_acti_combine, function(x) x[[1]]) %>% 
  do.call(rbind,.) %>% 
  mutate(Group = factor(Group, levels = c("Crossing", "Crossing_back", "Last_crossing")))

pairwise.t.test(dat_cell_acti$cell_acti, dat_cell_acti$Group, paired = T, alternative = "less")


## with Fatih's suggestion to calculate the signed d------
back_crossing_frame <- tibble(ID_mouse = c("m3", "m7", "m17", "m18", "m855", "m857"), Frame_d3= c(2292, 498, 610, 535,1044,1770),
                              Frame_d7 = c(1003, 1632, 2113,1923,1572,1892))
c_miniscope_matlab_ft <- function(file_trace) {
  ## import and format the data
  dat_trace1 <- raveio::read_mat(file_trace)
  ID <- str_extract(file_trace, regex("m\\d+"))
  dat_trace1[[2]][1] <- ifelse(ID=="m857", 1694, dat_trace1[[2]][1])
  num_compare <- c(4, 8)
  t_crossing <- c(1, 5)
  t_crossing_back <- back_crossing_frame %>% 
    dplyr::filter(ID_mouse == ID) %>% 
    dplyr::select(Frame_d3, Frame_d7) %>% 
    unlist() %>% 
    unname()
  
  
  dat_stim_trace <- vector(mode = "list", length = length(num_compare))
  
  
  for (i in seq_along(num_compare)) {
    
    dat_trace <- dat_trace1[[num_compare[i]]] %>% 
      as_tibble() %>% 
      apply(., 2, scale)
    
    ## number of rows to be binned
    n <- 2 # 0.05*2=0.1
    t1_p <- dat_trace1$crossing[t_crossing[i]]
    t1_p_back <- t_crossing_back[i]
    
    ## extract cell activity when they cross the border
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
    
    dat_stim_combine <- full_join(dat_stim_cross_sta, dat_stim_back_sta, by = 'name') %>% 
      mutate(sd_pool = sqrt(sd_cross^2 + sd_back^2)/sqrt(2)) %>% 
      mutate(sd_pool = ifelse(sd_pool < 1e-06, 1e-06, sd_pool)) %>% 
      mutate(d = (mean_cross - mean_back)/sd_pool) %>% 
      dplyr::select(name, d) %>% 
      mutate(ID = ID)
    
    dat_stim_trace[[i]] <- dat_stim_combine
    
    
  }
  return(dat_stim_trace)
}

mouse_file <- as.list(list.files("~cchen/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/FD-processed", full.names = T))

dat_cell_trace <- mapply(c_miniscope_matlab_ft, mouse_file, SIMPLIFY = F)

dat_cell_d_pre <- lapply(dat_cell_trace, function(x)  x[[1]]) %>% 
  do.call(rbind,.) %>% 
  mutate(Group = "Pre")

dat_cell_d_combine <- lapply(dat_cell_trace, function(x)  x[[2]]) %>% 
  do.call(rbind,.) %>% 
  mutate(Group = "Post") %>% 
  rbind(dat_cell_d_pre, .) %>% 
  ddply(., .(ID, Group), summarise,n=length(d),mean=mean(d),sd=sd(d),se=sd(d)/sqrt(length(d))) %>% 
  mutate(Group = factor(Group, levels = c("Pre", "Post")))


p_d_combine <- ggplot(dat_cell_d_combine, aes(Group, mean, group = Group))+
  geom_boxplot(outlier.shape = NA)+
  geom_line(aes(group=ID), colour="gray90")+
  geom_jitter(aes(colour = Group, shape = Group),width = 0.2,  size=2)+
  labs(x="", y="d'", title = "All cells, p = 0.015")+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  # scale_y_continuous(limits = c(-0.5, 0.6))+
  theme(legend.title = element_blank(), legend.position = "none")

wilcox.test(mean~Group, dat_cell_d_combine, paired = T, alternative = "less")

## compare by cells
p_d_combine_cell <- lapply(dat_cell_trace, function(x)  x[[2]]) %>% 
  do.call(rbind,.) %>% 
  mutate(Group = "Post") %>% 
  rbind(dat_cell_d_pre, .) %>% 
  mutate(Group = factor(Group, levels = c("Pre", "Post"))) %>% 
  ggplot(., aes(Group, d, group = Group))+
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(aes(colour = Group, shape = Group),width = 0.2,  size=2)+
  labs(x="", y="d'", title = "All cells")+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  # scale_y_continuous(limits = c(-0.5, 0.6))+
  theme(legend.title = element_blank(), legend.position = "none")

p_d_all <- plot_grid(p_d_combine, p_d_combine_cell, nrow = 1)

setwd("~cchen/Documents/neuroscience/Pn\ project/Figure/PDF/")
cairo_pdf("p_d_all.pdf", width = 120/25.6, height = 60/25.6, family = "Arial")
p_d_all
dev.off()

## test for individual mouse
t_d_mice <- lapply(dat_cell_trace, function(x)  x[[2]]) %>% 
  do.call(rbind,.) %>% 
  mutate(Group = "Post") %>% 
  rbind(dat_cell_d_pre, .) %>% 
  mutate(Group = factor(Group, levels = c("Pre", "Post"))) %>% 
  filter(ID == "m855") %>% 
  wilcox.test(d~Group,.)

## cross-day aligned neurons-----
back_crossing_frame <- tibble(ID_mouse = c("m3", "m7", "m17", "m18", "m855", "m857"), Frame_d3= c(2292, 498, 610, 535,1044,1770),
                              Frame_d7 = c(1003, 1632, 2113,1923,1572,1892))
c_miniscope_matlab_ft <- function(file_trace) {
  ## import and format the data
  dat_trace1 <- raveio::read_mat(file_trace)
  ID <- str_extract(file_trace, regex("m\\d+"))
  dat_trace1[[2]][1] <- ifelse(ID=="m857", 1694, dat_trace1[[2]][1])
  num_compare <- c(4, 8)
  t_crossing <- c(1, 5)
  t_crossing_back <- back_crossing_frame %>% 
    dplyr::filter(ID_mouse == ID) %>% 
    dplyr::select(Frame_d3, Frame_d7) %>% 
    unlist() %>% 
    unname()
  
  cross_ID <- dat_trace1$global_map %>% 
    as_tibble() %>% 
    select(V1, V5) %>% 
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
    t1_p_back <- t_crossing_back[i]
    
    ## extract cell activity when they cross the border
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
    
    dat_stim_combine <- full_join(dat_stim_cross_sta, dat_stim_back_sta, by = 'name') %>% 
      mutate(sd_pool = sqrt(sd_cross^2 + sd_back^2)/sqrt(2)) %>% 
      mutate(sd_pool = ifelse(sd_pool < 1e-06, 1e-06, sd_pool)) %>% 
      mutate(d = (mean_cross - mean_back)/sd_pool) %>% 
      dplyr::select(name, d) %>% 
      mutate(ID = ID)
    
    dat_stim_trace[[i]] <- dat_stim_combine
    
    
  }
  return(dat_stim_trace)
}

mouse_file <- as.list(list.files("~cchen/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/FD-processed", full.names = T))

dat_cell_trace <- mapply(c_miniscope_matlab_ft, mouse_file, SIMPLIFY = F)

dat_cell_d_pre <- lapply(dat_cell_trace, function(x)  x[[1]]) %>% 
  do.call(rbind,.) %>% 
  mutate(Group = "Pre")

dat_cell_d_combine <- lapply(dat_cell_trace, function(x)  x[[2]]) %>% 
  do.call(rbind,.) %>% 
  mutate(Group = "Post") %>% 
  rbind(dat_cell_d_pre, .) %>% 
  ddply(., .(ID, Group), summarise,n=length(d),mean=mean(d),sd=sd(d),se=sd(d)/sqrt(length(d))) %>% 
  mutate(Group = factor(Group, levels = c("Pre", "Post")))


p_d_combine <- ggplot(dat_cell_d_combine, aes(Group, mean, group = Group))+
  geom_boxplot(outlier.shape = NA)+
  geom_line(aes(group=ID), colour="gray90")+
  geom_jitter(aes(colour = Group, shape = Group),width = 0.2,  size=2)+
  labs(x="", y="d'", title = "Aligned cells, p = 0.015")+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  # scale_y_continuous(limits = c(-0.5, 0.6))+
  theme(legend.title = element_blank(), legend.position = "none")

wilcox.test(mean~Group, dat_cell_d_combine, paired = T, alternative = "less")

## compare by cells
p_d_combine_cell <- lapply(dat_cell_trace, function(x)  x[[2]]) %>% 
  do.call(rbind,.) %>% 
  mutate(Group = "Post") %>% 
  rbind(dat_cell_d_pre, .) %>% 
  mutate(Group = factor(Group, levels = c("Pre", "Post"))) %>% 
  ggplot(., aes(Group, d, group = Group))+
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(aes(colour = Group, shape = Group),width = 0.2,  size=2)+
  labs(x="", y="d'",title = "Aligned cells, p= 0.059" )+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  # scale_y_continuous(limits = c(-0.5, 0.6))+
  theme(legend.title = element_blank(), legend.position = "none")

t_d_combine_cell <- lapply(dat_cell_trace, function(x)  x[[2]]) %>% 
  do.call(rbind,.) %>% 
  mutate(Group = "Post") %>% 
  rbind(dat_cell_d_pre, .) %>% 
  mutate(Group = factor(Group, levels = c("Pre", "Post"))) %>% 
  wilcox.test(d~Group, ., paired = T, alternative = "less")

p_d_align <- plot_grid(p_d_combine, p_d_combine_cell, nrow = 1)

setwd("~cchen/Documents/neuroscience/Pn\ project/Figure/PDF/")
cairo_pdf("p_d_align.pdf", width = 120/25.6, height = 60/25.6, family = "Arial")
p_d_align
dev.off()


