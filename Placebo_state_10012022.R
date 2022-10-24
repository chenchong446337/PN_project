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

dat_anti_sta <- 

dat_cell_acti <- lapply(dat_acti_combine, function(x) x[[1]]) %>% 
  do.call(rbind,.) %>% 
  mutate(Group = factor(Group, levels = c("Pre", "Cond", "Post")))

pairwise.t.test(dat_cell_acti$cell_acti, dat_cell_acti$Group, paired = T, alternative = "greater")

setwd("~cchen/Documents/neuroscience/Pn\ project/Figure/PDF/")
cairo_pdf("p_dat_anti.pdf", width = 50/25.6, height = 60/25.6, family = "Arial")
p_dat_anti
dev.off()

## compare all activity by all cells
dat_cell_pre_trace <- lapply(dat_cell_trace, function(x) x[[1]]) %>% 
  do.call(cbind,.) %>% 
  add_column(Group = "Pre")

dat_cell_cond_trace <- lapply(dat_cell_trace, function(x) x[[2]]) %>% 
  do.call(cbind,.) %>% 
  add_column(Group = "Cond")

p_cell_acti <- lapply(dat_cell_trace, function(x) x[[3]]) %>% 
  do.call(cbind,.) %>% 
  add_column(Group = "Post") %>% 
  rbind(dat_cell_pre_trace, dat_cell_cond_trace,.) %>%
  mutate(Group = factor(Group, levels = c("Pre", "Cond", "Post"))) %>% 
  pivot_longer(-Group) %>% 
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
  scale_y_continuous(limits = c(-1, 6))+
  theme(legend.title = element_blank(), legend.position = "none")

p_acti_combine <- plot_grid( p_cell_acti,p_dat_anti, nrow = 1)
setwd("~cchen/Documents/neuroscience/Pn\ project/Figure/PDF/")
cairo_pdf("p_acti_combine.pdf", width = 80/25.6, height = 60/25.6, family = "Arial")
p_acti_combine
dev.off()

t_cell_acti <- lapply(dat_cell_trace, function(x) x[[3]]) %>% 
  do.call(cbind,.) %>% 
  add_column(Group = "Post") %>% 
  rbind(dat_cell_pre_trace, dat_cell_cond_trace,.) %>%
  mutate(Group = factor(Group, levels = c("Pre", "Cond", "Post"))) %>% 
  pivot_longer(-Group) %>% 
  ddply(., .(name, Group), summarise, mean_acti = mean(value)) %>% 
  aov(mean_acti~Group,.)

summary(t_cell_acti)
TukeyHSD(t_cell_acti, ordered = T)

pairwise.t.test(t_cell_acti$mean_acti, t_cell_acti$Group, paired = T, alternative = "greater")


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
                              Last_frame= c(1756, 2777, 2536, 2238, 3062, 2306)) ## m3, last one 3942

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
  labs(x="", y="Population activity (∆F/F)")+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  scale_y_continuous(limits = c(-1, 6))+
  theme(legend.title = element_blank(), legend.position = "none")

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
TukeyHSD(t_cell_firing)


p_d7_combin <- plot_grid(p_cell_firing, p_dat_anti_cross, nrow = 1)
setwd("~cchen/Documents/neuroscience/Pn\ project/Figure/PDF/")
cairo_pdf("p_dat_firing_d7.pdf", width = 80/25.6, height = 60/25.6, family = "Arial" )
p_d7_combin
dev.off()

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
      mutate(d2 = d^2) %>% 
      dplyr::select(name, d2) %>% 
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

## for d^2-----
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
      mutate(d2 = d^2) %>% 
      dplyr::select(name, d2) %>% 
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
  ddply(., .(ID, Group), summarise,n=length(d2),mean=mean(d2),sd=sd(d2),se=sd(d2)/sqrt(length(d2))) %>% 
  mutate(Group = factor(Group, levels = c("Pre", "Post")))


p_d_combine <- ggplot(dat_cell_d_combine, aes(Group, mean, group = Group))+
  geom_boxplot(outlier.shape = NA)+
  geom_line(aes(group=ID), colour="gray90")+
  geom_jitter(aes(colour = Group, shape = Group),width = 0.2,  size=2)+
  labs(x="", y="(d')^2")+
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
  ggplot(., aes(Group, d2, group = Group))+
  geom_violin()+
  geom_jitter(aes(colour = Group, shape = Group),width = 0.2,  size=2, alpha = 0.3)+
  labs(x="", y="(d')^2")+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  # scale_y_continuous(limits = c(-0.5, 0.6))+
  theme(legend.title = element_blank(), legend.position = "none")

t_combine_cells <- lapply(dat_cell_trace, function(x)  x[[2]]) %>% 
  do.call(rbind,.) %>% 
  mutate(Group = "Post") %>% 
  rbind(dat_cell_d_pre, .) %>% 
  mutate(Group = factor(Group, levels = c("Pre", "Post"))) %>% 
  wilcox.test(d2~Group,.)
  

p_d_all <- plot_grid(p_d_combine, p_d_combine_cell, nrow = 1)

setwd("~cchen/Documents/neuroscience/Pn\ project/Figure/PDF/")
cairo_pdf("p_d_all.pdf", width = 90/25.6, height = 60/25.6, family = "Arial")
p_d_all
dev.off()

## test for individual mouse
t_d_mice <- lapply(dat_cell_trace, function(x)  x[[2]]) %>% 
  do.call(rbind,.) %>% 
  mutate(Group = "Post") %>% 
  rbind(dat_cell_d_pre, .) %>% 
  mutate(Group = factor(Group, levels = c("Pre", "Post"))) %>% 
  filter(ID == "m857") %>% 
  wilcox.test(d2~Group,.)

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
      mutate(d2 = d^2) %>% 
      dplyr::select(name, d2) %>% 
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
  ddply(., .(ID, Group), summarise,n=length(d2),mean=mean(d2),sd=sd(d2),se=sd(d2)/sqrt(length(d2))) %>% 
  mutate(Group = factor(Group, levels = c("Pre", "Post")))


p_d_combine <- ggplot(dat_cell_d_combine, aes(Group, mean, group = Group))+
  geom_boxplot(outlier.shape = NA)+
  geom_line(aes(group=ID), colour="gray90")+
  geom_jitter(aes(colour = Group, shape = Group),width = 0.2,  size=2)+
  scale_color_manual(values=c("#8491B4FF", "#3C5488FF"))+
  labs(x="", y="(d')^2")+
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
  ggplot(., aes(Group, d2, group = Group))+
  geom_violin( )+
  geom_jitter(aes(colour = Group, shape = Group),width = 0.2,  size=2, alpha = 0.3 )+
  scale_color_manual(values=c("#8491B4FF", "#3C5488FF"))+
  labs(x="", y="(d')^2" )+
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
  wilcox.test(d2~Group, ., paired = T, alternative = "less")

p_d_align <- plot_grid(p_d_combine, p_d_combine_cell, nrow = 1)

setwd("~cchen/Documents/neuroscience/Pn\ project/Figure/PDF/")
cairo_pdf("p_d_align.pdf", width = 90/25.6, height = 60/25.6, family = "Arial")
p_d_align
dev.off()

## test for individual mouse
t_d_mice <- lapply(dat_cell_trace, function(x)  x[[2]]) %>% 
  do.call(rbind,.) %>% 
  mutate(Group = "Post") %>% 
  rbind(dat_cell_d_pre, .) %>% 
  mutate(Group = factor(Group, levels = c("Pre", "Post"))) %>% 
  filter(ID == "m3") %>% 
  wilcox.test(d2~Group, ., paired = T)


## for permutation test----
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
  
  cell_port <- rep(0, 2)
  
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
      as_tibble() 
    
    dat_stim_back_sta <- dat_stim_back %>% 
      as_tibble() 
    
    dat_stim_combine <- rbind(dat_stim_cross_sta, dat_stim_back_sta) 
    
    cc_d_fun <- function(cell_acti, cell_treat) {
      cell_mean_diff <- mean(cell_acti[cell_treat=="Cross"])-mean(cell_acti[cell_treat=="Cross_back"])
      sd_pool <-  sqrt(sd(cell_acti[cell_treat=="Cross"])^2 + sd(cell_acti[cell_treat=="Cross_back"])^2)/sqrt(2)
      sd_pool <- ifelse(sd_pool < 1e-06, 1e-06, sd_pool)
      real_cell_d <- (cell_mean_diff/sd_pool)^2
      return(real_cell_d)
      
    }
    
    cell_coding <- rep(0, ncol(dat_stim_combine))
    ## for each neuron
    for (j in c(1:ncol(dat_stim_combine))) {
      cell_acti <- unname(unlist(dat_stim_combine[,j]))
      cell_treat <- rep(c("Cross", "Cross_back"), each = length(comp_time))
      
      real_d <- cc_d_fun(cell_acti, cell_treat)
      
      cell_d_sh <- rep(0, 10000)
      
      for (k in c(1:10000)) {
        cell_treat_sh <- sample(cell_treat)
        cell_d_sh[k] <- cc_d_fun(cell_acti, cell_treat_sh)
        
      }
      
      big_value <- sum( cell_d_sh > real_d )
      
      cell_coding[j] <- ifelse(big_value > 100, "N", "Y")
      
    }
    
    cell_port[i] <- length(cell_coding[cell_coding=="Y"])/ncol(dat_stim_combine)
    
    
  }
  return(cell_port)
}



mouse_file <- as.list(list.files("~cchen/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/FD-processed", full.names = T))

dat_cell_trace <- mapply(c_miniscope_matlab_ft, mouse_file)

p_cell_coding <- dat_cell_trace %>% 
  t() %>% 
  as_tibble() %>% 
  rename(Pre = V1, Post = V2) %>% 
  add_column(ID = sort(c("m3", "m7", "m17", "m18", "m855", "m857"))) %>% 
  pivot_longer(-ID) %>% 
  mutate(value = value *100) %>% 
  mutate(name = factor(name, levels = c("Pre", "Post"))) %>% 
  ggplot(., aes(name, value, group= name))+
  geom_boxplot(outlier.shape = NA)+
  geom_line(aes(group=ID), colour="gray90")+
  geom_jitter(aes(colour = name, shape = name),width = 0.2,  size=2)+
  labs(x="", y="Proportion of coding cells (%)", title = "All cells, p = 0.03")+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  # scale_y_continuous(limits = c(-0.5, 0.6))+
  theme(legend.title = element_blank(), legend.position = "none")

t_cell_coding <- dat_cell_trace %>% 
  t() %>% 
  as_tibble() %>% 
  rename(Pre = V1, Post = V2) %>% 
  add_column(ID = sort(c("m3", "m7", "m17", "m18", "m855", "m857"))) %>% 
  pivot_longer(-ID) %>% 
  mutate(value = value *100) %>% 
  mutate(name = factor(name, levels = c("Pre", "Post"))) %>% 
  wilcox.test(value~name, ., paired = T, alternative = "less" )

setwd("~cchen/Documents/neuroscience/Pn\ project/Figure/PDF/")
cairo_pdf("p_cell_coding_all.pdf", width = 60/25.6, height = 60/25.6, family = "Arial")
p_cell_coding
dev.off()


## For coding cell with all aligned neurons-----

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
  
  cell_port <- rep(0, 2)
  
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
      as_tibble() 
    
    dat_stim_back_sta <- dat_stim_back %>% 
      as_tibble() 
    
    dat_stim_combine <- rbind(dat_stim_cross_sta, dat_stim_back_sta) 
    
    cc_d_fun <- function(cell_acti, cell_treat) {
      cell_mean_diff <- mean(cell_acti[cell_treat=="Cross"])-mean(cell_acti[cell_treat=="Cross_back"])
      sd_pool <-  sqrt(sd(cell_acti[cell_treat=="Cross"])^2 + sd(cell_acti[cell_treat=="Cross_back"])^2)/sqrt(2)
      sd_pool <- ifelse(sd_pool < 1e-06, 1e-06, sd_pool)
      real_cell_d <- (cell_mean_diff/sd_pool)^2
      return(real_cell_d)
      
    }
    
    cell_coding <- rep(0, ncol(dat_stim_combine))
    ## for each neuron
    for (j in c(1:ncol(dat_stim_combine))) {
      cell_acti <- unname(unlist(dat_stim_combine[,j]))
      cell_treat <- rep(c("Cross", "Cross_back"), each = length(comp_time))
      
      real_d <- cc_d_fun(cell_acti, cell_treat)
      
      cell_d_sh <- rep(0, 10000)
      
      for (k in c(1:10000)) {
        cell_treat_sh <- sample(cell_treat)
        cell_d_sh[k] <- cc_d_fun(cell_acti, cell_treat_sh)
        
      }
      
      big_value <- sum( cell_d_sh > real_d )
      
      cell_coding[j] <- ifelse(big_value > 100, "N", "Y")
      
    }
    
    cell_port[i] <- length(cell_coding[cell_coding=="Y"])/ncol(dat_stim_combine)
    
    
  }
  return(cell_port)
}



mouse_file <- as.list(list.files("~cchen/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/FD-processed", full.names = T))

dat_cell_trace <- mapply(c_miniscope_matlab_ft, mouse_file)

p_cell_coding <- dat_cell_trace %>% 
  t() %>% 
  as_tibble() %>% 
  rename(Pre = V1, Post = V2) %>% 
  add_column(ID = sort(c("m3", "m7", "m17", "m18", "m855", "m857"))) %>% 
  pivot_longer(-ID) %>% 
  mutate(value = value *100) %>% 
  mutate(name = factor(name, levels = c("Pre", "Post"))) %>% 
  ggplot(., aes(name, value, group= name))+
  geom_boxplot(outlier.shape = NA)+
  geom_line(aes(group=ID), colour="gray90")+
  geom_jitter(aes(colour = name, shape = name),width = 0.2,  size=2)+
  labs(x="", y="Proportion of coding cells (%)", title = "Aligned cells, p = 0.05")+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  # scale_y_continuous(limits = c(-0.5, 0.6))+
  theme(legend.title = element_blank(), legend.position = "none")

t_cell_coding <- dat_cell_trace %>% 
  t() %>% 
  as_tibble() %>% 
  rename(Pre = V1, Post = V2) %>% 
  add_column(ID = sort(c("m3", "m7", "m17", "m18", "m855", "m857"))) %>% 
  pivot_longer(-ID) %>% 
  mutate(value = value *100) %>% 
  mutate(name = factor(name, levels = c("Pre", "Post"))) %>% 
  wilcox.test(value~name, ., paired = T, alternative = "less" )

setwd("~cchen/Documents/neuroscience/Pn\ project/Figure/PDF/")
cairo_pdf("p_cell_coding_align.pdf", width = 60/25.6, height = 60/25.6, family = "Arial")
p_cell_coding
dev.off()


## with raw data------
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
      select(all_of(global_cell)) 
      # apply(., 2, function(x) x/max(x))
    
    ## number of rows to be binned
    n <- 2 # 0.05*2=0.1
    t1_p <- dat_trace1$crossing[t_crossing[i]]
    ## extract cell activity when they cross the border
    #dat_stim1 <- dat_trace[(t1_p-40):(t1_p+140-1),] 
    dat_stim1 <- dat_trace[(t1_p-40):(t1_p+40-1),] 
    
    dat_stim <- aggregate(dat_stim1,list(rep(1:(nrow(dat_stim1)%/%n+1),each=n,len=nrow(dat_stim1))),mean)[-1]

    
    
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
  # scale_y_continuous(limits = c(0.05, 0.25))+
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

## compare peak value before and after conditioning-----
c_miniscope_matlab_ft <- function(file_trace) {
  ## import and format the data
  dat_trace1 <- raveio::read_mat(file_trace)
  ID <- str_extract(file_trace, regex("m\\d+"))
  dat_trace1[[2]][1] <- ifelse(ID=="m857", 1694, dat_trace1[[2]][1])
  num_compare <- c(4, 7, 8)
  t_crossing <- c(1, 4, 5)
  
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
      apply(., 2, scale) %>% 
      apply(., 2, max)
    
    dat_stim_trace[[i]] <- dat_trace
  }
  return(dat_stim_trace)
}


## analyze for all mouse 
mouse_file <- as.list(list.files("~cchen/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/FD-processed", full.names = T))

dat_cell_trace <- mapply(c_miniscope_matlab_ft, mouse_file, SIMPLIFY = F)

dat_cell_peak_pre <- lapply(dat_cell_trace, function(x) x[[1]]) %>% 
  do.call(c,.) %>% 
  unname()

dat_cell_peak_cond <- lapply(dat_cell_trace, function(x) x[[2]]) %>% 
  do.call(c,.) %>% 
  unname()

dat_cell_peak_post <- lapply(dat_cell_trace, function(x) x[[3]]) %>% 
  do.call(c,.) %>% 
  unname()

dat_peak <- tibble(pre = dat_cell_peak_pre, cond = dat_cell_peak_cond,post = dat_cell_peak_post)
  
p_pre_con <- ggplot(dat_peak, aes(cond, pre))+
  geom_point()+
  labs(x="Peak value during Cond (s.d.)", y="Peak value during Pre (s.d.)")+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  scale_y_continuous(limits = c(2,12))+
  scale_x_continuous(limits = c(2, 12))+
  theme(legend.position = 'none')+
  geom_abline(intercept = 0, slope = 1, size = 0.5) 

setwd("~cchen/Documents/neuroscience/Pn\ project/Figure/PDF/")
cairo_pdf("p_pre_con.pdf", width = 65/25.6, height = 65/25.6, family = "Arial")
p_pre_con
dev.off()

p_pre_post <- ggplot(dat_peak, aes(post, pre))+
  geom_point()+
  labs(x="Peak value during Post (s.d.)", y="Peak value during Pre (s.d)")+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  scale_y_continuous(limits = c(2,14))+
  scale_x_continuous(limits = c(2, 14))+
  theme(legend.position = 'none')+
  geom_abline(intercept = 0, slope = 1, size = 0.5) 

setwd("~cchen/Documents/neuroscience/Pn\ project/Figure/PDF/")
cairo_pdf("p_pre_post.pdf", width = 65/25.6, height = 65/25.6, family = "Arial")
p_pre_post
dev.off()

p_cond_post <- ggplot(dat_peak, aes(post, cond))+
  geom_point()+
  labs(x="Peak value during Post (s.d.)", y="Peak value during Pre (s.d.)")+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  scale_y_continuous(limits = c(2,14))+
  scale_x_continuous(limits = c(2, 14))+
  theme(legend.position = 'none')+
  geom_abline(intercept = 0, slope = 1, size = 0.5) 

setwd("~cchen/Documents/neuroscience/Pn\ project/Figure/PDF/")
cairo_pdf("p_dat_anti.pdf", width = 50/25.6, height = 65/25.6, family = "Arial")
p_dat_anti
dev.off()
  

cor.test(dat_peak$pre, dat_peak$cond)
cor.test(dat_peak$pre, dat_peak$post)
cor.test(dat_peak$cond, dat_peak$post)


## raw data mean and sd------
c_miniscope_matlab_ft <- function(file_trace) {
  ## import and format the data
  dat_trace1 <- raveio::read_mat(file_trace)
  ID <- str_extract(file_trace, regex("m\\d+"))
  dat_trace1[[2]][1] <- ifelse(ID=="m857", 1694, dat_trace1[[2]][1])
  num_compare <- c(4, 7, 8)
  t_crossing <- c(1, 4, 5)
  
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
      #apply(., 2, scale) %>% 
      apply(., 2, sd)
    
    dat_stim_trace[[i]] <- dat_trace
  }
  return(dat_stim_trace)
}


## analyze for all mouse 
mouse_file <- as.list(list.files("~cchen/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/FD-processed", full.names = T))

dat_cell_trace <- mapply(c_miniscope_matlab_ft, mouse_file, SIMPLIFY = F)

dat_cell_peak_pre <- lapply(dat_cell_trace, function(x) x[[1]]) %>% 
  do.call(c,.) %>% 
  unname()

dat_cell_peak_cond <- lapply(dat_cell_trace, function(x) x[[2]]) %>% 
  do.call(c,.) %>% 
  unname()

dat_cell_peak_post <- lapply(dat_cell_trace, function(x) x[[3]]) %>% 
  do.call(c,.) %>% 
  unname()

dat_peak <- tibble(pre = dat_cell_peak_pre, cond = dat_cell_peak_cond,post = dat_cell_peak_post)

p_pre_con <- ggplot(dat_peak, aes(cond, pre))+
  geom_point()+
  labs(x="sd during Cond", y="sd during Pre")+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  scale_y_continuous(limits = c(0,0.25))+
  scale_x_continuous(limits = c(0, 0.25))+
  theme(legend.position = 'none')+
  geom_abline(intercept = 0, slope = 1, size = 0.5) 

setwd("~cchen/Documents/neuroscience/Pn\ project/Figure/PDF/")
cairo_pdf("p_pre_con.pdf", width = 65/25.6, height = 65/25.6, family = "Arial")
p_pre_con
dev.off()

p_pre_post <- ggplot(dat_peak, aes(post, pre))+
  geom_point()+
  labs(x="sd during Post", y="sd during Pre")+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  scale_y_continuous(limits = c(0,0.25))+
  scale_x_continuous(limits = c(0, 0.25))+
  theme(legend.position = 'none')+
  geom_abline(intercept = 0, slope = 1, size = 0.5) 

setwd("~cchen/Documents/neuroscience/Pn\ project/Figure/PDF/")
cairo_pdf("p_pre_post.pdf", width = 65/25.6, height = 65/25.6, family = "Arial")
p_pre_post
dev.off()

p_cond_post <- ggplot(dat_peak, aes(post, cond))+
  geom_point()+
  labs(x="Peak value during Post (s.d.)", y="Peak value during Pre (s.d.)")+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  scale_y_continuous(limits = c(0,0.18))+
  scale_x_continuous(limits = c(0, 0.18))+
  theme(legend.position = 'none')+
  geom_abline(intercept = 0, slope = 1, size = 0.5) 

setwd("~cchen/Documents/neuroscience/Pn\ project/Figure/PDF/")
cairo_pdf("p_dat_anti.pdf", width = 50/25.6, height = 65/25.6, family = "Arial")
p_dat_anti
dev.off()


cor.test(dat_peak$pre, dat_peak$cond)
cor.test(dat_peak$pre, dat_peak$post)
cor.test(dat_peak$cond, dat_peak$post)

## for d^2 on D7-----
back_crossing_frame <- tibble(ID = c("m3", "m7", "m17", "m18", "m855", "m857"), 
                              Frame_back = c(1003, 1632, 2113,1923,1572,1892), 
                              Frame_back_back = c(1209, 1915, 2536, 2239, 1784, 1153*2  ),
                              Last_frame= c(3271,4078, 2536, 2239, 1776*2, 1153*2),
                              Last_frame_back = c(3781, 4159, 2770, 2674, 1832*2, 1422*2))

c_miniscope_matlab_ft <- function(file_trace) {
  ## import and format the data
  dat_trace1 <- raveio::read_mat(file_trace)
  ID_mouse <- str_extract(file_trace, regex("m\\d+"))
  
  t_crossing_back <- back_crossing_frame %>% 
    filter(ID== ID_mouse )
  t_crossing <- t_crossing_back %>% 
    add_column(cross = dat_trace1[[2]][5], .before = "Frame_back" ) %>% 
    mutate(cross_back = Frame_back, .before = "Frame_back") %>% 
    select(-ID)
  
  dat_stim_trace <- vector(mode = "list", length = 3)
  
  
  for (i in seq_along(1:3)) {
    dat_trace <- dat_trace1[[8]] %>% 
      as_tibble() %>% 
      apply(., 2, scale)
    
    ## number of rows to be binned
    n <- 2 # 0.05*2=0.1
    t1_p <- unname(unlist(t_crossing[,i*2 -1]))
    t1_p_back <- unname(unlist(t_crossing[,i*2]))
    
    ## extract cell activity when they cross the border
    dat_stim1 <- dat_trace[(t1_p-40):(t1_p+40-1),] 
    
    dat_stim_cross <- aggregate(dat_stim1,list(rep(1:(nrow(dat_stim1)%/%n+1),each=n,len=nrow(dat_stim1))),mean)[-1]
    colnames(dat_stim_cross) <- str_c(ID_mouse,"Cell", 1: ncol(dat_stim_cross))
    
    ## extract cell activity when cross back
    dat_stim2 <- dat_trace[(t1_p_back-40):(t1_p_back+40-1),] 
    
    dat_stim_back <- aggregate(dat_stim2,list(rep(1:(nrow(dat_stim2)%/%n+1),each=n,len=nrow(dat_stim2))),mean)[-1]
    colnames(dat_stim_back) <- str_c(ID_mouse,"Cell", 1: ncol(dat_stim_back))
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
      mutate(ID = ID_mouse)
    
    dat_stim_trace[[i]] <- dat_stim_combine
    
    
  }
  return(dat_stim_trace)
}

mouse_file <- as.list(list.files("~cchen/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/FD-processed", full.names = T))

dat_cell_trace <- mapply(c_miniscope_matlab_ft, mouse_file, SIMPLIFY = F)

dat_cell_d_1st_back <- lapply(dat_cell_trace, function(x)  x[[1]]) %>% 
  do.call(rbind,.) %>% 
  mutate(Group = "1")

dat_cell_d_1st_last <- lapply(dat_cell_trace, function(x)  x[[2]]) %>% 
  do.call(rbind,.) %>% 
  mutate(Group = "2")

dat_cell_d_combine <- lapply(dat_cell_trace, function(x)  x[[3]]) %>% 
  do.call(rbind,.) %>% 
  mutate(Group = "3") %>% 
  rbind(dat_cell_d_1st_back, dat_cell_d_1st_last, .) %>% 
  ddply(., .(ID, Group), summarise,n=length(d2),mean=mean(d2),sd=sd(d2),se=sd(d2)/sqrt(length(d2))) %>% 
  mutate(Group = factor(Group, levels = c("1", "2", "3")))


p_d_combine <- dat_cell_d_combine %>% 
  #dplyr::filter(Group != "3") %>% 
  ggplot(., aes(Group, mean, group = Group))+
  geom_boxplot(outlier.shape = NA)+
  geom_line(aes(group=ID), colour="gray90")+
  geom_jitter(aes(colour = Group, shape = Group),width = 0.2,  size=2)+
  labs(x="", y="(d')^2")+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  scale_y_continuous(limits = c(0, 80))+
  theme(legend.title = element_blank(), legend.position = "none")



t_cell_d_combine <- lapply(dat_cell_trace, function(x)  x[[3]]) %>% 
  do.call(rbind,.) %>% 
  mutate(Group = "3") %>% 
  rbind(dat_cell_d_1st_back, dat_cell_d_1st_last, .) %>% 
  ddply(., .(ID, Group), summarise,n=length(d2),mean=mean(d2),sd=sd(d2),se=sd(d2)/sqrt(length(d2))) %>% 
  mutate(Group = factor(Group, levels = c("1", "2", "3"))) %>% 
  aov(mean~Group, .)

summary(t_cell_d_combine)
TukeyHSD(t_cell_d_combine)

dat_cell_d <- lapply(dat_cell_trace, function(x)  x[[3]]) %>% 
  do.call(rbind,.) %>% 
  mutate(Group = "3") %>% 
  rbind(dat_cell_d_1st_back, dat_cell_d_1st_last, .) %>% 
  ddply(., .(ID, Group), summarise,n=length(d2),mean=mean(d2),sd=sd(d2),se=sd(d2)/sqrt(length(d2))) %>% 
  mutate(Group = factor(Group, levels = c("1", "2", "3")))

pairwise.t.test(dat_cell_d$mean, dat_cell_d$Group, paired = T, alternative = "less")


## compare by cells
p_d_combine_cell <- lapply(dat_cell_trace, function(x)  x[[3]]) %>% 
  do.call(rbind,.) %>% 
  mutate(Group = "3") %>% 
  rbind(dat_cell_d_1st_back, dat_cell_d_1st_last, .) %>%  
  mutate(Group = factor(Group, levels = c("1", "2", "3"))) %>% 
  ggplot(., aes(Group, d2, group = Group))+
  geom_violin()+
  geom_jitter(aes(colour = Group, shape = Group),width = 0.2,  size=2, alpha = 0.3)+
  labs(x="", y="(d')^2")+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  scale_y_continuous(limits = c(-5, 600))+
  theme(legend.title = element_blank(), legend.position = "none")

t_combine_cells <- lapply(dat_cell_trace, function(x)  x[[3]]) %>% 
  do.call(rbind,.) %>% 
  mutate(Group = "3") %>% 
  rbind(dat_cell_d_1st_back, dat_cell_d_1st_last, .) %>%  
  mutate(Group = factor(Group, levels = c("1", "2", "3"))) %>% 
  aov(d2~Group,.)
summary(t_combine_cells)
TukeyHSD(t_combine_cells)


p_d_all <- plot_grid(p_d_combine, p_d_combine_cell, nrow = 1)

setwd("~cchen/Documents/neuroscience/Pn\ project/Figure/PDF/")
cairo_pdf("p_d_all.pdf", width = 100/25.6, height = 60/25.6, family = "Arial")
p_d_all
dev.off()

## test for individual mouse
t_d_mice <- lapply(dat_cell_trace, function(x)  x[[2]]) %>% 
  do.call(rbind,.) %>% 
  mutate(Group = "Post") %>% 
  rbind(dat_cell_d_pre, .) %>% 
  mutate(Group = factor(Group, levels = c("Pre", "Post"))) %>% 
  filter(ID == "m857") %>% 
  wilcox.test(d2~Group,.)

## barplot the number of cells-----
mouse_ID <- c("m3", "m7", "m17", "m18", "m855", "m857")
n_cells <- c(82, 62, 36, 45, 90, 42)
mean(n_cells)
se=sd(n_cells)/sqrt(length(n_cells))

p_number <- tibble(num = n_cells, Group = "Num") %>% 
  ggplot(., aes(x = Group, y = num))+
  geom_boxplot()+
  geom_jitter(width = 0.2, shape=1)+
  labs(x="", y="No. of cells")+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  scale_y_continuous(limits = c(0, 100))+
  theme(legend.position = "none")

setwd("~cchen/Documents/neuroscience/Pn\ project/Figure/PDF/")
cairo_pdf("p_number.pdf", width = 30/25.6, height = 60/25.6, family = "Arial")
p_number
dev.off()

## with new data-----
c_miniscope_matlab_ft <- function(file_trace) {
  ## import and format the data
  dat_trace1 <- raveio::read_mat(file_trace)
  ID <- str_extract(file_trace, regex("m\\d+"))
  while (length(dat_trace1)< 8) {
    dat_trace1 <-lapply(rapply(dat_trace1, enquote, how="unlist"), eval)
    ref <- list(c(0,0))
    dat_trace1 <- append(ref, dat_trace1)
    names(dat_trace1)[3] = "global_map"
    
  }
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
TukeyHSD(p_test)


dat_cell_acti <- lapply(dat_acti_combine, function(x) x[[1]]) %>% 
  do.call(rbind,.) %>% 
  mutate(Group = factor(Group, levels = c("Pre", "Cond", "Post")))

pairwise.t.test(dat_cell_acti$cell_acti, dat_cell_acti$Group, paired = T, alternative = "greater")

setwd("~cchen/Documents/neuroscience/Pn\ project/Figure/PDF/")
cairo_pdf("p_dat_anti.pdf", width = 50/25.6, height = 60/25.6, family = "Arial")
p_dat_anti
dev.off()

## compare all activity by all cells
dat_cell_pre_trace <- lapply(dat_cell_trace, function(x) x[[1]]) %>% 
  do.call(cbind,.) %>% 
  add_column(Group = "Pre")

dat_cell_cond_trace <- lapply(dat_cell_trace, function(x) x[[2]]) %>% 
  do.call(cbind,.) %>% 
  add_column(Group = "Cond")

p_cell_acti <- lapply(dat_cell_trace, function(x) x[[3]]) %>% 
  do.call(cbind,.) %>% 
  add_column(Group = "Post") %>% 
  rbind(dat_cell_pre_trace, dat_cell_cond_trace,.) %>%
  mutate(Group = factor(Group, levels = c("Pre", "Cond", "Post"))) %>% 
  pivot_longer(-Group) %>% 
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
  scale_y_continuous(limits = c(-1, 6))+
  theme(legend.title = element_blank(), legend.position = "none")

p_acti_combine <- plot_grid( p_cell_acti,p_dat_anti, nrow = 1)
setwd("~cchen/Documents/neuroscience/Pn\ project/Figure/PDF/")
cairo_pdf("p_acti_combine.pdf", width = 80/25.6, height = 60/25.6, family = "Arial")
p_acti_combine
dev.off()

t_cell_acti <- lapply(dat_cell_trace, function(x) x[[3]]) %>% 
  do.call(cbind,.) %>% 
  add_column(Group = "Post") %>% 
  rbind(dat_cell_pre_trace, dat_cell_cond_trace,.) %>%
  mutate(Group = factor(Group, levels = c("Pre", "Cond", "Post"))) %>% 
  pivot_longer(-Group) %>% 
  ddply(., .(name, Group), summarise, mean_acti = mean(value)) %>% 
  aov(mean_acti~Group,.)

summary(t_cell_acti)
TukeyHSD(t_cell_acti, ordered = T)

pairwise.t.test(t_cell_acti$mean_acti, t_cell_acti$Group, paired = T, alternative = "greater")

## plot the cell number across days-----
num_pre <- c(74, 55, 35, 45, 89, 42)
num_cond <- c(68, 53, 36, 41, 87, 29)
num_post <- c(82, 62, 27, 42, 82, 24)
num_aligned = c(48, 36, 19, 25, 58, 19)

dat_num <- tibble(Pre = num_pre, Cond = num_cond, Post = num_post, Align = num_aligned) %>% 
  add_column(ID = c("m3", "m7", "m17", "m18", "m855", "m857")) %>%
  pivot_longer(-ID) %>% 
  mutate(name = factor(name, levels = c("Pre", "Cond", "Post", "Align"))) %>% 
  ggplot(., aes(name, value, group = name))+
  geom_boxplot(outlier.shape = NA)+
  #geom_line(aes(group=ID), colour="gray90")+
  geom_jitter(aes(colour = name, shape = name),width = 0.2,  size=2)+
  scale_color_manual(values=c("#8491B4FF", "#00A087FF", "#3C5488FF", "violetred3"))+
  labs(x="", y="No. of cells")+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  scale_y_continuous(limits = c(10, 100))+
  theme(legend.title = element_blank(), legend.position = "none")

dat_num_sta <- tibble(Pre = num_pre, Cond = num_cond, Post = num_post, Align = num_aligned) %>% 
  add_column(ID = c("m3", "m7", "m17", "m18", "m855", "m857")) %>%
  pivot_longer(-ID) %>% 
  mutate(name = factor(name, levels = c("Pre", "Cond", "Post", "Align"))) %>% 
  ddply(., .(name), summarise,n=length(value),mean=mean(value),sd=sd(value),se=sd(value)/sqrt(length(value)))
  
setwd("~cchen/Documents/neuroscience/Pn\ project/Figure/PDF/")
cairo_pdf("p_number.pdf", width = 55/25.6, height = 65/25.6, family = "Arial")
dat_num
dev.off()

## corss-day aligned neuron number
p_number <- tibble(num = c(48, 36, 19, 25, 58, 19), Group = "Num") %>% 
  ggplot(., aes(x = Group, y = num))+
  geom_boxplot()+
  geom_jitter(width = 0.2, shape=1)+
  labs(x="", y="No. of cross-day aligned cells")+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  scale_y_continuous(limits = c(0, 100))+
  theme(legend.position = "none")

setwd("~cchen/Documents/neuroscience/Pn\ project/Figure/PDF/")
cairo_pdf("p_number_cross.pdf", width = 40/25.6, height = 65/25.6, family = "Arial")
p_number
dev.off()
