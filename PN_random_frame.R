## randomize the frame of crossing to exclude overfitting issue

library(svMisc)

dat_cell_trace_pre_random <- vector(mode = "list", 100)
dat_cell_trace_test_random <- vector(mode = "list", 100)
random_trails <- c(1:100)
for ( i in seq_along(random_trails)){
  
  progress(i, progress.bar = TRUE)
  Sys.sleep (0.01)
  
  frame_random <- sample(c(50:3600), 12)
  
  path_trace_m3 <- as.list(list.files(path ="~cchen/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/matlab_analysis/m3/trace_days/", pattern = "*.xlsx", full.names = T ))[c(3,7)]
  t_stim_m3 <- list(t_stim_m3_d3 = frame_random[1], t_stim_m3_d7 = frame_random[2])
  dat_trace_m3 <- mapply(c_miniscope_matlab, path_trace_m3, t_stim_m3, SIMPLIFY = F)
  
  # for m7
  path_trace_m7 <- as.list(list.files(path ="~cchen/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/matlab_analysis/m7/trace_days/", pattern = "*.xlsx", full.names = T ))[c(3,7)]
  t_stim_m7 <- list(t_stim_m7_d3 = frame_random[3], t_stim_m7_d7 = frame_random[4])
  dat_trace_m7 <- mapply(c_miniscope_matlab, path_trace_m7, t_stim_m7,SIMPLIFY = F)
  
  # for m17
  path_trace_m17 <- as.list(list.files(path ="~cchen/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/matlab_analysis/m17/trace_days/", pattern = "*.xlsx", full.names = T ))[c(3,7)]
  t_stim_m17 <- list(t_stim_m17_d3 = frame_random[5], t_stim_m17_d7 = sample(frame_random[frame_random<2900], 1) )
  dat_trace_m17 <- mapply(c_miniscope_matlab, path_trace_m17, t_stim_m17, SIMPLIFY = F)
  
  # for m18
  path_trace_m18 <- as.list(list.files(path ="~cchen/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/matlab_analysis/m18/trace_days/", pattern = "*.xlsx", full.names = T ))[c(3,7)]
  t_stim_m18 <- list(t_stim_m18_d3 = frame_random[7], t_stim_m18_d7 = sample(frame_random[frame_random<2800], 1))
  dat_trace_m18 <- mapply(c_miniscope_matlab, path_trace_m18, t_stim_m18, SIMPLIFY = F)
  
  # for m855
  path_trace_m855 <- as.list(list.files(path ="~cchen/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/matlab_analysis/m855/trace_days", pattern = "*.xlsx", full.names = T ))[c(3,7)]
  t_stim_m855 <- list(t_stim_m855_d3 = frame_random[9], t_stim_m855_d7 = frame_random[10])
  dat_trace_m855 <- mapply(c_miniscope_matlab, path_trace_m855, t_stim_m855, SIMPLIFY = F)
  
  # for m857
  path_trace_m857 <- as.list(list.files(path ="~cchen/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/matlab_analysis/m857/trace_days", pattern = "*.xlsx", full.names = T ))[c(1, 4)]
  t_stim_m857 <- list(t_stim_m857_d3 = frame_random[11], t_stim_m857_d7 = frame_random[12])
  dat_trace_m857 <- mapply(c_miniscope_matlab, path_trace_m857, t_stim_m857, SIMPLIFY = F)
  
  mouse_ID <- c("m3", "m7", "m17", "m18", "m855", "m857")
  
  
  dat_cell_trace <- vector(mode = "list", 2)
  for (j in 1:2) {
    dat_cell_trace[[j]] <- cbind(dat_trace_m3[[j]], dat_trace_m7[[j]], dat_trace_m17[[j]], dat_trace_m18[[j]], dat_trace_m855[[j]], dat_trace_m857[[j]])
    
  }
  
  dat_cell_trace_re <- vector(mode = "list", 2)
  
  for (k in 1:length(dat_cell_trace)) {
    dat_cell_trace_d <- dat_cell_trace[[k]]
    cell_id <- sprintf("Cell%s",seq(1:ncol(dat_cell_trace_d)))
    colnames(dat_cell_trace_d)<- cell_id
    
    dat_cell_trace_d_cluster <- unname(kmeans(t(dat_cell_trace_d), centers = 3, nstart = 50)[[1]])
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
    
    rep_time <- c(ncol(dat_trace_m3[[k]]), ncol(dat_trace_m7[[k]]), ncol(dat_trace_m17[[k]]), ncol(dat_trace_m18[[k]]), ncol(dat_trace_m855[[k]]), ncol(dat_trace_m857[[k]]))
    dat_cell_trace_d_re <- mutate(dat_cell_trace_d_re, Group= factor(Group, levels = c("Excited", "Neutral", "Inhibited"))) %>% 
      mutate(., ID = rep(mouse_ID, rep_time*length(stim_time)) )
    
    dat_cell_trace_re[[k]] <- dat_cell_trace_d_re
  }
  
  dat_cell_trace_pre_random[[i]] <-  ddply(dat_cell_trace_re[[1]], .(ID,Time, Group), summarise,value=mean(value)) 
  dat_cell_trace_test_random[[i]] <- ddply(dat_cell_trace_re[[2]], .(ID,Time, Group), summarise,value=mean(value))

  
}

dat_cell_trace_pre_random_sta <- do.call(rbind, dat_cell_trace_pre_random) %>% 
  ddply(., .(ID, Time, Group), summarise, value1=mean(value)) %>% 
  ddply(., .(ID, Time), summarise,mean_pre=sum(value1)) 
  

dat_cell_trace_test_random_sta <- do.call(rbind, dat_cell_trace_test_random) %>% 
  ddply(., .(ID, Time, Group), summarise, value1=mean(value)) %>% 
  ddply(., .(ID, Time), summarise,mean_test=sum(value1))


dat_cell_trace_pre_sta <-  ddply(dat_cell_trace_re[[1]], .(ID,Time, Group), summarise,value1=mean(value)) %>% 
  ddply(., .(ID, Time), summarise,mean_pre1=sum(value1)) %>% 
  mutate(mean_pre = dat_cell_trace_pre_random_sta$mean_pre) %>% 
  mutate(value = mean_pre1 - mean_pre) %>% 
  ddply(., .(Time), summarise,n=length(value),mean=mean(value),sd=sd(value),se=sd(value)/sqrt(length(value))) %>% 
  mutate(Day = "Pre")


dat_cell_trace_test_sta <- ddply(dat_cell_trace_re[[2]], .(ID,Time, Group), summarise,value1=mean(value))%>% 
  ddply(., .(ID, Time), summarise,mean_test1=sum(value1)) %>% 
  mutate(mean_test = dat_cell_trace_test_random_sta$mean_test) %>% 
  mutate(value = mean_test1 - mean_test) %>% 
  ddply(., .(Time), summarise,n=length(value),mean=mean(value),sd=sd(value),se=sd(value)/sqrt(length(value))) %>% 
  mutate(Day = "Test")


p_trace_compare <- rbind(dat_cell_trace_pre_sta, dat_cell_trace_test_sta) %>% 
  ggplot(., aes(Time, mean, colour=Day))+
  geom_line()+
  geom_ribbon(aes(ymin=mean-se, ymax=mean+se, fill=Day), alpha=0.2, linetype=0)+
  scale_colour_manual(values=c("deepskyblue4", "indianred"))+
  scale_fill_manual(values=c("deepskyblue4", "indianred"))+
  labs(x="Time relative to crossing (s)", y="z score", title = "After subtraction")+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  geom_vline(xintercept = 0, col= 'red', linetype=2)+
  theme(legend.title = element_blank(), legend.position = c(0.2, 0.8))
  



p_pre <- ggplot(dat_cell_trace_pre_random_sta, aes(Time, mean_pre,  colour = ID))+
  geom_line()+
  labs(x="Time relative to crossing (s)", y="z score", title = "Averaged from 100 trails: Pre")+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  geom_vline(xintercept = 0, col= 'red', linetype=2)

p_test <- ggplot(dat_cell_trace_test_random_sta, aes(Time, mean_test,  colour = ID))+
  geom_line()+
  labs(x="Time relative to crossing (s)", y="z score", title = "Averaged from 100 trails: Test")+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  geom_vline(xintercept = 0, col= 'red', linetype=2)
  

## random crossing with data from Fatih, placebo state---------
c_miniscope_matlab_ft <- function(file_trace) {
  ## import and format the data
  dat_trace1 <- raveio::read_mat(file_trace)
  ID <- str_extract(file_trace, regex("m\\d+"))
  num_compare <- c(4,7, 8)
  length_pre <- dat_trace1[[4]] %>% 
    nrow()
  t_crossing_pre <- sample(c(40:(length_pre-40)), 100)
  
  length_cond <- dat_trace1[[7]] %>% 
    nrow()
  t_crossing_cond <- sample(c(40:(length_cond-40)), 100)
  
  length_post <- dat_trace1[[8]] %>% 
    nrow()
  t_crossing_post <- sample(c(40:(length_post-40)), 100)
  
  t_crossing <- list(t_crossing_pre, t_crossing_cond, t_crossing_post)
  
  cross_ID <- dat_trace1$global_map %>% 
    as_tibble() %>% 
    select(V1, V4, V5) %>% 
    mutate_all(na_if, 0) %>% 
    drop_na()
  
  dat_stim_trace <- rep(0, length(num_compare))
  
  
  for (i in seq_along(num_compare)) {
    global_cell <- pull(cross_ID[,i])
    dat_trace <- dat_trace1[[num_compare[i]]] %>% 
      as_tibble() %>% 
      select(all_of(global_cell)) %>% 
      apply(., 2, scale)
    
    ## number of rows to be binned
    n <- 2 # 0.05*2=0.1
    
    t_crossing_day <- t_crossing[[i]]
    t_crossing_mean <- rep(0, 100)
    
    for (j in seq_along(t_crossing_day)){
      t1_p <- t_crossing_day[j]
      ## extract cell activity when they cross the border
      #dat_stim1 <- dat_trace[(t1_p-40):(t1_p+140-1),] 
      dat_stim1 <- dat_trace[(t1_p-40):(t1_p+40-1),] 
      
      dat_stim <- aggregate(dat_stim1,list(rep(1:(nrow(dat_stim1)%/%n+1),each=n,len=nrow(dat_stim1))),mean)[-1] %>% 
        apply(., 2, mean) %>% 
        mean()
     
     t_crossing_mean[j] <- dat_stim    
      
      
    }

    dat_stim_trace[i] <- mean(t_crossing_mean)

  }
  return(dat_stim_trace)
}

mouse_file <- as.list(list.files("~cchen/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/FD-processed", full.names = T))

dat_cell_trace <- mapply(c_miniscope_matlab_ft, mouse_file, SIMPLIFY = F)

p_random <- dat_cell_trace %>% 
  do.call(rbind,.) %>% 
  as_tibble() %>% 
  rename(Pre =V1, Cond= V2, Post = V3) %>% 
  mutate(ID = sort(c("m3", "m7", "m17", "m18", "m855", "m857"))) %>% 
  pivot_longer(-ID) %>%
  mutate(name = factor(name, levels = c("Pre", "Cond", "Post"))) %>% 
  ggplot(., aes(name, value, color = name))+
  geom_boxplot(outlier.shape = NA)+
  geom_line(aes(group=ID), colour="gray90")+
  geom_jitter(aes(colour = name, shape = name),width = 0.2,  size=2)+
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

setwd("~cchen/Documents/neuroscience/Pn\ project/Figure/PDF/")
cairo_pdf("p_random_crossings.pdf", width = 50/25.6, height = 65/25.6, family = "Arial")
p_random
dev.off()

## for random crossing control---------
# 1. for the d^2

c_miniscope_matlab_ft <- function(file_trace) {
  ## import and format the data
  dat_trace1 <- raveio::read_mat(file_trace)
  ID <- str_extract(file_trace, regex("m\\d+"))
  num_compare <- c(4, 8)
  length_pre <- dat_trace1[[4]] %>% 
    nrow()
  t_crossing_pre <- tibble(cross = sample(c(40:(length_pre-40)), 100), back = sample(c(40:(length_pre-40)), 100))
  
  length_post <- dat_trace1[[8]] %>% 
    nrow()
  t_crossing_post <- tibble(cross = sample(c(40:(length_post-40)), 100), back = sample(c(40:(length_post-40)), 100))
  
  t_crossing <- list(t_crossing_pre,  t_crossing_post)
  
  cross_ID <- dat_trace1$global_map %>% 
    as_tibble() %>% 
    select(V1, V5) %>% 
    mutate_all(na_if, 0) %>% 
    drop_na()
  
  group <- c("Pre", "Post")
  
  dat_stim_trace <- rep(0, length(num_compare))
  dat_stim_trace_cell <- vector(mode = "list", length = length(num_compare))

  
  for (i in seq_along(num_compare)) {
    
    global_cell <- pull(cross_ID[,i])
    dat_trace <- dat_trace1[[num_compare[i]]] %>% 
      as_tibble() %>% 
      select(all_of(global_cell)) %>% 
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

mouse_file <- as.list(list.files("~cchen/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/FD-processed", full.names = T))

dat_cell_trace <- mapply(c_miniscope_matlab_ft, mouse_file, SIMPLIFY = F)



p_d_combine <- lapply(dat_cell_trace, function(x)  x[[1]]) %>% 
  do.call(rbind,.) %>% 
  as_tibble() %>% 
  rename(Pre = V1, Post = V2) %>% 
  mutate(ID = sort(c("m3", "m7", "m17", "m18", "m855", "m857"))) %>% 
  pivot_longer(-ID) %>% 
  mutate(name = factor(name, levels = c("Pre", "Post"))) %>% 
  ggplot(., aes(name, value, group = name))+
  geom_boxplot(outlier.shape = NA)+
  geom_line(aes(group=ID), colour="gray90")+
  geom_jitter(aes(colour = name, shape = name),width = 0.2,  size=2)+
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
  scale_y_continuous(limits = c(0, 60))+
  theme(legend.title = element_blank(), legend.position = "none")

wilcox.test(mean~Group, dat_cell_d_combine, paired = T, alternative = "less")

## compare by cells

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
  scale_y_continuous(limits = c(-1, 600))+
  theme(legend.title = element_blank(), legend.position = "none")

t_d_combine_cell <- lapply(dat_cell_trace, function(x)  x[[2]]) %>% 
  do.call(rbind,.) %>% 
  mutate(Group = "Post") %>% 
  rbind(dat_cell_d_pre, .) %>% 
  mutate(Group = factor(Group, levels = c("Pre", "Post"))) %>% 
  wilcox.test(d2~Group, ., paired = T, alternative = "less")

p_d_align <- plot_grid(p_d_combine, p_d_combine_cell, nrow = 1,rel_widths = c(0.8, 1) )

setwd("~cchen/Documents/neuroscience/Pn\ project/Figure/PDF/")
cairo_pdf("p_d_align.pdf", width = 90/25.6, height = 60/25.6, family = "Arial")
p_d_align
dev.off()

## compare crossing, back-crossing and last crossing---------


c_miniscope_matlab_ft <- function(file_trace) {
  ## import and format the data
  dat_trace1 <- raveio::read_mat(file_trace)
  
  
  length_post <- dat_trace1[[8]] %>% 
    nrow()
  t_crossing_post <- sample(c(40:(length_post-40)), 100)
  t_crossing_back <- sample(c(40:(length_post-40)), 100)
  t_crossing_last <- sample(c(40:(length_post-40)), 100)
  
  t_crossing <- tibble(cross = t_crossing_post, cross_back = t_crossing_back, last_cross = t_crossing_last)
  
  
  
  
  dat_stim_trace <- rep(0, 3)
  
  
  for (i in c(1:3)) {
    dat_trace <- dat_trace1[[8]] %>% 
      as_tibble() %>% 
      apply(., 2, scale)
    
    n <- 2 # 0.05*2=0.1
    
    mean_cell_acti <- rep(0, 100)
    
    for (j in c(1:100)) {
      t1_p <- t_crossing[j,i] %>% 
        unlist()
      
      dat_stim1 <- dat_trace[(t1_p-40):(t1_p+40-1),] 
      
      dat_stim <- aggregate(dat_stim1,list(rep(1:(nrow(dat_stim1)%/%n+1),each=n,len=nrow(dat_stim1))),mean)[-1]
      # apply(., 2, scale) %>% 
      # as_tibble() %>% 
      # replace(is.na(.), 0)
      
      mean_cell_acti[j] <- dat_stim %>% 
        apply(., 2, mean) %>% 
        mean()
      
    }
    
    dat_stim_trace[i] <- mean_cell_acti %>% 
      mean()
    
  }
  return(dat_stim_trace)
}


## analyze for all mouse 
mouse_file <- as.list(list.files("~cchen/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/FD-processed", full.names = T))

dat_cell_trace <- mapply(c_miniscope_matlab_ft, mouse_file, SIMPLIFY = F)

p_dat_anti_cross <- dat_cell_trace%>% 
  do.call(rbind,.) %>% 
  as_tibble() %>% 
  rename(Crossing = V1, Crossing_back = V2, Last_crossing = V3) %>% 
  mutate(ID = sort(c("m3", "m7", "m17", "m18", "m855", "m857"))) %>% 
  pivot_longer(-ID) %>% 
  mutate(name = factor(name, levels = c("Crossing", "Crossing_back", "Last_crossing"))) %>% 
  ggplot(., aes(name, value, group = name))+
  geom_boxplot(outlier.shape = NA)+
  geom_line(aes(group=ID), colour="gray90")+
  geom_jitter(aes(colour = name, shape = name),width = 0.2,  size=2)+
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


## d^2 cross-day aligned neurons-----
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
  ## for the spike frequency file
  dat_spike <- list.files("~cchen/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/Extract/Events_0.5/", pattern = ID, full.names = T)[c(1,3)]
  
  
  
  dat_stim_trace <- vector(mode = "list", length = length(num_compare))
  
  
  for (i in seq_along(num_compare)) {
    
    global_cell <- pull(cross_ID[,i])
    dat_trace <- read.xlsx(dat_spike[i], colNames = F) %>% 
      as_tibble() %>% 
      select(all_of(global_cell))
    
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
  scale_shape_manual(values=c(16,15))+
  labs(x="", y="(d')^2")+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  scale_y_continuous(limits = c(0, 8))+
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
  # scale_y_continuous(limits = c(-0.5, 0.6))+
  theme(legend.title = element_blank(), legend.position = "none")

t_d_combine_cell <- lapply(dat_cell_trace, function(x)  x[[2]]) %>% 
  do.call(rbind,.) %>% 
  mutate(Group = "Post") %>% 
  rbind(dat_cell_d_pre, .) %>% 
  mutate(Group = factor(Group, levels = c("Pre", "Post"))) %>% 
  wilcox.test(d2~Group, ., paired = T, alternative = "less")

p_d_align <- plot_grid(p_d_combine, p_d_combine_cell, nrow = 1, rel_widths = c(0.8,1))

setwd("~cchen/Documents/neuroscience/Pn\ project/Figure/PDF/")
cairo_pdf("p_d_align.pdf", width = 90/25.6, height = 60/25.6, family = "Arial")
p_d_align
dev.off()

## random crossing with data from Fatih, placebo state, dor firing frequencey---------
c_miniscope_matlab_ft <- function(file_trace) {
  ## import and format the data
  dat_trace1 <- raveio::read_mat(file_trace)
  ID <- str_extract(file_trace, regex("m\\d+"))
  num_compare <- c(4,7, 8)
  length_pre <- dat_trace1[[4]] %>% 
    nrow()
  t_crossing_pre <- sample(c(40:(length_pre-40)), 100)
  
  length_cond <- dat_trace1[[7]] %>% 
    nrow()
  t_crossing_cond <- sample(c(40:(length_cond-40)), 100)
  
  length_post <- dat_trace1[[8]] %>% 
    nrow()
  t_crossing_post <- sample(c(40:(length_post-40)), 100)
  
  t_crossing <- list(t_crossing_pre, t_crossing_cond, t_crossing_post)
  
  cross_ID <- dat_trace1$global_map %>% 
    as_tibble() %>% 
    select(V1, V4, V5) %>% 
    mutate_all(na_if, 0) %>% 
    drop_na()
  dat_spike <- list.files("~cchen/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/Extract/Events_0.5/", pattern = ID, full.names = T)
  
  
  dat_stim_trace <- rep(0, length(num_compare))
  
  
  for (i in seq_along(num_compare)) {
    global_cell <- pull(cross_ID[,i])
    dat_event <- read.xlsx(dat_spike[i], colNames = F) %>% 
      as_tibble() %>% 
      select(all_of(global_cell))
    
    ## number of rows to be binned
    n <- 2 # 0.05*2=0.1
    
    t_crossing_day <- t_crossing[[i]]
    t_crossing_mean <- rep(0, 100)
    
    for (j in seq_along(t_crossing_day)){
      t1_p <- t_crossing_day[j]
      ## extract cell activity when they cross the border
      #dat_stim1 <- dat_trace[(t1_p-40):(t1_p+140-1),] 
      dat_stim1 <- dat_event[(t1_p-40):(t1_p+40-1),] 
      
      dat_stim <- aggregate(dat_stim1,list(rep(1:(nrow(dat_stim1)%/%n+1),each=n,len=nrow(dat_stim1))),mean)[-1] %>% 
        apply(., 2, mean) %>% 
        mean()
      
      t_crossing_mean[j] <- dat_stim    
      
      
    }
    
    dat_stim_trace[i] <- mean(t_crossing_mean)
    
  }
  return(dat_stim_trace)
}

mouse_file <- as.list(list.files("~cchen/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/FD-processed", full.names = T))

dat_cell_trace <- mapply(c_miniscope_matlab_ft, mouse_file, SIMPLIFY = F)

p_random <- dat_cell_trace %>% 
  do.call(rbind,.) %>% 
  as_tibble() %>% 
  rename(Pre =V1, Cond= V2, Post = V3) %>% 
  mutate(ID = sort(c("m3", "m7", "m17", "m18", "m855", "m857"))) %>% 
  pivot_longer(-ID) %>%
  mutate(name = factor(name, levels = c("Pre", "Cond", "Post"))) %>% 
  ggplot(., aes(name, value, color = name))+
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
  scale_y_continuous(limits = c(-0.1, 0.6))+
  theme(legend.title = element_blank(), legend.position = "none")

setwd("~cchen/Documents/neuroscience/Pn\ project/Figure/PDF/")
cairo_pdf("p_random_crossings.pdf", width = 50/25.6, height = 65/25.6, family = "Arial")
p_random
dev.off()



## for d^2 with firing event----
c_miniscope_matlab_ft <- function(file_trace) {
  ## import and format the data
  dat_trace1 <- raveio::read_mat(file_trace)
  ID <- str_extract(file_trace, regex("m\\d+"))
  num_compare <- c(4, 8)
  length_pre <- dat_trace1[[4]] %>% 
    nrow()
  t_crossing_pre <- tibble(cross = sample(c(40:(length_pre-40)), 100), back = sample(c(40:(length_pre-40)), 100))
  
  length_post <- dat_trace1[[8]] %>% 
    nrow()
  t_crossing_post <- tibble(cross = sample(c(40:(length_post-40)), 100), back = sample(c(40:(length_post-40)), 100))
  
  t_crossing <- list(t_crossing_pre,  t_crossing_post)
  
  cross_ID <- dat_trace1$global_map %>% 
    as_tibble() %>% 
    select(V1, V5) %>% 
    mutate_all(na_if, 0) %>% 
    drop_na()
  
  group <- c("Pre", "Post")
  dat_spike <- list.files("~cchen/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/Extract/Events_0.5/", pattern = ID, full.names = T)[c(1,3)]
  
  
  dat_stim_trace <- rep(0, length(num_compare))
  dat_stim_trace_cell <- vector(mode = "list", length = length(num_compare))
  
  
  for (i in seq_along(num_compare)) {
    
    global_cell <- pull(cross_ID[,i])
    dat_trace <- read.xlsx(dat_spike[i], colNames = F) %>% 
      as_tibble() %>% 
      select(all_of(global_cell))
    
    
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

mouse_file <- as.list(list.files("~cchen/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/FD-processed", full.names = T))

dat_cell_trace <- mapply(c_miniscope_matlab_ft, mouse_file, SIMPLIFY = F)



p_d_combine <- lapply(dat_cell_trace, function(x)  x[[1]]) %>% 
  do.call(rbind,.) %>% 
  as_tibble() %>% 
  rename(Pre = V1, Post = V2) %>% 
  mutate(ID = sort(c("m3", "m7", "m17", "m18", "m855", "m857"))) %>% 
  pivot_longer(-ID) %>% 
  mutate(name = factor(name, levels = c("Pre", "Post"))) %>% 
  ggplot(., aes(name, value, group = name))+
  geom_boxplot(outlier.shape = NA)+
  geom_line(aes(group=ID), colour="gray90")+
  geom_jitter(aes(colour = name, shape = name),width = 0.2,  size=2)+
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
  scale_y_continuous(limits = c(0, 8))+
  theme(legend.title = element_blank(), legend.position = "none")

wilcox.test(mean~Group, dat_cell_d_combine, paired = T, alternative = "less")
