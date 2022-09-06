## this script is for analzing Fatih's analysis

library(raveio)

c_miniscope_matlab_ft <- function(file_trace) {
  ## import and format the data
  
  ID <- str_extract(file_trace, regex("m\\d+"))
  dat_trace1 <- raveio::read_mat(file_trace)
  num_compare <- c(4, 8)
  t_crossing <- c(1,5)
  
  dat_stim_trace <- vector(mode = "list", length = length(num_compare))
  
  
  for (i in seq_along(num_compare)) {
    dat_trace <- dat_trace1[[num_compare[i]]] %>% 
      apply(., 2, scale)
    
    ## divide ca2+ trace within -2 and 6s before and after stimuls, take -2 to 0 as baseline
    stim_time<- seq(-2, 6.5, by=0.5)
    ## number of rows to be binned
    n <- 10 # 0.05*10=0.5
    t1_p <- dat_trace1$crossing[t_crossing[i]]
    ## extract cell activity when they cross the border
    dat_stim1 <- dat_trace[(t1_p-40):(t1_p+140-1),] 
    dat_stim1 <- aggregate(dat_stim1,list(rep(1:(nrow(dat_stim1)%/%n+1),each=n,len=nrow(dat_stim1))),mean)[-1]
    dat_stim <- dat_stim1[1:5,] %>%  ## baseline as -2 to 0
      colMeans(., na.rm = T) %>% 
      sweep(dat_stim1, 2, ., FUN = "-")
    
    colnames(dat_stim) <- str_c(ID,"Cell", 1: ncol(dat_stim))
    dat_stim_trace[[i]] <- dat_stim
    }
    return(dat_stim_trace)
  }

## analyze for all mouse 
mouse_file <- as.list(list.files("~cchen/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/FD-processed", full.names = T))

dat_cell_trace <- mapply(c_miniscope_matlab_ft, mouse_file, SIMPLIFY = F)

# ## 
# dat_trace_m3 <- c_miniscope_matlab_ft("~cchen/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/FD-processed/m3.mat")
# dat_trace_m7 <- c_miniscope_matlab_ft("~cchen/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/FD-processed/m7.mat")
# dat_trace_m17 <- c_miniscope_matlab_ft("~cchen/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/FD-processed/m17.mat")
# dat_trace_m18 <- c_miniscope_matlab_ft("~cchen/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/FD-processed/m18.mat")
# dat_trace_m855 <- c_miniscope_matlab_ft("~cchen/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/FD-processed/m855.mat")
# dat_trace_m857 <- c_miniscope_matlab_ft("~cchen/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/FD-processed/m857.mat")

mouse_ID <- sort(c("m3", "m7", "m17", "m18", "m855", "m857"))

dat_cell_trace_re <- vector(mode = "list", 2)
for (i in c(1,2)) {
  dat_cell_trace_d <- map(dat_cell_trace, i) %>% 
    do.call(cbind,.)
  cell_id <- sprintf("Cell%s",seq(1:ncol(dat_cell_trace_d)))
  colnames(dat_cell_trace_d)<- cell_id
  
  dat_cell_trace_d_cluster <- unname(kmeans(t(dat_cell_trace_d), centers = 3)[[1]])
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
  
  rep_time <- map(dat_cell_trace, i) %>% 
    mapply(ncol, .)
  dat_cell_trace_d_re <- mutate(dat_cell_trace_d_re, Group= factor(Group, levels = c("Excited", "Neutral", "Inhibited"))) %>% 
    mutate(., ID = rep(mouse_ID, rep_time*length(stim_time)) )
  
  dat_cell_trace_re[[i]] <- dat_cell_trace_d_re
}


heat_m_ID <- "m855"
score_rang1 <- dat_cell_trace_re[[1]] 

score_range <- dat_cell_trace_re[[2]] %>% 
  rbind(., score_rang1) %>% 
  .$value %>% 
  range()
group_day <- c("Pre", "Test")

for (i in c(1: 2)) {
  dat_trace <- dat_cell_trace_re[[i]] 
  
  dat_trace_sta <- dat_trace %>% 
    ddply(., .(variable, Group), summarise,mean=mean(value), sum=sum(value))
  dat_trace_sta <- dat_trace_sta[order(dat_trace_sta[,'mean']),]
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


## for statistical test

c_miniscope_matlab_ft <- function(file_trace) {
  ## import and format the data
  dat_trace1 <- raveio::read_mat(file_trace)
  num_compare <- c(4, 8)
  t_crossing <- c(1,5)
  
  dat_stim_trace <- vector(mode = "list", length = length(num_compare))
  
  
  for (i in seq_along(num_compare)) {
    dat_trace <- dat_trace1[[4]] %>% 
      apply(., 2, scale)
    
    ## divide ca2+ trace within -2 and 6s before and after stimuls, take -2 to 0 as baseline
    stim_time<- seq(-2, 6.5, by=0.5)
    ## number of rows to be binned
    n <- 10 # 0.05*10=0.5
    t1_p <- dat_trace1$crossing[t_crossing[i]]
    ## extract cell activity when they cross the border
    dat_stim1 <- dat_trace[(t1_p-40):(t1_p+140-1),] 
    dat_stim1 <- aggregate(dat_stim1,list(rep(1:(nrow(dat_stim1)%/%n+1),each=n,len=nrow(dat_stim1))),mean)[-1]
    dat_stim <- dat_stim1[1:5,] %>%  ## baseline as -2 to 0
      colMeans(., na.rm = T) %>% 
      sweep(dat_stim1, 2, ., FUN = "-")
    
    dat_stim_trace[[i]] <- dat_stim
  }
  return(dat_stim_trace)
}


## cluster all the cells based on their activiry during pre, cond and test days-----
c_miniscope_matlab_ft <- function(file_trace) {
  ## import and format the data
  dat_trace1 <- raveio::read_mat(file_trace)
  ID <- str_extract(file_trace, regex("m\\d+"))
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
    
    ## divide ca2+ trace within -2 and 6s before and after stimuls, take -2 to 0 as baseline
    stim_time<- seq(-2, 6.5, by=0.5)
    ## number of rows to be binned
    n <- 5 # 0.05*10=0.5
    t1_p <- dat_trace1$crossing[t_crossing[i]]
    ## extract cell activity when they cross the border
    #dat_stim1 <- dat_trace[(t1_p-40):(t1_p+140-1),] 
    dat_stim1 <- dat_trace[(t1_p-40):(t1_p+140-1),] 
    
    dat_stim1 <- aggregate(dat_stim1,list(rep(1:(nrow(dat_stim1)%/%n+1),each=n,len=nrow(dat_stim1))),mean)[-1]
    dat_stim <- dat_stim1[1:10,] %>%  ## baseline as -2 to 0
      colMeans(., na.rm = T) %>% 
      sweep(dat_stim1, 2, ., FUN = "-")
    
    colnames(dat_stim) <- str_c(ID,"Cell", 1: ncol(dat_stim))
    dat_stim_trace[[i]] <- dat_stim
  }
  return(dat_stim_trace)
}


## analyze for all mouse 
mouse_file <- as.list(list.files("~cchen/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/FD-processed", full.names = T))

dat_cell_trace <- mapply(c_miniscope_matlab_ft, mouse_file, SIMPLIFY = F)

dat_cell_trace_combine <- lapply(dat_cell_trace, function(x)  do.call(rbind, x)) %>% 
  do.call(cbind,.)
                         
                      
library(factoextra)
fviz_nbclust(t(dat_cell_trace_combine), kmeans, method = "silhouette")

dat_cell_trace_cluster <- kmeans(t(dat_cell_trace_combine), centers = 6, iter.max = 10)[[1]]

cluster1_name <- names(dat_cell_trace_cluster[dat_cell_trace_cluster ==6])
dat_cell_average_cluster1 <- dat_cell_trace_combine %>% 
  select(all_of(cluster1_name)) %>%
  rowMeans() 

dat_cluster_trace <- c()

for (i in c(1:6)){
  cluster1_name <- names(dat_cell_trace_cluster[dat_cell_trace_cluster ==i])
  dat_cell_value <- dat_cell_trace_combine %>% 
    select(all_of(cluster1_name)) %>% 
    as_tibble() %>% 
    mutate(Cluster = str_c("C", i),Group = rep(c("Pre",  "Cond", "Test"), each = 36), Time = rep(seq(-2, 6.75, by=0.25), 3)) %>% 
    pivot_longer(-c(Time, Cluster, Group))
  
  dat_cluster_trace <- rbind(dat_cluster_trace, dat_cell_value)
  
}



setwd("~cchen/Documents/neuroscience/Pn\ project/Figure/PDF/")
cairo_pdf("p_cluster.pdf", width = 120/25.6, height = 70/25.6, family = "Arial")
p_cluster
dev.off() 

## heatmap plot of all the neurons
# order the neurons by their clustering group and activity
dat_cell_cluster_order <- dat_cell_trace_combine %>% 
  slice(1:36) %>% 
  as_tibble() %>% 
  add_column(Time = seq(-2, 6.75, 0.25)) %>% 
  pivot_longer(-Time) %>% 
  add_column(Cluster = rep(str_c("C", dat_cell_trace_cluster), 36)) %>% 
  ddply(.,.(Cluster), summarise, mean = mean(value)) %>% 
  arrange(mean)
  
  
dat_cell_order <- dat_cell_trace_combine %>% 
  slice(1:36) %>% 
  as_tibble() %>% 
  add_column(Time = seq(-2, 6.75, 0.25)) %>% 
  pivot_longer(-Time) %>% 
  add_column(Cluster = rep(str_c("C", dat_cell_trace_cluster), 36)) %>% 
  ddply(.,.(name, Cluster), summarise, mean = mean(value)) %>% 
  mutate(Cluster = factor(Cluster, levels = dat_cell_cluster_order$Cluster)) %>% 
  group_by(Cluster) %>% 
  arrange(mean, .by_group = T) %>% 
  select(name, Cluster)

p_cluster <- dat_cluster_trace %>% 
  ddply(., .(Time, Cluster, Group), summarise,n=length(value),mean=mean(value),sd=sd(value),se=sd(value)/sqrt(length(value))) %>%
  mutate(Cluster = factor(Cluster, levels = rev(dat_cell_cluster_order$Cluster))) %>% 
  ggplot(., aes(Time, mean, color= Group))+
  geom_line()+
  geom_ribbon(aes(ymin=mean-se, ymax=mean+se, fill=Group), alpha=0.1, linetype=0)+
  facet_wrap(~Cluster, nrow = 2)+
  theme_void()+
  theme(legend.title = element_blank())+
  annotate(x=c(5,5,7,7), y=c(-2),"path")+
  annotate(x=c(7), y=c(-2, -2, 0, 0),"path")


setwd("~cchen/Documents/neuroscience/Pn\ project/Figure/PDF/")
cairo_pdf("p_cluster.pdf", width = 90/25.6, height = 90/25.6, family = "Arial")
p_cluster
dev.off() 


## heatmap plot of three groups
score_range <- dat_cell_trace_combine %>% 
  range()


dat_rect <- dat_cell_order %>% 
  ddply(., .(Cluster), summarise, n = length(name)) %>% 
  add_column(xmin = 7, xmax = 7.2) %>% 
  mutate(ymax = cumsum(n) ) %>% 
  mutate(ymin = c(1, lag(ymax)[-1]))

p_pro_cluster<- ggplot(dat_rect, aes(x="", y=n, fill=Cluster)) +
  geom_bar(stat="identity", width=1, color="white") +
  coord_polar("y", start=0) +
  theme_void()

setwd("~cchen/Documents/neuroscience/Pn\ project/Figure/PDF/")
cairo_pdf("p_pro_cluster.pdf", width = 60/25.6, height = 60/25.6, family = "Arial")
p_pro_cluster
dev.off() 

p_heat_com <- dat_cell_trace_combine %>% 
  as_tibble() %>% 
  add_column(Time = rep(seq(-2, 6.75, 0.25), 3))%>%
  add_column(Group = rep(c("Pre", "Cond", "Post"), each = 36)) %>% 
  pivot_longer(-c(Time, Group)) %>% 
  mutate(name = factor(name, levels = dat_cell_order$name)) %>% 
  mutate(Group = factor(Group, levels = c("Pre", "Cond", "Post"))) %>% 
  merge(.,dat_cell_order, by = "name") %>% 
  ggplot(., aes(Time, name,fill= value))+ 
  geom_tile(height=2)+
  facet_grid(col= vars(Group))+
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
  theme(legend.position = "none") +
  annotate("rect", xmin = 7, xmax = 7.2, ymin = dat_rect$ymin, ymax = dat_rect$ymax-1, alpha = .2)

  


setwd("~cchen/Documents/neuroscience/Pn\ project/Figure/PDF/")
cairo_pdf("p_heat.pdf", width = 120/25.6, height = 70/25.6, family = "Arial")
p_heat_com
dev.off()

## for heat bar
setwd("~cchen/Documents/neuroscience/Pn\ project/Figure/PDF/")
cairo_pdf("p_heat_bar.pdf", width = 60/25.6, height = 60/25.6, family = "Arial")
p_heat
dev.off()

## compare 2s before and 2s after-------
c_miniscope_matlab_ft <- function(file_trace) {
  ## import and format the data
  dat_trace1 <- raveio::read_mat(file_trace)
  ID <- str_extract(file_trace, regex("m\\d+"))
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

    
    colnames(dat_stim) <- str_c(ID,"Cell", 1: ncol(dat_stim))
    dat_stim_trace[[i]] <- dat_stim
  }
  return(dat_stim_trace)
}


## analyze for all mouse 
mouse_file <- as.list(list.files("~cchen/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/FD-processed", full.names = T))

dat_cell_trace <- mapply(c_miniscope_matlab_ft, mouse_file, SIMPLIFY = F)

dat_cell_trace_combine <- lapply(dat_cell_trace, function(x)  do.call(rbind, x)) %>% 
  do.call(cbind,.)


library(factoextra)
fviz_nbclust(t(dat_cell_trace_combine), kmeans, method = "silhouette")
k_num <- 3

dat_cell_trace_cluster <- kmeans(t(dat_cell_trace_combine), centers = k_num, iter.max = 10)[[1]]

cluster1_name <- names(dat_cell_trace_cluster[dat_cell_trace_cluster ==k_num])
dat_cell_average_cluster1 <- dat_cell_trace_combine %>% 
  select(all_of(cluster1_name)) %>%
  rowMeans() 

dat_cluster_trace <- c()

for (i in c(1:k_num)){
  cluster1_name <- names(dat_cell_trace_cluster[dat_cell_trace_cluster ==i])
  dat_cell_value <- dat_cell_trace_combine %>% 
    select(all_of(cluster1_name)) %>% 
    as_tibble() %>% 
    mutate(Cluster = str_c("C", i),Group = rep(c("Pre",  "Cond", "Test"), each = 40), Time = rep(seq(-2, 1.9, by=0.1), 3)) %>% 
    pivot_longer(-c(Time, Cluster, Group))
  
  dat_cluster_trace <- rbind(dat_cluster_trace, dat_cell_value)
  
}



## heatmap plot of all the neurons
# order the neurons by their clustering group and activity
dat_cell_cluster_order <- dat_cell_trace_combine %>% 
  slice(1:40) %>% 
  as_tibble() %>% 
  add_column(Time = seq(-2, 1.9, 0.1)) %>% 
  pivot_longer(-Time) %>% 
  add_column(Cluster = rep(str_c("C", dat_cell_trace_cluster), 40)) %>% 
  ddply(.,.(Cluster), summarise, mean = mean(value)) %>% 
  arrange(mean)


dat_cell_order <- dat_cell_trace_combine %>% 
  slice(1:40) %>% 
  as_tibble() %>% 
  add_column(Time = seq(-2, 1.9, 0.1)) %>% 
  pivot_longer(-Time) %>% 
  add_column(Cluster = rep(str_c("C", dat_cell_trace_cluster), 40)) %>% 
  ddply(.,.(name, Cluster), summarise, mean = mean(value)) %>% 
  mutate(Cluster = factor(Cluster, levels = dat_cell_cluster_order$Cluster)) %>% 
  group_by(Cluster) %>% 
  arrange(mean, .by_group = T) %>% 
  select(name, Cluster)

p_cluster <- dat_cluster_trace %>% 
  ddply(., .(Time, Cluster, Group), summarise,n=length(value),mean=mean(value),sd=sd(value),se=sd(value)/sqrt(length(value))) %>%
  mutate(Cluster = factor(Cluster, levels = rev(dat_cell_cluster_order$Cluster))) %>% 
  ggplot(., aes(Time, mean, color= Group))+
  geom_line()+
  geom_ribbon(aes(ymin=mean-se, ymax=mean+se, fill=Group), alpha=0.1, linetype=0)+
  facet_wrap(~Cluster, nrow = 2)+
  theme_void()+
  theme(legend.title = element_blank())+
  annotate(x=c(5,5,7,7), y=c(-2),"path")+
  annotate(x=c(7), y=c(-2, -2, 0, 0),"path")


setwd("~cchen/Documents/neuroscience/Pn\ project/Figure/PDF/")
cairo_pdf("p_cluster.pdf", width = 90/25.6, height = 90/25.6, family = "Arial")
p_cluster
dev.off() 


## heatmap plot of three groups
score_range <- dat_cell_trace_combine %>% 
  range()


dat_rect <- dat_cell_order %>% 
  ddply(., .(Cluster), summarise, n = length(name)) %>% 
  add_column(xmin = 2, xmax = 2.2) %>% 
  mutate(ymax = cumsum(n) ) %>% 
  mutate(ymin = c(1, lag(ymax)[-1]))

p_pro_cluster<- ggplot(dat_rect, aes(x="", y=n, fill=Cluster)) +
  geom_bar(stat="identity", width=1, color="white") +
  coord_polar("y", start=0) +
  theme_void()

setwd("~cchen/Documents/neuroscience/Pn\ project/Figure/PDF/")
cairo_pdf("p_pro_cluster.pdf", width = 60/25.6, height = 60/25.6, family = "Arial")
p_pro_cluster
dev.off() 

p_heat_com <- dat_cell_trace_combine %>% 
  as_tibble() %>% 
  add_column(Time = rep(seq(-2, 1.9, 0.1), 3))%>%
  add_column(Group = rep(c("Pre", "Cond", "Post"), each = 40)) %>% 
  pivot_longer(-c(Time, Group)) %>% 
  mutate(name = factor(name, levels = dat_cell_order$name)) %>% 
  mutate(Group = factor(Group, levels = c("Pre", "Cond", "Post"))) %>% 
  merge(.,dat_cell_order, by = "name") %>% 
  ggplot(., aes(Time, name,fill= value))+ 
  geom_tile(height=2)+
  facet_grid(col= vars(Group))+
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
  theme(legend.position = "none") +
  annotate("rect", xmin = 2, xmax = 2.2, ymin = dat_rect$ymin, ymax = dat_rect$ymax-1, alpha = .2)




setwd("~cchen/Documents/neuroscience/Pn\ project/Figure/PDF/")
cairo_pdf("p_heat.pdf", width = 120/25.6, height = 70/25.6, family = "Arial")
p_heat_com
dev.off()

## compare 2s before and 2s after of cond and post-------
c_miniscope_matlab_ft <- function(file_trace) {
  ## import and format the data
  dat_trace1 <- raveio::read_mat(file_trace)
  ID <- str_extract(file_trace, regex("m\\d+"))
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

## compare the activity cross days

comp_time <- seq(-2, 1.9, by=0.1)
leng_time <- length(comp_time)
dat_cell_trace_combine_d6 <- lapply(dat_cell_trace, function(x)  x[[2]]) %>% 
  do.call(cbind,.) %>% 
  mutate(Time = comp_time) %>% 
  mutate(Group = "Cond")

dat_cell_trace_combine_d3 <- lapply(dat_cell_trace, function(x)  x[[1]]) %>% 
  do.call(cbind,.) %>% 
  mutate(Time = comp_time) %>% 
  mutate(Group = "Pre")

dat_cell_trace_combine_d7 <- lapply(dat_cell_trace, function(x)  x[[3]]) %>% 
  do.call(cbind,.) %>% 
  mutate(Time = comp_time) %>% 
  mutate(Group = "Post")

## venn plot the cells show more activity
dat_trace_combine <- bind_rows(dat_cell_trace_combine_d3, dat_cell_trace_combine_d6, dat_cell_trace_combine_d7) %>% 
  pivot_longer(-c(Time, Group)) 

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

dat_cell_compare <- dat_trace_combine %>% 
  ddply(., .(name, Group), summarise, mean_acti = mean(value)) %>% 
  pivot_wider(., names_from = Group, values_from = mean_acti) %>% 
  mutate(diff1 = Cond - Pre) %>% 
  mutate(diff2 = Post - Pre) %>% 
  mutate(diff3 = Post - Cond) 

Pre_name_group <- dat_cell_compare %>% 
  dplyr::select(name) %>% 
  mutate(Group = ifelse(name %in% ratio_pre_con$name, 2,1)) %>% 
  mutate(Group = ifelse(name %in% ratio_pre_post$name & Group ==1, 4, Group)) %>% 
  mutate(Group = ifelse(name %in% ratio_pre_post$name & Group ==2, 3, Group)) %>% 
  mutate(Group = str_c("G", Group)) 



Cond_name_group <- dat_cell_compare %>% 
  dplyr::select(name) %>% 
  mutate(Group = ifelse(name %in% ratio_pre_con$name, 2,6)) %>% 
  mutate(Group = ifelse(name %in% ratio_con_post$name & Group ==6, 5, Group)) %>% 
  mutate(Group = ifelse(name %in% ratio_con_post$name & Group ==2, 3, Group)) %>% 
  mutate(Group = str_c("G", Group)) 


Pre_name_group3 <- Pre_name_group %>% 
  filter(Group == "G3")
Pre_name_group4 <- Pre_name_group %>% 
  filter(Group == "G4")
Cond_name_group5 <- Cond_name_group %>% 
  filter(Group == "G5")

Post_name_group <- dat_cell_compare %>% 
  dplyr::select(name) %>% 
  mutate(Group = "G7") %>% 
  mutate(Group = ifelse(name %in% Pre_name_group3$name, "G3",Group)) %>% 
  mutate(Group = ifelse(name %in% Pre_name_group4$name, "G4", Group)) %>% 
  mutate(Group = ifelse(name %in% Cond_name_group5$name, "G5", Group)) %>% 
  mutate(name = str_c(name, Group))
Pre_name_group <- Pre_name_group %>% 
  mutate(name = str_c(name, Group))

Cond_name_group <- Cond_name_group %>% 
  mutate(name = str_c(name, Group))



library("ggvenn")
dat_ven <- list(Pre = Pre_name_group$name, Cond = Cond_name_group$name, Post = Post_name_group$name )
ggvenn::ggvenn(dat_ven,  c("Pre", "Cond", "Post"), fill_color = c("#8491B4FF", "#00A087FF", "blue"), stroke_color = "black")





p_cell_combine <- bind_rows(dat_cell_trace_combine_d3, dat_cell_trace_combine_d6, dat_cell_trace_combine_d7) %>% 
  pivot_longer(-c(Time, Group)) %>% 
  ddply(., .(name, Group), summarise, mean_acti = mean(value)) %>% 
  ddply(., .(Group), summarise, mean = mean(mean_acti), se=sd(mean_acti)/sqrt(length(mean_acti))) %>% 
  mutate(Group = factor(Group, levels = c("Pre", "Cond", "Post"))) %>%
  ggplot(., aes(Group, mean, fill = Group))+
  geom_bar(stat="identity", color="black", position=position_dodge()) +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.2,position=position_dodge(.9))



## take day6 as the reference for clustering
comp_time <- seq(-2, 1.9, by=0.1)
leng_time <- length(comp_time)
dat_cell_trace_combine_d6 <- lapply(dat_cell_trace, function(x)  x[[2]]) %>% 
  do.call(cbind,.)

library(factoextra)
fviz_nbclust(t(dat_cell_trace_combine_d6), kmeans, method = "silhouette")
k_num <- 2

dat_cell_trace_cluster <- kmeans(t(dat_cell_trace_combine_d6), centers = k_num, iter.max = 10)[[1]]

dat_cell_trace_combine_d6_cluster <- dat_cell_trace_combine_d6 %>% 
  as_tibble() %>% 
  add_column(Time = comp_time) %>% 
  pivot_longer(-Time) %>% 
  add_column(Cluster = rep(dat_cell_trace_cluster, leng_time))

## order the cells based on cluster and activity
dat_cell_cluster_order <- dat_cell_trace_combine_d6_cluster %>% 
  ddply(.,.(Cluster), summarise, mean = mean(value)) %>% 
  arrange(mean)

dat_cell_order <- dat_cell_trace_combine_d6_cluster %>% 
  ddply(.,.(name, Cluster), summarise, mean = mean(value)) %>% 
  mutate(Cluster = factor(Cluster, levels = dat_cell_cluster_order$Cluster)) %>% 
  group_by(Cluster) %>% 
  arrange(mean, .by_group = T) %>% 
  select(name, Cluster)

dat_cell_trace_combine_d5_cluster <- lapply(dat_cell_trace, function(x)  x[[1]]) %>% 
  do.call(cbind,.) %>% 
  as_tibble() %>% 
  add_column(Time = comp_time) %>% 
  pivot_longer(-Time) %>% 
  add_column(Cluster = rep(dat_cell_trace_cluster, leng_time))

dat_cell_trace_combine_cluster <- lapply(dat_cell_trace, function(x)  x[[3]]) %>% 
  do.call(cbind,.) %>% 
  as_tibble() %>% 
  add_column(Time = comp_time) %>% 
  pivot_longer(-Time) %>% 
  add_column(Cluster = rep(dat_cell_trace_cluster, leng_time)) %>% 
  bind_rows(dat_cell_trace_combine_d5_cluster, dat_cell_trace_combine_d6_cluster,.) %>% 
  add_column(Day = rep(c("Pre", "Cond", "Post"), each = nrow(dat_cell_trace_combine_d5_cluster))) %>% 
  mutate(Day = factor(Day, levels = c("Pre", "Cond", "Post"))) %>% 
  mutate(name = factor(name, levels = dat_cell_order$name))


## plot the trace
p_cluster <- dat_cell_trace_combine_cluster %>% 
  ddply(., .(Time, Cluster, Day), summarise,n=length(value),mean_acti=mean(value)) %>%
  ddply(., .(Cluster, Day), summarise, mean = mean(mean_acti), sd=sd(mean_acti),se=sd(mean_acti)/sqrt(length(mean_acti))) %>% 
  ggplot(., aes(Day, mean, group = Day,fill= Day))+
  geom_bar(stat="identity", color="black", position=position_dodge()) +
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.2,position=position_dodge(.9))+ 
  facet_wrap(~Cluster)

p_cluster_total <- dat_cell_trace_combine_cluster %>% 
  ddply(., .(name, Day), summarise,n=length(value),mean_acti=mean(value)) %>%
  ddply(., .( Day), summarise,mean=mean(mean_acti)) %>%
  ggplot(., aes(Day, mean, group = Day,fill= Day))+
  geom_boxplot() 

## heatmap plot
score_range <- dat_cell_trace_combine_cluster$value %>% 
  range()

p_heat_com <- dat_cell_trace_combine_cluster %>% 
  ggplot(., aes(Time, name,fill= value))+ 
  geom_tile(height=2)+
  facet_grid(col= vars(Day))+
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



## analyze the calcium imaging data including all neurons--------
c_miniscope_matlab_ft <- function(file_trace) {
  ## import and format the data
  dat_trace1 <- raveio::read_mat(file_trace)
  ID <- str_extract(file_trace, regex("m\\d+"))
  num_compare <- c(4,7, 8)
  t_crossing <- c(1,4, 5)
  

  
  dat_stim_trace <- vector(mode = "list", length = length(num_compare))
  
  
  for (i in seq_along(num_compare)) {
    dat_trace <- dat_trace1[[num_compare[i]]] %>% 
      as_tibble() %>% 
      apply(., 2, scale)
    
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



## take day6 as the reference
comp_time <- seq(-2, 1.9, by=0.1)
leng_time <- length(comp_time)



dat_cell_trace_combine_d6 <- lapply(dat_cell_trace, function(x)  x[[2]]) %>% 
  do.call(cbind,.) %>% 
  as_tibble() %>% 
  add_column(Time = comp_time) %>% 
  add_column(Group = "Cond") %>% 
  pivot_longer(-c(Time, Group))

dat_cell_cluster_order_d6 <- dat_cell_trace_combine_d6 %>% 
  ddply(.,.(name), summarise, mean = mean(value)) %>% 
  arrange(mean)


## for D3
dat_cell_trace_combine_d3_pca <- lapply(dat_cell_trace, function(x)  x[[1]]) %>% 
  do.call(cbind,.) %>% 
  t() %>% 
  prcomp(., scale = F, center = T)

summary(dat_cell_trace_combine_d3_pca)

dat_cell_trace_combine_d6_pca <- lapply(dat_cell_trace, function(x)  x[[2]]) %>% 
  do.call(cbind,.) %>% 
  t() %>% 
  prcomp(., scale = F, center = T)

summary(dat_cell_trace_combine_d6_pca)

dat_cell_trace_combine_d7_pca <- lapply(dat_cell_trace, function(x)  x[[3]]) %>% 
  do.call(cbind,.) %>% 
  t() %>% 
  prcomp(., scale = F, center = T)

summary(dat_cell_trace_combine_d7_pca)


dat_cell_trace_combine_d3 <- lapply(dat_cell_trace, function(x)  x[[1]]) %>% 
  do.call(cbind,.) %>% 
  as_tibble() %>% 
  add_column(Time = comp_time) %>% 
  add_column(Group = "Pre") %>% 
  pivot_longer(-c(Time, Group))

dat_cell_cluster_order_d3 <- dat_cell_trace_combine_d3 %>% 
  ddply(.,.(name), summarise, mean = mean(value)) %>% 
  arrange(mean)

## for d7
dat_cell_trace_combine_d7 <- lapply(dat_cell_trace, function(x)  x[[3]]) %>% 
  do.call(cbind,.) %>% 
  as_tibble() %>% 
  add_column(Time = comp_time) %>% 
  add_column(Group = "Post") %>% 
  pivot_longer(-c(Time, Group))

dat_cell_cluster_order_d7 <- dat_cell_trace_combine_d7 %>% 
  ddply(.,.(name), summarise, mean = mean(value)) %>% 
  arrange(mean)

score_range <- range(rbind(dat_cell_trace_combine_d3, dat_cell_trace_combine_d6, dat_cell_trace_combine_d7)$value)

p_dat_cell_d6 <- dat_cell_trace_combine_d6 %>% 
  mutate(name = factor(name, levels = dat_cell_cluster_order_d6$name)) %>% 
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
  
p_dat_cell_d3 <- dat_cell_trace_combine_d3 %>% 
  mutate(name = factor(name, levels = dat_cell_cluster_order_d3$name)) %>% 
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

p_dat_cell_d7 <- dat_cell_trace_combine_d7 %>% 
  mutate(name = factor(name, levels = dat_cell_cluster_order_d7$name)) %>% 
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
  
  
p_heat <- plot_grid(p_dat_cell_d3, p_dat_cell_d6, p_dat_cell_d7, ncol = 3)
  
  
## calculate population activity

p_trace <- rbind(dat_cell_trace_combine_d3, dat_cell_trace_combine_d6, dat_cell_trace_combine_d7) %>% 
  ddply(., .(Time, Group),summarise,n=length(value),mean_acti=mean(value),sd=sd(value),se=sd(value)/sqrt(length(value))) %>% 
  ddply(., .(Group), summarise,mean = sum(mean_acti)) %>% 
  mutate(Group = factor(Group, levels = c("Pre", "Cond", "Post"))) %>% 
  ggplot(., aes(Group,  mean, fill = Group))+
  geom_bar(stat="identity", color="black", position=position_dodge())


p_trace <- rbind(dat_cell_trace_combine_d3, dat_cell_trace_combine_d6, dat_cell_trace_combine_d7) %>% 
  ddply(., .(name, Group),summarise,n=length(value),mean_acti=mean(value),sd=sd(value),se=sd(value)/sqrt(length(value))) %>% 
  ddply(., .(Group), summarise,mean = mean(mean_acti), se=sd(mean_acti)/sqrt(length(mean_acti))) %>% 
  mutate(Group = factor(Group, levels = c("Pre", "Cond", "Post"))) %>% 
  ggplot(., aes(Group,  mean, fill = Group))+
  geom_bar(stat="identity", color="black", position=position_dodge())
  # geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.2,position=position_dodge(.9))

  

dat_cell_trace_combine_d6_cluster <- dat_cell_trace_combine_d6 %>% 
  as_tibble() %>% 
  add_column(Time = comp_time) %>% 
  pivot_longer(-Time) %>% 
  add_column(Cluster = rep(dat_cell_trace_cluster, leng_time))

## order the cells based on cluster and activity
dat_cell_cluster_order <- dat_cell_trace_combine_d6_cluster %>% 
  ddply(.,.(Cluster), summarise, mean = mean(value)) %>% 
  arrange(mean)

dat_cell_order <- dat_cell_trace_combine_d6_cluster %>% 
  ddply(.,.(name, Cluster), summarise, mean = mean(value)) %>% 
  mutate(Cluster = factor(Cluster, levels = dat_cell_cluster_order$Cluster)) %>% 
  group_by(Cluster) %>% 
  arrange(mean, .by_group = T) %>% 
  select(name, Cluster)

dat_cell_trace_combine_d5_cluster <- lapply(dat_cell_trace, function(x)  x[[1]]) %>% 
  do.call(cbind,.) %>% 
  as_tibble() %>% 
  add_column(Time = comp_time) %>% 
  pivot_longer(-Time) %>% 
  add_column(Cluster = rep(dat_cell_trace_cluster, leng_time))

dat_cell_trace_combine_cluster <- lapply(dat_cell_trace, function(x)  x[[3]]) %>% 
  do.call(cbind,.) %>% 
  as_tibble() %>% 
  add_column(Time = comp_time) %>% 
  pivot_longer(-Time) %>% 
  add_column(Cluster = rep(dat_cell_trace_cluster, leng_time)) %>% 
  bind_rows(dat_cell_trace_combine_d5_cluster, dat_cell_trace_combine_d6_cluster,.) %>% 
  add_column(Day = rep(c("Pre", "Cond", "Post"), each = nrow(dat_cell_trace_combine_d5_cluster))) %>% 
  mutate(Day = factor(Day, levels = c("Pre", "Cond", "Post"))) %>% 
  mutate(name = factor(name, levels = dat_cell_order$name))


## plot the trace
p_cluster <- dat_cell_trace_combine_cluster %>% 
  ddply(., .(Time, Cluster, Day), summarise,n=length(value),mean_acti=mean(value)) %>%
  ddply(., .(Cluster, Day), summarise, mean = mean(mean_acti), sd=sd(mean_acti),se=sd(mean_acti)/sqrt(length(mean_acti))) %>% 
  ggplot(., aes(Day, mean, group = Day,fill= Day))+
  geom_bar(stat="identity", color="black", position=position_dodge()) +
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.2,position=position_dodge(.9))+ 
  facet_wrap(~Cluster)


ggplot(., aes(Day, mean, group = Day,fill= Day))+
  geom_boxplot() +
  geom_jitter()+
  facet_wrap(~Cluster)

p_cluster_total <- dat_cell_trace_combine_cluster %>% 
  ddply(., .(name, Day), summarise,n=length(value),mean_acti=mean(value)) %>%
  ddply(., .( Day), summarise,mean=mean(mean_acti)) %>%
  ggplot(., aes(Day, mean, group = Day,fill= Day))+
  geom_boxplot() 

## heatmap plot
score_range <- dat_cell_trace_combine %>% 
  range()

p_heat_com <- dat_cell_trace_combine_cluster %>% 
  ggplot(., aes(Time, name,fill= value))+ 
  geom_tile(height=2)+
  facet_grid(col= vars(Day))+
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


## analyze the data mouse by mouse-------
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
comp_time <- seq(-2, 1.9, by=0.1)

cc_heat_function <- function(dat_trace){
  m_id <- dat_trace[[1]] %>% 
    colnames() %>%
    .[1] %>% 
    str_extract(., regex("m\\d+"))
  
  score_range <- range(unlist(dat_trace))
  
  for (i in 1: length(dat_trace)){
    dat_cell_trace <- dat_trace[[i]]%>% 
      as_tibble() %>% 
      add_column(Time = comp_time) %>% 
      pivot_longer(-Time)
    
    dat_cell_order <-  dat_cell_trace %>% 
      ddply(.,.(name), summarise, mean = mean(value)) %>% 
      arrange(mean) %>% 
      select(name)
    
    p_heat <- dat_cell_trace %>% 
      mutate(name = factor(name, levels = dat_cell_order$name)) %>% 
      ggplot(., aes(Time, name,fill= value))+ 
      geom_tile(height=2)+
      scale_fill_gradientn(limits= score_range, colours = c("navy", "white", "red4"), values = rescale(c(score_range[1], 0, score_range[2])))+
      labs(x="Time relative to crossing (s)", y="Number of cells")+
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
    assign(str_c("p_heat", i), p_heat)
    
  }
  
  legend <- get_legend(p_heat1 + guides(color = guide_legend(nrow = 1)) + theme(legend.position = "right"))
  
  p_combine <- plot_grid(p_heat1, p_heat2, p_heat3, ncol = 3, labels = m_id) %>% 
    plot_grid(., legend, nrow = 1,rel_widths = c(3, 0.4))
  
  return(p_combine)
}

heat_map <- mapply(cc_heat_function, dat_cell_trace, SIMPLIFY = F)

p_heat_m3 <- heat_map[[3]]
setwd("~cchen/Documents/neuroscience/Pn\ project/Figure/PDF/")
cairo_pdf("p_heat_m3.pdf", width = 140/25.6, height = 60/25.6, family = "Arial")
p_heat_m3
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
  
  return(list(dat_trace_acti, dat_ratio))
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
  #scale_y_continuous(limits = c(-0.3, 0.6))+
  theme(legend.title = element_blank(), legend.position = "none")

## PCA analysis

cc_pca_fun <- function(dat_trace){
  dat_variance <- mapply(function(x) prcomp(t(x), scale. = F, center = T), dat_trace, SIMPLIFY = F) %>% 
    lapply(., function(x) x[[1]]) %>% 
    mapply(function(x) x^2, ., SIMPLIFY = F) %>% 
    mapply(function(x) x/sum(x), .) %>% 
    as_tibble() %>% 
    rename(Pre= V1, Cond = V2, Post = V3) %>% 
    slice(1:10) %>% 
    add_column(Group = str_c("PC", 1:10))
  
  return(dat_variance)
}


dat_pca_mouse <- mapply(cc_pca_fun, dat_cell_trace, SIMPLIFY = F) %>% 
  do.call(rbind,.) %>% 
  as_tibble() %>% 
  filter(Group == "PC1") %>% 
  add_column(ID = sort(c("m3", "m7", "m17", "m18", "m855", "m857"))) %>%
  dplyr::select(-Group) %>% 
  pivot_longer(-ID) %>% 
  mutate(name = factor(name, levels = c("Pre", "Cond", "Post"))) %>% 
  ggplot(., aes(name, value, group = name))+
  geom_boxplot(outlier.shape = NA)+
  geom_line(aes(group=ID), colour="gray90")+
  geom_jitter(aes(colour = name, shape = name),width = 0.2,  size=2)+
  scale_color_manual(values=c("#8491B4FF", "#00A087FF", "#3C5488FF"))+
  labs(x="", y="Variance explained by PC1")+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  theme(legend.title = element_blank(), legend.position = "none")
  
## correaltion analysis between cell population activity and behavior
dat_anti_pain <- read.xlsx("~cchen/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/matlab_analysis/behavior_tracking/pn_anti_pain_behavior.xlsx", sheet = 2) %>% 
  mutate(Latency_licking = (First_licking - First_crossing)/Frame_rate,
         Latency_rearing = (First_rearing - First_crossing)/Frame_rate) 

dat_ca_behavior <- lapply(dat_acti_combine, function(x) x[[1]]) %>% 
  do.call(rbind,.) %>% 
  mutate(Group = factor(Group, levels = c("Pre", "Cond", "Post"))) %>% 
  filter(Group =="Post") %>% 
  merge(dat_anti_pain, ., by="ID") %>% 
  ggplot(., aes(cell_acti, Latency_licking))+
  geom_point()+
  geom_smooth(method = "lm")


## for d3, d4 ,d5 -------
c_miniscope_matlab_ft <- function(file_trace) {
  ## import and format the data
  dat_trace1 <- raveio::read_mat(file_trace)
  ID <- str_extract(file_trace, regex("m\\d+"))
  dat_trace1[[2]][1] <- ifelse(ID=="m857", 1694, dat_trace1[[2]][1])
  
  num_compare <- c(4,5)
  t_crossing <- c(1,2)
  
  
  
  cross_ID <- dat_trace1$global_map %>% 
    as_tibble() %>% 
    select(V1, V2) %>% 
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
comp_time <- seq(-2, 1.9, by=0.1)

cc_heat_function <- function(dat_trace){
  m_id <- dat_trace[[1]] %>% 
    colnames() %>%
    .[1] %>% 
    str_extract(., regex("m\\d+"))
  
  score_range <- range(unlist(dat_trace))
  
  for (i in 1: length(dat_trace)){
    dat_cell_trace <- dat_trace[[i]]%>% 
      as_tibble() %>% 
      add_column(Time = comp_time) %>% 
      pivot_longer(-Time)
    
    dat_cell_order <-  dat_cell_trace %>% 
      ddply(.,.(name), summarise, mean = mean(value)) %>% 
      arrange(mean) %>% 
      select(name)
    
    p_heat <- dat_cell_trace %>% 
      mutate(name = factor(name, levels = dat_cell_order$name)) %>% 
      ggplot(., aes(Time, name,fill= value))+ 
      geom_tile(height=2)+
      scale_fill_gradientn(limits= score_range, colours = c("navy", "white", "red4"), values = rescale(c(score_range[1], 0, score_range[2])))+
      labs(x="Time relative to crossing (s)", y="Number of cells")+
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
    assign(str_c("p_heat", i), p_heat)
    
  }
  
  p_combine <- plot_grid(p_heat1, p_heat2, ncol = 2, labels = m_id)
  return(p_combine)
}

heat_map <- mapply(cc_heat_function, dat_cell_trace, SIMPLIFY = F)

## calculate cell activity change mouse by mouse
cc_active_fun <- function (dat_trace){
  m_id <- dat_trace[[1]] %>% 
    colnames() %>%
    .[1] %>% 
    str_extract(., regex("m\\d+"))
  
  dat_trace_combine <- dat_trace %>% 
    do.call(rbind,.) %>% 
    as_tibble() %>% 
    add_column(Time = rep(comp_time, 2)) %>% 
    add_column(Group = rep(c("Pre", "Cond"), each = length(comp_time))) %>% 
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
  
  
  dat_ratio <- nrow(ratio_pre_con)/ncol(dat_trace[[1]])
  
  return(list(dat_trace_acti, dat_ratio))
}


dat_acti_combine <- mapply(cc_active_fun, dat_cell_trace, SIMPLIFY = F)

p_dat_anti <- lapply(dat_acti_combine, function(x) x[[1]]) %>% 
  do.call(rbind,.) %>% 
  mutate(Group = factor(Group, levels = c("Pre", "Cond"))) %>% 
  ggplot(., aes(Group, cell_acti, color = Group))+
  geom_boxplot(outlier.shape = NA)+
  geom_line(aes(group = ID), color="gray90")+
  geom_jitter(width = 0.2, shape=1)+
  #scale_colour_manual(values=c("deepskyblue4", "indianred"))+
  labs(x="", y="Population activity (∆F/F)")+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  scale_y_continuous(limits = c(-0.2, 0.4))+
  theme(legend.title = element_blank(), legend.position = "none")

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

cc_heat_function <- function(dat_trace){
  m_id <- dat_trace[[1]] %>% 
    colnames() %>%
    .[1] %>% 
    str_extract(., regex("m\\d+"))
  
  score_range <- range(unlist(dat_trace))
  
  for (i in 1: length(dat_trace)){
    dat_cell_trace <- dat_trace[[i]]%>% 
      as_tibble() %>% 
      add_column(Time = comp_time) %>% 
      pivot_longer(-Time)
    
    dat_cell_order <-  dat_cell_trace %>% 
      ddply(.,.(name), summarise, mean = mean(value)) %>% 
      arrange(mean) %>% 
      select(name)
    
    p_heat <- dat_cell_trace %>% 
      mutate(name = factor(name, levels = dat_cell_order$name)) %>% 
      ggplot(., aes(Time, name,fill= value))+ 
      geom_tile(height=2)+
      scale_fill_gradientn(limits= score_range, colours = c("navy", "white", "red4"), values = rescale(c(score_range[1], 0, score_range[2])))+
      labs(x="Time relative to crossing (s)", y="Number of cells")+
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
    assign(str_c("p_heat", i), p_heat)
    
  }
  
  legend <- get_legend(p_heat1 + guides(color = guide_legend(nrow = 1)) + theme(legend.position = "right"))
  
  p_combine <- plot_grid(p_heat1, p_heat2, p_heat3, ncol = 3, labels = m_id) %>% 
    plot_grid(., legend, nrow = 1,rel_widths = c(3, 0.4))
  
  return(p_combine)
}

heat_map <- mapply(cc_heat_function, dat_cell_trace, SIMPLIFY = F)

p_heat_m3 <- heat_map[[3]]
setwd("~cchen/Documents/neuroscience/Pn\ project/Figure/PDF/")
cairo_pdf("p_heat_m3.pdf", width = 140/25.6, height = 60/25.6, family = "Arial")
p_heat_m3
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

p_dat_anti <- lapply(dat_acti_combine, function(x) x[[1]]) %>% 
  do.call(rbind,.) %>% 
  mutate(Group = factor(Group, levels = c("Crossing", "Crossing_back", "Last_crossing"))) %>% 
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
  scale_y_continuous(limits = c(-0.5, 0.6))+
  theme(legend.title = element_blank(), legend.position = "none")

## ## for d3, d5, d6, d7 and latency to first crossing -------
c_miniscope_matlab_ft <- function(file_trace) {
  ## import and format the data
  dat_trace1 <- raveio::read_mat(file_trace)
  ID <- str_extract(file_trace, regex("m\\d+"))
  dat_trace1[[2]][1] <- ifelse(ID=="m857", 1694, dat_trace1[[2]][1])
  
  num_compare <- c(4,6:8)
  t_crossing <- c(1,3:5)
  
  
  
  cross_ID <- dat_trace1$global_map %>% 
    as_tibble() %>% 
    select(V1, V3, V4, V5) %>% 
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
  
  comp_time <- seq(-2, 1.9, by=0.1)
  
  dat_trace_combine <- dat_stim_trace %>% 
    do.call(rbind,.) %>% 
    as_tibble() %>% 
    add_column(Time = rep(comp_time, 4)) %>% 
    add_column(Group = rep(c("D3", "D5", "D6", "D7"), each = length(comp_time))) %>% 
    pivot_longer(-c(Time, Group)) %>% 
    ddply(., .(name, Group), summarise, mean_acti = mean(value)) %>% 
    ddply(., .(Group), summarise, cell_acti = mean(mean_acti)) %>% 
    add_column(ID = ID) %>% 
    add_column(Latency_cross = dat_trace1$crossing[t_crossing]/20)
  
    return(dat_trace_combine)
}


## analyze for all mouse 
mouse_file <- as.list(list.files("~cchen/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/FD-processed", full.names = T))

dat_cell_trace <- mapply(c_miniscope_matlab_ft, mouse_file, SIMPLIFY = F)

p_ca_crossing <- dat_cell_trace %>% 
  do.call(rbind,.) %>% 
  as_tibble() %>% 
  ggplot(., aes(cell_acti, Latency_cross))+
  geom_point()+
  geom_smooth(method = "lm")+
  labs(x="Mean activity of rACC-Pn neurons (s.d.)", y="Latency of 1st crossing (s)")+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))

setwd("~cchen/Documents/neuroscience/Pn\ project/Figure/PDF/")
cairo_pdf("p_ca_crossing.pdf", width = 90/25.6, height = 60/25.6, family = "Arial")
p_ca_crossing
dev.off()

p_cor <-  dat_cell_trace %>% 
  do.call(rbind,.) %>% 
  as_tibble()

cor.test(p_cor$cell_acti, p_cor$Latency_cross)  

  
  
  
  
  