## Using PCA to analyze the neural population trajectory during crossing with data from Fatih

### test diff time window-----------
c_miniscope_matlab_ft <- function(file_trace) {
  ## import and format the data
  ID <- str_extract(file_trace, regex("m\\d+"))
  
  dat_trace1 <- raveio::read_mat(file_trace)
  num_compare <- c(4,7, 8)
  t_crossing <- c(1,4, 5)
  
  cross_ID <- dat_trace1$global_map %>% 
    as_tibble() %>% 
    dplyr::select(V1, V4, V5) %>% 
    mutate_all(na_if, 0) %>% 
    drop_na()
  
  dat_stim_trace <- vector(mode = "list", length = length(num_compare))
  
  
  for (i in seq_along(num_compare)) {
    global_cell <- pull(cross_ID[,i])
    dat_trace <- dat_trace1[[num_compare[i]]] %>% 
      as_tibble() %>% 
      dplyr::select(all_of(global_cell)) %>% 
      apply(., 2, scale)
    
    n <- 5 # 0.05*10=0.5
    t1_p <- dat_trace1$crossing[t_crossing[i]]
    ## extract cell activity when they cross the border
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

mouse_file <- as.list(list.files("~cchen/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/FD-processed", full.names = T))

dat_cell_trace <- mapply(c_miniscope_matlab_ft, mouse_file, SIMPLIFY = F)

## trajector function
cc_trajectory <- function(x){
  ## varaince explained by PCA first two components
  
  dat_variance <- mapply(function(x) prcomp(t(x), scale. = F, center = T), x, SIMPLIFY = F) %>% 
    lapply(., function(x) x[[1]]) %>% 
    mapply(function(x) x^2, ., SIMPLIFY = F) %>% 
    mapply(function(x) x/sum(x), .) %>% 
    as_tibble() %>% 
    rename(Pre= V1, Cond = V2, Post = V3) %>% 
    slice(1:10) %>% 
    add_column(Group = str_c("PC", 1:10))
    
  

  dat_total <- x %>% 
    do.call(rbind, .) %>% 
    base::as.matrix()
  
  dat_total_eigen <- dat_total %>% 
    prcomp(., scale. = F, center = T) %>% 
    .$rotation 
  
  dat_trajector <- dat_total %*% dat_total_eigen %>% 
    as_tibble() %>% 
    mutate(Group = c(rep(c("Pre","Cond", "Post"), each = 36))) %>% 
    mutate(Time = rep(seq(-2, 6.75, by=0.25), 3))
  
  
  dat_trajector_start <- dat_trajector %>% 
    filter(Time == -2.0 | Time == 0 | Time == 6.75) %>% 
    mutate(Catlog = rep(c("Begin","Crossing", "End"), 3 ))
  
  
  p_2d_trajector <- dat_trajector %>% 
    ggplot(., aes(x = PC1, y = PC2, group = Group,color = Group))+
    geom_line(size = 1)+
    geom_point(data = dat_trajector_start,aes(x=PC1, y = PC2, group = Group,  color = Group, shape = Catlog, size = 1) )
  
  dat_dis <- rep(0, 3) 
  dat_group <- c("Pre", "Cond", "Post")
  for (i in seq_along(dat_group)) {
    dat <- dat_trajector %>% 
      filter(Group ==dat_group[i]) %>% 
      dplyr::select(PC1, PC2) %>% 
      dist(., diag = T, upper = T) %>% 
      as.matrix() %>% 
      .[,1] %>% 
      sum
    dat_dis[i] <- dat
  }
  
  return(list(p_2d_trajector, dat_dis, dat_variance))
}

dat_trajector_combine <- mapply(cc_trajectory, dat_cell_trace, SIMPLIFY = F)


## for trajector distance-----
dat_trajector_dis <- lapply(dat_trajector_combine, function(x) x[[2]]) %>% 
  do.call(rbind,.) %>% 
  as_tibble() %>% 
  rename(Pre = V1, Cond = V2, Test = V3) %>% 
  mutate(ID = sort(c("m3", "m7", "m17", "m18", "m855", "m857"))) %>% 
  pivot_longer(-ID) %>% 
  mutate(name = factor(name, levels = c("Pre", "Cond", "Test")))

p_trajector_dis <- dat_trajector_dis %>% 
  ggplot(., aes(name, value, color = name))+
  geom_boxplot(outlier.shape = NA)+
  geom_jitter()+
  labs(x="PC1", y="PC2")+
  theme_classic()
  

t_trajectory <- dat_trajector_dis %>% 
  friedman_test(value ~ name | ID)

## for variance explained by pcs-----
dat_variance_pc <- lapply(dat_trajector_combine, function(x) x[[3]]) %>% 
  do.call(rbind,.) %>% 
  as_tibble() %>% 
  pivot_longer(-Group) %>% 
  mutate(name = factor(name, levels = c("Pre", "Cond", "Post"))) %>% 
  mutate(Group = factor(Group, levels = c(str_c("PC", 1:10)))) %>% 
  ddply(., .(name, Group), summarise,n=length(value),mean=mean(value),sd=sd(value),se=sd(value)/sqrt(length(value))) %>% 
  ggplot(., aes(Group, mean, color = Group))+
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=0.1) +
  geom_point(aes(colour = Group, shape = Group),size=2)+
  geom_line()+
  labs(x="Variance ", y="Time spend in chammber 2 (%)")+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  theme(legend.position = 'none')

dat_variance_pc <- lapply(dat_trajector_combine, function(x) x[[3]]) %>% 
  do.call(rbind,.) %>% 
  as_tibble() %>% 
  pivot_longer(-Group) %>% 
  mutate(name = factor(name, levels = c("Pre", "Cond", "Post"))) %>% 
  filter(Group == "PC1") %>% 
  ggplot(., aes(name, value, color = name))+
  geom_boxplot(outlier.shape = NA)+
  geom_jitter()+
  labs(x=" ", y="Variance explained by PC1")+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  theme(legend.position = 'none')


## compare the crossing and crossing back--------
mouse_file <- as.list(list.files("~cchen/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/FD-processed", full.names = T))

t_stim_d7 <- list(m3=  c(346, 1003), m7 = c(236, 1632), m17 = c(157, 2113), m18 = c(124,1923), m855 = c(84,1572), m857 = c(334,1892))
t_stim_d7<-t_stim_d7[sort(names(t_stim_d7))]


c_miniscope_matlab_ft_back <- function(file_trace, t_stim) {
  ## import and format the data
  dat_trace1 <- raveio::read_mat(file_trace)
  num_compare <- c(8)

  
  dat_stim_trace <- vector(mode = "list", length = 2)
  
  
  for (i in seq_along(t_stim)) {
    dat_trace <- dat_trace1[[num_compare]] %>% 
      as_tibble() %>% 
      apply(., 2, scale)
    
    n <- 5 # 0.05*10=0.5
    t1_p <- t_stim[i]
    ## extract cell activity when they cross the border
    dat_stim1 <- dat_trace[(t1_p-40):(t1_p+140-1),] 
    
    dat_stim1 <- aggregate(dat_stim1,list(rep(1:(nrow(dat_stim1)%/%n+1),each=n,len=nrow(dat_stim1))),mean)[-1]
    dat_stim <- dat_stim1[1:10,] %>%  ## baseline as -2 to 0
      colMeans(., na.rm = T) %>% 
      sweep(dat_stim1, 2, ., FUN = "-")
    
    colnames(dat_stim) <- str_c("Cell", 1: ncol(dat_stim))
    dat_stim_trace[[i]] <- dat_stim
  }
  return(dat_stim_trace)
}

dat_trace_d7 <- mapply(c_miniscope_matlab_ft_back, mouse_file, t_stim_d7, SIMPLIFY = F) 

cc_trajectory <- function(x){
  ## varaince explained by PCA first two components
  
  dat_vrance_PC1 <- mapply(function(x) prcomp(t(x), scale. = F, center = T), x, SIMPLIFY = F) %>% 
    lapply(., function(x) x[[1]]) %>% 
    mapply(function(x) x^2, ., SIMPLIFY = F) %>% 
    mapply(function(x) x[1]/sum(x), .)
  
  dat_vrance_PC2 <- mapply(function(x) prcomp(t(x), scale. = F, center = T), x, SIMPLIFY = F) %>% 
    lapply(., function(x) x[[1]]) %>% 
    mapply(function(x) x^2, ., SIMPLIFY = F) %>% 
    mapply(function(x) x[2]/sum(x), .)
  dat_virance <- tibble(PC1 = dat_vrance_PC1, PC2= dat_vrance_PC2)
  
  dat_total <- x %>% 
    do.call(rbind, .) %>% 
    as.matrix()
  
  dat_total_eigen <- dat_total %>% 
    prcomp(., scale. = F, center = T) %>% 
    .$rotation 
  

  dat_trajector <- dat_total %*% dat_total_eigen %>% 
    as_tibble() %>% 
    mutate(Group = c(rep(c("Cross","Cross_back"), each = 36))) %>% 
    mutate(Time = rep(seq(-2, 6.75, by=0.25), 2)) 
  
  dat_trajector_start <- dat_trajector %>% 
    filter(Time == -2.0 | Time == 0 | Time == 6.75) %>% 
    mutate(Catlog = rep(c("Begin","Crossing", "End"), 2 ))
  

  
  p_2d_trajector <- dat_trajector %>% 
    ggplot(., aes(x = PC1, y = PC2, group = Group,color = Group))+
    geom_line(size = 1)+
    geom_point(size = 1, shape =1)+
    geom_point(data = dat_trajector_start,aes(x=PC1, y = PC2, group = Group,  color = Group, shape = Catlog, size = 1) )
  
  dat_dis <- rep(0, 2) 
  dat_group <- c("Cross","Cross_back")
  for (i in seq_along(dat_group)) {
    dat <- dat_trajector %>% 
      filter(Group ==dat_group[i]) %>% 
      select(PC1, PC2) %>% 
      dist(., diag = T, upper = T) %>% 
      as.matrix() %>% 
      .[,1] %>% 
      sum
    dat_dis[i] <- dat
  }
  
  return(list(p_2d_trajector, dat_dis, dat_virance))
}

dat_trajector_combine_d7 <- mapply(cc_trajectory, dat_trace_d7, SIMPLIFY = F)


## for trajector distance
dat_trajector_dis <- lapply(dat_trajector_combine_d7, function(x) x[[2]]) %>% 
  do.call(rbind,.) %>% 
  as_tibble() %>% 
  rename(Cross = V1, Cross_back = V2) %>% 
  mutate(ID = sort(c("m3", "m7", "m17", "m18", "m855", "m857"))) %>% 
  pivot_longer(-ID) %>% 
  mutate(name = factor(name, levels = c("Cross","Cross_back")))

p_trajector_dis <- dat_trajector_dis %>% 
  ggplot(., aes(name, value, color = name))+
  geom_boxplot(outlier.shape = NA)+
  geom_jitter()+
  labs(x="PC1", y="PC2")+
  theme_classic()


t_trajectory <- dat_trajector_dis %>% 
  friedman_test(value ~ name | ID)

## for variance explained by pcs-----
dat_variance_pc <- lapply(dat_trajector_combine_d7, function(x) x[[3]]) %>% 
  do.call(rbind,.) %>% 
  as_tibble() %>% 
  mutate(Group = rep(c("Cross","Cross_back"), 6)) %>% 
  mutate(ID = rep(sort(c("m3", "m7", "m17", "m18", "m855", "m857")), each = 2)) %>% 
  mutate(PC1_2 = PC1 + PC2) %>% 
  mutate(Group = factor(Group, levels = c("Cross","Cross_back"))) %>% 
  ggplot(., aes(Group, PC1, color = Group))+
  geom_boxplot()

## for compare 1st and last crossing------
t_stim_d7_2nd <- list(m3=  c(c(346, 1756)), m7 = c(236, 2777), m17 = c(157, 2536), m18 = c(124, 2238), m855 = c(82, 3062), m857 = c(167*2, 1153*2))
t_stim_d7_2nd<-t_stim_d7_2nd[sort(names(t_stim_d7_2nd))]

dat_trace_d7_2nd <- mapply(c_miniscope_matlab_ft_back, mouse_file, t_stim_d7_2nd, SIMPLIFY = F) 


dat_trajector_combine_d7_2nd <- mapply(cc_trajectory, dat_trace_d7_2nd, SIMPLIFY = F)


## for trajector distance
dat_trajector_dis <- lapply(dat_trajector_combine_d7_2nd, function(x) x[[2]]) %>% 
  do.call(rbind,.) %>% 
  as_tibble() %>% 
  rename(Cross = V1, Cross_back = V2) %>% 
  mutate(ID = sort(c("m3", "m7", "m17", "m18", "m855", "m857"))) %>% 
  pivot_longer(-ID) %>% 
  mutate(name = factor(name, levels = c("Cross","Cross_back")))

p_trajector_dis <- dat_trajector_dis %>% 
  ggplot(., aes(name, value, color = name))+
  geom_boxplot(outlier.shape = NA)+
  geom_jitter()+
  labs(x="PC1", y="PC2")+
  theme_classic()


t_trajectory <- dat_trajector_dis %>% 
  friedman_test(value ~ name | ID)

## for variance explained by pcs-----
dat_variance_pc <- lapply(dat_trajector_combine_d7_2nd, function(x) x[[3]]) %>% 
  do.call(rbind,.) %>% 
  as_tibble() %>% 
  mutate(Group = rep(c("Cross","Cross_back"), 6)) %>% 
  mutate(ID = rep(sort(c("m3", "m7", "m17", "m18", "m855", "m857")), each = 2)) %>% 
  mutate(PC1_2 = PC1 + PC2) %>% 
  mutate(Group = factor(Group, levels = c("Cross","Cross_back"))) %>% 
  ggplot(., aes(Group, PC1, color = Group))+
  geom_boxplot()


## trajectory similarity analysis-----
c_miniscope_matlab_ft <- function(file_trace) {
  ## import and format the data
  ID <- str_extract(file_trace, regex("m\\d+"))
  
  dat_trace1 <- raveio::read_mat(file_trace)
  num_compare <- c(4,7, 8)
  t_crossing <- dat_trace1$crossing[c(1,4, 5)]
  
  cross_ID <- dat_trace1$global_map %>% 
    as_tibble() %>% 
    dplyr::select(V1, V4, V5) %>% 
    mutate_all(na_if, 0) %>% 
    drop_na()
  
  dat_stim_trace <- vector(mode = "list", length = length(num_compare))
  
  
  for (i in seq_along(num_compare)) {
    global_cell <- pull(cross_ID[,i])
    dat_trace <- dat_trace1[[num_compare[i]]] %>% 
      as_tibble() %>% 
      dplyr::select(all_of(global_cell)) %>% 
      apply(., 2, scale)
    
    
    colnames(dat_trace) <- str_c(ID,"Cell", 1: ncol(dat_trace))
    dat_stim_trace[[i]] <- dat_trace
  }
  
  time_crossing_frame <- c((t_crossing[1] - 40): (t_crossing[1] + 80-1), nrow(dat_stim_trace[[1]]) + ((t_crossing[2] -40):(t_crossing[2]  + 80-1)), (nrow(dat_stim_trace[[1]]) + nrow(dat_stim_trace[[2]])) + ((t_crossing[3] -40):(t_crossing[3]+ 80-1)))
  
  dat_total <- dat_stim_trace %>% 
    do.call(rbind, .) %>% 
    as_tibble() %>% 
    as.matrix()
  
  dat_total_eigen <- dat_total %>% 
    prcomp(., scale. = F, center = T) 
  
  
  compare_time <- 120
  
  dat_trajector <- dat_total_eigen$x %>% 
    as_tibble() %>% 
    dplyr::slice(all_of(time_crossing_frame)) %>% 
    mutate(Group = c(rep(c("Pre","Cond", "Post"), each = compare_time))) %>% 
    add_column(Time = rep(c(1:compare_time), 3))
  
  p_2d_trajector <- dat_trajector %>% 
    ggplot(., aes(x = PC1, y = PC2, group = Group,color = Group))+
    geom_line(size = 1)
  
  ## calculat the cosine similarity between two PCA space vectors
  
  dat_pre_pc <- dat_trajector %>% 
    filter(Group == "Pre")
  dat_cond_pc <- dat_trajector %>% 
    filter(Group == "Cond")
  dat_post_pc <- dat_trajector %>% 
    filter(Group == "Post")
  
  cosine_pre_cond <- rep(0, compare_time)
  
  for (i in c(1:compare_time)){
    x <- dat_pre_pc %>%
      dplyr::select(PC1, PC2) %>% 
      slice(i) %>% 
      unlist(., use.names=FALSE)
    
    y <- dat_cond_pc %>%
      dplyr::select(PC1, PC2) %>% 
      slice(i) %>% unlist(., use.names=FALSE)
    
    cosine_pre_cond[i] <- cosine(x, y)
  }
  
  cosine_pre_post <- rep(0, compare_time)
  
  for (i in c(1:compare_time)){
    x <- dat_pre_pc %>%
      dplyr::select(PC1, PC2) %>% 
      slice(i) %>% 
      unlist(., use.names=FALSE)
    
    y <- dat_post_pc %>%
      dplyr::select(PC1, PC2) %>% 
      slice(i) %>% unlist(., use.names=FALSE)
    
    cosine_pre_post[i] <- cosine(x, y)
  }
  
  cosine_cond_post <- rep(0, compare_time)
  
  for (i in c(1:compare_time)){
    x <- dat_cond_pc %>%
      dplyr::select(PC1, PC2) %>% 
      slice(i) %>% 
      unlist(., use.names=FALSE)
    
    y <- dat_post_pc %>%
      dplyr::select(PC1, PC2) %>% 
      slice(i) %>% unlist(., use.names=FALSE)
    
    cosine_cond_post[i] <- cosine(x, y)
  }
  
  dat_similarity <- tibble(value = c(cosine_pre_cond, cosine_pre_post, cosine_cond_post), Time= rep(seq(1:compare_time), 3)) %>% 
    add_column(Group = rep(c("Pre_cond", "Pre_post", "Cond_post"), each = compare_time)) %>% 
    ggplot(., aes(Time, value, col = Group))+
    geom_line()+
    geom_vline(xintercept = 40, col = "black")+
    labs(title = ID)
  
  similarity_mean <- tibble(value = c(cosine_pre_cond, cosine_pre_post, cosine_cond_post), Time= rep(seq(1:compare_time), 3)) %>% 
    add_column(Group = rep(c("Pre_cond", "Pre_post", "Cond_post"), each = compare_time)) %>% 
    ddply(., .(Group), summarise, mean = mean(value))
  
  return(similarity_mean)
}


mouse_file <- as.list(list.files("~cchen/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/FD-processed", full.names = T))

dat_cell_trace <- mapply(c_miniscope_matlab_ft, mouse_file, SIMPLIFY = F)

dat_similarity <- dat_cell_trace %>% 
  do.call(rbind,.) %>% 
  ggplot(., aes(Group, mean))+
  geom_boxplot()

