## for PCA analysis Kushal suggested-----


cc_globalID_fun <- function(path_ID, path_trace, t_stim){
  Global_ID <- read.csv(path_ID, header = F) %>% 
    na_if(., 0) %>% 
    drop_na()
  
  
  ## read.xlsx with no header and remove the first col
  cc_read.xlsx<- function(x){
    dat <- read.xlsx(x, colNames = F) %>% 
      select(., -1)
    dat
  }
  
  ## do z score of the whole trace
  dat_trace <- mapply(cc_read.xlsx, path_trace,SIMPLIFY = F) %>% 
    mapply(function(x) apply(x, 2, scale), ., SIMPLIFY = F)
  
  stim_time<- seq(-2, 6.5, by=0.5)
  ## number of rows to be binned
  n <- 10 # 0.05*10=0.5
  
  ## create a 2 days list and put the activation of each cell in it
  cc_trace_pick <- function(dat_trace, cell_pick){
    dat_trace_pick <- t(dat_trace[cell_pick,])
    dat_trace_pick
  }
  
  dat_cell_trace_day <- mapply(cc_trace_pick, dat_trace, as.list(Global_ID))
  
  cc_trace_extract <- function(dat_cell_trace, t_stim_day){
    dat_stim <- vector(mode = "list", length = length(t_stim_day))
    for (i in seq_along(t_stim_day)){
      t1_p <- t_stim_day[i]
      dat_stim1 <- dat_cell_trace[(t1_p-40):(t1_p+140-1),]
      dat_stim1 <- aggregate(dat_stim1,list(rep(1:(nrow(dat_stim1)%/%n+1),each=n,len=nrow(dat_stim1))),mean)[-1]
      dat_stim[[i]] <- dat_stim1[1:5,] %>%  ## baseline as -2 to 0
        colMeans(., na.rm = T) %>% 
        sweep(dat_stim1, 2, ., FUN = "-")
    }
    ## average the trace by cell number
    if (length(t_stim_day)>1){
      dat_cell_trace_average <- aaply(laply(dat_stim, as.matrix), c(2, 3), function(x) mean(x, na.rm=T))
      
    } else {
      dat_cell_trace_average <- data.matrix(dat_stim[[1]])
      
    }
    return(dat_cell_trace_average)
  }
  dat_cell_trace_day_extract <- mapply(cc_trace_extract, dat_cell_trace_day, t_stim, SIMPLIFY = F)
  
  return(dat_cell_trace_day_extract)
}

global_ID_m3 <- "~cchen/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/matlab_analysis/m3/Global_id/m3_anti_pin_global_ID.csv"
path_trace_m3 <- as.list(list.files(path ="~cchen/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/matlab_analysis/m3/trace_days/", pattern = "*.xlsx", full.names = T ))[c(3, 6,7)]
t_stim_m3 <- list(t_stim_m3_d3 = c(1885),t_stim_m3_d6 = c(132), t_stim_m3_d7 = c(346))
dat_global_trace_m3 <- cc_globalID_fun(global_ID_m3, path_trace_m3, t_stim_m3)



global_ID_m857 <- "~cchen/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/matlab_analysis/m857/Global_ID/m857_global_ID.csv"
path_trace_m857 <- as.list(list.files(path ="~cchen/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/matlab_analysis/m857/trace_days", pattern = "*.xlsx", full.names = T ))[c(1, 4, 5)]
t_stim_m857 <- list(t_stim_m857_d5 = c(392) *2, t_stim_m857_d6 = c(723) *2, t_stim_m857_d7 = c(167)*2)
dat_global_trace_m857 <- cc_globalID_fun(global_ID_m857, path_trace_m857, t_stim_m857)




## for corss-day aligned cells---

c_miniscope_matlab_ft <- function(file_trace) {
  ## import and format the data
  dat_trace1 <- raveio::read_mat(file_trace)
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
    n <- 10 # 0.05*10=0.5
    t1_p <- dat_trace1$crossing[t_crossing[i]]
    ## extract cell activity when they cross the border
    #dat_stim1 <- dat_trace[(t1_p-40):(t1_p+140-1),] 
    dat_stim1 <- dat_trace[(t1_p-40):(t1_p+140-1),] 
    
    dat_stim1 <- aggregate(dat_stim1,list(rep(1:(nrow(dat_stim1)%/%n+1),each=n,len=nrow(dat_stim1))),mean)[-1]
    dat_stim <- dat_stim1[1:5,] %>%  ## baseline as -2 to 0
      colMeans(., na.rm = T) %>% 
      sweep(dat_stim1, 2, ., FUN = "-")
    
    colnames(dat_stim) <- str_c("Cell", 1: ncol(dat_stim))
    dat_stim_trace[[i]] <- dat_stim
  }
  return(dat_stim_trace)
}



library(factoextra)

dat_trace_m3 <- c_miniscope_matlab_ft("~cchen/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/FD-processed/m3.mat")

pca_m3_pre_t <- prcomp(t(dat_trace_m3[[1]]), scale = F, center = T)
summary(pca_m3_pre_t)
biplot(pca_m3_pre_t, scale = 0)
pca_m3_cond_t <- prcomp(dat_trace_m3[[2]], scale = F, center = T)
summary(pca_m3_cond_t)

pca_m3_post_t <- prcomp(dat_trace_m3[[3]], scale = F, center = T)
summary(pca_m3_post_t)

dat_pca <- tibble(value = c(pca_m3_pre_t$center, pca_m3_cond_t$center,pca_m3_post_t$center), Time = rep(c(1:18), 3),Group = c(rep(c("Pre","Cond", "Post"), each = 18))) %>% 
  ggplot(., aes(Time, value, color = Group))+
  geom_line()

## trajector plot-------
library(plotly)
dat_m3_total <- mapply(function(x) prcomp(x,scale = F, center = T), dat_trace_m3, SIMPLIFY = F) %>% 
  mapply(function(x) x[[5]],., SIMPLIFY = F) %>% 
  do.call(rbind, .) %>% 
  as_tibble()  %>% 
  mutate(Group = c(rep(c("Pre","Cond", "Post"), each = 18))) %>% 
  mutate(Time = rep(c(1:18), 3)) %>% 
  mutate(Time = rep(seq(-2, 6.5, by=0.5), 3))



plot_ly(dat_m3_total, x = ~PC1, y = ~PC2, z = ~PC3, type = 'scatter3d', mode = 'lines', color = ~Group)

 ggplot(dat_m3_total, aes(Time,PC1, color = Group))+
  geom_line()


## based on the nature paper-----

 dat_m3_total <- dat_trace_m3 %>% 
   do.call(rbind, .) %>% 
   as.matrix()

 dat_m3_total_eigen <- dat_m3_total %>% 
   prcomp(., scale. = F, center = T) %>% 
   .$rotation 
 
 dat_m3_trajector <- dat_m3_total %*% dat_m3_total_eigen %>% 
   as_tibble() %>% 
   mutate(Group = c(rep(c("Pre","Cond", "Post"), each = 18))) %>% 
   mutate(Time = rep(seq(-2, 6.5, by=0.5), 3))
 
   
 dat_m3_trajector_start <- dat_m3_trajector %>% 
   filter(Time == -2.0 | Time == 0 | Time == 6.5) %>% 
   mutate(Catlog = rep(c("Begin","Crossing", "End"), 3 ))
 
 
 p_2d_trajector <- dat_m3_trajector %>% 
   ggplot(., aes(x = PC1, y = PC2, group = Group,color = Group))+
   geom_line(size = 1)+
   geom_point(data = dat_m3_trajector_start,aes(x=PC1, y = PC2, group = Group,  color = Group, shape = Catlog, size = 1) )
   
dat_dis <- rep(0, 3) 
dat_group <- c("Pre", "Cond", "Post")
for (i in seq_along(dat_group)) {
  dat <- dat_m3_trajector %>% 
    filter(Group ==dat_group[i]) %>% 
    select(PC1, PC2) %>% 
    dist(., diag = T, upper = T) %>% 
    as.matrix() %>% 
    .[,1] %>% 
    sum
  dat_dis[i] <- dat
}

 
 
 plot_ly(dat_m3_trajector, x = ~PC1, y = ~PC2, z = ~PC3, type = 'scatter3d', mode = 'lines', color = ~Group)
 
 
 ## for m7
 dat_m7_total <- dat_trace_m7 %>% 
   do.call(rbind, .) %>% 
   as.matrix()
 
 dat_m7_total_eigen <- dat_m7_total %>% 
   prcomp(., scale. = F, center = T) %>% 
   .$rotation 
 
 dat_m7_trajector <- dat_m7_total %*% dat_m7_total_eigen %>% 
   as_tibble() %>% 
   mutate(Group = c(rep(c("Pre","Cond", "Post"), each = 18))) %>% 
   mutate(Time = rep(seq(-2, 6.5, by=0.5), 3))
 
 
 dat_m7_trajector_start <- dat_m7_trajector %>% 
   filter(Time == -2.0 | Time == 0 | Time == 6.5) %>% 
   mutate(Catlog = rep(c("Begin","Crossing", "End"), 3 ))
 
 
 p_2d_trajector <- dat_m7_trajector %>% 
   ggplot(., aes(x = PC1, y = PC2, group = Group,color = Group))+
   geom_line(size = 1)+
   geom_point(data = dat_m7_trajector_start,aes(x=PC1, y = PC2, group = Group,  color = Group, shape = Catlog, size = 1) )
 
 
 plot_ly(dat_m7_trajector, x = ~PC1, y = ~PC2, z = ~PC3, type = 'scatter3d', mode = 'lines', color = ~Group)
 


dat_trace_m7 <- c_miniscope_matlab_ft("~cchen/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/FD-processed/m7.mat")

pca_m7_pre_t <- prcomp(dat_trace_m7[[1]], scale = F, center = T)
summary(pca_m7_pre_t)
pca_m7_cond_t <- prcomp(t(dat_trace_m7[[2]]), scale = F, center = T)
summary(pca_m7_cond_t)

pca_m7_post_t <- prcomp(t(dat_trace_m7[[3]]), scale = F, center = T)
summary(pca_m7_post_t)


dat_m7_total <- mapply(function(x) prcomp(x,scale = F, center = F), dat_trace_m7, SIMPLIFY = F) %>% 
  mapply(function(x) x[[5]],., SIMPLIFY = F) %>% 
  do.call(rbind, .) %>% 
  as_tibble()  %>% 
  mutate(Group = c(rep(c("Pre","Cond", "Post"), each = 18))) %>% 
  mutate(Time = rep(c(1:18), 3)) %>% 
  mutate(Time = rep(seq(-2, 6.5, by=0.5), 3))



plot_ly(dat_m7_total, x = ~PC1, y = ~PC2, z = ~PC3, type = 'scatter3d', mode = 'lines', color = ~Group)

ggplot(dat_m7_total, aes(Time,PC1, color = Group))+
  geom_line()



dat_trace_m17 <- c_miniscope_matlab_ft("~cchen/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/FD-processed/m17.mat")

pca_m17_pre_t <- prcomp(t(dat_trace_m17[[1]]), scale = F)

pca_m17_cond_t <- prcomp(t(dat_trace_m17[[2]]), scale = F)

pca_m17_post_t <- prcomp(t(dat_trace_m17[[3]]), scale = F)

dat_m17_total <- mapply(function(x) prcomp(x,scale = F, center = F), dat_trace_m17, SIMPLIFY = F) %>% 
  mapply(function(x) x[[5]],., SIMPLIFY = F) %>% 
  do.call(rbind, .) %>% 
  as_tibble()  %>% 
  mutate(Group = c(rep(c("Pre","Cond", "Post"), each = 18))) %>% 
  mutate(Time = rep(c(1:18), 3)) %>% 
  mutate(Time = rep(seq(-2, 6.5, by=0.5), 3))



plot_ly(dat_m17_total, x = ~PC1, y = ~PC2, z = ~PC3, type = 'scatter3d', mode = 'lines', color = ~Group)

ggplot(dat_m17_total, aes(Time,PC1, color = Group))+
  geom_line()


dat_trace_m857 <- c_miniscope_matlab_ft("~cchen/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/FD-processed/m857.mat")

pca_m857_pre_t <- prcomp(t(dat_trace_m857[[1]]), scale = F)

pca_m857_cond_t <- prcomp(t(dat_trace_m857[[2]]), scale = F)

pca_m857_post_t <- prcomp(t(dat_trace_m857[[3]]), scale = F)


dat_m857_total <- mapply(function(x) prcomp(x,scale = F, center = F), dat_trace_m857, SIMPLIFY = F) %>% 
  mapply(function(x) x[[5]],., SIMPLIFY = F) %>% 
  do.call(rbind, .) %>% 
  as_tibble()  %>% 
  mutate(Group = c(rep(c("Pre","Cond", "Post"), each = 18))) %>% 
  mutate(Time = rep(c(1:18), 3)) %>% 
  mutate(Time = rep(seq(-2, 6.5, by=0.5), 3))



plot_ly(dat_m857_total, x = ~PC1, y = ~PC2, z = ~PC3, type = 'scatter3d', mode = 'lines', color = ~Group)

ggplot(dat_m857_total, aes(Time,PC1, color = Group))+
  geom_line()

## m855
dat_trace_m855 <- c_miniscope_matlab_ft("~cchen/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/FD-processed/m855.mat")

pca_m855_pre_t <- prcomp(t(dat_trace_m855[[1]]), scale = F)

pca_m855_cond_t <- prcomp(t(dat_trace_m855[[2]]), scale = F)

pca_m855_post_t <- prcomp(t(dat_trace_m855[[3]]), scale = F)


dat_m855_total <- mapply(function(x) prcomp(x,scale = F, center = F), dat_trace_m855, SIMPLIFY = F) %>% 
  mapply(function(x) x[[5]],., SIMPLIFY = F) %>% 
  do.call(rbind, .) %>% 
  as_tibble()  %>% 
  mutate(Group = c(rep(c("Pre","Cond", "Post"), each = 18))) %>% 
  mutate(Time = rep(c(1:18), 3)) %>% 
  mutate(Time = rep(seq(-2, 6.5, by=0.5), 3))



plot_ly(dat_m855_total, x = ~PC1, y = ~PC2, z = ~PC3, type = 'scatter3d', mode = 'lines', color = ~Group)

ggplot(dat_m855_total, aes(Time,PC1, color = Group))+
  geom_line()



pca_m3_pre_t <- prcomp(t(dat_global_trace_m3[[1]]), scale = F)

pca_m3_cond_t <- prcomp(t(dat_global_trace_m3[[2]]), scale = F)

pca_m3_post_t <- prcomp(t(dat_global_trace_m3[[3]]), scale = F)


dat_pca <- tibble(value = c(pca_m3_pre_t$center, pca_m3_cond_t$center,pca_m3_post_t$center), Time = rep(c(1:18), 3),Group = c(rep(c("Pre","Cond", "Post"), each = 18))) %>% 
  ggplot(., aes(Time, value, color = Group))+
  geom_line()


pca_m7_pre <- prcomp(dat_global_trace_m7[[1]], scale = F)

pca_m7_pre_t <- prcomp(t(dat_global_trace_m7[[1]]), scale = F)


pca_m7_post_t <- prcomp(t(dat_global_trace_m7[[2]]), scale = F)


dat_pca <- tibble(value = c(pca_m7_pre_t$center, pca_m7_post_t$center), Time = rep(c(1:18), 2),Group = c(rep(c("Pre", "Post"), each = 18))) %>% 
  ggplot(., aes(Time, value, color = Group))+
  geom_line()


pca_m17_pre <- prcomp(dat_global_trace_m17[[1]], scale = F)

pca_m17_pre_t <- prcomp(t(dat_global_trace_m17[[1]]), scale = F)


pca_m17_post_t <- prcomp(t(dat_global_trace_m17[[2]]), scale = F)


dat_pca <- tibble(value = c(pca_m17_pre_t$center, pca_m17_post_t$center), Time = rep(c(1:18), 2),Group = c(rep(c("Pre", "Post"), each = 18))) %>% 
  ggplot(., aes(Time, value, color = Group))+
  geom_line()

pca_m18_pre_t <- prcomp(t(dat_global_trace_m18[[1]]), scale = F)


pca_m18_post_t <- prcomp(t(dat_global_trace_m18[[2]]), scale = F)


dat_pca <- tibble(value = c(pca_m18_pre_t$center, pca_m18_post_t$center), Time = rep(c(1:18), 2),Group = c(rep(c("Pre", "Post"), each = 18))) %>% 
  ggplot(., aes(Time, value, color = Group))+
  geom_line()


pca_m855_pre_t <- prcomp(t(dat_global_trace_m855[[1]]), scale = F)


pca_m855_post_t <- prcomp(t(dat_global_trace_m855[[2]]), scale = F)


dat_pca <- tibble(value = c(pca_m855_pre_t$center, pca_m855_post_t$center), Time = rep(c(1:18), 2),Group = c(rep(c("Pre", "Post"), each = 18))) %>% 
  ggplot(., aes(Time, value, color = Group))+
  geom_line()


pca_m857_pre_t <- prcomp(t(dat_global_trace_m857[[1]]), scale = F)


pca_m857_post_t <- prcomp(t(dat_global_trace_m857[[2]]), scale = F)


dat_pca <- tibble(value = c(pca_m857_pre_t$center, pca_m857_post_t$center), Time = rep(c(1:18), 2),Group = c(rep(c("Pre", "Post"), each = 18))) %>% 
  ggplot(., aes(Time, value, color = Group))+
  geom_line()

