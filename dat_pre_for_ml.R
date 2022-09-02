

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
    time_window <- (t_crossing[i]-40):(t_crossing[i]+100-1)
    
    dat_trace <- dat_trace1[[num_compare[i]]] %>% 
      as_tibble() %>% 
      dplyr::select(all_of(global_cell)) %>% 
      apply(., 2, scale) %>% 
      as_tibble() %>% 
      slice(all_of(time_window))
    
    
    colnames(dat_trace) <- str_c(ID,"Cell", 1: ncol(dat_trace))
    dat_stim_trace[[i]] <- dat_trace
  }
  
  dat_stim_trace_pc <- dat_stim_trace %>% 
    do.call(rbind, .) %>% 
    prcomp(., scale. = F, center = T) %>% 
    .$x %>% 
    as_tibble() %>% 
    dplyr::select(c(1:10)) %>% 
    split(.,rep(1:length(num_compare),each=140))
    
    
  return(dat_stim_trace_pc)
}

mouse_file <- as.list(list.files("~cchen/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/FD-processed", full.names = T))

dat_cell_trace <- mapply(c_miniscope_matlab_ft, mouse_file, SIMPLIFY = F)


## pre data as non-placebo and cond data as placebo
library(abind)
dat_pre <- lapply(dat_cell_trace, function (x) x[[1]]) %>% 
  unlist() %>% 
  as.numeric() %>% 
  array(., dim=c(140, 10, 6))

dat_train_array <- lapply(dat_cell_trace, function (x) x[[2]]) %>% 
  unlist() %>% 
  as.numeric() %>% 
  array(., dim=c(140, 10, 6)) %>% 
  abind(dat_pre, ., along = 3) %>% 
  aperm(., c(3,1,2))

## 0 for non-placebo, 1 for placebo
dat_train_array_label <- matrix(rep(c(0,1), each = 6))

dat_test <- lapply(dat_cell_trace, function (x) x[[3]]) %>% 
  unlist() %>% 
  as.numeric() %>% 
  array(., dim=c(140, 10, 6)) %>% 
  aperm(., c(3,1,2))

dat_test_y <- matrix(rep(1, 6))


dat_placebo <- list(dat_train = dat_train_array, dat_train_y = dat_train_array_label, dat_test = dat_test, dat_test_y= dat_test_y)


library(jsonlite)
write_json(dat_placebo, "dat_placebo.json")

## use PCA to reduce the dimension-------

crossing_info <- read.xlsx("~cchen/Documents/neuroscience/Pn\ project/Pn_project/Frame_crossing_eachday.xlsx", sheet = 1) %>% 
  as_tibble()

c_miniscope_matlab_ft <- function(file_trace) {
  ## import and format the data
  ID_mouse <- str_extract(file_trace, regex("m\\d+"))
  
  dat_trace1 <- raveio::read_mat(file_trace) 
  
  dat_trace_filter <- dat_trace1 %>% 
    .[c(4, 6, 7, 8)]

  ## nrow of each day
  day_row <- mapply(nrow, dat_trace_filter) %>% 
    unlist()

  cross_ID <- dat_trace1$global_map %>% 
    as_tibble() %>% 
    dplyr::select(V1, V3, V4, V5) %>% 
    mutate_all(na_if, 0) %>% 
    drop_na()
  
  ## select cross day aligned neurons
  for (i in 1: length(dat_trace_filter)) {
    global_cell <- pull(cross_ID[,i])
    
    dat <- dat_trace_filter[[i]] %>% 
      as_tibble() %>% 
      dplyr::select(all_of(global_cell))
    
    colnames(dat) <- str_c(ID_mouse,"Cell", 1: ncol(dat))
    
    dat_trace_filter[[i]] <- dat
    
  }

  ## do PCA analysis
  dat_trace_filter_pca <- dat_trace_filter %>% 
    do.call(rbind,.) %>% 
    prcomp() %>% 
    .$x  %>% 
    as_tibble() %>% 
    dplyr::select(c(1:10)) %>% 
    split(.,rep(1:length(dat_trace_filter), day_row))
  
  t_crossing_ID <- crossing_info %>% 
    filter(ID == ID_mouse) %>% 
    mutate(Day = as.factor(Day))
  
  day_crossing <- levels(t_crossing_ID$Day)
  
  dat_stim_trace <- vector(mode = "list", length = length(dat_trace_filter_pca))
  
  t_compare <- 40 # 2s *20 Hz
  
  for (j in 1: length(dat_trace_filter_pca)) {
    
    t_crossing <- t_crossing_ID %>% 
      dplyr::filter(Day == day_crossing[j]) %>% 
      select(Frame) %>% 
      unlist()
    
    dat_list <- vector(mode = "list", length = length(t_crossing))
    
    if (length(t_crossing) == 1){
      time_window <- (t_crossing-40):(t_crossing+t_compare-1)
      
      dat_trace <- dat_trace_filter_pca[[j]] %>% 
        slice(all_of(time_window)) 
      
      dat_list[[1]] <- dat_trace
      
    } else {
      for (k in seq_along(t_crossing)) {
        time_window <- (t_crossing[k]-40):(t_crossing[k]+t_compare-1)
        
        dat_trace <- dat_trace_filter_pca[[j]] %>% 
          slice(all_of(time_window)) 
        
        dat_list[[k]] <- dat_trace
      }
    }
    
    dat_stim_trace[[j]] <- dat_list

  }
  return(dat_stim_trace)
}


mouse_file <- as.list(list.files("~cchen/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/FD-processed", full.names = T))

dat_cell_trace <- mapply(c_miniscope_matlab_ft, mouse_file, SIMPLIFY = F)


## pre data as non-placebo and cond data as placebo
t_compare_length <- 80
library(abind)

length_pre <- lapply(dat_cell_trace, function (x) x[[1]]) %>% 
  unlist(., recursive = FALSE) %>% 
  length()

dat_pre <- lapply(dat_cell_trace, function (x) x[[1]]) %>% 
  unlist(., recursive = FALSE) %>% 
  unlist() %>% 
  as.numeric() %>% 
  array(., dim=c(t_compare_length, 10, length_pre))

## with day5
length_placebo <- lapply(dat_cell_trace, function (x) x[[2]]) %>% 
  c(., lapply(dat_cell_trace, function (x) x[[3]])) %>% 
  unlist(., recursive = FALSE) %>% 
  length()

dat_train_array <- lapply(dat_cell_trace, function (x) x[[2]]) %>% 
  c(., lapply(dat_cell_trace, function (x) x[[3]])) %>% 
  unlist(., recursive = FALSE) %>% 
  unlist() %>% 
  as.numeric() %>% 
  array(., dim=c(t_compare_length, 10, length_placebo)) %>% 
  abind(dat_pre, ., along = 3) %>% 
  aperm(., c(3,1,2))


## 0 for non-placebo, 1 for placebo
dat_train_array_label <- matrix(c(rep(0, length_pre), rep(1, length_placebo)))

dat_test <- lapply(dat_cell_trace, function (x) x[[4]]) %>% 
  unlist() %>% 
  as.numeric() %>% 
  array(., dim=c(t_compare_length, 10, 6)) %>% 
  aperm(., c(3,1,2))

dat_test_y <- matrix(rep(1, 6))


dat_placebo <- list(dat_train = dat_train_array, dat_train_y = dat_train_array_label, dat_test = dat_test, dat_test_y= dat_test_y)


library(jsonlite)
write_json(dat_placebo, "dat_placebo.json")

