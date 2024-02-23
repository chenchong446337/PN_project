c_miniscope_matlab <- function(file_trace, stim_trace) {
 
  
  dat_trace1 <- raveio::read_mat(file_trace) %>% 
    .[[1]] %>% 
    as_tibble() %>% 
    apply(., 2, scale) %>% 
    as_tibble()
 
  dat_stim <- read.xlsx(stim_trace, sheet = 1) %>% 
    as_tibble() %>% 
    mutate(Type = as.factor(Type))
  
  stim_type <- levels(dat_stim$Type)
  
  stim_time<- seq(-0.95, 2, by=0.05)
  ## number of rows to be binned
  n <- 2 # 0.05*2=0.1
  
  dat_trace_stim <- vector(mode = "list", length = length(stim_type))
  
  
  for (i in seq_along(stim_type)) {
    stim_type1 <- stim_type[i]
    
    t_stim <- dat_stim %>% 
      filter(Type ==stim_type1) %>% 
      select(Frame) %>% 
      unlist() %>% 
      unname()*2
    
    if(stim_type1 == "Pin_left"){
      t_stim = t_stim + 17877
    } else {
      t_stim = t_stim
    }
    
    dat_trace_average <- vector(mode = "list", length = length(t_stim))
    for (j in seq_along(t_stim)) {
      t1_p <- (t_stim[j]-20+1):(t_stim[j]+40)
      
      dat_stim1 <- dat_trace1 %>% 
        slice(t1_p)
      
      #dat_stim1 <- aggregate(dat_stim1,list(rep(1:(nrow(dat_stim1)%/%n+1),each=n,len=nrow(dat_stim1))),mean)[-1]
      
      dat_trace_average[[j]] <- dat_stim1[1:20,] %>%  ## baseline as -2 to 0
        colMeans(., na.rm = T) %>% 
        sweep(dat_stim1, 2, ., FUN = "-") %>% 
        as_tibble() %>% 
        add_column(Time = stim_time,.before = "V1")
      
    }
    dat_trace_stim[[i]] <- dat_trace_average %>% 
      do.call(rbind,.) %>% 
      pivot_longer(-Time) %>% 
      ddply(.,.(Time, name), summarise, value1 = mean(value)) %>% 
      add_column(Group= stim_type[i])
    
    
  }

  
  return(dat_trace_stim)
  
}

file_trace <- "~cchen/Documents/neuroscience/Pn\ project/Data_analysis/CB_PC_pain/CB_PC_Pain_m625.mat"

stim_trace <- "~cchen/Documents/neuroscience/Pn\ project/Data_analysis/CB_PC_pain/CB_PC_pain.xlsx"


dat_pain_m625 <- c_miniscope_matlab(file_trace, stim_trace) %>% 
  do.call(rbind,.) %>% 
  as_tibble() %>% 
  mutate(Group = as.factor(Group)) 


p_pain_trace <-dat_pain_m625 %>% 
  ddply(., .(Time, Group), summarise,n=length(value1),value=mean(value1),sd=sd(value1),se=sd(value1)/sqrt(length(value1))) %>% 
  ggplot(., aes(Time, value, color= Group))+
  geom_line()

group_stim <- levels(dat_pain_m625$Group)

score_range <- dat_pain_m625 %>% 
  select(value1) %>% 
  range()

for (i in seq_along(group_stim)) {
  dat_plot_order <- dat_pain_m625 %>% 
    filter(Group == group_stim[3]) %>% 
    ddply(.,.(name), summarise, value = mean(value1)) %>% 
    arrange(value)
  
  p_heat <-  dat_pain_m625 %>% 
    filter(Group == group_stim[i]) %>% 
    mutate(name = factor(name, levels =   dat_plot_order$name)) %>% 
    ggplot(., aes(Time, name,fill= value1))+ 
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
  assign(str_c("p_heat_", group_stim[i]), p_heat)
  
}

p_heat <- plot_grid(p_heat_Non, p_heat_Hot, p_heat_Pin, p_heat_Pin_left, nrow  = 1)

## for ca2+ events---------
c_miniscope_matlab <- function(file_trace, stim_trace) {
  
  
  dat_trace1 <- read.xlsx(file_trace, colNames = F) %>% 
    as_tibble() %>% 
    apply(., 2, scale) %>% 
    as_tibble()
  
  dat_stim <- read.xlsx(stim_trace, sheet = 1) %>% 
    as_tibble() %>% 
    mutate(Type = as.factor(Type))
  
  stim_type <- levels(dat_stim$Type)
  
  stim_time<- seq(-0.95, 2, by=0.05)
  ## number of rows to be binned
  n <- 2 # 0.05*2=0.1
  
  dat_trace_stim <- vector(mode = "list", length = length(stim_type))
  
  
  for (i in seq_along(stim_type)) {
    stim_type1 <- stim_type[i]
    
    t_stim <- dat_stim %>% 
      filter(Type ==stim_type1) %>% 
      select(Frame) %>% 
      unlist() %>% 
      unname()*2
    
    if(stim_type1 == "Pin_left"){
      t_stim = t_stim + 17877
    } else {
      t_stim = t_stim
    }
    
    dat_trace_average <- vector(mode = "list", length = length(t_stim))
    for (j in seq_along(t_stim)) {
      t1_p <- (t_stim[j]-20+1):(t_stim[j]+40)
      
      dat_stim1 <- dat_trace1 %>% 
        slice(t1_p)
      
      #dat_stim1 <- aggregate(dat_stim1,list(rep(1:(nrow(dat_stim1)%/%n+1),each=n,len=nrow(dat_stim1))),mean)[-1]
      
      dat_trace_average[[j]] <- dat_stim1[1:20,] %>%  ## baseline as -2 to 0
        colMeans(., na.rm = T) %>% 
        sweep(dat_stim1, 2, ., FUN = "-") %>% 
        as_tibble() %>% 
        add_column(Time = stim_time,.before = "X1")
      
    }
    dat_trace_stim[[i]] <- dat_trace_average %>% 
      do.call(rbind,.) %>% 
      pivot_longer(-Time) %>% 
      ddply(.,.(Time, name), summarise, value1 = mean(value)) %>% 
      add_column(Group= stim_type[i])
    
    
  }
  
  
  return(dat_trace_stim)
  
}

file_trace <- "~cchen/Documents/neuroscience/Pn\ project/Data_analysis/CB_PC_pain/CB_PC_Pain_m625_spike.xlsx"

stim_trace <- "~cchen/Documents/neuroscience/Pn\ project/Data_analysis/CB_PC_pain/CB_PC_pain.xlsx"


dat_pain_m625 <- c_miniscope_matlab(file_trace, stim_trace) %>% 
  do.call(rbind,.) %>% 
  as_tibble() %>% 
  mutate(Group = as.factor(Group)) 


p_pain_trace <-dat_pain_m625 %>% 
  ddply(., .(Time, Group), summarise,n=length(value1),value=mean(value1),sd=sd(value1),se=sd(value1)/sqrt(length(value1))) %>% 
  ggplot(., aes(Time, value, color= Group))+
  geom_line()

group_stim <- levels(dat_pain_m625$Group)

score_range <- dat_pain_m625 %>% 
  select(value1) %>% 
  range()

for (i in seq_along(group_stim)) {
  dat_plot_order <- dat_pain_m625 %>% 
    filter(Group == group_stim[3]) %>% 
    ddply(.,.(name), summarise, value = mean(value1)) %>% 
    arrange(value)
  
  p_heat <-  dat_pain_m625 %>% 
    filter(Group == group_stim[i]) %>% 
    mutate(name = factor(name, levels =   dat_plot_order$name)) %>% 
    ggplot(., aes(Time, name,fill= value1))+ 
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
  assign(str_c("p_heat_", group_stim[i]), p_heat)
  
}

p_heat <- plot_grid(p_heat_Non, p_heat_Hot, p_heat_Pin, p_heat_Pin_left, nrow  = 1)
