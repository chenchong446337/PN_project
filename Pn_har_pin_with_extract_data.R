## function for analysis-----
c_miniscope_matlab <- function(file_trace) {
  ## import and format the data
  ID <- str_extract(file_trace, regex("m\\d+"))
  t_stim <- stim_t_list[[ID]] %>% 
    unlist() %>% 
    unname()
  
  dat_trace1 <- raveio::read_mat(file_trace) %>% 
    .[[1]] %>% 
    as_tibble() %>% 
    mutate(across(.fns = scale))
  
  stim_time<- seq(-0.95, 2, by=0.05)
  ## number of rows to be binned
  n <- 2 # 0.05*2=0.1
  
  dat_trace_average <- vector(mode = "list", length = length(t_stim))
  for (i in seq_along(t_stim)) {
    t1_p <- (t_stim[i]-20+1):(t_stim[i]+40)
    
    dat_stim1 <- dat_trace1 %>% 
      slice(t1_p)
    
    #dat_stim1 <- aggregate(dat_stim1,list(rep(1:(nrow(dat_stim1)%/%n+1),each=n,len=nrow(dat_stim1))),mean)[-1]
    
    dat_trace_average[[i]] <- dat_stim1[1:20,] %>%  ## baseline as -2 to 0
      colMeans(., na.rm = T) %>% 
      sweep(dat_stim1, 2, ., FUN = "-") %>% 
      as_tibble() %>% 
      add_column(Time = stim_time,.before = "V1")
    
  }
  
  dat_trace_average1 <- dat_trace_average %>% 
    do.call(rbind,.) %>% 
    pivot_longer(-Time) %>% 
    ddply(.,.(Time, name), summarise, value1 = mean(value)) %>% 
    add_column(ID = ID)
  
  return(dat_trace_average1)
  
}
## for pinprick-----
file_trace_pin <- as.list(list.files("~cchen/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/Pin_extract/pinprick/Pin/", full.names = T))

stim_t_list <- list(m3 = c(337,665,916,1303,1541,1767,2032,2277,2496,2824),
                    m7 = c(493, 1323, 1867, 2415, 3035, 3556, 4079, 4519, 5020, 5517) ,
                    m16 = c(115, 414, 641,965, 1272, 1597, 1879, 2130, 2438,2679)*2,
                    m17 = c(262, 497, 803, 1057, 1404, 1676, 2020, 2248, 2543, 2805 )*2,
                    m18 = c(139, 681, 1012, 1365, 1663, 1900, 2110, 2324, 2566, 2769)*2)

dat_pin_trace <- mapply(c_miniscope_matlab, file_trace_pin, SIMPLIFY = F) %>% 
  do.call(rbind,.) %>% 
  add_column(Group = "Pin", .before = "Time")
rm(stim_t_list)
rm(file_trace_pin)

## for light touch as control------
file_trace_pin_ctrl <- as.list(list.files("~cchen/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/Pin_extract/pinprick/ctrl/", full.names = T))

stim_t_list <- list(m3 = c(134, 528, 820, 1084, 1261, 1493, 1729, 1964, 2142,2367)*2,
                    m7 = c(195,594, 901, 1163, 1536, 1818, 2137, 2450, 2624, 3008 )*2 ,
                    m16 = c(192, 498, 779, 1068, 1359, 1624, 1918, 2209, 2592, 2858)*2 ,
                    m17 = c(216, 605, 982, 1198, 1604, 1829, 2098, 2348, 2520, 2743)*2,
                    m18 = c(225, 465, 883, 1385, 1561, 1735, 1897, 2085, 2370, 2691 ) *2)

dat_pin_ctrl_trace <- mapply(c_miniscope_matlab, file_trace_pin_ctrl, SIMPLIFY = F) %>% 
  do.call(rbind,.) %>% 
  add_column(Group = "Touch", .before = "Time")
rm(stim_t_list)
rm(file_trace_pin_ctrl)


## for har------
file_trace_har <- as.list(list.files("~cchen/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/Pin_extract/Har/", pattern = ".mat",full.names = T))

stim_t_list <- read.xlsx("~cchen/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/har_06112020/2020_06_01_har_events.xlsx") %>% 
  select(Mice, Frame) %>% 
  group_by(Mice) %>% 
  group_split(., .keep = F) 
names(stim_t_list) <- sort(c("m3", "m7", "m16", "m17", "m18"))
  


dat_har_trace <- mapply(c_miniscope_matlab, file_trace_har, SIMPLIFY = F) %>% 
  do.call(rbind,.) %>% 
  add_column(Group = "Har", .before = "Time")
rm(stim_t_list)
rm(file_trace_har)

## for pan expectation----
file_trace_exp_pain <- as.list(list.files("~cchen/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/Pin_extract/Expectation_pain/", pattern = ".mat",full.names = T))

dat_stim_type <- read.csv("~cchen/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/pinprick_pain_exp/random_pinprick.csv", row.names = 1)
stim_t_list <- read.xlsx("~cchen/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/pinprick_pain_exp/pin_event.xlsx") %>% 
  .[dat_stim_type =="Pain"] %>% 
  matrix(., nrow = 10) %>%
  as_tibble()%>% 
  mutate(across(.fns = ~ .x *2)) %>% 
  map(unlist)

names(stim_t_list) <- c("m3", "m7", "m16", "m17", "m18")



dat_exp_pain_trace <- mapply(c_miniscope_matlab, file_trace_exp_pain, SIMPLIFY = F) %>% 
  do.call(rbind,.) %>% 
  add_column(Group = "Exp_pain", .before = "Time")
rm(stim_t_list)

## for expect touch
stim_t_list <- read.xlsx("~cchen/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/pinprick_pain_exp/pin_event.xlsx") %>% 
  .[dat_stim_type !="Pain"] %>% 
  matrix(., nrow = 10) %>%
  as_tibble()%>% 
  mutate(across(.fns = ~ .x *2)) %>% 
  map(unlist)

names(stim_t_list) <- c("m3", "m7", "m16", "m17", "m18")



dat_exp_touch_trace <- mapply(c_miniscope_matlab, file_trace_exp_pain, SIMPLIFY = F) %>% 
  do.call(rbind,.) %>% 
  add_column(Group = "Exp_touch", .before = "Time")
rm(stim_t_list)
rm(file_trace_exp_pain)


## combine data together for heatmap

dat_pain_combine <- bind_rows(dat_pin_trace, dat_pin_ctrl_trace, dat_har_trace, dat_exp_touch_trace, dat_exp_pain_trace)
rm(dat_pin_trace, dat_pin_ctrl_trace, dat_har_trace, dat_exp_touch_trace, dat_exp_pain_trace)

## heatmap plot
score_range <- dat_pain_combine %>% 
  filter(Group != "Exp_pain") %>% 
  filter(Group !="Exp_touch") %>% 
  select(value1) %>% 
  range()

group_stim <- c("Touch", "Pin", "Har")

for (i in seq_along(group_stim)) {
  dat_plot_order <- dat_pain_combine %>% 
    filter(Group == group_stim[i]) %>% 
    mutate(order = str_c(ID, name)) %>% 
    ddply(.,.(order), summarise, value = mean(value1)) %>% 
    arrange(value)
  
  p_heat <-  dat_pain_combine %>% 
    filter(Group == group_stim[i]) %>% 
    mutate(order = str_c(ID, name)) %>% 
    mutate(order = factor(order, levels =   dat_plot_order$order)) %>% 
    ggplot(., aes(Time, order,fill= value1))+ 
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

p_heat <- plot_grid(p_heat_Touch, p_heat_Pin, p_heat_Har, nrow  = 1)
rm(p_heat_Touch, p_heat_Pin, p_heat_Har, p_heat_Exp_touch, p_heat_Exp_pain)

setwd("~cchen/Documents/neuroscience/Pn\ project/Figure/PDF/")
cairo_pdf("p_heat.pdf", width = 120/25.6, height = 90/25.6, family = "Arial")
p_heat
dev.off()


## plot the traces
p_com_trace <- dat_pain_combine %>% 
  filter(Group != "Exp_pain") %>% 
  filter(Group !="Exp_touch") %>% 
  ddply(., .(Group, Time, ID), summarise, value = mean(value1)) %>% 
  ddply(., .(Group, Time), summarise,mean=mean(value),sd=sd(value),se=sd(value)/sqrt(length(value))) %>%  
  mutate(Group = factor(Group, levels = c("Touch", "Pin","Har"))) %>% 
  ggplot(., aes(Time, mean, color = Group))+
  geom_line()+
  geom_ribbon(aes(ymin=mean-se, ymax=mean+se, fill=Group), alpha=0.1, linetype=0)+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  scale_y_continuous(limits = c(-0.1, 0.1)) +
  theme(legend.title = element_blank())

setwd("~cchen/Documents/neuroscience/Pn\ project/Figure/PDF/")
cairo_pdf("p_com_trace.pdf", width = 110/25.6, height = 65/25.6, family = "Arial")
p_com_trace
dev.off()

## compare pin and expect pain-----
p_com_trace <- dat_pain_combine %>% 
  filter(Group == "Exp_pain"|Group =="Pin" ) %>% 
  ddply(., .(Group, Time, ID), summarise, value = mean(value1)) %>% 
  ddply(., .(Group, Time), summarise,mean=mean(value),sd=sd(value),se=sd(value)/sqrt(length(value))) %>%  
  mutate(Group = factor(Group, levels = c("Pin","Exp_pain"))) %>% 
  ggplot(., aes(Time, mean, color = Group))+
  geom_line()+
  geom_ribbon(aes(ymin=mean-se, ymax=mean+se, fill=Group), alpha=0.1, linetype=0)+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  scale_y_continuous(limits = c(-0.1, 0.1)) +
  theme(legend.title = element_blank())
## compare the cell activity during expect touch and expect pain------

score_range <- dat_pain_combine %>% 
  filter(Group == "Exp_pain"| Group =="Exp_touch") %>% 
  select(value1) %>% 
  range()
group_stim <- c("Exp_pain", "Exp_touch")

for (i in seq_along(group_stim)) {
  dat_plot_order <- dat_pain_combine %>% 
    filter(Group == group_stim[1]) %>% 
    mutate(order = str_c(ID, name)) %>% 
    ddply(.,.(order), summarise, value = mean(value1)) %>% 
    arrange(value)
  
  p_heat <-  dat_pain_combine %>% 
    filter(Group == group_stim[i]) %>% 
    mutate(order = str_c(ID, name)) %>% 
    mutate(order = factor(order, levels =   dat_plot_order$order)) %>% 
    ggplot(., aes(Time, order,fill= value1))+ 
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

p_exp_com <- plot_grid(p_heat_Exp_pain, p_heat_Exp_touch, nrow = 1)



### for cross-stim-aligned neurons-------

# 1. arrange the data
## trace order: Expect, har, pin, touch

## time of stimulation
stim_t_list_pin <- list(m3 = c(337,665,916,1303,1541,1767,2032,2277,2496,2824),
                        m7 = c(493, 1323, 1867, 2415, 3035, 3556, 4079, 4519, 5020, 5517) ,
                        m16 = c(115, 414, 641,965, 1272, 1597, 1879, 2130, 2438,2679)*2,
                        m17 = c(262, 497, 803, 1057, 1404, 1676, 2020, 2248, 2543, 2805 )*2,
                        m18 = c(139, 681, 1012, 1365, 1663, 1900, 2110, 2324, 2566, 2769)*2)

stim_t_list_touch <- list(m3 = c(134, 528, 820, 1084, 1261, 1493, 1729, 1964, 2142,2367)*2,
                          m7 = c(195,594, 901, 1163, 1536, 1818, 2137, 2450, 2624, 3008 )*2 ,
                          m16 = c(192, 498, 779, 1068, 1359, 1624, 1918, 2209, 2592, 2858)*2 ,
                          m17 = c(216, 605, 982, 1198, 1604, 1829, 2098, 2348, 2520, 2743)*2,
                          m18 = c(225, 465, 883, 1385, 1561, 1735, 1897, 2085, 2370, 2691 ) *2)
stim_t_list_har <- read.xlsx("~cchen/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/har_06112020/2020_06_01_har_events.xlsx") %>% 
  select(Mice, Frame) %>% 
  group_by(Mice) %>% 
  group_split(., .keep = F) 
names(stim_t_list_har) <- sort(c("m3", "m7", "m16", "m17", "m18"))

dat_stim_type <- read.csv("~cchen/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/pinprick_pain_exp/random_pinprick.csv", row.names = 1)
stim_t_list_expect_pain <- read.xlsx("~cchen/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/pinprick_pain_exp/pin_event.xlsx") %>% 
  .[dat_stim_type =="Pain"] %>% 
  matrix(., nrow = 10) %>%
  as_tibble()%>% 
  mutate(across(.fns = ~ .x *2)) %>% 
  map(unlist)

names(stim_t_list_expect_pain) <- c("m3", "m7", "m16", "m17", "m18")

stim_t_list_expect_touch <- read.xlsx("~cchen/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/pinprick_pain_exp/pin_event.xlsx") %>% 
  .[dat_stim_type !="Pain"] %>% 
  matrix(., nrow = 10) %>%
  as_tibble()%>% 
  mutate(across(.fns = ~ .x *2)) %>% 
  map(unlist)

names(stim_t_list_expect_touch) <- c("m3", "m7", "m16", "m17", "m18")

c_extract_align <- function(file_trace){
  ID <- str_extract(file_trace, regex("m\\d+"))
  
  dat_trace <-   raveio::read_mat(file_trace)
  
  global_ID <- dat_trace[[1]] %>% 
    as_tibble() %>% 
    mutate_all(na_if, 0) %>% 
    drop_na()
  
  t_expect_touch <- stim_t_list_expect_touch[ID] %>% 
    unlist() %>% 
    unname()
  
  t_expect_pain <- stim_t_list_expect_pain[ID] %>% 
  unlist() %>% 
  unname()
  
  t_har <- stim_t_list_har[ID] %>% 
    unlist() %>% 
    unname()
  
  t_pin <- stim_t_list_pin[ID] %>% 
    unlist() %>% 
    unname()
  
  t_touch <- stim_t_list_touch[ID] %>% 
    unlist() %>% 
    unname()
  
  t_stim <- list(t_expect_touch, t_expect_pain, t_har, t_pin, t_touch)
  
  trace_num <- c(1,1,2:4)
  stim_time<- seq(-0.95, 2, by=0.05)
  
  dat_stim_average <- vector(mode = "list", length = 5)

  for (i in 1:length(trace_num)) {
    global_cell <- pull(global_ID[,trace_num[i]])
    
    dat_stim_trace <- dat_trace[[2]][[trace_num[i]]] %>% 
      as_tibble() %>% 
      dplyr::select(all_of(global_cell)) %>% 
      mutate(across(.fns = scale))
    
    stim_t <- t_stim[[i]]
    
    dat_trace_t <- vector(mode = "list", length = length(stim_t))
    for (j in seq_along(stim_t)){
      t1_p <- (stim_t[j]-20+1):(stim_t[j]+40)
      
      dat_stim1 <- dat_stim_trace %>% 
        slice(t1_p)
      
      #dat_stim1 <- aggregate(dat_stim1,list(rep(1:(nrow(dat_stim1)%/%n+1),each=n,len=nrow(dat_stim1))),mean)[-1]
      
      dat_trace_t[[j]] <- dat_stim1[1:20,] %>%  ## baseline as -2 to 0
        colMeans(., na.rm = T) %>% 
        sweep(dat_stim1, 2, ., FUN = "-") %>% 
        as_tibble() %>% 
        add_column(Time = stim_time)
      
    }
    dat_stim_average[[i]] <- dat_trace_t %>% 
      do.call(rbind,.) %>% 
      pivot_longer(-Time) %>% 
      ddply(., .(Time, name), summarise, value1 = mean(value)) %>% 
      add_column(ID = ID)
    
  }
  
  return(dat_stim_average)
  }

path_trace_file <- as.list(list.files("~cchen/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/Pin_extract/stim_combine/", full.names = T))

dat_trace_combine <- mapply(c_extract_align, path_trace_file, SIMPLIFY = F)
rm(list=setdiff(ls(), "dat_trace_combine"))


## split by stim
stim <- c("Expect_touch", "Expect_pain", "Har", "Pin", "touch")

dat_trace_combine_re <- NULL
for (i in seq_along(stim)) {
  dat <- dat_trace_combine %>% 
    lapply(., function(x) x[[i]]) %>% 
    do.call(rbind,.) %>% 
    add_column(Group = stim[i])
  
dat_trace_combine_re <- rbind(dat_trace_combine_re, dat)  
  
}
rm(list=setdiff(ls(), "dat_trace_combine_re"))

## heatmap plot
score_range <- dat_trace_combine_re %>% 
  # filter(Group != "Exp_pain") %>% 
  # filter(Group !="Exp_touch") %>% 
  select(value1) %>% 
  range()

#group_stim <- c("Expect_touch", "Expect_pain", "Har", "Pin", "touch")
group_stim <- c( "Har", "Pin", "touch")
for (i in seq_along(group_stim)) {
  dat_plot_order <- dat_trace_combine_re %>% 
    filter(Group == group_stim[i]) %>% 
    mutate(order = str_c(ID, name)) %>% 
    ddply(.,.(order), summarise, value = mean(value1)) %>% 
    arrange(value)
  
  p_heat <-  dat_trace_combine_re %>% 
    filter(Group == group_stim[i]) %>% 
    mutate(order = str_c(ID, name)) %>% 
    mutate(order = factor(order, levels =   dat_plot_order$order)) %>% 
    ggplot(., aes(Time, order,fill= value1))+ 
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


#p_heat <- plot_grid(p_heat_touch, p_heat_Pin, p_heat_Har, p_heat_Expect_touch, p_heat_Expect_pain,nrow  = 1)
p_heat <- plot_grid(p_heat_touch, p_heat_Pin, p_heat_Har,nrow  = 1)

setwd("~cchen/Documents/neuroscience/Pn\ project/Figure/PDF/")
cairo_pdf("p_heat.pdf", width = 90/25.6, height = 90/25.6, family = "Arial")
p_heat
dev.off()
rm(p_heat_touch, p_heat_Pin, p_heat_Har, p_heat_Expect_touch, p_heat_Expect_pain)

## plot the traces
p_com_trace <- dat_trace_combine_re %>% 
  filter(Group != "Expect_pain") %>% 
  filter(Group !="Expect_touch") %>% 
  ddply(., .(Group, Time, ID), summarise, value = mean(value1)) %>% 
  ddply(., .(Group, Time), summarise,mean=mean(value),sd=sd(value),se=sd(value)/sqrt(length(value))) %>%  
  mutate(Group = factor(Group, levels = c("touch",  "Pin" , "Har"))) %>% 
  ggplot(., aes(Time, mean, color = Group))+
  geom_line()+
  geom_ribbon(aes(ymin=mean-se, ymax=mean+se, fill=Group), alpha=0.1, linetype=0)+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  scale_y_continuous(limits = c(-0.1, 0.2)) +
  theme(legend.title = element_blank())

setwd("~cchen/Documents/neuroscience/Pn\ project/Figure/PDF/")
cairo_pdf("p_com_trace.pdf", width = 110/25.6, height = 65/25.6, family = "Arial")
p_com_trace
dev.off()


## compare pin and expect pain-----
p_com_trace <- dat_trace_combine_re %>% 
  filter(Group == "Expect_touch" | Group =="touch" ) %>% 
  ddply(., .(Group, Time, ID), summarise, value = mean(value1)) %>% 
  ddply(., .(Group, Time), summarise,mean=mean(value),sd=sd(value),se=sd(value)/sqrt(length(value))) %>%  
  #mutate(Group = factor(Group, levels = c("Pin","Expect_pain"))) %>% 
  ggplot(., aes(Time, mean, color = Group))+
  geom_line()+
  geom_ribbon(aes(ymin=mean-se, ymax=mean+se, fill=Group), alpha=0.1, linetype=0)+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  scale_y_continuous(limits = c(-0.1, 0.1)) +
  theme(legend.title = element_blank())


## plot the trace to show stim response
c_extract_align <- function(file_trace){
  ID <- str_extract(file_trace, regex("m\\d+"))
  
  dat_trace <-   raveio::read_mat(file_trace)
  
  global_ID <- dat_trace[[1]] %>% 
    as_tibble() %>% 
    mutate_all(na_if, 0) %>% 
    drop_na() %>% 
    slice(all_of(sample(1:34, 10)))
  
  
  t_har <- stim_t_list_har[ID] %>% 
    unlist() %>% 
    unname()
  
  t_pin <- stim_t_list_pin[ID] %>% 
    unlist() %>% 
    unname()
  
  t_touch <- stim_t_list_touch[ID] %>% 
    unlist() %>% 
    unname()

  ## for har  

  har_cells <- pull(global_ID[,2])
  p_trace_har <- dat_trace[[2]][[2]] %>% 
    as_tibble() %>% 
    dplyr::select(all_of(har_cells)) %>% 
    mutate(across(.fns = scale)) %>% 
    add_column(Time = 1: nrow(.)) %>% 
    pivot_longer(-Time) %>% 
    ggplot(., aes(Time, value, color = name))+
    geom_line()+
    #scale_y_continuous(limits = plot_range)+
    facet_grid(rows = vars(name))+
    theme_void()+
    theme(legend.position = "none")+
    theme(strip.text.y = element_blank()) +
    geom_vline(xintercept = c(t_har))
  
  setwd("~cchen/Documents/neuroscience/Pn\ project/Figure/PDF/")
  cairo_pdf("p_trace_har.pdf", width = 90/25.6, height = 90/25.6, family = "Arial")
  p_trace_har
  dev.off()
  

  
  for (i in 1:length(trace_num)) {
    global_cell <- pull(global_ID[,trace_num[i]])
    
    dat_stim_trace <- dat_trace[[2]][[trace_num[i]]] %>% 
      as_tibble() %>% 
      dplyr::select(all_of(global_cell)) %>% 
      mutate(across(.fns = scale))
    
    stim_t <- t_stim[[i]]
    
    dat_trace_t <- vector(mode = "list", length = length(stim_t))
    for (j in seq_along(stim_t)){
      t1_p <- (stim_t[j]-20+1):(stim_t[j]+40)
      
      dat_stim1 <- dat_stim_trace %>% 
        slice(t1_p)
      
      #dat_stim1 <- aggregate(dat_stim1,list(rep(1:(nrow(dat_stim1)%/%n+1),each=n,len=nrow(dat_stim1))),mean)[-1]
      
      dat_trace_t[[j]] <- dat_stim1[1:20,] %>%  ## baseline as -2 to 0
        colMeans(., na.rm = T) %>% 
        sweep(dat_stim1, 2, ., FUN = "-") %>% 
        as_tibble() %>% 
        add_column(Time = stim_time)
      
    }
    dat_stim_average[[i]] <- dat_trace_t %>% 
      do.call(rbind,.) %>% 
      pivot_longer(-Time) %>% 
      ddply(., .(Time, name), summarise, value1 = mean(value)) %>% 
      add_column(ID = ID)
    
  }
  
  return(dat_stim_average)
}


