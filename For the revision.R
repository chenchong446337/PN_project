## for the revision, 05152023


## 1. correlation analysis of neuron activity and locomotion of animals------

cc_neuron_cor_fun <- function(ca_file){
  dat_trace1 <- raveio::read_mat(ca_file) [c(4, 8)]
  mouse_ID <- str_extract(ca_file, regex("m\\d+"))
  
  ## import tracking file
  
  tracking_file <- str_c("~cchen/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/behavior_tracking/", mouse_ID) %>% 
    list.files(., full.names = T)
  
  
  if (length(tracking_file) > 2) {
    tracking_file <- tracking_file[c(3, 7)]
  } else {
    tracking_file <- tracking_file
  }
  
  
  ## correlation analysis from D4 to D7
 
  day_r <- rep(0, 2)
  day_r_shuffle <- rep(0, 2)
  for (i in c(1, 2)){
    day_velocity <- read.csv(tracking_file[i], header = F, stringsAsFactors = F) %>% 
      slice(-(1:3)) %>% 
      select(1:3) %>% 
      mutate_all(as.numeric) %>% 
      mutate(vx = (V2 - lag(V2)), vy = (V3 - lag(V3)),
             velocity = sqrt(vx^2 + vy^2)) %>% 
      mutate(velocity_mean = zoo::rollmean(velocity, k = 5, na.pad = TRUE))
    
    day_ca <- dat_trace1[[i]] %>%
      as_tibble() %>% 
      apply(., 2, scale) %>% 
      rowMeans() 
    
    day_ca_shuffle <- sample(day_ca)
    
    length_com <- min(length(day_velocity$velocity_mean), length(day_ca))
    
    day_r[i]<- cor(day_velocity$velocity_mean[1:length_com], day_ca[1:length_com], use = "complete.obs")
    day_r_shuffle[i]<- cor(day_velocity$velocity_mean[1:length_com], day_ca_shuffle[1:length_com], use = "complete.obs")
    
  }
  
  dat_r <- c(day_r, day_r_shuffle)
  return(dat_r)
}


mouse_file <- as.list(list.files("~cchen/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/FD-processed", full.names = T))


dat_ca_loc_r <- mapply(cc_neuron_cor_fun, mouse_file, SIMPLIFY = F) 


p_ca_loc_r <- do.call(cbind, dat_ca_loc_r) %>% 
  as_tibble() %>% 
  mutate(Group = c("Pre", "Post", "Pre", "Post")) %>% 
  mutate(Type = c("Real", "Real", "Shuffle", "Shuffle")) %>% 
  mutate(Group = factor(Group, levels = c("Pre", "Post"))) %>% 
  mutate(Type = factor(Type, levels = c("Real", "Shuffle"))) %>%
  filter(Group != "Post") %>% 
  pivot_longer(-c(Group, Type)) %>% 
  ggplot(., aes( Type, value, fill= Type))+
  geom_boxplot( outlier.shape = NA)+
  geom_point(aes(colour = Group,shape = Type),position=position_jitterdodge(jitter.width = .2), size = 2)+
  scale_shape_manual(values=c(20, 18))+
  scale_color_manual(values=c("darkcyan", "indianred"))+
  labs(x="", y="Correlation coefficient (r)")+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  scale_y_continuous(limits = c(-0.2, 0.2), expand = c(0, 0))+
  theme(legend.position = "none")

t_ca_loc <- do.call(cbind, dat_ca_loc_r) %>% 
  as_tibble() %>% 
  mutate(Group = c("Pre", "Post", "Pre", "Post")) %>% 
  mutate(Type = c("Real", "Real", "Shuffle", "Shuffle")) %>% 
  mutate(Group = factor(Group, levels = c("Pre", "Post"))) %>% 
  mutate(Type = factor(Type, levels = c("Real", "Shuffle"))) %>% 
  pivot_longer(-c(Group, Type)) %>% 
  aov(value~Type,.)

summary(t_ca_loc)

setwd("~cchen/Documents/neuroscience/Pn\ project/Figure/PDF/")
cairo_pdf("p_ca_loc.pdf", width = 40/25.6, height = 60/25.6, family = "Arial")
p_ca_loc_r
dev.off()

p_ca_loc_r_com <- do.call(cbind, dat_ca_loc_r) %>% 
  as_tibble() %>% 
  mutate(Group = c("Pre", "Post", "Pre", "Post")) %>% 
  mutate(Type = c("Real", "Real", "Shuffle", "Shuffle")) %>% 
  pivot_longer(-c(Group, Type)) %>% 
  ggplot(., aes(Type, value, fill= Type))+
  geom_boxplot( outlier.shape = NA)+
  geom_jitter(position=position_jitterdodge(), shape=1)


## check the ca2+ activity during fast moving

cc_neuron_cor_fun1 <- function(ca_file){
  dat_trace1 <- raveio::read_mat(ca_file) [c(4, 8)]
  mouse_ID <- str_extract(ca_file, regex("m\\d+"))
  
  ## import tracking file
  
  tracking_file <- str_c("~cchen/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/behavior_tracking/", mouse_ID) %>% 
    list.files(., full.names = T)
  
  
  if (length(tracking_file) > 2) {
    tracking_file <- tracking_file[c(3, 7)]
  } else {
    tracking_file <- tracking_file
  }
  
  
  ## correlation analysis from D4 to D7
  mouse_ID_high <- c("m855", "m857")
  r_pix_mm <- 165 / 387*20/100 # to cm/s
  r_pix_mm1 <- 165 / 387*10/100
  
  
  mean_ca_peak <- rep(0, 20)
  veloc_peak <- rep(0, 20)
  for (i in c(1, 2)){
    
    r_convert <- ifelse(mouse_ID %in% mouse_ID_high, r_pix_mm1, r_pix_mm)
    day_velocity <- read.csv(tracking_file[i], header = F, stringsAsFactors = F) %>% 
      slice(-(1:3)) %>% 
      select(1:3) %>% 
      mutate_all(as.numeric) %>% 
      mutate(vx = (V2 - lag(V2)), vy = (V3 - lag(V3)),
             velocity = sqrt(vx^2 + vy^2)) %>% 
      mutate(velocity = velocity * r_convert) %>% 
      mutate(velocity_mean = zoo::rollmean(velocity, k = 5, na.pad = TRUE))
    
    day_ca <- dat_trace1[[i]] %>%
      as_tibble() %>% 
      apply(., 2, scale) %>% 
      rowMeans() 
    
    ## mouse ID have high recorded frequency
    
    
    if (mouse_ID %in% mouse_ID_high) {
      day_ca <- day_ca[seq(1, length(day_ca), by = 2)]
    } else {
      day_ca <- day_ca
    }
    
    ## find peaks of velocity trace
    peak_indices <- findpeaks(day_velocity$velocity_mean, nups = 2, minpeakdistance = 10,npeaks = 10)
    
    # veloc_peak[((i-1)*10 + 1): (i*10)] <- peak_indices[,1]
    # plot(day_velocity$velocity_mean, type = "l")
    # points(peak_indices[,2], peak_indices[,1], col = "red")
    
    ## calculate the activity of neurons during peak velocity
    mean_ca_peak1 <- rep(0, 10)
    mean_veloc_peak <- rep(0, 10)
    for (j in c(1:10)) {
      peak_time <- (peak_indices[,2][j] -9) : (peak_indices[,2][j] + 10)
      
      ## mean Ca2+ activity during this period
      
      mean_ca_peak1[j] <- mean(day_ca[peak_time], na.rm = T)
      mean_veloc_peak[j] <- mean(day_velocity$velocity_mean[peak_time], na.rm = T)
      
    }
    
    mean_ca_peak[((i-1)*10 + 1): (i*10)] <- mean_ca_peak1
    veloc_peak[((i-1)*10 + 1): (i*10)] <- mean_veloc_peak
    
  }
  
  dat_ca_veloc <- tibble(veloc = veloc_peak, Acti = mean_ca_peak, ID = mouse_ID)
  return(dat_ca_veloc)
}


mouse_file <- as.list(list.files("~cchen/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/FD-processed", full.names = T))


dat_ca_loc_corr <- mapply(cc_neuron_cor_fun1, mouse_file, SIMPLIFY = F) 


p_dat_ca_corr <- dat_ca_loc_corr %>% 
  do.call(rbind,.) %>% 
  as_tibble() %>% 
  mutate(Group = rep(rep(c("Pre", "Post"), each= 10), 6)) %>% 
  filter(Group != "Post") %>% 
  ggplot(., aes(veloc, Acti))+
  geom_point()+
  geom_smooth(method = "lm", se = T)+
  labs(x = "Velocity (mm/s)", y = "Ca2+ activity (z-score)") +
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))

dat_dat_ca_corr_t <- dat_ca_loc_corr %>% 
  do.call(rbind,.) %>% 
  as_tibble() %>% 
  mutate(Group = rep(rep(c("Pre", "Post"), each= 10), 6)) %>% 
  filter(Group =="Pre")

cor.test(dat_dat_ca_corr_t$veloc, dat_dat_ca_corr_t$Acti)

setwd("~cchen/Documents/neuroscience/Pn\ project/Figure/PDF/")
cairo_pdf("p_dat_ca_cor.pdf", width = 80/25.6, height = 60/25.6, family = "Arial")
p_dat_ca_corr
dev.off()

## Calculate the mean velocity during peak Ca activity


cc_neuron_cor_fun2 <- function(ca_file){
  dat_trace1 <- raveio::read_mat(ca_file) [c(4, 8)]
  mouse_ID <- str_extract(ca_file, regex("m\\d+"))
  
  ## import tracking file
  
  tracking_file <- str_c("~cchen/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/behavior_tracking/", mouse_ID) %>% 
    list.files(., full.names = T)
  
  
  if (length(tracking_file) > 2) {
    tracking_file <- tracking_file[c(3, 7)]
  } else {
    tracking_file <- tracking_file
  }
  
  
  ## correlation analysis from D4 to D7
  
  mean_ca_peak <- rep(0, 20)
  veloc_peak <- rep(0, 20)
  for (i in c(1, 2)){
    day_velocity <- read.csv(tracking_file[i], header = F, stringsAsFactors = F) %>% 
      slice(-(1:3)) %>% 
      select(1:3) %>% 
      mutate_all(as.numeric) %>% 
      mutate(vx = (V2 - lag(V2)), vy = (V3 - lag(V3)),
             velocity = sqrt(vx^2 + vy^2)) %>% 
      mutate(velocity_mean = zoo::rollmean(velocity, k = 5, na.pad = TRUE))
    
    day_ca <- dat_trace1[[i]] %>% 
      rowMeans() %>% 
      zoo::rollmean(., k = 5, na.pad = TRUE)
    
    ## mouse ID have high recorded frequency
    mouse_ID_high <- c("m855", "m857")
    
    if (mouse_ID %in% mouse_ID_high) {
      day_ca <- day_ca[seq(1, length(day_ca), by = 2)]
    } else {
      day_ca <- day_ca
    }
    
    
    ## find peaks of velocity trace
    peak_indices <- findpeaks(day_ca, nups = 2, minpeakdistance = 100,npeaks = 10)
    
    #plot(day_ca, type = "l")
    # points(peak_indices[,2], peak_indices[,1], col = "red")
    
    ## calculate the activity of neurons during peak velocity
    mean_ca_peak1 <- rep(0, 10)
    mean_veloc_peak <- rep(0, 10)
    for (j in c(1:10)) {
      peak_time <- (peak_indices[,2][j] -9) : (peak_indices[,2][j] + 10)
      peak_time <- peak_time[peak_time>0]
      
      ## mean Ca2+ activity during this period
      
      mean_ca_peak1[j] <- mean(day_ca[peak_time], na.rm = T)
      mean_veloc_peak[j] <- mean(day_velocity$velocity_mean[peak_time], na.rm = T)
      
      
    }
    
    mean_ca_peak[((i-1)*10 + 1): (i*10)] <- mean_ca_peak1
    veloc_peak[((i-1)*10 + 1): (i*10)] <- mean_veloc_peak
  }
  
  dat_ca_veloc <- tibble(veloc = veloc_peak, Acti = mean_ca_peak, ID = mouse_ID)
  return(dat_ca_veloc)
}


mouse_file <- as.list(list.files("~cchen/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/FD-processed", full.names = T))


dat_ca_loc_corr <- mapply(cc_neuron_cor_fun2, mouse_file, SIMPLIFY = F) 


p_dat_ca_corr1 <- dat_ca_loc_corr %>% 
  do.call(rbind,.) %>% 
  as_tibble() %>% 
  mutate(Group = rep(rep(c("Pre", "Post"), each= 10), 6)) %>% 
  filter(Group != "Post") %>% 
  filter(veloc < 90) %>% 
  ggplot(., aes(veloc, Acti))+
  geom_point()+
  geom_smooth(method = "lm", se = T)+
  labs(x = "Velocity (mm/s)", y = "Ca2+ activity (z-score)") +
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))

cor.test(p_dat_ca_corr1$veloc, p_dat_ca_corr1$Acti)

setwd("~cchen/Documents/neuroscience/Pn\ project/Figure/PDF/")
cairo_pdf("p_dat_ca_cor1.pdf", width = 80/25.6, height = 60/25.6, family = "Arial")
p_dat_ca_corr1
dev.off()


## 2. check Ca2+ activity of rACC neurons during pain behavior-----

dat_anti_pain <- read.xlsx("~cchen/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/matlab_analysis/behavior_tracking/pn_anti_pain_behavior.xlsx") 


c_miniscope_pain <- function(ca_file) {
  
  dat_trace1 <- raveio::read_mat(ca_file) [[8]]
  mouse_ID <- str_extract(ca_file, regex("m\\d+"))
  
  mouse_ID_high <- c("m855", "m857")
  
  t_stim_licking <- dat_anti_pain %>% 
    filter(ID == mouse_ID) %>% 
    select(First_licking) %>% 
    unlist() %>% 
    unname()
  
  t_stim_rearing <- dat_anti_pain %>% 
    filter(ID == mouse_ID) %>% 
    select(First_rearing) %>% 
    unlist() %>% 
    unname()
  
  t_stim_licking <- ifelse(mouse_ID %in% mouse_ID_high, t_stim_licking*2, t_stim_licking)
  
  t_stim_rearing <- ifelse(mouse_ID %in% mouse_ID_high, t_stim_rearing*2, t_stim_rearing)
  
  dat_trace <- dat_trace1 %>%
    as_tibble() %>% 
    apply(., 2, scale)
  
  ## import and format the data

  
  

  stim_time<- seq(-2, 1.9, by=0.1)
  ## number of rows to be binned
  n <- 2 # 0.05*10=0.5
  
  ## extract cell activity when they cross the border

  dat_stim_licking1 <- dat_trace[(t_stim_licking-40):(t_stim_licking+40-1),] 
  dat_stim_licking1 <- aggregate(dat_stim_licking1,list(rep(1:(nrow(dat_stim_licking1)%/%n+1),each=n,len=nrow(dat_stim_licking1))),mean)[-1]
    
  dat_stim_licking <- dat_stim_licking1[1:5,] %>%  ## baseline as -2 to 0
    colMeans(., na.rm = T) %>% 
    sweep(dat_stim_licking1, 2, ., FUN = "-")
  
  
  dat_stim_rearing1 <- dat_trace[(t_stim_rearing-40):(t_stim_rearing+40-1),] 
  dat_stim_rearing1 <- aggregate(dat_stim_rearing1,list(rep(1:(nrow(dat_stim_rearing1)%/%n+1),each=n,len=nrow(dat_stim_rearing1))),mean)[-1]
  
  dat_stim_rearing <- dat_stim_rearing1[1:10,] %>%  ## baseline as -2 to 0
    colMeans(., na.rm = T) %>% 
    sweep(dat_stim_rearing1, 2, ., FUN = "-")
  
  dat_cell_trace_average <- rbind(dat_stim_licking, dat_stim_rearing) %>% 
    as_tibble() 
  
  colnames(dat_cell_trace_average) <- str_c(mouse_ID,"Cell", 1: ncol(dat_cell_trace_average))
  
  
  return(dat_cell_trace_average)
  
}


mouse_file <- as.list(list.files("~cchen/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/FD-processed", full.names = T))

stim_time <- seq(-2, 1.9, by=0.1)
dat_ca_licking <- mapply(c_miniscope_pain, mouse_file, SIMPLIFY = F) 

dat_pain<- dat_ca_licking %>% 
  do.call(cbind,.) %>% 
  as_tibble() %>% 
  add_column(Time = rep(stim_time, 2)) %>% 
  add_column(Type = rep(c("Licking", "Rearing"), each = 40))
  



score_range <- dat_pain %>% 
  select(-c(Time, Type)) %>% 
  range()

cell_order_licking <-  dat_pain %>% 
  pivot_longer(-c(Time, Type)) %>% 
  ddply(.,.(Type, name), summarise, mean = mean(value)) %>% 
  filter(Type =="Licking") %>% 
  arrange(mean)

cell_order_rearing <-  dat_pain %>% 
  pivot_longer(-c(Time, Type)) %>% 
  ddply(.,.(Type, name), summarise, mean = mean(value)) %>% 
  filter(Type =="Rearing") %>% 
  arrange(mean)
  
p_heat_licking <- dat_pain %>% 
  pivot_longer(-c(Time, Type)) %>%
  filter(Type == "Licking") %>% 
  mutate(name = factor(name, levels = cell_order_licking$name)) %>% 
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


p_heat_rearing <- dat_pain %>% 
  pivot_longer(-c(Time, Type)) %>%
  filter(Type == "Rearing") %>% 
  mutate(name = factor(name, levels = cell_order_rearing$name)) %>% 
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

p_heat_com <- plot_grid(p_heat_licking, p_heat_rearing, nrow = 1)

setwd("~cchen/Documents/neuroscience/Pn\ project/Figure/PDF/")
cairo_pdf("p_heat_com.pdf", width = 100/25.6, height = 60/25.6, family = "Arial")
p_heat_com
dev.off()

p_trace_pain <- dat_pain %>% 
  pivot_longer(-c(Time, Type)) %>% 
  ddply(.,.(Time, Type), summarise,n=length(value),mean=mean(value),sd=sd(value),se=sd(value)/sqrt(length(value))) %>% 
  ggplot(., aes(Time, mean, color = Type))+
  geom_line()+
  geom_ribbon(aes(ymin=mean-se, ymax=mean+se), alpha=0.1, linetype=0)+
  labs(x="Time (s)", y="Fluorescence (s.d.)")+
  geom_vline(xintercept = 0, col= 'red', linetype=2)+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  scale_y_continuous(limits = c(-0.1, 0.2))+
  theme(legend.title = element_blank(),legend.position = c(0.2, 0.8))
  
setwd("~cchen/Documents/neuroscience/Pn\ project/Figure/PDF/")
cairo_pdf("p_trace_pain.pdf", width = 90/25.6, height = 60/25.6, family = "Arial")
p_trace_pain
dev.off()


## for behavioral bouts------
file_path <- list.files("~cchen/Documents/neuroscience/Pn\ project/Data_analysis/For\ revision/mouse\ behavior\ bouts", full.names = T)

dat_bouts <- vector(mode = "list", length = length(file_path))
for (i in c(1:length(file_path))) {
  path_file <- file_path[[i]]
  
  Group1 <- read.xlsx(path_file) %>% 
    as_tibble() %>% 
    slice(34) %>% 
    pull(2)
  
    
  df <- read.xlsx(path_file) %>% 
    as_tibble() %>% 
    drop_na() %>%
    { colnames(.) <- .[1, ]; . } %>%
    slice(-1, -2)
  
  if (any("Jumping" %in% df$Behavior)) {
    trimmed_df <- df %>%
      slice(1:(which(Behavior == "Jumping")[2]))
  } else {
    # In this example, we're keeping the entire dataset if "Jumping" doesn't exist
    trimmed_df <- df
  }  
  
  dat_bouts[[i]] <- df %>% 
    mutate(Group = Group1) %>% 
    mutate(ID = str_c("m", i))
  
}

df <- dat_bouts %>% 
  do.call(rbind,.) %>% 
  select(-c("Recording time","Subject" )) %>% 
  rename(Time = "Trial time") %>% 
  mutate(Time = as.numeric(Time)) %>% 
  mutate(Time = Time *3) %>% 
  group_by(ID) %>%
  mutate(adjust_time = first(Time[Behavior == "Crossing" & Event == "state start"])) %>%
  ungroup() %>%
  mutate(Time = Time - adjust_time) %>%
  select(-adjust_time) %>% 
  filter(Time >= 0)




# Extract 'start' and 'stop' events
df_start <- df[df$Event == "state start",]
df_stop <- df[df$Event == "state stop",]

df_duration <- data.frame(Time.x = df_start$Time,
                          Time.y = df_stop$Time,
                          Behavior = df_start$Behavior,
                          Group = df_start$Group,
                          ID = df_start$ID) %>% 
  as_tibble() %>% 
  mutate(Duration = Time.y - Time.x)




# Create the raster plot using the earlier provided function
plot_behavior <- function(group_name) {
  data_group <- df_duration[df_duration$Group == group_name, ]
  ggplot(data_group, aes(xmin = Time.x, xmax = Time.y, ymin = as.numeric(as.factor(ID)), ymax = as.numeric(as.factor(ID)) + 1)) +
    geom_rect(aes(fill = Behavior)) +
    scale_fill_manual(values = c("Crossing" = "steelblue", "Rearing" = "darkolivegreen", "Licking" = "coral", "Jumping"= "slateblue")) +
    labs(x = "Time", y = "ID") + 
    theme(axis.line.x = element_line(),
          axis.line.y = element_line(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank(),
          axis.title = element_text(family = "Arial", size = 12, face = "plain"),
          axis.text.y = element_blank())+
    scale_x_continuous(limits = c(0, 180))
}

# Plot for 'Cond' group
p_con <- plot_behavior("Cond")

p_control <- plot_behavior("Ctrl")


p_bouts <- plot_grid(p_control, p_con, ncol = 1)

## calculate the frequecy of each nocifensive behavior during early phase
df_licking <- df %>% 
  filter(Time < 45) %>% 
  filter(Behavior == "Licking") %>% 
  filter(Event =="state start") %>% 
  ddply(., .(Group, ID), summarise, n = length(Event))

df_licking_time <- df %>% 
  filter(Event =="state stop") %>% 
  group_by(Group, ID) %>%
  summarise(Duration = max(Time)) %>% 
  ungroup() %>% 
  mutate(Duration = 60)

p_licking_freq <- right_join(df_licking, df_licking_time) %>% 
  as_tibble() %>% 
  replace(is.na(.), 0) %>% 
  mutate(Group = factor(Group, levels= c("Ctrl", "Cond"))) %>% 
  mutate(n = n*60/45) %>% 
  ggplot(., aes(Group, n, group = Group))+
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(aes(colour= Group,shape = Group),width = 0.2, size = 2)+
  scale_color_manual(values=c("darkcyan", "indianred"))+
  scale_shape_manual(values=c(20, 18))+
  labs(x="", y="Freq. of licking (#/min)")+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  scale_y_continuous(limits = c(-1, 10))+
  theme(legend.position = "none")

t_licking_freq <- right_join(df_licking, df_licking_time) %>% 
  as_tibble() %>% 
  replace(is.na(.), 0) %>% 
  mutate(Group = factor(Group, levels= c("Ctrl", "Cond"))) %>% 
  wilcox.test(n~Group,.)


## for rearing
df_rearing <- df %>% 
  filter(Time < 45) %>% 
  filter(Behavior == "Rearing") %>% 
  filter(Event =="state start") %>% 
  ddply(., .(Group, ID), summarise, n = length(Event))

df_rearing_time <- df %>% 
  filter(Event =="state stop") %>% 
  group_by(Group, ID) %>%
  summarise(Duration = max(Time)) %>% 
  ungroup() %>% 
  mutate(Duration = 60)

p_rearing_freq <- right_join(df_rearing, df_rearing_time) %>% 
  as_tibble() %>% 
  replace(is.na(.), 0) %>% 
  mutate(Group = factor(Group, levels= c("Ctrl", "Cond"))) %>% 
  ggplot(., aes(Group, n*60/45, group = Group))+
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(aes(colour= Group,shape = Group),width = 0.2, size = 2)+
  scale_color_manual(values=c("darkcyan", "indianred"))+
  scale_shape_manual(values=c(20, 18))+
  labs(x="", y="Freq. of rearing (#/min)")+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  scale_y_continuous(limits = c(-1, 10))+
  theme(legend.position = "none")

t_licking_freq <- right_join(df_rearing, df_rearing_time) %>% 
  as_tibble() %>% 
  replace(is.na(.), 0) %>% 
  mutate(Group = factor(Group, levels= c("Ctrl", "Cond"))) %>% 
  wilcox.test(n~Group,.)


p_PAC <- plot_grid(p_licking_freq, p_rearing_freq, nrow = 1)

setwd("~cchen/Documents/neuroscience/Pn\ project/Figure/PDF/")
cairo_pdf("p_PAC.pdf", width = 70/25.6, height = 60/25.6, family = "Arial")
p_PAC
dev.off()


## calculate the duration of each nocifensive behavior during early phase
df_licking_duration <- df_duration %>% 
  filter(Time.y <= 45) %>% 
  filter(Behavior == "Licking") %>% 
  ddply(., .(Group, ID), summarise, t_duration = sum(Duration))

df_licking_time <- df %>% 
  filter(Event =="state stop") %>% 
  group_by(Group, ID) %>%
  summarise(Duration1 = max(Time)) %>% 
  ungroup() 

p_licking_duration <- right_join(df_licking_duration, df_licking_time) %>% 
  as_tibble() %>% 
  replace(is.na(.), 0) %>% 
  mutate(Group = factor(Group, levels= c("Ctrl", "Cond"))) %>% 
  ggplot(., aes(Group, t_duration, group = Group))+
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(aes(colour= Group,shape = Group),width = 0.2, size = 2)+
  scale_color_manual(values=c("darkcyan", "indianred"))+
  scale_shape_manual(values=c(20, 18))+
  labs(x="", y="Duration of licking (s)")+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  scale_y_continuous(limits = c(-1, 8))+
  theme(legend.position = "none")

t_licking_freq <- right_join(df_licking_duration, df_licking_time) %>% 
  as_tibble() %>% 
  replace(is.na(.), 0) %>% 
  mutate(Group = factor(Group, levels= c("Ctrl", "Cond"))) %>% 
  wilcox.test(t_duration~Group,.)


## for rearing
df_rearing_duration <- df_duration %>% 
  filter(Time.y <= 45) %>% 
  filter(Behavior == "Rearing") %>% 
  ddply(., .(Group, ID), summarise, t_duration = sum(Duration))

df_licking_time <- df %>% 
  filter(Event =="state stop") %>% 
  group_by(Group, ID) %>%
  summarise(Duration1 = max(Time)) %>% 
  ungroup() 

p_rearing_duration <- right_join(df_rearing_duration, df_licking_time) %>% 
  as_tibble() %>% 
  replace(is.na(.), 0) %>% 
  mutate(Group = factor(Group, levels= c("Ctrl", "Cond"))) %>% 
  ggplot(., aes(Group, t_duration, group = Group))+
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(aes(colour= Group,shape = Group),width = 0.2, size = 2)+
  scale_color_manual(values=c("darkcyan", "indianred"))+
  scale_shape_manual(values=c(20, 18))+
  labs(x="", y="Duration of rearing (s)")+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  scale_y_continuous(limits = c(-1, 17))+
  theme(legend.position = "none")

t_duration_freq <- right_join(df_rearing_duration, df_licking_time) %>% 
  as_tibble() %>% 
  replace(is.na(.), 0) %>% 
  mutate(Group = factor(Group, levels= c("Ctrl", "Cond"))) %>% 
  wilcox.test(t_duration~Group,.)


p_PAC <- plot_grid(p_licking_duration, p_rearing_duration, nrow = 1)
setwd("~cchen/Documents/neuroscience/Pn\ project/Figure/PDF/")
cairo_pdf("p_PAC.pdf", width = 60/25.6, height = 60/25.6, family = "Arial")
p_PAC
dev.off()
## for total duration
df_licking_duration <- df_duration %>% 
  filter(Behavior == "Licking") %>% 
  ddply(., .(Group, ID), summarise, t_duration = sum(Duration))


p_licking_duration <- df_licking_duration %>% 
  mutate(Group = factor(Group, levels= c("Ctrl", "Cond"))) %>% 
  ggplot(., aes(Group, t_duration, group = Group))+
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(aes(colour= Group,shape = Group),width = 0.2, size = 2)+
  scale_color_manual(values=c("darkcyan", "indianred"))+
  scale_shape_manual(values=c(20, 18))+
  labs(x="", y="Duration of licking (s)")+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  scale_y_continuous(limits = c(-1, 15))+
  theme(legend.position = "none")

t_licking_freq <- right_join(df_licking_duration, df_licking_time) %>% 
  as_tibble() %>% 
  replace(is.na(.), 0) %>% 
  mutate(Group = factor(Group, levels= c("Ctrl", "Cond"))) %>% 
  wilcox.test(t_duration~Group,.)


## for rearing
df_rearing_duration <- df_duration %>% 
  filter(Time.y <= 45) %>% 
  filter(Behavior == "Rearing") %>% 
  ddply(., .(Group, ID), summarise, t_duration = sum(Duration))

df_licking_time <- df %>% 
  filter(Event =="state stop") %>% 
  group_by(Group, ID) %>%
  summarise(Duration1 = max(Time)) %>% 
  ungroup() 

p_rearing_duration <- right_join(df_rearing_duration, df_licking_time) %>% 
  as_tibble() %>% 
  replace(is.na(.), 0) %>% 
  mutate(Group = factor(Group, levels= c("Ctrl", "Cond"))) %>% 
  ggplot(., aes(Group, t_duration, group = Group))+
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(aes(colour= Group,shape = Group),width = 0.2, size = 2)+
  scale_color_manual(values=c("darkcyan", "indianred"))+
  scale_shape_manual(values=c(20, 18))+
  labs(x="", y="Duration of rearing (s)")+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  scale_y_continuous(limits = c(-1, 17))+
  theme(legend.position = "none")

t_duration_freq <- right_join(df_rearing_duration, df_licking_time) %>% 
  as_tibble() %>% 
  replace(is.na(.), 0) %>% 
  mutate(Group = factor(Group, levels= c("Ctrl", "Cond"))) %>% 
  wilcox.test(t_duration~Group,.)


p_PAC <- plot_grid(p_licking_duration, p_rearing_duration, nrow = 1)
setwd("~cchen/Documents/neuroscience/Pn\ project/Figure/PDF/")
cairo_pdf("p_PAC.pdf", width = 80/25.6, height = 60/25.6, family = "Arial")
p_PAC
dev.off()

## calculate the frequecy of each nocifensive behavior during whole behavior
df_licking <- df %>% 
  filter(Behavior == "Licking") %>% 
  filter(Event =="state start") %>% 
  ddply(., .(Group, ID), summarise, n = length(Event))

df_licking_time <- df %>% 
  filter(Event =="state stop") %>% 
  group_by(Group, ID) %>%
  summarise(Duration = max(Time)) %>% 
  ungroup() 

p_licking_freq <- right_join(df_licking, df_licking_time) %>% 
  as_tibble() %>% 
  replace(is.na(.), 0) %>% 
  mutate(Freq = n/Duration) %>% 
  mutate(Freq = Freq*60) %>% 
  mutate(Group = factor(Group, levels= c("Ctrl", "Cond"))) %>% 
  ggplot(., aes(Group, Freq, group = Group))+
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(aes(colour= Group,shape = Group),width = 0.2, size = 2)+
  scale_color_manual(values=c("darkcyan", "indianred"))+
  scale_shape_manual(values=c(20, 18))+
  labs(x="", y="Freq. of licking (#/min)")+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  scale_y_continuous(limits = c(0, 6))+
  theme(legend.position = "none")

t_licking_freq <- right_join(df_licking, df_licking_time) %>% 
  as_tibble() %>% 
  replace(is.na(.), 0) %>% 
  mutate(Freq = n/Duration) %>% 
  mutate(Freq = Freq*60) %>% 
  mutate(Group = factor(Group, levels= c("Ctrl", "Cond"))) %>% 
  wilcox.test(Freq~Group,.)

## for rearing
df_rearing <- df %>% 
  filter(Behavior == "Rearing") %>% 
  filter(Event =="state start") %>% 
  ddply(., .(Group, ID), summarise, n = length(Event))

df_rearing_time <- df %>% 
  filter(Event =="state stop") %>% 
  group_by(Group, ID) %>%
  summarise(Duration = max(Time)) %>% 
  ungroup() 

p_rearing_freq <- right_join(df_rearing, df_rearing_time) %>% 
  as_tibble() %>% 
  replace(is.na(.), 0) %>% 
  mutate(Freq = n/Duration) %>% 
  mutate(Freq = Freq*60) %>% 
  mutate(Group = factor(Group, levels= c("Ctrl", "Cond"))) %>% 
  ggplot(., aes(Group, Freq, group = Group))+
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(aes(colour= Group,shape = Group),width = 0.2, size = 2)+
  scale_color_manual(values=c("darkcyan", "indianred"))+
  scale_shape_manual(values=c(20, 18))+
  labs(x="", y="Freq. of rearing (#/min)")+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  scale_y_continuous(limits = c(0, 10))+
  theme(legend.position = "none")

t_licking_freq <- right_join(df_rearing, df_rearing_time) %>% 
  as_tibble() %>% 
  replace(is.na(.), 0) %>% 
  mutate(Freq = n/Duration) %>% 
  mutate(Freq = Freq*60) %>% 
  mutate(Group = factor(Group, levels= c("Ctrl", "Cond"))) %>% 
  wilcox.test(Freq~Group,.)


p_PAC <- plot_grid(p_licking_freq, p_rearing_freq, nrow = 1)

setwd("~cchen/Documents/neuroscience/Pn\ project/Figure/PDF/")
cairo_pdf("p_PAC.pdf", width = 80/25.6, height = 60/25.6, family = "Arial")
p_PAC
dev.off()
## for the cummulative plot
p_cum_licking <- df %>% 
  filter(Behavior == "Licking") %>% 
  filter(Event =="state start") %>% 
  mutate(Group = factor(Group, levels= c("Ctrl", "Cond"))) %>% 
  ggplot(., aes(Time, color = Group)) + stat_ecdf(geom = "step")+
  scale_color_manual(values=c("darkcyan", "indianred"))+
  labs(x="Time of licking (s)", y="Cumulative probability ")+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  theme(legend.position = "none")

t_cum_test <-  df %>% 
  filter(Behavior == "Licking") %>% 
  filter(Event =="state start") %>% 
  mutate(Group = factor(Group, levels= c("Ctrl", "Cond"))) %>% 
  ks.test(Time~Group,.)


p_cum_rearing <- df %>% 
  filter(Behavior == "Rearing") %>% 
  filter(Event =="state start") %>% 
  mutate(Group = factor(Group, levels= c("Ctrl", "Cond"))) %>% 
  ggplot(., aes(Time, color = Group)) + stat_ecdf(geom = "step")+
  scale_color_manual(values=c("darkcyan", "indianred"))+
  labs(x="Time of rearing (s)", y="Cumulative probability ")+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  theme(legend.position = c(0.2, 0.8), legend.title = element_blank())

t_cum_rearing <- df %>% 
  filter(Behavior == "Rearing") %>% 
  filter(Event =="state start") %>% 
  mutate(Group = factor(Group, levels= c("Ctrl", "Cond"))) %>% 
  ks.test(Time~Group,.)

p_cum_pac <- plot_grid(p_cum_licking, p_cum_rearing, nrow = 1)
setwd("~cchen/Documents/neuroscience/Pn\ project/Figure/PDF/")
cairo_pdf("p_cum_pac.pdf", width = 120/25.6, height = 60/25.6, family = "Arial")
p_cum_pac
dev.off()
## For mechanical and thermal pain after conditioning-----
p_mecha_con <- read.xlsx("~cchen/Documents/neuroscience/Pn\ project/Data_analysis/For\ revision/Pn_revision.xlsx", sheet = 'VonFrey') %>% 
  mutate(Group = factor(Group, levels = c('Ctrl', 'Cond'))) %>% 
  ggplot(., aes(Group,Threshold, fill=Group))+
  geom_boxplot(outlier.shape = NA)+
  geom_point(aes(colour = Group,shape = Group),position=position_jitterdodge(jitter.width = .5), size = 2)+
  scale_fill_manual(values = c("white", "white" ))+
  scale_shape_manual(values=c(20, 18))+
  scale_color_manual(values=c("darkcyan", "indianred"))+
  labs(x="", y="Pain threshold (g)")+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  theme(legend.title = element_blank(), legend.position = "none")

p_thermal_con <- read.xlsx("~cchen/Documents/neuroscience/Pn\ project/Data_analysis/For\ revision/Pn_revision.xlsx", sheet = 'Har') %>% 
  mutate(Group = factor(Group, levels = c('Ctrl', 'Cond'))) %>% 
  ddply(.,.(Group, ID), summarise, latency1 = mean(Latency)) %>% 
  ggplot(., aes(Group,latency1, fill=Group))+
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(aes(colour = Group,shape = Group),position=position_jitterdodge(jitter.width = 0.5), size = 2)+
  scale_fill_manual(values = c("white", "white" ))+
  scale_shape_manual(values=c(20, 18))+
  scale_color_manual(values=c("darkcyan", "indianred"))+
  labs(x="", y="Paw withdrawal latency (s)")+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  theme(legend.title = element_blank(), legend.position = "none")

p_pain_cond <- plot_grid(p_mecha_con, p_thermal_con, nrow = 1)

setwd("~cchen/Documents/neuroscience/Pn\ project/Figure/PDF/")
cairo_pdf("p_pain_cond.pdf", width = 80/25.6, height = 60/25.6, family = "Arial")
p_pain_cond
dev.off()


## for female mice during PAC------
p_ratio <- read.xlsx("~cchen/Documents/neuroscience/Pn\ project/Data_analysis/For\ revision/Pn_revision.xlsx", sheet = 'Female') %>% 
  select(Group, ID, Time_hot_d3, Time_hot_d7) %>% 
  pivot_longer(-c(Group, ID)) %>% 
  mutate(name = factor(name, levels=c("Time_hot_d3", "Time_hot_d7"))) %>% 
  mutate(Group = factor(Group, levels = c("Ctrl", "Cond"))) %>% 
  ggplot(., aes(interaction(name, Group), value, fill= name))+
  geom_boxplot( outlier.shape = NA)+
  geom_line(aes(group=interaction(ID, Group)), colour="gray90")+
  geom_point(aes(colour = Group,shape = Group),position=position_jitterdodge(jitter.width = .2))+
  scale_fill_manual(values = c("white", "white"))+
  scale_shape_manual(values=c(20, 18))+
  scale_color_manual(values=c("darkcyan", "indianred"))+
  labs(x="", y="Time spend in 2nd plate (%)")+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  scale_y_continuous(limits = c(0, 100), expand = c(0, 0))+
  theme(legend.position = 'none', legend.title = element_blank())

p_return_latency <- read.xlsx("~cchen/Documents/neuroscience/Pn\ project/Data_analysis/For\ revision/Pn_revision.xlsx", sheet = 'Female') %>% 
  mutate(Latency = (Crossing_back - Latency_1st_crossing)/10) %>% 
  mutate(Latency1 = Crossing_back_d3/10) %>% 
  select(Group, ID, Latency, Latency1) %>% 
  pivot_longer(-c(Group, ID)) %>% 
  mutate(name = factor(name, levels=c("Latency1", "Latency"))) %>% 
  mutate(Group = factor(Group, levels = c("Ctrl", "Cond"))) %>% 
  ggplot(., aes(interaction(name, Group), value, fill= name))+
  geom_boxplot( outlier.shape = NA)+
  geom_line(aes(group=interaction(ID, Group)), colour="gray90")+
  geom_point(aes(colour = Group,shape = Group),position=position_jitterdodge(jitter.width = .2))+
  scale_fill_manual(values = c("white", "white"))+
  scale_shape_manual(values=c(20, 18))+
  scale_color_manual(values=c("darkcyan", "indianred"))+
  labs(x="", y="Latency of crossing back (s)")+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  theme(legend.position = 'none', legend.title = element_blank())

p_female <- plot_grid(p_return_latency, p_ratio, nrow = 1)
  
setwd("~cchen/Documents/neuroscience/Pn\ project/Figure/PDF/")
cairo_pdf("p_female.pdf", width = 110/25.6, height = 60/25.6, family = "Arial")
p_female
dev.off()

p_female_licking <- read.xlsx("~cchen/Documents/neuroscience/Pn\ project/Data_analysis/For\ revision/Pn_revision.xlsx", sheet = 'Female') %>% 
  select(Group, Latency_1st_crossing, First_licking) %>% 
  mutate(Latency = (First_licking - Latency_1st_crossing)/10) %>% 
  mutate(Group = factor(Group, levels = c("Ctrl", "Cond"))) %>% 
  ggplot(., aes(Group,Latency, fill=Group))+
  geom_boxplot(outlier.shape = NA)+
  geom_point(aes(colour = Group,shape = Group),position=position_jitterdodge(jitter.width = .5), size = 2)+
  scale_fill_manual(values = c("white", "white" ))+
  scale_shape_manual(values=c(20, 18))+
  scale_color_manual(values=c("darkcyan", "indianred"))+
  labs(x="", y="Latency of 1st licking (s)")+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  theme(legend.title = element_blank(), legend.position = "none")

t_female_licking <- read.xlsx("~cchen/Documents/neuroscience/Pn\ project/Data_analysis/For\ revision/Pn_revision.xlsx", sheet = 'Female') %>% 
  select(Group, Latency_1st_crossing, First_licking) %>% 
  mutate(Latency = (First_licking - Latency_1st_crossing)/10) %>% 
  mutate(Group = factor(Group, levels = c("Ctrl", "Cond"))) %>%
  wilcox.test(Latency~Group, .)

p_female_rearing <- read.xlsx("~cchen/Documents/neuroscience/Pn\ project/Data_analysis/For\ revision/Pn_revision.xlsx", sheet = 'Female') %>% 
  select(Group, Latency_1st_crossing, First_rearing) %>% 
  mutate(Latency = (First_rearing - Latency_1st_crossing)/10) %>% 
  mutate(Group = factor(Group, levels = c("Ctrl", "Cond"))) %>% 
  ggplot(., aes(Group,Latency, fill=Group))+
  geom_boxplot(outlier.shape = NA)+
  geom_point(aes(colour = Group,shape = Group),position=position_jitterdodge(jitter.width = .5), size = 2)+
  scale_fill_manual(values = c("white", "white" ))+
  scale_shape_manual(values=c(20, 18))+
  scale_color_manual(values=c("darkcyan", "indianred"))+
  labs(x="", y="Latency of 1st rearing (s)")+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  theme(legend.title = element_blank(), legend.position = "none")

t_female_rearing <- read.xlsx("~cchen/Documents/neuroscience/Pn\ project/Data_analysis/For\ revision/Pn_revision.xlsx", sheet = 'Female') %>% 
  select(Group, Latency_1st_crossing, First_rearing) %>% 
  mutate(Latency = (First_rearing - Latency_1st_crossing)/10) %>% 
  mutate(Group = factor(Group, levels = c("Ctrl", "Cond"))) %>% 
  wilcox.test(Latency~Group, .)

  
p_female_jumping <- read.xlsx("~cchen/Documents/neuroscience/Pn\ project/Data_analysis/For\ revision/Pn_revision.xlsx", sheet = 'Female') %>% 
  select(Group, Latency_1st_crossing, Jumping) %>% 
  mutate(Latency = (Jumping - Latency_1st_crossing)/10) %>% 
  mutate(Group = factor(Group, levels = c("Ctrl", "Cond"))) %>% 
  ggplot(., aes(Group,Latency, fill=Group))+
  geom_boxplot(outlier.shape = NA)+
  geom_point(aes(colour = Group,shape = Group),position=position_jitterdodge(jitter.width = .5), size = 2)+
  scale_fill_manual(values = c("white", "white" ))+
  scale_shape_manual(values=c(20, 18))+
  scale_color_manual(values=c("darkcyan", "indianred"))+
  labs(x="", y="Latency of 1st jumping (s)")+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  theme(legend.title = element_blank(), legend.position = "none")
  

t_female_jumping <- read.xlsx("~cchen/Documents/neuroscience/Pn\ project/Data_analysis/For\ revision/Pn_revision.xlsx", sheet = 'Female') %>% 
  select(Group, Latency_1st_crossing, Jumping) %>% 
  mutate(Latency = (Jumping - Latency_1st_crossing)/10) %>% 
  mutate(Group = factor(Group, levels = c("Ctrl", "Cond"))) %>% 
  wilcox.test(Latency~Group, .)


p_female_pain <- plot_grid(p_female_licking, p_female_rearing, p_female_jumping, nrow = 1)
  
setwd("~cchen/Documents/neuroscience/Pn\ project/Figure/PDF/")
cairo_pdf("p_female_pain.pdf", width = 100/25.6, height = 60/25.6, family = "Arial")
p_female_pain
dev.off()  
  
## Mearing the duration of PAC-induced pain relief-----

## compare the latency of crossing back
p_duration_crossing_back_d3 <- read.xlsx("~cchen/Documents/neuroscience/Pn\ project/Data_analysis/For\ revision/Pn_revision.xlsx", sheet = 'PAC_anti_duration') %>% 
  select(Group, Duration, Latency_1st_crossing, Crossing_back) %>%   
  mutate(Latency = (Crossing_back - Latency_1st_crossing)/10) %>% 
  filter(Duration == 3) %>% 
  mutate(Group = factor(Group, levels = c("Ctrl", "Cond"))) %>% 
  ggplot(., aes(Group, Latency, fill= Group))+
  geom_boxplot( outlier.shape = NA)+
  geom_point(aes(colour = Group,shape = Group),position=position_jitterdodge(jitter.width = .2))+
  scale_fill_manual(values = c("white", "white"))+
  scale_shape_manual(values=c(20, 18))+
  scale_color_manual(values=c("darkcyan", "indianred"))+
  labs(x="", y="Latency of crossing back (s)")+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  theme(legend.position = 'none', legend.title = element_blank())

p_duration_crossing_back_d7 <- read.xlsx("~cchen/Documents/neuroscience/Pn\ project/Data_analysis/For\ revision/Pn_revision.xlsx", sheet = 'PAC_anti_duration') %>% 
  select(Group, Duration, Latency_1st_crossing, Crossing_back) %>%   
  mutate(Latency = (Crossing_back - Latency_1st_crossing)/10) %>% 
  filter(Duration == 7) %>% 
  mutate(Group = factor(Group, levels = c("Ctrl", "Cond"))) %>% 
  ggplot(., aes(Group, Latency, fill= Group))+
  geom_boxplot( outlier.shape = NA)+
  geom_point(aes(colour = Group,shape = Group),position=position_jitterdodge(jitter.width = .2))+
  scale_fill_manual(values = c("white", "white"))+
  scale_shape_manual(values=c(20, 18))+
  scale_color_manual(values=c("darkcyan", "indianred"))+
  labs(x="", y="Latency of crossing back (s)")+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  theme(legend.position = 'none', legend.title = element_blank())
  

p_duration_plate_d3 <- read.xlsx("~cchen/Documents/neuroscience/Pn\ project/Data_analysis/For\ revision/Pn_revision.xlsx", sheet = 'PAC_anti_duration') %>% 
  select(Group, Duration, Time_hot) %>%   
  filter(Duration == 3) %>% 
  mutate(Group = factor(Group, levels = c("Ctrl", "Cond"))) %>% 
  ggplot(., aes(Group, Time_hot, fill= Group))+
  geom_boxplot( outlier.shape = NA)+
  geom_point(aes(colour = Group,shape = Group),position=position_jitterdodge(jitter.width = .2))+
  scale_fill_manual(values = c("white", "white"))+
  scale_shape_manual(values=c(20, 18))+
  scale_color_manual(values=c("darkcyan", "indianred"))+
  labs(x="", y="Time spend in 2nd plate (%)")+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  theme(legend.position = 'none', legend.title = element_blank())

p_duration_plate_d7 <- read.xlsx("~cchen/Documents/neuroscience/Pn\ project/Data_analysis/For\ revision/Pn_revision.xlsx", sheet = 'PAC_anti_duration') %>% 
  select(Group, Duration, Time_hot) %>%   
  filter(Duration == 7) %>% 
  mutate(Group = factor(Group, levels = c("Ctrl", "Cond"))) %>% 
  ggplot(., aes(Group, Time_hot, fill= Group))+
  geom_boxplot( outlier.shape = NA)+
  geom_point(aes(colour = Group,shape = Group),position=position_jitterdodge(jitter.width = .2))+
  scale_fill_manual(values = c("white", "white"))+
  scale_shape_manual(values=c(20, 18))+
  scale_color_manual(values=c("darkcyan", "indianred"))+
  labs(x="", y="Time spend in 2nd plate (%)")+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  theme(legend.position = 'none', legend.title = element_blank())

p_duration_licking_d3 <- read.xlsx("~cchen/Documents/neuroscience/Pn\ project/Data_analysis/For\ revision/Pn_revision.xlsx", sheet = 'PAC_anti_duration') %>% 
  select(Group, Duration, Latency_1st_crossing, First_licking) %>%   
  mutate(Latency = (First_licking - Latency_1st_crossing)/10) %>% 
  filter(Duration == 3) %>% 
  mutate(Group = factor(Group, levels = c("Ctrl", "Cond"))) %>% 
  ggplot(., aes(Group, Latency, fill= Group))+
  geom_boxplot( outlier.shape = NA)+
  geom_point(aes(colour = Group,shape = Group),position=position_jitterdodge(jitter.width = .2))+
  scale_fill_manual(values = c("white", "white"))+
  scale_shape_manual(values=c(20, 18))+
  scale_color_manual(values=c("darkcyan", "indianred"))+
  labs(x="", y="Latency of 1st licking (s)")+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  scale_y_continuous(limits = c(5, 90))+
  theme(legend.position = 'none', legend.title = element_blank())

p_duration_licking_d7 <- read.xlsx("~cchen/Documents/neuroscience/Pn\ project/Data_analysis/For\ revision/Pn_revision.xlsx", sheet = 'PAC_anti_duration') %>% 
  select(Group, Duration, Latency_1st_crossing, First_licking) %>%   
  mutate(Latency = (First_licking - Latency_1st_crossing)/10) %>% 
  filter(Duration == 7) %>% 
  mutate(Group = factor(Group, levels = c("Ctrl", "Cond"))) %>% 
  ggplot(., aes(Group, Latency, fill= Group))+
  geom_boxplot( outlier.shape = NA)+
  geom_point(aes(colour = Group,shape = Group),position=position_jitterdodge(jitter.width = .2))+
  scale_fill_manual(values = c("white", "white"))+
  scale_shape_manual(values=c(20, 18))+
  scale_color_manual(values=c("darkcyan", "indianred"))+
  labs(x="", y="Latency of 1st licking (s)")+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  scale_y_continuous(limits = c(5, 90))+
  theme(legend.position = 'none', legend.title = element_blank())


p_duration_rearing_d3 <- read.xlsx("~cchen/Documents/neuroscience/Pn\ project/Data_analysis/For\ revision/Pn_revision.xlsx", sheet = 'PAC_anti_duration') %>% 
  select(Group, Duration, Latency_1st_crossing, First_rearing) %>%   
  mutate(Latency = (First_rearing - Latency_1st_crossing)/10) %>% 
  filter(Duration == 3) %>% 
  mutate(Group = factor(Group, levels = c("Ctrl", "Cond"))) %>% 
  ggplot(., aes(Group, Latency, fill= Group))+
  geom_boxplot( outlier.shape = NA)+
  geom_point(aes(colour = Group,shape = Group),position=position_jitterdodge(jitter.width = .2))+
  scale_fill_manual(values = c("white", "white"))+
  scale_shape_manual(values=c(20, 18))+
  scale_color_manual(values=c("darkcyan", "indianred"))+
  labs(x="", y="Latency of 1st rearing (s)")+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  scale_y_continuous(limits = c(0, 180))+
  theme(legend.position = 'none', legend.title = element_blank())


p_duration_rearing_d7 <- read.xlsx("~cchen/Documents/neuroscience/Pn\ project/Data_analysis/For\ revision/Pn_revision.xlsx", sheet = 'PAC_anti_duration') %>% 
  select(Group, Duration, Latency_1st_crossing, First_rearing) %>%   
  mutate(Latency = (First_rearing - Latency_1st_crossing)/10) %>% 
  filter(Duration == 7) %>% 
  mutate(Group = factor(Group, levels = c("Ctrl", "Cond"))) %>% 
  ggplot(., aes(Group, Latency, fill= Group))+
  geom_boxplot( outlier.shape = NA)+
  geom_point(aes(colour = Group,shape = Group),position=position_jitterdodge(jitter.width = .2))+
  scale_fill_manual(values = c("white", "white"))+
  scale_shape_manual(values=c(20, 18))+
  scale_color_manual(values=c("darkcyan", "indianred"))+
  labs(x="", y="Latency of 1st rearing (s)")+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  scale_y_continuous(limits = c(0, 180))+
  theme(legend.position = 'none', legend.title = element_blank())

p_duration_jumping_d3 <- read.xlsx("~cchen/Documents/neuroscience/Pn\ project/Data_analysis/For\ revision/Pn_revision.xlsx", sheet = 'PAC_anti_duration') %>% 
  select(Group, Duration, Latency_1st_crossing, Jumping) %>%   
  mutate(Latency = (Jumping - Latency_1st_crossing)/10) %>% 
  filter(Duration == 3) %>% 
  mutate(Group = factor(Group, levels = c("Ctrl", "Cond"))) %>% 
  ggplot(., aes(Group, Latency, fill= Group))+
  geom_boxplot( outlier.shape = NA)+
  geom_point(aes(colour = Group,shape = Group),position=position_jitterdodge(jitter.width = .2))+
  scale_fill_manual(values = c("white", "white"))+
  scale_shape_manual(values=c(20, 18))+
  scale_color_manual(values=c("darkcyan", "indianred"))+
  labs(x="", y="Latency of 1st jumping (s)")+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  scale_y_continuous(limits = c(50, 200))+
  theme(legend.position = 'none', legend.title = element_blank())

p_duration_jumping_d7 <- read.xlsx("~cchen/Documents/neuroscience/Pn\ project/Data_analysis/For\ revision/Pn_revision.xlsx", sheet = 'PAC_anti_duration') %>% 
  select(Group, Duration, Latency_1st_crossing, Jumping) %>%   
  mutate(Latency = (Jumping - Latency_1st_crossing)/10) %>% 
  filter(Duration == 7) %>% 
  mutate(Group = factor(Group, levels = c("Ctrl", "Cond"))) %>% 
  ggplot(., aes(Group, Latency, fill= Group))+
  geom_boxplot( outlier.shape = NA)+
  geom_point(aes(colour = Group,shape = Group),position=position_jitterdodge(jitter.width = .2))+
  scale_fill_manual(values = c("white", "white"))+
  scale_shape_manual(values=c(20, 18))+
  scale_color_manual(values=c("darkcyan", "indianred"))+
  labs(x="", y="Latency of 1st jumping (s)")+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  scale_y_continuous(limits = c(50, 200))+
  theme(legend.position = 'none', legend.title = element_blank())

p_duration_pain_d3 <- plot_grid(p_duration_crossing_back_d3, p_duration_plate_d3,p_duration_licking_d3, p_duration_rearing_d3, p_duration_jumping_d3, nrow = 1)

setwd("~cchen/Documents/neuroscience/Pn\ project/Figure/PDF/")
cairo_pdf("p_duration_pain_d3.pdf", width = 150/25.6, height = 60/25.6, family = "Arial")
p_duration_pain_d3
dev.off()  

## mean + error bar to save space
dat_cross1 <- dat_anti_wt %>% 
  filter(variable=="latency2") %>% 
  filter(Day == "Test") %>% 
  add_column(Duration = 1) %>% 
  select(Group, Duration, value) %>% 
  rename(Latency = "value") %>% 
  mutate(Group = ifelse(Group=="Cond.", "Cond", "Ctrl"))


p_duration_crossing_back <- read.xlsx("~cchen/Documents/neuroscience/Pn\ project/Data_analysis/For\ revision/Pn_revision.xlsx", sheet = 'PAC_anti_duration') %>% 
  select(Group, Duration, Latency_1st_crossing, Crossing_back) %>%   
  mutate(Latency = (Crossing_back - Latency_1st_crossing)/10) %>% 
  select(Group, Duration, Latency) %>% 
  bind_rows(., dat_cross1) %>% 
  drop_na() %>% 
  mutate(Duration = factor(Duration, levels = c(1, 3, 7))) %>% 
  mutate(Group = factor(Group, levels = c("Ctrl", "Cond"))) %>% 
  ddply(., .(Group, Duration), summarise, mean = mean(Latency), se=sd(Latency)/sqrt(length(Latency))) %>% 
  ggplot(., aes(Duration, mean, group = Group, color =Group))+
  geom_line(color = "gray90")+
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se, color = Group), width=.3) +
  geom_point(aes(colour =Group,shape =Group), size = 3)+
  scale_shape_manual(values=c(20, 18))+
  scale_color_manual(values=c("darkcyan", "indianred"))+
  labs(x="", y="Latency of crossing back (s)")+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  scale_y_continuous(limits = c(10, 100))+
  theme(legend.position = 'none', legend.title = element_blank())


dat_anti_plate1 <- dat_anti_wt %>% 
  filter(variable=="ratio" ) %>% 
  mutate(value=value*100) %>% 
  filter(Day == "Test") %>% 
  add_column(Duration = 1) %>% 
  select(Group, Duration, value) %>% 
  rename(Time_hot = "value") %>% 
  mutate(Group = ifelse(Group=="Cond.", "Cond", "Ctrl"))
  
p_duration_plate <- read.xlsx("~cchen/Documents/neuroscience/Pn\ project/Data_analysis/For\ revision/Pn_revision.xlsx", sheet = 'PAC_anti_duration') %>% 
  select(Group, Duration, Time_hot) %>%   
  bind_rows(., dat_anti_plate1) %>% 
  drop_na() %>% 
  mutate(Duration = factor(Duration, levels = c(1, 3, 7))) %>% 
  mutate(Group = factor(Group, levels = c("Ctrl", "Cond"))) %>% 
  ddply(., .(Group, Duration), summarise, mean = mean(Time_hot), se=sd(Time_hot)/sqrt(length(Time_hot))) %>% 
  ggplot(., aes(Duration, mean, group = Group, color =Group))+
  geom_line(color = "gray90")+
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se, color = Group), width=.3) +
  geom_point(aes(colour =Group,shape =Group), size = 3)+
  scale_shape_manual(values=c(20, 18))+
  scale_color_manual(values=c("darkcyan", "indianred"))+
  labs(x="", y="Time spend in 2nd plate (%)")+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  scale_y_continuous(limits = c(20, 100))+
  theme(legend.position = 'none', legend.title = element_blank())

p_duration_time <- plot_grid(p_duration_crossing_back, p_duration_plate, nrow = 1)

##data from previous experiment
dat_wt_manual <- read.xlsx("~cchen/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/wt_ctrl_07042020/wt_ctrl_07042020_manualy.xlsx") %>% 
  mutate('Rearing'= Total_rearing/Total_frame*600, 'First_rearing'=First_rearing*0.1, 
         'Guarding'= (Guarding - Frame_1st_crossing)*0.1, 'Acceleration'= (Acceleration -Frame_1st_crossing )*0.1) %>% 
  pivot_longer(-c(Group, ID), names_to = "variable", values_to = 'value') %>% 
  mutate(Group=factor(Group, levels=c("Ctrl", "Cond.")))

dat_licking1 <- dat_wt_manual %>% 
  filter(variable=="First_licking") %>% 
  add_column(Duration = 1) %>% 
  select(Group, Duration, value) %>% 
  rename(Latency = "value") %>% 
  mutate(Group = ifelse(Group=="Cond.", "Cond", "Ctrl"))
  

p_duration_licking <- read.xlsx("~cchen/Documents/neuroscience/Pn\ project/Data_analysis/For\ revision/Pn_revision.xlsx", sheet = 'PAC_anti_duration') %>% 
  select(Group, Duration, Latency_1st_crossing, First_licking) %>%   
  mutate(Latency = (First_licking - Latency_1st_crossing)/10) %>% 
  select(Group, Duration, Latency) %>% 
  bind_rows(., dat_licking1) %>% 
  drop_na() %>% 
  mutate(Duration = factor(Duration, levels = c(1, 3, 7))) %>% 
  mutate(Group = factor(Group, levels = c("Ctrl", "Cond"))) %>% 
  ddply(., .(Group, Duration), summarise, mean = mean(Latency), se=sd(Latency)/sqrt(length(Latency))) %>% 
  ggplot(., aes(Duration, mean, group = Group, color =Group))+
  geom_line(color = "gray90")+
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se, color = Group), width=.3) +
  geom_point(aes(colour =Group,shape =Group), size = 3)+
  scale_shape_manual(values=c(20, 18))+
  scale_color_manual(values=c("darkcyan", "indianred"))+
  labs(x="", y="Latency of 1st licking (s)")+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  scale_y_continuous(limits = c(5, 90))+
  theme(legend.position = 'none', legend.title = element_blank())

dat_rearing1 <- dat_wt_manual %>% 
  filter(variable=="First_rearing") %>% 
  add_column(Duration = 1) %>% 
  select(Group, Duration, value) %>% 
  rename(Latency = "value") %>% 
  mutate(Group = ifelse(Group=="Cond.", "Cond", "Ctrl"))

p_duration_rearing <- read.xlsx("~cchen/Documents/neuroscience/Pn\ project/Data_analysis/For\ revision/Pn_revision.xlsx", sheet = 'PAC_anti_duration') %>% 
  select(Group, Duration, Latency_1st_crossing, First_rearing) %>%   
  mutate(Latency = (First_rearing - Latency_1st_crossing)/10) %>% 
  select(Group, Duration, Latency) %>% 
  bind_rows(., dat_licking1) %>% 
  drop_na() %>% 
  mutate(Duration = factor(Duration, levels = c(1, 3, 7))) %>% 
  mutate(Group = factor(Group, levels = c("Ctrl", "Cond"))) %>% 
  ddply(., .(Group, Duration), summarise, mean = mean(Latency), se=sd(Latency)/sqrt(length(Latency))) %>% 
  ggplot(., aes(Duration, mean, group = Group, color =Group))+
  geom_line(color = "gray90")+
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se, color = Group), width=.3) +
  geom_point(aes(colour =Group,shape =Group), size = 3)+
  scale_shape_manual(values=c(20, 18))+
  scale_color_manual(values=c("darkcyan", "indianred"))+
  labs(x="", y="Latency of 1st rearing (s)")+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  scale_y_continuous(limits = c(0, 100))+
  theme(legend.position = 'none', legend.title = element_blank())

dat_jump1 <- dat_wt_manual %>% 
  filter(variable=="Jump") %>% 
  add_column(Duration = 1) %>% 
  select(Group, Duration, value) %>% 
  mutate(value = value/10) %>% 
  rename(Latency = "value") %>% 
  mutate(Group = ifelse(Group=="Cond.", "Cond", "Ctrl"))

p_duration_jumping <- read.xlsx("~cchen/Documents/neuroscience/Pn\ project/Data_analysis/For\ revision/Pn_revision.xlsx", sheet = 'PAC_anti_duration') %>% 
  select(Group, Duration, Latency_1st_crossing, Jumping) %>%   
  mutate(Latency = (Jumping - Latency_1st_crossing)/10) %>% 
  select(Group, Duration, Latency) %>% 
  bind_rows(., dat_jump1) %>% 
  drop_na() %>% 
  mutate(Duration = factor(Duration, levels = c(1, 3, 7))) %>% 
  mutate(Group = factor(Group, levels = c("Ctrl", "Cond"))) %>% 
  ddply(., .(Group, Duration), summarise, mean = mean(Latency), se=sd(Latency)/sqrt(length(Latency))) %>% 
  ggplot(., aes(Duration, mean, group = Group, color =Group))+
  geom_line(color = "gray90")+
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se, color = Group), width=.3) +
  geom_point(aes(colour =Group,shape =Group), size = 3)+
  scale_shape_manual(values=c(20, 18))+
  scale_color_manual(values=c("darkcyan", "indianred"))+
  labs(x="", y="Latency of 1st jumping (s)")+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  scale_y_continuous(limits = c(50, 200))+
  theme(legend.position = 'none', legend.title = element_blank())

p_duration_pain <- plot_grid(p_duration_crossing_back, p_duration_plate, p_duration_licking, p_duration_rearing, p_duration_jumping,nrow = 1)

setwd("~cchen/Documents/neuroscience/Pn\ project/Figure/PDF/")
cairo_pdf("p_duration_pain.pdf", width = 230/25.6, height = 60/25.6, family = "Arial")
p_duration_pain
dev.off()  

## Examine the role of OPRD1+ neurons during PAC-conditioning-----
p_oprd_opto_crossing <- read.xlsx("~cchen/Documents/neuroscience/Pn\ project/Data_analysis/For\ revision/Pn_revision.xlsx", sheet = 'PAC_oprd_cond_anti') %>% 
  select(Group, Latency_1st_crossing) %>%   
  mutate(Latency = (Latency_1st_crossing)/10) %>% 
  mutate(Group = factor(Group, levels = c("eYFP", "eNpHR"))) %>% 
  ggplot(., aes(Group, Latency, fill= Group))+
  geom_boxplot( outlier.shape = NA)+
  geom_point(aes(colour = Group,shape = Group),position=position_jitterdodge(jitter.width = .2), size = 2)+
  scale_fill_manual(values = c("white", "white"))+
  scale_shape_manual(values=c(20, 18))+
  scale_color_manual(values=c("darkolivegreen3", "dodgerblue4"))+
  labs(x="", y="Latency of crossing (s)")+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  theme(legend.position = 'none', legend.title = element_blank())

t_oprd_opto_crossing <- read.xlsx("~cchen/Documents/neuroscience/Pn\ project/Data_analysis/For\ revision/Pn_revision.xlsx", sheet = 'PAC_oprd_cond_anti') %>% 
  select(Group, Latency_1st_crossing) %>%   
  mutate(Latency = (Latency_1st_crossing)/10) %>% 
  mutate(Group = factor(Group, levels = c("eYFP", "eNpHR"))) %>% 
  wilcox.test(Latency~Group,.)

p_oprd_opto_crossing_back <- read.xlsx("~cchen/Documents/neuroscience/Pn\ project/Data_analysis/For\ revision/Pn_revision.xlsx", sheet = 'PAC_oprd_cond_anti') %>% 
  select(Group, Latency_1st_crossing, Crossing_back) %>%   
  mutate(Latency = (Crossing_back - Latency_1st_crossing)/10) %>% 
  mutate(Group = factor(Group, levels = c("eYFP", "eNpHR"))) %>% 
  ggplot(., aes(Group, Latency, fill= Group))+
  geom_boxplot( outlier.shape = NA)+
  geom_point(aes(colour = Group,shape = Group),position=position_jitterdodge(jitter.width = .2), size = 2)+
  scale_fill_manual(values = c("white", "white"))+
  scale_shape_manual(values=c(20, 18))+
  scale_color_manual(values=c("darkolivegreen3", "dodgerblue4"))+
  labs(x="", y="Latency of crossing back (s)")+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  theme(legend.position = 'none', legend.title = element_blank())

t_oprd_opto_crossing_back <- read.xlsx("~cchen/Documents/neuroscience/Pn\ project/Data_analysis/For\ revision/Pn_revision.xlsx", sheet = 'PAC_oprd_cond_anti') %>% 
  select(Group, Latency_1st_crossing, Crossing_back) %>%   
  mutate(Latency = (Crossing_back - Latency_1st_crossing)/10) %>% 
  mutate(Group = factor(Group, levels = c("eYFP", "eNpHR"))) %>% 
  wilcox.test(Latency~Group,.)

p_oprd_opto_plate <- read.xlsx("~cchen/Documents/neuroscience/Pn\ project/Data_analysis/For\ revision/Pn_revision.xlsx", sheet = 'PAC_oprd_cond_anti') %>% 
  select(Group, Time_hot) %>%   
  mutate(Group = factor(Group, levels = c("eYFP", "eNpHR"))) %>% 
  ggplot(., aes(Group, Time_hot, fill= Group))+
  geom_boxplot( outlier.shape = NA)+
  geom_point(aes(colour = Group,shape = Group),position=position_jitterdodge(jitter.width = .2), size = 2)+
  scale_fill_manual(values = c("white", "white"))+
  scale_shape_manual(values=c(20, 18))+
  scale_color_manual(values=c("darkolivegreen3", "dodgerblue4"))+
  labs(x="", y="Time spend in 2nd plate (%)")+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  theme(legend.position = 'none', legend.title = element_blank())

p_oprd_opto_licking <- read.xlsx("~cchen/Documents/neuroscience/Pn\ project/Data_analysis/For\ revision/Pn_revision.xlsx", sheet = 'PAC_oprd_cond_anti') %>% 
  select(Group, Latency_1st_crossing, First_licking) %>%   
  mutate(Latency = (First_licking - Latency_1st_crossing)/10) %>% 
  mutate(Group = factor(Group, levels = c("eYFP", "eNpHR"))) %>% 
  ggplot(., aes(Group, Latency, fill= Group))+
  geom_boxplot( outlier.shape = NA)+
  geom_point(aes(colour = Group,shape = Group),position=position_jitterdodge(jitter.width = .2), size = 2)+
  scale_fill_manual(values = c("white", "white"))+
  scale_shape_manual(values=c(20, 18))+
  scale_color_manual(values=c("darkolivegreen3", "dodgerblue4"))+
  labs(x="", y="Latency of 1st licking (s)")+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  theme(legend.position = 'none', legend.title = element_blank())

t_oprd_opto_licking <- read.xlsx("~cchen/Documents/neuroscience/Pn\ project/Data_analysis/For\ revision/Pn_revision.xlsx", sheet = 'PAC_oprd_cond_anti') %>% 
  select(Group, Latency_1st_crossing, First_licking) %>%   
  mutate(Latency = (First_licking - Latency_1st_crossing)/10) %>% 
  mutate(Group = factor(Group, levels = c("eYFP", "eNpHR"))) %>% 
  wilcox.test(Latency~ Group,.)
  

p_oprd_opto_rearing <- read.xlsx("~cchen/Documents/neuroscience/Pn\ project/Data_analysis/For\ revision/Pn_revision.xlsx", sheet = 'PAC_oprd_cond_anti') %>% 
  select(Group, Latency_1st_crossing, First_rearing) %>%   
  mutate(Latency = (First_rearing - Latency_1st_crossing)/10) %>% 
  mutate(Group = factor(Group, levels = c("eYFP", "eNpHR"))) %>% 
  ggplot(., aes(Group, Latency, fill= Group))+
  geom_boxplot( outlier.shape = NA)+
  geom_point(aes(colour = Group,shape = Group),position=position_jitterdodge(jitter.width = .2), size = 2)+
  scale_fill_manual(values = c("white", "white"))+
  scale_shape_manual(values=c(20, 18))+
  scale_color_manual(values=c("darkolivegreen3", "dodgerblue4"))+
  labs(x="", y="Latency of 1st rearing (s)")+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  theme(legend.position = 'none', legend.title = element_blank())

p_oprd_opto_jumping <- read.xlsx("~cchen/Documents/neuroscience/Pn\ project/Data_analysis/For\ revision/Pn_revision.xlsx", sheet = 'PAC_oprd_cond_anti') %>% 
  select(Group, Latency_1st_crossing, Jumping) %>%   
  mutate(Latency = (Jumping - Latency_1st_crossing)/10) %>% 
  mutate(Group = factor(Group, levels = c("eYFP", "eNpHR"))) %>% 
  ggplot(., aes(Group, Latency, fill= Group))+
  geom_boxplot( outlier.shape = NA)+
  geom_point(aes(colour = Group,shape = Group),position=position_jitterdodge(jitter.width = .2), size = 2)+
  scale_fill_manual(values = c("white", "white"))+
  scale_shape_manual(values=c(20, 18))+
  scale_color_manual(values=c("darkolivegreen3", "dodgerblue4"))+
  labs(x="", y="Latency of 1st jumping (s)")+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  theme(legend.position = 'none', legend.title = element_blank())


t_oprd_opto_jumping <- read.xlsx("~cchen/Documents/neuroscience/Pn\ project/Data_analysis/For\ revision/Pn_revision.xlsx", sheet = 'PAC_oprd_cond_anti') %>% 
  select(Group, Latency_1st_crossing, Jumping) %>%   
  mutate(Latency = (Jumping - Latency_1st_crossing)/10) %>% 
  mutate(Group = factor(Group, levels = c("eYFP", "eNpHR"))) %>% 
  wilcox.test(Latency~Group,.)

p_opto_oprd_conditioning <- plot_grid(p_oprd_opto_crossing, p_oprd_opto_crossing_back, p_oprd_opto_plate, nrow = 1)

p_opto_oprd_conditioning1 <- plot_grid(p_oprd_opto_licking, p_oprd_opto_rearing, p_oprd_opto_jumping, nrow = 1)


setwd("~cchen/Documents/neuroscience/Pn\ project/Figure/PDF/")
cairo_pdf("p_opto_oprd_conditniong.pdf", width = 100/25.6, height = 60/25.6, family = "Arial")
p_opto_oprd_conditioning
dev.off()


setwd("~cchen/Documents/neuroscience/Pn\ project/Figure/PDF/")
cairo_pdf("p_opto_oprd_conditniong1.pdf", width = 100/25.6, height = 60/25.6, family = "Arial")
p_opto_oprd_conditioning1
dev.off()

## DOR antagonist during PAC post test-----

p_DOR_ratio <- read.xlsx("~cchen/Documents/neuroscience/Pn\ project/Data_analysis/For\ revision/Pn_revision.xlsx", sheet = 'PAC_DOR_anti') %>% 
  select(Group, ID, Time_hot_d3, Time_hot_d7) %>% 
  pivot_longer(-c(Group, ID)) %>% 
  mutate(name = factor(name, levels=c("Time_hot_d3", "Time_hot_d7"))) %>% 
  mutate(Group = factor(Group, levels = c("Saline", "Dor", "Dor1", "Mor", "Mor1"))) %>% 
  ggplot(., aes(interaction(name, Group), value, fill= name))+
  geom_boxplot( outlier.shape = NA)+
  geom_line(aes(group=interaction(ID, Group)), colour="gray90")+
  geom_point(aes(colour = Group,shape = Group),position=position_jitterdodge(jitter.width = .2))+
  scale_fill_manual(values = c("white", "white"))+
  scale_shape_manual(values=c(20, 18))+
  scale_color_manual(values=c("darkcyan", "indianred"))+
  labs(x="", y="Time spend in 2nd plate (%)")+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  scale_y_continuous(limits = c(0, 100), expand = c(0, 0))+
  theme(legend.position = 'none', legend.title = element_blank())

p_DOR_return_latency <- read.xlsx("~cchen/Documents/neuroscience/Pn\ project/Data_analysis/For\ revision/Pn_revision.xlsx", sheet = 'PAC_DOR_anti') %>% 
  mutate(Latency = (Crossing_back - Latency_1st_crossing)/10) %>% 
  mutate(Latency1 = Crossing_back_d3/10) %>% 
  select(Group, ID, Latency, Latency1) %>% 
  pivot_longer(-c(Group, ID)) %>% 
  mutate(name = factor(name, levels=c("Latency1", "Latency"))) %>% 
  mutate(Group = factor(Group, levels = c("Saline", "Dor", "Dor1", "Mor", "Mor1"))) %>% 
  ggplot(., aes(interaction(name, Group), value, fill= name))+
  geom_boxplot( outlier.shape = NA)+
  geom_line(aes(group=interaction(ID, Group)), colour="gray90")+
  geom_point(aes(colour = Group,shape = Group),position=position_jitterdodge(jitter.width = .2))+
  scale_fill_manual(values = c("white", "white"))+
  scale_shape_manual(values=c(20, 18))+
  scale_color_manual(values=c("darkcyan", "indianred"))+
  labs(x="", y="Latency of crossing back (s)")+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  theme(legend.position = 'none', legend.title = element_blank())

p_female <- plot_grid(p_return_latency, p_ratio, nrow = 1)

setwd("~cchen/Documents/neuroscience/Pn\ project/Figure/PDF/")
cairo_pdf("p_female.pdf", width = 110/25.6, height = 60/25.6, family = "Arial")
p_female
dev.off()

p_DOR_licking <- read.xlsx("~cchen/Documents/neuroscience/Pn\ project/Data_analysis/For\ revision/Pn_revision.xlsx", sheet = 'PAC_DOR_anti') %>% 
  select(Group, Latency_1st_crossing, First_licking) %>% 
  mutate(Latency = (First_licking - Latency_1st_crossing)/10) %>% 
  mutate(Group = factor(Group, levels = c("Saline", "Dor", "Dor1", "Mor", "Mor1"))) %>% 
  ggplot(., aes(Group,Latency, fill=Group))+
  geom_boxplot(outlier.shape = NA)+
  geom_point(aes(colour = Group,shape = Group),position=position_jitterdodge(jitter.width = .5), size = 2)+
  scale_fill_manual(values = c("white", "white", "white", "white", "white"))+
  labs(x="", y="Latency of 1st licking (s)")+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  scale_y_continuous(limits = c(5, 80))+
  theme(legend.title = element_blank(), legend.position = "none")

t_DOR_licking <- read.xlsx("~cchen/Documents/neuroscience/Pn\ project/Data_analysis/For\ revision/Pn_revision.xlsx", sheet = 'PAC_DOR_anti') %>% 
  select(Group, Latency_1st_crossing, First_licking) %>% 
  mutate(Latency = (First_licking - Latency_1st_crossing)/10) %>% 
  mutate(Group = factor(Group, levels = c("Saline", "Dor", "Dor1", "Mor", "Mor1"))) %>% 
  ddply(., .(Group), summarise, mean= mean(Latency),se=sd(Latency)/sqrt(length(Latency)))
  aov(Latency~Group, .)

summary(t_DOR_licking)
TukeyHSD(t_DOR_licking)

p_DOR_rearing <- read.xlsx("~cchen/Documents/neuroscience/Pn\ project/Data_analysis/For\ revision/Pn_revision.xlsx", sheet = 'PAC_DOR_anti') %>% 
  select(Group, Latency_1st_crossing, First_rearing) %>% 
  mutate(Latency = (First_rearing - Latency_1st_crossing)/10) %>% 
  mutate(Group = factor(Group, levels = c("Saline", "Dor", "Dor1", "Mor", "Mor1"))) %>% 
  ggplot(., aes(Group,Latency, fill=Group))+
  geom_boxplot(outlier.shape = NA)+
  geom_point(aes(colour = Group,shape = Group),position=position_jitterdodge(jitter.width = .5), size = 2)+
  scale_fill_manual(values = c("white", "white", "white", "white", "white"))+
  labs(x="", y="Latency of 1st rearing (s)")+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  scale_y_continuous(limits = c(0, 130))+
  theme(legend.title = element_blank(), legend.position = "none")

t_DOR_rearing <- read.xlsx("~cchen/Documents/neuroscience/Pn\ project/Data_analysis/For\ revision/Pn_revision.xlsx", sheet = 'PAC_DOR_anti') %>% 
  select(Group, Latency_1st_crossing, First_rearing) %>% 
  mutate(Latency = (First_rearing - Latency_1st_crossing)/10) %>% 
  mutate(Group = factor(Group, levels = c("Saline", "Dor", "Dor1", "Mor", "Mor1"))) %>% 
  ddply(., .(Group), summarise, mean= mean(Latency),se=sd(Latency)/sqrt(length(Latency)))

  aov(Latency~Group, .)


summary(t_DOR_rearing)
TukeyHSD(t_DOR_rearing)

p_DOR_jumping <- read.xlsx("~cchen/Documents/neuroscience/Pn\ project/Data_analysis/For\ revision/Pn_revision.xlsx", sheet = 'PAC_DOR_anti') %>% 
  select(Group, Latency_1st_crossing, Jumping) %>% 
  mutate(Latency = (Jumping - Latency_1st_crossing)/10) %>% 
  mutate(Group = factor(Group, levels = c("Saline", "Dor", "Dor1", "Mor", "Mor1"))) %>% 
  ggplot(., aes(Group,Latency, fill=Group))+
  geom_boxplot(outlier.shape = NA)+
  geom_point(aes(colour = Group,shape = Group),position=position_jitterdodge(jitter.width = .5), size = 2)+
  scale_fill_manual(values = c("white", "white", "white", "white", "white"))+
  labs(x="", y="Latency of 1st jumping (s)")+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  scale_y_continuous(limits = c(20, 220))+
  theme(legend.title = element_blank(), legend.position = "none")


t_DOR_jumping <- read.xlsx("~cchen/Documents/neuroscience/Pn\ project/Data_analysis/For\ revision/Pn_revision.xlsx", sheet = 'PAC_DOR_anti') %>% 
  select(Group, Latency_1st_crossing, Jumping) %>% 
  mutate(Latency = (Jumping - Latency_1st_crossing)/10) %>% 
  mutate(Group = factor(Group, levels = c("Saline", "Dor", "Dor1", "Mor", "Mor1"))) %>% 
  ddply(., .(Group), summarise, mean= mean(Latency),se=sd(Latency)/sqrt(length(Latency)))

  aov(Latency~Group, .)


summary(t_DOR_jumping)
TukeyHSD(t_DOR_jumping)

p_DOR_pain <- plot_grid(p_DOR_licking, p_DOR_rearing, p_DOR_jumping, nrow = 1)

setwd("~cchen/Documents/neuroscience/Pn\ project/Figure/PDF/")
cairo_pdf("p_DOR_pain.pdf", width = 180/25.6, height = 65/25.6, family = "Arial")
p_DOR_pain
dev.off()  



## for noloxone injection during conditioning phase------
## combine the data from WT cond mice without injection, in the 'Pn_behavior_analysis' file


p_latency2 <- read.xlsx("~cchen/Documents/neuroscience/Pn\ project/Data_analysis/For\ revision/Pn_revision.xlsx", sheet = 'PAC_nalox_cond') %>% 
  drop_na() %>% 
  mutate(Latency2 = (Crossing_back - Latency_1st_crossing)/10) %>% 
  mutate(Latency1 = Crossing_back_d3/10) %>% 
  select(Group, ID, Latency1, Latency2) %>% 
  pivot_longer(-c(Group, ID)) %>% 
  rename(Day= name) %>% 
  mutate(Group = factor(Group, levels = c("Saline", "Nalox"))) %>% 
  ggplot(., aes(interaction(Day, Group), value, fill= Day))+
  geom_boxplot( outlier.shape = NA)+
  geom_line(aes(group=interaction(ID, Group)), colour="gray90")+
  geom_point(aes(colour = Group,shape = Group),position=position_jitterdodge(jitter.width = .4))+
  scale_fill_manual(values = c( "white", "white"))+
  scale_shape_manual(values=c( 17, 18))+
  scale_color_manual(values=c("darkorange", "deepskyblue3"))+
  labs(x="", y="Latency of crossing back (s)")+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  scale_y_continuous(limits = c(0, 150))+
  theme(legend.position = "none")

t_latency2 <- read.xlsx("~cchen/Documents/neuroscience/Pn\ project/Data_analysis/For\ revision/Pn_revision.xlsx", sheet = 'PAC_nalox_cond') %>% 
  drop_na() %>% 
  mutate(Latency2 = (Crossing_back - Latency_1st_crossing)/10) %>% 
  mutate(Latency1 = Crossing_back_d3/10) %>% 
  select(Group, ID, Latency1, Latency2) %>% 
  pivot_longer(-c(Group, ID)) %>% 
  filter(name =="Latency2") %>% 
  wilcox.test(value~ Group,.)



summary(t_latency2)
TukeyHSD(t_latency2, which = "Group")  
pairwise.t.test(t_latency2$value, t_latency2$name,p.adjust.method = "BH")


p_ratio <- read.xlsx("~cchen/Documents/neuroscience/Pn\ project/Data_analysis/For\ revision/Pn_revision.xlsx", sheet = 'PAC_nalox_cond') %>% 
  select(Group, ID, Time_hot, Time_hot3) %>% 
  pivot_longer(-c(Group, ID)) %>% 
  rename(Day= name) %>% 
  mutate(Day = factor(Day, levels = c("Time_hot3", "Time_hot"))) %>% 
  mutate(Group = factor(Group, levels = c("Saline", "Nalox"))) %>% 
  ggplot(., aes(interaction(Day, Group), value, fill= Day))+
  geom_boxplot( outlier.shape = NA)+
  geom_line(aes(group=interaction(ID, Group)), colour="gray90")+
  geom_point(aes(colour = Group,shape = Group),position=position_jitterdodge(jitter.width = .4))+
  scale_fill_manual(values = c("white", "white"))+
  scale_shape_manual(values=c( 17, 18))+
  scale_color_manual(values=c("darkorange", "deepskyblue3"))+
  labs(x="", y="Time spent in 2nd plate (%)")+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  scale_y_continuous(limits = c(0, 100), expand = c(0, 0))+
  theme(legend.position = "none")


t_ratio <- read.xlsx("~cchen/Documents/neuroscience/Pn\ project/Data_analysis/For\ revision/Pn_revision.xlsx", sheet = 'PAC_nalox_cond') %>% 
  select(Group, ID, Time_hot, Time_hot3) %>% 
  pivot_longer(-c(Group, ID)) %>% 
  filter(name =="Time_hot") %>% 
  wilcox.test(value~ Group,.)

summary(t_ratio)
TukeyHSD(t_ratio, which = "Group")  
pairwise.t.test(t_latency2$value, t_latency2$name,p.adjust.method = "BH")


p_licking <- read.xlsx("~cchen/Documents/neuroscience/Pn\ project/Data_analysis/For\ revision/Pn_revision.xlsx", sheet = 'PAC_nalox_cond') %>% 
  select(Group, Latency_1st_crossing, First_licking) %>% 
  mutate(Latency = (First_licking - Latency_1st_crossing)/10) %>% 
  dplyr::select(Group, Latency) %>% 
  mutate(Group = factor(Group, levels = c("Saline", "Nalox"))) %>% 
  ggplot(., aes(Group,Latency, fill=Group))+
  geom_boxplot(outlier.shape = NA)+
  geom_point(aes(colour = Group,shape = Group),position=position_jitterdodge(jitter.width = .4))+
  scale_fill_manual(values = c("white", "white", "white"))+
  scale_shape_manual(values=c( 17, 18))+
  scale_color_manual(values=c("darkorange", "deepskyblue3"))+
  labs(x="", y="Latency to 1st licking (s)")+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  scale_y_continuous(limits = c(0, 70), expand = c(0,0))+
  theme(legend.position = "none")

t_licking <- read.xlsx("~cchen/Documents/neuroscience/Pn\ project/Data_analysis/For\ revision/Pn_revision.xlsx", sheet = 'PAC_nalox_cond') %>% 
  select(Group, Latency_1st_crossing, First_licking) %>% 
  mutate(Latency = (First_licking - Latency_1st_crossing)/10) %>% 
  dplyr::select(Group, Latency) %>% 
  wilcox.test(Latency~Group, .)

dat_licking_sta <- read.xlsx("~cchen/Documents/neuroscience/Pn\ project/Data_analysis/For\ revision/Pn_revision.xlsx", sheet = 'PAC_nalox_cond') %>% 
  select(Group, Latency_1st_crossing, First_licking) %>% 
  mutate(Latency = (First_licking - Latency_1st_crossing)/10) %>% 
  dplyr::select(Group, Latency) %>% 
  bind_rows(., dat_licking_wt) %>% 
  ddply(.,.(Group), summarise, mean=mean(Latency), se=sd(Latency)/sqrt(length(Latency)))
  
  
  
summary(t_licking)
TukeyHSD(t_licking)  
t_licking <- read.xlsx("~cchen/Documents/neuroscience/Pn\ project/Data_analysis/For\ revision/Pn_revision.xlsx", sheet = 'PAC_nalox_cond') %>% 
  select(Group, Latency_1st_crossing, First_licking) %>% 
  mutate(Latency = (First_licking - Latency_1st_crossing)/10) %>% 
  mutate(Group = factor(Group, levels = c("Saline", "Nalox"))) %>% 
  wilcox.test(Latency~ Group, .)



p_rearing <- read.xlsx("~cchen/Documents/neuroscience/Pn\ project/Data_analysis/For\ revision/Pn_revision.xlsx", sheet = 'PAC_nalox_cond') %>% 
  select(Group, Latency_1st_crossing, First_rearing) %>% 
  mutate(Latency = (First_rearing - Latency_1st_crossing)/10) %>% 
  dplyr::select(Group, Latency) %>% 
  mutate(Group = factor(Group, levels = c("Saline", "Nalox"))) %>% 
  ggplot(., aes(Group,Latency, fill=Group))+
  geom_boxplot(outlier.shape = NA)+
  geom_point(aes(colour = Group,shape = Group),position=position_jitterdodge(jitter.width = .4))+
  scale_fill_manual(values = c("white", "white", "white"))+
  scale_shape_manual(values=c( 17, 18))+
  scale_color_manual(values=c("darkorange", "deepskyblue3"))+
  labs(x="", y="Latency to 1st rearing (s)")+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  scale_y_continuous(limits = c(0, 150), expand = c(0,0))+
  theme(legend.position = "none")


t_rearing <- read.xlsx("~cchen/Documents/neuroscience/Pn\ project/Data_analysis/For\ revision/Pn_revision.xlsx", sheet = 'PAC_nalox_cond') %>% 
  select(Group, Latency_1st_crossing, First_rearing) %>% 
  mutate(Latency = (First_rearing - Latency_1st_crossing)/10) %>% 
  dplyr::select(Group, Latency) %>% 
  wilcox.test(Latency~Group, .) 
summary(t_rearing)
TukeyHSD(t_rearing) 

dat_rearing_sta <- read.xlsx("~cchen/Documents/neuroscience/Pn\ project/Data_analysis/For\ revision/Pn_revision.xlsx", sheet = 'PAC_nalox_cond') %>% 
  select(Group, Latency_1st_crossing, First_rearing) %>% 
  mutate(Latency = (First_rearing - Latency_1st_crossing)/10) %>% 
  dplyr::select(Group, Latency) %>% 
  bind_rows(., dat_rearing_wt) %>% 
  ddply(.,.(Group), summarise, mean=mean(Latency), se=sd(Latency)/sqrt(length(Latency)))


t_rearing <- read.xlsx("~cchen/Documents/neuroscience/Pn\ project/Data_analysis/For\ revision/Pn_revision.xlsx", sheet = 'PAC_nalox_cond') %>% 
  select(Group, Latency_1st_crossing, First_rearing) %>% 
  mutate(Latency = (First_rearing - Latency_1st_crossing)/10) %>% 
  wilcox.test(Latency~Group,.)

  

p_jumping <- read.xlsx("~cchen/Documents/neuroscience/Pn\ project/Data_analysis/For\ revision/Pn_revision.xlsx", sheet = 'PAC_nalox_cond') %>% 
  select(Group, Latency_1st_crossing, Jumping) %>% 
  mutate(Latency = (Jumping - Latency_1st_crossing)/10) %>% 
  dplyr::select(Group, Latency) %>% 
  mutate(Group = factor(Group, levels = c("Saline", "Nalox"))) %>% 
  ggplot(., aes(Group,Latency, fill=Group))+
  geom_boxplot(outlier.shape = NA)+
  geom_point(aes(colour = Group,shape = Group),position=position_jitterdodge(jitter.width = .4))+
  scale_fill_manual(values = c("white", "white", "white"))+
  scale_shape_manual(values=c( 17, 18))+
  scale_color_manual(values=c("darkorange", "deepskyblue3"))+
  labs(x="", y="Latency to 1st jumping (s)")+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  scale_y_continuous(limits = c(50, 200))+
  theme(legend.position = "none")

t_jumping <- read.xlsx("~cchen/Documents/neuroscience/Pn\ project/Data_analysis/For\ revision/Pn_revision.xlsx", sheet = 'PAC_nalox_cond') %>% 
  select(Group, Latency_1st_crossing, Jumping) %>% 
  mutate(Latency = (Jumping - Latency_1st_crossing)/10) %>% 
  dplyr::select(Group, Latency) %>% 
  wilcox.test(Latency~Group, .) 

summary(t_jumping)
TukeyHSD(t_jumping) 


dat_jumping_sta <- read.xlsx("~cchen/Documents/neuroscience/Pn\ project/Data_analysis/For\ revision/Pn_revision.xlsx", sheet = 'PAC_nalox_cond') %>% 
  select(Group, Latency_1st_crossing, Jumping) %>% 
  mutate(Latency = (Jumping - Latency_1st_crossing)/10) %>% 
  dplyr::select(Group, Latency) %>% 
  bind_rows(., dat_jumping_wt) %>% 
  ddply(.,.(Group), summarise, mean=mean(Latency), se=sd(Latency)/sqrt(length(Latency)))

p_nalox_con1 <- plot_grid(p_latency2, p_ratio, nrow = 1)
p_nalox_con2 <- plot_grid( p_licking, p_rearing, p_jumping, nrow = 1)


setwd("~cchen/Documents/neuroscience/Pn\ project/Figure/PDF/")
cairo_pdf("p_nalox_con1.pdf", width = 100/25.6, height = 60/25.6, family = "Arial")
p_nalox_con1
dev.off()


setwd("~cchen/Documents/neuroscience/Pn\ project/Figure/PDF/")
cairo_pdf("p_nalox_con2.pdf", width = 105/25.6, height = 60/25.6, family = "Arial")
p_nalox_con2
dev.off()


## without Cond group
dat_wt_cond_latency <- dat_anti_wt %>% 
  dplyr::filter(variable=="latency2" & Group == "Cond.") %>% 
  mutate(Day = ifelse(Day == "Pre", "Latency1", "Latency2")) %>% 
  select(-variable) %>% 
  as_tibble() 


p_latency2 <- read.xlsx("~cchen/Documents/neuroscience/Pn\ project/Data_analysis/For\ revision/Pn_revision.xlsx", sheet = 'PAC_nalox_cond') %>% 
  drop_na() %>% 
  mutate(Latency2 = (Crossing_back - Latency_1st_crossing)/10) %>% 
  mutate(Latency1 = Crossing_back_d3/10) %>% 
  select(Group, ID, Latency1, Latency2) %>% 
  pivot_longer(-c(Group, ID)) %>% 
  rename(Day= name) %>% 
  bind_rows(., dat_wt_cond_latency) %>% 
  mutate(Group = factor(Group, levels = c("Cond.","Saline", "Nalox"))) %>% 
  filter(Group!= "Nalox") %>% 
  ggplot(., aes(interaction(Day, Group), value, fill= Day))+
  geom_boxplot( outlier.shape = NA)+
  geom_line(aes(group=interaction(ID, Group)), colour="gray90")+
  geom_point(aes(colour = Group,shape = Group),position=position_jitterdodge(jitter.width = .4))+
  scale_fill_manual(values = c("white", "white"))+
  scale_shape_manual(values=c(16, 17, 18))+
  scale_color_manual(values=c("indianred","darkorange"))+
  labs(x="", y="Latency of crossing back (s)")+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  scale_y_continuous(limits = c(0, 150))+
  theme(legend.position = "none")

t_latency2 <- read.xlsx("~cchen/Documents/neuroscience/Pn\ project/Data_analysis/For\ revision/Pn_revision.xlsx", sheet = 'PAC_nalox_cond') %>% 
  drop_na() %>% 
  mutate(Latency2 = (Crossing_back - Latency_1st_crossing)/10) %>% 
  mutate(Latency1 = Crossing_back_d3/10) %>% 
  select(Group, ID, Latency1, Latency2) %>% 
  pivot_longer(-c(Group, ID)) %>% 
  rename(Day= name) %>% 
  bind_rows(., dat_wt_cond_latency) %>% 
  aov(value ~ Day + Group,.)


summary(t_latency2)
TukeyHSD(t_latency2, which = "Group")  
pairwise.t.test(t_latency2$value, t_latency2$name,p.adjust.method = "BH")

dat_wt_cond_ratio <- dat_anti_wt %>% 
  filter(variable=="ratio" & Group =="Cond.") %>% 
  mutate(value=value*100) %>% 
  mutate(Day = ifelse(Day == "Pre", "Time_hot3", "Time_hot")) %>% 
  select(-variable) 


p_ratio <- read.xlsx("~cchen/Documents/neuroscience/Pn\ project/Data_analysis/For\ revision/Pn_revision.xlsx", sheet = 'PAC_nalox_cond') %>% 
  select(Group, ID, Time_hot, Time_hot3) %>% 
  pivot_longer(-c(Group, ID)) %>% 
  rename(Day= name) %>% 
  bind_rows(., dat_wt_cond_ratio) %>% 
  mutate(Day = factor(Day, levels = c("Time_hot3", "Time_hot"))) %>% 
  mutate(Group = factor(Group, levels = c("Cond.","Saline", "Nalox"))) %>% 
  filter(Group!= "Nalox") %>% 
  ggplot(., aes(interaction(Day, Group), value, fill= Day))+
  geom_boxplot( outlier.shape = NA)+
  geom_line(aes(group=interaction(ID, Group)), colour="gray90")+
  geom_point(aes(colour = Group,shape = Group),position=position_jitterdodge(jitter.width = .4))+
  scale_fill_manual(values = c("white", "white"))+
  scale_shape_manual(values=c(16, 17, 18))+
  scale_color_manual(values=c("indianred","darkorange"))+
  labs(x="", y="Time spent in 2nd plate (%)")+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  scale_y_continuous(limits = c(0, 100), expand = c(0, 0))+
  theme(legend.position = "none")


t_ratio <- read.xlsx("~cchen/Documents/neuroscience/Pn\ project/Data_analysis/For\ revision/Pn_revision.xlsx", sheet = 'PAC_nalox_cond') %>% 
  select(Group, ID, Time_hot, Time_hot3) %>% 
  pivot_longer(-c(Group, ID)) %>% 
  rename(Day= name) %>% 
  bind_rows(., dat_wt_cond_ratio) %>% 
  mutate(Day = factor(Day, levels = c("Time_hot3", "Time_hot"))) %>% 
  mutate(Group = factor(Group, levels = c("Cond.","Saline", "Nalox"))) %>% 
  aov(value ~ Day + Group,.)

summary(t_ratio)
TukeyHSD(t_ratio, which = "Group")  
pairwise.t.test(t_latency2$value, t_latency2$name,p.adjust.method = "BH")

dat_licking_wt <- dat_wt_manual %>% 
  filter(variable=="First_licking") %>% 
  dplyr::filter(Group =="Cond.") %>% 
  dplyr::select(Group, value) %>% 
  rename(Latency = value)

p_licking <- read.xlsx("~cchen/Documents/neuroscience/Pn\ project/Data_analysis/For\ revision/Pn_revision.xlsx", sheet = 'PAC_nalox_cond') %>% 
  select(Group, Latency_1st_crossing, First_licking) %>% 
  mutate(Latency = (First_licking - Latency_1st_crossing)/10) %>% 
  dplyr::select(Group, Latency) %>% 
  bind_rows(., dat_licking_wt) %>% 
  mutate(Group = factor(Group, levels = c("Cond.","Saline", "Nalox"))) %>% 
  filter(Group!= "Nalox") %>% 
  ggplot(., aes(Group,Latency, fill=Group))+
  geom_boxplot(outlier.shape = NA)+
  geom_point(aes(colour = Group,shape = Group),position=position_jitterdodge(jitter.width = .4))+
  scale_fill_manual(values = c("white", "white"))+
  scale_shape_manual(values=c(16, 17, 18))+
  scale_color_manual(values=c("indianred","darkorange"))+
  labs(x="", y="Latency to 1st licking (s)")+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  scale_y_continuous(limits = c(0, 70), expand = c(0,0))+
  theme(legend.position = "none")

t_licking <- read.xlsx("~cchen/Documents/neuroscience/Pn\ project/Data_analysis/For\ revision/Pn_revision.xlsx", sheet = 'PAC_nalox_cond') %>% 
  select(Group, Latency_1st_crossing, First_licking) %>% 
  mutate(Latency = (First_licking - Latency_1st_crossing)/10) %>% 
  dplyr::select(Group, Latency) %>% 
  bind_rows(., dat_licking_wt) %>% 
  aov(Latency~Group, .)

dat_licking_sta <- read.xlsx("~cchen/Documents/neuroscience/Pn\ project/Data_analysis/For\ revision/Pn_revision.xlsx", sheet = 'PAC_nalox_cond') %>% 
  select(Group, Latency_1st_crossing, First_licking) %>% 
  mutate(Latency = (First_licking - Latency_1st_crossing)/10) %>% 
  dplyr::select(Group, Latency) %>% 
  bind_rows(., dat_licking_wt) %>% 
  ddply(.,.(Group), summarise, mean=mean(Latency), se=sd(Latency)/sqrt(length(Latency)))

p_licking <- read.xlsx("~cchen/Documents/neuroscience/Pn\ project/Data_analysis/For\ revision/Pn_revision.xlsx", sheet = 'PAC_nalox_cond') %>% 
  select(Group, Latency_1st_crossing, First_licking) %>% 
  mutate(Latency = (First_licking - Latency_1st_crossing)/10) %>% 
  dplyr::select(Group, Latency) %>% 
  bind_rows(., dat_licking_wt) %>% 
  mutate(Group = factor(Group, levels = c("Cond.","Saline", "Nalox"))) %>% 
  filter(Group!= "Nalox") %>% 
  wilcox.test(Latency~Group,.)


summary(t_licking)
TukeyHSD(t_licking)  
t_licking <- read.xlsx("~cchen/Documents/neuroscience/Pn\ project/Data_analysis/For\ revision/Pn_revision.xlsx", sheet = 'PAC_nalox_cond') %>% 
  select(Group, Latency_1st_crossing, First_licking) %>% 
  mutate(Latency = (First_licking - Latency_1st_crossing)/10) %>% 
  mutate(Group = factor(Group, levels = c("Saline", "Nalox"))) %>% 
  wilcox.test(Latency~ Group, .)


dat_rearing_wt <- dat_wt_manual %>% 
  filter(variable=="First_rearing") %>% 
  dplyr::filter(Group =="Cond.") %>% 
  dplyr::select(Group, value) %>% 
  rename(Latency = value)

p_rearing <- read.xlsx("~cchen/Documents/neuroscience/Pn\ project/Data_analysis/For\ revision/Pn_revision.xlsx", sheet = 'PAC_nalox_cond') %>% 
  select(Group, Latency_1st_crossing, First_rearing) %>% 
  mutate(Latency = (First_rearing - Latency_1st_crossing)/10) %>% 
  dplyr::select(Group, Latency) %>% 
  bind_rows(., dat_rearing_wt) %>% 
  mutate(Group = factor(Group, levels = c("Cond.","Saline", "Nalox"))) %>% 
  filter(Group!= "Nalox") %>% 
  ggplot(., aes(Group,Latency, fill=Group))+
  geom_boxplot(outlier.shape = NA)+
  geom_point(aes(colour = Group,shape = Group),position=position_jitterdodge(jitter.width = .4))+
  scale_fill_manual(values = c("white", "white"))+
  scale_shape_manual(values=c(16, 17, 18))+
  scale_color_manual(values=c("indianred","darkorange"))+
  labs(x="", y="Latency to 1st rearing (s)")+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  scale_y_continuous(limits = c(0, 150), expand = c(0,0))+
  theme(legend.position = "none")


t_rearing <- read.xlsx("~cchen/Documents/neuroscience/Pn\ project/Data_analysis/For\ revision/Pn_revision.xlsx", sheet = 'PAC_nalox_cond') %>% 
  select(Group, Latency_1st_crossing, First_rearing) %>% 
  mutate(Latency = (First_rearing - Latency_1st_crossing)/10) %>% 
  dplyr::select(Group, Latency) %>% 
  bind_rows(., dat_rearing_wt) %>% 
  aov(Latency~Group, .) 
summary(t_rearing)
TukeyHSD(t_rearing) 

dat_rearing_sta <- read.xlsx("~cchen/Documents/neuroscience/Pn\ project/Data_analysis/For\ revision/Pn_revision.xlsx", sheet = 'PAC_nalox_cond') %>% 
  select(Group, Latency_1st_crossing, First_rearing) %>% 
  mutate(Latency = (First_rearing - Latency_1st_crossing)/10) %>% 
  dplyr::select(Group, Latency) %>% 
  bind_rows(., dat_rearing_wt) %>% 
  ddply(.,.(Group), summarise, mean=mean(Latency), se=sd(Latency)/sqrt(length(Latency)))


p_rearing <- read.xlsx("~cchen/Documents/neuroscience/Pn\ project/Data_analysis/For\ revision/Pn_revision.xlsx", sheet = 'PAC_nalox_cond') %>% 
  select(Group, Latency_1st_crossing, First_rearing) %>% 
  mutate(Latency = (First_rearing - Latency_1st_crossing)/10) %>% 
  dplyr::select(Group, Latency) %>% 
  bind_rows(., dat_rearing_wt) %>% 
  mutate(Group = factor(Group, levels = c("Cond.","Saline", "Nalox"))) %>% 
  filter(Group!= "Nalox")  %>% 
  wilcox.test(Latency~Group,.)

dat_jumping_wt <- dat_wt_manual %>% 
  filter(variable=="Jump") %>% 
  dplyr::filter(Group =="Cond.") %>% 
  dplyr::select(Group, value) %>% 
  rename(Latency = value) %>% 
  mutate(Latency = Latency/10)


p_jumping <- read.xlsx("~cchen/Documents/neuroscience/Pn\ project/Data_analysis/For\ revision/Pn_revision.xlsx", sheet = 'PAC_nalox_cond') %>% 
  select(Group, Latency_1st_crossing, Jumping) %>% 
  mutate(Latency = (Jumping - Latency_1st_crossing)/10) %>% 
  dplyr::select(Group, Latency) %>% 
  bind_rows(., dat_jumping_wt) %>% 
  mutate(Group = factor(Group, levels = c("Cond.","Saline", "Nalox"))) %>% 
  filter(Group!= "Nalox") %>% 
  ggplot(., aes(Group,Latency, fill=Group))+
  geom_boxplot(outlier.shape = NA)+
  geom_point(aes(colour = Group,shape = Group),position=position_jitterdodge(jitter.width = .4))+
  scale_fill_manual(values = c("white", "white"))+
  scale_shape_manual(values=c(16, 17, 18))+
  scale_color_manual(values=c("indianred","darkorange"))+
  labs(x="", y="Latency to 1st jumping (s)")+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  scale_y_continuous(limits = c(50, 200))+
  theme(legend.position = "none")

t_jumping <- read.xlsx("~cchen/Documents/neuroscience/Pn\ project/Data_analysis/For\ revision/Pn_revision.xlsx", sheet = 'PAC_nalox_cond') %>% 
  select(Group, Latency_1st_crossing, Jumping) %>% 
  mutate(Latency = (Jumping - Latency_1st_crossing)/10) %>% 
  dplyr::select(Group, Latency) %>% 
  bind_rows(., dat_jumping_wt) %>% 
  aov(Latency~Group, .) 
summary(t_jumping)
TukeyHSD(t_jumping) 

t_jumping <- read.xlsx("~cchen/Documents/neuroscience/Pn\ project/Data_analysis/For\ revision/Pn_revision.xlsx", sheet = 'PAC_nalox_cond') %>% 
  select(Group, Latency_1st_crossing, Jumping) %>% 
  mutate(Latency = (Jumping - Latency_1st_crossing)/10) %>% 
  dplyr::select(Group, Latency) %>% 
  bind_rows(., dat_jumping_wt) %>% 
  mutate(Group = factor(Group, levels = c("Cond.","Saline", "Nalox"))) %>% 
  filter(Group!= "Nalox") %>% 
  wilcox.test(Latency~Group,.)

dat_jumping_sta <- read.xlsx("~cchen/Documents/neuroscience/Pn\ project/Data_analysis/For\ revision/Pn_revision.xlsx", sheet = 'PAC_nalox_cond') %>% 
  select(Group, Latency_1st_crossing, Jumping) %>% 
  mutate(Latency = (Jumping - Latency_1st_crossing)/10) %>% 
  dplyr::select(Group, Latency) %>% 
  bind_rows(., dat_jumping_wt) %>% 
  ddply(.,.(Group), summarise, mean=mean(Latency), se=sd(Latency)/sqrt(length(Latency)))

p_nalox_con1 <- plot_grid(p_latency2, p_ratio, nrow = 1)
p_nalox_con2 <- plot_grid( p_licking, p_rearing, p_jumping, nrow = 1)


setwd("~cchen/Documents/neuroscience/Pn\ project/Figure/PDF/")
cairo_pdf("p_nalox_con1.pdf", width = 80/25.6, height = 60/25.6, family = "Arial")
p_nalox_con1
dev.off()


setwd("~cchen/Documents/neuroscience/Pn\ project/Figure/PDF/")
cairo_pdf("p_nalox_con2.pdf", width = 105/25.6, height = 60/25.6, family = "Arial")
p_nalox_con2
dev.off()
## compare ratio and latency when both plates are at 30 degree-----

p_latency_30 <- read.xlsx("~cchen/Documents/neuroscience/Pn\ project/Data_analysis/For\ revision/Pn_revision.xlsx", sheet = 'Post_30') %>% 
  mutate(Latency1 = Latency/10) %>% 
  mutate(Group = factor(Group, levels = c("Ctrl", "Cond"))) %>% 
  ggplot(., aes(Group,Latency1, fill=Group))+
  geom_boxplot(outlier.shape = NA)+
  geom_point(aes(colour = Group,shape = Group),position=position_jitterdodge(jitter.width = .5), size = 2)+
  scale_fill_manual(values = c("white", "white" ))+
  scale_shape_manual(values=c(20, 18))+
  scale_color_manual(values=c("darkcyan", "indianred"))+
  labs(x="", y="Latency of crossing back (s)")+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  scale_y_continuous(limits = c(0, 80))+
  theme(legend.position = "none")

t_latency_30 <- read.xlsx("~cchen/Documents/neuroscience/Pn\ project/Data_analysis/For\ revision/Pn_revision.xlsx", sheet = 'Post_30') %>% 
  mutate(Latency1 = Latency/10) %>% 
  wilcox.test(Latency1~ Group,.)
# ddply(.,.(Group), summarise, mean = mean(Latency1), se=sd(Latency1)/sqrt(length(Latency1))) 


p_ratio_30 <- read.xlsx("~cchen/Documents/neuroscience/Pn\ project/Data_analysis/For\ revision/Pn_revision.xlsx", sheet = 'Post_30') %>% 
  mutate(Group = factor(Group, levels = c("Ctrl", "Cond"))) %>% 
  ggplot(., aes(Group,Ratio, fill=Group))+
  geom_boxplot(outlier.shape = NA)+
  geom_point(aes(colour = Group,shape = Group),position=position_jitterdodge(jitter.width = .5), size = 2)+
  scale_fill_manual(values = c("white", "white" ))+
  scale_shape_manual(values=c(20, 18))+
  scale_color_manual(values=c("darkcyan", "indianred"))+
  labs(x="", y="Time spend in chamber 2 (%)")+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  scale_y_continuous(limits = c(0, 100))+
  theme(legend.position = "none")

t_ratio_30 <- read.xlsx("~cchen/Documents/neuroscience/Pn\ project/Data_analysis/For\ revision/Pn_revision.xlsx", sheet = 'Post_30') %>% 
  wilcox.test(Ratio ~ Group,.)
  
p_pain_30 <- plot_grid(p_latency_30, p_ratio_30)

setwd("~cchen/Documents/neuroscience/Pn\ project/Figure/PDF/")
cairo_pdf("p_pain_30.pdf", width = 60/25.6, height = 60/25.6, family = "Arial")
p_pain_30
dev.off()


## time profile of the ca2+ traces------
comp_time <- seq(-2, 1.9, by=0.1)
group_day <- c("Pre", "Cond", "Post")
score_range <- range(dat_cell_trace)

dat_mean_trace <- vector(mode = "list", length = 3)
for(i in c(1:3)){

  
  dat_mean_trace[[i]] <- lapply(dat_cell_trace, function(x)  x[[i]]) %>% 
    do.call(cbind,.) %>% 
    mutate(Time = comp_time) %>% 
    pivot_longer(-Time) %>% 
    mutate(Group = group_day[i])
    
  
}

dat_trace_mean_com <- dat_mean_trace %>% 
  do.call(rbind,.) %>% 
  ddply(.,.(Group, Time), summarise, mean= mean(value),se=sd(value)/sqrt(length(value))) %>% 
  ggplot(., aes(Time, mean, color= Group))+
  geom_line()+
  geom_ribbon(aes(ymin=mean-se, ymax=mean+se, fill=Group), alpha=0.2, linetype=0)+
  scale_color_manual(values=c("#8491B4FF", "#00A087FF", "#3C5488FF"))+
  labs(x="Time (s)", y="Ca2+ event rate (Hz)")+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  theme(legend.title = element_blank())

## plot the activity of rACC-Pn neurons 10s after crossing-----
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
    dat_stim1 <- dat_trace[(t1_p-40):(t1_p+200-1),] 
    
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
comp_time <- seq(-2, 9.9, by=0.1)
group_day <- c("Pre", "Cond", "Post")


dat_mean_trace <- vector(mode = "list", length = 3)
for(i in c(1:3)){
  
  
  dat_mean_trace[[i]] <- lapply(dat_cell_trace, function(x)  x[[i]]) %>% 
    do.call(cbind,.) %>% 
    mutate(Time = comp_time) %>% 
    pivot_longer(-Time) %>% 
    mutate(Group = group_day[i])
  
  
}

dat_trace_mean_com <- dat_mean_trace %>% 
  do.call(rbind,.) %>% 
  ddply(.,.(Group, Time), summarise, mean= mean(value),se=sd(value)/sqrt(length(value))) %>% 
  ggplot(., aes(Time, mean, color= Group))+
  geom_line()+
  geom_ribbon(aes(ymin=mean-se, ymax=mean+se, fill=Group), alpha=0.2, linetype=0)+
  scale_color_manual(values=c("#8491B4FF", "#00A087FF", "#3C5488FF"))+
  labs(x="Time (s)", y="Ca2+ event rate (Hz)")+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  theme(legend.title = element_blank(), legend.position = c(0.8, 0.9))

setwd("~cchen/Documents/neuroscience/Pn\ project/Figure/PDF/")
cairo_pdf("p_trace_mean_com.pdf", width = 90/25.6, height = 60/25.6, family = "Arial")
dat_trace_mean_com
dev.off()

## for histology quantification------
p_TRAP <- read.xlsx("~cchen/Documents/neuroscience/Pn\ project/Data_analysis/For\ revision/Pn_revision.xlsx", sheet = "TRAPed_neurons") %>% 
  as_tibble() %>% 
  mutate(Group = factor(Group, levels= c("Ctrl", "Cond"))) %>% 
  ggplot(., aes(Group, Ratio, Group = Group, fill = Group))+
  geom_boxplot(outlier.shape = NA, )+
  geom_point(aes(color = Group, shape = Group),position=position_jitterdodge(jitter.width = .2))+
  scale_fill_manual(values = c("white", "white"))+
  scale_shape_manual(values=c(20, 18))+
  scale_color_manual(values=c("darkcyan", "indianred"))+
  labs(x="", y="TRAPed neurons (#/mm^2)")+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  scale_y_continuous(limits = c(1, 40), expand = c(0, 0))+
  theme(legend.position = "none")

setwd("~cchen/Documents/neuroscience/Pn\ project/Figure/PDF/")
cairo_pdf("p_trap.pdf", width = 25/25.6, height = 60/25.6, family = "Arial")
p_TRAP
dev.off()

t_TRAP <- read.xlsx("~cchen/Documents/neuroscience/Pn\ project/Data_analysis/For\ revision/Pn_revision.xlsx", sheet = "TRAPed_neurons") %>% 
  as_tibble() %>% 
  mutate(Group = factor(Group, levels= c("Ctrl", "Cond"))) %>% 
  wilcox.test(Ratio~Group, .)






p_density <- read.xlsx("~cchen/Documents/neuroscience/Pn\ project/Data_analysis/For\ revision/Pn_revision.xlsx", sheet = "TRAP2_tracing_quantification_te") %>% 
  as_tibble() %>% 
  mutate(Group = factor(Group, levels= c("Ctrl", "Cond"))) %>% 
  mutate(Region = factor(Region, levels = c("CP", "ZI","MD", "VPM", "PN"))) %>% 
  ggplot(., aes(Region, Density,  Group = Group, fill = Group))+
  geom_boxplot(outlier.shape = NA, )+
  geom_point(aes(color = Group, shape = Group),position=position_jitterdodge(jitter.width = .2))+
  scale_fill_manual(values = c("white", "white"))+
  scale_shape_manual(values=c(20, 18))+
  scale_color_manual(values=c("darkcyan", "indianred"))+
  labs(x="", y="Axon Terminals (#/mm^2)")+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  scale_y_continuous(limits = c(1, 4000), expand = c(0, 0))+
  theme(legend.position = "none")


setwd("~cchen/Documents/neuroscience/Pn\ project/Figure/PDF/")
cairo_pdf("p_density.pdf", width = 50/25.6, height = 60/25.6, family = "Arial")
p_density
dev.off()

t_density <- read.xlsx("~cchen/Documents/neuroscience/Pn\ project/Data_analysis/For\ revision/Pn_revision.xlsx", sheet = "TRAP2_tracing_quantification_te") %>% 
  as_tibble() %>% 
  mutate(Group = factor(Group, levels= c("Ctrl", "Cond"))) %>% 
  mutate(Region = factor(Region, levels = c("CP", "ZI","MD", "VPM", "PN"))) %>% 
  aov(Density ~ Group * Area , data = .)
summary(t_density)
TukeyHSD(t_density)
  
## compare the activity of rACC-Pn neurons during rearing and licking-----
c_miniscope_matlab_ft <- function(file_trace, t_event) {
  ## import and format the data
  
  ID <- str_extract(file_trace, regex("m\\d+"))
  dat_trace1 <- raveio::read_mat(file_trace)
  num_compare <- c( 8)
  
  dat_trace <- dat_trace1[[8]] %>% 
      apply(., 2, scale)
    
    ## divide ca2+ trace within -2 and 6s before and after stimuls, take -2 to 0 as baseline
    stim_time<- seq(-2, 6.5, by=0.5)
    ## number of rows to be binned
    n <- 10 # 0.05*10=0.5
    t1_p <- t_event
    ## extract cell activity when they cross the border
    dat_stim1 <- dat_trace[(t1_p-40):(t1_p+140-1),] 
    dat_stim1 <- aggregate(dat_stim1,list(rep(1:(nrow(dat_stim1)%/%n+1),each=n,len=nrow(dat_stim1))),mean)[-1]
    dat_stim <- dat_stim1[1:5,] %>%  ## baseline as -2 to 0
      colMeans(., na.rm = T) %>% 
      sweep(dat_stim1, 2, ., FUN = "-")
    
    colnames(dat_stim) <- str_c(ID,"Cell", 1: ncol(dat_stim))

  return(dat_stim)
}

## analyze for all mouse 
mouse_file <- as.list(list.files("~cchen/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/FD-processed", full.names = T))
t_event_licking <- as.list(c(880, 431, 866, 2177, 350, 726))
t_event_rearing <- as.list(c(2372, 2158, 1646, 791, 984, 1138))

dat_cell_trace_licking <- mapply(c_miniscope_matlab_ft, mouse_file,t_event_licking, SIMPLIFY = F)


mouse_ID <- sort(c("m3", "m7", "m17", "m18", "m855", "m857"))

dat_cell_trace_d <- dat_cell_trace_licking %>% 
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
  
  rep_time <- dat_cell_trace_licking %>% 
    mapply(ncol, .)
  dat_cell_trace_d_re <- mutate(dat_cell_trace_d_re, Group= factor(Group, levels = c("Excited", "Neutral", "Inhibited"))) %>% 
    mutate(., ID = rep(mouse_ID, rep_time*length(stim_time)) )
  
  dat_cell_trace_re_licking <- dat_cell_trace_d_re

  

  
  score_range <- dat_cell_trace_re_licking %>% 
    .$value %>% 
    range()


    dat_trace <- dat_cell_trace_re_licking
    
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

    
    p_com_trace_licking <- dat_cell_trace_re_licking %>% 
      ddply(., .(Time, ID), summarise, value = mean(value)) %>% 
      ddply(., .(Time), summarise,mean=mean(value),sd=sd(value),se=sd(value)/sqrt(length(value))) %>%  
      ggplot(., aes(Time, mean))+
      geom_line()+
      geom_ribbon(aes(ymin=mean-se, ymax=mean+se), alpha=0.1, linetype=0)+
      theme(axis.line.x = element_line(),
            axis.line.y = element_line(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.border = element_blank(),
            panel.background = element_blank(),
            axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
      scale_y_continuous(limits = c(-0.1, 0.5)) +
      theme(legend.title = element_blank())
    
    
    
    dat_com_trace_licking <- dat_cell_trace_licking %>% 
      mapply(rowMeans, ., SIMPLIFY = F) %>% 
      do.call(cbind,.) %>% 
      as_tibble() %>% 
      rename()
    
    colnames(dat_com_trace_licking) <- mouse_ID
    
    p_com_trace_licking <- dat_com_trace_licking %>% 
      mutate(Time = stim_time) %>% 
      pivot_longer(-Time) %>% 
      mutate(Group = "Licking")
    
    ## for rearing
    
    dat_cell_trace_rearing <- mapply(c_miniscope_matlab_ft, mouse_file,t_event_rearing, SIMPLIFY = F)
    
    dat_com_trace_rearing <- dat_cell_trace_rearing %>% 
      mapply(rowMeans, ., SIMPLIFY = F) %>% 
      do.call(cbind,.) %>% 
      as_tibble() %>% 
      rename()
    
    colnames(dat_com_trace_rearing) <- mouse_ID
    
    p_com_trace_rearing_licking <- dat_com_trace_rearing %>% 
      mutate(Time = stim_time) %>% 
      pivot_longer(-Time) %>% 
      mutate(Group = "Rearing") %>% 
      bind_rows(., p_com_trace_licking) %>% 
      ddply(., .(Time, Group), summarise,mean=mean(value),sd=sd(value),se=sd(value)/sqrt(length(value))) %>%  
      ggplot(., aes(Time, mean, colour = Group))+
      geom_line()+
      geom_ribbon(aes(ymin=mean-se, ymax=mean+se), alpha=0.1, linetype=0)+
      theme(axis.line.x = element_line(),
            axis.line.y = element_line(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.border = element_blank(),
            panel.background = element_blank(),
            axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
      scale_y_continuous(limits = c(-0.2, 0.5)) +
      theme(legend.title = element_blank())
    
## Examine the OPRD1 neurons receive ACC input-----
    
    p_ratio_acc_oprd <- read.xlsx("~cchen/Documents/neuroscience/Pn\ project/Data_analysis/For\ revision/Pn_revision.xlsx", sheet = 'PAC_ACC_OPRD_OPTO') %>% 
      select(Group, ID, Time_hot_d3, Time_hot_d7) %>% 
      pivot_longer(-c(Group, ID)) %>% 
      mutate(name = factor(name, levels=c("Time_hot_d3", "Time_hot_d7"))) %>% 
      mutate(Group = factor(Group, levels = c("eYFP", "eNpHR"))) %>% 
      ggplot(., aes(interaction(name, Group), value, fill= name))+
      geom_boxplot( outlier.shape = NA)+
      geom_line(aes(group=interaction(ID, Group)), colour="gray90")+
      geom_point(aes(colour = Group,shape = Group),position=position_jitterdodge(jitter.width = .2))+
      scale_fill_manual(values = c("white", "white"))+
      scale_shape_manual(values=c(20, 18))+
      scale_color_manual(values=c("darkcyan", "indianred"))+
      labs(x="", y="Time spend in 2nd plate (%)")+
      theme(axis.line.x = element_line(),
            axis.line.y = element_line(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.border = element_blank(),
            panel.background = element_blank(),
            axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
      scale_y_continuous(limits = c(0, 100), expand = c(0, 0))+
      theme(legend.position = 'none', legend.title = element_blank())
    
    p_return_latency <- read.xlsx("~cchen/Documents/neuroscience/Pn\ project/Data_analysis/For\ revision/Pn_revision.xlsx", sheet = 'PAC_ACC_OPRD_OPTO') %>% 
      drop_na() %>% 
      mutate(Latency = (Crossing_back - Latency_1st_crossing)/10) %>% 
      mutate(Latency1 = Crossing_back_d3/10) %>% 
      select(Group, ID, Latency, Latency1) %>% 
      pivot_longer(-c(Group, ID)) %>% 
      mutate(name = factor(name, levels=c("Latency1", "Latency"))) %>% 
      mutate(Group = factor(Group, levels = c("eYFP", "eNpHR"))) %>% 
      ggplot(., aes(interaction(name, Group), value, fill= name))+
      geom_boxplot( outlier.shape = NA)+
      geom_line(aes(group=interaction(ID, Group)), colour="gray90")+
      geom_point(aes(colour = Group,shape = Group),position=position_jitterdodge(jitter.width = .2))+
      scale_fill_manual(values = c("white", "white"))+
      scale_shape_manual(values=c(20, 18))+
      scale_color_manual(values=c("darkcyan", "indianred"))+
      labs(x="", y="Latency of crossing back (s)")+
      theme(axis.line.x = element_line(),
            axis.line.y = element_line(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.border = element_blank(),
            panel.background = element_blank(),
            axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
      theme(legend.position = 'none', legend.title = element_blank())
    
    p_ACC_oprd_com <- plot_grid(p_ratio_acc_oprd, p_return_latency, nrow = 1)

    setwd("~cchen/Documents/neuroscience/Pn\ project/Figure/PDF/")
    cairo_pdf("p_ACC_oprd_com.pdf", width = 90/25.6, height =60/25.6, family = "Arial")
    p_ACC_oprd_com
    dev.off()   
    
    p_ACC_OPRD_licking <- read.xlsx("~cchen/Documents/neuroscience/Pn\ project/Data_analysis/For\ revision/Pn_revision.xlsx", sheet = 'PAC_ACC_OPRD_OPTO') %>% 
      select(Group, Latency_1st_crossing, First_licking) %>% 
      mutate(Latency = (First_licking - Latency_1st_crossing)/10) %>% 
      mutate(Group = factor(Group, levels = c("eYFP", "eNpHR"))) %>% 
      ggplot(., aes(Group,Latency, fill=Group))+
      geom_boxplot(outlier.shape = NA)+
      geom_point(aes(colour = Group,shape = Group),position=position_jitterdodge(jitter.width = .5), size = 2)+
      scale_fill_manual(values = c("white", "white" ))+
      scale_shape_manual(values=c(20, 18))+
      scale_color_manual(values=c("darkcyan", "indianred"))+
      labs(x="", y="Latency of 1st licking (s)")+
      theme(axis.line.x = element_line(),
            axis.line.y = element_line(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.border = element_blank(),
            panel.background = element_blank(),
            axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
      theme(legend.title = element_blank(), legend.position = "none")
    
    t_ACC_OPRD_licking <- read.xlsx("~cchen/Documents/neuroscience/Pn\ project/Data_analysis/For\ revision/Pn_revision.xlsx", sheet = 'PAC_ACC_OPRD_OPTO') %>% 
      select(Group, Latency_1st_crossing, First_licking) %>% 
      mutate(Latency = (First_licking - Latency_1st_crossing)/10) %>% 
      mutate(Group = factor(Group, levels = c("eYFP", "eNpHR"))) %>% 
      wilcox.test(Latency~Group, .)
    
    p_ACC_OPRD_rearing <- read.xlsx("~cchen/Documents/neuroscience/Pn\ project/Data_analysis/For\ revision/Pn_revision.xlsx", sheet = 'PAC_ACC_OPRD_OPTO') %>% 
      select(Group, Latency_1st_crossing, First_rearing) %>% 
      mutate(Latency = (First_rearing - Latency_1st_crossing)/10) %>% 
      mutate(Group = factor(Group, levels = c("eYFP", "eNpHR"))) %>% 
      ggplot(., aes(Group,Latency, fill=Group))+
      geom_boxplot(outlier.shape = NA)+
      geom_point(aes(colour = Group,shape = Group),position=position_jitterdodge(jitter.width = .5), size = 2)+
      scale_fill_manual(values = c("white", "white" ))+
      scale_shape_manual(values=c(20, 18))+
      scale_color_manual(values=c("darkcyan", "indianred"))+
      labs(x="", y="Latency of 1st rearing (s)")+
      theme(axis.line.x = element_line(),
            axis.line.y = element_line(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.border = element_blank(),
            panel.background = element_blank(),
            axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
      theme(legend.title = element_blank(), legend.position = "none")
    
    t_ACC_OPRD_rearing <- read.xlsx("~cchen/Documents/neuroscience/Pn\ project/Data_analysis/For\ revision/Pn_revision.xlsx", sheet = 'PAC_ACC_OPRD_OPTO') %>% 
      select(Group, Latency_1st_crossing, First_rearing) %>% 
      mutate(Latency = (First_rearing - Latency_1st_crossing)/10) %>% 
      mutate(Group = factor(Group, levels = c("eYFP", "eNpHR"))) %>% 
      wilcox.test(Latency~Group, .)
    
    
    p_ACC_OPRD_jumping <- read.xlsx("~cchen/Documents/neuroscience/Pn\ project/Data_analysis/For\ revision/Pn_revision.xlsx", sheet = 'PAC_ACC_OPRD_OPTO') %>% 
      select(Group, Latency_1st_crossing, Jumping) %>% 
      mutate(Latency = (Jumping - Latency_1st_crossing)/10) %>% 
      mutate(Group = factor(Group, levels = c("eYFP", "eNpHR"))) %>% 
      ggplot(., aes(Group,Latency, fill=Group))+
      geom_boxplot(outlier.shape = NA)+
      geom_point(aes(colour = Group,shape = Group),position=position_jitterdodge(jitter.width = .5), size = 2)+
      scale_fill_manual(values = c("white", "white" ))+
      scale_shape_manual(values=c(20, 18))+
      scale_color_manual(values=c("darkcyan", "indianred"))+
      labs(x="", y="Latency of 1st jumping (s)")+
      theme(axis.line.x = element_line(),
            axis.line.y = element_line(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.border = element_blank(),
            panel.background = element_blank(),
            axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
      theme(legend.title = element_blank(), legend.position = "none")
    
    
    t_ACC_OPRD_jumping <- read.xlsx("~cchen/Documents/neuroscience/Pn\ project/Data_analysis/For\ revision/Pn_revision.xlsx", sheet = 'PAC_ACC_OPRD_OPTO') %>% 
      select(Group, Latency_1st_crossing, Jumping) %>% 
      mutate(Latency = (Jumping - Latency_1st_crossing)/10) %>% 
      mutate(Group = factor(Group, levels = c("eYFP", "eNpHR"))) %>% 
      wilcox.test(Latency~Group, .)
    
    
    p_ACC_OPRD_pain <- plot_grid(p_ACC_OPRD_licking, p_ACC_OPRD_rearing, p_ACC_OPRD_jumping, nrow = 1)
    
    setwd("~cchen/Documents/neuroscience/Pn\ project/Figure/PDF/")
    cairo_pdf("p_ACC_OPRD_pain.pdf", width = 100/25.6, height = 60/25.6, family = "Arial")
    p_ACC_OPRD_pain
    dev.off()  
    
    
    ## For AAV1 wt to determine the Pn neurons receving ACC inputs----
    
    p_AAV1 <- read.xlsx("~cchen/Documents/neuroscience/Pn\ project/Data_analysis/For\ revision/Pn_revision.xlsx", sheet = 'AAV1') %>% 
      mutate(Ratio = Ratio* 100) %>% 
      mutate(DummyGroup = "All Data") %>%  # Add a dummy grouping variable
      ggplot(aes(x = DummyGroup, y = Ratio)) +  # Use the dummy grouping variable for x
      geom_boxplot(fill = "white") +
      geom_jitter(width = 0.2, size = 1.5, color = "black")+
      labs(x="", y="eGFP/DAPI (%)")+
      theme(axis.line.x = element_line(),
            axis.line.y = element_line(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.border = element_blank(),
            panel.background = element_blank(),
            axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
      scale_y_continuous(limits = c(20, 100))+
      theme(legend.title = element_blank(), legend.position = "none")

    setwd("~cchen/Documents/neuroscience/Pn\ project/Figure/PDF/")
    cairo_pdf("p_AAV1.pdf", width = 30/25.6, height = 70/25.6, family = "Arial")
    p_AAV1
    dev.off()  
    

    
    
## For formalin injection-----
    p_formalin <- read.xlsx("~cchen/Documents/neuroscience/Pn\ project/Data_analysis/For\ revision/Formalin_duration.xlsx", sheet= 1) %>% 
      as_tibble() %>% 
      mutate(Group = factor(Group, levels = c('Ctrl', 'Cond'))) %>% 
      ggplot(., aes(Group, attending, col= Group))+
      geom_boxplot(outlier.shape = NA)+
      geom_point(aes(colour = Group,shape = Group),position=position_jitterdodge(jitter.width = .5), size = 2)+
      scale_fill_manual(values = c("white", "white" ))+
      scale_shape_manual(values=c(20, 18))+
      scale_color_manual(values=c("darkcyan", "indianred"))+
      labs(x="", y="Attending duration (s)")+
      theme(axis.line.x = element_line(),
            axis.line.y = element_line(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.border = element_blank(),
            panel.background = element_blank(),
            axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
      scale_y_continuous(limits = c(0, 400))+
      theme(legend.title = element_blank(), legend.position = "none")
      
    setwd("~cchen/Documents/neuroscience/Pn\ project/Figure/PDF/")
    cairo_pdf("p_formalin.pdf", width = 42/25.6, height = 60/25.6, family = "Arial")
    p_formalin
    dev.off()  
  
    p_formalin_duration <- read.xlsx("~cchen/Documents/neuroscience/Pn\ project/Data_analysis/For\ revision/Formalin_duration.xlsx", sheet= 2) %>% 
      as_tibble() %>% 
      mutate(Group = factor(Group, levels = c('Ctrl', 'Cond'))) %>% 
      filter(Phase =="Early") %>% 
      ggplot(., aes(Group, Attending, col= Group))+
      geom_boxplot(outlier.shape = NA)+
      geom_point(aes(colour = Group,shape = Group),position=position_jitterdodge(jitter.width = .5), size = 2)+
      scale_fill_manual(values = c("white", "white" ))+
      scale_shape_manual(values=c(20, 18))+
      scale_color_manual(values=c("darkcyan", "indianred"))+
      labs(x="", y="Attending duration (s)")+
      theme(axis.line.x = element_line(),
            axis.line.y = element_line(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.border = element_blank(),
            panel.background = element_blank(),
            axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
      scale_y_continuous(limits = c(0, 120))+
      theme(legend.title = element_blank(), legend.position = "none")
    
    p_formalin_ratio <- read.xlsx("~cchen/Documents/neuroscience/Pn\ project/Data_analysis/For\ revision/Formalin_duration.xlsx", sheet= 2) %>% 
      as_tibble() %>% 
      mutate(Group = factor(Group, levels = c('Ctrl', 'Cond'))) %>% 
      filter(Phase =="Early") %>% 
      ggplot(., aes(Group, Ratio, col= Group))+
      geom_boxplot(outlier.shape = NA)+
      geom_point(aes(colour = Group,shape = Group),position=position_jitterdodge(jitter.width = .5), size = 2)+
      scale_fill_manual(values = c("white", "white" ))+
      scale_shape_manual(values=c(20, 18))+
      scale_color_manual(values=c("darkcyan", "indianred"))+
      labs(x="", y="Time spent in 2nd plate (%)")+
      theme(axis.line.x = element_line(),
            axis.line.y = element_line(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.border = element_blank(),
            panel.background = element_blank(),
            axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
      scale_y_continuous(limits = c(30, 120))+
      theme(legend.title = element_blank(), legend.position = "none")

    t_formalin_ratio <- read.xlsx("~cchen/Documents/neuroscience/Pn\ project/Data_analysis/For\ revision/Formalin_duration.xlsx", sheet= 2) %>% 
      as_tibble() %>% 
      mutate(Group = factor(Group, levels = c('Ctrl', 'Cond'))) %>% 
      filter(Phase =="Early") %>% 
      wilcox.test(Attending~Group,.)
    
    
    p_formalin_duration1 <- read.xlsx("~cchen/Documents/neuroscience/Pn\ project/Data_analysis/For\ revision/Formalin_duration.xlsx", sheet= 2) %>% 
      as_tibble() %>% 
      mutate(Group = factor(Group, levels = c('Ctrl', 'Cond'))) %>% 
      filter(Phase =="Total") %>% 
      ggplot(., aes(Group, Attending, col= Group))+
      geom_boxplot(outlier.shape = NA)+
      geom_point(aes(colour = Group,shape = Group),position=position_jitterdodge(jitter.width = .5), size = 2)+
      scale_fill_manual(values = c("white", "white" ))+
      scale_shape_manual(values=c(20, 18))+
      scale_color_manual(values=c("darkcyan", "indianred"))+
      labs(x="", y="Attending duration (s)")+
      theme(axis.line.x = element_line(),
            axis.line.y = element_line(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.border = element_blank(),
            panel.background = element_blank(),
            axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
      scale_y_continuous(limits = c(100, 400))+
      theme(legend.title = element_blank(), legend.position = "none")
    
    p_formalin_ratio1 <- read.xlsx("~cchen/Documents/neuroscience/Pn\ project/Data_analysis/For\ revision/Formalin_duration.xlsx", sheet= 2) %>% 
      as_tibble() %>% 
      mutate(Group = factor(Group, levels = c('Ctrl', 'Cond'))) %>% 
      filter(Phase =="Total") %>% 
      ggplot(., aes(Group, Ratio, col= Group))+
      geom_boxplot(outlier.shape = NA)+
      geom_point(aes(colour = Group,shape = Group),position=position_jitterdodge(jitter.width = .5), size = 2)+
      scale_fill_manual(values = c("white", "white" ))+
      scale_shape_manual(values=c(20, 18))+
      scale_color_manual(values=c("darkcyan", "indianred"))+
      labs(x="", y="Time spent in 2nd plate (%)")+
      theme(axis.line.x = element_line(),
            axis.line.y = element_line(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.border = element_blank(),
            panel.background = element_blank(),
            axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
      scale_y_continuous(limits = c(30, 120))+
      theme(legend.title = element_blank(), legend.position = "none")
    
    p_formalin_total <- plot_grid(p_formalin_ratio1, p_formalin_duration1, nrow = 1)
    
    setwd("~cchen/Documents/neuroscience/Pn\ project/Figure/PDF/")
    cairo_pdf("p_formalin_total.pdf", width = 60/25.6, height = 60/25.6, family = "Arial")
    p_formalin_total
    dev.off()
    
    
    p_formalin_early <- plot_grid(p_formalin_ratio, p_formalin_duration, nrow = 1)
    
    setwd("~cchen/Documents/neuroscience/Pn\ project/Figure/PDF/")
    cairo_pdf("p_formalin_early.pdf", width = 60/25.6, height = 60/25.6, family = "Arial")
    p_formalin_early
    dev.off()
## check Ca2+ activity during conditioning-------
    c_miniscope_matlab_ft <- function(file_trace) {
      ## import and format the data
      dat_trace1 <- raveio::read_mat(file_trace)
      ID <- str_extract(file_trace, regex("m\\d+"))
      dat_trace1[[2]][1] <- ifelse(ID=="m857", 1694, dat_trace1[[2]][1])
      num_compare <- c(4, 6, 7)
      t_crossing <- c(1, 3, 4)
      
      cross_ID <- dat_trace1$global_map %>% 
        as_tibble() %>% 
        select(V1, V3, V4, V5) %>% 
        mutate_all(na_if, 0) %>%
        drop_na() %>% 
        select(V1, V3, V4)
        
      
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
    group_day <- c("Pre","D5", "D6")
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
    
    p_heat <- plot_grid(p_heat_Pre,p_heat_D5, p_heat_D6, nrow  = 1)
    
    setwd("~cchen/Documents/neuroscience/Pn\ project/Figure/PDF/")
    cairo_pdf("p_heat.pdf", width = 140/25.6, height = 60/25.6, family = "Arial")
    p_heat
    dev.off()
    
    ## calcilate cell activity by individual cells
    group_day <- c("Pre","D5", "D6")
    
    dat_cell_trace_combine <- vector(mode = "list", length  = length(group_day))
    
    for (i in seq_along(group_day)){
      dat_cell_trace1 <- dat_cell_trace %>% 
        lapply(., function(x)  x[[i]]) %>% 
        do.call(cbind,.) %>% 
        add_column(Group = group_day[i]) %>% 
        add_column(Time = comp_time) %>% 
        pivot_longer(-c(Time, Group)) %>% 
        mutate(ID = str_extract(name, "m\\d+")) 
        
      dat_cell_trace_combine[[i]] <- dat_cell_trace1
    }

    p_cell_cond <- dat_cell_trace_combine %>% 
      do.call(rbind,.) %>% 
      ddply(., .(ID, Group, name), summarise, value1 = mean(value)) %>% 
      mutate(Group = factor(Group, levels = c("Pre","D5", "D6"))) %>% 
      ddply(., .(ID, Group), summarise, mean = mean(value1), se=sd(value1)/sqrt(length(value1))) %>% 
      ggplot(., aes(Group, mean, group = Group)) +
      geom_boxplot() +
      geom_jitter(aes(color = Group),position = position_dodge(0.1), size = 2)+
      labs(x="", y="Mean activity\n during first crossing (s.d.)")+
      theme(axis.line.x = element_line(),
            axis.line.y = element_line(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.border = element_blank(),
            panel.background = element_blank(),
            axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
      scale_y_continuous(limits = c(-0.15, 0.15))+
      theme(legend.title = element_blank(), legend.position = "none")
    
    p_cell_cond <- dat_cell_trace_combine %>% 
      do.call(rbind,.) %>% 
      ddply(., .(Time, ID, Group), summarise, value1 = mean(value)) %>% 
      mutate(Group = factor(Group, levels = c("Pre","D5", "D6"))) %>% 
      ddply(., .(Time, Group), summarise, mean = mean(value1), se=sd(value1)/sqrt(length(value1))) %>% 
      ggplot(., aes(Time, mean, group = Group, color = Group)) +
      geom_line()+
      #geom_ribbon(aes(ymin = mean - se, ymax = mean + se, fill = Group), alpha = 0.2)
      
      labs(x="", y="Mean activity\n during first crossing (s.d.)")+
      theme(axis.line.x = element_line(),
            axis.line.y = element_line(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.border = element_blank(),
            panel.background = element_blank(),
            axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
      scale_y_continuous(limits = c(-0.2, 0.15))+
      theme(legend.title = element_blank(), legend.position = "none")

    
    
    p_dat_anti <- dat_cell_trace_combine %>% 
      do.call(rbind,.) %>% 
      mutate(Group = factor(Group, levels = c("Pre","D5", "D6"))) %>% 
      ddply(.,.(ID, Group, name), summarise, cell_acti = mean(value)) %>% 
      ddply(.,.(ID, Group), summarise, cell_acti = mean(cell_acti)) %>% 
      ddply(., .(Group), summarise, mean = mean(cell_acti), se = sd(cell_acti) / sqrt(length(cell_acti))) %>% 
      ggplot(., aes(Group, mean, group = Group)) +
      geom_errorbar(aes(ymin = mean - se, ymax = mean + se), width = .1, color='black',position = position_dodge(0.1)) +
      geom_line(aes(group = 1), color='black') +  # Use a common group for all data points
      geom_point(aes(color = Group),position = position_dodge(0.1), size = 2)+
      labs(x="", y="Mean activity\n during first crossing (s.d.)")+
      theme(axis.line.x = element_line(),
            axis.line.y = element_line(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.border = element_blank(),
            panel.background = element_blank(),
            axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
      scale_y_continuous(limits = c(-0.2, 0.2))+
      theme(legend.title = element_blank(), legend.position = "none")
    
    setwd("~cchen/Documents/neuroscience/Pn\ project/Figure/PDF/")
    cairo_pdf("p_dat_anti.pdf", width = 60/25.6, height = 60/25.6, family = "Arial")
    p_dat_anti
    dev.off()    
    
    p_dat_anti <- dat_cell_trace_combine %>% 
      do.call(rbind,.) %>% 
      mutate(Group = factor(Group, levels = c("Pre","D5", "D6"))) %>% 
      ddply(.,.(ID, Group, name), summarise, cell_acti = mean(value)) %>% 
      ddply(.,.(ID, Group), summarise, cell_acti = mean(cell_acti)) %>% 
      filter(Group !="D5") %>% 
      wilcox.test(cell_acti~Group, ., paired = T, alternative = "l")
    
    
    ## calculate cell activity change mouse by mouse (cell becomes more active)
    cc_active_fun <- function (dat_trace){
      m_id <- dat_trace[[1]] %>% 
        colnames() %>%
        .[1] %>% 
        str_extract(., regex("m\\d+"))
      
      dat_trace_combine <- dat_trace %>% 
        do.call(rbind,.) %>% 
        as_tibble() %>% 
        add_column(Time = rep(comp_time, 4)) %>% 
        add_column(Group = rep(c("Pre","D4","D5", "D6"), each = length(comp_time))) %>% 
        pivot_longer(-c(Time, Group)) 
      
      dat_trace_acti <- dat_trace_combine %>% 
        ddply(., .(name, Group), summarise, mean_acti = mean(value)) %>% 
        ddply(., .(Group), summarise, cell_acti = mean(mean_acti)) %>% 
        add_column(ID = m_id)
      
      
      return(dat_trace_acti)
    }
    
    
    dat_acti_combine <- mapply(cc_active_fun, dat_cell_trace, SIMPLIFY = F)
    
    p_dat_anti <- dat_acti_combine %>% 
      do.call(rbind, .) %>% 
      mutate(Group = factor(Group, levels = c("Pre", "D5", "D6"))) %>% 
      ddply(., .(Group), summarise, mean = mean(cell_acti), se = sd(cell_acti) / sqrt(length(cell_acti))) %>% 
      ggplot(aes(Group, mean, group = Group)) +
      geom_errorbar(aes(ymin = mean - se, ymax = mean + se), width = .1, color='black',position = position_dodge(0.1)) +
      geom_line(aes(group = 1), color='black') +  # Use a common group for all data points
      geom_point(aes(color = Group),position = position_dodge(0.1), size = 2)+
      labs(x="", y="Mean activity\n during first crossing (s.d.)")+
      theme(axis.line.x = element_line(),
            axis.line.y = element_line(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.border = element_blank(),
            panel.background = element_blank(),
            axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
      scale_y_continuous(limits = c(-0.15, 0.15))+
      theme(legend.title = element_blank(), legend.position = "none")
    
    setwd("~cchen/Documents/neuroscience/Pn\ project/Figure/PDF/")
    cairo_pdf("p_dat_anti2.pdf", width = 60/25.6, height = 60/25.6, family = "Arial")
    p_dat_anti
    dev.off()    
    
    
    
    t_dat_anti <- dat_acti_combine %>% 
      do.call(rbind, .) %>% 
      mutate(Group = factor(Group, levels = c("Pre", "D5", "D6"))) %>% 
      aov(cell_acti~ Group, .)
    
    summary(t_dat_anti)
    
    
    ## plot D3, D5, and D6
    p_dat_anti <- dat_acti_combine %>% 
      do.call(rbind,.) %>% 
      mutate(Group = factor(Group, levels = c("Pre", "D4", "D5", "D6", "Post"))) %>% 
      dplyr::filter(!Group %in% c("D4", "Post")) %>% 
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
      scale_y_continuous(limits = c(-0.3, 0.6))+
      theme(legend.title = element_blank(), legend.position = "none")
    
    ## plot the trace to show Ca2+ activity during crossing
    cc_active_fun1 <- function (dat_trace){
      m_id <- dat_trace[[1]] %>% 
        colnames() %>%
        .[1] %>% 
        str_extract(., regex("m\\d+"))
      
      dat_trace_combine <- dat_trace %>% 
        do.call(rbind,.) %>% 
        as_tibble() %>% 
        add_column(Time = rep(comp_time, 4)) %>% 
        add_column(Group = rep(c("Pre", "D4","D5", "D6", "Post"), each = length(comp_time))) %>% 
        pivot_longer(-c(Time, Group)) %>% 
        ddply(.,.(Time, Group), summarise, mean = mean(value)) %>% 
        add_column(ID = m_id)
      
      return(dat_trace_combine)
    }
    
    
    dat_acti_combine1 <- mapply(cc_active_fun1, dat_cell_trace, SIMPLIFY = F)
    
    p_dat_anti1 <- dat_acti_combine1 %>% 
      do.call(rbind,.) %>% 
      filter(!Group %in% c("D4", "Post")) %>% 
      ddply(., .(Time, Group), summarise, mean1 = mean(mean), se=sd(mean)/sqrt(length(mean))) %>% 
      mutate(Group = factor(Group, levels = c("Pre", "D5", "D6"))) %>% 
      ggplot(., aes(Time, mean1, group = Group, color = Group))+
      geom_line()+
      # geom_ribbon(aes(ymin = mean1 - se, ymax = mean1 + se, fill = Group), alpha = 0.2)
      labs(x="Time (s)", y="Mean activity \nduring first crossing (s.d.)")+
      theme(axis.line.x = element_line(),
            axis.line.y = element_line(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.border = element_blank(),
            panel.background = element_blank(),
            axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
      scale_y_continuous(limits = c(-0.2, 0.1))+
      theme(legend.title = element_blank(), legend.position = c(0.1, 0.8))
    
    setwd("~cchen/Documents/neuroscience/Pn\ project/Figure/PDF/")
    cairo_pdf("p_dat_anti1.pdf", width = 90/25.6, height = 60/25.6, family = "Arial")
    p_dat_anti1
    dev.off()    
    
## Heatmap plot to show mouse behavior during block------
    dat_block_cond <- read.xlsx("~cchen/Documents/neuroscience/Pn\ project/Data_analysis/For\ revision/Block_tracking/PN_WT_naloxTrial61.xlsx") %>% 
      as_tibble() %>% 
      filter(In.zone ==1) %>% 
      select(X.center, Y.center) %>% 
      ggplot(., aes(x = X.center, y = Y.center))+
      stat_density2d(aes(fill = after_stat(level), alpha = after_stat(level)), geom = "polygon") +
      scale_fill_viridis_c() +  # Optional: Change to your preferred color scale
      theme_void() +
      guides(alpha = "none") 
      
      
 
      
      
    dat_block_ctrl <- read.xlsx("~cchen/Documents/neuroscience/Pn\ project/Data_analysis/For\ revision/Block_tracking/PN_WT_nalox_condTrial62.xlsx") %>% 
      as_tibble() %>% 
      filter(In.zone ==1) %>% 
      select(X.center, Y.center) %>% 
      ggplot(., aes(x = X.center, y = Y.center))+
      stat_density2d(aes(fill = after_stat(level), alpha = after_stat(level)), geom = "polygon") +
      scale_fill_viridis_c() +  # Optional: Change to your preferred color scale
      theme_void() +
      guides(alpha = "none") 

    p_block <- plot_grid(dat_block_ctrl, dat_block_cond, nrow = 1)      
    
    
    setwd("~cchen/Documents/neuroscience/Pn\ project/Figure/PDF/")
    cairo_pdf("p_block.pdf", width = 90/25.6, height = 60/25.6, family = "Arial")
    p_block
    dev.off()    
 
    # plot the visiting frequency
    p_door_visiting <- read.xlsx("~cchen/Documents/neuroscience/Pn\ project/Data_analysis/For\ revision/Pn_revision.xlsx", sheet = 'PAC_block') %>% 
      as_tibble() %>% 
      select(Group, Freq_back) %>% 
      mutate(Group = factor(Group, levels= c("Ctrl", "Cond"))) %>% 
      ggplot(., aes(x = Group, y = Freq_back , color = Group))+
      geom_boxplot(outlier.shape = NA)+
      geom_jitter(position=position_jitterdodge(), shape=1)+
      scale_color_manual(values=c("darkcyan", "indianred"))+
      labs(x="", y="Freq. of door visiting (#/min)")+
      theme(axis.line.x = element_line(),
            axis.line.y = element_line(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.border = element_blank(),
            panel.background = element_blank(),
            axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
      scale_y_continuous(limits = c(0, 8), expand = c(0, 0))+
      theme(legend.position = "none")
      
    p_door_visiting <- read.xlsx("~cchen/Documents/neuroscience/Pn\ project/Data_analysis/For\ revision/Pn_revision.xlsx", sheet = 'PAC_block') %>% 
      as_tibble() %>% 
      select(Group, Freq_back) %>% 
      mutate(Group = factor(Group, levels= c("Ctrl", "Cond"))) %>% 
      ggplot(., aes(x = Group, y = Freq_back , color = Group))+
      geom_boxplot(outlier.shape = NA)+
      geom_jitter(position=position_jitterdodge(), shape=1)+
      scale_color_manual(values=c("darkcyan", "indianred"))+
      labs(x="", y="Freq. of door visiting (#/min)")+
      theme(axis.line.x = element_line(),
            axis.line.y = element_line(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.border = element_blank(),
            panel.background = element_blank(),
            axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
      scale_y_continuous(limits = c(0, 8), expand = c(0, 0))+
      theme(legend.position = "none")
    
    t_door_visiting <- read.xlsx("~cchen/Documents/neuroscience/Pn\ project/Data_analysis/For\ revision/Pn_revision.xlsx", sheet = 'PAC_block') %>% 
      as_tibble() %>% 
      select(Group, Freq_back) %>% 
      wilcox.test(Freq_back~Group,.)
    
    p_freq_rearing <- read.xlsx("~cchen/Documents/neuroscience/Pn\ project/Data_analysis/For\ revision/Pn_revision.xlsx", sheet = 'PAC_block') %>% 
      as_tibble() %>% 
      select(Group, Freq_rearing) %>% 
      mutate(Group = factor(Group, levels= c("Ctrl", "Cond"))) %>% 
      ggplot(., aes(x = Group, y = Freq_rearing, color = Group))+
      geom_boxplot(outlier.shape = NA)+
      geom_jitter(position=position_jitterdodge(), shape=1)+
      scale_color_manual(values=c("darkcyan", "indianred"))+
      labs(x="", y="Freq. of rearing (#/min)")+
      theme(axis.line.x = element_line(),
            axis.line.y = element_line(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.border = element_blank(),
            panel.background = element_blank(),
            axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
      scale_y_continuous(limits = c(0, 12), expand = c(0, 0))+
      theme(legend.position = "none")
    
    t_freq_rearing <- read.xlsx("~cchen/Documents/neuroscience/Pn\ project/Data_analysis/For\ revision/Pn_revision.xlsx", sheet = 'PAC_block') %>% 
      as_tibble() %>% 
      select(Group, Freq_rearing) %>% 
      wilcox.test(Freq_rearing~Group,.)
    
    p_block <- plot_grid(p_door_visiting, p_freq_rearing, nrow = 1)
      
    setwd("~cchen/Documents/neuroscience/Pn\ project/Figure/PDF/")
    cairo_pdf("p_block.pdf", width = 70/25.6, height = 60/25.6, family = "Arial")
    p_block
    dev.off()    
    
  ## clustering analysis of rACC-Pn neurons------
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
    
    dat_cell_trace_combine <- lapply(dat_cell_trace, function(x)  do.call(rbind, x)) %>% 
      do.call(cbind,.)
    
    library(factoextra)
    setwd("~cchen/Documents/neuroscience/Pn\ project/Figure/PDF/")
    cairo_pdf("p_cluster_opto.pdf", width = 80/25.6, height = 60/25.6, family = "Arial")
    fviz_nbclust(t(dat_cell_trace_combine), kmeans, method = "silhouette")
    dev.off()   
    fviz_nbclust(t(dat_cell_trace_combine), kmeans, method = "silhouette")
    dat_cell_trace_cluster <- kmeans(t(dat_cell_trace_combine), centers = 3, iter.max = 10)[[1]]
    
    cluster1_name <- names(dat_cell_trace_cluster[dat_cell_trace_cluster ==2])
    
    dat_cell_average_cluster1 <- dat_cell_trace_combine %>% 
      select(all_of(cluster1_name)) %>%
      rowMeans() 
    
    dat_cluster_trace <- c()
    
    for (i in c(1:3)){
      cluster1_name <- names(dat_cell_trace_cluster[dat_cell_trace_cluster ==i])
      dat_cell_value <- dat_cell_trace_combine %>% 
        select(all_of(cluster1_name)) %>% 
        as_tibble() %>% 
        mutate(Cluster = str_c("C", i),Group = rep(c("Pre",  "Cond", "Post"), each = 40), Time = rep(seq(-2, 1.9, by=0.1), 3)) %>% 
        pivot_longer(-c(Time, Cluster, Group))
      
      dat_cluster_trace <- rbind(dat_cluster_trace, dat_cell_value)
      
    }

    
    dat_cell_cluster_order <- dat_cell_trace_combine %>% 
      slice(1:40) %>% 
      as_tibble() %>% 
      add_column(Time = seq(-2, 1.9, by=0.1)) %>% 
      pivot_longer(-Time) %>% 
      add_column(Cluster = rep(str_c("C", dat_cell_trace_cluster), 40)) %>% 
      ddply(.,.(Cluster), summarise, mean = mean(value)) %>% 
      arrange(mean)
    
    
    dat_cell_order <- dat_cell_trace_combine %>% 
      slice(1:40) %>% 
      as_tibble() %>% 
      add_column(Time = seq(-2, 1.9, by=0.1)) %>% 
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
      mutate(Group = factor(Group, levels = c("Pre", "Cond", "Post"))) %>% 
      ggplot(., aes(Time, mean, group = Cluster,color= Cluster))+
      geom_line()+
      geom_ribbon(aes(ymin=mean-se, ymax=mean+se, fill=Cluster), alpha=0.1, linetype=0)+
      scale_color_manual(values=c("#8491B4FF", "#00A087FF", "#3C5488FF"))+
      facet_wrap(~Group, nrow = 1)+
      theme_bw()+
      theme(legend.title = element_blank())+
      annotate(x=c(1.5,1.5,2,2), y=c(-1),"path")+
      annotate(x=c(2), y=c(-1, -1, 0, 0),"path")    
    
    setwd("~cchen/Documents/neuroscience/Pn\ project/Figure/PDF/")
    cairo_pdf("p_trace_cluster.pdf", width = 120/25.6, height = 60/25.6, family = "Arial")
    p_cluster
    dev.off()    
    
    ## heatmap plot of three groups
    score_range <- dat_cell_trace_combine %>% 
      range()
    
    
    dat_rect <- dat_cell_order %>% 
      ddply(., .(Cluster), summarise, n = length(name)) %>% 
      add_column(xmin = 2, xmax = 2.1) %>% 
      mutate(ymax = cumsum(n) ) %>% 
      mutate(ymin = c(1, lag(ymax)[-1]))
    
    p_pro_cluster<- ggplot(dat_rect, aes(x="", y=n, fill=Cluster)) +
      geom_bar(stat="identity", width=1, color="white") +
      coord_polar("y", start=0) +
      theme_void()
    
    
    p_heat_com <- dat_cell_trace_combine %>% 
      as_tibble() %>% 
      add_column(Time = rep(seq(-2, 1.9, by=0.1), 3))%>%
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
      annotate("rect", xmin = 2, xmax = 2.1, ymin = dat_rect$ymin, ymax = dat_rect$ymax-1, alpha = .2)

    
    setwd("~cchen/Documents/neuroscience/Pn\ project/Figure/PDF/")
    cairo_pdf("p_heat_com_cluster.pdf", width = 120/25.6, height = 60/25.6, family = "Arial")
    p_heat_com
    dev.off()    
    
    
## clustering neurons for each day-----

    dat_cell_trace_combine1 <- vector("list", 3)
    
    # Loop through each of the three sub-sublists
    for (i in 1:3) {
      dat_cell_trace_combine1[[i]] <- do.call(cbind, lapply(dat_cell_trace, function(x) x[[i]]))
    }
    
    # dat_cell_trace_combine1 now has three sublists, each a combination of the respective sublists from the original
    

    library(factoextra)
    
    cairo_pdf("p_cluster_opto.pdf", width = 80/25.6, height = 60/25.6, family = "Arial")
    fviz_nbclust(t(dat_cell_trace_combine1[[3]]), kmeans, method = "silhouette")
    dev.off()   
    
    c_cluster_fun <- function(dataset){
      
      k = 2
      dat_cell_trace_cluster <- kmeans(t(dataset), centers = k, iter.max = 10)[[1]]
      stim_time<- seq(-2, 1.9, by=0.1)
      
      dat_cell_trace_d_re<- dataset %>% 
        mutate(., Time= stim_time) %>% 
        melt(., id.vars='Time') %>% 
        mutate(., Group= rep(dat_cell_trace_cluster, each=length(stim_time)) )
      
      ## sort data by the value in each group
      dat_cell_d_sort <- as.numeric(names(sort(tapply(dat_cell_trace_d_re$value, dat_cell_trace_d_re$Group, mean), decreasing = T)))
      
      dat_cell_trace_d_re$Group[dat_cell_trace_d_re$Group == dat_cell_d_sort[1]] ="Excited"
      dat_cell_trace_d_re$Group[dat_cell_trace_d_re$Group == dat_cell_d_sort[2]] ="Neutral"
      dat_cell_trace_d_re$Group[dat_cell_trace_d_re$Group == dat_cell_d_sort[3]] ="Inhibited"
      

      return(dat_cell_trace_d_re)
 
    }
    
    group_day <- c("Pre", "Cond", "Test")
    
    dat_cell_trace_combine_cluster <- mapply(c_cluster_fun, dat_cell_trace_combine1, SIMPLIFY = F) %>% 
      do.call(rbind,.) %>% 
      mutate(Day = rep(group_day, each = nrow(.)/3))
    
    
    
    
    score_range <- range(dat_cell_trace_combine_cluster$value)
    
    group_day <- c("Pre", "Cond", "Test")
    
    for (i in c(1: 3)) {
      dat_trace <- dat_cell_trace_combine_cluster %>% 
        filter(Day == group_day[i])
      
      dat_rect <- dat_trace %>% 
        ddply(., .(Group), summarise, n = length(variable)/40) %>% 
        add_column(xmin = 2, xmax = 2.1) %>% 
        mutate(ymax = cumsum(n) ) %>% 
        mutate(ymin = c(1, lag(ymax)[-1]))
      
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
        theme(legend.position = "none")+
        annotate("rect", xmin = 2, xmax = 2.1, ymin = dat_rect$ymin, ymax = dat_rect$ymax-1, alpha = .2)
      
      assign(str_c("p_heat_", group_day[i]), p_heat)
    }
    
    ## combine the heat plot
    p_heat_com <- plot_grid(p_heat_Pre, p_heat_Cond, p_heat_Test, nrow = 1)
    
    setwd("~cchen/Documents/neuroscience/Pn\ project/Figure/PDF/")
    cairo_pdf("p_heat_com_cluster.pdf", width = 120/25.6, height = 60/25.6, family = "Arial")
    p_heat_com
    dev.off()    
    
    p_trace <- dat_cell_trace_combine_cluster %>% 
      ddply(.,.(Time, Day, Group), summarise,n=length(value),mean=mean(value),sd=sd(value),se=sd(value)/sqrt(length(value))) %>% 
      mutate(Day = factor(Day, levels = c("Pre", "Cond", "Test"))) %>% 
      mutate(Group = factor(Group, levels = c("Excited", "Neutral", "Inhibited"))) %>% 
      ggplot(., aes(Time, mean, colour=Group))+
      geom_line()+
      geom_ribbon(aes(ymin=mean-se, ymax=mean+se, fill=Group), alpha=0.1, linetype=0)+
      facet_grid(cols = vars(Day))+
      scale_colour_manual(values=c( "darkred", "navy"))+
      scale_fill_manual(values=c( "darkred", "navy"))+
      labs(x="Time relative to crossing (s)", y="Z score")+
      theme(axis.line.x = element_line(),
            axis.line.y = element_line(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.border = element_blank(),
            panel.background = element_blank(),
            axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
      geom_vline(xintercept = 0, col= 'red', linetype=2)+
      theme(legend.title = element_blank(), legend.position = 'none')
    
    setwd("~cchen/Documents/neuroscience/Pn\ project/Figure/PDF/")
    cairo_pdf("p_trace.pdf", width = 90/25.6, height = 60/25.6, family = "Arial")
    p_trace
    dev.off() 
  ## for clustering analysis of the mechanical pain data------
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
        mutate(name = str_c(ID, name)) %>% 
        pivot_wider(names_from = name, values_from = value1, values_fill = list(value1 = NA)) %>% 
        select(-Time)
      
      return(dat_trace_average1)
      
    }
    
    
    
    file_trace_pin <- as.list(list.files("~cchen/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/Pin_extract/pinprick/Pin/", full.names = T))
    
    stim_t_list <- list(m3 = c(337,665,916,1303,1541,1767,2032,2277,2496,2824),
                        m7 = c(493, 1323, 1867, 2415, 3035, 3556, 4079, 4519, 5020, 5517) ,
                        m16 = c(115, 414, 641,965, 1272, 1597, 1879, 2130, 2438,2679)*2,
                        m17 = c(262, 497, 803, 1057, 1404, 1676, 2020, 2248, 2543, 2805 )*2,
                        m18 = c(139, 681, 1012, 1365, 1663, 1900, 2110, 2324, 2566, 2769)*2)
    
    dat_pin_trace <- mapply(c_miniscope_matlab, file_trace_pin, SIMPLIFY = F) %>% 
      do.call(cbind,.)
      
    setwd("~cchen/Documents/neuroscience/Pn\ project/Figure/PDF/")
    cairo_pdf("p_cluster_opto.pdf", width = 80/25.6, height = 60/25.6, family = "Arial")
    fviz_nbclust(t(dat_pin_trace), kmeans, method = "silhouette")
    dev.off()   
    fviz_nbclust(t(dat_pin_trace), kmeans, method = "silhouette")

    dat_cell_trace_cluster <- kmeans(t(dat_pin_trace), centers = 2, iter.max = 10)[[1]]
    stim_time <- seq(-0.95, 2, by=0.05)
    
    dat_cell_pin_value <- dat_pin_trace %>% 
      as_tibble() %>% 
      mutate(Time = stim_time) %>% 
      pivot_longer(-Time) %>% 
      mutate(Cluster = rep(dat_cell_trace_cluster, 60)) %>% 
      mutate(Group = "Pin")
      
    mean_values <- dat_cell_pin_value %>%
      group_by(Cluster) %>%
      summarise(mean_value = mean(value, na.rm = TRUE)) %>%
      arrange(desc(mean_value))
    
    # Determine which Cluster has the larger mean
    larger_mean_cluster <- mean_values$Cluster[1]
    
    # Assign "Excitation" or "Inhibition" to Clusters
    dat_cell_pin_value <- dat_cell_pin_value %>%
      mutate(Cluster = ifelse(Cluster == larger_mean_cluster, "Excitation", "Inhibition"))
    
    ## for light touch as control------
    file_trace_touch <- as.list(list.files("~cchen/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/Pin_extract/pinprick/ctrl/", full.names = T))
    
    stim_t_list <- list(m3 = c(134, 528, 820, 1084, 1261, 1493, 1729, 1964, 2142,2367)*2,
                        m7 = c(195,594, 901, 1163, 1536, 1818, 2137, 2450, 2624, 3008 )*2 ,
                        m16 = c(192, 498, 779, 1068, 1359, 1624, 1918, 2209, 2592, 2858)*2 ,
                        m17 = c(216, 605, 982, 1198, 1604, 1829, 2098, 2348, 2520, 2743)*2,
                        m18 = c(225, 465, 883, 1385, 1561, 1735, 1897, 2085, 2370, 2691 ) *2)
    
    dat_touch_trace <- mapply(c_miniscope_matlab, file_trace_touch, SIMPLIFY = F) %>% 
      do.call(cbind,.)
    
    
    fviz_nbclust(t(dat_touch_trace), kmeans, method = "silhouette")
    
    dat_cell_trace_cluster <- kmeans(t(dat_touch_trace), centers = 2, iter.max = 10)[[1]]
    stim_time <- seq(-0.95, 2, by=0.05)
    
    dat_cell_touch_value <- dat_touch_trace %>% 
      as_tibble() %>% 
      mutate(Time = stim_time) %>% 
      pivot_longer(-Time) %>% 
      mutate(Cluster = rep(dat_cell_trace_cluster, 60)) %>% 
      mutate(Group = "Touch")
    
    mean_values <- dat_cell_touch_value %>%
      group_by(Cluster) %>%
      summarise(mean_value = mean(value, na.rm = TRUE)) %>%
      arrange(desc(mean_value))
    
    # Determine which Cluster has the larger mean
    larger_mean_cluster <- mean_values$Cluster[1]
    
    # Assign "Excitation" or "Inhibition" to Clusters
    dat_cell_touch_value <- dat_cell_touch_value %>%
      mutate(Cluster = ifelse(Cluster == larger_mean_cluster, "Excitation", "Inhibition"))
    
  ## for hargraveas
    file_trace_har <- as.list(list.files("~cchen/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/Pin_extract/Har/", pattern = ".mat",full.names = T))
    
    stim_t_list <- read.xlsx("~cchen/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/har_06112020/2020_06_01_har_events.xlsx") %>% 
      select(Mice, Frame) %>% 
      group_by(Mice) %>% 
      group_split(., .keep = F) 
    names(stim_t_list) <- sort(c("m3", "m7", "m16", "m17", "m18"))
    
    dat_har_trace <- mapply(c_miniscope_matlab, file_trace_har, SIMPLIFY = F) %>% 
      do.call(cbind,.)
    
    
    fviz_nbclust(t(dat_har_trace), kmeans, method = "silhouette")
    
    dat_cell_trace_cluster <- kmeans(t(dat_har_trace), centers = 2, iter.max = 10)[[1]]
    stim_time <- seq(-0.95, 2, by=0.05)
    
    dat_cell_har_value <- dat_har_trace %>% 
      as_tibble() %>% 
      mutate(Time = stim_time) %>% 
      pivot_longer(-Time) %>% 
      mutate(Cluster = rep(dat_cell_trace_cluster, 60)) %>% 
      mutate(Group = "Har")
    
    # Calculate mean values for each Cluster
    mean_values <- dat_cell_har_value %>%
      group_by(Cluster) %>%
      summarise(mean_value = mean(value, na.rm = TRUE)) %>%
      arrange(desc(mean_value))
    
    # Determine which Cluster has the larger mean
    larger_mean_cluster <- mean_values$Cluster[1]
    
    # Assign "Excitation" or "Inhibition" to Clusters
    dat_cell_har_value <- dat_cell_har_value %>%
      mutate(Cluster = ifelse(Cluster == larger_mean_cluster, "Excitation", "Inhibition"))
    
    # The Cluster column now contains "Excitation" or "Inhibition" based on the mean values
    
    
    p_cluster <- dat_cell_touch_value %>% 
      bind_rows(., dat_cell_pin_value, dat_cell_har_value) %>% 
      ddply(., .(Time, Cluster, Group), summarise,n=length(value),mean=mean(value),sd=sd(value),se=sd(value)/sqrt(length(value))) %>%
      mutate(Cluster = as.factor(Cluster)) %>% 
      mutate(Group = factor(Group, levels = c("Touch", "Pin", "Har"))) %>% 
      ggplot(., aes(Time, mean, group = Group,color= Group))+
      geom_line()+
      geom_ribbon(aes(ymin=mean-se, ymax=mean+se, fill=Group), alpha=0.1, linetype=0)+
      #scale_color_manual(values=c("#8491B4FF", "#00A087FF", "#3C5488FF"))+
      facet_wrap(~Cluster, nrow = 1)+
      theme_bw()+
      theme(legend.title = element_blank())
    
    t_cluster_excitation <- dat_cell_touch_value %>% 
      bind_rows(., dat_cell_pin_value, dat_cell_har_value) %>% 
      ddply(., .(Cluster, Group, name), summarise,value=mean(value)) %>%
      filter(Cluster == "Excitation") %>% 
      aov(value~Group, .)
    
    summary(t_cluster_excitation)
    TukeyHSD(t_cluster_excitation)
    
    p_cluster_inhibition <- dat_cell_touch_value %>% 
      bind_rows(., dat_cell_pin_value, dat_cell_har_value) %>% 
      ddply(., .(Cluster, Group, name), summarise,value=mean(value)) %>%
      mutate(Group = factor(Group, levels = c("Touch", "Pin","Har"))) %>% 
      ggplot(., aes(Group, value, color = Group))+
      geom_violin()+
      geom_jitter(width = 0.2, alpha = 0.5)+
      facet_wrap(~Cluster, nrow = 1)+
      theme(axis.line.x = element_line(),
            axis.line.y = element_line(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.border = element_blank(),
            panel.background = element_blank(),
            axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
      scale_y_continuous(limits = c(-1.5,1.5))+
      theme(legend.title = element_blank(), legend.position = "none")
      
    setwd("~cchen/Documents/neuroscience/Pn\ project/Figure/PDF/")
    cairo_pdf("p_cluster_inhibition.pdf", width = 80/25.6, height = 60/25.6, family = "Arial")
    p_cluster_inhibition
    dev.off()  
    
    t_cluster_inhibition <- dat_cell_touch_value %>% 
      bind_rows(., dat_cell_pin_value, dat_cell_har_value) %>% 
      ddply(., .(Cluster, Group, name), summarise,value=mean(value)) %>%
      filter(Cluster == "Inhibition") %>% 
      aov(value~Group, .)
    
    summary(t_cluster_inhibition)
    TukeyHSD(t_cluster_inhibition)

      mutate(Cluster = as.factor(Cluster)) %>% 
      mutate(Group = factor(Group, levels = c("Touch", "Pin", "Har"))) %>% 
    
    setwd("~cchen/Documents/neuroscience/Pn\ project/Figure/PDF/")
    cairo_pdf("p_trace_cluster.pdf", width = 80/25.6, height = 60/25.6, family = "Arial")
    p_cluster
    dev.off()    
    
    ## heatmap plot of three groups
    
  
    
    score_range <- dat_cell_touch_value %>% 
      bind_rows(dat_cell_pin_value, dat_cell_har_value) %>% 
      as_tibble() %>% 
      summarise(range = range(value, na.rm = TRUE)) %>%
      pull(range)
    
    group_day <- c("Touch", "Pin", "Har")
    
    for (i in c(1: 3)) {
      dat_trace <- dat_cell_touch_value %>% 
        bind_rows(dat_cell_pin_value, dat_cell_har_value) %>% 
        filter(Group == group_day[i])
      
      dat_rect <- dat_trace %>% 
        ddply(., .(Cluster), summarise, n = length(name)/60) %>% 
        add_column(xmin = 2, xmax = 2.1) %>% 
        mutate(ymax = cumsum(n) ) %>% 
        mutate(ymin = c(1, lag(ymax)[-1]))
      
      dat_cell_cluster_order <- dat_trace %>% 
        mutate(Cluster = factor(Cluster, levels= c("Excitation", "Inhibition"))) %>% 
        group_by(Cluster, name) %>%
        summarise(mean = mean(value), .groups = 'drop') %>%
        arrange(mean)
      


      p_heat <- dat_trace %>% 
        mutate(name = factor(name, levels = dat_cell_cluster_order$name)) %>% 
        mutate(Cluster = factor(Cluster, levels= c("Excitation", "Inhibition"))) %>% 
        ggplot(., aes(Time, name,fill= value))+ 
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
        theme(legend.position = "none")+
        annotate("rect", xmin = 2, xmax = 2.1, ymin = dat_rect$ymin, ymax = dat_rect$ymax-1, alpha = .2)
      
      assign(str_c("p_heat_", group_day[i]), p_heat)
    }
    
    ## combine the heat plot
    p_heat_com <- plot_grid(p_heat_Touch, p_heat_Pin, p_heat_Har, nrow = 1)
    
    setwd("~cchen/Documents/neuroscience/Pn\ project/Figure/PDF/")
    cairo_pdf("p_heat.pdf", width = 120/25.6, height = 90/25.6, family = "Arial")
    p_heat_com
    dev.off()
    
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
      
      dat_stim_average <- vector(mode = "list", length = length(t_stim))
      
      for (i in 1:length(trace_num)) {
        global_cell <- pull(global_ID[,trace_num[i]])
        
        col_name <- str_c("C", c(1: length(global_cell)))
        dat_stim_trace <- dat_trace[[2]][[trace_num[i]]] %>% 
          as_tibble() %>% 
          dplyr::select(all_of(global_cell)) %>% 
          mutate(across(.fns = scale)) %>% 
          setNames(col_name)
        

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
          ddply(.,.(Time, name), summarise, value1 = mean(value)) %>% 
          mutate(name = str_c(ID, name)) %>% 
          pivot_wider(names_from = name, values_from = value1, values_fill = list(value1 = NA)) %>% 
          select(-Time)
        
      }
      
      return(dat_stim_average)
    }
    
    path_trace_file <- as.list(list.files("~cchen/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/Pin_extract/stim_combine/", full.names = T))
    
    dat_trace_combine <- mapply(c_extract_align, path_trace_file, SIMPLIFY = F)
    rm(list=setdiff(ls(), "dat_trace_combine")) 
    
    dat_cell_trace_combine1 <- vector("list", 5)
    
    # Loop through each of the three sub-sublists
    for (i in 1:5) {
      # Extract the i-th sub-sublist from each main sublist and merge them
      dat_cell_trace_combine1[[i]] <- do.call(cbind, lapply(dat_trace_combine, function(x) x[[i]]))
      
    }
    
    
    stim <- c("Expect_touch", "Expect_pain", "Har", "Pin", "touch")
    
    dat_cell_trace_combine1 <- dat_cell_trace_combine1 %>% 
      do.call(rbind,.) %>% 
      mutate(Group = rep(stim, each = 60)) %>% 
      filter(!Group%in%c("Expect_touch", "Expect_pain")) 
    
    stim_select <- c( "Har", "Pin", "touch")
    
    dat_cell_trace_combine_cluster <- vector(mode = "list", length = length(stim_select))
    for (i in seq_along(stim_select)) {
      dat_cell_trace <- dat_cell_trace_combine1 %>% 
        filter(Group == stim_select[i]) %>% 
        select(-Group)
      
      dat_cell_trace_cluster <- kmeans(t(dat_cell_trace), centers = 2, iter.max = 10)[[1]]
      stim_time <- seq(-0.95, 2, by=0.05)
      
      dat_cell_value <- dat_cell_trace %>% 
        as_tibble() %>% 
        mutate(Time = stim_time) %>% 
        pivot_longer(-Time) %>% 
        mutate(Cluster = rep(dat_cell_trace_cluster, 60)) %>% 
        mutate(Group = stim_select[i])
      
      mean_values <- dat_cell_value %>%
        group_by(Cluster) %>%
        summarise(mean_value = mean(value, na.rm = TRUE)) %>%
        arrange(desc(mean_value))
      
      larger_mean_cluster <- mean_values$Cluster[1]
      
      # Assign "Excitation" or "Inhibition" to Clusters
      dat_cell_value <- dat_cell_value %>%
        mutate(Cluster = ifelse(Cluster == larger_mean_cluster, "Excitation", "Inhibition"))
      
      dat_cell_trace_combine_cluster[[i]] <- dat_cell_value
      
    }
    
    
    
    p_cluster <- dat_cell_trace_combine_cluster %>% 
      do.call(rbind,.) %>% 
      ddply(., .(Time, Cluster, Group), summarise,n=length(value),mean=mean(value),sd=sd(value),se=sd(value)/sqrt(length(value))) %>%
      mutate(Cluster = as.factor(Cluster)) %>% 
      mutate(Group = factor(Group, levels = c("touch", "Pin", "Har"))) %>% 
      ggplot(., aes(Time, mean, group = Group,color= Group))+
      geom_line()+
      geom_ribbon(aes(ymin=mean-se, ymax=mean+se, fill=Group), alpha=0.1, linetype=0)+
      #scale_color_manual(values=c("#8491B4FF", "#00A087FF", "#3C5488FF"))+
      facet_wrap(~Cluster, nrow = 1)+
      theme_bw()+
      theme(legend.title = element_blank())
    
    t_cluster_excitation <- dat_cell_trace_combine_cluster %>% 
      do.call(rbind,.) %>% 
      ddply(., .(Cluster, Group, name), summarise,value=mean(value)) %>%
      filter(Cluster == "Excitation") %>% 
      aov(value~Group, .)
    
    summary(t_cluster_excitation)
    TukeyHSD(t_cluster_excitation)
    
    p_cluster_inhibition <- dat_cell_trace_combine_cluster %>% 
      do.call(rbind,.)  %>% 
      ddply(., .(Cluster, Group, name), summarise,value=mean(value)) %>%
      mutate(Group = factor(Group, levels = c("touch", "Pin","Har"))) %>% 
      ggplot(., aes(Group, value, color = Group))+
      geom_violin()+
      geom_jitter(width = 0.2, alpha = 0.5)+
      facet_wrap(~Cluster, nrow = 1)+
      theme(axis.line.x = element_line(),
            axis.line.y = element_line(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.border = element_blank(),
            panel.background = element_blank(),
            axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
      scale_y_continuous(limits = c(-0.8,1))+
      theme(legend.title = element_blank(), legend.position = "none")
    
    setwd("~cchen/Documents/neuroscience/Pn\ project/Figure/PDF/")
    cairo_pdf("p_cluster_inhibition.pdf", width = 80/25.6, height = 60/25.6, family = "Arial")
    p_cluster_inhibition
    dev.off()  
    
    t_cluster_inhibition <- dat_cell_trace_combine_cluster %>% 
      do.call(rbind,.)  %>% 
      ddply(., .(Cluster, Group, name), summarise,value=mean(value)) %>%
      filter(Cluster == "Inhibition") %>% 
      aov(value~Group, .)
    
    summary(t_cluster_inhibition)
    TukeyHSD(t_cluster_inhibition)
    
 
      
    setwd("~cchen/Documents/neuroscience/Pn\ project/Figure/PDF/")
    cairo_pdf("p_trace_cluster.pdf", width = 80/25.6, height = 60/25.6, family = "Arial")
    p_cluster
    dev.off()    
    
    ## heatmap plot of three groups
    
    
    
    score_range <- dat_cell_trace_combine_cluster %>% 
      do.call(rbind,.)  %>% 
      as_tibble() %>% 
      summarise(range = range(value, na.rm = TRUE)) %>%
      pull(range)
    
    group_day <- c("touch", "Pin", "Har")
    
    for (i in c(1: 3)) {
      dat_trace <-  dat_cell_trace_combine_cluster %>% 
        do.call(rbind,.)  %>% 
        filter(Group == group_day[i])
      
      dat_rect <- dat_trace %>% 
        ddply(., .(Cluster), summarise, n = length(name)/60) %>% 
        add_column(xmin = 2, xmax = 2.1) %>% 
        mutate(ymax = cumsum(n) ) %>% 
        mutate(ymin = c(1, lag(ymax)[-1]))
      
      dat_cell_cluster_order <- dat_trace %>% 
        mutate(Cluster = factor(Cluster, levels= c("Excitation", "Inhibition"))) %>% 
        group_by(Cluster, name) %>%
        summarise(mean = mean(value), .groups = 'drop') %>%
        arrange(mean)
      
      
      
      p_heat <- dat_trace %>% 
        mutate(name = factor(name, levels = dat_cell_cluster_order$name)) %>% 
        mutate(Cluster = factor(Cluster, levels= c("Excitation", "Inhibition"))) %>% 
        ggplot(., aes(Time, name,fill= value))+ 
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
        theme(legend.position = "none")+
        annotate("rect", xmin = 2, xmax = 2.1, ymin = dat_rect$ymin, ymax = dat_rect$ymax-1, alpha = .2)
      
      assign(str_c("p_heat_", group_day[i]), p_heat)
    }
    
    ## combine the heat plot
    p_heat_com <- plot_grid(p_heat_touch, p_heat_Pin, p_heat_Har, nrow = 1)
    
    
    setwd("~cchen/Documents/neuroscience/Pn\ project/Figure/PDF/")
    cairo_pdf("p_heat.pdf", width = 120/25.6, height = 90/25.6, family = "Arial")
    p_heat_com
    dev.off()

    
## compare the locomotion and walk distance between mice with and without miniscope implantation-------
    cc_anti_behavior <- function(path_ID, frame_rate, length_pix) {
      dat_anti <- read.csv(path_ID, header = F, stringsAsFactors = F) %>% 
        .[-c(1:3), ] %>% 
        apply(., 2, as.numeric) %>%   ## only select the position of head
        as.data.frame() %>% 
        mutate(V2 = V2- V8, V3 = V3- V9, V5 = V5-V8, V6 = V6 -V9) %>%  ##normalize from FIJI
        select(str_c("V", 1:6)) 
      
      colnames(dat_anti)<- c("Frame", "Head_x","Head_y", "Prob", "Tail_x", "Tail_y")
      
      ## remove the frame for the first few frame
      dat_anti <- dat_anti[seq(1, nrow(dat_anti), frame_rate/10), ] %>% 
        subset(., Prob >0.5 | Frame> 60) %>% 
        slice(1:1600)
      
      ## change the dim from pix to mm (the dim of one chammer is 165mm)
      # value 387 for anti experiment 20200309
      # value 447 for miniscope 20191209
      r_pix_mm <- 165 / length_pix
      dat_anti <- dat_anti %>% 
        mutate(Head_x = abs(Head_x * r_pix_mm), Head_y = abs(Head_y * r_pix_mm)) %>% 
        mutate(Tail_x = abs(Tail_x * r_pix_mm), Tail_y = abs(Tail_y * r_pix_mm))
      
      chamber_div <- 165.5 ## the lenght of two chambels and the gap between
      
      ## assign each frame in hot as 1 (48) and nor plate(30) as 2
      dat_anti$chamber <- ifelse(dat_anti$Head_y < chamber_div, "Hot","Nor")
      ratio_time <- unname(by(dat_anti[,3], dat_anti$chamber, length)[2]/nrow(dat_anti))
      
      p_anti <- ggplot(dat_anti, aes(Head_y, Head_x, colour=chamber))+
        geom_point(shape=1)+
        theme_void()+
        theme(legend.position = "none")
      
      ## calculate how many times the mouse cross the border, and latency
      digit <- ifelse(dat_anti$Head_y< chamber_div, 1,2)
      
      for (i in c(2: (length(digit)-20))) {
        digit[i]<- ifelse(digit[i]!= digit[i-1] & digit[i+20] == digit[i], digit[i], digit[i-1] )
      }
      
      dat_anti$digit<- digit
      cross_digit <- diff(digit)
      cross_hot_nor <- which(cross_digit==1)+1
      
      cross_nor_hot <- which(cross_digit==-1)+1
      
      total_cross <- length(c(cross_hot_nor, cross_nor_hot))
      # calculate the latency of fist moving 
      # frame_rate <- 10 ##10 Hz
      cross_latency_1st  <- ifelse(digit[1]==2,2, cross_hot_nor[1]/10)
      cross_latency_1st <- ifelse(cross_latency_1st <2, 2, cross_latency_1st)
      
      # latency of cross back from nor to hot
      cross_latency_2nd <- ifelse(length(cross_nor_hot)==0, nrow(dat_anti)/10 - cross_latency_1st, 
                                  cross_nor_hot[1]/10 - cross_latency_1st)
      
      
      ## calculate the distance and speed of movement
      x_diff <- c(0, diff(dat_anti$Head_x))
      y_diff <- c(0, diff(dat_anti$Head_y))
      dis_frame <- sqrt(x_diff^2 + y_diff^2)
      move_speed <- dis_frame*10 
      
      move_speed_t <- zoo(move_speed) %>% 
        rollapply(., width = 10, by = 10, FUN = mean, na.rm = TRUE, align = "left") %>% 
        as.vector()
      peak_threshold <- mean(move_speed_t) +  sd(move_speed_t)
      peak_speed <- findpeaks(move_speed_t, nups = 2,minpeakdistance = 2, minpeakheight = peak_threshold) %>% 
        .[,2]*10 %>% 
        sort()
      latency_accerlation <- (peak_speed[which(peak_speed> cross_hot_nor[1])][1] - cross_hot_nor[1])/10
      # total distance, normalization
      total_distance <- round(sum(dis_frame), digits = 2)/1000  #unit to m
      
      # speed in different chambers
      speed_hot <- mean(move_speed[which(dat_anti$chamber == "Hot")]) 
      speed_nor <- mean(move_speed[which(dat_anti$chamber == "Nor")])
      speed_compare <- speed_nor/speed_hot
      # average speed (3s) before crossing the chambers
      speed_to_nor_1st <- ifelse(cross_hot_nor[1]< 20, mean(move_speed[0:cross_hot_nor[1]]) ,mean(move_speed[(cross_hot_nor[1]-20):cross_hot_nor[1]]))
      speed_to_hot_1st <- ifelse(length(cross_nor_hot)==0, NA, mean(move_speed[(cross_nor_hot[1]-20):cross_nor_hot[1]]))
      
      dat_anti_result <- tibble(ratio=ratio_time, latency1= cross_latency_1st, latency2 = cross_latency_2nd,
                                total_distance = total_distance, speed1 = speed_to_nor_1st, speed2= speed_to_hot_1st, total_cross= total_cross,
                                speed_compare=speed_hot, latency_accerlation= latency_accerlation)
      return(dat_anti_result)
    }
    
  ## for miniscope implanted mice

    path_trace_m3 <- list.files( "~cchen/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/matlab_analysis/behavior_tracking/m3", pattern = "*.csv", full.names = T )[c(1)]
    dat_behavior_m3 <- mapply(cc_anti_behavior, path_trace_m3, frame_rate = 20, length_pix=478,SIMPLIFY = F) %>% 
      do.call(rbind,.)
    
    # for m7

    path_trace_m7 <- list.files("~cchen/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/matlab_analysis/behavior_tracking/m7", pattern = "*.csv", full.names = T )[c(1)]
    dat_behavior_m7 <- mapply(cc_anti_behavior, path_trace_m7, frame_rate = 20, length_pix=478,SIMPLIFY = F) %>% 
      do.call(rbind,.)
    
    
    # for m17
    path_trace_m17 <- list.files("~cchen/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/matlab_analysis/behavior_tracking/m17", pattern = "*.csv", full.names = T )[c(1)]
    dat_behavior_m17 <- mapply(cc_anti_behavior, path_trace_m17, frame_rate = 20, length_pix=387,SIMPLIFY = F) %>% 
      do.call(rbind,.)
    
    # for m18

    path_trace_m18 <- list.files("~cchen/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/matlab_analysis/behavior_tracking/m18",  pattern = "*.csv", full.names = T )[c(1)]
    dat_behavior_m18 <- mapply(cc_anti_behavior, path_trace_m18, frame_rate = 20, length_pix=387,SIMPLIFY = F) %>% 
      do.call(rbind,.)
    
    # for m855

    path_trace_m855 <- list.files("~cchen/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/matlab_analysis/behavior_tracking/m855", pattern = "*.csv", full.names = T )[c(1)]
    dat_behavior_m855 <- mapply(cc_anti_behavior, path_trace_m855, frame_rate = 10, length_pix=364,SIMPLIFY = F) %>% 
      do.call(rbind,.)
    
    # for m857

    
    path_trace_m857 <- list.files("~cchen/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/matlab_analysis/behavior_tracking/m857", pattern = "*.csv", full.names = T )[c(1)]
    dat_behavior_m857 <- mapply(cc_anti_behavior, path_trace_m855, frame_rate = 10, length_pix=364,SIMPLIFY = F) %>% 
      do.call(rbind,.)

    
    dat_anti_mini <- rbind(dat_behavior_m3, dat_behavior_m7, dat_behavior_m17, dat_behavior_m18, dat_behavior_m855, dat_behavior_m857) %>% 
      mutate(Group = "Mini.") %>% 
      select(Group, speed_compare, total_distance)
    ## for mice without implantation
    dat_anti_ctrl <- vector(mode = "list", 7)
    dat_anti_con <- vector(mode = "list", 7)
    length_pix <- c(363,362,363,364,362,364,365)
    day_group <- c("D1","D2", "Pre", "Rm","Rm","Rm", "Test")
    for(i in c(1:7)){
      path_ctrl <- str_c("~cchen/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/wt_ctrl_07042020/wt_ctrl/T",i)
      path_con <-  str_c("~cchen/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/wt_ctrl_07042020/wt_con/T",i)
      dat_anti_ctrl[[i]]<- list.files(path_ctrl,pattern = ".csv", full.names = T ) %>% 
        as.list() %>% 
        mapply(cc_anti_behavior, ., frame_rate = 10, length_pix=length_pix[i], SIMPLIFY = F) %>% 
        do.call(rbind,.) %>% 
        add_column(Day = day_group[i], Group="Ctrl", ID=str_c("m",1:9), .before = 'ratio')
      
      dat_anti_con[[i]]<- list.files(path_con,pattern = ".csv", full.names = T ) %>% 
        as.list() %>% 
        mapply(cc_anti_behavior, ., frame_rate = 10, length_pix=length_pix[i], SIMPLIFY = F) %>% 
        do.call(rbind,.) %>% 
        add_column(Day = day_group[i], Group="Cond.", ID=str_c("m",1:10), .before = 'ratio')
    }
    
    dat_anti_compare <- do.call(rbind, dat_anti_ctrl) %>% 
      rbind(., do.call(rbind, dat_anti_con)) %>% 
      filter(Day == "D1") %>% 
      mutate(Group = "Ctrl") %>% 
      select(Group, speed_compare, total_distance) %>% 
      bind_rows(dat_anti_mini)  

    t_distance <- dat_anti_compare %>% 
      t.test(total_distance~Group, .)
    t_speed <- dat_anti_compare %>% 
      t.test(speed_compare~Group, .)
    
    p_distance <- dat_anti_compare %>%
      mutate(Group = factor(Group, levels = c("Ctrl", "Mini."))) %>% 
      ggplot(., aes(Group, total_distance, Group= Group))+
      geom_boxplot(outlier.shape = NA)+
      geom_jitter(aes(colour = Group, shape = Group),width = 0.2,  size=2)+
      scale_color_manual(values=c("darkcyan", "indianred"))+
      labs(x="", y="Total distance (m)")+
      theme(axis.line.x = element_line(),
            axis.line.y = element_line(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.border = element_blank(),
            panel.background = element_blank(),
            axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
      scale_y_continuous(limits = c(6, 16))+
      theme(legend.title = element_blank(), legend.position = "none")

        
    p_speed <- dat_anti_compare %>%
      mutate(Group = factor(Group, levels = c("Ctrl", "Mini."))) %>% 
      ggplot(., aes(Group, speed_compare, Group= Group))+
      geom_boxplot(outlier.shape = NA)+
      geom_jitter(aes(colour = Group, shape = Group),width = 0.2,  size=2)+
      scale_color_manual(values=c("darkcyan", "indianred"))+
      labs(x="", y="Moving speed (mm/s)")+
      theme(axis.line.x = element_line(),
            axis.line.y = element_line(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.border = element_blank(),
            panel.background = element_blank(),
            axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
      scale_y_continuous(limits = c(20, 120))+
      theme(legend.title = element_blank(), legend.position = "none")
    
    p_mini_compare <- plot_grid(p_distance, p_speed, nrow = 1)
    
    setwd("~cchen/Documents/neuroscience/Pn\ project/Figure/PDF/")
    cairo_pdf("p_mini_compare.pdf", width = 60/25.6, height = 60/25.6, family = "Arial")
    p_mini_compare
    dev.off()

    
    ## correlation analysis for the rACC-Pn neuron activity-----
    # run the script from 'Placebo_state_10012022.R'
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
      dat_cor_trace <- vector(mode = "list", length = length(num_compare))
      
      
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
        
        ## for correlation analysis 
        
        res_cor <- cor(dat_stim, method = "pearson")
        res_cor[res_cor==1]<- NA
        dat_cor_trace[[i]] <- res_cor
      }
      return(list(dat_stim_trace, dat_cor_trace))
    }
    
    mouse_file <- as.list(list.files("~cchen/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/FD-processed", full.names = T))
    
    dat_cell_trace <- mapply(c_miniscope_matlab_ft, mouse_file, SIMPLIFY = F)
    
    group_day <- c("Pre", "Cond", "Post")
    
    dat_cor_combine <- c()
    for (i in seq_along(group_day)){
      dat_cell_cor <- lapply(dat_cell_trace, function(x)  x[[2]]) %>% 
        lapply(., function(x)  x[[i]]) %>% 
        unlist() %>% 
        .[!is.na(.)] %>% 
        as_tibble() %>% 
        add_column(Group = group_day[i])
      
      dat_cor_combine <- bind_rows(dat_cor_combine, dat_cell_cor)
      
    }
    
    dat_cor_pac <- dat_cor_combine %>% 
      mutate(Group = factor(Group, levels = c("Pre", "Cond", "Post"))) %>% 
      filter(Group != 'Cond') %>% 
      ggplot(., aes(Group, value, cor = Group))+
      geom_violin()+
      geom_jitter(aes(colour = Group, shape = Group),width = 0.1,  size=2, alpha= 0.1)+
      scale_color_manual(values=c("#8491B4FF",  "#3C5488FF"))+
      labs(x="", y="C.C. of paired neurons")+
      theme(axis.line.x = element_line(),
            axis.line.y = element_line(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.border = element_blank(),
            panel.background = element_blank(),
            axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
      theme(legend.title = element_blank(), legend.position = "none")
    
    
    t_cor_pac <- dat_cor_combine %>% 
      aov(value~Group,.)
    
    summary(t_cor_pac)
    
    TukeyHSD(t_cor_pac)
    
    ## for frequency of ca2+ events---
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
      
      ## for the spike frequency file
      dat_spike <- list.files("~cchen/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/Extract/Events_0.5/", pattern = ID, full.names = T)
      
      dat_stim_trace <- vector(mode = "list", length = length(num_compare))
      dat_event_cor <- vector(mode = "list", length = length(num_compare))
      
      dat_event_day <- rep(0, 3)
      
      for (i in seq_along(num_compare)) {
        global_cell <- pull(cross_ID[,i])
        dat_trace <- dat_trace1[[num_compare[i]]] %>% 
          as_tibble() %>% 
          select(all_of(global_cell)) %>% 
          apply(., 2, scale)
        
        dat_event <- read.xlsx(dat_spike[i], colNames = F) %>% 
          as_tibble() %>% 
          select(all_of(global_cell))
        
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
        
        ## calculate the event rate
        dat_event1 <- dat_event[(t1_p-40):(t1_p+40-1),] %>% 
          apply(., 2, mean) %>% 
          mean()
        
        dat_event_day[i] <- dat_event1
        
        ## correlation analysis
        res_cor <- cor(dat_stim, method = "pearson")
        res_cor[res_cor==1]<- NA
        dat_event_cor[[i]] <- res_cor
        
      }
      return(dat_event_cor)
    }
    
    
    
    ## analyze for all mouse 
    mouse_file <- as.list(list.files("~cchen/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/FD-processed", full.names = T))
    
    dat_cell_trace <- mapply(c_miniscope_matlab_ft, mouse_file, SIMPLIFY = F)
    
    group_day <- c("Pre", "Cond", "Post")
    
    dat_event_cor_combine <- c()
    for (i in seq_along(group_day)){
      dat_cell_cor <- dat_cell_trace %>% 
        lapply(., function(x)  x[[i]]) %>% 
        unlist() %>% 
        .[!is.na(.)] %>% 
        as_tibble() %>% 
        add_column(Group = group_day[i])
      
      dat_event_cor_combine <- bind_rows(dat_event_cor_combine, dat_cell_cor)
      
    }
    
    dat_cor_pac1 <- dat_event_cor_combine %>% 
      mutate(Group = factor(Group, levels = c("Pre", "Cond", "Post"))) %>% 
      filter(Group != 'Cond') %>% 
      ggplot(., aes(Group, value, cor = Group))+
      geom_violin()+
      geom_jitter(aes(colour = Group, shape = Group),width = 0.1,  size=2, alpha= 0.2)+
      scale_color_manual(values=c("#8491B4FF", "#3C5488FF"))+
      labs(x="", y="C.C. of paired neurons")+
      theme(axis.line.x = element_line(),
            axis.line.y = element_line(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.border = element_blank(),
            panel.background = element_blank(),
            axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
      theme(legend.title = element_blank(), legend.position = "none")
    
    t_cor_pac <- dat_event_cor_combine %>% 
      aov(value~Group,.)
    
    summary(t_cor_pac)
    
    TukeyHSD(t_cor_pac)
    
    setwd("~cchen/Documents/neuroscience/Pn\ project/Figure/PDF/")
    cairo_pdf("p_dat_cor_pac.pdf", width = 50/25.6, height = 60/25.6, family = "Arial")
    dat_cor_pac1
    dev.off()
    