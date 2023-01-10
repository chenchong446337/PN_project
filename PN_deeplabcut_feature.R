## for analyze deeplabcut behavior, tsne
library(smoother)
library(Rtsne)
library(feather)
library(factoextra)


cc_feature_fun <- function(file_path){
  dat_behavior <- read.csv(file_path, header = T, skip = 2) %>% 
    as_tibble() %>% 
    select(!starts_with("likelihood")) 
  # 1. for velocity
  m_velocity <- c(NA, sqrt((diff(dat_behavior$x))^2 + (diff(dat_behavior$y))^2)/0.1) %>% 
    smth.gaussian(., window = 10) # 1s smooth window
  
  ## 2. for acceleration
  m_accerleration <- c(NA, NA, diff(sqrt((diff(dat_behavior$x))^2 + (diff(dat_behavior$y))^2)/0.1)/0.1) %>% 
    smth.gaussian(., window = 10) # 1s smooth window
  
  
  
  ## 3. distance to reference
  dist_refernce <- sqrt((dat_behavior$x - mean(dat_behavior$x.5))^2 + (dat_behavior$y - mean(dat_behavior$y.5))^2) %>% 
    smth.gaussian(., window = 10) # 1s smooth window
  
  # direction to the reference
  direct_refernce <- atan2(dat_behavior$y - dat_behavior$y.5, dat_behavior$x - dat_behavior$x.5) %>% 
    smth.gaussian(., window = 10)
  
  ## 4. distance to the center
  # provide the side information
  point_right <- which(dat_behavior$x > mean(dat_behavior$x.6))
  dist_center <- sqrt((dat_behavior$x - mean(dat_behavior$x.6))^2 + (dat_behavior$y - mean(dat_behavior$y.6))^2)
  dist_center[point_right] = -dist_center[point_right]
  dist_center <-  smth.gaussian(dist_center, window = 10) # 1s smooth window
  
  ## 5. for the body direction
  body_direct <- atan2(dat_behavior$y, dat_behavior$x) %>% 
    smth.gaussian(., window = 10) # 1s smooth window
  
  ## 6. For direction to the door
  body_direction_door <- atan2(dat_behavior$y - mean(dat_behavior$y.6), dat_behavior$x - mean(dat_behavior$x.6)) %>% 
    smth.gaussian(., window = 10) # 1s smooth window
  
  ## 7. for body strench
  body_length <- sqrt((dat_behavior$x - dat_behavior$x.3)^2 + (dat_behavior$y - dat_behavior$y.3)^2) %>% 
    smth.gaussian(., window = 10) # 1s smooth window
  
  ## 8. head turn right
  head_turn_right <- sqrt((dat_behavior$x.1 - dat_behavior$x.3)^2 + (dat_behavior$y.1 - dat_behavior$y.3)^2) %>% 
    smth.gaussian(., window = 10) # 1s smooth window
  
  ## 9. head turn left
  head_turn_left <- sqrt((dat_behavior$x.2 - dat_behavior$x.3)^2 + (dat_behavior$y.2 - dat_behavior$y.3)^2)
  
  ## 10. tail direction
  tail_direct <- atan2(dat_behavior$y.4 - dat_behavior$y.3, dat_behavior$x.4 - dat_behavior$x.3) %>% 
    smth.gaussian(., window = 10)
  
  ## 11. distance of tail to head
  dist_tail_head <- sqrt((dat_behavior$x - dat_behavior$x.4)^2 + (dat_behavior$y - dat_behavior$y.4)^2) %>% 
    smth.gaussian(., window = 10)
  
  ## 12. change direction to the door
  body_direction_door_change <- c(NA, diff(atan2(dat_behavior$y - mean(dat_behavior$y.6), dat_behavior$x - mean(dat_behavior$x.6)))/0.1) %>% 
    smth.gaussian(., window = 10) # 1s smooth window
  
  ## make a tibble 
  
  # transforme data for t-sne analysis
  cc_unlist_fun <- function(x){
    dat<- unlist(x) %>%
      unname()
  }
  
  num_bin <- 10 # 1s
  num_trail <- floor(length(m_velocity)/num_bin)
  dat_feature <- tibble(Velocity = m_velocity, Accerleration = m_accerleration,dist_refernce, direct_refernce, dist_center,
                        body_direct, body_direction_door, body_length, head_turn_right, head_turn_left, tail_direct, dist_tail_head, body_direction_door_change) %>% 
    slice(1:(num_bin*num_trail)) %>% 
    mutate_each(funs(scale)) %>% 
    add_column(Trail = rep(1:num_trail, each = num_bin)) %>% 
    group_by(Trail) %>% 
    group_split(., .keep = F) %>% 
    mapply(cc_unlist_fun, ., SIMPLIFY = F) %>% 
    do.call(rbind,.) %>% 
    as_tibble() %>%
    drop_na() %>% 
    mutate(Trail = c(1: nrow(.)), .before = V1)

  dat_feature_sta <- tibble(Velocity = m_velocity, Accerleration = m_accerleration,dist_refernce, direct_refernce, dist_center,
                        body_direct, body_direction_door, body_length, head_turn_right, head_turn_left, tail_direct, dist_tail_head, body_direction_door_change) %>% 
    slice(1:(num_bin*num_trail)) %>% 
    mutate_each(funs(scale)) %>% 
    add_column(Trail = rep(1:num_trail, each = num_bin)) %>% 
    pivot_longer(-Trail) %>% 
    ddply(., .(Trail, name), summarise, mean_value = mean(value, na.rm = F )) %>% 
    pivot_wider(names_from = name, values_from = mean_value) %>% 
    drop_na() %>% 
    mutate(Trail = dat_feature$Trail)
  
  return(list(dat_feature, dat_feature_sta))
  
}

## for ctrl d5 and D6
file_path_ctrl <- list.files("~cchen/Documents/neuroscience/Pn\ project/Data_analysis/Deeplabcut_behavior/Ctrl/", full.names = T, recursive = T )

  

dat_feature_ctrl <- vector(mode = "list", length = length(file_path_ctrl))
dat_feature_ctrl_sta <- vector(mode = "list", length = length(file_path_ctrl))

for (i in seq_along(file_path_ctrl)) {
  test_day <- str_extract(file_path_ctrl[i], regex("D\\d+"))
  dat_feature_ctrl[[i]] <- cc_feature_fun(file_path_ctrl[[i]])[[1]] %>% 
    mutate(Trail = str_c("m", i,"_", Trail)) %>% 
    mutate(Day = test_day, .before = Trail) %>% 
    drop_na()
  
  dat_feature_ctrl_sta[[i]] <- cc_feature_fun(file_path_ctrl[[i]])[[2]] %>% 
    mutate(Trail = str_c("m", i,"_", Trail)) %>% 
    mutate(Day = test_day, .before = Trail)
  
}

dat_feature_ctrl <- dat_feature_ctrl %>% 
  do.call(rbind,.) %>% 
  as_tibble() %>% 
  add_column(Group = "Ctrl", .before = "Day")

dat_feature_ctrl_sta <- dat_feature_ctrl_sta %>% 
  do.call(rbind,.) %>% 
  as_tibble() %>% 
  add_column(Group = "Ctrl", .before = "Day")



file_path_cond <- list.files("~cchen/Documents/neuroscience/Pn\ project/Data_analysis/Deeplabcut_behavior/Cond/", full.names = T , recursive = T) 
  

dat_feature_cond <- vector(mode = "list", length = length(file_path_cond))
dat_feature_cond_sta <- vector(mode = "list", length = length(file_path_cond))

for (i in seq_along(file_path_cond)) {
  test_day <- str_extract(file_path_cond[i], regex("D\\d+"))
  
  dat_feature_cond[[i]] <- cc_feature_fun(file_path_cond[[i]])[[1]] %>% 
    mutate(Trail = str_c("m", i,"_", Trail)) %>% 
    mutate(Day = test_day, .before = Trail) %>% 
    drop_na()
  
  dat_feature_cond_sta[[i]] <- cc_feature_fun(file_path_cond[[i]])[[2]] %>% 
    mutate(Trail = str_c("m", i,"_", Trail)) %>% 
    mutate(Day = test_day, .before = Trail)
  
}

dat_feature_cond <- dat_feature_cond %>% 
  do.call(rbind,.) %>% 
  as_tibble() %>% 
  add_column(Group = "Cond", .before = "Day")

dat_feature_cond_sta <- dat_feature_cond_sta %>% 
  do.call(rbind,.) %>% 
  as_tibble() %>% 
  add_column(Group = "Cond", .before = "Day")

## combine data to do clustering
dat_behavior_com <- full_join(dat_feature_ctrl, dat_feature_cond) 
rm(dat_feature_cond)
rm(dat_feature_ctrl)

dat_behavior_com_sta <- full_join(dat_feature_ctrl_sta, dat_feature_cond_sta) 
rm(dat_feature_cond_sta)
rm(dat_feature_ctrl_sta)


dat_behavior_com_tsne <- dat_behavior_com %>% 
  select(-c(Group, Day, Trail)) %>% 
  as.matrix() %>% 
  Rtsne(., check_duplicates = FALSE, pca = FALSE, perplexity=50, theta=0.5, dims=2) %>% 
  .$Y

dat_feature_tsne <- dat_behavior_com_tsne %>% 
  as_tibble() %>% 
  bind_cols(dat_behavior_com[,1:3], .) %>% 
  bind_cols(., dat_behavior_com_sta[,-c(1:3)])

rm(dat_behavior_com_tsne)

#saveRDS(dat_feature_tsne, file = "~cchen/Documents/neuroscience/Pn\ project/Data_analysis/Deeplabcut_behavior/data/dat_feature_tsne.rds")
#saveRDS(dat_behavior_com, file = "~cchen/Documents/neuroscience/Pn\ project/Data_analysis/Deeplabcut_behavior/data/dat_behavior_com.rds")


### start from here-----
dat_feature_tsne <- readRDS("~cchen/Documents/neuroscience/Pn\ project/Data_analysis/Deeplabcut_behavior/data/dat_feature_tsne.rds")

# dat_behavior_com <- readRDS("~cchen/Documents/neuroscience/Pn\ project/Data_analysis/Deeplabcut_behavior/data/dat_behavior_com.rds")


## plot to see which parameter is important for the clustering
name_col <- colnames(dat_feature_tsne)[6:18]

for (i in name_col) {
  plx <- ggplot(dat_feature_tsne, aes(V1, V2, shape = Group))+
    geom_point(aes_string(colour = i))+
    scale_colour_gradientn(colours = terrain.colors(10))+
    theme_bw()
  print(plx)
  
}


plx <- dat_feature_tsne %>% 
  filter(Day != "D7") %>% 
  ggplot(., aes(V1, V2, color = Group))+
  geom_point()+
  theme_bw()

## plot the figure to show distance reference
p_dist_reference <- ggplot(dat_feature_tsne_opto, aes(V1, V2))+
  geom_point(aes(colour = dist_refernce))+
  scale_color_gradientn(colours = c("navy", "white", "red4"))+
  theme_void()+
  theme(legend.position = "none")

p_dist_center <- ggplot(dat_feature_tsne, aes(V1, V2))+
  geom_point(aes(colour = dist_center))+
  scale_color_gradientn(colours = c("navy", "white", "red4"))+
  theme_void()+
  theme(legend.position = "none")

dat_feature_d7 <- dat_feature_tsne %>% 
  filter(Day != "D7") %>% 
  mutate(Trail1 = as.numeric(str_extract(Trail, '(?<=_)\\d+'))) %>% 
  filter(Trail1 <=120) %>% 
  dplyr::select("V1", "V2", "Group") %>% 
  mutate(Group = factor(Group, levels = c("Ctrl", "Cond")))

p_feature_d7 <- dat_feature_d7 %>% 
  ggplot(., aes(V1, V2, group = Group))+
  stat_density_2d(aes(fill = after_stat(density)), geom = "raster", contour = FALSE)+
  scale_fill_viridis_c()+
  facet_grid(cols = vars(Group))+
  theme_void()+
  theme(legend.title = element_blank())

setwd("~cchen/Documents/neuroscience/Pn\ project/Figure/PDF/")
cairo_pdf("p_feature_d7.pdf", width = 140/25.6, height = 60/25.6, family = "Arial")
p_feature_d7
dev.off()

W <- owin( range(dat_feature_tsne$V1), range(dat_feature_tsne$V2) )


dat_feature_ppp <- ppp(x = dat_feature_d7$V1, y = dat_feature_d7$V2, window = W)

marks(dat_feature_ppp) <- as.factor(dat_feature_d7$Group)

is.multitype(dat_feature_ppp)
ripK <- Kcross(dat_feature_ppp)

p_ripk <- ripK %>% 
  as_tibble() %>% 
  pivot_longer(-r) %>% 
  filter(name =="theo"| name == "iso") %>% 
  mutate(name = factor(name, levels = c("theo", "iso"))) %>% 
  ggplot(., aes(r, value, color = name))+
  geom_line()+
  scale_color_manual(values = c("gray", "black"))+
  labs(x = "r", y = "K (r)")+
  theme_classic()+
  theme(legend.position = c(0.2, 0.8))

setwd("~cchen/Documents/neuroscience/Pn\ project/Figure/PDF/")
cairo_pdf("p_ripk.pdf", width = 80/25.6, height = 60/25.6, family = "Arial")
p_ripk
dev.off()

## calculate the density of points
rang_x_y <- c(-45, 48, -40, 42)

dat_feature_tsne_density <- MASS::kde2d(dat_feature_tsne$V1, dat_feature_tsne$V2, n= 200,lims = rang_x_y)

dat_feature_tsne_density$z <- ifelse(dat_feature_tsne_density$z < 0.00002, 0, dat_feature_tsne_density$z)
rk <- raster::raster(dat_feature_tsne_density) %>% 
  terra::rast()
values(rk ) <- values(rk)*10000
terra::plot(rk, legend = F, axes= F)

library(factoextra)

dat_feature_tsne %>% 
  dplyr::select(-c(1:5)) %>% 
  as.matrix() %>% 
  fviz_nbclust(., kmeans, method = "silhouette")+ labs(subtitle = "Silhouette method")

cluster_num <- 2
set.seed(42)
k_cluster <- dat_feature_tsne %>% 
  dplyr::select(-c(1:5)) %>% 
  as.matrix() %>% 
  kmeans(., centers = cluster_num, nstart = 50, iter.max = 500)


k_cluster_centroid <- k_cluster$centers %>% 
  as_tibble()

p_tsne_density <- dat_feature_tsne %>% 
  add_column(cluster = k_cluster$cluster) %>% 
  mutate(cluster = as.factor(cluster)) %>% 
  ggplot(., aes(x = V1, y = V2, col = cluster))+
  geom_point(alpha = 0.4)+
  geom_text(data = k_cluster_centroid, mapping = aes(x = V1, y = V2,label = 1:cluster_num),
            color = "black", size = 4) 
  

dat_feature_tsne_sta1 <- dat_feature_tsne %>% 
  add_column(cluster = k_cluster$cluster) %>% 
  mutate(cluster = as.factor(cluster)) %>% 
  mutate(ID = str_extract(Trail, regex("m\\d+"))) %>% 
  ddply(., .(Group, ID), summarise, n_total = length(cluster)) 

dat_feature_tsne_sta <- dat_feature_tsne %>% 
  add_column(cluster = k_cluster$cluster) %>% 
  mutate(cluster = as.factor(cluster)) %>% 
  mutate(ID = str_extract(Trail, regex("m\\d+"))) %>% 
  mutate(Trail1 = as.numeric(str_extract(Trail, '(?<=_)\\d+'))) %>% 
  filter(Day =="D7") %>% 
  filter(Trail1 < 60) %>% 
  ddply(., .(Group, cluster, ID), summarise, n = length(cluster)) %>% 
  right_join(., dat_feature_tsne_sta1, by = c("Group")) %>% 
  mutate(ratio = n/n_total) %>% 
  ddply(., .(Group, cluster), summarise, mean=mean(ratio), sd=sd(ratio),se=sd(ratio)/sqrt(length(ratio))) %>%  
  ggplot(., aes(x = cluster, y = mean, fill = Group, group = Group))+
  geom_bar(stat="identity", color="black", position=position_dodge()) +
  geom_errorbar(aes(ymin=mean, ymax=mean+sd), width=.2,
                position=position_dodge(.9))+
  scale_y_continuous(limits = c(0, 0.3))+
  coord_flip()
rm(dat_feature_tsne_sta1)

t_feature_tsne_sta <- dat_feature_tsne %>% 
  add_column(cluster = k_cluster$cluster) %>% 
  mutate(cluster = as.factor(cluster)) %>% 
  mutate(ID = str_extract(Trail, regex("m\\d+"))) %>% 
  filter(Day =="D7") %>% 
  ddply(., .(Group, ID,cluster), summarise, n = length(cluster)) %>% 
  right_join(., dat_feature_tsne_sta1, by = c("Group", "ID")) %>% 
  mutate(ratio = n/n_total) %>% 
  aov(ratio ~ Group + cluster,.)
  
summary(t_feature_tsne_sta)

p_tsne_density <- dat_feature_tsne %>% 
  add_column(cluster = k_cluster$cluster) %>% 
  mutate(cluster = as.factor(cluster)) %>% 
  ggplot(., aes(x = V1, y = V2, col = Group))+
  geom_point(alpha = 0.4)
 
## check the behavior of each cluster
dat_behavior <- dat_feature_tsne %>% 
  add_column(cluster = k_cluster$cluster) %>% 
  select(Group, Day, Trail, cluster) %>% 
  filter(Day =="D7")



## import feature extraction from python for further analysis---

dat_feature_ctrl_d56 <- read.csv("~cchen/Documents/neuroscience/Pn\ project/Data_analysis/Deeplabcut_behavior/data/dat_feature_ctrl_d56_feature.csv")

dat_feature_ctrl_d56_pca <- dat_feature_ctrl_d56 %>%
  as_tibble() %>% 
  select(-X) %>% 
  Filter(function(x) sd(x) != 0, .) %>% 
  as.matrix() %>% 
  prcomp(., scale = T, center = T) %>% 
  .$x %>% 
  .[, 1:500]

rm(dat_feature_ctrl_d56)
  
tsne_ctrl_d56 <- Rtsne(dat_feature_ctrl_d56_pca, check_duplicates = FALSE, pca = FALSE, perplexity=50, theta=0.5, dims=2)$Y


## for cond
dat_feature_cond_d56 <- read.csv("~cchen/Documents/neuroscience/Pn\ project/Data_analysis/Deeplabcut_behavior/data/dat_feature_cond_d56_feature.csv")

dat_feature_cond_d56_pca <- dat_feature_cond_d56 %>%
  as_tibble() %>% 
  select(-X) %>% 
  Filter(function(x) sd(x) != 0, .) %>% 
  as.matrix() %>% 
  prcomp(., scale = T, center = T) %>% 
  .$x %>% 
  .[, 1:500]

rm(dat_feature_cond_d56)

tsne_cond_d56 <- Rtsne(dat_feature_cond_d56_pca, check_duplicates = FALSE, pca = FALSE, perplexity=50, theta=0.5, dims=2)$Y

## combine data together

dat_behavior_d56 <- rbind(tsne_ctrl_d56, tsne_cond_d56) %>% 
  as_tibble() %>% 
  add_column(Group = c(rep("Ctrl", nrow(tsne_ctrl_d56)), rep("Cond", nrow(tsne_cond_d56))))

# tsne_cluster <- kmeans(tsne, 5, iter.max = 10, nstart = 1)$cluster

dat_tsne_d56 <- dat_behavior_d56 %>% 
  as_tibble() %>%
  ggplot(., aes(V1, V2, col = Group))+
  geom_point()


## for D7----
file_path_ctrl <- list.files("~cchen/Documents/neuroscience/Pn\ project/Data_analysis/Deeplabcut_behavior/Ctrl/D7", full.names = T )



dat_feature_ctrl_d7 <- vector(mode = "list", length = length(file_path_ctrl))

for (i in seq_along(file_path_ctrl)) {
  dat_feature_ctrl_d7[[i]] <- cc_feature_fun(file_path_ctrl[[i]]) %>% 
    mutate(Trail = str_c("m", i,"_", Trail)) %>% 
    drop_na()
  
}

dat_feature_ctrl_d7 <- dat_feature_ctrl_d7 %>% 
  do.call(rbind,.)


write.csv(dat_feature_ctrl_d7,"~cchen/Documents/neuroscience/Pn\ project/Data_analysis/Deeplabcut_behavior/data/dat_feature_ctrl_d7.csv", row.names = FALSE)


## for cond d5 and D6
file_path_cond <- list.files("~cchen/Documents/neuroscience/Pn\ project/Data_analysis/Deeplabcut_behavior/Cond/D7", full.names = T )


dat_feature_cond_d7 <- vector(mode = "list", length = length(file_path_cond))

for (i in seq_along(file_path_cond)) {
  dat_feature_cond_d7[[i]] <- cc_feature_fun(file_path_cond[[i]]) %>% 
    mutate(Trail = str_c("m", i,"_", Trail)) %>% 
    drop_na()
  
}

dat_feature_cond_d7 <- dat_feature_cond_d7 %>% 
  do.call(rbind,.)


write.csv(dat_feature_cond_d7,"~cchen/Documents/neuroscience/Pn\ project/Data_analysis/Deeplabcut_behavior/data/dat_feature_cond_d7.csv", row.names = FALSE)



## import feature extraction from python for further analysis---

dat_feature_ctrl_d7 <- read.csv("~cchen/Documents/neuroscience/Pn\ project/Data_analysis/Deeplabcut_behavior/data/dat_feature_ctrl_d7_feature.csv")

dat_feature_ctrl_d7_pca <- dat_feature_ctrl_d7 %>%
  as_tibble() %>% 
  select(-X) %>% 
  Filter(function(x) sd(x) != 0, .) %>% 
  as.matrix() %>% 
  prcomp(., scale = T, center = T) %>% 
  .$x %>% 
  .[, 1:500]

rm(dat_feature_ctrl_d7)

tsne_ctrl_d7 <- Rtsne(dat_feature_ctrl_d7_pca, check_duplicates = FALSE, pca = FALSE, perplexity=50, theta=0.5, dims=2)$Y


## for cond
dat_feature_cond_d7 <- read.csv("~cchen/Documents/neuroscience/Pn\ project/Data_analysis/Deeplabcut_behavior/data/dat_feature_cond_d7_feature.csv")

dat_feature_cond_d7_pca <- dat_feature_cond_d7 %>%
  as_tibble() %>% 
  select(-X) %>% 
  Filter(function(x) sd(x) != 0, .) %>% 
  as.matrix() %>% 
  prcomp(., scale = T, center = T) %>% 
  .$x %>% 
  .[, 1:500]

rm(dat_feature_cond_d7)

tsne_cond_d7 <- Rtsne(dat_feature_cond_d7_pca, check_duplicates = FALSE, pca = FALSE, perplexity=50, theta=0.5, dims=2)$Y

# get velocity information
ctrl_d7_velocity <- dat_feature_ctrl_d7 %>% 
  filter(name =="body_length") %>% 
  ddply(., .(Trail), summarise, mean = mean(value)) %>% 
  select(mean) %>% 
  unlist()


cond_d7_velocity <- dat_feature_cond_d7 %>% 
  filter(name =="body_length") %>% 
  ddply(., .(Trail), summarise, mean = mean(value)) %>% 
  select(mean) %>% 
  unlist()
## combine data together
cluster_optimal <- rbind(tsne_ctrl_d7, tsne_cond_d7) %>% 
  fviz_nbclust(., kmeans, method = "silhouette")+
  labs(subtitle = "Silhouette method")

tsne_cluster <- rbind(tsne_ctrl_d7, tsne_cond_d7) %>% 
  kmeans(., centers = 3) %>% 
  .$cluster

dat_behavior_d7 <- rbind(tsne_ctrl_d7, tsne_cond_d7) %>% 
  as_tibble() %>% 
  add_column(Group = c(rep("Ctrl", nrow(tsne_ctrl_d7)), rep("Cond", nrow(tsne_cond_d7)))) %>% 
  add_column(Velocity = c(ctrl_d7_velocity, cond_d7_velocity)) %>% 
  add_column(Cluster = tsne_cluster)

# tsne_cluster <- kmeans(tsne, 5, iter.max = 10, nstart = 1)$cluster

dat_tsne_d7 <- dat_behavior_d7 %>% 
  as_tibble() %>%
  ggplot(., aes(V1, V2))+
  stat_density_2d(aes(fill = ..density..), geom = "raster", contour = FALSE)+
  geom_point()
  scale_colour_gradientn(colours = terrain.colors(10))

  
  
  geom_hex(bins = 70) +
  scale_fill_continuous(type = "viridis") +
  theme_bw()

  
  
  geom_point(aes(shape = Group))+
  stat_density2d()
  scale_colour_gradientn(colours = terrain.colors(10))


dat_behavior_d7_sta <- rbind(tsne_ctrl_d7, tsne_cond_d7) %>% 
  as_tibble() %>% 
  add_column(Group = c(rep("Ctrl", nrow(tsne_ctrl_d7)), rep("Cond", nrow(tsne_cond_d7)))) %>% 
  add_column(Velocity = c(ctrl_d7_velocity, cond_d7_velocity)) %>% 
  add_column(Cluster = tsne_cluster) %>% 
  ddply(., .(Cluster, Group),summarise, n= length(V1))



## for opto control the activity of rACC-Pn----
file_path_opto_ctrl <- list.files("~cchen/Documents/neuroscience/Pn\ project/Data_analysis/Deeplabcut_behavior/Opto_anti/eYFP/", full.names = T, recursive = T )

dat_feature_opto_ctrl <- vector(mode = "list", length = length(file_path_opto_ctrl))
dat_feature_opto_ctrl_sta <- vector(mode = "list", length = length(file_path_opto_ctrl))

for (i in seq_along(file_path_opto_ctrl)) {
  test_day <- str_extract(file_path_opto_ctrl[i], regex("D\\d+"))
  dat_feature_opto_ctrl[[i]] <- cc_feature_fun(file_path_opto_ctrl[[i]])[[1]] %>% 
    mutate(Trail = str_c("m", i,"_", Trail)) %>% 
    mutate(Day = test_day, .before = Trail) %>% 
    drop_na()
  
  dat_feature_opto_ctrl_sta[[i]] <- cc_feature_fun(file_path_opto_ctrl[[i]])[[2]] %>% 
    mutate(Trail = str_c("m", i,"_", Trail)) %>% 
    mutate(Day = test_day, .before = Trail)
  
}

dat_feature_opto_ctrl <- dat_feature_opto_ctrl %>% 
  do.call(rbind,.) %>% 
  as_tibble() %>% 
  add_column(Group = "eYFP", .before = "Day")

dat_feature_opto_ctrl_sta <- dat_feature_opto_ctrl_sta %>% 
  do.call(rbind,.) %>% 
  as_tibble() %>% 
  add_column(Group = "eYFP", .before = "Day")



file_path_opto_chr <- list.files("~cchen/Documents/neuroscience/Pn\ project/Data_analysis/Deeplabcut_behavior/Opto_anti/Chr2/", full.names = T , recursive = T) 


dat_feature_opto_chr <- vector(mode = "list", length = length(file_path_opto_chr))
dat_feature_opto_chr_sta <- vector(mode = "list", length = length(file_path_opto_chr))

for (i in seq_along(file_path_opto_chr)) {
  test_day <- str_extract(file_path_opto_chr[i], regex("D\\d+"))
  
  dat_feature_opto_chr[[i]] <- cc_feature_fun(file_path_opto_chr[[i]])[[1]] %>% 
    mutate(Trail = str_c("m", i,"_", Trail)) %>% 
    mutate(Day = test_day, .before = Trail) %>% 
    drop_na()
  
  dat_feature_opto_chr_sta[[i]] <- cc_feature_fun(file_path_opto_chr[[i]])[[2]] %>% 
    mutate(Trail = str_c("m", i,"_", Trail)) %>% 
    mutate(Day = test_day, .before = Trail)
  
}

dat_feature_opto_chr <- dat_feature_opto_chr %>% 
  do.call(rbind,.) %>% 
  as_tibble() %>% 
  add_column(Group = "Chr2", .before = "Day")

dat_feature_opto_chr_sta <- dat_feature_opto_chr_sta %>% 
  do.call(rbind,.) %>% 
  as_tibble() %>% 
  add_column(Group = "Chr2", .before = "Day")


file_path_opto_enphr <- list.files("~cchen/Documents/neuroscience/Pn\ project/Data_analysis/Deeplabcut_behavior/Opto_anti/eNPhR/", full.names = T , recursive = T) 


dat_feature_opto_enphr <- vector(mode = "list", length = length(file_path_opto_enphr))
dat_feature_opto_enphr_sta <- vector(mode = "list", length = length(file_path_opto_enphr))

for (i in seq_along(file_path_opto_enphr)) {
  test_day <- str_extract(file_path_opto_enphr[i], regex("D\\d+"))
  
  dat_feature_opto_enphr[[i]] <- cc_feature_fun(file_path_opto_enphr[[i]])[[1]] %>% 
    mutate(Trail = str_c("m", i,"_", Trail)) %>% 
    mutate(Day = test_day, .before = Trail) %>% 
    drop_na()
  
  dat_feature_opto_enphr_sta[[i]] <- cc_feature_fun(file_path_opto_enphr[[i]])[[2]] %>% 
    mutate(Trail = str_c("m", i,"_", Trail)) %>% 
    mutate(Day = test_day, .before = Trail)
  
}

dat_feature_opto_enphr <- dat_feature_opto_enphr %>% 
  do.call(rbind,.) %>% 
  as_tibble() %>% 
  add_column(Group = "eNpHR", .before = "Day")

dat_feature_opto_enphr_sta <- dat_feature_opto_enphr_sta %>% 
  do.call(rbind,.) %>% 
  as_tibble() %>% 
  add_column(Group = "eNpHR", .before = "Day")


## combine data to do clustering
dat_behavior_com_opto <- full_join(dat_feature_opto_ctrl, dat_feature_opto_chr) %>% 
  full_join(., dat_feature_opto_enphr)
  
rm(dat_feature_opto_ctrl)
rm(dat_feature_opto_chr)
rm(dat_feature_opto_enphr)

dat_behavior_com_sta <- full_join(dat_feature_opto_ctrl_sta, dat_feature_opto_chr_sta) %>% 
  full_join(., dat_feature_opto_enphr_sta)
rm(dat_feature_opto_ctrl_sta)
rm(dat_feature_opto_chr_sta)
rm(dat_feature_opto_enphr_sta)


dat_behavior_opto_com_tsne <- dat_behavior_com_opto %>% 
  select(-c(Group, Day, Trail)) %>% 
  as.matrix() %>% 
  Rtsne(., check_duplicates = FALSE, pca = FALSE, perplexity=50, theta=0.5, dims=2) %>% 
  .$Y

dat_feature_tsne_opto <- dat_behavior_opto_com_tsne %>% 
  as_tibble() %>% 
  bind_cols(dat_behavior_com_opto[,1:3], .) %>% 
  bind_cols(., dat_behavior_com_sta[,-c(1:3)])

rm(dat_behavior_com_tsne)


#saveRDS(dat_feature_tsne_opto, file = "~cchen/Documents/neuroscience/Pn\ project/Data_analysis/Deeplabcut_behavior/data/dat_feature_tsne_opto.rds")
dat_feature_tsne_opto <- readRDS("~cchen/Documents/neuroscience/Pn\ project/Data_analysis/Deeplabcut_behavior/data/dat_feature_tsne_opto.rds")

dat_feature_d7_opto <- dat_feature_tsne_opto %>% 
  filter(Day != "D7") %>% 
  mutate(Trail1 = as.numeric(str_extract(Trail, '(?<=_)\\d+'))) %>% 
  filter(Trail1 <=120) %>% 
  dplyr::select("V1", "V2", "Group") %>% 
  mutate(Group = factor(Group, levels = c("eYFP", "Chr2", "eNpHR")))

p_feature_d7_opto <- dat_feature_d7_opto %>% 
  ggplot(., aes(V1, V2, group = Group))+
  stat_density_2d(aes(fill = after_stat(density)), geom = "raster", contour = FALSE)+
  scale_fill_viridis_c()+
  facet_grid(cols = vars(Group))+
  theme_void()+
  theme(legend.title = element_blank())


W <- owin( range(dat_feature_tsne_opto$V1), range(dat_feature_tsne_opto$V2) )

dat_feature_d7_opto_ctrl_chr <- dat_feature_d7_opto %>% 
  filter(Group != "eNpHR")
dat_feature_ppp <- ppp(x = dat_feature_d7_opto_ctrl_chr$V1, y = dat_feature_d7_opto_ctrl_chr$V2, window = W)

marks(dat_feature_ppp) <- as.factor(dat_feature_d7_opto_ctrl_chr$Group)

is.multitype(dat_feature_ppp)
ripK <- Kcross(dat_feature_ppp)
plot(ripK)
p_ripk <- ripK %>% 
  as_tibble() %>% 
  pivot_longer(-r) %>% 
  filter(name =="theo"| name == "iso") %>% 
  mutate(name = factor(name, levels = c("theo", "iso"))) %>% 
  ggplot(., aes(r, value, color = name))+
  geom_line()+
  scale_color_manual(values = c("gray", "black"))+
  labs(x = "r", y = "K (r)")+
  theme_classic()+
  theme(legend.position = c(0.2, 0.8))

