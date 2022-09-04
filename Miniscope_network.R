## build neuronal network based on their Ca2+ activity

## import librarys
library(igraph)


# compare the neuronal network on day7 of mouse m3, cross vs crossing back-----
c_miniscope_matlab_d7 <- function(ID_trace, t_stim) {
  ## import and format the data
  dat_trace <- read.xlsx(ID_trace, colNames = F, rowNames = F)
  ## analyze the trace
  dat_trace<- dat_trace[dat_trace[,1]==1, ]
  dat_trace<- as.data.frame(t(dat_trace[,-1]))
  rownames(dat_trace) <- NULL
  colnames(dat_trace)<- NULL
  ## do z score of the whole trace
  #dat_trace <- apply(dat_trace, 2, scale)
  
  ## divide ca2+ trace within -5 and 9s before and after stimuls, take -5 to -3 as baseline
  dat_stim <- vector(mode = "list", length = length(t_stim))
  stim_time<- seq(-2, 6.5, by=0.5)
  ## number of rows to be binned
  n <- 10 # 0.05*10=0.5
  
  ## extract cell activity when they cross the border
  for (i in seq_along(t_stim)){
    t1_p <- t_stim[i]
    dat_stim1 <- dat_trace[(t1_p-40):(t1_p+140-1),] 
    # dat_stim1 <- aggregate(dat_stim1,list(rep(1:(nrow(dat_stim1)%/%n+1),each=n,len=nrow(dat_stim1))),mean)[-1]
    # dat_stim1_base <- dat_stim1[1:5,] 
    # dat_stim1_base_mean <- colMeans(dat_stim1_base, na.rm = T)
    # dat_stim1_nor1 <- sweep(dat_stim1, 2, dat_stim1_base_mean, FUN = "-")
    # rownames(dat_stim1_nor1) <- NULL
    # colnames(dat_stim1_nor1)<- NULL
    # dat_stim[[i]] <- dat_stim1_nor1
    rownames(dat_stim1) <- NULL
    colnames(dat_stim1)<- NULL
    dat_stim[[i]] <- dat_stim1
  }
  
  return(dat_stim)
  
}
## for m3
path_trace_m3 <- "~cchen/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/matlab_analysis/m3/trace_days/m3_trace_d7.xlsx"
t_stim_m3_d7 <- c(346, 1206)

dat_stim_m3_d7 <- c_miniscope_matlab_d7(path_trace_m3, t_stim_m3_d7)

## for m7
path_trace_m7 <- "~cchen/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/matlab_analysis/m7/trace_days/m7_trace_d7.xlsx"
t_stim_m7_d7 <- c(236, 1914)

dat_stim_m7_d7 <- c_miniscope_matlab_d7(path_trace_m7, t_stim_m7_d7)

## for m17
path_trace_m17 <- "~cchen/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/matlab_analysis/m17/trace_days/m17_trace_d7.xlsx"
t_stim_m17_d7 <- c(157, 2113)

dat_stim_m17_d7 <- c_miniscope_matlab_d7(path_trace_m17, t_stim_m17_d7)

## for m18
path_trace_m18 <- "~cchen/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/matlab_analysis/m18/trace_days/m18_trace_d7.xlsx"
t_stim_m18_d7 <- c(124,2238)

dat_stim_m18_d7 <- c_miniscope_matlab_d7(path_trace_m18, t_stim_m18_d7)
## for m855
path_trace_m855 <- "~cchen/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/matlab_analysis/m855/trace_days/m855_trace_d7.xlsx"
t_stim_m855_d7 <- c(82,1572)

dat_stim_m855_d7 <- c_miniscope_matlab_d7(path_trace_m855, t_stim_m855_d7)

## for m857
path_trace_m857 <- "~cchen/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/matlab_analysis/m857/trace_days/m857_trace_d7.xlsx"
t_stim_m857_d7 <- c(334,1892)

dat_stim_m857_d7 <- c_miniscope_matlab_d7(path_trace_m857, t_stim_m857_d7)


## function of changing neuron activity to graph and and measure some properties------

cc_activity_graph <- function(dat_trace){
  
  # for cross
  colnames(dat_trace[[1]]) <- str_c("Cell", 1:ncol(dat_trace[[1]]))
  adj_cor_cross <- dat_trace[[1]] %>% 
    cor() 
  
  
  adj_cor_cross[upper.tri(adj_cor_cross)] <- 42
  adj_cor_cross1 <- adj_cor_cross %>% 
    melt() %>% 
    as_tibble() %>% 
    rename(From = Var1, To = Var2, Cor = value) %>% 
    filter(Cor != 42) %>% 
    filter(From != To) %>% 
    filter(Cor>0.3 | Cor < (-0.3)) %>% 
    mutate(Cor = rescale(Cor, to = c(0,1)))
  
  ## for crossback
  colnames(dat_trace[[2]]) <- str_c("Cell", 1:ncol(dat_trace[[2]]))
  adj_cor_crossback <- dat_trace[[2]] %>% 
    cor() 
  
  
  adj_cor_crossback[upper.tri(adj_cor_crossback)] <- 42
  adj_cor_crossback1 <- adj_cor_crossback %>% 
    melt() %>% 
    as_tibble() %>% 
    rename(From = Var1, To = Var2, Cor = value) %>% 
    filter(Cor != 42) %>% 
    filter(From != To) %>% 
    filter(Cor>0.3 | Cor < (-0.3)) %>% 
    mutate(Cor = rescale(Cor, to = c(0,1)))
  ## change the martix to network
  net_cross <- graph.data.frame(adj_cor_cross1, directed = FALSE)
  E(net_cross)$weight <- adj_cor_cross1$Cor
  
  net_crossback <- graph.data.frame(adj_cor_crossback1, directed = FALSE)
  E(net_crossback)$weight <- adj_cor_crossback1$Cor
  
  ##1. for graph density
  density_cross <- edge_density(net_cross)
  density_crossback <- edge_density(net_crossback)
  density <- c(density_cross, density_crossback)
  
  ##2. min distance
  dis_cross<- mean_distance(net_cross, directed = F, unconnected = F)
  dis_crossback<- mean_distance(net_crossback, directed = F, unconnected = F)
  
  dis <- c(dis_cross, dis_crossback)
  ##3. degree of each vertex
  d_net_cross <- tibble::enframe(igraph::degree(net_cross, mode = "total",  normalized = T)) %>% 
    rename(ID = name,d_cross = value) %>% 
    as_tibble()
  d_net_crossback <- tibble::enframe(igraph::degree(net_crossback, mode = "total", normalized = T)) %>% 
    rename(ID = name, d_crossback = value) %>% 
    as_tibble()
  
  dat_degree <- full_join(d_net_cross, d_net_crossback, by = "ID")
  
  ##4. clustering coefficiency of each vertex
  c_net_cross <- tibble::enframe(igraph::transitivity(net_cross,type = "local", isolates="NaN")) %>% 
    rename(ID = name,cluster_cross = value)
  c_net_crossback <- tibble::enframe(igraph::transitivity(net_crossback,type = "local", isolates="NaN")) %>% 
    rename(ID = name,cluster_crossback = value)
  
  dat_cluster <- full_join(c_net_cross, c_net_crossback, by = "ID")
  
  ##5. clustering analysis
  # cluster_cross <- cluster_walktrap(net_cross)
  # cluster_crossback <- cluster_walktrap(net_crossback)
  # 
  # modularity(cluster_cross)

  
  
  return(list(density, dis, dat_degree, dat_cluster, net_cross, net_crossback))
}
## change the ca2+ activity to adjecent martix based on correlation analysis
# for cross

mouse_ID <- c("m3", "m7", "m17", "m18", "m855", "m857")
dat_density <- vector(mode= "list", length = length(mouse_ID))
dat_dis <- vector(mode= "list", length = length(mouse_ID))
dat_degree <- vector(mode= "list", length = length(mouse_ID))
dat_cluster <- vector(mode= "list", length = length(mouse_ID))

for (i in seq_along(mouse_ID)) {
  dat_trace <- get(str_c("dat_stim_", mouse_ID[i], "_d7"))
  graph_prop <- cc_activity_graph(dat_trace)
  dat_density[[i]] <- graph_prop[[1]]
  dat_dis[[i]] <- graph_prop[[2]]
  dat_degree[[i]] <- graph_prop[[3]]
  dat_cluster[[i]] <- graph_prop[[4]]
}

## compare the density

p_density_com <- tibble(density = unlist(dat_density), ID = rep(mouse_ID, each=2), Group = rep(c("Cross", "Crossback"), length= length(mouse_ID)*2)) %>% 
  ggplot(., aes(Group, density, color = Group))+
  geom_boxplot()+
  geom_jitter()

t_density_com <- tibble(density = unlist(dat_density), ID = rep(mouse_ID, each=2), Group = rep(c("Cross", "Crossback"), length= length(mouse_ID)*2)) %>% 
  t.test(density~Group, ., alternative = c("greater"), paired = T)
  ## compare the dis
p_dis_com <- tibble(dis = unlist(dat_dis), ID = rep(mouse_ID, each=2), Group = rep(c("Cross", "Crossback"), length= length(mouse_ID)*2)) %>% 
  ggplot(., aes(Group, dis, color = Group))+
  geom_boxplot()+
  geom_jitter()

## compare the degree
p_degree_com <- do.call(rbind, dat_degree) %>% 
  pivot_longer(-ID) %>% 
  drop_na() %>% 
  ggplot(., aes(value, color=name))+
  stat_ecdf(geom = "step")+
  theme_classic()


## compare dat_cluster
dat_cluster_com <- do.call(rbind, dat_cluster) %>% 
  pivot_longer(-ID) %>% 
  drop_na() %>% 
  ggplot(., aes(value, color=name))+
  stat_ecdf(geom = "step")+
  theme_classic()












## plot the graph--------
colnames(dat_stim_m3_d7[[1]]) <- str_c("Cell", 1:ncol(dat_stim_m3_d7[[1]]))
adj_cor_cross<- dat_stim_m3_d7[[1]] %>% 
  cor() 

adj_cor_cross[upper.tri(adj_cor_cross)] <- 42
adj_cor_cross1 <- adj_cor_cross %>% 
  melt() %>% 
  as_tibble() %>% 
  rename(From = Var1, To = Var2, Cor = value) %>% 
  filter(Cor != 42) %>% 
  filter(From != To) %>% 
  filter(Cor>0.3 | Cor < (-0.3)) %>% 
  mutate(Cor = rescale(Cor, to = c(0,1)))

## for crossback
colnames(dat_stim_m3_d7[[2]]) <- str_c("Cell", 1:ncol(dat_stim_m3_d7[[2]]))
adj_cor_crossback <- dat_stim_m3_d7[[2]] %>% 
  cor() 


adj_cor_crossback[upper.tri(adj_cor_crossback)] <- 42
adj_cor_crossback1 <- adj_cor_crossback %>% 
  melt() %>% 
  as_tibble() %>% 
  rename(From = Var1, To = Var2, Cor = value) %>% 
  filter(Cor != 42) %>% 
  filter(From != To) %>% 
  filter(Cor>0.3 | Cor < (-0.3)) %>% 
  mutate(Cor = rescale(Cor, to = c(0,1)))
## change the martix to network

net_cross <- graph.data.frame(adj_cor_cross1, directed = FALSE)
E(net_cross)$weight <- adj_cor_cross1$Cor

coord_cross <- layout_nicely(net_cross)
coord_cross <- layout_with_fr(net_cross)



par(mar=c(0,0,0,0)+0.2)
plot(net_cross, layout= coord_cross,vertex.size=5, edge.width = E(net_cross)$weight)

par(mar=c(5,4,4,2)+0.1)

## for crossback network
net_crossback <- graph.data.frame(adj_cor_crossback1, directed = FALSE)
E(net_crossback)$weight <- adj_cor_crossback1$Cor

coord_crossback <- layout_nicely(net_crossback)
coord_crossback <- layout_with_fr(net_crossback)

par(mar=c(0,0,0,0)+0.2)
plot(net_crossback, layout= coord_crossback,vertex.size=5, edge.width = E(net_crossback)$weight)

par(mar=c(5,4,4,2)+0.1)

## compare the density
edge_density(net_cross)
edge_density(net_crossback)

## compare the min distance
mean_distance(net_cross, directed = F, unconnected = F)
mean_distance(net_crossback, directed = F, unconnected = F)

## histogram plot of the degree

d_net_cross <- tibble::enframe(igraph::degree(net_cross, mode = "total",  normalized = T)) %>% 
  rename(ID = name,d_cross = value)
d_net_crossback <- tibble::enframe(igraph::degree(net_crossback, mode = "total", normalized = T)) %>% 
  rename(d_crossback = value)

dat_d_cross <- cbind(d_net_cross, d_net_crossback[,2]) %>% 
  pivot_longer(-ID) %>% 
  ggplot(., aes(value, color=name))+
  stat_ecdf(geom = "step")+
  theme_classic()

## histogram plot of the clustering coefficient, compare with other three group

c_net_cross <- tibble::enframe(igraph::transitivity(net_cross,type = "local", isolates="NaN")) %>% 
  rename(ID = name,cluster_cross = value)
c_net_crossback <- tibble::enframe(igraph::transitivity(net_crossback,type = "local", isolates="NaN")) %>% 
  rename(ID = name,cluster_crossback = value)

##cummulative plot

dat_net_cluster <- cbind(c_net_cross, c_net_crossback[2]) %>% 
  pivot_longer(-ID) %>% 
  filter(value!='NaN' & value!=0) %>% 
  ggplot(., aes(value, color=name))+
  stat_ecdf(geom = "step")+
  theme_classic()


## compare pre and post conditiong using graph theory------
c_miniscope_matlab <- function(ID_trace, t_stim) {
  ## import and format the data
  dat_trace <- read.xlsx(ID_trace, colNames = F, rowNames = F) %>% 
    filter(., X1 ==1) %>% 
    select(., -c(1,2)) %>% 
    t() %>% 
    as.data.frame()
  
  rownames(dat_trace) <- NULL
  colnames(dat_trace)<- NULL
  ## do z score of the whole trace
  #dat_trace <- apply(dat_trace, 2, scale)
  
  
  ## extract cell activity when they cross the border-3, 15
  dat_stim <- dat_trace[(t_stim-60):(t_stim+200-1),] 

  return(dat_stim)
  
}
## Extract trace from each mice------
# for m3
path_trace_m3 <- as.list(list.files(path ="~cchen/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/matlab_analysis/m3/trace_days/", pattern = "*.xlsx", full.names = T ))[c(3,7)]
t_stim_m3 <- list(t_stim_m3_d3 = c(1885), t_stim_m3_d7 = c(346))
dat_trace_m3 <- mapply(c_miniscope_matlab, path_trace_m3, t_stim_m3, SIMPLIFY = F)

# for m7
path_trace_m7 <- as.list(list.files(path ="~cchen/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/matlab_analysis/m7/trace_days/", pattern = "*.xlsx", full.names = T ))[c(3,7)]
t_stim_m7 <- list(t_stim_m7_d3 = c(360), t_stim_m7_d7 = c(236))
dat_trace_m7 <- mapply(c_miniscope_matlab, path_trace_m7, t_stim_m7,SIMPLIFY = F)

# for m17
path_trace_m17 <- as.list(list.files(path ="~cchen/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/matlab_analysis/m17/trace_days/", pattern = "*.xlsx", full.names = T ))[c(3,7)]
t_stim_m17 <- list(t_stim_m17_d3 = c(437), t_stim_m17_d7 = c(157) )
dat_trace_m17 <- mapply(c_miniscope_matlab, path_trace_m17, t_stim_m17, SIMPLIFY = F)

# for m18
path_trace_m18 <- as.list(list.files(path ="~cchen/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/matlab_analysis/m18/trace_days/", pattern = "*.xlsx", full.names = T ))[c(3,7)]
t_stim_m18 <- list(t_stim_m18_d3 = c(784), t_stim_m18_d7 = c(124))
dat_trace_m18 <- mapply(c_miniscope_matlab, path_trace_m18, t_stim_m18, SIMPLIFY = F)

# for m855
path_trace_m855 <- as.list(list.files(path ="~cchen/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/matlab_analysis/m855/trace_days", pattern = "*.xlsx", full.names = T ))[c(3,7)]
t_stim_m855 <- list(t_stim_m855_d3 = c(392) *2, t_stim_m855_d7 = c(41)*2)
dat_trace_m855 <- mapply(c_miniscope_matlab, path_trace_m855, t_stim_m855, SIMPLIFY = F)

# for m857
path_trace_m857 <- as.list(list.files(path ="~cchen/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/matlab_analysis/m857/trace_days", pattern = "*.xlsx", full.names = T ))[c(1, 4)]
t_stim_m857 <- list(t_stim_m857_d3 = c(847) *2, t_stim_m857_d7 = c(167)*2)
dat_trace_m857 <- mapply(c_miniscope_matlab, path_trace_m857, t_stim_m857, SIMPLIFY = F)


cc_activity_graph <- function(dat_trace){
  
  # for cross
  colnames(dat_trace[[1]]) <- str_c("Cell", 1:ncol(dat_trace[[1]]))
  adj_cor_cross <- dat_trace[[1]] %>% 
    cor() 
  
  
  adj_cor_cross[upper.tri(adj_cor_cross)] <- 42
  adj_cor_cross1 <- adj_cor_cross %>% 
    melt() %>% 
    as_tibble() %>% 
    rename(From = Var1, To = Var2, Cor = value) %>% 
    filter(Cor != 42) %>% 
    filter(From != To) %>% 
    filter(Cor>0.3 | Cor < (-0.3)) %>% 
    mutate(Cor = rescale(Cor, to = c(0,1)))
  
  ## for crossback
  colnames(dat_trace[[2]]) <- str_c("Cell", 1:ncol(dat_trace[[2]]))
  adj_cor_crossback <- dat_trace[[2]] %>% 
    cor() 
  
  
  adj_cor_crossback[upper.tri(adj_cor_crossback)] <- 42
  adj_cor_crossback1 <- adj_cor_crossback %>% 
    melt() %>% 
    as_tibble() %>% 
    rename(From = Var1, To = Var2, Cor = value) %>% 
    filter(Cor != 42) %>% 
    filter(From != To) %>% 
    filter(Cor>0.3 | Cor < (-0.3)) %>% 
    mutate(Cor = rescale(Cor, to = c(0,1)))
  ## change the martix to network
  net_cross <- graph.data.frame(adj_cor_cross1, directed = FALSE)
  E(net_cross)$weight <- adj_cor_cross1$Cor
  
  net_crossback <- graph.data.frame(adj_cor_crossback1, directed = FALSE)
  E(net_crossback)$weight <- adj_cor_crossback1$Cor
  
  ##1. for graph density
  density_cross <- edge_density(net_cross)
  density_crossback <- edge_density(net_crossback)
  density <- c(density_cross, density_crossback)
  
  ##2. min distance
  dis_cross<- mean_distance(net_cross, directed = F, unconnected = F)
  dis_crossback<- mean_distance(net_crossback, directed = F, unconnected = F)
  
  dis <- c(dis_cross, dis_crossback)
  ##3. degree of each vertex
  d_net_cross <- tibble::enframe(igraph::degree(net_cross, mode = "total",  normalized = T)) %>% 
    rename(ID = name, pre = value) %>% 
    as_tibble()
  d_net_crossback <- tibble::enframe(igraph::degree(net_crossback, mode = "total", normalized = T)) %>% 
    rename(ID = name, post = value) %>% 
    as_tibble()
  
  dat_degree <- full_join(d_net_cross, d_net_crossback, by = "ID")
  
  ##4. clustering coefficiency of each vertex
  c_net_cross <- tibble::enframe(igraph::transitivity(net_cross,type = "local", isolates="NaN")) %>% 
    rename(ID = name, pre = value)
  c_net_crossback <- tibble::enframe(igraph::transitivity(net_crossback,type = "local", isolates="NaN")) %>% 
    rename(ID = name, post = value)
  
  dat_cluster <- full_join(c_net_cross, c_net_crossback, by = "ID")
  
  ##5. clustering analysis
  # cluster_cross <- cluster_walktrap(net_cross)
  # cluster_crossback <- cluster_walktrap(net_crossback)
  # 
  # modularity(cluster_cross)
  
  
  
  return(list(density, dis, dat_degree, dat_cluster))
}

# get graph properties of pre and post crossing

mouse_ID <- c("m3", "m7", "m17", "m18", "m855", "m857")
dat_density <- vector(mode= "list", length = length(mouse_ID))
dat_dis <- vector(mode= "list", length = length(mouse_ID))
dat_degree <- vector(mode= "list", length = length(mouse_ID))
dat_cluster <- vector(mode= "list", length = length(mouse_ID))
dat_net_pre <- vector(mode= "list", length = length(mouse_ID))
dat_net_post <- vector(mode= "list", length = length(mouse_ID))

for (i in seq_along(mouse_ID)) {
  dat_trace <- get(str_c("dat_trace_", mouse_ID[i]))
  graph_prop <- cc_activity_graph(dat_trace)
  dat_density[[i]] <- graph_prop[[1]]
  dat_dis[[i]] <- graph_prop[[2]]
  dat_degree[[i]] <- graph_prop[[3]]
  dat_cluster[[i]] <- graph_prop[[4]]
  dat_net_pre[[i]] <- graph_prop[[5]]
  dat_net_post[[i]] <- graph_prop[[6]]
  
}

## compare the density

p_density_com <- tibble(density = unlist(dat_density), ID = rep(mouse_ID, each=2), Group = rep(c("Pre", "Post"), length= length(mouse_ID)*2)) %>% 
  ggplot(., aes(Group, density, color = Group))+
  geom_boxplot()+
  geom_jitter()

t_density_com <- tibble(density = unlist(dat_density), ID = rep(mouse_ID, each=2), Group = rep(c("Pre", "Post"), length= length(mouse_ID)*2)) %>% 
  wilcox.test(density~Group, .)
## compare the dis
p_dis_com <- tibble(dis = unlist(dat_dis), ID = rep(mouse_ID, each=2), Group = rep(c("Pre", "Post"), length= length(mouse_ID)*2)) %>% 
  ggplot(., aes(Group, dis, color = Group))+
  geom_boxplot()+
  geom_jitter()

## compare the degree
p_degree_com <- do.call(rbind, dat_degree) %>% 
  pivot_longer(-ID) %>% 
  drop_na() %>% 
  mutate(name = factor(name, levels = c("d_cross",  "d_crossback"))) %>% 
  ggplot(., aes(value, color=name))+
  stat_ecdf(geom = "step")+
  scale_colour_manual(values=c("deepskyblue4", "indianred"))+
  labs(x="Degree", y="Probability")+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  theme(legend.title = element_blank(), legend.position = c(0.4, 0.8))


## compare dat_cluster
p_cluster_com <- do.call(rbind, dat_cluster) %>% 
  pivot_longer(-ID) %>% 
  drop_na() %>% 
  mutate(name = factor(name, levels = c("cluster_cross",  "cluster_crossback"))) %>% 
  ggplot(., aes(value, color=name))+
  stat_ecdf(geom = "step")+
  scale_colour_manual(values=c("deepskyblue4", "indianred"))+
  labs(x="Clustering coefficient", y="Probability")+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  theme(legend.title = element_blank(), legend.position = c(0.4, 0.8))

t_cluster <- do.call(rbind, dat_cluster) %>% 
  pivot_longer(-ID) %>% 
  drop_na() %>%
  as_tibble() %>% 
  wilcox.test(value ~ name, .)
summary(t_cluster)

p_graph_compare <- plot_grid(p_degree_com, p_cluster_com, nrow = 1)

setwd("~cchen/Documents/neuroscience/Pn\ project/Figure/PDF/")
cairo_pdf("p_graph_com.pdf", width = 120/25.6, height = 60/25.6, family = "Arial")
p_graph_compare
dev.off()

## plot the graph of pre and post conditioning
net_pre_m3 <- dat_net_pre[[1]]
coord_cross <- layout_nicely(net_pre_m3)

par(mar=c(0,0,0,0)+0.2)
plot(net_pre_m3, layout=layout.fruchterman.reingold,vertex.label=NA,vertex.size=4, vertex.color= "deepskyblue4", edge.width = E(net_pre_m3)$weight*0.5)

par(mar=c(5,4,4,2)+0.1)

## for m3 post
net_post_m3 <- dat_net_post[[1]]
coord_cross <- layout_nicely(net_post_m3)

par(mar=c(0,0,0,0)+0.2)
plot(net_post_m3, layout=layout.fruchterman.reingold, vertex.size=4, vertex.label=NA, vertex.color= "indianred", edge.width = E(net_post_m3)$weight*0.5)

par(mar=c(5,4,4,2)+0.1)
