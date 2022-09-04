## for the nvue system
dat_cell_red <- read.xlsx("~cchen2/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/nvue/nvue_cells/cc_nvue_red_cell.xlsx") %>%
  as_tibble() %>% 
  mutate_at(vars(-("Time")),scale) %>% 
  mutate(mean = rowMeans(select(.,starts_with("C")))) %>% 
  select(Time, mean)

dat_cell_green <- read.xlsx("~cchen2/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/nvue/nvue_cells/cc_nvue_green_cell.xlsx") %>%
  as_tibble() %>% 
  mutate_at(vars(-("Time")),scale) %>% 
  mutate(mean = rowMeans(select(.,starts_with("C")))) %>% 
  select(Time, mean)

dat_cell_nvue <- rbind(dat_cell_green, dat_cell_red) %>% 
  mutate(Group = c(rep("IT", nrow(dat_cell_green)), rep("PT", nrow(dat_cell_red)))) %>% 
  ggplot(., aes(Time, mean, colour = Group))+
  geom_line()+
  labs(x = "Time (s)", y = "z-Score")+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  theme(legend.title =  element_blank())
  
dat_cell_green_cells <- read.xlsx("~cchen2/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/nvue/nvue_cells/cc_nvue_green_cell.xlsx") %>%
  as_tibble() %>% 
  mutate_at(vars(-("Time")),scale) %>% 
  mutate(mean = rowMeans(select(.,starts_with("C")))) %>% 
  select(Time, C00, C01, C02, C03) %>% 
  pivot_longer(-Time, names_to = "Cell") %>% 
  mutate(Group = "IT")

dat_cell_red_cells <- read.xlsx("~cchen2/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/nvue/nvue_cells/cc_nvue_red_cell.xlsx") %>%
  as_tibble() %>% 
  mutate_at(vars(-("Time")),scale) %>% 
  mutate(mean = rowMeans(select(.,starts_with("C")))) %>% 
  select(Time, C00, C01, C02, C03) %>% 
  pivot_longer(-Time, names_to = "Cell") %>% 
  mutate(Group = "PT")

plot_cells_compare <- rbind(dat_cell_green_cells, dat_cell_red_cells) %>% 
  ggplot(., aes(Time, value))+
  facet_grid(Cell ~ Group)+
  geom_line(aes(colour = Group))+
  theme_void()

## for m18
dat_cell_red <- read.csv("~cchen2/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/nvue/IT_PT/m18/m18_d1_RedTrace.csv", skip = 1) %>% 
  as_tibble() %>% 
  select(!starts_with("r")) %>% 
  rename(Time = names(.)[1]) %>% 
  mutate_at(vars(-("Time")),scale) %>% 
  mutate(mean = rowMeans(select(.,starts_with("a")))) %>% 
  select(Time, mean)



dat_cell_green <- read.csv("~cchen2/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/nvue/IT_PT/m18/m18_d1_GreenTrace.csv", skip = 1) %>% 
  as_tibble() %>% 
  select(!starts_with("r")) %>% 
  rename(Time = names(.)[1]) %>% 
  mutate_at(vars(-("Time")),scale) %>% 
  mutate(mean = rowMeans(select(.,starts_with("a")))) %>% 
  select(Time, mean)

## Plot the trace to show
dat_cell_green_trace <- read.csv("~cchen2/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/nvue/IT_PT/m18/m18_d1_GreenTrace.csv", skip = 1) %>% 
  as_tibble() %>% 
  select(c(1,5, 7, 15, 35,36)) %>% 
  rename(Time = names(.)[1]) %>% 
  pivot_longer(-(Time), names_to = "Cell") %>% 
  ggplot(., aes(Time, value, color = Cell))+
  facet_grid(rows = vars(Cell))+
  geom_line()+
  theme_void()+
  theme(legend.position = 'none')
  

cor(dat_cell_red$mean, dat_cell_green$mean)


dat_cell_nvue_d1 <- rbind(dat_cell_green, dat_cell_red) %>% 
  mutate(Group = c(rep("IT", nrow(dat_cell_green)), rep("PT", nrow(dat_cell_red)))) %>% 
  mutate(Group = factor(Group, levels = c("IT", "PT"))) %>% 
  ggplot(., aes(Time, mean, colour = Group))+
  geom_line()+
  scale_color_manual(values=c("darkcyan", "indianred"))+
  labs(x = "Time (s)", y = "z-Score", title = "D1")+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  theme(legend.title =  element_blank())

print(ccf(dat_cell_green$mean, dat_cell_red$mean))

## for day 2
dat_cell_red <- read.csv("~cchen2/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/nvue/IT_PT/m18/m18_d2_redtrace.csv", skip = 1) %>% 
  as_tibble() %>% 
  select(!starts_with("r")) %>% 
  rename(Time = names(.)[1]) %>% 
  mutate_at(vars(-("Time")),scale) %>% 
  mutate(mean = rowMeans(select(.,starts_with("a")))) %>% 
  select(Time, mean)



dat_cell_green <- read.csv("~cchen2/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/nvue/IT_PT/m18/m18_d2_Greentrace.csv", skip = 1) %>% 
  as_tibble() %>% 
  select(!starts_with("r")) %>% 
  rename(Time = names(.)[1]) %>% 
  mutate_at(vars(-("Time")),scale) %>% 
  mutate(mean = rowMeans(select(.,starts_with("a")))) %>% 
  select(Time, mean)

cor(dat_cell_red$mean, dat_cell_green$mean)


dat_cell_nvue_d2 <- rbind(dat_cell_green, dat_cell_red) %>% 
  mutate(Group = c(rep("IT", nrow(dat_cell_green)), rep("PT", nrow(dat_cell_red)))) %>% 
  ggplot(., aes(Time, mean, colour = Group))+
  geom_line()+
  labs(x = "Time (s)", y = "z-Score", title = "D2")+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  theme(legend.title =  element_blank())

print(ccf(dat_cell_green$mean, dat_cell_red$mean))

## for day 3
dat_cell_red <- read.csv("~cchen2/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/nvue/IT_PT/m18/m18_d3_redtrace.csv", skip = 1) %>% 
  as_tibble() %>% 
  select(!starts_with("r")) %>% 
  rename(Time = names(.)[1]) %>% 
  mutate_at(vars(-("Time")),scale) %>% 
  mutate(mean = rowMeans(select(.,starts_with("a")))) %>% 
  select(Time, mean)



dat_cell_green <- read.csv("~cchen2/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/nvue/IT_PT/m18/m18_d3_Greentrace.csv", skip = 1) %>% 
  as_tibble() %>% 
  select(!starts_with("r")) %>% 
  rename(Time = names(.)[1]) %>% 
  mutate_at(vars(-("Time")),scale) %>% 
  mutate(mean = rowMeans(select(.,starts_with("a")))) %>% 
  select(Time, mean)

cor(dat_cell_red$mean, dat_cell_green$mean)


dat_cell_nvue_d3 <- rbind(dat_cell_green, dat_cell_red) %>% 
  mutate(Group = c(rep("IT", nrow(dat_cell_green)), rep("PT", nrow(dat_cell_red)))) %>% 
  ggplot(., aes(Time, mean, colour = Group))+
  geom_line()+
  labs(x = "Time (s)", y = "z-Score", title = "D3")+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  theme(legend.title =  element_blank())

r_cross <- ccf(dat_cell_red$mean,  dat_cell_green$mean, plot = F)

p_auto_correlation <- tibble(lag = r_cross$lag, r = r_cross$acf) %>% 
  ggplot(., aes(lag, r))+
  geom_line()

## for day 4
dat_cell_red <- read.csv("~cchen2/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/nvue/IT_PT/m18/m18_d4_redtrace.csv", skip = 1) %>% 
  as_tibble() %>% 
  select(!starts_with("r")) %>% 
  rename(Time = names(.)[1]) %>% 
  mutate_at(vars(-("Time")),scale) %>% 
  mutate(mean = rowMeans(select(.,starts_with("a")))) %>% 
  select(Time, mean)



dat_cell_green <- read.csv("~cchen2/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/nvue/IT_PT/m18/m18_d4_Greentrace.csv", skip = 1) %>% 
  as_tibble() %>% 
  select(!starts_with("r")) %>% 
  rename(Time = names(.)[1]) %>% 
  mutate_at(vars(-("Time")),scale) %>% 
  mutate(mean = rowMeans(select(.,starts_with("a")))) %>% 
  select(Time, mean)

cor(dat_cell_red$mean, dat_cell_green$mean)


dat_cell_nvue_d4 <- rbind(dat_cell_green, dat_cell_red) %>% 
  mutate(Group = c(rep("IT", nrow(dat_cell_green)), rep("PT", nrow(dat_cell_red)))) %>% 
  ggplot(., aes(Time, mean, colour = Group))+
  geom_line()+
  labs(x = "Time (s)", y = "z-Score", title = "D4")+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  theme(legend.title =  element_blank())

r_cross <- ccf(dat_cell_red$mean,  dat_cell_green$mean, plot = F)

p_auto_correlation <- tibble(lag = r_cross$lag, r = r_cross$acf) %>% 
  ggplot(., aes(lag, r))+
  geom_line()

## for day 5
dat_cell_red <- read.csv("~cchen2/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/nvue/IT_PT/m18/m18_d5_redtrace.csv", skip = 1) %>% 
  as_tibble() %>% 
  select(!starts_with("r")) %>% 
  rename(Time = names(.)[1]) %>% 
  mutate_at(vars(-("Time")),scale) %>% 
  mutate(mean = rowMeans(select(.,starts_with("a")))) %>% 
  select(Time, mean)



dat_cell_green <- read.csv("~cchen2/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/nvue/IT_PT/m18/m18_d5_Greentrace.csv", skip = 1) %>% 
  as_tibble() %>% 
  select(!starts_with("r")) %>% 
  rename(Time = names(.)[1]) %>% 
  mutate_at(vars(-("Time")),scale) %>% 
  mutate(mean = rowMeans(select(.,starts_with("a")))) %>% 
  select(Time, mean)

cor(dat_cell_red$mean, dat_cell_green$mean)


dat_cell_nvue_d5 <- rbind(dat_cell_green, dat_cell_red) %>% 
  mutate(Group = c(rep("IT", nrow(dat_cell_green)), rep("PT", nrow(dat_cell_red)))) %>% 
  ggplot(., aes(Time, mean, colour = Group))+
  geom_line()+
  labs(x = "Time (s)", y = "z-Score", title = "D5")+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  theme(legend.title =  element_blank())

r_cross <- ccf(dat_cell_red$mean,  dat_cell_green$mean, plot = F)

p_auto_correlation <- tibble(lag = r_cross$lag, r = r_cross$acf) %>% 
  ggplot(., aes(lag, r))+
  geom_line()

## correlation between velocity and calcium transit--
dat_velocity <- read.xlsx("~cchen2/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/nvue/IT_PT/m18_d5_velocity.xlsx") %>% 
  as_tibble() %>% 
  select(Time, Velocity) %>% 
  slice(which(row_number() %% 2 == 1)) %>% 
  mutate_at(vars(-("Time")),scale) %>%
  mutate(Green = dat_cell_green$mean[1:1322], red = dat_cell_red$mean[1:1322])
  
cor.test(dat_velocity$Velocity, dat_velocity$Green)
cor.test(dat_velocity$Velocity, dat_velocity$red)
print(ccf(dat_velocity$Velocity[1:1322], dat_velocity$Green[1:1322]))

p_velocity <- dat_velocity %>% 
  pivot_longer(-(Time)) %>% 
  ggplot(., aes(Time, value, colour = name))+
  geom_line()

## for day 6
dat_cell_red <- read.csv("~cchen2/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/nvue/IT_PT/m18/m18_d6_redtrace.csv", skip = 1) %>% 
  as_tibble() %>% 
  select(!starts_with("r")) %>% 
  rename(Time = names(.)[1]) %>% 
  mutate_at(vars(-("Time")),scale) %>% 
  mutate(mean = rowMeans(select(.,starts_with("a")))) %>% 
  select(Time, mean)%>% 
  slice(1:(n()-1))



dat_cell_green <- read.csv("~cchen2/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/nvue/IT_PT/m18/m18_d6_Greentrace.csv", skip = 1) %>% 
  as_tibble() %>% 
  select(!starts_with("r")) %>% 
  rename(Time = names(.)[1]) %>% 
  mutate_at(vars(-("Time")),scale) %>% 
  mutate(mean = rowMeans(select(.,starts_with("a")))) %>% 
  select(Time, mean) 
cor(dat_cell_red$mean, dat_cell_green$mean)


dat_cell_nvue_d6 <- rbind(dat_cell_green, dat_cell_red) %>% 
  mutate(Group = c(rep("IT", nrow(dat_cell_green)), rep("PT", nrow(dat_cell_red)))) %>% 
  ggplot(., aes(Time, mean, colour = Group))+
  geom_line()+
  labs(x = "Time (s)", y = "z-Score", title = "D6")+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  theme(legend.title =  element_blank())

r_cross <- ccf(dat_cell_red$mean,  dat_cell_green$mean, plot = F)

p_auto_correlation <- tibble(lag = r_cross$lag, r = r_cross$acf) %>% 
  ggplot(., aes(lag, r))+
  geom_line()

## correlation between velocity and calcium transit--
dat_velocity <- read.xlsx("~cchen2/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/nvue/IT_PT/m18_d6_velocity.xlsx") %>% 
  as_tibble() %>% 
  select(Time, Velocity) %>% 
  slice(which(row_number() %% 2 == 1)) %>% 
  mutate_at(vars(-("Time")),scale) %>%
  mutate(Green = dat_cell_green$mean[1:1371], red = dat_cell_red$mean[1:1371])

cor.test(dat_velocity$Velocity, dat_velocity$Green)
cor.test(dat_velocity$Velocity, dat_velocity$red)
print(ccf(dat_velocity$Velocity[1:1371], dat_velocity$Green[1:1371]))
p_velocity <- dat_velocity %>% 
  pivot_longer(-(Time)) %>% 
  mutate(name = factor(name, levels = c("Green", "red", "Velocity"))) %>% 
  ggplot(., aes(Time, value, colour = name))+
  geom_line()+
  scale_color_manual(values=c("darkcyan", "indianred", "blue"))+
  labs(x = "Time (s)", y = "z-Score")+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  theme(legend.title =  element_blank())


## analyze the trace when mouse crossing border
time_crossing <- which.min(abs(dat_cell_green$Time - 217.2))

dat_m18_d6 <- rbind(dat_cell_green[(time_crossing-40):(time_crossing+80),], dat_cell_red[(time_crossing-40):(time_crossing+80),]) %>% 
  as_tibble() %>% 
  mutate(Group = rep(c("Green", "red"), each = 121 )) %>% 
  ggplot(., aes(Time, mean, colour = Group))+
  geom_line()

## for day 7
dat_cell_red <- read.csv("~cchen2/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/nvue/IT_PT/m18/m18_d7_redtrace.csv", skip = 1) %>% 
  as_tibble() %>% 
  select(!starts_with("r")) %>% 
  rename(Time = names(.)[1]) %>% 
  mutate_at(vars(-("Time")),scale) %>% 
  mutate(mean = rowMeans(select(.,starts_with("a")))) %>% 
  select(Time, mean)%>% 
  slice(1:(n()-1))



dat_cell_green <- read.csv("~cchen2/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/nvue/IT_PT/m18/m18_d7_Greentrace.csv", skip = 1) %>% 
  as_tibble() %>% 
  select(!starts_with("r")) %>% 
  rename(Time = names(.)[1]) %>% 
  mutate_at(vars(-("Time")),scale) %>% 
  mutate(mean = rowMeans(select(.,starts_with("a")))) %>% 
  select(Time, mean) 

cor(dat_cell_red$mean, dat_cell_green$mean)


dat_cell_nvue_d7 <- rbind(dat_cell_green, dat_cell_red) %>% 
  mutate(Group = c(rep("IT", nrow(dat_cell_green)), rep("PT", nrow(dat_cell_red)))) %>% 
  ggplot(., aes(Time, mean, colour = Group))+
  geom_line()+
  labs(x = "Time (s)", y = "z-Score", title = "D7")+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  theme(legend.title =  element_blank())

r_cross <- ccf(dat_cell_red$mean,  dat_cell_green$mean, plot = F)

p_auto_correlation <- tibble(lag = r_cross$lag, r = r_cross$acf) %>% 
  ggplot(., aes(lag, r))+
  geom_line()

p_m18 <- plot_grid(dat_cell_nvue_d1, dat_cell_nvue_d2, dat_cell_nvue_d3, dat_cell_nvue_d4, dat_cell_nvue_d5,dat_cell_nvue_d6,dat_cell_nvue_d7, ncol = 3)


## for m19  

dat_cell_red <- read.csv("~cchen2/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/nvue/IT_PT/m19/m19_d1_RedTraces_11_43.csv", skip = 1) %>% 
  as_tibble() %>% 
  select(!starts_with("r")) %>% 
  rename(Time = names(.)[1]) %>% 
  mutate_at(vars(-("Time")),scale) %>% 
  mutate(mean = rowMeans(select(.,starts_with("a")))) %>% 
  select(Time, mean)



dat_cell_green <- read.csv("~cchen2/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/nvue/IT_PT/m19/m19_d1_GreenTraces_11_43.csv", skip = 1) %>% 
  as_tibble() %>% 
  select(!starts_with("r")) %>% 
  rename(Time = names(.)[1]) %>% 
  mutate_at(vars(-("Time")),scale) %>% 
  mutate(mean = rowMeans(select(.,starts_with("a")))) %>% 
  select(Time, mean)

cor(dat_cell_red$mean, dat_cell_green$mean)
dat_cell_nvue <- rbind(dat_cell_green, dat_cell_red) %>% 
  mutate(Group = c(rep("IT", nrow(dat_cell_green)), rep("PT", nrow(dat_cell_red)))) %>% 
  mutate(Group = factor(Group, levels = c("IT", "PT"))) %>% 
  ggplot(., aes(Time, mean, colour = Group))+
  geom_line()+
  scale_color_manual(values=c("darkcyan", "indianred"))+
  labs(x = "Time (s)", y = "z-Score")+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  theme(legend.title =  element_blank())

## for day 2
dat_cell_red <- read.csv("~cchen2/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/nvue/IT_PT/m19/m19_d2_Redtrace.csv", skip = 1) %>% 
  as_tibble() %>% 
  select(!starts_with("r")) %>% 
  rename(Time = names(.)[1]) %>% 
  mutate_at(vars(-("Time")),scale) %>% 
  mutate(mean = rowMeans(select(.,starts_with("a")))) %>% 
  select(Time, mean)



dat_cell_green <- read.csv("~cchen2/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/nvue/IT_PT/m19/m19_d2_Greentrace.csv", skip = 1) %>% 
  as_tibble() %>% 
  select(!starts_with("r")) %>% 
  rename(Time = names(.)[1]) %>% 
  mutate_at(vars(-("Time")),scale) %>% 
  mutate(mean = rowMeans(select(.,starts_with("a")))) %>% 
  select(Time, mean)

cor(dat_cell_red$mean, dat_cell_green$mean)


dat_cell_nvue <- rbind(dat_cell_green, dat_cell_red) %>% 
  mutate(Group = c(rep("IT", nrow(dat_cell_green)), rep("PT", nrow(dat_cell_red)))) %>% 
  ggplot(., aes(Time, mean, colour = Group))+
  geom_line()+
  labs(x = "Time (s)", y = "z-Score")+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  theme(legend.title =  element_blank())

## for day 3
dat_cell_red <- read.csv("~cchen2/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/nvue/IT_PT/m19/m19_d3_Redtrace.csv", skip = 1) %>% 
  as_tibble() %>% 
  select(!starts_with("r")) %>% 
  rename(Time = names(.)[1]) %>% 
  mutate_at(vars(-("Time")),scale) %>% 
  mutate(mean = rowMeans(select(.,starts_with("a")))) %>% 
  select(Time, mean)



dat_cell_green <- read.csv("~cchen2/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/nvue/IT_PT/m19/m19_d3_Greentrace.csv", skip = 1) %>% 
  as_tibble() %>% 
  select(!starts_with("r")) %>% 
  rename(Time = names(.)[1]) %>% 
  mutate_at(vars(-("Time")),scale) %>% 
  mutate(mean = rowMeans(select(.,starts_with("a")))) %>% 
  select(Time, mean)

cor(dat_cell_red$mean, dat_cell_green$mean)


dat_cell_nvue <- rbind(dat_cell_green, dat_cell_red) %>% 
  mutate(Group = c(rep("IT", nrow(dat_cell_green)), rep("PT", nrow(dat_cell_red)))) %>% 
  ggplot(., aes(Time, mean, colour = Group))+
  geom_line()+
  labs(x = "Time (s)", y = "z-Score")+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  theme(legend.title =  element_blank())
## for day 4
dat_cell_red <- read.csv("~cchen2/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/nvue/IT_PT/m19/m19_d4_Redtrace.csv", skip = 1) %>% 
  as_tibble() %>% 
  select(!starts_with("r")) %>% 
  rename(Time = names(.)[1]) %>% 
  mutate_at(vars(-("Time")),scale) %>% 
  mutate(mean = rowMeans(select(.,starts_with("a")))) %>% 
  select(Time, mean)



dat_cell_green <- read.csv("~cchen2/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/nvue/IT_PT/m19/m19_d4_Greentrace.csv", skip = 1) %>% 
  as_tibble() %>% 
  select(!starts_with("r")) %>% 
  rename(Time = names(.)[1]) %>% 
  mutate_at(vars(-("Time")),scale) %>% 
  mutate(mean = rowMeans(select(.,starts_with("a")))) %>% 
  select(Time, mean)

cor(dat_cell_red$mean, dat_cell_green$mean)


dat_cell_nvue <- rbind(dat_cell_green, dat_cell_red) %>% 
  mutate(Group = c(rep("IT", nrow(dat_cell_green)), rep("PT", nrow(dat_cell_red)))) %>% 
  ggplot(., aes(Time, mean, colour = Group))+
  geom_line()+
  labs(x = "Time (s)", y = "z-Score")+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  theme(legend.title =  element_blank())

## for day 5
dat_cell_red <- read.csv("~cchen2/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/nvue/IT_PT/m19/m19_d5_Redtrace.csv", skip = 1) %>% 
  as_tibble() %>% 
  select(!starts_with("r")) %>% 
  rename(Time = names(.)[1]) %>% 
  mutate_at(vars(-("Time")),scale) %>% 
  mutate(mean = rowMeans(select(.,starts_with("a")))) %>% 
  select(Time, mean)



dat_cell_green <- read.csv("~cchen2/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/nvue/IT_PT/m19/m19_d5_Greentrace.csv", skip = 1) %>% 
  as_tibble() %>% 
  select(!starts_with("r")) %>% 
  rename(Time = names(.)[1]) %>% 
  mutate_at(vars(-("Time")),scale) %>% 
  mutate(mean = rowMeans(select(.,starts_with("a")))) %>% 
  select(Time, mean)

cor(dat_cell_red$mean, dat_cell_green$mean)


dat_cell_nvue <- rbind(dat_cell_green, dat_cell_red) %>% 
  mutate(Group = c(rep("IT", nrow(dat_cell_green)), rep("PT", nrow(dat_cell_red)))) %>% 
  ggplot(., aes(Time, mean, colour = Group))+
  geom_line()+
  labs(x = "Time (s)", y = "z-Score")+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  theme(legend.title =  element_blank())

## for day 6
dat_cell_red <- read.csv("~cchen2/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/nvue/IT_PT/m19/m19_d6_Redtrace.csv", skip = 1) %>% 
  as_tibble() %>% 
  select(!starts_with("r")) %>% 
  rename(Time = names(.)[1]) %>% 
  mutate_at(vars(-("Time")),scale) %>% 
  mutate(mean = rowMeans(select(.,starts_with("a")))) %>% 
  select(Time, mean)



dat_cell_green <- read.csv("~cchen2/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/nvue/IT_PT/m19/m19_d6_Greentrace.csv", skip = 1) %>% 
  as_tibble() %>% 
  select(!starts_with("r")) %>% 
  rename(Time = names(.)[1]) %>% 
  mutate_at(vars(-("Time")),scale) %>% 
  mutate(mean = rowMeans(select(.,starts_with("a")))) %>% 
  select(Time, mean)

cor(dat_cell_red$mean, dat_cell_green$mean)


dat_cell_nvue_d6 <- rbind(dat_cell_green, dat_cell_red) %>% 
  mutate(Group = c(rep("IT", nrow(dat_cell_green)), rep("PT", nrow(dat_cell_red)))) %>% 
  ggplot(., aes(Time, mean, colour = Group))+
  geom_line()+
  labs(x = "Time (s)", y = "z-Score")+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  theme(legend.title =  element_blank())

## for day 7
dat_cell_red <- read.csv("~cchen2/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/nvue/IT_PT/m19/m19_d7_Redtrace.csv", skip = 1) %>% 
  as_tibble() %>% 
  select(!starts_with("r")) %>% 
  rename(Time = names(.)[1]) %>% 
  mutate_at(vars(-("Time")),scale) %>% 
  mutate(mean = rowMeans(select(.,starts_with("a")))) %>% 
  select(Time, mean)



dat_cell_green <- read.csv("~cchen2/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/nvue/IT_PT/m19/m19_d7_Greentrace.csv", skip = 1) %>% 
  as_tibble() %>% 
  select(!starts_with("r")) %>% 
  rename(Time = names(.)[1]) %>% 
  mutate_at(vars(-("Time")),scale) %>% 
  mutate(mean = rowMeans(select(.,starts_with("a")))) %>% 
  select(Time, mean)

cor(dat_cell_red$mean, dat_cell_green$mean)


dat_cell_nvue_d7 <- rbind(dat_cell_green, dat_cell_red) %>% 
  mutate(Group = c(rep("IT", nrow(dat_cell_green)), rep("PT", nrow(dat_cell_red)))) %>% 
  ggplot(., aes(Time, mean, colour = Group))+
  geom_line()+
  labs(x = "Time (s)", y = "z-Score")+
  theme(axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(family = "Arial",size = 12, face ="plain"))+
  theme(legend.title =  element_blank())


## analyze inter- and intra-correlation (Green vs Green and Red vs Red and Red vs Green)----
# correlation of 2 s before the point and 5 s after event time point
file_path <- "~cchen2/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/nvue/IT_PT/m19/"
cc_cell_cor <- function(file_path){
  file <- list.files(file_path, pattern = ".csv", full.names = T)
  ## for the green trace
  
  file_green_event <- file[grep("Green", file)] %>% 
    .[grep("event", .)]
  file_green_trace <- file[grep("Green", file)] %>% 
    setdiff(., file_green_event)
  
  ## for the red trace
  file_red_event <- file[grep("Red", file)] %>% 
    .[grep("event", .)]
  file_red_trace <- file[grep("Red", file)] %>% 
    setdiff(., file_red_event)
  
  ## get the cell ID which has events be detected
  green_cell_ID <- unique(read_csv(file_green_event)$`Cell Name`)
  red_cell_ID <- unique(read_csv(file_red_event)$`Cell Name`)
  
  ## keep cells which have detected events
  dat_cell_green <- read.csv(file_green_trace) %>% 
    as_tibble() %>% 
    slice(-1) %>% 
    rename(Time = X) %>% 
    select(Time, green_cell_ID) %>% 
    mutate_if(is.character,as.numeric) %>% 
    mutate_at(vars(-("Time")),scale)
  
  dat_cell_red <- read.csv(file_red_trace) %>% 
    as_tibble() %>% 
    slice(-1) %>% 
    rename(Time = X) %>% 
    select(Time, red_cell_ID) %>% 
    mutate_if(is.character,as.numeric) %>% 
    mutate_at(vars(-("Time")),scale)


  ## correlation analysis of one cell between others in Green cells group
  # split traces based on time of event
  dat_cor_inter <- c()
  dat_cor_intra <- c()
  for (i in seq_along(red_cell_ID)){
    time_point <- read_csv(file_red_event) %>% 
      as_tibble() %>% 
      rename(Time = 'Time (s)', Cell_ID = 'Cell Name' ) %>% 
      filter(Cell_ID == red_cell_ID[i]) %>% 
      .$Time
    # time point to spilt data
    
    for (j in seq_along(time_point)){
      slice_index <- which(dat_cell_red$Time == time_point[j])
      # correlation of cells between each other cells
      cor_value_red <- dat_cell_red %>% 
        slice((slice_index-20):(slice_index + 30)) %>% 
        mutate_if(is.character,as.numeric) %>% 
        select(-Time) %>% 
        as.matrix() %>% 
        cor %>% 
        .[,i] %>% 
        unname()
      
      cor_value_red <- cor_value_red[cor_value_red!= 1]
      
      dat_cor_inter <- c(dat_cor_inter, cor_value_red)
      
      ## correlation analysis of cell between cells in green group
      cell_ID_trace <- dat_cell_red %>% 
        slice((slice_index-10):(slice_index + 25)) %>% 
        mutate_if(is.character,as.numeric) %>% 
        select(red_cell_ID[i]) %>% 
        rename_with(., tolower)
      
      cor_value_green <- dat_cell_green %>% 
        slice((slice_index-10):(slice_index + 25)) %>% 
        mutate_if(is.character,as.numeric) %>% 
        select(-Time) %>% 
        add_column(cell_ID_trace, .before = 1) %>% 
        as.matrix() %>% 
        cor %>% 
        .[,i] %>% 
        unname() 
        
      cor_value_green <- cor_value_red[cor_value_green!= 1]
      dat_cor_intra <- c(dat_cor_intra, cor_value_green)
        
      
    }
  }
}


## test the activation of PT neurons on IT neuron ----

## for 50hz
trace_info <- read_csv("~cchen2/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/nvue/IT_PT/IT_PT_opto/stimulation_info.csv") %>% 
  rename(Time = names(.)[1], Channel = names(.)[2] ) %>% 
  filter(Channel == "OG-LED") 

## get the time when the stimulaion start
stim_time <- trace_info$Time[which(diff(trace_info$Time) >10) + 1]


## import the trace and truncated by stim time
dat_trace_PT_stim <- read.csv("~cchen2/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/nvue/IT_PT/IT_PT_opto/pn_it_gcamp_pt_opto_50hz_test_cell_trace.csv", skip = 1) %>% 
  as_tibble() %>% 
  select(!starts_with("r")) %>% 
  rename(Time = names(.)[1]) %>% 
  # mutate_at(vars(-("Time")),scale) %>% 
  mutate(mean = rowMeans(select(.,starts_with("a")))) %>% 
  select(Time, mean)

stim_time_trace <- sapply(stim_time, function(x) which.min(abs(dat_trace_PT_stim$Time - x)))

dat_trace_PT_stim_group_50hz <- sapply(stim_time_trace, function(x) dat_trace_PT_stim$mean[(x-20): (x+100)] - mean(dat_trace_PT_stim$mean[(x-20):x])) %>% 
  as_tibble() %>% 
  mutate(Time = seq(-2, 10, 0.1)) %>% 
  select_if(~ !any(is.na(.))) %>% 
  mutate(mean = rowMeans(select(.,starts_with("V")))) %>% 
  select(Time, mean)

## for 20hz
trace_info <- read_csv("~cchen2/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/nvue/IT_PT/IT_PT_opto/stimulation_info_20hz.csv") %>% 
  rename(Time = names(.)[1], Channel = names(.)[2] ) %>% 
  filter(Channel == "OG-LED") 

## get the time when the stimulaion start
stim_time <- trace_info$Time[which(diff(trace_info$Time) >10) + 1]


## import the trace and truncated by stim time
dat_trace_PT_stim <- read.csv("~cchen2/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/nvue/IT_PT/IT_PT_opto/pn_it_gcamp_pt_opto_20hz_test_cell_trace.csv", skip = 1) %>% 
  as_tibble() %>% 
  select(!starts_with("r")) %>% 
  rename(Time = names(.)[1]) %>% 
  # mutate_at(vars(-("Time")),scale) %>% 
  mutate(mean = rowMeans(select(.,starts_with("a")))) %>% 
  select(Time, mean)

## plot the mean of the trace
p_IT_trace <- ggplot(dat_trace_PT_stim, aes(Time, mean))+
  geom_line()
stim_time_trace <- sapply(stim_time, function(x) which.min(abs(dat_trace_PT_stim$Time - x)))

dat_trace_PT_stim_group_20hz <- sapply(stim_time_trace, function(x) dat_trace_PT_stim$mean[(x-20): (x+100)] - mean(dat_trace_PT_stim$mean[(x-20):x])) %>% 
  as_tibble() %>% 
  mutate(Time = seq(-2, 10, 0.1)) %>% 
  select_if(~ !any(is.na(.))) %>% 
  mutate(mean = rowMeans(select(.,starts_with("V")))) %>% 
  select(Time, mean)

## for 10hz
trace_info <- read_csv("~cchen2/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/nvue/IT_PT/IT_PT_opto/stimulation_info_10hz.csv") %>% 
  rename(Time = names(.)[1], Channel = names(.)[2] ) %>% 
  filter(Channel == "OG-LED") 

## get the time when the stimulaion start
stim_time <- trace_info$Time[which(diff(trace_info$Time) >10) + 1]


## import the trace and truncated by stim time
dat_trace_PT_stim <- read.csv("~cchen2/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/nvue/IT_PT/IT_PT_opto/pn_it_gcamp_pt_opto_10hz_test_cell_trace.csv", skip = 1) %>% 
  as_tibble() %>% 
  select(!starts_with("r")) %>% 
  rename(Time = names(.)[1]) %>% 
  # mutate_at(vars(-("Time")),scale) %>% 
  mutate(mean = rowMeans(select(.,starts_with("a")))) %>% 
  select(Time, mean)

stim_time_trace <- sapply(stim_time, function(x) which.min(abs(dat_trace_PT_stim$Time - x)))

dat_trace_PT_stim_group_10hz <- sapply(stim_time_trace, function(x) dat_trace_PT_stim$mean[(x-20): (x+100)] - mean(dat_trace_PT_stim$mean[(x-20):x])) %>% 
  as_tibble() %>% 
  mutate(Time = seq(-2, 10, 0.1)) %>% 
  select_if(~ !any(is.na(.))) %>% 
  mutate(mean = rowMeans(select(.,starts_with("V")))) %>% 
  select(Time, mean)

p_PT_stim_trace <- rbind(dat_trace_PT_stim_group_10hz,dat_trace_PT_stim_group_20hz, dat_trace_PT_stim_group_50hz) %>% 
  mutate(Group = rep(c("10hz","20hz", "50hz"), each = nrow(dat_trace_PT_stim_group_20hz))) %>% 
  ggplot(., aes(Time, mean, color = Group))+
  geom_line()+
  theme_classic()+
  labs(x = "Time (s)", y = "z-score")+
  geom_hline(yintercept = 0, linetype = 2,color = "red")+
  theme(legend.title = element_blank())

setwd("~cchen2/Documents/neuroscience/Pn\ project/Figure/PDF/")
cairo_pdf("p_PT_stim.pdf", width = 110/25.6, height = 60/25.6, family = "Arial")
p_PT_stim_trace
dev.off()  
  

