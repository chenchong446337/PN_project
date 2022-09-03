## to analyze the correlation of IT and PT neurons during OP exploring

## compare the cell activity of D1

## For green IT neurons
dat_cell_green_d1 <- read.csv("~cchen/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/nvue/IT_PT/m18/m18_d3_GreenTrace.csv", skip = 1) %>% 
  as_tibble() %>% 
  select(!starts_with("r")) %>% 
  rename(Time = names(.)[1])%>%
  mutate_at(vars(-("Time")),scale) %>% 
  pivot_longer(-Time) %>% 
  mutate(Group = "IT") %>% 
  mutate(name = as.factor(name))


  
dat_cell_red_d1 <- read.csv("~cchen/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/nvue/IT_PT/m18/m18_d3_RedTrace.csv", skip = 1) %>% 
  as_tibble() %>% 
  select(!starts_with("r")) %>% 
  rename(Time = names(.)[1])%>%
  mutate_at(vars(-("Time")),scale) %>% 
  pivot_longer(-Time) %>% 
  mutate(Group = "PT") %>% 
  mutate(name = as.factor(name))


p_trace_green <- dat_cell_green_d1 %>% 
  filter(name %in% sample(levels(dat_cell_green_d1$name), 10)) %>% 
  ggplot(., aes(x = Time, y = value))+
  geom_line(color ="green")+
  facet_grid(rows = vars(name))+
  theme_void()+
  theme(strip.text.y = element_blank())
  

p_trace_red <- dat_cell_red_d1 %>% 
  filter(name %in% sample(levels(dat_cell_red_d1$name), 10)) %>% 
  ggplot(., aes(x = Time, y = value))+
  geom_line(color ="red")+
  facet_grid(rows = vars(name))+
  theme_void()+
  theme(strip.text.y = element_blank())

p_trace <- plot_grid(p_trace_green, p_trace_red, nrow = 2)

## to see how cells activity changes

dat_cell_green_d1 <- read.csv("~cchen/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/nvue/IT_PT/m18/m18_d1_GreenTrace.csv", skip = 1) %>% 
  as_tibble() %>% 
  select(!starts_with("r")) %>% 
  rename(Time = names(.)[1])%>%
  mutate_at(vars(-("Time")),scale) %>% 
  pivot_longer(-Time) %>% 
  mutate(Group = "d1") %>% 
  mutate(name = as.factor(name)) %>% 
  filter(name %in% sample(levels(name), 10))

dat_cell_green_d2 <- read.csv("~cchen/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/nvue/IT_PT/m18/m18_d2_GreenTrace.csv", skip = 1) %>% 
  as_tibble() %>% 
  select(!starts_with("r")) %>% 
  rename(Time = names(.)[1])%>%
  mutate_at(vars(-("Time")),scale) %>% 
  pivot_longer(-Time) %>% 
  mutate(Group = "d2") %>% 
  mutate(name = as.factor(name)) %>% 
  filter(name %in% sample(levels(name), 10)) 

dat_cell_green_d3 <- read.csv("~cchen/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/nvue/IT_PT/m18/m18_d3_GreenTrace.csv", skip = 1) %>% 
  as_tibble() %>% 
  select(!starts_with("r")) %>% 
  rename(Time = names(.)[1])%>%
  mutate_at(vars(-("Time")),scale) %>% 
  pivot_longer(-Time) %>% 
  mutate(Group = "d3") %>% 
  mutate(name = as.factor(name)) %>% 
  filter(name %in% sample(levels(name), 10))

p_green_trace1 <- dat_cell_green_d1 %>% 
  ggplot(., aes(x = Time, y = value))+
  geom_line(color ="black")+
  facet_grid(rows = vars(name))+
  theme_void()+
  theme(strip.text.y = element_blank())


p_green_trace2 <- dat_cell_green_d2 %>% 
  ggplot(., aes(x = Time, y = value))+
  geom_line(color ="red")+
  facet_grid(rows = vars(name))+
  theme_void()+
  theme(strip.text.y = element_blank())

p_green_trace3 <- dat_cell_green_d3 %>% 
  ggplot(., aes(x = Time, y = value))+
  geom_line(color ="green")+
  facet_grid(rows = vars(name))+
  theme_void()+
  theme(strip.text.y = element_blank())

p_it_trace <- plot_grid(p_green_trace1, p_green_trace2, p_green_trace3, nrow = 1)

## for the red cells
dat_cell_red_d1 <- read.csv("~cchen/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/nvue/IT_PT/m18/m18_d1_RedTrace.csv", skip = 1) %>% 
  as_tibble() %>% 
  select(!starts_with("r")) %>% 
  rename(Time = names(.)[1])%>%
  mutate_at(vars(-("Time")),scale) %>% 
  pivot_longer(-Time) %>% 
  mutate(Group = "d1") %>% 
  mutate(name = as.factor(name)) %>% 
  filter(name %in% sample(levels(name), 10))

dat_cell_red_d2 <- read.csv("~cchen/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/nvue/IT_PT/m18/m18_d2_RedTrace.csv", skip = 1) %>% 
  as_tibble() %>% 
  select(!starts_with("r")) %>% 
  rename(Time = names(.)[1])%>%
  mutate_at(vars(-("Time")),scale) %>% 
  pivot_longer(-Time) %>% 
  mutate(Group = "d2") %>% 
  mutate(name = as.factor(name)) %>% 
  filter(name %in% sample(levels(name), 10)) 

dat_cell_red_d3 <- read.csv("~cchen/Documents/neuroscience/Pn\ project/Data_analysis/miniscope/nvue/IT_PT/m18/m18_d3_RedTrace.csv", skip = 1) %>% 
  as_tibble() %>% 
  select(!starts_with("r")) %>% 
  rename(Time = names(.)[1])%>%
  mutate_at(vars(-("Time")),scale) %>% 
  pivot_longer(-Time) %>% 
  mutate(Group = "d3") %>% 
  mutate(name = as.factor(name)) %>% 
  filter(name %in% sample(levels(name), 10))

p_red_trace1 <- dat_cell_red_d1 %>% 
  ggplot(., aes(x = Time, y = value))+
  geom_line(color ="black")+
  facet_grid(rows = vars(name))+
  theme_void()+
  theme(strip.text.y = element_blank())


p_red_trace2 <- dat_cell_red_d2 %>% 
  ggplot(., aes(x = Time, y = value))+
  geom_line(color ="red")+
  facet_grid(rows = vars(name))+
  theme_void()+
  theme(strip.text.y = element_blank())

p_red_trace3 <- dat_cell_red_d3 %>% 
  ggplot(., aes(x = Time, y = value))+
  geom_line(color ="blue")+
  facet_grid(rows = vars(name))+
  theme_void()+
  theme(strip.text.y = element_blank())

p_pt_trace <- plot_grid(p_red_trace1, p_red_trace2, p_red_trace3, nrow = 1)

p_trace <- plot_grid(p_pt_trace, p_it_trace, nrow = 2)

