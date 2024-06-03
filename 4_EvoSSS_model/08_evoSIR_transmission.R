rm(list = ls())
library(deSolve)   
library(ggplot2)   
library(dplyr)
library(tidyr)
library(scales)
library(ggnewscale)
# Define the model
viral_model <- function(t, state, parameters) {
  V1 <- state[1]
  V2 <- state[2]
  
  r1 <- parameters["r1"]
  r2 <- parameters["r2"]
  K <- parameters["K"]
  alpha12 <- parameters["alpha12"]
  alpha21 <- parameters["alpha21"]
  mu <- parameters["mu"]
  

  dV1dt <- r1 * V1 * (1 - (V1 + alpha21 * V2) / K) - mu * V1
  dV2dt <- r2 * V2 * (1 - (V2 + alpha12 * V1) / K) - mu * V1
  
  return(list(c(dV1dt, dV2dt)))
}

withinhost_fun = function(param_sets){
  # Initial state and time sequence
  # state <- c(V1 = 1, V2 = 1)
  times <- seq(0, 72, by = 0.1)
  
  # Run simulations for each parameter set
  results <- lapply(seq(nrow(param_sets)), function(i) {
    params <- unlist(param_sets[i, ])
    out <- ode(y = c(params['V1'],params['V2']), times = times, 
               func = viral_model, parms = params)
    out_df <- as.data.frame(out)
    out_df$group <- params['group']
    return(out_df)
  })
  combined_results = bind_rows(results)
  
  return(combined_results)
  
}

transform_data = function(combined_results){
  long_data <- pivot_longer(combined_results, cols = c("V1", "V2"), 
                            names_to = "Strain", values_to = "Population")
  long_data$Population_label = long_data$Population/max(long_data$Population)
  long_data$group = factor(long_data$group, levels = unique(long_data$group))
  return(long_data) 
}
ratio_fun = function(combined_results){
  data_ratio <- combined_results %>%
    mutate(ratio = V1/(V1+V2), group = group) %>%
    select(time, group, ratio)
  return(data_ratio)
}
# Parameters
load(file = 'evoSIR.rdata')
fit = fitlist[[4]]
posterior_samples <- rstan::extract(fit)

c = posterior_samples$beta1/posterior_samples$beta2
deltar = -log(c)

group = c(0.02,0.05,0.1,0.2)
param_sets <- expand.grid(r1 = 0.2, 
                          r2 = 0.2 + group, 
                          K = 5000, 
                          alpha12 = 2, 
                          alpha21 = 0, 
                          mu = 0,
                          V1 = 1, V2 = 1)
param_sets['group'] = group

long_data_list = list()
data_ratio_list = list()

for (i in 1:10) {
  combined_results = withinhost_fun(param_sets)
  long_data = transform_data(combined_results)
  data_ratio = ratio_fun(combined_results)
  long_data_list[[i]] = long_data
  data_ratio_list[[i]] = data_ratio
  param_sets$V1= data_ratio$ratio[data_ratio == 24]*2
  param_sets$V2= (1-data_ratio$ratio[data_ratio == 24])*2
}


data_ratio_combined = data.frame()
long_data_combined = data.frame()
for (i in 1:10) {
  data_ratio = data_ratio_list[[i]] 
  y = data_ratio[data_ratio$time <= 24,]
  y$time = (i-1)*24+y$time
  data_ratio_combined = rbind(data_ratio_combined, y)
  
  long_data = long_data_list[[i]]
  y = long_data[long_data$time <= 24,]
  y$time = (i-1)*24+y$time
  long_data_combined = rbind(long_data_combined, y)
}

v = hue_pal()(3)[1]
values1 = c(alpha(v, 0.9), alpha(v, 0.6), 
            alpha(v, 0.4), alpha(v, 0.2))
v = hue_pal()(3)[3]
values2 = c(alpha(v, 0.9), alpha(v, 0.6), 
            alpha(v, 0.4), alpha(v, 0.2))
v = '#2A41AF'
values = c(alpha(v, 0.9), alpha(v, 0.6), 
           alpha(v, 0.4), alpha(v, 0.2))

p = ggplot() +
  geom_point(data = long_data_combined,
            aes(x = time, y = Population_label,
                group = interaction(group, Strain),
                color = interaction(group, Strain)),
            size = 0.2) +
  scale_color_manual(values = c(values1, values2)) +
  new_scale_color() + 
  geom_line(data = data_ratio_combined[data_ratio_combined$time %in% seq(0,24*10,24),], 
            aes(x = time, y = ratio, color = factor(group))) +
  geom_point(data = data_ratio_combined[data_ratio_combined$time %in% seq(0,24*10,24),], 
             aes(x = time, y = ratio, color = factor(group))) +
  scale_color_manual(values = values, 
                     name = expression(Delta*r)) +
  labs(y = "Population (%)", x = "Transmission", color = "r2 value") +
  theme_bw() +
  theme(legend.position = "none") +
  scale_x_continuous(breaks = seq(0,24*10,24),
                     labels = 0:10) + 
  scale_y_continuous(limits = c(0, 0.5), 
                     breaks = c(0,0.25,0.5), 
                     labels = c(0,25,50),
                     sec.axis = sec_axis(
                       transform = ~ ., 
                       name = 'Ratio (%)',
                       breaks = c(0,0.25,0.5), 
                       labels = c(0,25,50)
                     ))

p
pdf(paste0("Output/withinhost_transmission.pdf"), width = 2.4, height = 1.6)
print(p)
dev.off()


getplot = function(long_data){


  # Plotting
  p1 = ggplot(long_data, aes(x = time, y = Population_label)) +
    geom_line(data = long_data[long_data$Strain == 'V1',], 
              aes(group = group, color = group)) +
    scale_color_manual(values = values1,
                       name = 'A') +
    new_scale_color() +
    geom_line(data = long_data[long_data$Strain == 'V2',], 
              aes(group = group, color = group)) +
    labs(y = "", x = "", color = "Strain") +
    scale_color_manual(values = values2,
                       name = 'B') +
    theme_bw() +
    theme(legend.position = "none") +
    scale_x_continuous(breaks = c(0,24,48,72)) +
    scale_y_continuous(n.breaks = 3, labels = c(0,50,100))
    
  return(p1)
}


pdf("Output/withinhost_transmission_2.pdf", 
    width = 1.5, height = 1)
for(i in c(1,2,5,10)){
  print(getplot(long_data_list[[i]]))
}
dev.off()


