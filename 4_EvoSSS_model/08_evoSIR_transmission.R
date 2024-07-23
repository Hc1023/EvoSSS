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




h_function = function(h){
  long_data_list = list()
  data_ratio_list = list()
  group = c(0.02,0.04,0.1,0.2)
  param_sets <- expand.grid(r1 = 0.2, 
                            r2 = 0.2 + group, 
                            K = 200, 
                            alpha12 = 1, 
                            alpha21 = 0, 
                            mu = 0,
                            V1 = 1, V2 = 1)
  param_sets['group'] = group
  for (i in 1:20) {
    combined_results = withinhost_fun(param_sets)
    long_data = transform_data(combined_results)
    data_ratio = ratio_fun(combined_results)
    long_data_list[[i]] = long_data
    data_ratio_list[[i]] = data_ratio
    param_sets$V1= data_ratio$ratio[data_ratio == h]*2
    param_sets$V2= (1-data_ratio$ratio[data_ratio ==h])*2
  }
  
  data_ratio_combined = data.frame()
  long_data_combined = data.frame()
  for (i in 1:20) {
    data_ratio = data_ratio_list[[i]] 
    y = data_ratio[data_ratio$time <= h,]
    y$time = (i-1)*h+y$time
    data_ratio_combined = rbind(data_ratio_combined, y)
    
    long_data = long_data_list[[i]]
    y = long_data[long_data$time <= h,]
    y$time = (i-1)*h+y$time
    long_data_combined = rbind(long_data_combined, y)
  }
  return(list(data_ratio_combined, long_data_list))
}

getplot1 = function(data_ratio_combined, h){
  values = c(alpha('#6A619F',0.9), alpha('#2A419F',0.8), 
             alpha('#4A71BF',0.7), alpha('#9AC1FF',0.7))
  breaks = c(h*c(0,1,2,3,5,10,20))
  p = ggplot() +
    geom_line(data = data_ratio_combined[data_ratio_combined$time %in% breaks,], 
              aes(x = time, y = ratio, color = factor(group))) +
    geom_point(data = data_ratio_combined[data_ratio_combined$time %in% breaks,], 
               aes(x = time, y = ratio, color = factor(group))) +
    scale_color_manual(values = values, 
                       name = expression(r[2]-r[1])) +
    labs(y = expression('Ratio of V'['1']), x = "Transmission cycle") +
    theme_bw() +
    theme(legend.position = 'right',
          legend.background = element_rect(color = NA, fill = NA),
          legend.key = element_blank(),
          legend.key.height = unit(0.2, 'cm'),
          legend.key.width = unit(0.5, 'cm'),
          panel.grid.minor = element_blank())+
    scale_x_continuous(breaks = breaks,
                       labels = breaks/h) + 
    scale_y_continuous(n.breaks = 3)
  p
  return(p)
}
h_vec = c(12,24,48)
results1 = h_function(h_vec[1])
results2 = h_function(h_vec[2])
results3 = h_function(h_vec[3])

p1 = getplot1(results1[[1]], h_vec[1])
p2 = getplot1(results2[[1]], h_vec[2])
p3 = getplot1(results3[[1]], h_vec[3])

pdf(paste0("Output/withinhost_transmission.pdf"), 
    width = 3, height = 1.5)
print(p1)
print(p2)
print(p3)
dev.off()


getplot2 = function(long_data){
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
    scale_y_continuous(n.breaks = 3, labels = c('0.0','0.5','1.0'))
    
  return(p1)
}

long_data_list = results1[[2]]

pdf("Output/withinhost_transmission_2.pdf", 
    width = 1.5, height = 1)
for(i in c(1,2,5,10)){
  print(getplot2(long_data_list[[i]]))
}
dev.off()


