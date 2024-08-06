rm(list = ls())
library(deSolve)   
library(ggplot2)   
library(dplyr)
library(tidyr)
library(scales)

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
  dV2dt <- r2 * V2 * (1 - (V2 + alpha12 * V1) / K) + mu * V1
  
  return(list(c(dV1dt, dV2dt)))
}


withinhost_fun = function(param_sets){
  # Initial state and time sequence
  state <- c(V1 = 1, V2 = 0)
  times <- seq(0, 100, by = 1)
  
  # Run simulations for each parameter set
  results <- lapply(seq(nrow(param_sets)), function(i) {
    params <- unlist(param_sets[i, ])
    out <- ode(y = state, times = times, func = viral_model, parms = params)
    out_df <- as.data.frame(out)
    out_df$mu <- params['mu']
    return(out_df)
  })
  
  # Combine all dataframes into one
  combined_results <- bind_rows(results)
  
  return(combined_results)
  
}

ratio_fun = function(combined_results){
  long_data_transformed <- combined_results %>%
    mutate(ratio = V1/(V1+V2), group = format(mu, digits = 2)) %>%
    select(time, group, ratio)
  data_ratio = long_data_transformed[long_data_transformed$time %in%
                                       c(0,24,48,72),]
  return(data_ratio)
}

getplot = function(combined_results){

  data_ratio = ratio_fun(combined_results)
  data_ratio$group = factor(data_ratio$group,
                            levels = as.character(format(mu_vec, digits = 2)))
  
  values = c(alpha('#6A619F',0.9), alpha('#2A419F',0.8), 
             alpha('#4A71BF',0.7), alpha('#9AC1FF',0.7))
  values = c(alpha('#2A419F',0.8), alpha('#9AC1FF',0.7))
  values = c(alpha('#6A619F',0.8), alpha('#AA712F',0.7))
  # show_col(values)
  p = ggplot(data_ratio, 
             aes(x = time, y = ratio, 
                 color = factor(group))) +
    geom_line() +
    geom_point() +
    scale_color_manual(values = values,
                       name = expression(mu)) +
    labs(y = "", x = "Time unit", color = "r2 value") +
    theme_bw() +
    theme(legend.position = "none") +
    scale_x_continuous(breaks = c(0,24,48,72)) +
    scale_y_continuous(limits = c(0, 1),
                       minor_breaks = seq(0 , 1, 0.25),
                       n.breaks = 3)

  return(p)
}
transform_data = function(combined_results){
  # Transform data for plotting
  long_data <- pivot_longer(combined_results, cols = c("V1", "V2"), names_to = "Strain", values_to = "Population")
  long_data$Population_label = long_data$Population/10^4
  long_data$group = long_data$r2 - long_data$r1
  long_data$group[long_data$Strain == 'V1'] = 'V1'
  long_data$group = factor(long_data$group, levels = unique(long_data$group))
  return(long_data) 
}

getplot2 = function(combined_results){
  values = c(hue_pal()(3)[1], hue_pal()(3)[3])
  long_data = pivot_longer(combined_results, cols = c("V1", "V2"), names_to = "Strain", values_to = "Population")
  # long_data = long_data[long_data$mu == mu_vec[2],]
  long_data$Population = long_data$Population/max(long_data$Population)
  v = c("#F8766D", "#619CFF")
  values = c(alpha(v,0.9), alpha(v,0.6),
             alpha(v,0.4), alpha(v,0.2))
  p = ggplot() + 
    geom_line(data = long_data, 
              aes(x = time, y = Population,
                  group = interaction(Strain, mu), color = interaction(Strain, mu))) +
    theme_bw() + 
    scale_color_manual(values = values) +
    xlab('Time unit') + 
    ylab('') +
    theme(legend.position = "none") +
    scale_x_continuous(breaks = c(0,24,48,72)) +
    scale_y_continuous(limits = c(0, 1),
                       minor_breaks = seq(0 , 1, 0.25),
                       breaks = c(0,0.5,1),
                       labels = c('0.0','0.5','1.0'))
  return(p)
}

Kbase = 200
mu = 10^(-6)
mu_vec = c(10^(-6), 10^(-3))
r = 0.2
param_sets1 <- expand.grid(r1 = r, r2 = r, K = 200, 
                           alpha12 = 1, alpha21 = 1, 
                           mu = mu_vec)
param_sets2 <- expand.grid(r1 = r, r2 = r + 0.04, K = 200, 
                          alpha12 = 1, alpha21 = 1, 
                          mu = mu_vec)
param_sets3 <- expand.grid(r1 = r, r2 = r + 0.4, K = 200, 
                          alpha12 = 1, alpha21 = 1, 
                          mu = mu_vec)
param_sets4 <- expand.grid(r1 = r, r2 = r, K = 200, 
                          alpha12 = 0, alpha21 = 1, 
                          mu = mu_vec)

pdf(paste0("Output/withinhost_evolution_mutation_rate.pdf"), width = 1.5, height = 1.2)
print(getplot(withinhost_fun(param_sets1)))
print(getplot(withinhost_fun(param_sets2)))
print(getplot(withinhost_fun(param_sets3)))
print(getplot(withinhost_fun(param_sets4)))
print(getplot2(withinhost_fun(param_sets1)))
print(getplot2(withinhost_fun(param_sets2)))
print(getplot2(withinhost_fun(param_sets3)))
print(getplot2(withinhost_fun(param_sets4)))
dev.off()

pdf(paste0("Output/withinhost_evolution_legend.pdf"), width = 3, height = 1.2)
{
  combined_results = withinhost_fun(param_sets1)
  data_ratio = ratio_fun(combined_results)
  data_ratio$group = factor(data_ratio$group,
                            levels = as.character(format(mu_vec, digits = 2)))
  
  values = c(alpha('#6A619F',0.8), alpha('#AA712F',0.7))
  # show_col(values)
  p = ggplot(data_ratio, 
             aes(x = time, y = ratio, 
                 color = factor(group))) +
    geom_line() +
    geom_point() +
    scale_color_manual(values = values,
                       name = expression(mu)) +
    labs(y = "", x = "Time unit", color = "r2 value") +
    theme_bw() +
    theme(legend.position = "top") +
    scale_x_continuous(breaks = c(0,24,48,72)) +
    scale_y_continuous(limits = c(0, 1),
                       minor_breaks = seq(0 , 1, 0.25),
                       n.breaks = 3)
  print(p)
}
dev.off()

param_sets2 <- expand.grid(r1 = r, r2 = r, K = 200, 
                          alpha12 = 0, alpha21 = 0, 
                          mu = mu_vec)
param_sets3 <- expand.grid(r1 = r, r2 = r, K = 200, 
                           alpha12 = 1, alpha21 = 0, 
                           mu = mu_vec)




r = 0.3
param_sets2 <- expand.grid(r1 = r, r2 = r, K = 200, 
                           alpha12 = 1, alpha21 = 1, 
                           mu = mu_vec)
results2 = withinhost_fun(param_sets2)
getplot(results2)
getplot2(results2)



param_sets2 <- expand.grid(r1 = r, r2 = r, K = 100000, 
                          alpha12 = 0, alpha21 = 1, 
                          mu = mu_vec)
param_sets3 <- expand.grid(r1 = r, r2 = r, K = 100000, 
                          alpha12 = 1, alpha21 = 0, 
                          mu = mu_vec)
param_list = list(param_sets1, param_sets2, param_sets3)

results_list1 = list()
for (i in 1:3) {
  results_list1[[i]] = withinhost_fun(param_list[[i]])
}

r = 0.4
param_sets1 <- expand.grid(r1 = r, r2 = r, K = 100000, 
                           alpha12 = 0, alpha21 = 0, 
                           mu = mu_vec)
param_sets2 <- expand.grid(r1 = r, r2 = r, K = 100000, 
                           alpha12 = 0, alpha21 = 2, 
                           mu = mu_vec)
param_sets3 <- expand.grid(r1 = r, r2 = r, K = 100000, 
                           alpha12 = 2, alpha21 = 0, 
                           mu = mu_vec)
param_list = list(param_sets1, param_sets2, param_sets3)

results_list2 = list()
for (i in 1:3) {
  results_list2[[i]] = withinhost_fun(param_list[[i]])
}


pdf(paste0("Output/withinhost_evolution1.pdf"), width = 1.5, height = 1.2)
results_list = results_list1
getplot(results_list[[1]])
getplot(results_list[[2]])
getplot(results_list[[3]])
getplot2(results_list[[1]])
getplot2(results_list[[2]])
getplot2(results_list[[3]])
dev.off()

pdf(paste0("Output/withinhost_evolution2.pdf"), width = 1.5, height = 1.2)
results_list = results_list2
getplot(results_list[[1]])
getplot(results_list[[2]])
getplot(results_list[[3]])
getplot2(results_list[[1]])
getplot2(results_list[[2]])
getplot2(results_list[[3]])
dev.off()


# SARS-CoV-2 mutation rate estimates of around 
# 1 × 10–6–2 × 10–6 mutations per nucleotide per replication cycle 

data_ratio = ratio_fun(results_list[[1]])
data_ratio$group = factor(data_ratio$group,
                          levels = as.character(format(mu_vec, digits = 2)))

getplot_legend = function(){

  values = c(alpha('#6A619F',0.9), alpha('#2A419F',0.8), 
             alpha('#4A71BF',0.7), alpha('#9AC1FF',0.7))

  p = ggplot(data_ratio, 
             aes(x = time, y = ratio, 
                 color = factor(group))) +
    geom_line() +
    geom_point() +
    scale_color_manual(values = values,
                       name = expression(mu),
                       label = c('1e-12','1e-9','1e-6','1e-3')) +
    labs(y = "", x = "Time unit") +
    theme_bw() +
    theme(legend.position = "top",
          legend.key.size = unit(0.4,'cm'),
          legend.key.width = unit(0.4,'cm'),
          legend.key = element_blank(),
          legend.background = element_blank(),
          panel.background = element_blank())
  
  return(p)
}
# v = c('#2A41AF', "#F8766D", "#619CFF")
# getplot_legend(v[3])

pdf(paste0("Output/withinhost_evolution_legend.pdf"), 
    width = 3, height = 1.4)
print(getplot_legend())
dev.off()
# 2*10^(-6)
