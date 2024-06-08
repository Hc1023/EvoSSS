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
  dV2dt <- r2 * V2 * (1 - (V2 + alpha12 * V1) / K) - mu * V1
  
  return(list(c(dV1dt, dV2dt)))
}

########### Try different settings ##########

withinhost_fun = function(param_sets){
  # Initial state and time sequence
  state <- c(V1 = 1, V2 = 1)
  times <- seq(0, 72, by = 0.1)
  
  # Run simulations for each parameter set
  results <- lapply(seq(nrow(param_sets)), function(i) {
    params <- unlist(param_sets[i, ])
    out <- ode(y = state, times = times, func = viral_model, parms = params)
    out_df <- as.data.frame(out)
    out_df$r1 <- params['r1']
    out_df$r2 <- params['r2']
    return(out_df)
  })
  
  # Combine all dataframes into one
  combined_results <- bind_rows(results)
  
  return(combined_results)
  
}

transform_data = function(combined_results){
  # Transform data for plotting
  long_data <- pivot_longer(combined_results, cols = c("V1", "V2"), names_to = "Strain", values_to = "Population")
  long_data$Population_label = long_data$Population/10^3
  long_data$group = long_data$r2 - long_data$r1
  long_data$group[long_data$Strain == 'V1'] = 'V1'
  long_data$group = factor(long_data$group, levels = unique(long_data$group))
  return(long_data) 
}

ratio_fun = function(combined_results){
  long_data_transformed <- combined_results %>%
    mutate(ratio = V1/(V1+V2), group = r2-r1) %>%
    select(time, group, ratio)
  data_ratio = long_data_transformed[long_data_transformed$time %in%
                                     c(0,24,48,72),]
  return(data_ratio)
}

getplot = function(param_sets){
  combined_results = withinhost_fun(param_sets)
  long_data = transform_data(combined_results)
  data_ratio = ratio_fun(combined_results)
  v = '#2A41AF'
  values = c(alpha(v, 0.9), alpha(v, 0.6), 
             alpha(v, 0.4), alpha(v, 0.2))
  # values = alpha(c(hue_pal()(3)[1], '#2A41AF', '#4D7AAF', '#619CAF', '#85BDAF'),0.7)
  values = c(hue_pal()(3)[1], values)
  # Plotting
  p1 = ggplot(long_data, aes(x = time, y = Population_label)) +
    geom_line(data = long_data, 
              aes(group = group, color = group)) +
    scale_color_manual(values = values, name = expression('Strain'/Delta*r)) +
    labs(y = "", x = "Time unit", color = "Strain") +
    theme_bw() +
    theme(legend.position = "none") +
    scale_x_continuous(breaks = c(0,24,48,72)) +
    scale_y_continuous(n.breaks = 3, labels = c(0,50,100))

  p2 = ggplot(data_ratio, aes(x = time, y = ratio, color = factor(group))) +
    geom_line() +
    geom_point() +
    scale_color_manual(values = values[-1], 
                       name = expression(Delta*r)) +
    labs(y = "", x = "Time unit", color = "r2 value") +
    theme_bw() +
    theme(legend.position = "none") +
    scale_x_continuous(breaks = c(0,24,48,72)) + 
    scale_y_continuous(limits = c(0, 1), 
                       minor_breaks = seq(0 , 1, 0.25),
                       breaks = c(0,0.5,1), labels = c(0,50,100))

  return(list(p1, p2))
}

# Parameter sets for multiple simulations

param_sets <- expand.grid(r1 = 0.2, 
                          r2 = 0.2 + c(0.02, 0.05, 0.1, 0.2), 
                          K = 5000, alpha12 = 2, alpha21 = 0, mu = 0)


plist1 = getplot(param_sets)
  
param_sets <- expand.grid(r1 = 0.4, 
                          r2 = 0.4 + c(0.02, 0.05, 0.1, 0.2), 
                          K = 5000, alpha12 = 2, alpha21 = 0, mu = 0)
plist2 = getplot(param_sets)

param_sets <- expand.grid(r1 = 0.4, 
                          r2 = 0.4 + c(0.02, 0.05, 0.1, 0.2), 
                          K = 5000, alpha12 = 4, alpha21 = 0, mu = 0)
plist3 = getplot(param_sets)
pdf(paste0("Output/withinhost_competition.pdf"), width = 1.6, height = 1.2)
print(plist1[[1]])
print(plist1[[2]])
print(plist2[[1]])
print(plist2[[2]])
print(plist3[[1]])
print(plist3[[2]])
dev.off()

combined_results = withinhost_fun(param_sets)
data_ratio = ratio_fun(combined_results)
v = '#2A41AF'
values = c(alpha(v, 0.9), alpha(v, 0.6), 
           alpha(v, 0.4), alpha(v, 0.2))
# values = alpha(c(hue_pal()(3)[1], '#2A41AF', '#4D7AAF', '#619CAF', '#85BDAF'),0.7)
values = c(hue_pal()(3)[1], values)
# values = alpha(c(hue_pal()(3)[1], '#2A41AF', '#4D7AAF', '#619CAF', '#85BDAF'),0.7)
p1 = ggplot(data_ratio, aes(x = time, y = ratio, color = factor(group))) +
  geom_line() +
  geom_point() +
  scale_color_manual(values = values[-1], 
                     name = expression(r[B]-r[A])) +
  labs(y = "", x = "Time unit", color = "r2 value") +
  theme_bw() +
  theme(legend.position = "right",
        legend.key.size = unit(0.4,'cm')) +
  scale_x_continuous(breaks = c(0,24,48,72)) + 
  scale_y_continuous(limits = c(0, 1), 
                     minor_breaks = seq(0 , 1, 0.25),
                     n.breaks = 3)

p1

pdf(paste0("Output/withinhost_competition_legend.pdf"), width = 1.8, height = 1.5)
print(p1)
dev.off()
