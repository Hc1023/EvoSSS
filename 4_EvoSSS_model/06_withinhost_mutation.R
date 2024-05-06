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
  times <- seq(0, 72, by = 0.1)
  
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


getplot = function(param_sets){
  combined_results = withinhost_fun(param_sets)
  data_ratio = ratio_fun(combined_results)
  data_ratio$group = factor(data_ratio$group,
                            levels = as.character(format(mu_vec, digits = 2)))
  
  values = alpha(c('#2A41AF', '#4D7AAF', '#619CAF'),0.7)
  
  
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
rbase = log2(exp(1))/3
Kbase_vec = c(2^20, 2^12)
Kbase = c(2^18)
mu = (2*10^(-6))/(1/rbase)
mu_vec = c(mu/2, mu, mu*10^3)
param_sets <- expand.grid(r1 = rbase, r2 = rbase, K = Kbase, 
                          alpha12 = 0, alpha21 = 0, 
                          mu = mu_vec)
p1 = getplot(param_sets)
param_sets <- expand.grid(r1 = rbase, r2 = rbase, K = Kbase, 
                          alpha12 = 0, alpha21 = 1, 
                          mu = mu_vec)
p2 = getplot(param_sets)
param_sets <- expand.grid(r1 = rbase, r2 = rbase, K = Kbase, 
                          alpha12 = 1, alpha21 = 0, 
                          mu = mu_vec)
p3 = getplot(param_sets)

pdf(paste0("Output/within_host_mutations.pdf"), 
    width = 2, height = 1.2)
print(p1)
print(p2)
print(p3)
dev.off()

# SARS-CoV-2 mutation rate estimates of around 
# 1 × 10–6–2 × 10–6 mutations per nucleotide per replication cycle 


# mu = 0.03-0.06/(1/rbase)




# Initial state
state <- c(V1 = 1, V2 = 0)

# Time sequence for the simulation
times <- seq(0, 72, by = 0.1)

# Solve the model
out <- ode(y = state, times = times, func = viral_model, parms = params)
out_df <- as.data.frame(out)
# Plotting
ggplot(data = out_df, aes(x = time)) +
  geom_line(aes(y = V1, color = "Strain 1")) +
  geom_line(aes(y = V2, color = "Strain 2")) +
  labs(y = "Population unit", x = "Time unit", color = "Strain") +
  theme_bw() +
  scale_x_continuous(breaks = c(0,24,48,72)) +
  scale_y_continuous(labels = label_scientific(digits = 3))


