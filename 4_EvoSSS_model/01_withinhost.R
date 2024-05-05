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

# Parameters
rbase = log2(exp(1))/3
Kbase = 2^12
Kbase = 2^20
params <- c(r1 = rbase, r2 = (rbase+0.05), K = Kbase, 
            alpha12 = 1, alpha21 = 0, 
            mu = 0)

# Initial state
state <- c(V1 = 1, V2 = 1)

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
  long_data$Population_label = long_data$Population/10^4
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
  
  values = alpha(c(hue_pal()(3)[1], '#2A41AF', '#4D7AAF', '#619CAF', '#85BDAF'),0.7)
  # Plotting
  p1 = ggplot(long_data, aes(x = time, y = Population_label)) +
    geom_line(data = long_data, 
              aes(group = group, color = group)) +
    scale_color_manual(values = values, name = expression('Strain'/Delta*r)) +
    labs(y = "", x = "Time unit", color = "Strain") +
    theme_bw() +
    theme(legend.position = "none") +
    scale_x_continuous(breaks = c(0,24,48,72)) +
    scale_y_continuous(n.breaks = 3)

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
                       n.breaks = 3)

  return(list(p1, p2))
}

# Parameter sets for multiple simulations
rbase = log2(exp(1))/3
Kbase_vec = c(2^20, 2^12)
plist_all = list()
for (i in  1:2) {
  Kbase = Kbase_vec[i]
  param_sets <- expand.grid(r1 = rbase, 
                            r2 = rbase + c(0.01, 0.02, 0.05, 0.1), 
                            K = Kbase, 
                            alpha12 = 1, 
                            alpha21 = 0, 
                            mu = 0)
  plist = getplot(param_sets)
  plist_all[[i*2-1]] = plist[[1]]
  plist_all[[i*2]] = plist[[2]]
  
}


pdf(paste0("Output/within_host.pdf"), width = 2, height = 1.2)
print(plist_all[[1]])
print(plist_all[[2]])
print(plist_all[[3]])
print(plist_all[[4]])
dev.off()



# Define a function for finding equilibrium, adjusted for rootSolve
viral_model_equilibrium <- function(x, params) {
  viral_model(t = 0, state = x, parameters = params)[[1]]
}

# Find equilibrium
equilibrium <- multiroot(f = viral_model_equilibrium, start = state, parameters = params)
print(equilibrium$root)  # Display equilibrium points

# Function to wrap viral_model_equilibrium for numDeriv
func_for_jacobian <- function(x) {
  viral_model_equilibrium(x, params)
}

# Compute Jacobian at equilibrium using numDeriv
J_at_eq <- jacobian(func = func_for_jacobian, x = equilibrium$root)
print(J_at_eq)  # Print Jacobian matrix

# Eigenvalues for stability analysis
eigenvalues <- eigen(J_at_eq)$values
print(eigenvalues)  # Print eigenvalues

# Check stability
if (all(Re(eigenvalues) < 0)) {
  cat("The system is stable at the equilibrium.\n")
} else {
  cat("The system is unstable at the equilibrium.\n")
}
