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
params <- c(r1 = 0.1, r2 = 0.2, K = 100, 
            alpha12 = 0.8, alpha21 = 0.2, 
            mu = 0)

# Initial state
state <- c(V1 = 10, V2 = 10)

# Time sequence for the simulation
times <- seq(0, 100, by = 1)

# Solve the model
out <- ode(y = state, times = times, func = viral_model, parms = params)
out_df <- as.data.frame(out)

# Plotting
ggplot(data = out_df, aes(x = time)) +
  geom_line(aes(y = V1, color = "Strain 1")) +
  geom_line(aes(y = V2, color = "Strain 2")) +
  labs(y = "Population unit", x = "Time unit", color = "Strain") +
  theme_minimal()




########### Try different settings ##########

# Parameter sets for multiple simulations
param_sets <- expand.grid(r1 = c(0.1), r2 = c(0.09, 0.1, 0.12, 0.15, 0.2), K = 100,
                          alpha12 = 0, alpha21 = 0, mu = 0)

# Initial state and time sequence
state <- c(V1 = 10, V2 = 10)
times <- seq(0, 100, by = 1)

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

# Transform data for plotting
long_data <- pivot_longer(combined_results, cols = c("V1", "V2"), names_to = "Strain", values_to = "Population")

# Plotting
values = alpha(c(hue_pal()(3)[1], hue_pal()(3)[3]),0.7)
# Define shapes for different r2 values
# shapes <- c("0.05" = 16, "0.1" = 17, "0.15" = 18, "0.2" = 19)  # Using different shape numbers
filtered_points <- long_data %>%
  filter(Strain == "V2", time %in% seq(0, 100, by = 5))

# Plotting
ggplot(long_data, aes(x = time, y = Population, color = Strain)) +
  geom_line(data = filter(long_data, Strain == "V1"), linewidth = 1.2) +
  geom_line(data = filter(long_data, Strain == "V2"), aes(group = factor(r2)), linewidth = 0.5) +
  geom_point(data = filtered_points, aes(shape = factor(r2)), size = 1) +
  scale_color_manual(values = values, labels = c("Strain 1", "Strain 2")) +
  # scale_shape_manual(values = shapes, name = "r2", labels = c("0.05", "0.1", "0.15", "0.2")) +
  labs(y = "Population unit", x = "Time unit", color = "Strain") +
  theme_bw() +
  theme(legend.position = "top")


# Calculate e^(V1 - V2) and prepare the data
long_data_transformed <- combined_results %>%
  mutate(ExpDiff = V1/V2) %>%
  select(time, r2, ExpDiff)

values <- c('#4D7ACC', '#5A8BD9', '#619CFF', '#7CB2FF', '#99C9FF')
values <- c('#2A41AF', '#4D7AAF', '#619CAF', '#85BDAF', '#A8D4AF')

# Plotting
ggplot(long_data_transformed, aes(x = time, y = ExpDiff, color = factor(r2))) +
  geom_line() +
  scale_color_manual(values = values, name = "r2") +
  labs(y = expression(beta[1]/beta[2]), x = "Time unit", color = "r2 value") +
  theme_bw() +
  theme(legend.position = "top")

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
