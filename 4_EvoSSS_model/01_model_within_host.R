rm(list = ls())
library(deSolve)   
library(ggplot2)   
# R1 = 2.65; R2 = 5.78
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
params <- c(r1 = 0.1, r2 = 0.11, K = 100, 
            alpha12 = 0.5, alpha21 = 0.5, 
            mu = 0)

# Initial state
state <- c(V1 = 10, V2 = 10)

# Time sequence for the simulation
times <- seq(0, 100, by = 0.1)

# Solve the model
out <- ode(y = state, times = times, func = viral_model, parms = params)
out_df <- as.data.frame(out)

# Plotting
ggplot(data = out_df, aes(x = time)) +
  geom_line(aes(y = V1, color = "Strain 1")) +
  geom_line(aes(y = V2, color = "Strain 2")) +
  labs(y = "Population unit", x = "Time unit", color = "Strain") +
  theme_minimal()





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
