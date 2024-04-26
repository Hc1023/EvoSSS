library(deSolve)  # For solving differential equations
library(rootSolve)  # For finding equilibrium points
library(numDeriv)
# Define the model
viral_model <- function(state, parameters) {
  V1 <- state[1]
  V2 <- state[2]
  
  r1 <- parameters["r1"]
  r2 <- parameters["r2"]
  K <- parameters["K"]
  alpha12 <- parameters["alpha12"]
  alpha21 <- parameters["alpha21"]
  mu <- parameters["mu"]
  
  dV1dt <- r1 * V1 * (1 - (V1 + alpha12 * V2) / K) - mu * V1
  dV2dt <- r2 * V2 * (1 - (V2 + alpha21 * V1) / K) + mu * V1
  
  return(c(dV1dt, dV2dt))
}

# Parameters
params <- c(r1 = 0.1, r2 = 0.2, K = 100, alpha12 = 0.1, alpha21 = 0.0, mu = 0.1)

# Initial state
state <- c(V1 = 80, V2 = 80)

# Find equilibrium
equilibrium <- multiroot(f = viral_model, start = state, parameters = params)
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
