rm(list = ls())
library(deSolve) 
library(rootSolve) 
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
  
  dV1dt <- r1 * V1 * (1 - (V1 + alpha21 * V2) / K) - mu * V1
  dV2dt <- r2 * V2 * (1 - (V2 + alpha12 * V1) / K) + mu * V1
  
  return(c(dV1dt, dV2dt))
}

state <- c(V1 = 1, V2 = 1)

# Time sequence for the simulation
times <- seq(0, 72, by = 0.1)
params1 = c(r1 = 0.4, 
           r2 = 0.45, 
           K = 100000, 
           alpha12 = 1, 
           alpha21 = 0, 
           mu = 0)
params2 = c(r1 = 0.6, 
           r2 = 0.65, 
           K = 100000, 
           alpha12 = 1, 
           alpha21 = 0, 
           mu = 0)
params3 = c(r1 = 0.4, 
            r2 = 0.45, 
            K = 100000, 
            alpha12 = 0, 
            alpha21 = 0, 
            mu = 0)
paramlist = list(params1, params2, params3)

equilibrium_fun = function(params){
  # Initial state
  state <- c(V1 = 100000, V2 = 100000)
  times <- seq(0, 72, by = 1)
  viral_model_equilibrium <- function(x, params) {
    viral_model(state = x, parameters = params)
  }
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
}
for (i in 1:3) {
  equilibrium_fun(paramlist[[i]])
}

# V1           V2 
# 1.000000e+05 1.906811e-01 
# [,1]          [,2]
# [1,] -4.000000e-01  0.000000e+00
# [2,] -8.580649e-07 -1.716126e-06
# [1] -4.000000e-01 -1.716126e-06
# The system is stable at the equilibrium.
# V1           V2 
# 1.000000e+05 9.545885e-02 
# [,1]          [,2]
# [1,] -6.000000e-01  0.000000e+00
# [2,] -6.204825e-07 -1.240963e-06
# [1] -6.000000e-01 -1.240963e-06
# The system is stable at the equilibrium.
# V1    V2 
# 1e+05 1e+05 
# [,1]  [,2]
# [1,] -0.4  0.00
# [2,]  0.0 -0.45
# [1] -0.40 -0.45
# The system is stable at the equilibrium.
