rm(list = ls())
library(rstan)
library(ggplot2)
library(scales)

realData_all <- read.csv("Covid19CasesWH.csv", row.names = 1)
CaseNum = realData_all$CaseNum

observed_cases = realData$CaseNum

retro = 40
observed_cases = CaseNum[-c(1:24)]

# Define data
N = sum(CaseNum)

stan_data <- list(
  N = N,
  S0 = N - 1,
  T =  length(observed_cases) + retro,
  cases = c(rep(0, retro), observed_cases)
)

# Define data
N = sum(CaseNum)
S0 <- N - I0
ndays <- length(observed_cases)

update_fun = function(pars, states_old){
  
  S <- states_old[1]
  I1 <- states_old[2]
  I2 <- states_old[3]
  
  beta = pars[3]
  gamma = pars[4]
  
  S_new = S - beta*S*I1/N - beta*S*I2/N
  I1_new = I1 + beta*S*I1/N - gamma*I1
  I2_new = I2 + beta*S*I2/N - gamma*I2
  Onset1 = beta*S*I1/N
  Onset2 = beta*S*I2/N
  
  return(c(S_new, I1_new, I2_new, Onset1, Onset2))
}


simu <- function(pars = c(10,5,0.3,0.15), 
                 states_old = c(S0,0,0), 
                 ndays = 100, f = update_fun) {
  # Initial conditions
  states_old = c(states_old, NA, NA)
  
  # Simulate the dynamics over ndays
  mycol <- c("time", "S", "I1", "I2", "Onset1", "Onset2")
  states_mat <- matrix(0, ndays, length(mycol))
  states_mat[, 1] <- 1:ndays
  states_mat[1, 2:length(mycol)] <- states_old
  
  colnames(states_mat) <- mycol
  t1_intro = pars[1]
  t2_intro = pars[2]
  
  for (t in 1:(ndays-1)) {
    if (t >= t1_intro && states_old['I1'] == 0) {
      states_old['I1'] = 1
    }
    if (t >= t2_intro && states_old['I2'] == 0) {
      states_old['I2'] = 1 
    }
    states_mat[t+1, -1] <- f(pars = pars, states_old = states_old)
    states_old <- states_mat[t+1, -1]
  }
  
  return(states_mat)
}


states_mat = simu()
observed_cases
ndays = 100
data_observed = data.frame(Day =  1:ndays,
                           Observed = c(rep(0, ndays- length(observed_cases)), observed_cases))

plot_data <- data.frame(
  Day = rep(1:ndays,3),
  Fitted = c(states_mat[,'Onset1'], states_mat[,'Onset2'],
             states_mat[,'Onset1']+states_mat[,'Onset2']),
  V = rep(c('V1','V2','All'), each = ndays)
)

ggplot() +
  geom_point(data = data_observed, 
             aes(x = Day, y = Observed), 
             size = 0.5, alpha = 0.5) +  # Observed data
  geom_line(data = plot_data, 
            aes(x = Day, group = V, 
                y = Fitted, color = V))


