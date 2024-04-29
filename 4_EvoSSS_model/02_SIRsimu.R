rm(list = ls())
library(rstan)
library(ggplot2)
load(file = 'SIR.rdata')

update_fun = function(pars, states_old){
  
  S = states_old[1]
  I = states_old[2]
  R = states_old[3]

  S_new = S - beta*S*I/N
  I_new = I + beta*S*I/N - gamma*I
  R_new = R + gamma*I
  Onset = beta*S*I/N
 
  return(c(S_new, I_new, R_new, Onset))
}

update_fun_stochastic = function(pars, states_old){
  
  S <- states_old[1]
  I <- states_old[2]
  R <- states_old[3]

  beta = pars[1]
  gamma = pars[2]

  
  pS_vec = c(beta*I/N, 1-beta*I/N)
  sample_S <- rmultinom(1, size = S, prob = pS_vec)
  pI_vec <- c(gamma, 1-gamma)
  sample_I <- rmultinom(1, size = I, prob = pI_vec)
  
  ## new values
  S_new <- sample_S[2]
  I_new <- sample_I[2] + sample_S[1] 
  R_new <- R + sample_I[1]
  Onset <- sample_S[1]
  
  return(c(S_new, I_new, R_new, Onset))
}

simu <- function(pars, S0, I0, N, ndays, f = update_fun) {
  # Initial conditions
  states_old = c(S0, I0, N - S0 - I0, NA)
  
  # Simulate the dynamics over T days
  mycol <- c("time", "S", "I", "R", "Onset")
  states_mat <- matrix(0, ndays, length(mycol))
  states_mat[, 1] <- 1:ndays
  states_mat[1, 2:length(mycol)] <- states_old
  
  colnames(states_mat) <- mycol
  
  for (t in 1:(ndays-1)) {
    states_mat[t+1, -1] <- f(pars = pars, states_old = states_old)
    states_old <- states_mat[t+1, -1]
  }
  
  return(states_mat)
}


realData_all <- read.csv("Covid19CasesWH.csv", row.names = 1)
realData <- realData_all[-c(1:24), ] 
observed_cases = realData$CaseNum

# Define data

N = sum(observed_cases)
I0 = sum(realData_all$CaseNum[22:24])

# Setting initial conditions
N <- N
I0 <- I0
S0 <- N - I0
ndays <- length(observed_cases)

posterior_samples <- extract(fit)
simu_Onset <- matrix(NA, nrow = 15000, ncol = ndays)
# Simulate the epidemic for each set of sampled parameters
for (i in 1:15000) {
  pars = c(posterior_samples$beta[i], posterior_samples$gamma[i])
  result <- simu(pars, S0, I0, N, ndays,
                 f = update_fun_stochastic)
  simu_Onset[i, ] <- result[,'Onset']
}

# Calculate the 2.5th and 97.5th percentiles for the confidence interval
ci_lower <- apply(simu_Onset, 2, quantile, probs = 0.025, na.rm = T)
ci_upper <- apply(simu_Onset, 2, quantile, probs = 0.975, na.rm = T)

start_date <- as.Date("2020-01-01") 
date_vector <- seq.Date(start_date, by = "day", length.out = ndays)

# Create data frame for plotting
plot_data <- data.frame(
  Day = 1:ndays,
  datep_vector = date_vector,
  Observed = observed_cases,
  Fitted = colMeans(simu_Onset),
  LowerCI = ci_lower,
  UpperCI = ci_upper
)


# Plot using ggplot2
ggplot(plot_data, aes(x = date_vector)) +
  geom_ribbon(aes(ymin = LowerCI, ymax = UpperCI), fill = '#ead5bf', alpha = 0.6) +  # Confidence interval
  geom_point(aes(y = Observed), color = '#241508', size = 1, alpha = 0.5, shape = 16) +  # Observed data
  geom_line(aes(y = Fitted), color = '#80553c', size = 0.8, alpha = 0.7) +  # Fitted line
  labs(x = "Date (2020)", y = "Cases") +
  scale_x_date(date_labels = "%b-%d", date_breaks = "1 month") +  
  theme_bw()





