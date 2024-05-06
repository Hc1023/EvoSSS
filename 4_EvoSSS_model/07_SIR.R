rm(list = ls())
library(rstan)
library(ggplot2)

if(F){
  # stan fitting
  realData_all <- read.csv("Covid19CasesWH.csv", row.names = 1)
  CaseNum = realData_all$CaseNum
  
  observed_cases = CaseNum[-c(1:24)]
  
  # Define data
  N = sum(CaseNum)
  I0 = sum(CaseNum[22:24])
  
  stan_data <- list(
    N = N,
    I0 = I0,
    S0 = N - I0,
    T = length(observed_cases),
    cases = observed_cases  
  )

  # Fit the model
  fit <- stan(file = 'SIR.stan', data = stan_data, 
              iter = 15000, chains = 4, warmup = 10000,
              verbose = TRUE)
  print(fit)
  traceplot(fit)
  save(fit, file = 'sir.rdata')
  return(fit)
}

update_fun_stochastic = function(pars, states_old){
  
  S <- states_old[1]
  I <- states_old[2]
  
  beta = pars[1]
  gamma = pars[2]
  
  pS_vec = c(beta*I/N, 1-beta*I/N)
  sample_S <- rmultinom(1, size = S, prob = pS_vec)
  pI_vec <- c(gamma, 1-gamma)
  sample_I <- rmultinom(1, size = I, prob = pI_vec)
  
  ## new values
  S_new <- sample_S[2]
  I_new <- sample_I[2] + sample_S[1] 
  
  Onset <- sample_S[1]
  return(c(S_new, I_new, Onset))
}

simu <- function(pars, states_old, ndays, f = update_fun_stochastic) {
  # Initial conditions
  states_old = c(states_old, NA)
  
  # Simulate the dynamics over ndays
  mycol <- c("time", "S", "I", "Onset")
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
S0 <- N - I0
ndays <- length(observed_cases)
load(file = 'sir.rdata')
posterior_samples <- extract(fit)

simu_Onset <- matrix(NA, nrow = 20000, ncol = ndays)
# Simulate the epidemic for each set of sampled parameters
for (i in 1:20000) {
  
  beta = posterior_samples$beta[i]
  gamma = posterior_samples$gamma[i]
  pars = c(beta = beta, gamma = gamma)
  states_old = c(S0,I0)
  result <- simu(pars, states_old, ndays,
                 f = update_fun_stochastic)
  simu_Onset[i, ] <- result[,'Onset']
}

start_date <- as.Date("2020-01-01") 
date_vector <- seq.Date(start_date, by = "day", length.out = ndays)

generate_plotdata = function(df){
  simu_Onset = df
  # Calculate the 2.5th and 97.5th percentiles for the confidence interval
  ci_lower <- apply(simu_Onset, 2, quantile, probs = 0.025, na.rm = T)
  ci_upper <- apply(simu_Onset, 2, quantile, probs = 0.975, na.rm = T)
  
  plot_data <- data.frame(
    Day = 1:ndays,
    date_vector = date_vector,
    Fitted = colMeans(simu_Onset),
    LowerCI = ci_lower,
    UpperCI = ci_upper
  )
  
  return(plot_data)
}

data_observed <- data.frame(
  Day = 1:ndays,
  date_vector = date_vector,
  Observed = observed_cases
)


plot_data = generate_plotdata(df = simu_Onset)

# Plot using ggplot2


p = ggplot() +
  geom_ribbon(data = plot_data, 
              aes(x = date_vector, ymin = LowerCI, ymax = UpperCI),
              fill = alpha("#aa85a6",0.3)) +  # Confidence interval
  geom_point(data = data_observed, 
             aes(x = date_vector, y = Observed), 
             size = 0.5, alpha = 0.5) +  # Observed data
  geom_line(data = plot_data, 
            aes(x = date_vector, y = Fitted),
            color = alpha("#aa85a6", 0.8)) +  
  labs(x = "Date (2020)", y = "Cases") +
  scale_x_date(date_labels = "%b-%d", date_breaks = "1 month") +  
  theme_bw() +
  theme(legend.background = element_blank(),
        legend.position = c(0.83,0.75))
p
pdf(paste0("Output/sirWH.pdf"), width = 2.5, height = 1.5)
print(p)
dev.off()

pdf(file = 'Output/sir_trace.pdf', 
    width = 6.3, height = 2)

traceplot(fit)

dev.off()




