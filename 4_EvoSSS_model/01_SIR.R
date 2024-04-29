library(rstan)
library(ggplot2)
realData_all <- read.csv("Covid19CasesWH.csv", row.names = 1)
realData <- realData_all[-c(1:24), ] 
observed_cases = realData$CaseNum

# Define data

N = sum(observed_cases)
I0 = sum(realData_all$CaseNum[22:24])

sir_data <- list(
  N = N,
  I0 = I0,
  S0 = N - I0,
  T = length(observed_cases),
  cases = observed_cases  
)

# Fit the model
fit <- stan(file = 'SIR.stan', data = sir_data, 
            iter = 15000, chains = 4, warmup = 10000,
            verbose = TRUE)

save(fit, file = 'SIR.rdata')
print(fit)

# Summary of parameters
print(summary(fit)$summary)

# Plot diagnostics
stan_diag(fit, chain = 0)
# stan_par(fit, par = c('beta'))
# stan_par(fit, par = c('gamma'))

# combine correlation and hist
stan_hist(fit)
# Generate trace plots for beta and gamma
traceplot(fit, pars = c("beta", "gamma"))






