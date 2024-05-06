rm(list = ls())
library(rstan)
library(ggplot2)

load(file = 'sir.rdata')
posterior_samples <- extract(fit)
beta = mean(posterior_samples$beta)
gamma = mean(posterior_samples$gamma)
stanfit = function(CaseNum, I0){
  realData_all <- read.csv("Covid19CasesWH.csv", row.names = 1)
  CaseNum = realData_all$CaseNum
  
  
  I0 = 0.5
  observed_cases = CaseNum[-c(1:24)]
  re = round(log(observed_cases[1]/I0/2)/(beta-gamma))
  ndays = length(observed_cases) + re
  # Define data
  N = sum(CaseNum)
 
  stan_data <- list(
    N = N,
    S0 = N - 1,
    T =  ndays,
    cases = c(rep(0, re), observed_cases),
    beta = beta,
    gamma = gamma,
    I0 = I0
  )
  
  # Fit the model
  fit <- stan(file = 'evoSIR_intro.stan', data = stan_data, 
              iter = 20000, chains = 1, warmup = 15000,
              verbose = TRUE)
  print(fit)
  traceplot(fit)
  return(fit)
}


realData_all <- read.csv("Covid19CasesWH.csv", row.names = 1)
CaseNum = realData_all$CaseNum

## sensitivity analysis
fitlist = list()
for (i in c(1:9)) {
  ratio_initial  = i/10
  fitlist[[i]] = stanfit(CaseNum, ratio_initial)
}

# print(fit)
save(fitlist, file = 'SIR.rdata')