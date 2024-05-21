rm(list = ls())
library(rstan)
library(ggplot2)

stanfit = function(CaseNum, ratio){
  observed_cases = CaseNum[-c(1:24)]
  
  # Define data
  N = sum(CaseNum)
  I10 = round(sum(CaseNum[22:24])*ratio)
  I20 = round(sum(CaseNum[22:24])*(1-ratio))
  
  stan_data <- list(
    N = N,
    I10 = I10,
    I20 = I20,
    S0 = N - I10 - I20,
    T = length(observed_cases),
    cases = observed_cases  
  )

  # Fit the model
  fit <- stan(file = 'evoSIR.stan', data = stan_data, 
              iter = 15000, chains = 4, warmup = 10000,
              verbose = TRUE)
  
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
save(fitlist, file = 'evoSIR.rdata')

# Generate trace plots for beta and gamma
pdf(file = 'Output/acrosshost_sensitivity.pdf', 
    width = 9, height = 2)

for (i in 1:9) {
  print(traceplot(fitlist[[i]]))
}

dev.off()



