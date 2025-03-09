rm(list = ls())
library(rstan)
library(ggplot2)

stanfit = function(CaseNum, ratio, chains_num){
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
  fit <- stan(file = 'F3D_evoSIR.stan', data = stan_data, 
              iter = 3000, chains = chains_num, 
              warmup = 2000, verbose = TRUE)
  
  return(fit)
}

# fit = stanfit(CaseNum,  0.3, chains_num = 1)

realData_all <- read.csv("evoSIR/F3C_Covid19CasesWH.csv", row.names = 1)
CaseNum = realData_all$CaseNum

if(F){
  ## sensitivity analysis
  fitlist = list()
  for (i in c(1:9)) {
    ratio_initial  = i/10
    fitlist[[i]] = stanfit(CaseNum, ratio_initial, 
                           chains_num = 4)
  }
  # print(fit)
  save(fitlist, file = 'evoSIR.rdata')
}
load('evoSIR/F3D_evoSIR.rdata')
# Generate trace plots for beta and gamma
pdf(file = 'Output/S7C_acrosshost_sensitivity.pdf', 
    width = 9, height = 2)

for (i in 1:9) {
  print(traceplot(fitlist[[i]]))
}

dev.off()



