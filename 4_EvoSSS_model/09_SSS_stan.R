rm(list = ls())
library(rstan)
library(ggplot2)

df = read.csv('../3_Epidemiological_analysis/Covid19CasesGISAID.csv')

df$Var1 = as.Date(df$Var1)
df = df[df$Var1 < as.Date('2021-12-01'),]
df = df[df$Mutations %in% c('Lineage A', 'Lineage B'),]


stan_data <- list(
  N = N,
  I10 = I10,
  I20 = I20,
  S0 = N - I10 - I20,
  T = length(observed_cases),
  cases = observed_cases  
)
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

