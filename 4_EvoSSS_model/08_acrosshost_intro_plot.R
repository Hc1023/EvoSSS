rm(list = ls())
library(rstan)
library(ggplot2)
library(scales)

realData_all <- read.csv("Covid19CasesWH.csv", row.names = 1)
CaseNum = realData_all$CaseNum
I0 = sum(CaseNum[22:24])

load(file = 'evoSIR_eqbeta.rdata')
posterior_samples <- extract(fit)

beta = posterior_samples$beta
gamma = posterior_samples$gamma
I10 = posterior_samples$I10
I20 = I0 - I10
hist(beta-gamma)


time_intro1 = log(I10)/(beta-gamma) # - log(ratio)
time_intro2 = log(I20)/(beta-gamma) # - log(ratio)

delta1 = time_intro2-time_intro1
hist(delta1)

load(file = 'evoSIR.rdata')
df2 = data.frame()
for (i in 1:9) {
  fit = fitlist[[i]]
  posterior_samples <- extract(fit)
  beta1 = posterior_samples$beta1
  beta2 = posterior_samples$beta2
  gamma = posterior_samples$gamma
  ratio = i/10
  I10 = round(sum(CaseNum[22:24])*ratio)
  I20 = round(sum(CaseNum[22:24])*(1-ratio))
  time_intro1 = log(I10)/(beta1-gamma) # - log(m)
  time_intro2 = log(I20)/(beta2-gamma) # - log(m)
  df = data.frame(x = i, y = time_intro2-time_intro1)
  df2 = rbind(df2, df)
}
library(ggpubr)
ggboxplot(df2, x = "x", y = "y")
