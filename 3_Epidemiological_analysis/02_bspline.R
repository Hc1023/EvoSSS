rm(list = ls())
library(dplyr)
library(ggplot2)
library(scales)
library(mgcv)
library(rstan)
df = read.csv('Covid19CasesGISAID.csv')
df$Var1 = as.Date(df$Var1)
df = df[df$Var1 < as.Date('2021-12-01'),]

bspline = function(df){
  ndays = nrow(df)
  time <- 1:ndays
  cases <- df$Freq
  # Fit a GAM with P-splines
  gam_model <- gam(cases ~ s(time, bs = "ps"), 
                   family = poisson(link = "log"))
  
  
  # Extracting model matrix
  X <- model.matrix(gam_model)
  
  # Preparing data for Stan
  stan_data <- list(
    N = length(cases),
    K = ncol(X),
    y = cases,
    X = X,
    mu = gam_model$coefficients
  )
  
  
  stan_code <- sprintf("
data {
  int<lower=1> N;
  int<lower=1> K;
  int y[N];
  matrix[N, K] X;
  vector[K] mu;
}
parameters {
  vector[K] beta;
}
model {
  for (k in 1:K) {
    beta[k] ~ normal(mu[k], 10);
  }
  y ~ poisson_log(X * beta);
}
")
  # Run the model
  # The default and preferred algorithm is "NUTS", which is the No-U-Turn sampler variant of Hamiltonian Monte Carlo (Hoffman and Gelman 2011, Betancourt 2017). 
  fit <- stan(model_code = stan_code, data = stan_data, 
              iter = 10000, chains = 3, verbose = T,
              control = list(max_treedepth = 20))
  
  return(fit)
}


df1 = df[df$Mutations == 'Lineage A',]
fit1 = bspline(df1)
df2 = df[df$Mutations == 'Lineage B',]
fit2 = bspline(df2)

save(fit1, fit2, file = 'grBspline.Rdata')

load('grBspline.Rdata')

evaluate_performance = function(fit){
  stan_est = extract(fit) 
  stan_est = as.matrix(fit)
  nchain = 3
  npar = 10
  performance = data.frame()
  for (i in 1:npar) {
    beta_est = stan_est[,i]
    beta_est_mat = matrix(beta_est, ncol=nchain)
    rhat = Rhat(beta_est_mat)
    essb = ess_bulk(beta_est_mat)
    esst = ess_tail(beta_est_mat)
    performance = rbind(performance, c(rhat, essb, esst))
  }
  colnames(performance) = c('rhat', 'essb', 'esst')
  return(performance)
}

performance1 = evaluate_performance(fit1)
performance2 = evaluate_performance(fit2)

traceplot(fit, pars = c("beta[2]"))


p = stan_diag(fit)
pdf(paste0("Output/bspline_diag.pdf"), width = 4, height = 5)
print(p)
dev.off()

p = traceplot(fit)
pdf(paste0("Output/bspline_performance.pdf"), width = 10, height = 4)
print(p)
dev.off()
