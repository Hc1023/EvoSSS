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


traceplot(fit, pars = c("beta[2]"))

getplot = function(df, fit){
  ndays = nrow(df)
  time = 1:ndays
  cases =  df$Freq
  # Fit a GAM with P-splines
  gam_model <- gam(cases ~ s(time, bs = "ps"), 
                   family = poisson(link = "log")
                   )
  # Extracting model matrix
  X <- model.matrix(gam_model)
  posterior_samples <- extract(fit)
  simu_Onset <- matrix(NA, nrow = 15000, ncol = ndays)
  for (i in 1:15000) {
    pars = posterior_samples$beta[i,]
    simu_Onset[i, ] <- exp(X %*% pars)  
  }
  ci_lower <- apply(simu_Onset, 2, quantile, probs = 0.025, na.rm = T)
  ci_upper <- apply(simu_Onset, 2, quantile, probs = 0.975, na.rm = T)
  start_date <- as.Date("2019-12-24") 
  date_vector <- seq.Date(start_date, by = "day", length.out = ndays)
  
  plot_data <- data.frame(
    Day = 1:ndays,
    date_vector = date_vector,
    Observed = cases,
    Fitted = colMeans(simu_Onset),
    LowerCI = ci_lower,
    UpperCI = ci_upper
  )
  return(plot_data)
}



posterior_samples <- extract(fit2)
pars<- apply(posterior_samples$beta, 2, mean)
fitted_values <- exp(X %*% pars)  
original_data <- data.frame(time = 1:ndays, cases = cases)
original_data$fitted_values <- fitted_values
p <- ggplot(original_data, aes(x = time)) +
  geom_line(aes(y = cases), colour = "blue", size = 1, alpha = 0.7, linetype = "dotted") +
  geom_line(aes(y = fitted_values), colour = "red", size = 1) +
  labs(title = "Comparison of Observed and Fitted Cases",
       x = "Time",
       y = "Number of Cases",
       caption = "Blue dotted line: Observed Data\nRed line: Fitted Values") +
  theme_minimal()

# Display the plot
print(p)


