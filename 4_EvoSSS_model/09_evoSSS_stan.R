rm(list = ls())
library(rstan)
library(ggplot2)
library(scales)
library(ggplot2)
library(scales)
library(ggnewscale)
library(tidyverse)

df = read.csv('../3_Epidemiological_analysis/Covid19CasesGISAID.csv')

df$Var1 = as.Date(df$Var1)
# df = df[df$Var1 < as.Date('2021-12-01'),]
df = df[df$Mutations %in% c('Lineage A', 'Lineage B'),]

update_fun = function(pars, states_old, N){
  
  S <- states_old[,1]
  I1 <- states_old[,2]
  I2 <- states_old[,3]
  
  beta1 = pars[1]
  beta2 = pars[2]
  gamma = pars[3]
  
  S_new = S - beta1*S*I1/N - beta2*S*I2/N
  I1_new = I1 + beta1*S*I1/N - gamma*I1
  I2_new = I2 + beta2*S*I2/N - gamma*I2
  Onset1 = beta1*S*I1/N
  Onset2 = beta2*S*I2/N
  
  return(cbind(S_new, I1_new, I2_new, Onset1, Onset2))
}

simu <- function(seed_mat_I1, seed_mat_I2, N, poolday, pars) {
  # Initial conditions
  I1_old = seed_mat_I1[1,]
  I2_old = seed_mat_I2[1,]
  S_old = N - I1_old - I2_old
  states_old = cbind(S_old, I1_old, I2_old)
  
  # Simulate the dynamics over 
  ndays = 200
  
  Onsets_mat <- matrix(0, ndays, 2)
  
  for (t in 2:ndays) {
    states_old = update_fun(pars = pars, states_old = states_old, N = N)
    if(t <= nrow(seed_mat_I1)){
      states_old[t,2:3] = states_old[t,2:3] + 
        c(seed_mat_I1[t,t], seed_mat_I2[t,t])
    }
    Onsets_mat[t,] = c(sum(states_old[,4]), sum(states_old[,5]))
  }
  
  return(Onsets_mat)
}

pars = c(0.379, 0.398, 0.157)
poolday = 30
# The initial cycle - epidemic outbreak
seed_vec = matrix(0,2,2)
seed_vec[,1] = round(34*c(0.4,0.6))
seed_mat_I1 = diag(seed_vec[1,])
seed_mat_I2 = diag(seed_vec[2,])
N = rep(32583, 2)
Onsets_mat_list = list()
Onsets_mat = simu(seed_mat_I1, seed_mat_I2, N, poolday, pars)
Onsets_mat_list[[1]] = Onsets_mat

observed1 = df[df$Mutations == 'Lineage A' & df$Var1 >= as.Date('2020-01-01'),]
observed2 = df[df$Mutations == 'Lineage B' & df$Var1 >= as.Date('2020-01-01'),]
observed_matrix = cbind(observed1$Freq, observed2$Freq)
# document_scale = max(Onsets_mat[,2])/max(df$Freq[df$Var1<as.Date('2020-02-01')])
observed_matrix = observed_matrix*28

fitlist = list()

for (n in 1:24) {
  Onsets_mat = Onsets_mat_list[[n]]
  Onset1 = Onsets_mat[poolday + 1:poolday, 1]
  Onset2 = Onsets_mat[poolday + 1:poolday, 2]
  
  mobility = rep(0.01,30) # Mobility: Control force
  Mobility_matrix = diag(mobility)
  
  seed_vec =  (Onset1 + Onset2) %*% Mobility_matrix %>% as.numeric()
  p = Onset1/(Onset1 + Onset2)
  
  seed_matrix = rbind(seed_vec * p, seed_vec * (1-p))
  seed_mat_I1 = diag(seed_matrix[1,])
  seed_mat_I2 = diag(seed_matrix[2,])
  
  nday = 100
  
  expected_matrix = observed_matrix[poolday*n+1:nday,] - Onsets_mat[poolday+1:nday,]
  
  expected_matrix[expected_matrix[,1]<0,1] = 0
  expected_matrix[expected_matrix[,2]<0,2] = 0
  expected_matrix = round(expected_matrix)
  if(n == 1){
    expected_matrix[(nday-30):nday,] = 0
  }else{
    expected_matrix[(nday-50):nday,] = 0
  }
  expected_matrix[1:30,] = 0
  
  stan_data <- list(
    poolday = poolday,
    nday = nday,
    expected_matrix = expected_matrix, 
    pars = pars,
    seed_mat_I1 = seed_mat_I1,
    seed_mat_I2 = seed_mat_I2,
    seed_vec = seed_vec
  )
  # Fit the model
  fit <- stan(file = 'evoSSS.stan', data = stan_data, 
              iter = 15000, chains = 4, warmup = 10000,
              verbose = TRUE)
  
  posterior = rstan::extract(fit)
  contact = mean(posterior$contact)
  N = seed_vec * contact + 1

  Onsets_mat = simu(seed_mat_I1, seed_mat_I2, N, poolday, pars)
  
  fitlist[[n]] = fit
  Onsets_mat_list[[n+1]] = Onsets_mat
}

save(fitlist, file = 'evoSSS_chain.rdata')
posterior = rstan::extract(fitlist[[1]])
posterior$contact[1:5000]


# Generate trace plots 

for (i in 1:length(fitlist)) {
  wd = paste0('Output/traceplot/evoSSS_', 
              as.character(i), '.png')
  png(file = wd, 
      width = 2.5, height = 2, 
      unit = 'in',
      res = 600)
  p = traceplot(fitlist[[i]]) + 
    ggtitle(as.character(i)) + ylab('') +
    theme(legend.position = 'none') + 
    scale_x_continuous(breaks = seq(10000,15000,2500),
                       labels = seq(0,5000,2500))
  print(p)
  dev.off()
}

dfposterior = data.frame()
for (i in 1:length(fitlist)) {
  fit = fitlist[[i]]
  posterior_samples = rstan::extract(fit)
  
  dfposterior[i,1] = i
  dfposterior$m[i] = mean(posterior_samples$contact)
  dfposterior$q1[i] = quantile(posterior_samples$contact, 0.025)
  dfposterior$q2[i] = quantile(posterior_samples$contact, 0.975)
  dfposterior$text[i] = paste0(sprintf('%.3f', dfposterior$m[i]),' (',
                            sprintf('%.3f',dfposterior$q1[i]),' ~ ',
                            sprintf('%.3f',dfposterior$q2[i]),')')
  
}

write.csv(dfposterior, file = 'Output/evoSSS_parameters.csv',
          row.names = F)


dfposterior
ggplot(dfposterior, aes(V1, m)) +
  geom_line() +
  geom_point(shape = 19, alpha = 0.6) +
  geom_errorbar(
    aes(ymin = q1, ymax = q2),
    width = 0.5
  ) + theme_bw() +
  scale_y_continuous(trans='log10') +
  annotation_logticks(sides = "l", linewidth = 0.1, alpha = 0.5) +
  xlab('Cycle') + ylab('Susceptible ratio')
# Add ratio of A to B
