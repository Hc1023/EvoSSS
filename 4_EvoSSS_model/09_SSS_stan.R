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
df = df[df$Var1 < as.Date('2021-12-01'),]
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
document_scale = max(Onsets_mat[,2])/max(df$Freq[df$Var1<as.Date('2020-02-01')])
observed_matrix = observed_matrix*28
n=3
for (n in 6:20) {
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
  
  # first fit poolday_fit
  
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
  fit <- stan(file = 'evoSSS_poolday.stan', data = stan_data, 
              iter = 15000, chains = 1, warmup = 10000,
              verbose = TRUE)
  fit
  posterior = rstan::extract(fit)
  contact = mean(posterior$contact)
  N = seed_vec * contact + 1
  # poolday*n + 1:100
  contact_vec[n] = contact
  Onsets_mat = simu(seed_mat_I1, seed_mat_I2, N, poolday, pars)
  Onsets_mat_list[[n+1]] = Onsets_mat
}



contact_vec_init = c(10000,100,100,150,280,260,260,260,260,260,260,
                     260,260,260,150,160,200,400,300,200,200,260)

contact_vec = c()
contact_vec = c(contact_vec, contact)


dfobserve = data.frame(
  x = rep(1:700,2),
  y = c(observed_matrix[,1], observed_matrix[,2]),
  group = c(rep('A', 700),rep('B', 700))
)

{
  Onsets_mat = Onsets_mat_list[[1]]
  for (n in 1:length(contact_vec)) {
    Onset1 = Onsets_mat[1:poolday+poolday,1]
    Onset2 = Onsets_mat[1:poolday+poolday,2]
    
    seed_vec =  (Onset1 + Onset2) %*% Mobility_matrix %>% as.numeric()
    p = Onset1/(Onset1 + Onset2)
    seed_matrix = rbind(seed_vec * p, seed_vec * (1-p))
    
    seed_mat_I1 = diag(seed_matrix[1,])
    seed_mat_I2 = diag(seed_matrix[2,])
    N = (seed_vec+1)*contact_vec[n]
    Onsets_mat = simu(seed_mat_I1, seed_mat_I2, N, poolday, pars)
    Onsets_mat_list[[n+1]] = Onsets_mat
  }
  
  dfplot_simu = data.frame()
  for (i in 1:(length(contact_vec)+1)) {
    Onsets_mat = Onsets_mat_list[[i]]
    n = i-1
    
    dfplot_simu1 = data.frame(x = rep(1:nrow(Onsets_mat)+poolday*n,2),
                              y = c(Onsets_mat[,1],Onsets_mat[,2]),
                              group = rep(paste0(c('A','B'),n), 
                                          each = nrow(Onsets_mat)),
                              color = rep(c('A','B'), 
                                          each = nrow(Onsets_mat)))
    dfplot_simu = rbind(dfplot_simu, dfplot_simu1)
  }
  
  df2 = dfplot_simu %>% group_by(x, color) %>% 
    summarise(y = sum(y)) %>%
    as.data.frame()
  values = c(hue_pal()(3)[1], hue_pal()(3)[3])
  
  
  ggplot() +
    geom_point(data = dfobserve, 
               aes(x = x, y = y, 
                   group = group, color = group)) +
    geom_line(data = dfplot_simu, 
              aes(x = x, y = y, 
                  group = group, color = color)) + 
    geom_line(data = df2, 
              aes(x = x, y = y, 
                  group = color, color = color)) + 
    scale_y_continuous(trans=scales::pseudo_log_trans(base=2),
                       breaks = c(0, 2^seq(0,20,2))) + 
    theme_bw() +
    scale_color_manual(name="Variant",
                       labels=c("A", "B"),
                       values = alpha(values, 0.6)) +
    scale_fill_manual(name="Variant",
                      labels=c("A", "B"),
                      values = alpha(values, 0.6))
  
}
