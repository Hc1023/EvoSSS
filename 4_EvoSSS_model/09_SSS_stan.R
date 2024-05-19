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
  ndays = dim*2
  
  I1_old = seed_mat_I1[1,]
  I2_old = seed_mat_I2[1,]
  S_old = N - I1_old - I2_old
  states_old = cbind(S_old, I1_old, I2_old)
  
  # Simulate the dynamics over ndays
  Onsets_mat <- matrix(0, poolday, 2)
  
  colnames(Onsets_mat) <- mycol
  for (t in 1:(ndays-1)) {
    states_new = update_fun(pars = pars, states_old = states_old, N = N)
    states_new[,2:3] = states_new[,2:3] + 
      cbind(seed_mat_I1[t+1,], seed_mat_I2[t+1,])
    Onsets_mat[t+1,-1] = c(sum(states_new[,4]), sum(states_new[,5]))
    states_old = states_new
  }
  
  return(Onsets_mat)
}



dim = 100
pars = c(0.379, 0.398, 0.157)
# The initial cycle - epidemic outbreak
seed_vec = rbind(rep(0,dim),rep(0,dim))
seed_vec[,1] = round(34*c(0.4,0.6))
seed_mat_I1 = rbind(diag(seed_vec[1,]), matrix(0, dim, dim))
seed_mat_I2 = rbind(diag(seed_vec[2,]), matrix(0, dim, dim))
N = rep(32583, dim)
Onsets_mat_list = list()
Onsets_mat = simu(seed_mat_I1, seed_mat_I2, N, dim, pars)
Onsets_mat_list[[1]] = Onsets_mat

n=1
Onset1 = Onsets_mat[poolday*n + 1:poolday, 2]
Onset2 = Onsets_mat[poolday*n + 1:poolday, 3]
mobility = rep(0.01,30) # Mobility: Control force
Mobility_matrix = diag(mobility)
seed_vec =  (Onset1 + Onset2) %*% Mobility_matrix
seed_vec = seed_vec[1,]
p = Onset1/(Onset1 + Onset2)
seed_matrix = rbind(seed_vec * p, seed_vec * (1-p))
seed_mat_I1 = diag(seed_matrix[1,])
seed_mat_I2 = diag(seed_matrix[2,])

observed1 = df[df$Mutations == 'Lineage A' & df$Var1 >= as.Date('2020-01-01'),]
observed2 = df[df$Mutations == 'Lineage B' & df$Var1 >= as.Date('2020-01-01'),]
observed_matrix = cbind(observed1$Freq, observed2$Freq)
document_scale = max(simu_matrix[,2])/max(df$Freq[df$Var1<as.Date('2020-02-01')])
observed_matrix = round(observed_matrix*document_scale)
expected_matrix = observed_matrix[poolday+1:poolday,] - Onsets_mat[poolday+1:poolday,2:3]
expected_matrix[expected_matrix[,1]<0,1] = 0
expected_matrix[expected_matrix[,2]<0,2] = 0
expected_matrix = round(expected_matrix)
contact_vec_init = c(10000,100,100,150,280,260,260,260,260,260,260,
                     260,260,260,150,160,200,400,300,200,200,260)
contact_init = contact_vec_init[1]

stan_data <- list(
  poolday = poolday,
  contact_init = contact_init,
  expected_matrix = expected_matrix, 
  pars = pars,
  seed_mat_I1 = seed_mat_I1,
  seed_mat_I2 = seed_mat_I2,
  seed_vec = seed_vec
)

# Fit the model
fit <- stan(file = 'evoSSS.stan', data = stan_data, 
            iter = 15000, chains = 1, warmup = 10000,
            verbose = TRUE)
# 
n = 1
poolday = 30
mobility = rep(0, dim)
mobility[seq(1,30,1)] = 0.01  # Mobility: Control force

contact_vec_init = c(10000,100,100,150,280,260,260,260,260,260,260,
                260,260,260,150,160,200,400,300,200,200,260)
contact_vec = contact_vec_init[1]
Onsets_mat = Onsets_mat_list[[1]]
cycle = length(contact_vec)
simu_days = dim * 2 + poolday *cycle
observed_days = poolday * (cycle+1)


simu_matrix = matrix(0, simu_days, 2)
simu_matrix[1:nrow(Onsets_mat_list[[1]]),] = Onsets_mat[,c(2,3)] 



{
  for (n in 1:length(contact_vec)) {
    # seeding_time = 1:dim + spacing*(n-1)
    Onset1 = Onsets_mat[1:dim+poolday,2]
    Onset2 = Onsets_mat[1:dim+poolday,3]
    
    seed_vec =  (Onset1 + Onset2) %*% Mobility_matrix
    seed_vec = seed_vec[1,]
    p = Onset1/(Onset1 + Onset2)
    seed_matrix = rbind(seed_vec * p, seed_vec * (1-p))
    # seed_matrix[,1] = 0
    
    seed_mat_I1 = rbind(diag(seed_matrix[1,]), matrix(0, dim, dim))
    seed_mat_I2 = rbind(diag(seed_matrix[2,]), matrix(0, dim, dim))
    N = (seed_vec+1)*contact_vec[n]
    Onsets_mat = simu(seed_mat_I1, seed_mat_I2, N, dim, pars)
    Onsets_mat[,1] = Onsets_mat[,1] + poolday*n
    Onsets_mat_list[[n+1]] = Onsets_mat
    simu_matrix[1:(2*dim) + poolday * n, ] = simu_matrix[1:(2*dim) + poolday * n, ] + Onsets_mat[,2:3]
  }
  
  simu_matrix = simu_matrix[1:700,]
  dfobserve = data.frame(
    x = rep(1:700,2),
    y = c(observed_matrix[,1], observed_matrix[,2]),
    group = c(rep('A', 700),rep('B', 700))
  )
  dfsimu = data.frame(
    x = rep(1:700,2),
    y = c(simu_matrix[,1], simu_matrix[,2]),
    group = c(rep('A', 700),rep('B', 700))
  )
  
  ggplot() +
    geom_point(data = dfobserve, 
               aes(x = x, y = y, 
                   group = group, color = group)) +
    geom_line(data = dfsimu, 
              aes(x = x, y = y, 
                  group = group, color = group)) + 
    scale_y_continuous(trans=scales::pseudo_log_trans(base=2),
                       breaks = c(0, 2^seq(0,14,2))) + 
    theme_bw()
  
  str(simu_matrix)
  str(observed_matrix)
  log_likelihood <- 0  # Initialize log-likelihood
  
  # Loop over each day
  for (i in 700 + 1:700) {
    # Calculate the log-likelihood for each variable in the matrix
    log_likelihood <- log_likelihood + 
      dpois(dfobserve[i, 2], lambda = dfsimu[i, 2]+1, log = TRUE)
  }
  log_likelihood 
  
  for (i in 1:observed_days) {
    
  }
  
  dfplot_simu = data.frame()
  for (i in 1:(length(contact_vec)+1)) {
    Onsets_mat = Onsets_mat_list[[i]]
    Onsets_mat = Onsets_mat[-1,]
    n = i-1
    
    dfplot_simu1 = data.frame(x = rep(Onsets_mat[,1],2),
                              y = c(Onsets_mat[,2],Onsets_mat[,3]),
                              group = rep(paste0(c('A','B'),n), each = nrow(Onsets_mat)),
                              color = rep(c('A','B'), each = nrow(Onsets_mat)))
    dfplot_simu = rbind(dfplot_simu, dfplot_simu1)
  }
  
  
  df1 = dfplot_simu
  df1$x = as.Date('2019-12-31') + df1$x
  
  
  document_scale = max(df1[df1$group == 'B0',2])/max(df$Freq[df$Var1<as.Date('2020-02-01')])
  
  df1[,2] = df1[,2]/document_scale
  
  df2 = df1 %>% group_by(x, color) %>% 
    summarise(y = sum(y)) %>%
    as.data.frame()
  values = c(hue_pal()(3)[1], hue_pal()(3)[3])
  
  ggplot() + 
    geom_point(data = df, 
               aes(x = Var1, y = Freq, colour = Mutations),
               size = 0.5) +
    scale_color_manual(name="Variant",
                       labels=c("A", "B"),
                       values = alpha(values, 0.6)) +
    scale_y_continuous(trans=scales::pseudo_log_trans(base=2),
                       breaks = c(0, 2^seq(0,14,2))) +
    new_scale_color() +
    new_scale_fill() +
    geom_line(data = df1, 
              aes(x, y, color = color, 
                  group = group), linewidth = 1) + 
    geom_line(data = df2, 
              aes(x, y, color = color, 
                  group = color), linewidth = 1) + 
    scale_color_manual(name="Variant",
                       labels=c("A", "B"),
                       values = alpha(values, 0.6)) +
    scale_fill_manual(name="Variant",
                      labels=c("A", "B"),
                      values = alpha(values, 0.6)) +
    scale_x_date(breaks = seq(as.Date('2020-01-01'), as.Date('2022-11-01'), by="6 months"),
                 minor_breaks = seq(as.Date('2019-12-01'), as.Date('2022-11-01'), by ='1 month'),
                 date_labels = "%y-%b") +
    xlab('Date (2019-2021)') +
    ylab('') + theme_bw() 
  
}
