rm(list = ls())
library(rstan)
library(ggplot2)
library(scales)
library(ggnewscale)
library(tidyverse)
library(dplyr)
library(RColorBrewer)

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

update_fun_stochastic = function(pars, states_old, N){
  
  beta1 = pars[1]
  beta2 = pars[2]
  gamma = pars[3]
  
  mat = data.frame()
  for (h in 1:nrow(states_old)) {
    S <- states_old[h,1]
    I1 <- states_old[h,2]
    I2 <- states_old[h,3]
    
    pS_vec = c(max(0,beta1*I1/N), max(0,beta2*I2/N), 
               max(0,1-beta1*I1/N-beta2*I2/N))
    sample_S <- rmultinom(1, size = S, prob = pS_vec)
    pI_vec <- c(gamma, 1-gamma)
    sample_I1 <- rmultinom(1, size = I1, prob = pI_vec)
    sample_I2 <- rmultinom(1, size = I2, prob = pI_vec)
    
    ## new values
    S_new <- sample_S[3]
    I1_new <- sample_I1[2] + sample_S[1] 
    I2_new <- sample_I2[2] + sample_S[2] 
    
    Onset1 <- sample_S[1]
    Onset2 <- sample_S[2]
    mat[h,1:5] = c(S_new, I1_new, I2_new, Onset1, Onset2)
  }
  
  return(mat)
}


simu <- function(seed_mat_I1, seed_mat_I2, N, poolday,
                 pars = c(0.379, 0.398, 0.157)) {
  
  # Initial conditions
  I1_old = seed_mat_I1[1,]
  I2_old = seed_mat_I2[1,]
  S_old = N - I1_old - I2_old
  states_old = cbind(S_old, I1_old, I2_old)
  
  # Simulate the dynamics over 
  ndays = 200
  
  Onsets_mat <- matrix(0, ndays, 2)
  
  for (t in 2:ndays) {
    states_old = update_fun_stochastic(pars = pars, states_old = states_old, N = N)
    if(t <= nrow(seed_mat_I1)){
      states_old[t,2:3] = states_old[t,2:3] + 
        c(seed_mat_I1[t,t], seed_mat_I2[t,t])
    }
    Onsets_mat[t,] = c(sum(states_old[,4]), sum(states_old[,5]))
  }
  
  return(Onsets_mat)
}


poolday = 30
# The initial cycle - epidemic outbreak
seed_vec = matrix(0,2,2)
seed_vec[,1] = round(34*c(0.3,0.7))
seed_mat_I1 = diag(seed_vec[1,])
seed_mat_I2 = diag(seed_vec[2,])
N = rep(32583, 2)
Onsets_mat_list = list()
Onsets_mat = simu(seed_mat_I1, seed_mat_I2, N, poolday, pars = c(0.38, 0.38, 0.157))
sum(Onsets_mat[,1])/(sum(Onsets_mat[,1])+sum(Onsets_mat[,2]))

Onsets_mat_list[[1]] = Onsets_mat

# document_scale = max(Onsets_mat[,2])/max(df$Freq[df$Var1<as.Date('2020-02-01')])
observed_matrix = observed_matrix*28

poolday = 30
nday = 100
expected_matrix = observed_matrix[1:nday,]
expected_matrix[expected_matrix < 0] = 0
expected_matrix = round(expected_matrix)
expected_matrix[(poolday*2+1):nday,] = 0
fexpect = data.frame(x = rep(1:nday,2), 
                     y = c(expected_matrix[1:nday,1],
                           expected_matrix[1:nday,2]),
                     group = factor(rep(c('A','B'), 
                                        each = nday),
                                    levels = c('A','B')))

fonset = data.frame(x = rep(1:nday,2), 
                    y = c(Onsets_mat[1:nday,1],
                          Onsets_mat[1:nday,2]),
                    group = factor(rep(c('A','B'), each = nday),
                                   levels = c('A','B')))


ggplot() +
  geom_point(data = fexpect, 
             aes(x = x, y = y, group = group, color = group)) +
  geom_line(data = fonset,
            aes(x = x, y = y, group = group, color = group))

