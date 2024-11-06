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

N = rep(32583, 2)

pc = c()
seeds = c(50,20,10,6,3)
nn = 1000
for (s in seeds) {
  print(s)
  for (i in 1:nn) {
    seed_vec = matrix(0,2,2)
    seed_vec[,1] = rmultinom(1, size = s, prob = c(0.3,0.7))
    seed_mat_I1 = diag(seed_vec[1,])
    seed_mat_I2 = diag(seed_vec[2,])
    Onsets_mat = simu(seed_mat_I1, seed_mat_I2, N, poolday, pars = c(0.38, 0.38, 0.157))
    p = sum(Onsets_mat[,1])/(sum(Onsets_mat[,1])+sum(Onsets_mat[,2]))
    pc = c(pc, p)
  }
  
}

dat1 = data.frame(pc = pc,
                  seed = rep(as.character(seeds), 
                             each = nn))
dat1$seed = factor(dat1$seed, levels = seeds)

pc = c()
seeds = c(50,20,10,6,3)
nn = 1000
for (s in seeds) {
  print(s)
  for (i in 1:nn) {
    seed_vec = matrix(0,2,2)
    seed_vec[,1] = rmultinom(1, size = s, prob = c(0.5,0.5))
    seed_mat_I1 = diag(seed_vec[1,])
    seed_mat_I2 = diag(seed_vec[2,])
    Onsets_mat = simu(seed_mat_I1, seed_mat_I2, N, poolday, pars = c(0.38, 0.38, 0.157))
    p = sum(Onsets_mat[,1])/(sum(Onsets_mat[,1])+sum(Onsets_mat[,2]))
    pc = c(pc, p)
  }
  
}

dat2 = data.frame(pc = pc,
                  seed = rep(as.character(seeds), 
                             each = nn))
dat2$seed = factor(dat2$seed, levels = seeds)

save(dat1, dat2, file = 'bottleneck.rdata')
values = c('#fbb365','#d3d667', 
           '#98d5b8','#87cad7','#b7b1ea')
values2 = c( '#cb8335','#a3a637',
   '#68a588','#579aa7','#8781ba')
library(ggpubr)

p1 = ggdensity(dat1, x = "pc", add = "mean",
               fill = 'seed', color = 'seed', rug = F) +
  ylab('Density') + xlab('') +
  geom_vline(xintercept = 0.3, linetype = "dashed") +
  scale_fill_manual(name = '',
                    values = values) +
  scale_color_manual(name = '',
                     values = values2) +
  theme(legend.position = 'right',
        legend.background = element_rect(color = NA, fill = NA),
        legend.key=element_blank(),
        legend.spacing.y = unit(0.05,'cm'),
        legend.key.size = unit(0.4,'cm'))
p1
p2 = ggdensity(dat2, x = "pc", add = "mean",
               fill = 'seed', color = 'seed', rug = F) +
  ylab('Density') + xlab('') +
  geom_vline(xintercept = 0.5, linetype = "dashed") +
  scale_fill_manual(name = '',
                    values = values) +
  scale_color_manual(name = '',
                     values = values2) +
  theme(legend.position = 'right',
        legend.background = element_rect(color = NA, fill = NA),
        legend.key=element_blank(),
        legend.spacing.y = unit(0.05,'cm'),
        legend.key.size = unit(0.4,'cm')) 

p2


pdf(file = paste0('Output/bottleneck.pdf'), width = 3, height = 1.8)
print(p1)
print(p2)
dev.off() 
