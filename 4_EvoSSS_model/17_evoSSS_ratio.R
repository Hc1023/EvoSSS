rm(list = ls())
library(rstan)
library(ggplot2)
library(scales)
library(ggplot2)
library(scales)
library(ggnewscale)
library(tidyverse)
library(dplyr)

load('evoSSS_chain.rdata')

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
  fit = fitlist[[n]]
  
  posterior = rstan::extract(fit)
  contact = mean(posterior$contact)
  N = seed_vec * contact + 1
  Onsets_mat = simu(seed_mat_I1, seed_mat_I2, N, poolday, pars)
  Onsets_mat_list[[n+1]] = Onsets_mat
}

# calibration
for (i in 1:length(Onsets_mat_list)) {
  Onsets_mat_list[[i]] = Onsets_mat_list[[i]]/28
}


dfplot_simu = data.frame()
for (i in 1:25) {
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

dfp_all = data.frame()
for (i in 0:5) {
  x1 = dfplot_simu[dfplot_simu$group == paste0('A',i),]
  x2 = dfplot_simu[dfplot_simu$group == paste0('B',i),]
  dfp = data.frame(p1 = x1$y/(x1$y + x2$y), cycle = i, x = poolday*(i+1))
  dfp_all = rbind(dfp_all, dfp)
}

library(ggdist)
dfp_all$cycle = factor(dfp_all$cycle)
dfp_all$x = dfp_all$x + as.Date('2019-12-31')
p = ggplot(dfp_all, aes(x = x, y = p1)) +
  stat_halfeye(aes(fill = cycle), color = 'black', alpha = 0.5,
               adjust = 1, justification = -0.22,
               .width = 0, width = 60,
               point_colour = NA,
               show.legend = 'none') +
  geom_boxplot(aes(fill=cycle), color = 'black', alpha = 0.5,
               linewidth = 0.12, 
               width = 20, outlier.shape = NA) + 
  ggsci::scale_color_nejm() +
  ggsci::scale_fill_nejm(name = 'Cycle') +
  theme_bw() +
  scale_x_date(breaks = seq(as.Date('2019-12-01'), as.Date('2022-11-01'), by="3 months"),
               minor_breaks = seq(as.Date('2019-12-01'), as.Date('2022-11-01'), by ='1 month'),
               date_labels = "%y-%b") +
  xlab('') + ylab(expression(p[A])) +
  theme(legend.key.size = unit(0.25,'cm'))

pdf(paste0("Output/evoSSS_ratio.pdf"), width = 3, height = 1.3)
print(p)
dev.off()
