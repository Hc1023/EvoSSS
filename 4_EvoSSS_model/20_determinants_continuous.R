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
  I3 <- states_old[,4]
  
  
  beta1 = pars[1]
  c = pars[2]
  beta2 = beta1 * c
  beta3 = beta2 * c
  gamma = pars[3]
  mu1 = 1e-6
  
  S_new = S - beta1*S*I1/N - beta2*S*I2/N - beta3*S*I3/N
  Onset1 = beta1*S*I1/N
  Onset2 = beta2*S*I2/N
  Onset3 = beta3*S*I2/N
  
  if(sum(Onset2) >= 0.2*sum(Onset1+Onset2)){
    mu2 = 1e-6
  }else{
    mu2 = 0
  }
  mu1=0
  mu2=0
  I1_new = I1 + beta1*S*I1/N - gamma*I1 - mu1*I1
  I2_new = I2 + beta2*S*I2/N - gamma*I2 + mu1*I1 - mu2*I2
  I3_new = I3 + beta3*S*I3/N - gamma*I3 + mu2*I2
  Onset1 = beta1*S*I1/N 
  Onset2 = beta2*S*I2/N + mu1*I1
  Onset3 = beta3*S*I3/N + mu2*I2
  
  mat = cbind(S_new, I1_new, I2_new, I3_new, 
              Onset1, Onset2, Onset3)
  return(mat)
}

simu <- function(seed_mat_I1, seed_mat_I2, seed_mat_I3,  N, poolday, pars) {
  # Initial conditions
  I1_old = seed_mat_I1[1,]
  I2_old = seed_mat_I2[1,]
  I3_old = seed_mat_I3[1,]
  S_old = N - I1_old - I2_old -I3_old
  states_old = cbind(S_old, I1_old, I2_old, I3_old)
  
  # Simulate the dynamics over 
  ndays = 200
  
  Onsets_mat <- matrix(0, ndays, 3)
  
  for (t in 2:ndays) {
    states_old = update_fun(pars = pars, states_old = states_old, N = N)
    if(t <= nrow(seed_mat_I1)){
      states_old[t,2:4] = states_old[t,2:4] + 
        c(seed_mat_I1[t,t], seed_mat_I2[t,t], seed_mat_I3[t,t])
    }
    Onsets_mat[t,] = c(sum(states_old[,5]), sum(states_old[,6]),
                       sum(states_old[,7]))
  }
  
  return(Onsets_mat)
}


b = 0.4
c = 1.12
pars = c(b, c, 0.157)

poolday = 30
# The initial cycle - epidemic outbreak
seed_vec = matrix(0,3,2)
seed_vec[,1] = round(34*c(1,0,0))
seed_mat_I1 = diag(seed_vec[1,])
seed_mat_I2 = diag(seed_vec[2,])
seed_mat_I3 = diag(seed_vec[3,])
N = rep(32583, 2)

Onsets_mat_list = list()
Onsets_mat = simu(seed_mat_I1, seed_mat_I2, seed_mat_I3, N, poolday, pars)
Onsets_mat_list[[1]] = Onsets_mat


for (n in 1:24) {
  Onsets_mat = Onsets_mat_list[[n]]
  
  Onset1 = Onsets_mat[poolday + 1:poolday, 1]
  Onset2 = Onsets_mat[poolday + 1:poolday, 2]
  Onset3 = Onsets_mat[poolday + 1:poolday, 3]
  
  mobility = rep(0.01,30) # Mobility: Control force
  Mobility_matrix = diag(mobility)
  
  seed_vec =  (Onset1 + Onset2 + Onset3) %*% Mobility_matrix %>% as.numeric()
  
  p1 = Onset1/(Onset1 + Onset2 + Onset3)
  p2 = Onset2/(Onset1 + Onset2 + Onset3)
  p3 = Onset3/(Onset1 + Onset2 + Onset3)
  seed_matrix = rbind(seed_vec * p1, seed_vec * p2, seed_vec *p3)
  seed_mat_I1 = diag(seed_matrix[1,])
  seed_mat_I2 = diag(seed_matrix[2,])
  seed_mat_I3 = diag(seed_matrix[3,])
  fit = fitlist[[n]]

  posterior = rstan::extract(fit)
  contact = mean(posterior$contact)
  N  = contact*seed_vec
  # N = seed_vec/34 * rep(32583, 30)+1
  Onsets_mat = simu(seed_mat_I1, seed_mat_I2, seed_mat_I3, N, poolday, pars)
  Onsets_mat_list[[n+1]] = Onsets_mat
}

# # calibration
# for (i in 1:length(Onsets_mat_list)) {
#   Onsets_mat_list[[i]] = Onsets_mat_list[[i]]/28
# }


dfplot_simu = data.frame()

for (i in 1:25) {
  Onsets_mat = Onsets_mat_list[[i]]
  n = i-1
  
  dfplot_simu1 = data.frame(x = rep(1:nrow(Onsets_mat)+poolday*n,3),
                            y = c(Onsets_mat[,1],Onsets_mat[,2], Onsets_mat[,3]),
                            group = rep(paste0(c('A','B','C'),n), 
                                        each = nrow(Onsets_mat)),
                            color = rep(c('A','B','C'), 
                                        each = nrow(Onsets_mat)))
  dfplot_simu = rbind(dfplot_simu, dfplot_simu1)
}

df2 = dfplot_simu %>% group_by(x, color) %>% 
  summarise(y = sum(y)) %>%
  as.data.frame()

# Normalizing proportions for each date
data <- df2 %>%
  group_by(x) %>%
  mutate(p = y/sum(y))


ggplot(data, aes(x = x, y = y, fill = color)) +
  geom_area() +
  # scale_y_continuous(trans='log10') +
  # coord_cartesian(ylim = c(2,max(data$y)^3)) +
  labs(x = "Date", y = "Proportion", fill = "Variant") +
  theme_minimal() +
  theme(legend.position = "right",
        plot.title = element_text(hjust = 0.5))


data2 = df2
data2[data2$color == 'A','y'] = df2[df2$color == 'A','y']*100

data2 <- data2 %>%
  group_by(x) %>%
  mutate(p = y/sum(y))

data3 = df2
data3[data3$color == 'A','y'] = df2[df2$color == 'A','y']*10000
data3 <- data3 %>%
  group_by(x) %>%
  mutate(p = y/sum(y))



dfdata = rbind(cbind(data[data$color == 'B',],group = '1'),
               cbind(data2[data$color == 'B',],group = '2'),
               cbind(data3[data$color == 'B',],group = '3'))
ggplot() +
  geom_area(data = data, aes(x = x, y = p, fill = color),
            position = 'fill', alpha = 0.2) +
  geom_area(data = data2, aes(x = x, y = p, fill = color),
            position = 'fill', alpha = 0.2) +
  geom_area(data = data3, aes(x = x, y = p, fill = color),
            position = 'fill', alpha = 0.2) +
  geom_line(data= dfdata, aes(x = x, y = p, 
                              group = group, color = group)) +
  theme_minimal() +
  theme(legend.position = "right",
        plot.title = element_text(hjust = 0.5))
