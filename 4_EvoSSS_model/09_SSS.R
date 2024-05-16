rm(list = ls())
library(rstan)
library(ggplot2)
library(scales)
library(ggplot2)
library(scales)
library(ggnewscale)

df = read.csv('../3_Epidemiological_analysis/Covid19CasesGISAID.csv')

df$Var1 = as.Date(df$Var1)
df = df[df$Var1 < as.Date('2021-12-01'),]
df = df[df$Mutations %in% c('Lineage A', 'Lineage B'),]

update_fun = function(pars, states_old){
  
  S <- states_old[,1]
  I1 <- states_old[,2]
  I2 <- states_old[,3]
  N = states_old[,4]
  beta1 = pars[1]
  beta2 = pars[2]
  gamma = pars[3]
  
  S_new = S - beta1*S*I1/N - beta2*S*I2/N
  I1_new = I1 + beta1*S*I1/N - gamma*I1
  I2_new = I2 + beta2*S*I2/N - gamma*I2
  Onset1 = beta1*S*I1/N
  Onset2 = beta2*S*I2/N
  
  return(data.frame(S_new, I1_new, I2_new, Onset1, Onset2))
}

simu <- function(seed_mat_list, dim, f = update_fun) {
  # Initial conditions
  seed_mat_I1 = seed_mat_list[[1]]
  seed_mat_I2 = seed_mat_list[[2]]
  N_mat = seed_mat_list[[3]]
  I1_old = seed_mat_I1[1,]
  I2_old = seed_mat_I2[1,]
  S_old = N_mat - I1_old - I2_old
  states_old = data.frame(S_old = S_old, I1_old = I1_old, 
                          I2_old = I2_old)
  # str(states_old)
  # Simulate the dynamics over ndays
  ndays = dim*2
  mycol <- c("time", "Onset1", "Onset2")
  Onsets_mat <- matrix(0, ndays, length(mycol))
  Onsets_mat[,1] = 1:ndays
  
  colnames(Onsets_mat) <- mycol
  pars = c(0.379, 0.398, 0.157)
  for (t in 1:(ndays-1)) {
    states_mat = f(pars = pars, states_old = as.matrix(cbind(states_old, N_mat = N_mat)))
    states_old = states_mat[1:3]
    states_old[,2:3] = states_old[,2:3] + cbind(seed_mat_I1[t+1,], seed_mat_I2[t+1,])
    
    Onsets_mat[t+1,-1] = c(sum(states_mat[4]), sum(states_mat[5]))
  }
  
  return(data.frame(Onsets_mat))
}

dim = 100
# The initial cycle - epidemic outbreak
seed_matrix = rbind(rep(0,dim),rep(0,dim))
seed_matrix[,1] = round(34*c(0.4,0.6))


seed_mat_list = list()
seed_mat_list[[1]] = rbind(diag(seed_matrix[1,]), matrix(0, dim, dim))
seed_mat_list[[2]] = rbind(diag(seed_matrix[2,]), matrix(0, dim, dim))
seed_mat_list[[3]] = rep(32583, dim)
Onsets_mat_list = list()
Onsets_mat = simu(seed_mat_list, dim)
Onsets_mat_list[[1]] = Onsets_mat
Onsets_mat = Onsets_mat_list[[1]]
# 
n = 1


poolday = 30
contact_vec = c(10^3,10^3)
for (n in 1:length(contact_vec)) {
  # seeding_time = 1:dim + spacing*(n-1)
  
  Onsets_mat_old = Onsets_mat
  Onset1 = Onsets_mat_old[1:dim+poolday,2]
  Onset2 = Onsets_mat_old[1:dim+poolday,3]
  I_last = Onset1 + Onset2
  mobility = rep(0, dim)
  mobility[seq(1,30,1)] = 0.01  # Mobility: Control force
  Mobility_matrix = diag(mobility)
  
  seed_vec =  I_last %*% Mobility_matrix
  seed_vec = seed_vec[1,]
  p = Onset1/(Onset1 + Onset2)
  seed_matrix = rbind(seed_vec * p, seed_vec * (1-p))
  # seed_matrix[,1] = 0
  seed_mat_list = list()
  seed_mat_list[[1]] = rbind(diag(seed_matrix[1,]), matrix(0, dim, dim))
  seed_mat_list[[2]] = rbind(diag(seed_matrix[2,]), matrix(0, dim, dim))
  seed_mat_list[[3]] = (seed_vec+1)*contact_vec[n]
  Onsets_mat = simu(seed_mat_list, dim)
  Onsets_mat[,1] = Onsets_mat[,1] + sum(poolday*n)
  Onsets_mat_list[[n+1]] = Onsets_mat

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
library(tidyverse)

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

