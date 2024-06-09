rm(list = ls())
library(rstan)
library(ggplot2)
library(scales)
library(ggnewscale)
library(tidyverse)
library(dplyr)


update_fun = function(pars, states_old, N){
  
  S <- states_old[,1]
  I1 <- states_old[,2]
  I2 <- states_old[,3]
  
  beta1 = pars[1]
  c = pars[2]
  beta2 = beta1 * c
  gamma = pars[3]
  mu1 = pars[4]
  
  S_new = S - beta1*S*I1/N - beta2*S*I2/N
  I1_new = I1 + beta1*S*I1/N - gamma*I1 - mu1*I1
  I2_new = I2 + beta2*S*I2/N - gamma*I2 + mu1*I1
  Onset1 = beta1*S*I1/N 
  Onset2 = beta2*S*I2/N + mu1*I1
  mat = cbind(S_new, I1_new, I2_new,
              Onset1, Onset2)
  return(mat)
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



determinant_fun = function(pars){
  
  poolday = 30
  # The initial cycle - epidemic outbreak
  variants_num = 2
  seed_vec = matrix(0,variants_num,2)
  seed_vec[,1] = round(34*c(1,0))
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
    
    p1 = Onset1/(Onset1 + Onset2)
    p2 = Onset2/(Onset1 + Onset2)
    
    seed_matrix = rbind(seed_vec * p1, seed_vec * p2)
    seed_mat_I1 = diag(seed_matrix[1,])
    seed_mat_I2 = diag(seed_matrix[2,])
    
    N = seed_vec/34 * rep(32583, 30)+1
    Onsets_mat = simu(seed_mat_I1, seed_mat_I2, N, poolday, pars)
    Onsets_mat_list[[n+1]] = Onsets_mat
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
  
  df2 = dfplot_simu %>% group_by(x, color) %>% 
    summarise(y = sum(y)) %>%
    as.data.frame()
  df2$x = df2$x #+ as.Date('2019-12-31')
  values = c(hue_pal()(3)[1], hue_pal()(3)[3])
  
  # Normalizing proportions for each date
  data <- df2 %>%
    group_by(x) %>%
    mutate(p = y/sum(y))
  
  return(data)
}
b = 0.4
c = 1
pars = c(b, c, 0.157,1e-6)
data = determinant_fun(pars)
b = 0.4
c = 1.05
pars = c(b, c, 0.157,1e-6)
data2 = determinant_fun(pars)
b = 0.4
c = 1.1
pars = c(b, c, 0.157,1e-6)
data3 = determinant_fun(pars)
b = 0.4
c = 1.2
pars = c(b, c, 0.157,1e-6)
data4 = determinant_fun(pars)



data$x = data$x + as.Date('2019-12-31')
data2$x = data2$x + as.Date('2019-12-31')
data3$x = data3$x + as.Date('2019-12-31')
data4$x = data4$x + as.Date('2019-12-31')
dfdata = rbind(cbind(data[data$color == 'B',],group = '1'),
               cbind(data2[data$color == 'B',],group = '2'),
               cbind(data3[data$color == 'B',],group = '3'),
               cbind(data4[data$color == 'B',],group = '4'))
values = rev(hue_pal()(4))
values2 = hue_pal()(2)

p = ggplot() +
  geom_area(data = data, aes(x = x, y = p, fill = color),
            position = 'fill') +
  geom_area(data = data2, aes(x = x, y = p, fill = color),
            position = 'fill') +
  geom_area(data = data3, aes(x = x, y = p, fill = color),
            position = 'fill') +
  geom_area(data = data4, aes(x = x, y = p, fill = color),
            position = 'fill') +
  geom_line(data= dfdata, 
            aes(x = x, y = p, 
                group = group, color = group),
            linewidth = 2) +
  scale_color_manual(values = alpha(values, 0.8)) +
  scale_fill_manual(values = alpha(values2, 0.2)) +
  theme_minimal() +
  theme(legend.position = "right",
        plot.title = element_text(hjust = 0.5),
        panel.grid.major = element_blank(),  # Remove major grid lines
        panel.grid.minor = element_blank(),  # Remove minor grid lines
        axis.ticks = element_blank()) +
  scale_y_continuous(breaks = c(0,0.5,1),
                     labels = c('0.0','0.5','1.0')) +
  coord_cartesian(xlim = c(as.Date('2020-01-31'), as.Date('2021-10-31')),
                  ylim = c(0.04,0.96)) +
  xlab('') + ylab('') 
p
pdf(paste0("Output/demerminants_tranmissibility.pdf"), width = 2.5, height = 1.2)
print(p)
dev.off()


b = 0.4
c = 1.1
mu = 1e-12
pars = c(b, c, 0.157,mu)
data = determinant_fun(pars)
mu = 1e-6
pars = c(b, c, 0.157,mu)
data2 = determinant_fun(pars)
mu = 2e-6
pars = c(b, c, 0.157,mu)
data3 = determinant_fun(pars)
mu = 2e-3
pars = c(b, c, 0.157,mu)
data4 = determinant_fun(pars)


data$x = data$x + as.Date('2019-12-31')
data2$x = data2$x + as.Date('2019-12-31')
data3$x = data3$x + as.Date('2019-12-31')
data4$x = data4$x + as.Date('2019-12-31')
dfdata = rbind(cbind(data[data$color == 'B',],group = '1'),
               cbind(data2[data$color == 'B',],group = '2'),
               cbind(data3[data$color == 'B',],group = '3'),
               cbind(data4[data$color == 'B',],group = '4'))
values = rev(hue_pal()(4))
values2 = hue_pal()(2)

p = ggplot() +
  geom_area(data = data, aes(x = x, y = p, fill = color),
            position = 'fill') +
  geom_area(data = data2, aes(x = x, y = p, fill = color),
            position = 'fill') +
  geom_area(data = data3, aes(x = x, y = p, fill = color),
            position = 'fill') +
  geom_area(data = data4, aes(x = x, y = p, fill = color),
            position = 'fill') +
  geom_line(data= dfdata, 
            aes(x = x, y = p, 
                group = group, color = group),
            linewidth = 2) +
  scale_color_manual(values = alpha(values, 0.8)) +
  scale_fill_manual(values = alpha(values2, 0.2)) +
  theme_minimal() +
  theme(legend.position = "right",
        panel.grid.major = element_blank(),  # Remove major grid lines
        panel.grid.minor = element_blank(),  # Remove minor grid lines
        axis.ticks = element_blank()) +
  scale_y_continuous(breaks = c(0,0.5,1),
                     labels = c('0.0','0.5','1.0')) +
  coord_cartesian(xlim = c(as.Date('2020-01-31'), as.Date('2021-10-31')),
                  ylim = c(0.04,0.96)) +
  xlab('') + ylab('') 
p

pdf(paste0("Output/demerminants_mutation.pdf"), width = 2.5, height = 1.2)
print(p)
dev.off()

ggplot(data, aes(x = x, y = y, fill = color)) +
  geom_area() +
  # scale_y_continuous(trans='log10') +
  # coord_cartesian(ylim = c(2,max(data$y)^3)) +
  labs(x = "Date", y = "Proportion", fill = "Variant") +
  theme_minimal() +
  theme(legend.position = "right",
        plot.title = element_text(hjust = 0.5))

ggplot(data, aes(x = x, y = p, fill = color)) +
  geom_area(position = 'fill') +
  theme_minimal() +
  theme(legend.position = "right",
        plot.title = element_text(hjust = 0.5))
