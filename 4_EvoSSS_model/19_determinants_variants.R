rm(list = ls())
library(rstan)
library(ggplot2)
library(scales)
library(ggnewscale)
library(tidyverse)
library(dplyr)

load('evoSSS_chain.rdata')


update_fun <- function(pars, states_old, N, n) {
  S <- states_old[,1]
  I <- states_old[,2:(n+1)]
  Onset <- states_old[,(n+2):(2*n+1)]
  
  
  betas = pars

  gamma = 0.157
  mu1 = 0.01
  mu = rep(mu1, n)
  
  # mu[1] = mu1
  
  S_new <- S
  for (i in 1:n) {
    S_new <- S_new - betas[i] * S * I[,i] / N
  }
  
  I_new <- matrix(0, nrow(states_old), n)
  I_new[,1] <- I[,1] + betas[1] * S * I[,1] / N - gamma * I[,1] - mu[1] * I[,1]
  for (i in 2:n) {
    I_new[,i] <- I[,i] + betas[i] * S * I[,i] / N - gamma * I[,i] + mu[i-1] * I[,i-1] - mu[i] * I[,i]
  }
  
  Onset_new <- matrix(0, nrow(states_old), n)
  Onset_new[,1] <- betas[1] * S * I_new[,1] / N - mu[1] * I[,1]
  
  for (i in 2:n) {
    Onset_new[,i] <- betas[i] * S * I_new[,i] / N + mu[i-1] * I[,i-1] - mu[i] * I[,i]
  }
  
  mat <- cbind(S_new, I_new, Onset_new)

  return(mat)
}

simu <- function(seed_mats, N, poolday, pars, n) {
  # Initial conditions
  I_old <- matrix(0, ncol=n, nrow=nrow(seed_mats[[1]]))
  S_old <- N - rowSums(I_old)
  
  
  for (i in 1:n) {
    I_old[,i] <- seed_mats[[i]][1,]
  }
  
  states_old <- cbind(S_old, I_old, matrix(0, ncol=n, nrow=nrow(seed_mats[[1]])))
  
  ndays <- 200
  Onsets_mat <- matrix(0, ndays, n)

  for (t in 2:ndays) {

    states_old <- update_fun(pars = pars, states_old = states_old, N = N, n = n)
    if (t <= nrow(seed_mats[[1]])) {
      for (i in 1:n) {
        states_old[t,(i+1):(i+1)] <- states_old[t,(i+1):(i+1)] + seed_mats[[i]][t,t]
      }
    }
    Onsets_mat[t,] <- colSums(states_old[,(n+2):(2*n+1)])
  }
  
  return(Onsets_mat)
}

determinant_fun = function(pars){
  
  poolday = 30
  # The initial cycle - epidemic outbreak
  n= length(pars)

  seed_vec = c(34,0)
  probs <- c(1,rep(0,n-1))
  
  seed_matrix <- matrix(0, nrow=n, ncol=length(seed_vec))
  for (i in 1:n) {
    seed_matrix[i,] <- seed_vec * probs[[i]]
  }
  seed_mats <- list()
  for (i in 1:n) {
    seed_mats[[i]] <- diag(seed_matrix[i,])
  }
  N = rep(32583, 2)
  Onsets_mat_list = list()
  b = 0.4
  c = 1.05
  mu = 0.01

  Onsets_mat = simu(seed_mats, N, poolday, pars, n)
  Onsets_mat_list[[1]] = Onsets_mat
  
  for (j in 1:24) {
    Onsets_mat <- Onsets_mat_list[[j]]

    # Generalized extraction of Onset columns for all variants
    Onsets <- list()
    for (i in 1:n) {
      Onsets[[i]] <- Onsets_mat[poolday + 1:poolday, i] + 1
    }
    
    # Mobility control force
    mobility <- rep(0.01, 30) 
    Mobility_matrix <- diag(mobility)
    
    # Create seed_vec considering all variants
    seed_vec <- rowSums(sapply(Onsets, function(Onset) Onset)) %*% Mobility_matrix %>% as.numeric()
    
    # Calculate probabilities for each variant
    probs <- lapply(Onsets, function(Onset) Onset / rowSums(do.call(cbind, Onsets)))
    seed_matrix <- matrix(0, nrow=n, ncol=length(seed_vec))
    for (i in 1:n) {
      seed_matrix[i,] <- seed_vec * probs[[i]]
    }
    
    # Create seed_matrices for each variant
    seed_mats <- list()
    for (i in 1:n) {
      seed_mats[[i]] <- diag(seed_matrix[i,])
    }
    
    fit = fitlist[[j]]
    posterior = rstan::extract(fit)
    contact = mean(posterior$contact)
    N  = contact*seed_vec
    print(contact)

    # N = seed_vec/34 * rep(32583, 30)+1

    Onsets_mat = simu(seed_mats, N, poolday, pars, n)
    Onsets_mat_list[[j+1]] = Onsets_mat
  }
  
  dfplot_simu = data.frame()
  for (i in 1:25) {
    Onsets_mat = Onsets_mat_list[[i]]
    k = i-1
    
    dfplot_simu1 = data.frame(x = rep(1:nrow(Onsets_mat)+poolday*k,n),
                              y = as.vector(Onsets_mat),
                              group = rep(paste0(1:n,'_',k), 
                                          each = nrow(Onsets_mat)),
                              color = rep(1:n, 
                                          each = nrow(Onsets_mat)))
    dfplot_simu = rbind(dfplot_simu, dfplot_simu1)
  }
  str(dfplot_simu)

  df2 = dfplot_simu %>% group_by(x, color) %>% 
    summarise(y = sum(y)) %>%
    as.data.frame()
  str(df2)

  # Normalizing proportions for each date
  data <- df2 %>%
    group_by(x) %>%
    mutate(p = y/sum(y)) %>%
    as.data.frame()
  
  str(data)
  levels = unique(data$color)
  levels[1:2] = c(2,1)
  data$color = factor(data$color,levels = levels)

  return(data)
}
pars = c(0.4,0.2,0.41)
data = determinant_fun(pars)
ggplot() +
  geom_area(data = data, 
            aes(x = x, y = p, fill = color),
            position = 'fill')
dfall = data %>% group_by(x) %>% summarise(yall = sum(y))
ggplot() +
  geom_line(data = data, aes(x = x, y = y/28, group = color, color= color)) +
  geom_line(data = dfall, aes(x = x, y = yall/28)) +
  scale_y_continuous(trans='log10') +
  coord_cartesian(ylim = c(2,max(data$y)/28)) +
  labs(x = "Date", y = "Proportion") +
  theme_minimal() +
  theme(legend.position = "right",
        plot.title = element_text(hjust = 0.5))

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
        plot.title = element_text(hjust = 0.5)) +
  scale_x_date(breaks = seq(as.Date('2020-01-01'), as.Date('2022-11-01'), by="6 months"),
               minor_breaks = seq(as.Date('2019-12-01'), as.Date('2022-11-01'), by ='1 month'),
               date_labels = "%y-%b") +
  scale_y_continuous(breaks = c(0,0.5,1),
                     labels = c('0.0','0.5','1.0')) +
  coord_cartesian(xlim = c(as.Date('2020-01-31'), as.Date('2021-10-31')),
                  ylim = c(0.04,0.96)) +
  xlab('') + ylab('') 

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
        plot.title = element_text(hjust = 0.5)) +
  scale_x_date(breaks = seq(as.Date('2020-01-01'), as.Date('2022-11-01'), by="6 months"),
               minor_breaks = seq(as.Date('2019-12-01'), as.Date('2022-11-01'), by ='1 month'),
               date_labels = "%y-%b") +
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
