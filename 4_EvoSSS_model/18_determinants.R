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
  
  S_new <- S
  for (i in 1:n) {
    S_new <- S_new - betas[i] * S * I[,i] / N
  }
  
  I_new <- matrix(0, nrow(states_old), n)
  I_new[,1] <- I[,1] + betas[1] * S * I[,1] / N - gamma * I[,1] 
  for (i in 2:n) {
    I_new[,i] <- I[,i] + betas[i] * S * I[,i] / N - gamma * I[,i] 
  }
  
  Onset_new <- matrix(0, nrow(states_old), n)
  Onset_new[,1] <- betas[1] * S * I_new[,1] / N
  
  for (i in 2:n) {
    Onset_new[,i] <- betas[i] * S * I_new[,i] / N 
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

determinant_fun = function(parsmat){
  
  poolday = 30
  # The initial cycle - epidemic outbreak
  n= length(parsmat[1,])
  
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
  
  Onsets_mat = simu(seed_mats, N, poolday, parsmat[1,], n)
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
  
    Onsets_mat = simu(seed_mats, N, poolday, parsmat[j+1,], n)
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
  
  df2 = dfplot_simu %>% group_by(x, color) %>% 
    summarise(y = sum(y)) %>%
    as.data.frame()

  
  # Normalizing proportions for each date
  data <- df2 %>%
    group_by(x) %>%
    mutate(p = y/sum(y)) %>%
    as.data.frame()
  
  levels = rev(unique(data$color))
  
  data$color = factor(data$color,levels = levels)
  
  return(data)
}
parsmat = matrix(0.3,nrow = 25, ncol = 3)
parsmat[10:25,2] = 0.36
parsmat[15:25,3] = 0.5
data = determinant_fun(parsmat)
ggplot() +
  geom_area(data = data, 
            aes(x = x, y = p, fill = color),
            position = 'fill')

dfall = data %>% group_by(x) %>% 
  summarise(yall = sum(y))

ggplot() +
  geom_line(data = data, aes(x = x, y = y/28, group = color, color= color)) +
  geom_line(data = dfall, aes(x = x, y = yall/28)) +
  scale_y_continuous(trans='log10') +
  coord_cartesian(ylim = c(2,max(data$y)/28)) +
  labs(x = "Date", y = "Proportion") +
  theme_minimal() +
  theme(legend.position = "right",
        plot.title = element_text(hjust = 0.5))

