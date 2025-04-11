rm(list = ls())
library(rstan)
library(ggplot2)
library(scales)
library(ggnewscale)
library(tidyverse)
library(dplyr)

df = read.csv('F6B_influenza.csv')
observed_matrix = data.frame(v1 = df$A.H1*df$Total,
                             v2 = df$A.H3*df$Total,
                             v3 = df$A.H5*df$Total,
                             v4 = df$B.Victoria*df$Total)
voc = c('A.H1','A.H3','A.H5','B.Victoria')
rownames(observed_matrix) = df$X

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
if(F){
  fitlist = list()
}

determinant_fun = function(cond = T, ifsimu  = T, n_simu = 1){
  # cond = F; ifsimu  = F; n_simu = 1
  n = 4
  poolday = 30
  nday = 100
  seed_matrix <- matrix(0, nrow=n, ncol=poolday)
  seed_matrix[,1] = rep(n,n)
  if(ifsimu){
    seed_matrix[,1] = rmultinom(1, n*n, rep(1/n,n))
  }
  seed_vec = colSums(seed_matrix)
  seed_mats <- list()
  
  for (i in 1:n) {
    seed_mats[[i]] <- diag(seed_matrix[i,])
  }
  
  if(cond){
    stan_data <- list(
      poolday = poolday,
      nday = nday,
      expected_total = unlist(round(observed_matrix[1,1:4])),
      seed_mat_I1 = seed_mats[[1]],
      seed_mat_I2 = seed_mats[[2]],
      seed_mat_I3 = seed_mats[[3]],
      seed_mat_I4 = seed_mats[[4]],
      seed_vec = seed_vec,
      gamma = 0.157,
      pars_last = c(200, rep(0.3,4))
    )
    # Fit the model
    fit <- stan(file = 'F6B_FLU.stan', data = stan_data, 
                iter = 3000, chains = 1, warmup = 2000,
                verbose = TRUE)
    fitlist[[1]] = fit
  }
  
  
  Onsets_mat_list = list()
  fit= fitlist[[1]]
  
  posterior = rstan::extract(fit)
  
  pars_last = c(mean(posterior$contact), 
                mean(posterior$beta1),
                mean(posterior$beta2), 
                mean(posterior$beta3),
                mean(posterior$beta4))
  if(ifsimu){
    pars_last = c(posterior$contact[n_simu], 
                  posterior$beta1[n_simu],
                  posterior$beta2[n_simu], 
                  posterior$beta3[n_simu],
                  posterior$beta4[n_simu])
  }
  Onsets_mat = simu(seed_mats, 
                    N = seed_vec * pars_last[1] + 1, 
                    poolday, pars = pars_last[-1], n)
  if(!ifsimu){
    fonset = data.frame(x = rep(1:nday,4), 
                        y = c(Onsets_mat[1:nday,1],
                              Onsets_mat[1:nday,2],
                              Onsets_mat[1:nday,3],
                              Onsets_mat[1:nday,4]),
                        group = factor(rep(voc, each = nday), 
                                       levels = voc))
    
    tmp = fonset %>% group_by(group) %>% 
      summarise(y = sum(y)) %>% 
      data.frame()
    tmp$x = c(55,60,65,70)
    
    observed_df = data.frame(y = unlist(observed_matrix[1,]))
    
    observed_df$x = c(25,30,35,40)
    observed_df$group = voc
    observed_df = rbind(observed_df, tmp)
    
    ggplot() +
      geom_line(data = fonset,
                aes(x = x, y = y, 
                    group = group, color = group))+
      geom_col(data = observed_df, 
               aes(x = x, y = y/30, 
                   group = group, 
                   color = group, fill = group), 
               width = 5)
    
  }
  
  Onsets_mat_list[[1]] = Onsets_mat
  
  for (j in 1:36) {
    if(!ifsimu)print(j)
    {
      Onsets_mat = Onsets_mat_list[[j]]
      
      # Generalized extraction of Onset columns for all variants
      Onsets <- list()
      for (i in 1:n) {
        Onsets[[i]] <- Onsets_mat[poolday + 1:poolday, i] + 1e-3
      }
      # Mobility control force
      mobility <- rep(1/30, 30) 
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
      
      if(cond){
        stan_data <- list(
          poolday = poolday,
          nday = nday,
          expected_total = unlist(round(observed_matrix[j+1,1:n])),
          seed_mat_I1 = seed_mats[[1]],
          seed_mat_I2 = seed_mats[[2]],
          seed_mat_I3 = seed_mats[[3]],
          seed_mat_I4 = seed_mats[[4]],
          seed_vec = seed_vec,
          gamma = 0.157,
          pars_last = pars_last
        )
        # Fit the model
        fit <- stan(file = 'F6B_FLU.stan', data = stan_data, 
                    iter = 3000, chains = 1, warmup = 2000,
                    verbose = TRUE)
        fitlist[[j+1]] = fit
      }
      
      fit = fitlist[[j+1]]
      posterior = rstan::extract(fit)
      pars_last = c(mean(posterior$contact), 
                    mean(posterior$beta1),
                    mean(posterior$beta2),
                    mean(posterior$beta3),
                    mean(posterior$beta4))
      if(ifsimu){
        pars_last = c(posterior$contact[n_simu], 
                      posterior$beta1[n_simu],
                      posterior$beta2[n_simu], 
                      posterior$beta3[n_simu],
                      posterior$beta4[n_simu])
      }
      Onsets_mat = simu(seed_mats, 
                        N = seed_vec * pars_last[1] + 1, 
                        poolday, pars = pars_last[-1], n)
      
      if(!ifsimu){
        fonset = data.frame(x = rep(1:nday,4), 
                            y = c(Onsets_mat[1:nday,1],
                                  Onsets_mat[1:nday,2],
                                  Onsets_mat[1:nday,3],
                                  Onsets_mat[1:nday,4]),
                            group = factor(rep(voc, each = nday), 
                                           levels = voc))
        
        
        tmp = fonset %>% group_by(group) %>% 
          summarise(y = sum(y)) %>% 
          data.frame()
        tmp$x = c(55,60,65,70)
        
        observed_df = data.frame(y = unlist(observed_matrix[j+1,]))
        
        observed_df$x = c(25,30,35,40)
        observed_df$group = voc
        observed_df = rbind(observed_df, tmp)
        
        ggplot() +
          geom_line(data = fonset,
                    aes(x = x, y = y, 
                        group = group, color = group))+
          geom_col(data = observed_df, 
                   aes(x = x, y = y/50, 
                       group = group, 
                       color = group, fill = group), 
                   width = 5)
      }
    }
    Onsets_mat_list[[j+1]] = Onsets_mat
  }
  
  return(Onsets_mat_list)
  # j = 35 last simulation -- 24-July
  
  
}

pred_sensitivity = function(s, Onsets_mat_list, ifsimu = F){
  n = 4
  poolday = 30
  nday = 100
  # s=1,2,3,4
  j = 36
  fit = fitlist[[j+1-s]]
  posterior = rstan::extract(fit)
  pars_last = c(mean(posterior$contact), 
                mean(posterior$beta1),
                mean(posterior$beta2),
                mean(posterior$beta3),
                mean(posterior$beta4))
  if(ifsimu){
    pars_last = c(posterior$contact[n_simu], 
                  posterior$beta1[n_simu],
                  posterior$beta2[n_simu], 
                  posterior$beta3[n_simu],
                  posterior$beta4[n_simu])
  }
  
  for (j in 37:46) {
    # print(j)
    {
      Onsets_mat = Onsets_mat_list[[j]]
      
      # Generalized extraction of Onset columns for all variants
      Onsets <- list()
      for (i in 1:n) {
        Onsets[[i]] <- Onsets_mat[poolday + 1:poolday, i] + 1e-3
      }
      # Mobility control force
      mobility <- rep(1/30, 30) 
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
      
      fit = fitlist[[j+1-12]]
      posterior = rstan::extract(fit)
      pars_last[1] = mean(posterior$contact)
      if(ifsimu){
        pars_last[1] = posterior$contact[n_simu]
      }
      Onsets_mat = simu(seed_mats, 
                        N = seed_vec * pars_last[1] + 1, 
                        poolday, pars = pars_last[-1], n)
      
    }
    Onsets_mat_list[[j+1]] = Onsets_mat
  }
  dfplot_simu = data.frame()
  for (i in 1:47) {
    Onsets_mat = Onsets_mat_list[[i]]
    k = i-1
    
    dfplot_simu1 = data.frame(x = rep(1:nrow(Onsets_mat)+poolday*k,n),
                              y = as.vector(Onsets_mat),
                              group = rep(paste0(1:n,'_',k), 
                                          each = nrow(Onsets_mat)),
                              color = rep(voc, 
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
  
  data$group = factor(data$color,levels = levels)
  # save(fitlist, Onsets_mat_list, file = 'flu.rdata')
  return(data)
}
load('F6B_fluH5.rdata')

dfob = data.frame(y = c(as.matrix(observed_matrix)),
                  x = rep(as.Date(rownames(observed_matrix)),
                          4),
                  group = rep(voc, each = nrow(observed_matrix)))

Onsets_mat_list0 = determinant_fun(cond = F, ifsimu  = F)

df = read.csv('F6B_influenza.csv')
df2 = data.frame(x = as.Date(rep(df$X, 3)),
                 y = c(observed_matrix$v1/rowSums(observed_matrix), 
                       (observed_matrix$v1 + observed_matrix$v2)/rowSums(observed_matrix),
                       (observed_matrix$v1 + observed_matrix$v2 + observed_matrix$v3)/rowSums(observed_matrix)),
                 group = rep(voc[1:3], each = nrow(df)))

values = c('#98afc7','#0041c2', "#ffa500",'#a37ca1')

datalist = list()
for (s in 0:5) {
  data = pred_sensitivity(s = s, Onsets_mat_list0)
  data$date = as.Date('2021-07-01') + data$x
  data$group = factor(data$group, levels = rev(voc))
  datalist[[s+1]] = data
}

pdf(paste0("Output/S12E_fluH5_pred2_prevalence.pdf"), width = 1.5, height = 1.2)
for (s in 0:5) {
  data = datalist[[s+1]]
  p2 =  ggplot() +
    geom_area(data = data, 
              aes(x = date, y = p, fill = group),
              position = 'fill') +
    scale_fill_manual(name="", breaks = rev(levels(data$group)),
                      values = alpha(values, 1)) +
    scale_x_date(breaks = seq(as.Date('2021-01-01'), as.Date('2025-05-01'), by="6 months"),
                 minor_breaks = seq(as.Date('2021-01-01'), as.Date('2025-05-01'), by ='1 month'),
                 date_labels = "%y-%b",
                 expand = c(0,0)) +
    scale_y_continuous(limits = c(0,1), expand = c(0,0),
                       breaks = c(0,0.5,1)) +
    ylab('Prevalence') + xlab('') +
    theme(legend.position = 'none',
          panel.background = element_blank(),
          plot.background = element_blank()) +
    coord_cartesian(xlim = c(as.Date('2024-07-01'), as.Date('2025-05-31')))
  
  print(p2)
}

dev.off()
