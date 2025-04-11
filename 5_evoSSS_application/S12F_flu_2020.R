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


determinant_fun = function(cond = T, ifsimu  = T, n_simu = 1, cutoff){
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
  
  Onsets_mat_list[[1]] = Onsets_mat
  
  for (j in 1:cutoff) {
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
      
    }
    Onsets_mat_list[[j+1]] = Onsets_mat
  }
  
  # j = 35 last simulation -- 24-July
  
  for (j in (cutoff+1):36) {
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
  
  # save(fitlist, Onsets_mat_list, file = 'flu.rdata')
  dfplot_simu = data.frame()
  for (i in 1:37) {
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
  return(data)
}

dfob = data.frame(y = c(as.matrix(observed_matrix)),
                  x = rep(as.Date(rownames(observed_matrix)),
                          4),
                  group = rep(voc, each = nrow(observed_matrix)))
dfob$group[dfob$group!=voc[4]] = 'A'
dfob = dfob %>% group_by(group, x) %>% summarise(y = sum(y))

if(F){
  df2_list = list()
  for (n_simu in 1:100) {
    print(n_simu)
    data= determinant_fun(cond = F, ifsimu  = T, n_simu = n_simu, cutoff = 25)
    df2_list[[n_simu]] = data$y
  }
  simu_Onset1 = data.frame(bind_cols(df2_list))
  
  df2_list = list()
  for (n_simu in 1:100) {
    print(n_simu)
    data= determinant_fun(cond = F, ifsimu  = T, n_simu = n_simu, cutoff = 13)
    df2_list[[n_simu]] = data$y
  }
  
  simu_Onset2 = data.frame(bind_cols(df2_list))
  
  save(simu_Onset1, simu_Onset2, file = 'flu_plot1.rdata')
  
}

load('F6B_fluH5.rdata')
load('S8C_flu_2020.rdata')
data = determinant_fun(cond = F, ifsimu  = F, n_simu = n_simu,
                       cutoff = 25) # 13,25)
getks = function(simu_Onset = simu_Onset1,
                 preddate = as.Date('2023-08-01'),
                 d1 = as.Date('2022-08-01'), d2 = as.Date('2024-05-31')){
  plot_data <- data.frame(
    x = data$x,
    group = data$color,
    Fitted = rowMeans(simu_Onset)
  )
  
  plot_data$date = as.Date('2021-07-01') + plot_data$x
  
  k = 30
  pd = 25
  
  plot_data$group[plot_data$group!=voc[4]] = 'A'
  plot_data = plot_data %>% group_by(group, date) %>% 
    summarise(Fitted = sum(Fitted))
  plot_data = plot_data[plot_data$date >= preddate &
                          plot_data$date <= d2,]
  plot_data$date2 = as.character(format(plot_data$date, "%y-%m"))
  plot_data2 = plot_data %>% group_by(group, date2) %>% 
    summarise(Fitted2 = sum(Fitted))
  dfob2 = dfob[dfob$x >= preddate &
                 dfob$x <= d2,]
  r = ks.test(dfob2$y, plot_data2$Fitted2)
  return(r)
  # ks.test(dfob$y, plot_data2$Fitted2)
  # >   ks.test(dfob$y, plot_data2$Fitted2)
  # 
  # Two-sample Kolmogorov-Smirnov test
  # 
  # data:  dfob$y and plot_data2$Fitted2
  # D = 0.081081, p-value = 0.9699
  # alternative hypothesis: two-sided
}

getks(simu_Onset = simu_Onset1,
      preddate = as.Date('2023-08-01'),
      d1 = as.Date('2022-08-01'), d2 = as.Date('2024-05-31'))
getks(simu_Onset2, 
      preddate = as.Date('2022-08-01'),
      d1 = as.Date('2021-08-01'), d2 = as.Date('2023-05-31'))
getks(simu_Onset2, 
      preddate = as.Date('2022-08-01'),
      d1 = as.Date('2021-08-01'), d2 = as.Date('2023-01-31'))


getplot = function(simu_Onset, preddate = as.Date('2023-08-01'),
                   d1 = as.Date('2022-08-01'), d2 = as.Date('2024-05-31')){
  values = c('#98aff7','#a37c91')
  # show_col(values)
  ci_lower <- apply(simu_Onset, 1, quantile, probs = 0.025, na.rm = T)
  ci_upper <- apply(simu_Onset, 1, quantile, probs = 0.975, na.rm = T)
  plot_data <- data.frame(
    x = data$x,
    group = data$color,
    Fitted = rowMeans(simu_Onset),
    LowerCI = ci_lower,
    UpperCI = ci_upper
  )
  
  plot_data$date = as.Date('2021-07-01') + plot_data$x
  
  k = 30
  pd = 25
  
  
  plot_data$group[plot_data$group!=voc[4]] = 'A'
  plot_data = plot_data %>% group_by(group, date) %>% 
    summarise(Fitted = sum(Fitted),
              LowerCI = sum(LowerCI),
              UpperCI = sum(UpperCI))
  p = ggplot() +
    geom_point(data = dfob,
               aes(x = x, y = y/k,
                   group = group, color = group),
               alpha = 0.8) +
    geom_ribbon(data = plot_data,
                aes(x = date, group = group,
                    ymin = LowerCI, ymax = UpperCI, fill = group)) +
    geom_line(data = plot_data[plot_data$date < preddate,], 
              aes(x = date, y = Fitted, group = group, color= group)) +
    geom_line(data = plot_data[plot_data$date >= preddate,], 
              aes(x = date, y = Fitted, group = group, color= group),
              linetype = 'dashed') +
    scale_color_manual(name="",
                       values = alpha(values, 0.7)) +
    scale_fill_manual(name="",
                      values = alpha(values, 0.3)) +
    theme_bw() +
    theme(legend.position = "none",
          plot.title = element_text(hjust = 0.5)) +
    scale_x_date(breaks = seq(as.Date('2021-01-01'), as.Date('2025-05-01'), by="1 year"),
                 minor_breaks = seq(as.Date('2021-01-01'), as.Date('2025-05-01'), by ='1 month'),
                 date_labels = "%y-%b",
                 expand = c(0,0)) +
    xlab('') + ylab('Cases / Monthly average') + 
    coord_cartesian(xlim = c(d1, d2),
                    ylim = c(1,800))
  
  return(p)
}

p1 = getplot(simu_Onset1)

p2 = getplot(simu_Onset2, preddate = as.Date('2022-08-01'),
        d1 = as.Date('2021-08-01'), d2 = as.Date('2023-05-31'))

pdf(paste0("Output/S12F_flu_2020.pdf"), width = 2.1, height = 1.45)
print(p1)
print(p2)
dev.off()


