rm(list = ls())
library(rstan)
library(ggplot2)
library(scales)
library(ggnewscale)
library(tidyverse)
library(dplyr)

df = read.csv('influenza.csv')
observed_matrix = data.frame(v1 = df$A.H1*df$Total,
                             v2 = df$A.H3*df$Total,
                             v3 = df$B.Victoria*df$Total)
voc = c('A.H1','A.H3','B.Victoria')
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

fitlist = list()
load('flu.rdata')

poolday = 30
nday = 100
determinant_fun = function(cond = T, ifsimu  = T, n_simu = 1){
  # cond = F; ifsimu  = F; n_simu = 1
  n = 3
  seed_matrix <- matrix(0, nrow=n, ncol=poolday)
  seed_matrix[,1] = rep(3,3)
  if(ifsimu){
    seed_matrix[,1] = rmultinom(1, 9, rep(1/3,3))
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
      expected_total = unlist(round(observed_matrix[1,1:3])),
      seed_mat_I1 = seed_mats[[1]],
      seed_mat_I2 = seed_mats[[2]],
      seed_mat_I3 = seed_mats[[3]],
      seed_vec = seed_vec,
      gamma = 0.157
    )
    # Fit the model
    fit <- stan(file = 'FLU.stan', data = stan_data, 
                iter = 2500, chains = 1, warmup = 2000,
                verbose = TRUE)
    fitlist[[1]] = fit
  }
  
  
  Onsets_mat_list = list()
  fit= fitlist[[1]]
  
  posterior = rstan::extract(fit)
  
  pars_last = c(mean(posterior$contact), 
                mean(posterior$beta1),
                mean(posterior$beta2), 
                mean(posterior$beta3))
  
  if(ifsimu){
    pars_last = c(posterior$contact[n_simu], 
                  posterior$beta1[n_simu],
                  posterior$beta2[n_simu], 
                  posterior$beta3[n_simu])
  }
  Onsets_mat = simu(seed_mats, 
                    N = seed_vec * pars_last[1] + 1, 
                    poolday, pars = pars_last[-1], n)
  
  Onsets_mat_list[[1]] = Onsets_mat
  
  for (j in 1:34) {
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
      
      if(cond){
        stan_data <- list(
          poolday = poolday,
          nday = nday,
          expected_total = unlist(round(observed_matrix[j+1,1:3])),
          seed_mat_I1 = seed_mats[[1]],
          seed_mat_I2 = seed_mats[[2]],
          seed_mat_I3 = seed_mats[[3]],
          seed_vec = seed_vec,
          gamma = 0.157
        )
        # Fit the model
        fit <- stan(file = 'FLU.stan', data = stan_data, 
                    iter = 2500, chains = 1, warmup = 2000,
                    verbose = TRUE)
        fitlist[[j+1]] = fit
      }
      
      fit = fitlist[[j+1]]
      posterior = rstan::extract(fit)
      pars_last = c(mean(posterior$contact), 
                    mean(posterior$beta1),
                    mean(posterior$beta2),
                    mean(posterior$beta3))
      if(ifsimu){
        pars_last = c(posterior$contact[n_simu], 
                      posterior$beta1[n_simu],
                      posterior$beta2[n_simu], 
                      posterior$beta3[n_simu])
      }
      Onsets_mat = simu(seed_mats, 
                        N = seed_vec * pars_last[1] + 1, 
                        poolday, pars = pars_last[-1], n)
      
    }
    Onsets_mat_list[[j+1]] = Onsets_mat
  }
  
  # j = 35 last simulation -- 24-July
  
  return(Onsets_mat_list)
}

Onsets_mat_list = determinant_fun(cond = F, ifsimu = F, n_simu = 1)
n = 4
voc = c("A.H1", "A.H3", "B.Victoria", "A.H5") 
j = 35
Onsets_mat = Onsets_mat_list[[j]]
Onsets_mat_list[[j]] = cbind(Onsets_mat,Onsets_mat[,3]*0.005)

predictfun = function(beta){
  for (j in 35:47) {
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
      pars_last = c(mean(posterior$contact), 
                    mean(posterior$beta1),
                    mean(posterior$beta2),
                    mean(posterior$beta3),
                    beta)
      
      Onsets_mat = simu(seed_mats, 
                        N = seed_vec * pars_last[1] + 1, 
                        poolday, pars = pars_last[-1], n)
      
    }
    Onsets_mat_list[[j+1]] = Onsets_mat
  }
  
  dfplot_simu = data.frame()
  for (i in 36:48) {
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
  
  
  data$group = factor(data$color, levels = voc)
  
  data$date = as.Date('2021-07-01') + data$x
  
  return(data)
}

beta_vec = c(0.3,0.32,0.34,0.36)
data_list = list()
for (i in 1:length(beta_vec)) {
  data = predictfun(beta = beta_vec[i])
  data_list[[i]] = data
}

df2all = data.frame()
for (i in 1:length(beta_vec)) {
  data = data_list[[i]]
  df2 = data[data$color == voc[4],][-1,]
  df2all = rbind(df2all, df2)
}
df2all$beta = rep(beta_vec, each = nrow(df2all)/length(beta_vec))
df2all$beta = factor(df2all$beta, levels = beta_vec)
values = rev(c("#cc7722","#ffa500","#fedc56","#fff700"))
show_col(values)
p = ggplot() +
  geom_line(data = df2all, 
            aes(x = date, y = p, group = beta, color = beta)) +
  scale_color_manual(values = values,
                     name = expression(beta),
                     labels = c('0.30','0.32','0.34','0.36')) + 
  theme_bw() +
  theme(panel.grid = element_blank(),
        panel.background = element_blank(),
        plot.background = element_blank(),
        legend.position = 'right',
        legend.background = element_rect(color = NA, fill = NA),
        legend.key.size = unit(0.2,'cm')) + 
  scale_x_date(breaks = seq(as.Date('2021-01-01'), as.Date('2025-11-01'), by="6 months"),
               minor_breaks = seq(as.Date('2021-01-01'), as.Date('2025-11-01'), by ='1 month'),
               date_labels = "%y-%b",
               expand = c(0,0)) +
  scale_y_continuous(breaks = c(0,0.5,1)) +
  xlab('') + ylab('') + 
  coord_cartesian(xlim = c(as.Date('2024-06-30'), as.Date('2025-05-31')),
                  ylim = c(0,1))
pdf(paste0("Output/flu_H5.pdf"), width = 2.1, height = 1)
print(p)
dev.off()


ggplot() +
  geom_line(data = data, 
            aes(x = date, y = y, group = group, color= group)) +
  coord_cartesian(ylim = c(2,max(data$y))) +
  labs(x = "Date", y = "Proportion") +
  theme_bw() +
  theme(legend.position = "right",
        plot.title = element_text(hjust = 0.5)) +
  scale_x_date(breaks = seq(as.Date('2021-01-01'), as.Date('2025-11-01'), by="6 months"),
               minor_breaks = seq(as.Date('2021-01-01'), as.Date('2025-11-01'), by ='1 month'),
               date_labels = "%y-%b",
               expand = c(0,0)) +
  xlab('') + ylab('Cases') + 
  coord_cartesian(xlim = c(as.Date('2024-06-30'), as.Date('2025-05-31')),
                  ylim = c(1,max(data$y)))

values2 = c('#98afc7','#0041c2', '#a37ca1', values[2])

plot_data2 = data_list[[4]]
plist = list()
for (i in 1:4) {
  plot_data2 = data_list[[i]]
  p = ggplot() +
    geom_area(data = plot_data2, 
              aes(x = date, y = p, fill = group),
              position = 'fill') +
    scale_fill_manual(name="", 
                      breaks = voc,
                      values = alpha(values2, 1)) +
    coord_cartesian(xlim = c(as.Date('2024-06-30'), as.Date('2025-05-31')),
                    ylim = c(0,1)) + 
    scale_x_date(breaks = seq(as.Date('2021-01-01'), as.Date('2025-05-01'), by="6 months"),
                 minor_breaks = seq(as.Date('2021-01-01'), as.Date('2025-05-01'), by ='1 month'),
                 date_labels = "%y-%b",
                 expand = c(0,0)) +
    scale_y_continuous(breaks = seq(0,1,0.5), expand = c(0, 0)) +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_rect(color = NA, fill = NA),
          legend.position = 'right',
          legend.background = element_rect(color = NA, fill = NA),
          legend.key = element_blank(),
          legend.key.size = unit(0.4, units = 'cm')) +
    xlab('') + ylab('Prevalence') 
  plist[[i]] = p
}

p2
pdf(paste0("Output/flu_H5_prevalence.pdf"), width = 2.5, height = 1.2)
for (i in 1:4) {
  print(plist[[i]])
}
dev.off()

parsmat = data.frame()
for (i in 1:length(fitlist)) {
  posterior = rstan::extract(fitlist[[i]])
  parsmat = rbind(parsmat, 
                  data.frame(Tcycle = i,
                             y = c(posterior$contact, 
                                   posterior$beta1,
                                   posterior$beta2,
                                   posterior$beta3),
                             group = c(rep('contact', length(posterior$contact)),
                                       rep('beta1', length(posterior$beta1)),
                                       rep('beta2', length(posterior$beta2)),
                                       rep('beta3', length(posterior$beta3)))
                             
                  ))
  
}

poolday = 30
parsmat$date = poolday*(parsmat$Tcycle-1) + as.Date('2021-07-01') 
df = parsmat[parsmat$group != 'contact', ]
df$date = as.Date(df$date)
pd = 8

medians <- df %>%
  group_by(date, group) %>%
  summarize(median_y = median(y))
p3 = ggplot() +
  geom_boxplot(data = df[df$group == 'beta1',], 
               aes(x = date-pd, y = y, group = date),
               color = alpha(values[1],0.7), 
               fill = alpha(values[1],0.3),
               width = 6,
               outlier.shape = NA) +
  geom_boxplot(data = df[df$group == 'beta2',], 
               aes(x = date, y = y, group = date),
               color = alpha(values[2],0.7), 
               fill = alpha(values[2],0.3),
               width = 6,
               outlier.shape = NA) +
  geom_boxplot(data = df[df$group == 'beta3',], 
               aes(x = date+pd, y = y, group = date),
               color = alpha(values[3],0.7), 
               fill = alpha(values[3],0.3),
               width = 6,
               outlier.shape = NA) +
  scale_x_date(breaks = seq(as.Date('2020-01-01'), as.Date('2024-11-01'), by="6 months"),
               minor_breaks = seq(as.Date('2019-12-01'), as.Date('2024-11-01'), by ='1 month'),
               date_labels = "%y-%b",
               expand = c(0, 0)) +
  geom_line(data = medians[medians$group == 'beta1',], 
            aes(x = date - pd, y = median_y), 
            color = values[1], 
            linewidth = 0.4) +
  geom_line(data = medians[medians$group == 'beta2',], 
            aes(x = date, y = median_y), 
            color = values[2], 
            linewidth = 0.4) +
  geom_line(data = medians[medians$group == 'beta3',], 
            aes(x = date + pd, y = median_y), 
            color = values[3], 
            linewidth = 0.4) +
  coord_cartesian(xlim = c(as.Date('2021-09-01'), 
                           as.Date('2024-05-31')),
                  ylim = c(0.24,0.37)) +
  scale_y_continuous(breaks = c(0.25,0.3,0.35)) +
  theme_bw() + xlab('') + ylab(expression(beta))
pdf(paste0("Output/flu_plot_beta.pdf"), width = 3.3, height = 1.2)
print(p3)
dev.off()



if(F){
  data = determinant_fun(cond = F, ifsimu  = T, n_simu = n_simu)
  data$date = as.Date('2021-07-01') + data$x
  dfob = data.frame(y = c(observed_matrix[,1],
                          observed_matrix[,2],
                          observed_matrix[,3]),
                    x = rep(as.Date(rownames(observed_matrix)),
                            3),
                    group = rep(voc, each = nrow(observed_matrix)))
  
  tmp = sapply(Onsets_mat_list, function(x){
    colSums(x)
  })
  dfsimu = data.frame(y = c(tmp[1,],tmp[2,],tmp[3,]),
                      group = rep(voc, each = length(Onsets_mat_list)),
                      x = rep(dfob$x[1:length(Onsets_mat_list)], 3))
  for (i in 1:4) {
    dfsimu = rbind(dfsimu, dfsimu)
  }
  dfsimu$y2 = dfsimu$y + runif(nrow(dfsimu), min = -5, max = 1000)
  
  dfsimu2 = dfsimu %>% group_by(x, group) %>% 
    summarise(y = mean(y2),
              y1 = quantile(y2, 0.025),
              y2 = quantile(y2, 0.975)) %>%
    as.data.frame()
  dfsimu2$group = factor(dfsimu2$group, levels = voc)
  ggplot() +
    geom_col(data = dfob,
             aes(x = x, y = y/30, 
                 group = group, fill = group),
             alpha = 0.8,
             position = position_dodge(width = 20)) +
    geom_point(data = dfsimu2,
               aes(x = x, y = y/30, color = group),
               alpha = 0.5,
               position = position_dodge(width = 20))+ 
    geom_errorbar(data = dfsimu2,
                  aes(x = x,
                      ymin = y1/30, ymax = y2/30,
                      color = group),
                  position = position_dodge(width = 20))+
    geom_line(data = data, 
              aes(x = date, y = y, group = group, color= group)) +
    coord_cartesian(ylim = c(2,max(data$y))) +
    labs(x = "Date", y = "Proportion") +
    theme_minimal() +
    theme(legend.position = "right",
          plot.title = element_text(hjust = 0.5)) +
    scale_x_date(breaks = seq(as.Date('2021-01-01'), as.Date('2024-11-01'), by="6 months"),
                 minor_breaks = seq(as.Date('2021-01-01'), as.Date('2024-11-01'), by ='1 month'),
                 date_labels = "%y-%b") +
    xlab('') + ylab('Cases') + 
    coord_cartesian(xlim = c(as.Date('2021-06-01'), as.Date('2024-05-31')),
                    ylim = c(1,max(data$y)))
  
  data$group = factor(data$group, levels = c('Delta','Alpha','D614G'))
  
  
  ggplot() +
    geom_area(data = data, 
              aes(x = date, y = p, fill = group),
              position = 'fill') +
    coord_cartesian(xlim = c(as.Date('2020-06-30'), as.Date('2021-10-31')),
                    ylim = c(0,1)) +
    scale_x_date(breaks = seq(as.Date('2020-01-01'), as.Date('2024-11-01'), by="6 months"),
                 minor_breaks = seq(as.Date('2019-12-01'), as.Date('2024-11-01'), by ='1 month'),
                 date_labels = "%y-%b",
                 expand = c(0, 0)) +
    scale_y_continuous(breaks = seq(0,1,0.5), expand = c(0, 0))
  
  parsmat = data.frame()
  for (i in 1:length(fitlist)) {
    posterior = rstan::extract(fitlist[[i]])
    parsmat = rbind(parsmat, 
                    data.frame(Tcycle = i,
                               y = c(posterior$contact, 
                                     posterior$beta1,
                                     posterior$beta2,
                                     posterior$beta3),
                               group = c(rep('contact', length(posterior$contact)),
                                         rep('beta1', length(posterior$beta1)),
                                         rep('beta2', length(posterior$beta2)),
                                         rep('beta3', length(posterior$beta3)))
                               
                    ))
    
  }
  parsmat$group = factor(parsmat$group)
  parsmat$Tcycle = factor(parsmat$Tcycle)
  ggplot(parsmat[parsmat$group != 'contact', ], 
         aes(x = Tcycle, y = y, color = group, fill = group)) +
    geom_boxplot(outlier.shape = NA) 
  
}


