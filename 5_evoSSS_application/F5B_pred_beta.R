rm(list = ls())
library(rstan)
library(ggplot2)
library(scales)
library(ggnewscale)
library(tidyverse)
library(dplyr)

df = read.csv('../3_Epidemiological_analysis/F3B_Covid19CasesGISAID.csv')
load('F5B_AB.rdata')

VOC = c("B", "A")
names(VOC) = c("Lineage B","Lineage A")
df = df[df$Mutations %in% names(VOC),]
df$V = NA
for (i in 1:length(VOC)) {
  df$V[df$Mutations == names(VOC)[i]] = VOC[i]
}
df$V = factor(df$V, levels = VOC[1:2])
df = df %>% group_by(Var1, V) %>%
  summarise(y = sum(Freq))
colnames(df)[1] = 'date'
df <- na.omit(df)
max(df$date)
as.Date('2019-12-31') + 700
full_dates <- as.Date('2019-12-31') + 1:800
# Create a dataframe with all dates
full_df <- expand.grid(date = full_dates, 
                       V = unique(df$V))
df$date = as.Date(df$date)
merged_df <- full_df %>%
  left_join(df, by = c("date", "V")) %>%
  mutate(count = ifelse(is.na(y), 0, y))
voc = unique(merged_df$V)

observed_matrix = data.frame(v1 = merged_df[merged_df$V == 'A', 'count'],
                             v2 = merged_df[merged_df$V == 'B', 'count'])
rownames(observed_matrix) = full_dates

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

determinant_fun = function(cond = T, ifsimu  = T, n_simu = 1, kb){
  
  n = 2
  poolday = 30
  nday = 100
  
  seed_matrix <- matrix(0, nrow=2, ncol=poolday)
  seed_matrix[,1] = 6.6*c(0.4,0.6) + 1e-3
  seed_vec = colSums(seed_matrix)
  seed_mats <- list()
  for (i in 1:2) {
    seed_mats[[i]] <- diag(seed_matrix[i,])
  }

  Onsets_mat_list = list()
  fit= fitlist[[1]]
  
  posterior = rstan::extract(fit)
  
  pars_last = c(mean(posterior$contact), 
                mean(posterior$beta1),
                mean(posterior$beta2))
  if(ifsimu){
    pars_last = c(posterior$contact[n_simu], 
                  posterior$beta1[n_simu],
                  posterior$beta2[n_simu])
  }
  Onsets_mat = simu(seed_mats, 
                    N = seed_vec * pars_last[1] + 1, 
                    poolday, pars = pars_last[-1], n)
  Onsets_mat_list[[1]] = Onsets_mat

  thresh_cycle = 2
  all_cycle = thresh_cycle + 5
  for (j in 1:thresh_cycle) {
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
      
      fit = fitlist[[j+1]]
      posterior = rstan::extract(fit)
      
      pars_last = c(mean(posterior$contact), 
                    mean(posterior$beta1),
                    mean(posterior$beta2))
      if(ifsimu){
        pars_last = c(posterior$contact[n_simu], 
                      posterior$beta1[n_simu],
                      posterior$beta2[n_simu])
      }
      Onsets_mat = simu(seed_mats, 
                        N = seed_vec * pars_last[1] + 1, 
                        poolday, pars = pars_last[-1], n)
      
    }
    
    Onsets_mat_list[[j+1]] = Onsets_mat
  }
  
  for (j in (thresh_cycle+1):all_cycle) {
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
      seed_matrix = sapply(1:length(seed_vec), function(x){
        rmultinom(1, ceiling(seed_vec[x]), c(probs[[1]][x],1-probs[[1]][x]))
      })
      
      # Create seed_matrices for each variant
      seed_mats <- list()
      
      for (i in 1:n) {
        seed_mats[[i]] <- diag(seed_matrix[i,])
      }
      # 
      # fit = fitlist[[j+1]]
      # posterior = rstan::extract(fit)
      
      pars_last = c(mean(posterior$contact), 
                    mean(posterior$beta2)*kb,
                    mean(posterior$beta2))
      if(ifsimu){
        pars_last = c(posterior$contact[n_simu], 
                      posterior$beta2[n_simu]*kb,
                      posterior$beta2[n_simu])
      }
      Onsets_mat = simu(seed_mats, 
                        N = seed_vec * pars_last[1] + 1, 
                        poolday, pars = pars_last[-1], n)
      
    }
    
    Onsets_mat_list[[j+1]] = Onsets_mat
  }
  
  dfplot_simu = data.frame()
  for (i in 1:(all_cycle+1)) {
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

if(F) fitlist = list()
# save(fitlist, file = 'AB.rdata')

getplot = function(kb, n){
  df2_list = list()
  for (n_simu in 1:n) {
    print(n_simu)
    data = determinant_fun(cond = F, ifsimu  = T, 
                           n_simu = n_simu, kb = kb)
    df2_list[[n_simu]] = data$y
  }
  simu_Onset = data.frame(bind_cols(df2_list))
  ci_lower <- apply(simu_Onset, 1, quantile, probs = 0.025, na.rm = T)
  ci_upper <- apply(simu_Onset, 1, quantile, probs = 0.975, na.rm = T)
  plot_data <- data.frame(
    x = data$x,
    V = data$color,
    # Observed = observed_cases,
    Fitted = rowMeans(simu_Onset),
    LowerCI = ci_lower,
    UpperCI = ci_upper
  )
  plot_data$date = plot_data$x + as.Date('2019-12-31')
  plot_data$group = 'A'
  plot_data$group[data$color == '2'] = 'B'
  plot_data$group = factor(plot_data$group, levels = c('A','B'))
  plot_data = plot_data[plot_data$x>1,]
  
  return(plot_data)
}

plotfun = function(plot_data1){
  plot_data = plot_data1[plot_data1$date <= as.Date('2020-01-01') + 90, ]
  plot_data_pred = plot_data1[plot_data1$date >= as.Date('2020-01-01') + 90, ]
  values = c(hue_pal()(3)[1], hue_pal()(3)[3])
  values2 = c('#d43f3b','#0041c2')
  # maxy = max(max(fexpect0$y[fexpect0$date < as.Date('2020-01-01') + 8*30]), 
  #            max(plot_data1$UpperCI[plot_data1$date <  as.Date('2020-01-01') + 8*30]))
  p = ggplot() +
    scale_y_continuous(trans='log10') +
    coord_cartesian(xlim = c(as.Date('2020-01-01'), as.Date('2020-01-01') + 8*30),
                    ylim = c(1, 8*10^3)) +
    geom_point(data = fexpect0, 
               aes(x = date, y = y, 
                   group = group, color = group),
               size = 0.4, shape = 16) +
    geom_ribbon(data = plot_data, 
                aes(x = date, group = group, 
                    ymin = LowerCI, ymax = UpperCI, 
                    fill = group)) +  # Confidence interval
    geom_line(data = plot_data, 
              aes(x = date, y = Fitted, 
                  group = group, color = group), 
              linewidth = 1) +
    scale_color_manual(name="",
                       values = alpha(values, 0.7),
                       breaks = names(values)) +
    scale_fill_manual(name="",
                      values = alpha(values, 0.3),
                      breaks = names(values)) +
    new_scale_color() +
    new_scale_fill() +
    geom_ribbon(data = plot_data_pred, 
                aes(x = date, group = group, 
                    ymin = LowerCI, ymax = UpperCI, 
                    fill = group)) + 
    geom_line(data = plot_data_pred, 
              aes(x = date, y = Fitted, 
                  group = group, color = group),
              linetype = 'dashed',
              linewidth = 1) +
    scale_color_manual(name="",
                       values = alpha(values2, 0.7),
                       breaks = names(values2)) +
    scale_fill_manual(name="",
                      values = alpha(values2, 0.3),
                      breaks = names(values2)) +
    labs(x = "Date", y = "Proportion") +
    theme_bw() +
    theme(legend.position = "right",
          plot.title = element_text(hjust = 0.5)) +
    scale_x_date(breaks = seq(as.Date('2020-01-01'), as.Date('2022-11-01'), by="3 months"),
                 minor_breaks = seq(as.Date('2019-12-01'), as.Date('2022-05-01'), by ='1 month'),
                 date_labels = "%y-%b", expand = c(0, 0)) +
    scale_y_continuous(trans='log10', 
                       breaks = c(1, 10, 100, 1000, 10000),
                       labels = c(expression(10^0), expression(10^1),
                                  expression(10^2), expression(10^3),
                                  expression(10^4))) +
    xlab('') + ylab('Cases') + 
    theme(legend.position = 'right',
          legend.background = element_rect(color = NA, fill = NA),
          legend.key = element_blank(),
          legend.key.size = unit(0.2, units = 'cm'),
          legend.key.width = unit(1, units = 'cm'))
  return(p)
}

fexpect0 = data.frame(y = c(observed_matrix$v1,observed_matrix$v2),
                      x = rep(1:nrow(observed_matrix),2),
                      group = rep(c('A', 'B'), 
                                  each = nrow(observed_matrix)))
fexpect0$date = as.Date('2019-12-31') + fexpect0$x
fexpect0$group = factor(fexpect0$group, levels = c('A','B'))

if(F){
  plot_data1 = getplot(kb = 1, n = 500)
  plot_data2 = getplot(kb = 1.1, n = 500)
  plot_data3 = getplot(kb = 0.9, n = 500)
  save(plot_data1, plot_data2, plot_data3, 
       file = 'F5B_pred_beta.rdata')
}
load('F5B_pred_beta.rdata')
p1 = plotfun(plot_data1)
p2 = plotfun(plot_data2)
p3 = plotfun(plot_data3)

pdf(paste0("Output/F5B_pred_beta.pdf"), width = 2, height = 1.8)
print(p1)
print(p2)
print(p3)
dev.off()


if(F){
  plot_data = plot_data1[plot_data1$date <= as.Date('2020-01-01') + 90, ]
  plot_data_pred = rbind(plot_data1[plot_data1$date >= as.Date('2020-01-01') + 90, ],
                         plot_data2[plot_data2$date >= as.Date('2020-01-01') + 90, ],
                         plot_data3[plot_data3$date >= as.Date('2020-01-01') + 90, ])
  plot_data_pred$group_line = rep(c(1,2,3), each = nrow(plot_data_pred)/3)
  plot_data_pred$group_pred = paste0(plot_data_pred$group,
                                     rep(c(1,2,3), each = nrow(plot_data_pred)/3))
  values = c(hue_pal()(3)[1], hue_pal()(3)[3],
             '#a70107', '#06068d',
             '#f8766d', '#619cff',
             '#ff2020', '#2020ff')
  values = c(hue_pal()(3)[1], hue_pal()(3)[3],
             '#ff2020', '#2020ff',
             '#f8766d', '#619cff',
             '#a70107', '#06068d')
  names(values) = c('A','B',
                    'A1','B1',
                    'A2','B2',
                    'A3','B3')
  
  p = ggplot() +
    geom_point(data = fexpect0, 
               aes(x = date, y = y, 
                   group = group, color = group),
               size = 0.4, shape = 16) +
    geom_ribbon(data = plot_data, 
                aes(x = date, group = group, 
                    ymin = LowerCI, ymax = UpperCI, fill = group)) +  # Confidence interval
    geom_line(data = plot_data, 
              aes(x = date, y = Fitted, 
                  group = group, color = group), linewidth = 1) +
    geom_ribbon(data = plot_data_pred, 
                aes(x = date, group = group_pred, 
                    ymin = LowerCI, ymax = UpperCI, fill = group_pred)) + 
    geom_line(data = plot_data_pred, 
              aes(x = date, y = Fitted, 
                  group = group_pred, color = group_pred,
                  linetype = factor(group_line)), 
              linewidth = 1) +
    scale_color_manual(name="",
                       values = alpha(values, 0.7),
                       breaks = names(values)) +
    scale_fill_manual(name="",
                      values = alpha(values, 0.3),
                      breaks = names(values)) +
    scale_linetype_manual(values = c('solid','dotted','dashed')) + 
    scale_y_continuous(trans='log10') +
    coord_cartesian(ylim = c(2, max(plot_data$Fitted))) +
    labs(x = "Date", y = "Proportion") +
    theme_bw() +
    theme(legend.position = "right",
          plot.title = element_text(hjust = 0.5)) +
    scale_x_date(breaks = seq(as.Date('2020-01-01'), as.Date('2022-11-01'), by="6 months"),
                 minor_breaks = seq(as.Date('2019-12-01'), as.Date('2022-05-01'), by ='1 month'),
                 date_labels = "%y-%b", expand = c(0, 0)) +
    scale_y_continuous(trans='log10', 
                       breaks = c(1, 10, 100, 1000, 10000),
                       labels = c(expression(10^0), expression(10^1),
                                  expression(10^2), expression(10^3),
                                  expression(10^4))) +
    xlab('') + ylab('Cases') + 
    coord_cartesian(xlim = c(as.Date('2020-01-01'), as.Date('2021-10-31')),
                    ylim = c(1, 4*10^4)) +
    theme(legend.position = 'right',
          legend.background = element_rect(color = NA, fill = NA),
          legend.key = element_blank(),
          legend.key.size = unit(0.2, units = 'cm'),
          legend.key.width = unit(1, units = 'cm'))
  
  
}






if(F){
  data = determinant_fun(cond = F, ifsimu  = F, n_simu = n_simu, kb = 1)
  data = determinant_fun(cond = F, ifsimu  = F, n_simu = n_simu, kb = 1.1)
  data$date = as.Date('2019-12-31') + data$x
  data$group = 'A'
  data$group[data$color == '2'] = 'B'
  data$group = factor(data$group, levels = c('A','B'))
  fexpect0 = data.frame(y = c(observed_matrix$v1,observed_matrix$v2,
                              observed_matrix$v3),
                        x = rep(1:nrow(observed_matrix),2),
                        group = rep(c('A','B'), 
                                    each = nrow(observed_matrix)))
  fexpect0$date = as.Date('2019-12-31') + fexpect0$x
  fexpect0$group = factor(fexpect0$group, levels = c('A','B'))
  
  
  ggplot() +
    geom_point(data = fexpect0, 
               aes(x = date, y = y, 
                   group = group, color = group)) +
    geom_line(data = data, 
              aes(x = date, y = y, group = group, color= group)) +
    scale_y_continuous(trans='log10') +
    coord_cartesian(ylim = c(2,max(data$y))) +
    labs(x = "Date", y = "Proportion") +
    theme_minimal() +
    theme(legend.position = "right",
          plot.title = element_text(hjust = 0.5)) +
    scale_x_date(breaks = seq(as.Date('2020-01-01'), as.Date('2022-11-01'), by="6 months"),
                 minor_breaks = seq(as.Date('2019-12-01'), as.Date('2022-11-01'), by ='1 month'),
                 date_labels = "%y-%b") +
    xlab('') + ylab('Cases') + 
    coord_cartesian(xlim = c(as.Date('2020-01-01'), as.Date('2021-10-31')),
                    ylim = c(1,2*10^4))
  

  library(ggplot2)
  library(dplyr)
  
  
  
}
