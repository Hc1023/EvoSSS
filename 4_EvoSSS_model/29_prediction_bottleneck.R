rm(list = ls())
library(rstan)
library(ggplot2)
library(scales)
library(ggnewscale)
library(tidyverse)
library(dplyr)
load('AB.rdata')
df = read.csv('../3_Epidemiological_analysis/Covid19CasesGISAID.csv')

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

update_fun = function(pars, states_old, N){
  
  S <- states_old[,1]
  I1 <- states_old[,2]
  I2 <- states_old[,3]
  
  beta1 = pars[1]
  beta2 = pars[2]
  gamma = 0.157
  
  S_new = S - beta1*S*I1/N - beta2*S*I2/N
  I1_new = I1 + beta1*S*I1/N - gamma*I1
  I2_new = I2 + beta2*S*I2/N - gamma*I2
  Onset1 = beta1*S*I1/N
  Onset2 = beta2*S*I2/N
  
  return(cbind(S_new, I1_new, I2_new, Onset1, Onset2))
}

update_fun_stochastic = function(pars, states_old, N){
  
  beta1 = pars[1]
  beta2 = pars[2]
  gamma = 0.157
  
  mat = data.frame()
  for (h in 1:nrow(states_old)) {
    S <- states_old[h,1]
    I1 <- states_old[h,2]
    I2 <- states_old[h,3]
    
    pS_vec = c(max(0,beta1*I1/N[h]), max(0,beta2*I2/N[h]), 
               max(0,1-beta1*I1/N[h]-beta2*I2/N[h]))
    sample_S <- rmultinom(1, size = S, prob = pS_vec)
    pI_vec <- c(gamma, 1-gamma)
    sample_I1 <- rmultinom(1, size = I1, prob = pI_vec)
    sample_I2 <- rmultinom(1, size = I2, prob = pI_vec)
    
    ## new values
    S_new <- sample_S[3]
    I1_new <- sample_I1[2] + sample_S[1] 
    I2_new <- sample_I2[2] + sample_S[2] 
    
    Onset1 <- sample_S[1]
    Onset2 <- sample_S[2]
    mat[h,1:5] = c(S_new, I1_new, I2_new, Onset1, Onset2)
  }
  
  return(mat)
}

simu <- function(seed_mat_I1, seed_mat_I2, N, poolday,
                 pars = c(0.379, 0.398, 0.157),
                 f = update_fun) {
  
  # Initial conditions
  I1_old = seed_mat_I1[1,]
  I2_old = seed_mat_I2[1,]
  S_old = N - I1_old - I2_old
  states_old = cbind(S_old, I1_old, I2_old)
  
  # Simulate the dynamics over 
  ndays = 200
  
  Onsets_mat <- matrix(0, ndays, 2)
  
  for (t in 2:ndays) {
    states_old = f(pars = pars, states_old = states_old, N = N)
    if(t <= nrow(seed_mat_I1)){
      states_old[t,2:3] = states_old[t,2:3] + 
        c(seed_mat_I1[t,t], seed_mat_I2[t,t])
    }
    Onsets_mat[t,] = c(sum(states_old[,4]), sum(states_old[,5]))
  }
  
  return(Onsets_mat)
}


determinant_fun = function(ifsimu  = T, n_simu = 1,
                           m, h, k_capacity, f){
  
  poolday = 30
  nday = 100
  n = 2
  seed_matrix <- matrix(0, nrow=2, ncol=poolday)
  seed_matrix[,1] = 6.6*c(0.4,0.6) + 1e-3
  seed_vec = colSums(seed_matrix)
  seed_mats <- list()
  
  for (i in 1:2) {
    seed_mats[[i]] <- diag(seed_matrix[i,])
  }
  
  if(!ifsimu){
    expected_matrix = observed_matrix[1:nday,]
    expected_matrix[expected_matrix < 0] = 0
    expected_matrix = round(expected_matrix)
    expected_matrix[(poolday*2+1):nday,] = 0
    fexpect = data.frame(x = rep(1:nday,2), 
                         y = c(expected_matrix[1:nday,1],
                               expected_matrix[1:nday,2]),
                         group = factor(rep(c('Non Omicron','Omicron'), 
                                            each = nday),
                                        levels = c('Non Omicron','Omicron')))
    
    
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
  
  Onsets_mat = simu(seed_mat_I1 = seed_mats[[1]], 
                    seed_mat_I2 = seed_mats[[2]], 
                    N = seed_vec * pars_last[1] + 1, 
                    poolday = poolday, 
                    pars = pars_last[-1],
                    f = update_fun)
  if(!ifsimu){
    fonset = data.frame(x = rep(1:nday,2), 
                        y = c(Onsets_mat[1:nday,1],
                              Onsets_mat[1:nday,2]),
                        group = factor(rep(c('Non Omicron','Omicron'), 
                                           each = nday),
                                       levels = c('Non Omicron','Omicron')))
    
    
    ggplot() +
      geom_point(data = fexpect, 
                 aes(x = x, y = y, group = group, color = group)) +
      geom_line(data = fonset,
                aes(x = x, y = y, group = group, color = group))
    
  }
  
  
  Onsets_mat_list[[1]] = Onsets_mat
  
  thresh_cycle = 2
  all_cycle = thresh_cycle + 5
  for (j in 1:all_cycle) {
    if(j == thresh_cycle+1){
      {
        
        Onsets_mat = Onsets_mat_list[[j]]
        
        # Generalized extraction of Onset columns for all variants
        Onsets <- list()
        for (i in 1:n) {
          Onsets[[i]] <- Onsets_mat[poolday + 1:h, i] + 1e-3
        }
        # Mobility control force
        mobility <- rep(m, h) 
        Mobility_matrix <- diag(mobility)
        
        # Create seed_vec considering all variants
        seed_vec <- rowSums(sapply(Onsets, function(Onset) Onset)) %*% Mobility_matrix %>% as.numeric()
        
        # Calculate probabilities for each variant
        probs <- lapply(Onsets, function(Onset) Onset / rowSums(do.call(cbind, Onsets)))
        seed_matrix <- matrix(0, nrow=n, ncol=length(seed_vec))
        
        
        if(ifsimu){
          seed_matrix = sapply(1:length(seed_vec), function(x){
            rmultinom(1, ceiling(seed_vec[x]), c(probs[[1]][x],1-probs[[1]][x]))
          })
        }else{
          for (i in 1:n) {
            seed_matrix[i,] <- seed_vec * probs[[i]]
          }
        }
        
        # Create seed_matrices for each variant
        seed_mats <- list()
        
        for (i in 1:n) {
          seed_mats[[i]] <- diag(seed_matrix[i,])
        }
        
        if(!ifsimu){
          expected_matrix = observed_matrix[poolday*j+1:nday,] - Onsets_mat[poolday+1:nday,] 
          expected_matrix[expected_matrix<0] = 0
          expected_matrix = round(expected_matrix)
          
          expected_matrix[(2*poolday + 1):nday,] = 0
          fexpect = data.frame(x = rep(1:nday,2), 
                               y = c(expected_matrix[1:nday,1],
                                     expected_matrix[1:nday,2]),
                               group = factor(rep(c('Non Omicron','Omicron'), 
                                                  each = nday),
                                              levels = c('Non Omicron','Omicron')))
          
        }
        
        fit = fitlist[[j]]
        posterior = rstan::extract(fit)
        
        pars_last = c(mean(posterior$contact)*k_capacity, 
                      mean(posterior$beta1),
                      mean(posterior$beta2))
        if(ifsimu){
          pars_last = c(posterior$contact[n_simu]*k_capacity, 
                        posterior$beta1[n_simu],
                        posterior$beta2[n_simu])
        }
        
        
        Onsets_mat = simu(seed_mats[[1]], seed_mats[[2]],
                          N = seed_vec * pars_last[1] + 1, 
                          poolday, pars = pars_last[-1],
                          f = f)
        if(!ifsimu){
          fonset = data.frame(x = rep(1:nday,2), 
                              y = c(Onsets_mat[1:nday,1],
                                    Onsets_mat[1:nday,2]),
                              group = factor(rep(c('Non Omicron','Omicron'), 
                                                 each = nday),
                                             levels = c('Non Omicron','Omicron')))
          
          ggplot() +
            geom_point(data = fexpect, 
                       aes(x = x, y = y, group = group, color = group)) +
            geom_line(data = fonset,
                      aes(x = x, y = y, group = group, color = group)) 
          
          
        }
      }
    }else{
      
      
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
      
      
      if(ifsimu){
        seed_matrix = sapply(1:length(seed_vec), function(x){
          rmultinom(1, ceiling(seed_vec[x]), c(probs[[1]][x],1-probs[[1]][x]))
        })
      }else{
        for (i in 1:n) {
          seed_matrix[i,] <- seed_vec * probs[[i]]
        }
      }
      
      # Create seed_matrices for each variant
      seed_mats <- list()
      
      for (i in 1:n) {
        seed_mats[[i]] <- diag(seed_matrix[i,])
      }
      
      if(!ifsimu){
        expected_matrix = observed_matrix[poolday*j+1:nday,] - Onsets_mat[poolday+1:nday,] 
        expected_matrix[expected_matrix<0] = 0
        expected_matrix = round(expected_matrix)
        expected_matrix[(2*poolday + 1):nday,] = 0
        fexpect = data.frame(x = rep(1:nday,2), 
                             y = c(expected_matrix[1:nday,1],
                                   expected_matrix[1:nday,2]),
                             group = factor(rep(c('Non Omicron','Omicron'), 
                                                each = nday),
                                            levels = c('Non Omicron','Omicron')))
        
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
      Onsets_mat = simu(seed_mats[[1]], seed_mats[[2]],
                        N = seed_vec * pars_last[1] + 1, 
                        poolday, pars = pars_last[-1],
                        f = update_fun)
      
      fonset = data.frame(x = rep(1:nday,2), 
                          y = c(Onsets_mat[1:nday,1],
                                Onsets_mat[1:nday,2]),
                          group = factor(rep(c('Non Omicron','Omicron'), 
                                             each = nday),
                                         levels = c('Non Omicron','Omicron')))
      
      ggplot() +
        geom_point(data = fexpect, 
                   aes(x = x, y = y, group = group, color = group)) +
        geom_line(data = fonset,
                  aes(x = x, y = y, group = group, color = group)) 
      
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

fexpect0 = data.frame(y = c(observed_matrix$v1,observed_matrix$v2),
                      x = rep(1:nrow(observed_matrix),2),
                      group = rep(c('A', 'B'), 
                                  each = nrow(observed_matrix)))
fexpect0$date = as.Date('2019-12-31') + fexpect0$x
fexpect0$group = factor(fexpect0$group, levels = c('A','B'))

getplot = function(m, h, n, k_capacity, f){
  df2_list = list()
  for (n_simu in 1:n) {
    print(n_simu)
    data = determinant_fun(ifsimu  = T, 
                           n_simu = n_simu,
                           m = m, h = h, 
                           k_capacity = k_capacity, 
                           f = f)
    
    df2_list[[n_simu]] = data$y
  }
  simu_Onset = data.frame(bind_cols(df2_list))
  ci_lower <- apply(simu_Onset, 1, quantile, probs = 0.025, na.rm = T)
  ci_upper <- apply(simu_Onset, 1, quantile, probs = 0.975, na.rm = T)
  plot_data <- data.frame(
    x = data$x,
    V = data$color,
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

plotfun2 = function(plot_data1){
  plot_data = plot_data1[plot_data1$date <= as.Date('2020-01-01') + 90, ]
  plot_data_pred = plot_data1[plot_data1$date >= as.Date('2020-01-01') + 90, ]
  values = c(hue_pal()(3)[1], hue_pal()(3)[3])
  values2 = c('#a70107', '#06068d')
  values2 = values
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
  p
  return(p)
}

plot_data1 = getplot(m = 1/200, h = 2, n = 500, 
                     k_capacity = 20, 
                     f = update_fun_stochastic)
plot_data2 = getplot(m = 1/1000, h = 2, n = 500, 
                     k_capacity = 100, 
                     f = update_fun_stochastic)
plot_data3 = getplot(m = 1/2000, h = 2, n = 500, 
                     k_capacity = 200, 
                     f = update_fun_stochastic)

save(plot_data1, plot_data2, plot_data3,
     file = 'mobility.rdata')
load(file = 'mobility.rdata')
p1 = plotfun2(plot_data1 = plot_data1)
p1
p2 = plotfun2(plot_data1 = plot_data2)
p2
p3 = plotfun2(plot_data1 = plot_data3)
p3
pdf(paste0("Output/prediction_mobility.pdf"), 
    width = 2, height = 1.8)
print(p1)
print(p2)
print(p3)
dev.off()



data = determinant_fun(ifsimu  = F, n_simu = n_simu,
                       m = 1/30, h = 30, k_capacity = 1, f = update_fun)
data = determinant_fun(ifsimu  = F, n_simu = n_simu,
                       m = 1/30, h = 30, k_capacity = 0.5, f = update_fun)
data = determinant_fun(ifsimu  = F, n_simu = n_simu,
                       m = 1/30, h = 30, k_capacity = 0.9, f = update_fun)

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
  coord_cartesian(xlim = c(as.Date('2020-01-01'), as.Date('2020-01-01')+30*7),
                  ylim = c(1,2*10^3))



