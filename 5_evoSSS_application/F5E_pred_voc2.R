rm(list = ls())
library(rstan)
library(ggplot2)
library(scales)
library(ggnewscale)
library(tidyverse)
library(dplyr)
library(RColorBrewer)

{
  df = read.csv('F5C_VOC_gisaid.csv')
  
  VOC = c("Omicron", "Non Omicron","Non Omicron","Non Omicron")
  names(VOC) = c("GRA","GK","GRY","G")
  df = df[df$GISAID_clade %in% names(VOC),]
  df$V = NA
  for (i in 1:length(VOC)) {
    df$V[df$GISAID_clade == names(VOC)[i]] = VOC[i]
  }
  df$V = factor(df$V, levels = VOC[1:2])
  df = df %>% group_by(date, V) %>%
    summarise(y = sum(count))
  
  df <- na.omit(df)
  max(df$date)
  as.Date('2020-12-31') + 700
  full_dates <- as.Date('2020-12-31') + 1:700
  # Create a dataframe with all dates
  full_df <- expand.grid(date = full_dates, 
                         V = unique(df$V))
  df$date = as.Date(df$date)
  merged_df <- full_df %>%
    left_join(df, by = c("date", "V")) %>%
    mutate(count = ifelse(is.na(y), 0, y))
  voc = unique(merged_df$V)
  
  observed_matrix = data.frame(v1 = merged_df[merged_df$V == 'Non Omicron', 'count'],
                               v2 = merged_df[merged_df$V == 'Omicron', 'count'])
  rownames(observed_matrix) = full_dates
  fexpect0 = data.frame(y = c(observed_matrix$v1,observed_matrix$v2),
                        x = rep(1:nrow(observed_matrix),2),
                        group = rep(c('Non Omicron', 'Omicron'), 
                                    each = nrow(observed_matrix)))
  fexpect0$date = as.Date('2020-12-31') + fexpect0$x
  fexpect0$group = factor(fexpect0$group, levels = c('Non Omicron','Omicron'))
  
}
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

determinant_fun = function(cond = F, ifsimu  = F, n_simu = 1, predj = 1){
  # cond = F; ifsimu  = F; n_simu = 1
  n = 2
  poolday = 30
  nday = 100
  
  # Create a dataframe with all dates
  d0 <- expand.grid(date = as.Date('2020-12-20') + 1:10, 
                    V = unique(df$V))
  d1 <- d0 %>%
    left_join(df, by = c("date", "V")) %>%
    mutate(y = ifelse(is.na(y), 0, y))
  d1 = d1 %>% group_by(V) %>% summarise(m = mean(y)) %>% as.data.frame()
  
  
  seed_matrix <- matrix(0, nrow=n, ncol=poolday)
  seed_matrix[,1] = c(d1[d1$V == 'Non Omicron','m'],
                      d1[d1$V == 'Omicron','m']) + 1e-1
  seed_vec = colSums(seed_matrix)
  seed_mats <- list()
  
  for (i in 1:n) {
    seed_mats[[i]] <- diag(seed_matrix[i,])
  }
  
  expected_matrix = observed_matrix[1:nday,]
  expected_matrix[expected_matrix < 0] = 0
  expected_matrix = round(expected_matrix)
  expected_matrix[(poolday*2+1):nday,] = 0
  
  Onsets_mat_list = list()
  fit= fitlist[[1]]
  
  posterior = rstan::extract(fit)
  
  if(!ifsimu){
    fexpect = data.frame(x = rep(1:nday,2), 
                         y = c(expected_matrix[1:nday,1],
                               expected_matrix[1:nday,2]),
                         group = factor(rep(c('Non Omicron','Omicron'), 
                                            each = nday),
                                        levels = c('Non Omicron','Omicron')))
    pars_last = c(mean(posterior$contact), 
                  mean(posterior$beta1),
                  mean(posterior$beta2))
    Onsets_mat = simu(seed_mats, 
                      N = seed_vec * pars_last[1] + 1, 
                      poolday, pars = pars_last[-1], n)
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
  
  if(ifsimu){
    pars_last = c(posterior$contact[n_simu], 
                  posterior$beta1[n_simu],
                  posterior$beta2[n_simu])
    Onsets_mat = simu(seed_mats, 
                      N = seed_vec * pars_last[1] + 1, 
                      poolday, pars = pars_last[-1], n)
  }
  
  Onsets_mat_list[[1]] = Onsets_mat
  # predj = 1

  for (j in 1:predj) {
    {
      
      Onsets_mat = Onsets_mat_list[[j]]
      
      # Generalized extraction of Onset columns for all variants
      Onsets <- list()
      for (i in 1:n) {
        Onsets[[i]] <- Onsets_mat[poolday + 1:poolday, i] + 1e-1
      }
      # Mobility control force
      mobility <- rep(1/30, 30) 
      Mobility_matrix <- diag(mobility)
      
      # Create seed_vec considering all variants
      seed_vec <- rowSums(sapply(Onsets, function(Onset) Onset)) %*% Mobility_matrix %>% as.numeric()
      
      # Calculate probabilities for each variant
      probs <- lapply(Onsets, function(Onset) Onset / rowSums(do.call(cbind, Onsets)))
      seed_matrix <- matrix(0, nrow=n, ncol=length(seed_vec))
      
      if(ifsimu & min(max(Onsets[[1]]), max(Onsets[[2]])) > 0.5){
        
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
      
      fit = fitlist[[j+1]]
      if(j == predj) fit = fitlist[[j]]
      posterior = rstan::extract(fit)
      
      if(!ifsimu){
        pars_last = c(mean(posterior$contact), 
                      mean(posterior$beta1),
                      mean(posterior$beta2))
        Onsets_mat = simu(seed_mats, 
                          N = seed_vec * pars_last[1] + 1, 
                          poolday, pars = pars_last[-1], n)
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
      
      if(ifsimu){
        pars_last = c(posterior$contact[n_simu], 
                      posterior$beta1[n_simu],
                      posterior$beta2[n_simu])
        Onsets_mat = simu(seed_mats, 
                          N = seed_vec * pars_last[1] + 1, 
                          poolday, pars = pars_last[-1], n)
      }
      
    }
    
    Onsets_mat_list[[j+1]] = Onsets_mat
  }
  
  # save(fitlist, file = 'voc2_chain4.rdata')
  dfplot_simu = data.frame()
  for (i in 1:(predj+1)) {
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

if(F){
  # fitlist = list()
  load('F5C_voc2.rdata')
}


dfall = data.frame()

for (j in 1:20) {
  print(j)
  data = determinant_fun(predj = j)
  data = data[data$x > (30*j) & data$x <= (30*(j+1)),]
  ob1 = fexpect0[fexpect0$x > (30*j) & fexpect0$x <= (30*(j+1)),]
  data$group = 'Non Omicron'
  data$group[data$color == '2'] = 'Omicron'
  merged_df <- ob1 %>%
    left_join(data, by = c("x", "group")) 
  merged_df$pred = j
  merged_df$x = merged_df$x - 30*j
  # dfall  = merged_df
  dfall = rbind(dfall, merged_df)
}

dfall = dfall[dfall$y.x>1,]
df1 = dfall[dfall$color == 1,]
fit1 <- lm(log10(y.y) ~ log10(y.x), data = df1)
slope1 <- coef(fit1)[2]
intercept1 <- coef(fit1)[1]
adj_r2_1 <- summary(fit1)$adj.r.squared
coord1 = log10(c(min(df1[, c('y.x','y.y')]), max(df1[, c('y.x','y.y')])))

df2 = dfall[dfall$color == 2,]
fit2 <- lm(log10(y.y) ~ log10(y.x), data = df2)
slope2 <- coef(fit2)[2]
intercept2 <- coef(fit2)[1]
adj_r2_2 <- summary(fit2)$adj.r.squared
coord2 = log10(c(min(df2[, c('y.x','y.y')]), max(df2[, c('y.x','y.y')])))

p1 = ggplot(data = df1, aes(x = log10(y.x), y =log10(y.y),
                       group = x, color =x)) +
  geom_point(alpha = 0.6) +
  geom_abline(slope = slope1, intercept = intercept1) +
  annotate("text", x = 3.4, y = 0.4, size = 3,
           label = paste0("y = ", round(slope1, 2), "x + ", round(intercept1, 2), #1.01x
                          "\nR² = ", round(adj_r2_1, 3))) +
  xlab(expression(log[10](Expected))) + ylab(expression(log[10](Predicted))) +
  theme_bw() +
  coord_cartesian(xlim = coord1, ylim = coord1) +
  scale_x_continuous(breaks = c(0,2,4)) + scale_y_continuous(breaks = c(0,2,4)) +
  theme(panel.grid.minor = element_blank(),
        legend.position = 'right') +
  scale_color_continuous(name = 'Day')

p2 = ggplot(data = df2, aes(x = log10(y.x), y =log10(y.y),
                       group = x, color =x)) +
  geom_point(alpha = 0.6) +
  geom_abline(slope = slope2, intercept = intercept2) +
  annotate("text", x = 2.5, y = -0.5, size = 3,
           label = paste0("y = x + ", round(intercept2, 2), # 1x
                          "\nR² = ", round(adj_r2_2, 3))) +
  xlab(expression(log[10](Expected))) + ylab(expression(log[10](Predicted))) +
  theme_bw() +
  coord_cartesian(xlim = coord2, ylim = coord2) +
  scale_x_continuous(breaks = c(0,2,4)) + scale_y_continuous(breaks = c(0,2,4)) +
  theme(panel.grid.minor = element_blank(),
        legend.position = 'right') +
  scale_color_continuous(name = 'Day')

pdf(paste0("Output/F5E_pred_voc2.pdf"), width = 3, height = 2)
print(p1)
print(p2)
dev.off()

