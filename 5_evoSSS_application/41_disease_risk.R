rm(list = ls())
library(rstan)
library(ggplot2)
library(scales)
library(ggnewscale)
library(tidyverse)
library(dplyr)
library(data.table)
library(lubridate)
library(ISOweek)
df1 = fread('VIW_FNT.csv')
df1[is.na(df1)] = 0

ob1 = df1[,c('ISO_WEEKSTARTDATE','INF_A','INF_B','INF_ALL')]
colnames(ob1)[1] = 'x'
ob1 = ob1 %>% 
  group_by(x) %>%
  summarise(yA = sum(INF_A),
            yB = sum(INF_B),
            yALL = sum(INF_ALL)) %>%
  as.data.frame()


dates = c(as.Date('2017-08-01'),as.Date('2021-08-01'))
{
  observed_matrix = ob1[ob1$x >= dates[1] & ob1$x < dates[2],-4] 
  observed_matrix$x = format(observed_matrix$x, "%Y-%m")
  observed_matrix = observed_matrix %>% group_by(x) %>%
    summarise(yA=sum(yA), yB=sum(yB)) %>% 
    as.data.frame()
  rownames(observed_matrix) = as.Date(paste0(unlist(observed_matrix$x),'-01')) + months(1)-1
  observed_matrix = observed_matrix[,-1]
  voc = c('A','B')
  colnames(observed_matrix) = voc
  
  dfob = data.frame(y = c(as.matrix(observed_matrix)),
                    x = rep(as.Date(rownames(observed_matrix)),2),
                    group = rep(voc, each = nrow(observed_matrix)))
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



fitlist = list()

determinant_fun = function(cond = T, ifsimu  = T, n_simu = 1){
  # cond = T; ifsimu  = F; n_simu = 1
  n = 2
  poolday = 30
  nday = 100
  seed_matrix <- matrix(0, nrow=n, ncol=poolday)
  seed_matrix[,1] = unlist(observed_matrix[1,]/30+1)
  
  if(ifsimu){
    seed_matrix[,1] = rmultinom(1, sum(seed_matrix[,1]), seed_matrix[,1]/sum(seed_matrix[,1]))
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
      expected_total = unlist(round(sqrt(observed_matrix[1,1:2]))),
      seed_mat_I1 = seed_mats[[1]],
      seed_mat_I2 = seed_mats[[2]],
      seed_vec = seed_vec,
      gamma = 0.157,
      pars_last = c(200, rep(0.3,2))
    )
    # Fit the model
    fit <- stan(file = 'flu_covid.stan', data = stan_data, 
                iter = 3000, chains = 1, warmup = 2000,
                verbose = TRUE)
    fitlist[[1]] = fit
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
  if(!ifsimu){
    fonset = data.frame(x = rep(1:nday,2), 
                        y = c(Onsets_mat[1:nday,1],
                              Onsets_mat[1:nday,2]),
                        group = factor(rep(voc, each = nday), 
                                       levels = voc))
    
    tmp = fonset %>% group_by(group) %>% 
      summarise(y = sum(y)) %>% 
      data.frame()
    tmp$x = c(55,60)
    
    observed_df = data.frame(y = unlist(observed_matrix[1,]))
    
    observed_df$x = c(25,30)
    observed_df$group = voc
    observed_df = rbind(observed_df, tmp)
    
    ggplot() +
      geom_line(data = fonset,
                aes(x = x, y = y+1, 
                    group = group, color = group))+
      geom_col(data = observed_df, 
               aes(x = x, y = y/30+1, 
                   group = group, 
                   color = group, fill = group), 
               width = 5) + scale_y_continuous(transform = 'log')
    
  }
  
  Onsets_mat_list[[1]] = Onsets_mat
  
  for (j in 1:47) {
    if(!ifsimu)print(j)
    {
      Onsets_mat = Onsets_mat_list[[j]]
      
      # Generalized extraction of Onset columns for all variants
      Onsets <- list()
      for (i in 1:n) {
        Onsets[[i]] <- Onsets_mat[poolday + 1:poolday, i] + 1
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
          expected_total = unlist(round(sqrt(observed_matrix[j+1,1:n]))),
          seed_mat_I1 = seed_mats[[1]],
          seed_mat_I2 = seed_mats[[2]],
          seed_vec = seed_vec,
          gamma = 0.157,
          pars_last = pars_last
        )
        # Fit the model
        fit <- stan(file = 'flu_covid.stan', data = stan_data, 
                    iter = 3000, chains = 1, warmup = 2000,
                    verbose = TRUE)
        fitlist[[j+1]] = fit
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
      
      if(!ifsimu){
        fonset = data.frame(x = rep(1:nday,2), 
                            y = c(Onsets_mat[1:nday,1],
                                  Onsets_mat[1:nday,2]),
                            group = factor(rep(voc, each = nday), 
                                           levels = voc))
        
        
        tmp = fonset %>% group_by(group) %>% 
          summarise(y = sum(y)) %>% 
          data.frame()
        tmp$x = c(55,60)
        
        
        observed_df = data.frame(y = unlist(observed_matrix[j+1,]))
        
        observed_df$x = c(25,30)
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
                   width = 5) +
          scale_y_continuous(transform = 'log')
      }
    }
    Onsets_mat_list[[j+1]] = Onsets_mat
  }
  
  dfplot_simu = data.frame()
  for (i in 1:48) {
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
  
  if(F){
    data$date = dates[1]-1 + data$x
    dfob = data.frame(y = c(observed_matrix[,1],
                            observed_matrix[,2]),
                      x = rep(as.Date(rownames(observed_matrix)),2),
                      group = rep(voc, each = nrow(observed_matrix)))
    
    
    ggplot() +
      geom_point(data = dfob,
                 aes(x = x, y = y/30+1, 
                     group = group, color = group, fill = group),
                 alpha = 0.8,
                 position = position_dodge(width = 20)) +
      geom_line(data = data, 
                aes(x = date, y = y+1, group = group, color = group)) +
      scale_y_continuous(transform = 'log') +
      theme(legend.position = "right",
            plot.title = element_text(hjust = 0.5)) +
      scale_x_date(breaks = seq(as.Date('2019-01-01'), as.Date('2024-11-01'), by="1 year"),
                   minor_breaks = seq(as.Date('2019-01-01'), as.Date('2024-11-01'), by ='2 months'),
                   date_labels = "%y-%b") +
      xlab('') + ylab('Cases')
    
  }
  return(data)
}

if(F){
  save(fitlist, file = 'disease_risk.rdata')
  save(fitlist, plot_data, file = 'disease_risk.rdata')
  load(file = 'disease_risk.rdata')
}


if(F){
  df2_list = list()
  
  for (n_simu in 1:1000) {
    print(n_simu)
    data = determinant_fun(cond = F, ifsimu  = T, n_simu = n_simu)
    df2_list[[n_simu]] = data$y
  }
  
  simu_Onset = data.frame(bind_cols(df2_list))
  ci_lower <- apply(simu_Onset, 1, quantile, probs = 0.025, na.rm = T)
  ci_upper <- apply(simu_Onset, 1, quantile, probs = 0.975, na.rm = T)
  plot_data <- data.frame(
    x = data$x,
    group = data$color,
    Fitted = rowMeans(simu_Onset),
    LowerCI = ci_lower,
    UpperCI = ci_upper
  )
  
  plot_data$date = dates[1]-1 + plot_data$x
  plot_data$group = factor(plot_data$group, levels = voc)
  
}

k = 30
pd = 25

values = c('#98aff7','#a37c91')
p = ggplot() +
  geom_point(data = dfob,
             aes(x = x, y = y/k,
                 group = group, color = group),
             alpha = 0.8) +
  geom_ribbon(data = plot_data,
              aes(x = date, group = group,
                  ymin = LowerCI, ymax = UpperCI, fill = group)) +
  geom_line(data = plot_data, 
            aes(x = date, y = Fitted, group = group, color= group)) +
  scale_color_manual(name="",
                     values = alpha(values, 0.7)) +
  scale_fill_manual(name="",
                    values = alpha(values, 0.3)) +
  theme_bw() +
  scale_x_date(breaks = seq(as.Date('2017-01-01'), as.Date('2024-11-01'), by="6 months"),
               minor_breaks = seq(as.Date('2017-01-01'), as.Date('2024-11-01'), by ='1 month'),
               date_labels = "%y-%b",
               expand = c(0, 0)) +
  xlab('') + ylab('Daily cases') +
  coord_cartesian(xlim = c(dates[1], dates[2]), ylim = c(0,20000)) +
  theme(legend.position = "none",
        panel.grid.minor = element_blank(),
        plot.margin = margin(0.2,0.6,0,0, "cm"))

pdf(paste0("Output/disease_risk.pdf"), width = 4.4, height = 1.6)
print(p)
dev.off()

parsmat = data.frame()
for (i in 2:length(fitlist)) {
  posterior = rstan::extract(fitlist[[i]])
  parsmat = rbind(parsmat, 
                  data.frame(Tcycle = i,
                             y = c(posterior$contact, 
                                   posterior$beta1,
                                   posterior$beta2),
                             group = rep(c('contact','beta1','beta2'), 
                                         each = length(posterior$contact))
                             
                  ))
  
}


poolday = 30
parsmat$date = poolday*(parsmat$Tcycle-1) + dates[1]

df = parsmat
df$date = as.Date(df$date)
pd = 4

medians <- df %>%
  group_by(date, group) %>%
  summarize(median_y = median(y))

p3 = ggplot() +
  geom_boxplot(data = df[df$group == 'contact',], 
               aes(x = date, y = y, group = date),
               color = alpha('red',0.7), 
               fill = alpha('red',0.3),
               width = 6, outlier.shape = NA) +
  geom_line(data = medians[medians$group == 'contact',], 
            aes(x = date, y = median_y), 
            color = 'red', 
            linewidth = 0.4) +
  scale_x_date(breaks = seq(as.Date('2017-01-01'), as.Date('2024-11-01'), by="6 months"),
               minor_breaks = seq(as.Date('2017-01-01'), as.Date('2024-11-01'), by ='1 month'),
               date_labels = "%y-%b",
               expand = c(0, 0)) +
  coord_cartesian(xlim = c(dates[1],dates[2])) +
  # scale_y_continuous(breaks = c(0.3,0.4)) +
  theme_bw() + xlab('') + ylab('Q') +
  theme(plot.margin = unit(c(0.3,0.6,0,0.3), "cm"),
        panel.grid.minor = element_blank())

p4 = ggplot() +
  geom_boxplot(data = df[df$group == 'beta1',], 
               aes(x = date-pd, y = y, group = date),
               color = alpha(values[1],0.7), 
               fill = alpha(values[1],0.3),
               width = 6,
               outlier.shape = NA) +
  geom_boxplot(data = df[df$group == 'beta2',], 
               aes(x = date+pd, y = y, group = date),
               color = alpha(values[2],0.7), 
               fill = alpha(values[2],0.3),
               width = 6,
               outlier.shape = NA) +
  scale_x_date(breaks = seq(as.Date('2018-01-01'), as.Date('2024-11-01'), by="6 months"),
               minor_breaks = seq(as.Date('2018-01-01'), as.Date('2024-11-01'), by ='1 month'),
               date_labels = "%y-%b",
               expand = c(0, 0)) +
  geom_line(data = medians[medians$group == 'beta1',], 
            aes(x = date - pd, y = median_y), 
            color = values[1], 
            linewidth = 0.4) +
  geom_line(data = medians[medians$group == 'beta2',], 
            aes(x = date + pd, y = median_y), 
            color = values[2], 
            linewidth = 0.4) +
  coord_cartesian(xlim = c(dates[1],dates[2])) +
  # scale_y_continuous(breaks = c(0.3,0.4)) +
  theme_bw() + xlab('') + ylab(expression(beta)) +
  theme(plot.margin = unit(c(0.3,0.6,0,0.3), "cm"),
        panel.grid.minor = element_blank())
pdf(paste0("Output/disease_risk_par.pdf"), width = 4.3, height = 1)
print(p3)
print(p4)
dev.off()


