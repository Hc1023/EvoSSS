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


determinant_fun = function(cond = T, ifsimu  = T, n_simu = 1){
  # cond = T; ifsimu  = F; n_simu = 1
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
  
  
  if(cond){
    stan_data <- list(
      poolday = poolday,
      nday = nday,
      expected_matrix = expected_matrix,
      seed_mat_I1 = seed_mats[[1]],
      seed_mat_I2 = seed_mats[[2]],
      seed_vec = seed_vec,
      gamma = 0.157,
      pars_last = c(200,0.3,0.3)
    )
    # Fit the model
    fit <- stan(file = 'VOC2.stan', data = stan_data, 
                iter = 3000, chains = 4, warmup = 2000,
                verbose = TRUE)
    fitlist[[1]] = fit
  }
  
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
  # cond = T
  # cond = F
  # 1:20
  for (j in 1:20) {
    if(cond) print(j)
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
      
      if(cond){
        stan_data <- list(
          poolday = poolday,
          nday = nday,
          expected_matrix = expected_matrix,
          seed_mat_I1 = seed_mats[[1]],
          seed_mat_I2 = seed_mats[[2]],
          seed_vec = seed_vec,
          gamma = 0.157,
          pars_last = pars_last
        )
        # Fit the model
        fit <- stan(file = 'F5C_VOC2.stan', data = stan_data, 
                    iter = 3000, chains = 4, warmup = 2000,
                    verbose = TRUE)
        fitlist[[j+1]] = fit
      }
      
      
      fit = fitlist[[j+1]]
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
  for (i in 1:21) {
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
  fitlist = list()
  save(fitlist, file = 'voc2.rdata')
  load('voc2.rdata')
}
if(F){
  
  df2_list = list()
  for (n_simu in 1:100) {
    print(n_simu)
    data = determinant_fun(cond = F, ifsimu  = T, n_simu = n_simu)
    df2_list[[n_simu]] = data$y
  }
  {
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
    plot_data$date = plot_data$x + as.Date('2020-12-31')
    plot_data$group = 'Non Omicron'
    plot_data$group[data$color == '2'] = 'Omicron'
    plot_data$group = factor(plot_data$group, levels = c('Non Omicron','Omicron'))
    save(fitlist, plot_data, file = 'voc2.rdata')
    
  }
  
}
load('F5C_voc2.rdata')
values = brewer.pal(n=4, name = "Spectral")[1:2]
values = rev(c(values[1], '#41424c'))
plot_data = plot_data[plot_data$x>1,]

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
  scale_color_manual(name="",
                     values = alpha(values, 0.7)) +
  scale_fill_manual(name="",
                    values = alpha(values, 0.3)) +
  coord_cartesian(ylim = c(2,max(plot_data$Fitted))) +
  theme_bw() +
  scale_x_date(breaks = seq(as.Date('2020-01-01'), as.Date('2022-11-01'), by="6 months"),
               minor_breaks = seq(as.Date('2019-12-01'), as.Date('2022-05-01'), by ='1 month'),
               date_labels = "%y-%b", expand = c(0, 0)) +
  xlab('') + ylab(expression(Cases~(10^4))) + 
  scale_y_continuous(breaks = 0:3*10^4, labels = 0:3) +
  coord_cartesian(xlim = c(as.Date('2021-01-01'), as.Date('2022-08-01')),
                  ylim = c(1,3.3*10^4)) +
  theme(legend.position = 'none',
        legend.background = element_rect(color = NA, fill = NA),
        legend.key = element_blank(),
        legend.key.size = unit(0.2, units = 'cm'),
        panel.grid.minor = element_blank())

pdf(paste0("Output/F5C_voc2.pdf"), width = 2.5, height = 1.8)
print(p)
dev.off()

iF(F){
  plot_data2 <- plot_data %>%
    group_by(x) %>%
    mutate(p = Fitted/sum(Fitted)) %>%
    as.data.frame()
  plot_data2$group = factor(plot_data2$group, 
                            levels = rev(levels(plot_data$group)))
  p2 =  ggplot() +
    geom_area(data = plot_data2, 
              aes(x = date, y = p, fill = group),
              position = 'fill') +
    scale_fill_manual(name="", breaks = rev(levels(plot_data2$group)),
                      values = alpha(values, 1)) +
    coord_cartesian(xlim = c(as.Date('2021-01-01'), as.Date('2022-08-01')),
                    ylim = c(0,1)) +
    scale_x_date(breaks = seq(as.Date('2020-01-01'), as.Date('2024-11-01'), by="6 months"),
                 minor_breaks = seq(as.Date('2019-12-01'), as.Date('2024-11-01'), by ='1 month'),
                 date_labels = "%y-%b",
                 expand = c(0, 0)) +
    scale_y_continuous(breaks = seq(0,1,0.5), expand = c(0, 0)) +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_rect(color = NA, fill = NA),
          legend.position = 'none',
          legend.background = element_rect(color = NA, fill = NA),
          legend.key = element_blank(),
          legend.key.size = unit(0.4, units = 'cm')) +
    xlab('') + ylab('Prevalence') 
  
  pdf(paste0("Output/voc2_plot_prevalence.pdf"), width = 2.5, height = 1.2)
  print(p2)
  dev.off()
}


parsmat = data.frame()
for (i in 1:length(fitlist)) {
  posterior = rstan::extract(fitlist[[i]])
  parsmat = rbind(parsmat, 
                  data.frame(Tcycle = i,
                             y = c(posterior$contact, 
                                   posterior$beta1,
                                   posterior$beta2),
                             group = c(rep('contact', length(posterior$contact)),
                                       rep('beta1', length(posterior$beta1)),
                                       rep('beta2', length(posterior$beta2)))
                             
                  ))
  
}

poolday = 30
parsmat$date = poolday*(parsmat$Tcycle-1) + as.Date('2020-12-31') + 1
df = parsmat
df$date = as.Date(df$date)
pd = 0

medians <- df %>%
  group_by(date, group) %>%
  summarize(median_y = median(y))

values = brewer.pal(n=4, name = "Spectral")[1:2]
values = rev(c(values[1], '#41424c'))
p3 = ggplot() +
  geom_boxplot(data = df[df$group == 'contact',], 
               aes(x = date-pd, y = y, group = date),
               color = alpha('red',0.7), 
               fill = alpha('red',0.3),
               width = 6,
               outlier.shape = NA) +
  scale_x_date(breaks = seq(as.Date('2020-01-01'), as.Date('2024-11-01'), by="6 months"),
               minor_breaks = seq(as.Date('2019-12-01'), as.Date('2024-11-01'), by ='1 month'),
               date_labels = "%y-%b",
               expand = c(0, 0)) +
  geom_line(data = medians[medians$group == 'contact',], 
            aes(x = date - pd, y = median_y), 
            color = 'red', 
            linewidth = 0.4) +
  coord_cartesian(xlim = c(as.Date('2021-01-01'), as.Date('2022-08-01'))) +
  # scale_y_continuous(breaks = c(0.1,0.3,0.5,0.7)) +
  theme_bw() + xlab('') + ylab('Q') + 
  theme(panel.grid.minor = element_blank())

p4 = ggplot() +
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
  coord_cartesian(xlim = c(as.Date('2021-01-01'), as.Date('2022-08-01'))) +
  scale_y_continuous(breaks = c(0.1,0.3,0.5,0.7)) +
  theme_bw() + xlab('') + ylab(expression(beta)) + 
  theme(panel.grid.minor = element_blank())
pdf(paste0("Output/F5C_voc2_beta.pdf"), width = 2.5, height = 1.2)
print(p3)
print(p4)
dev.off()



if(F){
  data = determinant_fun(cond = F, ifsimu  = F, n_simu = n_simu)
  data$date = as.Date('2020-12-31') + data$x
  data$group = 'Non Omicron'
  data$group[data$color == '2'] = 'Omicron'
  data$group = factor(data$group, levels = c('Non Omicron','Omicron'))
  fexpect0 = data.frame(y = c(observed_matrix$v1,observed_matrix$v2,
                              observed_matrix$v3),
                        x = rep(1:nrow(observed_matrix),2),
                        group = rep(c('Non Omicron','Omicron'), 
                                    each = nrow(observed_matrix)))
  fexpect0$date = as.Date('2020-12-31') + fexpect0$x
  fexpect0$group = factor(fexpect0$group, levels = c('Non Omicron','Omicron'))
  
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
    coord_cartesian(xlim = c(as.Date('2021-2-05'), as.Date('2022-10-01')),
                    ylim = c(1,2*10^4))

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

iF(F){
  pdf(file = 'Output/trace_voc2.pdf', 
      width = 9, height = 2)
  for (i in 1:21) {
    print(traceplot(fitlist[[i]]))
  }
  dev.off()
}


