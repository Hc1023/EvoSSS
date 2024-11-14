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

fitlist = list()
load('flu.rdata')
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
    fit <- stan(file = 'FLU.stan', data = stan_data, 
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
        fit <- stan(file = 'FLU.stan', data = stan_data, 
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
  
  # j = 35 last simulation -- 24-July
  
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
  
  # save(fitlist, Onsets_mat_list, file = 'flu.rdata')
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
  if(ifsimu){
    return(list(data, Onsets_mat_list))
  }
  return(data)
}

dfob = data.frame(y = c(as.matrix(observed_matrix)),
                  x = rep(as.Date(rownames(observed_matrix)),
                          4),
                  group = rep(voc, each = nrow(observed_matrix)))

df2_list = list()
dfsimu2 = data.frame()
for (n_simu in 1:1000) {
  print(n_simu)
  datalist = determinant_fun(cond = F, ifsimu  = T, n_simu = n_simu)
  data = datalist[[1]]
  df2_list[[n_simu]] = data$y
  Onsets_mat_list = datalist[[2]]
  tmp = sapply(Onsets_mat_list, function(x){
    colSums(x)
  })
  dfsimu = data.frame(y = c(t(tmp)),
                      group = rep(voc, each = length(Onsets_mat_list)),
                      x = rep(dfob$x[1:length(Onsets_mat_list)], 4))
  dfsimu2 = rbind(dfsimu2, dfsimu)
}

simu_Onset = data.frame(bind_cols(df2_list))
if(F){
  save(simu_Onset, dfsimu2, file = 'flu_plot.rdata')
  load('flu_plot.rdata')
  datalist = determinant_fun(cond = F, ifsimu  = T, n_simu = 1)
  data  = datalist[[1]]
}


ci_lower <- apply(simu_Onset, 1, quantile, probs = 0.025, na.rm = T)
ci_upper <- apply(simu_Onset, 1, quantile, probs = 0.975, na.rm = T)
plot_data <- data.frame(
  x = data$x,
  group = data$color,
  Fitted = rowMeans(simu_Onset),
  LowerCI = ci_lower,
  UpperCI = ci_upper
)

values = c('#98afc7','#0041c2', "#ffa500",'#a37ca1')

plot_data$date = as.Date('2021-07-01') + plot_data$x
plot_data$group = factor(plot_data$group, levels = voc)


dfsimu3 = dfsimu2 %>% group_by(x, group) %>% 
  summarise(ym = mean(y),
            y1 = quantile(y, 0.025),
            y2 = quantile(y, 0.975)) %>%
  as.data.frame()

dfsimu3$group = factor(dfsimu3$group, levels = voc)

k = 30
pd = 25
preddate = as.Date('2024-08-01')
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
        plot.title = element_text(hjust = 0.5),
        panel.grid.minor = element_blank()) +
  scale_x_date(breaks = seq(as.Date('2021-01-01'), as.Date('2025-05-01'), by="6 months"),
               minor_breaks = seq(as.Date('2021-01-01'), as.Date('2025-05-01'), by ='1 month'),
               date_labels = "%y-%b",
               expand = c(0,0)) +
  xlab('') + ylab('Daily cases') + 
  coord_cartesian(xlim = c(as.Date('2021-09-01'), as.Date('2025-05-31')),
                  ylim = c(1,max(plot_data$UpperCI)))

pdf(paste0("Output/flu_plot.pdf"), width = 4, height = 1.6)
print(p)
dev.off()

plot_data2 <- plot_data %>%
  group_by(x) %>%
  mutate(p = Fitted/sum(Fitted)) %>%
  as.data.frame()
plot_data2$group = factor(plot_data2$group, 
                          levels = rev(levels(plot_data$group)))


df = read.csv('influenza.csv')
df2 = data.frame(x = as.Date(rep(df$X, 3)),
  y = c(observed_matrix$v1/rowSums(observed_matrix), 
        (observed_matrix$v1 + observed_matrix$v2)/rowSums(observed_matrix),
        (observed_matrix$v1 + observed_matrix$v2 + observed_matrix$v3)/rowSums(observed_matrix)),
  group = rep(voc[1:3], each = nrow(df)))
p2 =  ggplot() +
  geom_area(data = plot_data2, 
            aes(x = date, y = p, fill = group),
            position = 'fill') +
  scale_fill_manual(name="", breaks = rev(levels(plot_data2$group)),
                    values = alpha(values, 1)) +
  scale_x_date(breaks = seq(as.Date('2021-01-01'), as.Date('2025-05-01'), by="6 months"),
               minor_breaks = seq(as.Date('2021-01-01'), as.Date('2025-05-01'), by ='1 month'),
               date_labels = "%y-%b",
               expand = c(0,0)) +
  scale_y_continuous(limits = c(0,1), expand = c(0,0),
                     breaks = c(0,0.5,1)) +
  ylab('Prevalence') + xlab('') +
  theme(legend.position = 'none') +
  coord_cartesian(xlim = c(as.Date('2021-09-01'), as.Date('2025-05-31')))
p2
pdf(paste0("Output/flu_plot_prevalence.pdf"), width = 4, height = 1.2)
print(p2)
dev.off()

parsmat = data.frame()
for (i in 1:length(fitlist)) {
  posterior = rstan::extract(fitlist[[i]])
  parsmat = rbind(parsmat, 
                  data.frame(Tcycle = i,
                             y = c(posterior$contact, 
                                   posterior$beta1,
                                   posterior$beta2,
                                   posterior$beta3,
                                   posterior$beta4),
                             group = rep(c('contact','beta1','beta2','beta3','beta4'), 
                                         each = length(posterior$contact))
                             
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
  geom_boxplot(data = df[df$group == 'beta4',], 
               aes(x = date+2*pd, y = y, group = date),
               color = alpha(values[4],0.7), 
               fill = alpha(values[4],0.3),
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
  geom_line(data = medians[medians$group == 'beta4',], 
            aes(x = date + 2*pd, y = median_y), 
            color = values[4], 
            linewidth = 0.4) +
  coord_cartesian(xlim = c(as.Date('2021-09-01'), 
                           as.Date('2024-07-06')),
                  ylim = c(0.24,0.37)) +
  scale_y_continuous(breaks = c(0.25,0.3,0.35)) +
  theme_bw() + xlab('') + ylab(expression(beta)) +
  theme(panel.grid.minor = element_blank())
pdf(paste0("Output/flu_plot_beta.pdf"), width = 3.3, height = 1.2)
print(p3)
dev.off()


df = parsmat[parsmat$group == 'contact', ]
df$date = as.Date(df$date)
df2 = df
df2$group = 'contact2' 
df2$date = df2$date + years(1)
df3 = df2
df3$group = 'contact3' 
df3$date = df2$date + years(1)
pd = 8

medians <- df %>%
  group_by(date, group) %>%
  summarize(median_y = median(y))

medians2 <- df2 %>%
  group_by(date, group) %>%
  summarize(median_y = median(y))

medians3 <- df3 %>%
  group_by(date, group) %>%
  summarize(median_y = median(y))

values = c('#98afc7','#0041c2', "#ffa500",'#a37ca1')
p4 = ggplot() +
  geom_boxplot(data = df, 
               aes(x = date, y = y, group = date),
               color = alpha('#3d3da1',0.9), 
               fill = alpha('#3d3da1',0.3),
               width = 6,
               outlier.shape = NA) +
  geom_boxplot(data = df2, 
               aes(x = date, y = y, group = date),
               color = alpha('#adadc1',0.9), 
               fill = alpha('#adadc1',0.3),
               width = 6,
               outlier.shape = NA) +
  geom_boxplot(data = df3, 
               aes(x = date, y = y, group = date),
               color = alpha('#ad8d81',0.9), 
               fill = alpha('#ad8d81',0.3),
               width = 6,
               outlier.shape = NA) +
  scale_x_date(breaks = seq(as.Date('2023-01-01'), as.Date('2025-12-01'), by="3 months"),
               minor_breaks = seq(as.Date('2023-01-01'), as.Date('2025-11-01'), by ='1 month'),
               date_labels = "%b",
               expand = c(0, 0)) +
  geom_line(data = medians, 
            aes(x = date, y = median_y), 
            color = '#3d3da1', 
            linewidth = 0.4) +
  geom_line(data = medians2, 
            aes(x = date, y = median_y), 
            color = '#adadc1', 
            linewidth = 0.4) +
  geom_line(data = medians3, 
            aes(x = date, y = median_y), 
            color = '#ad8d81', 
            linewidth = 0.4) +
  coord_cartesian(xlim = c(as.Date('2023-07-01'), 
                           as.Date('2024-06-30'))) +
  scale_y_continuous() +
  theme_bw() + xlab('') + ylab('')

pdf(paste0("Output/flu_plot_contact.pdf"), width = 1.65, height = 1.2)
print(p4)
dev.off()



if(F){
  data = determinant_fun(cond = F, ifsimu  = T, n_simu = n_simu)
  data$date = as.Date('2021-07-01') + data$x
  dfob = data.frame(y = c(observed_matrix[,1],
                          observed_matrix[,2],
                          observed_matrix[,3],
                          observed_matrix[,4]),
                    x = rep(as.Date(rownames(observed_matrix)),
                            4),
                    group = rep(voc, each = nrow(observed_matrix)))
  
  tmp = sapply(Onsets_mat_list, function(x){
    colSums(x)
  })
  dfsimu = data.frame(y = c(tmp[1,],tmp[2,],tmp[3,],tmp[4,]),
                      group = rep(voc, each = length(Onsets_mat_list)),
                      x = rep(dfob$x[1:length(Onsets_mat_list)], 4))
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


