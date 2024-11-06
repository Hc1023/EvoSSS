rm(list = ls())
library(rstan)
library(ggplot2)
library(scales)
library(ggnewscale)
library(tidyverse)
library(dplyr)
library(data.table)
library(lubridate)
df1 = fread('VIW_FNT.csv')
df1[is.na(df1)] = 0
ob1 = df1[,c('ISO_WEEKSTARTDATE','INF_ALL')]
ob1$x = as.character(format(ob1$ISO_WEEKSTARTDATE,'%Y-%m'))
ob1 = ob1%>% group_by(x) %>% summarise(y = sum(INF_ALL)) %>% as.data.frame()
ob1$x = as.Date(paste0(ob1$x,'-01'))

df2 = fread('WHO-COVID-19-global-data.csv')
df2[is.na(df2)] = 0
ob2 = df2[,c('Date_reported','New_cases')]
ob2$x = as.character(format(ob2$Date_reported,'%Y-%m'))
ob2 = ob2 %>% group_by(x) %>% summarise(y = sum(New_cases))
ob2$x = as.Date(paste0(ob2$x,'-01'))

values = c('#2159bf','#7b5378')
if(F){
  p1 = ggplot(ob1, aes(x = x, y = y)) + geom_line(color = '#8b6388')+
    geom_point(color = '#7b5378',size = 0.1) +
    scale_x_date(breaks = seq(as.Date('1990-01-01'), as.Date('2024-11-01'), by="5 years"),
                 minor_breaks = seq(as.Date('1990-01-01'), as.Date('2024-11-01'), by ='1 year'),
                 date_labels = "%Y-%b",
                 expand = c(0, 0)) + 
    xlab('') + ylab('') +
    theme_bw() 
  pdf(paste0("Output/fc_flu.pdf"), width = 4.5, height = 1.4)
  print(p1)
  dev.off()
  
 
  p2 = ggplot(ob2, aes(x = x, y = y)) + geom_line(color = '#3169cf') +
    geom_point(color = '#2159bf',size = 0.5) +
    scale_x_date(breaks = seq(as.Date('1990-01-01'), as.Date('2024-11-01'), by="1 year"),
                 minor_breaks = seq(as.Date('1990-01-01'), as.Date('2024-11-01'), by ='1 year'),
                 date_labels = "%y-%b",
                 expand = c(0, 0)) + 
    scale_y_continuous(breaks = c(0,3e7,6e7,9e7)) +
    xlab('') + ylab('') +
    theme_bw() 
  pdf(paste0("Output/fc_covid.pdf"), width = 2.8, height = 1.4)
  print(p2)
  dev.off()
}

start =  as.Date('2019-08-01')
observed_matrix = data.frame(ob1[ob1$x >start,])
observed_matrix[,3] = 0
for (i in 1:nrow(ob2)) {
  observed_matrix[observed_matrix$x == ob2$x[i],3] = ob2[i,'y']
}
rownames(observed_matrix) = unlist(observed_matrix$x)
observed_matrix = observed_matrix[,-1]
observed_matrix = observed_matrix[,c(2,1)]
voc = c('covid','flu')
colnames(observed_matrix) = voc
# > head(observed_matrix)
# covid    flu
# 2019-09-01     0  13553
# 2019-10-01     0  14855
# 2019-11-01     0  38204
# 2019-12-01     0 144701
# 2020-01-01  2032 177787
# 2020-02-01 76909 171979

x0 = c(0, unlist(ob1[ob1$x == start,'y']))

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
if(F){
  load('fc.rdata')
}
determinant_fun = function(cond = T, ifsimu  = T, n_simu = 1){
  # cond = T; ifsimu  = F; n_simu = 1
  n = 2
  poolday = 30
  nday = 100
  seed_matrix <- matrix(0, nrow=n, ncol=poolday)
  seed_matrix[,1] = unlist(x0/30+1)
  
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
  
  for (j in 1:59) {
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
  for (i in 1:60) {
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
                  x = rep(as.Date(rownames(observed_matrix)) + months(1) - 1,2),
                  group = rep(voc, each = nrow(observed_matrix)))

df2_list = list()

for (n_simu in 1:1000) {
  print(n_simu)
  data = determinant_fun(cond = F, ifsimu  = T, n_simu = n_simu)
  df2_list[[n_simu]] = data$y
}

simu_Onset = data.frame(bind_cols(df2_list))
if(F){
  save(fitlist, Onsets_mat_list, file = 'fc.rdata')
  save(simu_Onset, file = 'fc_plot.rdata')
  load('fc.rdata')
  load('fc_plot.rdata')
  data = determinant_fun(cond = F, ifsimu  = F, n_simu = n_simu)
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



plot_data$date = as.Date('2019-08-31') + plot_data$x
plot_data$group = factor(plot_data$group, levels = voc)


k = 30
pd = 25

dfob$y[dfob$group == 'flu'] = dfob$y[dfob$group == 'flu']*10^2
plot_data[plot_data$group == 'flu',3:5] = plot_data[plot_data$group == 'flu',3:5]*10^2
values = c('#6199ff','#7b5378')
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
  scale_x_date(breaks = seq(as.Date('2018-01-01'), as.Date('2024-11-01'), by="6 months"),
               minor_breaks = seq(as.Date('2018-01-01'), as.Date('2024-11-01'), by ='1 month'),
               date_labels = "%y-%b",
               expand = c(0, 0)) +
  scale_y_continuous('COVID-19 (*1e6)',
    sec.axis = sec_axis(~ . /100, name = "Influenza (*1e4)",
                        labels = c(0,2,4,6), breaks = c(0,2,4,6)*1e4),
    breaks = c(0,2,4,6)*1e6, labels = c(0,2,4,6)
  ) + 
  xlab('') +
  coord_cartesian(xlim = c(as.Date('2019-12-21'), as.Date('2024-07-01'))) +
  theme(legend.position = "none",
        axis.text.y.left = element_text(colour = values[1]),
        axis.text.y.right = element_text(colour = values[2]),
        axis.title.y.left = element_text(color = values[1]),
        axis.title.y.right = element_text(color = values[2]))

pdf(paste0("Output/fc_plot.pdf"), width = 4.5, height = 1.8)
print(p)
dev.off()

parsmat = data.frame()
for (i in 1:length(fitlist)) {
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
parsmat$date = poolday*(parsmat$Tcycle-1) + as.Date('2019-08-31') 

df = parsmat
df$date = as.Date(df$date)
pd = 4

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
  coord_cartesian(xlim = c(as.Date('2019-12-21'),
                           as.Date('2024-07-01'))) +
  # scale_y_continuous(breaks = c(0.3,0.4)) +
  theme_bw() + xlab('') + ylab(expression(beta)) +
  theme(plot.margin = unit(c(0.3,0.6,0,0), "cm"))
pdf(paste0("Output/fc_plot_beta.pdf"), width = 4.3, height = 1.2)
print(p3)
dev.off()
p4 = ggplot() +
  geom_boxplot(data = df[df$group == 'contact',], 
               aes(x = date, y = y, group = date),
               color = alpha('#3d3da1',0.9), 
               fill = alpha('#3d3da1',0.3),
               width = 6,
               outlier.shape = NA) +
  scale_x_date(breaks = seq(as.Date('2018-01-01'), as.Date('2024-11-01'), by="6 months"),
               minor_breaks = seq(as.Date('2018-01-01'), as.Date('2024-11-01'), by ='1 month'),
               date_labels = "%y-%b",
               expand = c(0, 0)) +
  geom_line(data = medians[medians$group == 'contact',], 
            aes(x = date, y = median_y), 
            color = '#3d3da1', 
            linewidth = 0.4) +
  coord_cartesian(xlim = c(as.Date('2019-12-21'),
                           as.Date('2024-07-01'))) +
  theme_bw() + xlab('') + ylab('') +
  theme(plot.margin = unit(c(0,0.6,0,0), "cm"))


p4
pdf(paste0("Output/fc_plot_Q.pdf"), width = 3, height = 1)
print(p4)
dev.off()


if(F){
  data = determinant_fun(cond = F, ifsimu  = T, n_simu = n_simu)
  data$date = as.Date('2019-08-31') + data$x
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


