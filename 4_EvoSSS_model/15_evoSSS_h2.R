rm(list = ls())
library(rstan)
library(ggplot2)
library(scales)
library(ggnewscale)
library(tidyverse)
library(dplyr)
library(RColorBrewer)

df = read.csv('../3_Epidemiological_analysis/Covid19CasesGISAID.csv')
load('AB_constant_beta.rdata')
pars = c(0.398, 0.398, 0.157)
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

fexpect0 = data.frame(y = c(observed_matrix$v1,observed_matrix$v2,
                            observed_matrix$v3),
                      x = rep(1:nrow(observed_matrix),2),
                      group = rep(c('A','B'), 
                                  each = nrow(observed_matrix)))
fexpect0$date = as.Date('2019-12-31') + fexpect0$x
fexpect0$group = factor(fexpect0$group, levels = c('A','B'))

update_fun = function(pars, states_old, N){
  
  S <- states_old[,1]
  I1 <- states_old[,2]
  I2 <- states_old[,3]
  
  beta1 = pars[1]
  beta2 = pars[2]
  gamma = pars[3]
  
  S_new = S - beta1*S*I1/N - beta2*S*I2/N
  I1_new = I1 + beta1*S*I1/N - gamma*I1
  I2_new = I2 + beta2*S*I2/N - gamma*I2
  Onset1 = beta1*S*I1/N
  Onset2 = beta2*S*I2/N
  
  return(cbind(S_new, I1_new, I2_new, Onset1, Onset2))
}

simu <- function(seed_mat_I1, seed_mat_I2, N, poolday) {
  
  # Initial conditions
  I1_old = seed_mat_I1[1,]
  I2_old = seed_mat_I2[1,]
  S_old = N - I1_old - I2_old
  states_old = cbind(S_old, I1_old, I2_old)
  
  # Simulate the dynamics over 
  ndays = 200
  
  Onsets_mat <- matrix(0, ndays, 2)
  
  for (t in 2:ndays) {
    states_old = update_fun(pars = pars, states_old = states_old, N = N)
    if(t <= nrow(seed_mat_I1)){
      states_old[t,2:3] = states_old[t,2:3] + 
        c(seed_mat_I1[t,t], seed_mat_I2[t,t])
    }
    Onsets_mat[t,] = c(sum(states_old[,4]), sum(states_old[,5]))
  }
  
  return(Onsets_mat)
}




determinant_fun = function(ifsimu  = T, n_simu = 1, require_dfplot_simu = F){
  poolday = 30
  # The initial cycle - epidemic outbreak
  seed_vec = matrix(0,2,2)
  seed_vec[,1] = round(34*c(0.4,0.6))
  seed_vec[,1] = rmultinom(1, size = 34, prob = c(0.4,0.6))
  seed_mat_I1 = diag(seed_vec[1,])
  seed_mat_I2 = diag(seed_vec[2,])
  N = rep(32583, 2)
  Onsets_mat_list = list()
  Onsets_mat = simu(seed_mat_I1, seed_mat_I2, N, poolday)
  Onsets_mat_list[[1]] = Onsets_mat
  poolday = 30
  nday = 100
  n = 2
  pars = c(0.398, 0.398, 0.157)
  
  
  for (j in 1:24) {
    # print(j)
    {
      
      Onsets_mat = Onsets_mat_list[[j]]
      
      # Generalized extraction of Onset columns for all variants
      Onsets <- list()
      for (i in 1:n) {
        Onsets[[i]] <- Onsets_mat[poolday + 1:poolday, i] + 1e-3
      }
      
      fit = fitlist[[j+1]]
      posterior = rstan::extract(fit)
      
      pars_last = c(mean(posterior$mobility), mean(posterior$contact))
      if(ifsimu){
        pars_last = c(posterior$mobility[n_simu],
                      posterior$contact[n_simu])
      }
      # Mobility control force
      mobility <- rep(pars_last[1], 30) 
      Mobility_matrix <- diag(mobility)
      seed_vec =  (Onsets[[1]] + Onsets[[2]]) %*% 
        Mobility_matrix %>% as.numeric()
      p = Onsets[[1]]/(Onsets[[1]] + Onsets[[2]])
      seed_matrix = sapply(1:length(p), function(x){
        rmultinom(1, ceiling(seed_vec[x]), c(p[x],1-p[x]))
      })
      
      # Create seed_matrices for each variant
      seed_mats <- list()
      
      for (i in 1:n) {
        seed_mats[[i]] <- diag(seed_matrix[i,])
      }
      N = seed_vec * pars_last[2] + 1
      
      Onsets_mat = simu(seed_mats[[1]], seed_mats[[2]], N, poolday)
    }
    
    Onsets_mat_list[[j+1]] = Onsets_mat
  }
  
  dfplot_simu = data.frame()
  for (i in 1:25) {
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
  
  data$color = factor(data$color, levels = levels)
  
  if(require_dfplot_simu){
    return(dfplot_simu)
  }
  return(data)
}

if(T){
  
  df2_list = list()
  for (n_simu in 1:100) {
    print(n_simu)
    data = determinant_fun(ifsimu  = T, n_simu = n_simu, require_dfplot_simu = F)
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
  plot_data$group[plot_data$V == '2'] = 'B'
  plot_data$group = factor(plot_data$group, levels = c('A','B'))
  plot_data = plot_data[plot_data$x>1,]
  plot_data_sampling = plot_data
  if(F){
    save(plot_data_sampling, file = 'h2.rdata')
    load('h2.rdata')
  }

  values = c(hue_pal()(3)[1], hue_pal()(3)[3])
  
  p = ggplot() +
    geom_point(data = fexpect0, 
               aes(x = date, y = y, 
                   group = group, color = group),
               size = 0.4, shape = 16) +
    geom_ribbon(data = plot_data_sampling, 
                aes(x = date, group = group, 
                    ymin = LowerCI/28, ymax = UpperCI/28, fill = group)) + 
    geom_line(data = plot_data_sampling, 
              aes(x = date, y = Fitted/28, 
                  group = group, color = group), linewidth = 1) +
    scale_color_manual(name="Variant",
                       values = alpha(values, 0.7)) +
    scale_fill_manual(name="Variant",
                      values = alpha(values, 0.3)) +
    scale_y_continuous(trans='log10',
                       breaks = c(1,10,100,1000,10000),
                       labels = c(expression(10^0),expression(10^1),
                                  expression(10^2), expression(10^3),
                                  expression(10^4))) +
    labs(x = "Date", y = "Proportion") +
    theme_bw() +
    theme(legend.position = "right",
          legend.key.size = unit(0.2,'cm'),
          legend.spacing = unit(0.0,'cm'),
          legend.text = element_text(size = 8),
          legend.title = element_text(size = 10,
                                      margin = margin(2, 0, 2, 0)),
          legend.margin = margin(0, 0, 0, 0),
          panel.grid.minor = element_blank()) +
    scale_x_date(breaks = seq(as.Date('2020-01-01'), as.Date('2022-11-01'), by="6 months"),
                 minor_breaks = seq(as.Date('2019-12-01'), as.Date('2022-11-01'), by ='1 month'),
                 date_labels = "%y-%b") +
    xlab('') + ylab('Cases') + 
    coord_cartesian(xlim = c(as.Date('2020-01-01'), as.Date('2021-10-31')),
                    ylim = c(1,2*10^4))
  
  p
  pdf(paste0("Output/evoSSS_h2.pdf"), width = 3, height = 1.5)
  print(p)
  dev.off()
}


