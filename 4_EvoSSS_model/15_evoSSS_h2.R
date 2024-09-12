rm(list = ls())
library(rstan)
library(ggplot2)
library(scales)
library(ggnewscale)
library(tidyverse)
library(dplyr)
library(RColorBrewer)

df = read.csv('../3_Epidemiological_analysis/Covid19CasesGISAID.csv')
load('h2_hotspots.rdata')
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
                             v2 = merged_df[merged_df$V == 'B', 'count']) * 28
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



fitlist = list()
determinant_fun = function(cond = F, ifsimu  = T, n_simu = 1, 
                           h = 2, seed = 6){
  # ifsimu  = F; cond = T; n_simu = 1
  
  poolday = 30
  # The initial cycle - epidemic outbreak
  seed_vec = matrix(0,2,2)
  seed_vec[,1] = round(34*c(0.3,0.7))
  seed_vec[,1] = rmultinom(1, size = 34, prob = c(0.3,0.7))
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
  pars_last = 10000
  for (j in 1:24) {
    if(cond) print(j)
    {

      Onsets_mat = Onsets_mat_list[[j]]
      
      # Generalized extraction of Onset columns for all variants
      Onsets <- list()
      for (i in 1:n) {
        Onsets[[i]] <- Onsets_mat[poolday + 1:poolday, i] + 1e-3
      }
      expected_matrix = observed_matrix[poolday*j+1:nday,] - Onsets_mat[poolday+1:nday,] 
      expected_matrix[expected_matrix<0] = 0
      expected_matrix = round(expected_matrix)
      expected_matrix[(2*poolday + 1):nday,] = 0
      
      fexpect = data.frame(x = rep(1:nday,2), 
                           y = c(expected_matrix[1:nday,1],
                                 expected_matrix[1:nday,2]),
                           group = factor(rep(c('A','B'), 
                                              each = nday),
                                          levels = c('A','B')))
      
      seed_matrix = matrix(0,2,poolday)
      seed_vec =  c(rep(seed,h), rep(0, poolday-h))
      p = Onsets[[1]]/(Onsets[[1]] + Onsets[[2]])
      
      if(ifsimu & j<=10){
        seed_matrix[,1:h] = sapply(1:h, function(x){
          rmultinom(1, ceiling(seed_vec[x]), c(p[x],1-p[x]))
        })
      }else{
        seed_matrix[,1:h] = sapply(1:h, function(x){
          seed_vec[x] * c(p[x],1-p[x])
        })
      }
      # Create seed_matrices for each variant
      seed_mats = list()
      
      for (i in 1:n) {
        seed_mats[[i]] <- diag(seed_matrix[i,])
      }
      if(cond){
        stan_data <- list(
          poolday = poolday,
          nday = nday,
          expected_matrix = expected_matrix,
          pars = pars,
          seed_mat_I1 = seed_mats[[1]],
          seed_mat_I2 = seed_mats[[2]],
          seed_vec = seed_vec,
          h = h,
          pars_last = pars_last
        )
        
        # Fit the model
        fit <- stan(file = 'h2.stan', data = stan_data, 
                    iter = 3000, chains = 1, warmup = 2000,
                    verbose = TRUE)
        fitlist[[j+1]] = fit
      }
      fit = fitlist[[j+1]]
      posterior = rstan::extract(fit)
    
      pars_last = mean(posterior$contact)
      if(ifsimu){
        pars_last = posterior$contact[n_simu]
      }
    
      Onsets_mat = simu(seed_mats[[1]], seed_mats[[2]], 
                        seed_vec * pars_last + 1, poolday)
      
      if(!ifsimu){
        fonset = data.frame(x = rep(1:nday,2), 
                            y = c(Onsets_mat[1:nday,1],
                                  Onsets_mat[1:nday,2]),
                            group = factor(rep(c('A','B'), each = nday),
                                           levels = c('A','B')))
        
        ggplot() +
          geom_point(data = fexpect, 
                     aes(x = x, y = y, group = group, color = group)) +
          geom_line(data = fonset,
                    aes(x = x, y = y, group = group, color = group)) 
      }
    }
    
    Onsets_mat_list[[j+1]] = Onsets_mat
  }
  
  if(cond == T){
    return(fitlist)
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
  
  return(data)
}



if(F){
  fitlist_s6_10 = determinant_fun(cond = T, ifsimu = F, h = 10, seed = 6)
  fitlist = fitlist_s6_10
  fitlist_s20_10 = determinant_fun(cond = T, ifsimu = F, h = 10, seed = 20)
  fitlist = fitlist_s20_10
  fitlist_s6_30 = determinant_fun(cond = T, ifsimu = F, h = 30, seed = 6)
  fitlist = fitlist_s6_30
  fitlist_s20_30 = determinant_fun(cond = T, ifsimu = F, h = 30, seed = 20)
  fitlist = fitlist_s20_30
  
  save(fitlist_s6_10, fitlist_s20_10, 
       fitlist_s6_30, fitlist_s20_30,
       file = 'h2_hotspots.rdata')
  load('h2_hotspots.rdata')
  
}
if(T){
  h = 10; seed = 6; fitlist = fitlist_s6_10
  h = 10; seed = 20; fitlist = fitlist_s20_10
  h = 30; seed = 6; fitlist = fitlist_s6_30
  h = 30; seed = 20; fitlist = fitlist_s20_30
  
  df2_list = list()
  for (n_simu in 1:1000) {
    print(n_simu)
    data = determinant_fun(cond = F, ifsimu  = T, n_simu = n_simu, 
                           h = h, seed = seed)
    df2_list[[n_simu]] = data$y
  }
  {
    simu_Onset = data.frame(bind_cols(df2_list))
    ci_lower1 <- apply(simu_Onset, 1, quantile, probs = 0.25, na.rm = T)
    ci_upper1 <- apply(simu_Onset, 1, quantile, probs = 0.75, na.rm = T)
    ci_lower2 <- apply(simu_Onset, 1, quantile, probs = 0.025, na.rm = T)
    ci_upper2 <- apply(simu_Onset, 1, quantile, probs = 0.975, na.rm = T)
    plot_data <- data.frame(
      x = data$x,
      V = data$color,
      # Observed = observed_cases,
      Fitted = rowMeans(simu_Onset),
      LowerCI1 = ci_lower1,
      UpperCI1 = ci_upper1,
      LowerCI2 = ci_lower2,
      UpperCI2 = ci_upper2
    )
    plot_data$date = plot_data$x + as.Date('2019-12-31')
    plot_data$group = 'A'
    plot_data$group[plot_data$V == '2'] = 'B'
    plot_data$group = factor(plot_data$group, levels = c('A','B'))
    plot_data = plot_data[plot_data$x>1,]
  }

  if(F){
    plot_data_s6_10 = plot_data
    plot_data_s20_10 =  plot_data
    plot_data_s6_30 = plot_data
    plot_data_s20_30 =  plot_data
    save(plot_data_s6_10, plot_data_s20_10, 
         plot_data_s6_30, plot_data_s20_30, 
         file = 'h2_hotspots_plot.rdata')
    load('h2_hotspots_plot.rdata')
    
    plot_data = plot_data_s6_10
    plot_data = plot_data_s20_10
    plot_data = plot_data_s6_30
    plot_data = plot_data_s20_30
  }

  values = c(hue_pal()(3)[1], hue_pal()(3)[3])
  getp = function(plot_data){
    p = ggplot() +
      geom_point(data = fexpect0, 
                 aes(x = date, y = y/28, 
                     group = group, color = group),
                 size = 0.1, shape = 16) +
      geom_ribbon(data = plot_data, 
                  aes(x = date, group = group, 
                      ymin = LowerCI2/28, ymax = UpperCI2/28, 
                      fill = group)) + 
      geom_ribbon(data = plot_data, 
                  aes(x = date, group = group, 
                      ymin = LowerCI1/28, ymax = UpperCI1/28, 
                      fill = group), alpha = 0.5) + 
      geom_line(data = plot_data, 
                aes(x = date, y = Fitted/28, 
                    group = group, color = group), linewidth = 0.2) +
      scale_color_manual(name="Variant",
                         values = alpha(values, 0.7)) +
      scale_fill_manual(name="Variant",
                        values = alpha(values, 0.3)) +
      scale_y_continuous(trans='log10',
                         breaks = c(1,100,10000),
                         labels = c(expression(10^0), expression(10^2), expression(10^4))) +
      labs(x = "Date", y = "Proportion") +
      theme_bw() +
      theme(legend.position = "none",
            legend.key.size = unit(0.2,'cm'),
            legend.spacing = unit(0.0,'cm'),
            legend.text = element_text(size = 8),
            legend.title = element_text(size = 10,
                                        margin = margin(2, 0, 2, 0)),
            legend.margin = margin(0, 0, 0, 0),
            panel.grid.minor = element_blank()) +
      scale_x_date(breaks = seq(as.Date('2020-01-01'), as.Date('2022-11-01'), by="12 months"),
                   minor_breaks = seq(as.Date('2019-12-01'), as.Date('2022-11-01'), by ='1 month'),
                   date_labels = "%y-%b") +
      xlab('') + ylab('Cases') + 
      coord_cartesian(xlim = c(as.Date('2020-01-01'), as.Date('2021-10-31')),
                      ylim = c(1,2*10^4))
    
    return(p)
  }
  p1 = getp(plot_data_s6_10)
  p2 = getp(plot_data_s20_10)
  p3 = getp(plot_data_s6_30)
  p4 = getp(plot_data_s20_30)

  pdf(paste0("Output/evoSSS_h2.pdf"), width = 1.6, height = 1.2)
  print(p1)
  print(p2)
  print(p3)
  print(p4)
  dev.off()
}


