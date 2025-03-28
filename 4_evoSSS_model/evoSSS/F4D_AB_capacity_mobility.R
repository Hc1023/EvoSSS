rm(list = ls())
library(rstan)
library(ggplot2)
library(scales)
library(ggnewscale)
library(tidyverse)
library(dplyr)
library(RColorBrewer)

df = read.csv('../3_Epidemiological_analysis/F1D_Covid19CasesGISAID.csv')

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
  pars = c(0.379, 0.398, 0.157)
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


poolday = 30
# The initial cycle - epidemic outbreak
seed_vec = matrix(0,2,2)
seed_vec[,1] = round(34*c(0.4,0.6))
seed_mat_I1 = diag(seed_vec[1,])
seed_mat_I2 = diag(seed_vec[2,])
N = rep(32583, 2)
Onsets_mat_list = list()
Onsets_mat = simu(seed_mat_I1, seed_mat_I2, N, poolday)
Onsets_mat_list[[1]] = Onsets_mat

# document_scale = max(Onsets_mat[,2])/max(df$Freq[df$Var1<as.Date('2020-02-01')])
observed_matrix = observed_matrix*28

poolday = 30
nday = 100
expected_matrix = observed_matrix[1:nday,]
expected_matrix[expected_matrix < 0] = 0
expected_matrix = round(expected_matrix)
expected_matrix[(poolday*2+1):nday,] = 0
fexpect = data.frame(x = rep(1:nday,2), 
                     y = c(expected_matrix[1:nday,1],
                           expected_matrix[1:nday,2]),
                     group = factor(rep(c('A','B'), 
                                        each = nday),
                                    levels = c('A','B')))

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

determinant_fun = function(cond = T, ifsimu  = T, n_simu = 1){
  
  n = 2
  poolday = 30
  nday = 100
  pars = c(0.379, 0.398, 0.157)
  # cond = T; ifsimu = F
  # cond = F
  pars_last = c(1/30, 200)
  for (j in 1:24) {
    if(!ifsimu) print(j)
    {
      
      Onsets_mat = Onsets_mat_list[[j]]
      
      # Generalized extraction of Onset columns for all variants
      Onsets <- list()
      for (i in 1:n) {
        Onsets[[i]] <- Onsets_mat[poolday + 1:poolday, i] + 1e-1
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
      
      if(cond){
        stan_data <- list(
          poolday = poolday,
          nday = nday,
          expected_matrix = expected_matrix,
          pars = pars,
          Onset1 = Onsets[[1]],
          Onset2 = Onsets[[2]],
          pars_last = pars_last
        )
        
        # Fit the model
        fit <- stan(file = 'evoSSS/F4D_capacity_mobility.stan', data = stan_data, 
                    iter = 3000, chains = 1, warmup = 2000,
                    verbose = TRUE)
        fitlist[[j]] = fit
      }
      
      
      fit = fitlist[[j]]
      posterior = rstan::extract(fit)
      
      pars_last = c(mean(posterior$mobility), mean(posterior$contact))
      if(ifsimu){
        pars_last = c(posterior$mobility[n_simu],
                      posterior$contact[n_simu])
      }
      # Mobility control force
      mobility <- rep(pars_last[1], 30) 
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
      N = seed_vec * pars_last[2] + 1
      
      Onsets_mat = simu(seed_mats[[1]], seed_mats[[2]], N, poolday)
      
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
  fitlist = list()
  save(fitlist, file = 'evoSSS/F4D_AB_mobility.rdata')
  load('evoSSS/F4D_AB_mobility.rdata')
}


######### plot dynamic ##############
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
  
  if(F){
    save(fitlist, simu_Onset, plot_data, file = 'F4D_evoSSS/AB_mobility.rdata')
    load('evoSSS/F4D_AB_mobility.rdata')
  }
}

values = c(hue_pal()(3)[1], hue_pal()(3)[3])
plot_data = plot_data[plot_data$x>1,]

p = ggplot() +
  geom_point(data = fexpect0, 
             aes(x = date, y = y, 
                 group = group, color = group),
             size = 0.4, shape = 16) +
  geom_ribbon(data = plot_data, 
              aes(x = date, group = group, 
                  ymin = LowerCI/28, ymax = UpperCI/28, fill = group)) +  # Confidence interval
  geom_line(data = plot_data, 
            aes(x = date, y = Fitted/28, 
                group = group, color = group), linewidth = 1) +
  scale_color_manual(name="",
                     values = alpha(values, 0.7)) +
  scale_fill_manual(name="",
                    values = alpha(values, 0.3)) +
  scale_y_continuous(trans='log10') +
  coord_cartesian(ylim = c(2,max(plot_data$Fitted))) +
  labs(x = "Date", y = "Proportion") +
  theme_bw() +
  scale_y_continuous(trans='log10', 
                     breaks = c(1,10,100,1000,10000),
                     labels = c(expression(10^0),expression(10^1),
                                expression(10^2), expression(10^3),
                                expression(10^4))) +
  xlab('') + ylab('Cases') + 
  scale_x_date(breaks = seq(as.Date('2020-01-01'), as.Date('2022-11-01'), by="6 months"),
               minor_breaks = seq(as.Date('2019-12-01'), as.Date('2022-05-01'), by ='1 month'),
               date_labels = "%y-%b", expand = c(0, 0)) +
  coord_cartesian(xlim = c(as.Date('2020-01-01'), as.Date('2021-10-31')),
                  ylim = c(1,4*10^4)) +
  theme(legend.position = 'none',
        panel.grid.minor = element_blank())

pdf(paste0("Output/F4D_AB_mobility.pdf"), width = 2.5, height = 1.6)
print(p)
dev.off()


######### parameter dynamic ##############
parsmat = data.frame()
for (i in 1:22) {
  posterior = rstan::extract(fitlist[[i]])
  parsmat = rbind(parsmat, 
                  data.frame(Tcycle = i,
                             y = c(posterior$mobility, posterior$contact),
                             group = c(rep('mobility', length(posterior$contact)),
                                       rep('contact', length(posterior$contact)))
                             
                  ))
  
}

poolday = 30
parsmat$date = as.Date(poolday*(parsmat$Tcycle) + 1 + as.Date('2019-12-31'))

medians <- parsmat %>%
  group_by(date, group) %>%
  summarize(median_y = median(y))

values = c('blue', 'red')
p1 = ggplot() +
  geom_boxplot(data = parsmat[parsmat$group == 'mobility',], 
               aes(x = date, y = y, group = date),
               color = alpha(values[1],0.7), 
               fill = alpha(values[1],0.3),
               width = 6,
               outlier.shape = NA) +
  geom_line(data = medians[medians$group == 'mobility',], 
            aes(x = date, y = median_y), 
            color = values[1], 
            linewidth = 0.4) +
  scale_x_date(breaks = seq(as.Date('2020-01-01'), as.Date('2022-11-01'), by="6 months"),
               minor_breaks = seq(as.Date('2019-12-01'), as.Date('2022-05-01'), by ='1 month'),
               date_labels = "%y-%b", expand = c(0, 0)) +
  coord_cartesian(xlim = c(as.Date('2020-01-01'), as.Date('2021-10-31'))) +
  theme_bw() + xlab('') + ylab('m') + 
  theme(panel.grid.minor = element_blank())

# > median(parsmat[parsmat$group == 'mobility','y'])
# [1] 0.002326484

p2 = ggplot() +
  geom_boxplot(data = parsmat[parsmat$group == 'contact',], 
               aes(x = date, y = log2(y), group = date),
               color = alpha(values[2],0.7), 
               fill = alpha(values[2],0.3),
               width = 6,
               outlier.shape = NA) +
  geom_line(data = medians[medians$group == 'contact',], 
            aes(x = date, y = log2(median_y)), 
            color = values[2], 
            linewidth = 0.4) +
  scale_x_date(breaks = seq(as.Date('2020-01-01'), as.Date('2022-11-01'), by="6 months"),
               minor_breaks = seq(as.Date('2019-12-01'), as.Date('2022-05-01'), by ='1 month'),
               date_labels = "%y-%b", expand = c(0, 0)) +
  coord_cartesian(xlim = c(as.Date('2020-01-01'), as.Date('2021-10-31'))) +
  theme_bw() + xlab('') + ylab('Q') + 
  scale_y_continuous(breaks = log2(c(500, 1000, 2000, 4000)),
                     labels = c(500, 1000, 2000, 4000)) +
  theme(panel.grid.minor = element_blank())

pdf(paste0("Output/F4D_AB_mobility_par.pdf"), width = 2.6, height = 1)
print(p1)
print(p2)
dev.off()


######### one trial ##############
if(F){
  data = determinant_fun(cond = F, ifsimu  = F, n_simu = n_simu)
  data$date = as.Date('2019-12-31') + data$x
  data$group = 'A'
  data$group[data$color == '2'] = 'B'
  data$group = factor(data$group, levels = c('A','B'))
  
  ggplot() +
    geom_point(data = fexpect0, 
               aes(x = date, y = y/28, 
                   group = group, color = group)) +
    geom_line(data = data, 
              aes(x = date, y = y/28, group = group, color= group)) +
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
}

######### tranceplot ##############
if(F){
  for (i in 1:length(fitlist)) {
    wd = paste0('Output/traceplot/AB_mobility_', as.character(i), '.png')
    png(file = wd, width = 2.5, height = 2, unit = 'in', res = 600)
    p = traceplot(fitlist[[i]]) + 
      ggtitle(as.character(i)) + ylab('') +
      theme(legend.position = 'none') + 
      scale_x_continuous(breaks = seq(2000,3000,500),
                         labels = seq(0,1000,500))
    print(p)
    dev.off()
  }
  
}

