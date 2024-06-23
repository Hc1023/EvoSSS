rm(list = ls())
library(rstan)
library(ggplot2)
library(scales)
library(ggnewscale)
library(tidyverse)
library(dplyr)
library(RColorBrewer)
df = read.csv('../3_Epidemiological_analysis/VOC_gisaid.csv')

VOC = c("Delta","Alpha","D614G","B 19A","A")
names(VOC) = c("GK","GRY","G","S","L")
df = df[df$GISAID_clade %in% names(VOC),]
df$V = NA
for (i in 1:length(VOC)) {
  df$V[df$GISAID_clade == names(VOC)[i]] = VOC[i]
}
df$V = factor(df$V, levels = VOC)
dfall = df %>% group_by(date) %>%
  summarise(yall = sum(count))

df <- na.omit(df)

full_dates <- as.Date('2020-07-01') + 0:800
# Create a dataframe with all dates
full_df <- expand.grid(date = full_dates, 
                       V = unique(df$V))
df$date = as.Date(df$date)
merged_df <- full_df %>%
  left_join(df, by = c("date", "V")) %>%
  mutate(count = ifelse(is.na(count), 0, count))
voc = unique(merged_df$V)
observed_matrix = data.frame(v1 = merged_df[merged_df$V == 'D614G', 'count'],
                             v2 = merged_df[merged_df$V == 'Alpha', 'count'],
                             v3 = merged_df[merged_df$V == 'Delta', 'count'])
rownames(observed_matrix) = full_dates


fexpect0 = data.frame(x = rep(1:length(full_dates),3), 
                      y = c(observed_matrix[,1],
                            observed_matrix[,2],
                            observed_matrix[,3]),
                      group = factor(rep(c('D614G','Alpha','Delta'), 
                                         each = length(full_dates)),
                                     levels = c('D614G','Alpha','Delta')))
ggplot(data = fexpect0, 
       aes(x = x, y = y, 
           group = group, color = group)) +
  geom_line() +
  scale_y_continuous(trans = 'log10')

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

load('voc_fitlist.rdata')

determinant_simu = function(n_simu){
  poolday = 30
  # The initial cycle - epidemic outbreak
  n= 3
  restrict_coef = 0.7
  
  seed_vec = c(34,rep(0, poolday-1))
  probs <- c(1,1e-3,1e-3)
  
  seed_matrix <- matrix(0, nrow=n, ncol=length(seed_vec))
  for (i in 1:n) {
    seed_matrix[i,] <- seed_vec * probs[[i]]
  }
  seed_mats <- list()
  
  for (i in 1:n) {
    seed_mats[[i]] <- diag(seed_matrix[i,])
  }
  
  nday = 100
  
  Onsets_mat_list = list()
  fit= fitlist[[1]]
  posterior = rstan::extract(fit)
  pars_last = c(posterior$contact[n_simu], posterior$beta1[n_simu],
                posterior$beta2[n_simu], posterior$beta3[n_simu])
  Onsets_mat = simu(seed_mats, 
                    N = seed_vec * pars_last[1] + 1, 
                    poolday, pars = c(pars_last[2], pars_last[3], pars_last[4]), n)
  
  Onsets_mat_list[[1]] = Onsets_mat
  for (j in 1:24) {
    # print(j)
    Onsets_mat = Onsets_mat_list[[j]]
    
    # Generalized extraction of Onset columns for all variants
    Onsets <- list()
    for (i in 1:n) {
      Onsets[[i]] <- Onsets_mat[poolday + 1:poolday, i] + 1e-10
    }
    
    # Mobility control force
    mobility <- rep(0.1, 30) 
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
    pars_last = c(posterior$contact[n_simu], posterior$beta1[n_simu],
                  posterior$beta2[n_simu], posterior$beta3[n_simu])
    Onsets_mat = simu(seed_mats, 
                      N = seed_vec * pars_last[1] + 1, 
                      poolday, pars = c(pars_last[2], pars_last[3], pars_last[4]), n)
    
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
  
  data$color = factor(data$color,levels = levels)
  return(data)
}

df2_list = list()
for (n_simu in 1:100) {
  print(n_simu)
  data = determinant_simu(n_simu)
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

plot_data$date = plot_data$x + as.Date('2020-07-01')
plot_data$group = 'D614G'
plot_data$group[data$color == '2'] = 'Alpha'
plot_data$group[data$color == '3'] = 'Delta'
save(simu_Onset, plot_data, file = 'determinants_plot.rdata')

# dfall = plot_data %>% group_by(date) %>% 
#   summarise(yall = sum(y))

fexpect0$date = as.Date('2020-07-01') + fexpect0$x
plot_data$group = factor(plot_data$group, levels = c('D614G','Alpha','Delta'))
brews = brewer.pal(n=4, name = "Spectral")
values = c(brews[4:2])
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
  scale_y_continuous(trans='log10') +
  coord_cartesian(ylim = c(2,max(data$y))) +
  labs(x = "Date", y = "Proportion") +
  theme_bw() +
  theme(legend.position = "right",
        plot.title = element_text(hjust = 0.5)) +
  scale_x_date(breaks = seq(as.Date('2020-01-01'), as.Date('2022-11-01'), by="6 months"),
               minor_breaks = seq(as.Date('2019-12-01'), as.Date('2022-11-01'), by ='1 month'),
               date_labels = "%y-%b") +
  scale_y_continuous(trans='log10', 
                     breaks = c(1,10,100,1000,10000),
                     labels = c(expression(10^0),expression(10^1),
                                expression(10^2), expression(10^3),
                                expression(10^4))) +
  xlab('') + ylab('Cases') + 
  coord_cartesian(xlim = c(as.Date('2020-06-30'), as.Date('2021-10-31')),
                  ylim = c(1,2*10^4)) +
  theme(legend.position = c(0.2,0.85),
        legend.background = element_rect(color = NA, fill = NA),
        legend.key = element_blank(),
        legend.key.size = unit(0.2, units = 'cm'))

pdf(paste0("Output/voc_plot.pdf"), width = 2.5, height = 1.8)
print(p)
dev.off()

# data_plot$group = factor(data$group, levels = c('Delta','Alpha','D614G'))
data2 = plot_data[,c('Fitted','date','group')]
data2 <- data2 %>%
  group_by(date) %>%
  mutate(p = Fitted/sum(Fitted)) %>%
  as.data.frame()
data2$group = factor(data2$group, levels = c('Delta','Alpha','D614G'))
p2= ggplot() +
  geom_area(data = data2, 
            aes(x = date, y = p, fill = group),
            position = 'fill') +
  scale_fill_manual(name="", breaks = rev(levels(data2$group)),
                    values = alpha(values, 1)) +
  coord_cartesian(xlim = c(as.Date('2020-07-01'), as.Date('2021-10-31')),
                  ylim = c(0,1)) +
  scale_x_date(breaks = seq(as.Date('2020-01-01'), as.Date('2024-11-01'), by="6 months"),
               minor_breaks = seq(as.Date('2019-12-01'), as.Date('2024-11-01'), by ='1 month'),
               date_labels = "%y-%b",
               expand = c(0, 0)) +
  scale_y_continuous(breaks = seq(0,1,0.5), expand = c(0, 0)) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = 'right',
        legend.background = element_rect(color = NA, fill = NA),
        legend.key = element_blank(),
        legend.key.size = unit(0.4, units = 'cm')) +
  xlab('') + ylab('Prevalence') 

pdf(paste0("Output/voc_prevalence.pdf"), width = 3, height = 1.2)
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
                                   posterior$beta3),
                             group = c(rep('contact', length(posterior$contact)),
                                       rep('beta1', length(posterior$beta1)),
                                       rep('beta2', length(posterior$beta2)),
                                       rep('beta3', length(posterior$beta3)))
                             
                  ))
  
}
poolday = 30
parsmat$date = poolday*parsmat$Tcycle + as.Date('2020-07-01')
parsmat$group = factor(parsmat$group)
parsmat$date = factor(parsmat$date)

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
  coord_cartesian(xlim = c(as.Date('2020-07-01'), as.Date('2021-10-31')),
                  ylim = c(0.15,0.3)) +
  theme_bw() + xlab('') + ylab(expression(beta))
pdf(paste0("Output/voc_beta.pdf"), width = 2.5, height = 1.2)
print(p3)
dev.off()
