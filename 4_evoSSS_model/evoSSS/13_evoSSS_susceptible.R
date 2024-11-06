rm(list = ls())
library(rstan)
library(ggplot2)
library(scales)
library(ggplot2)
library(scales)
library(ggnewscale)
library(tidyverse)
library(dplyr)

load('evoSSS_chain.rdata')

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

simu <- function(seed_mat_I1, seed_mat_I2, N, poolday, pars) {
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
options(dplyr.summarise.inform=F)
df2_list = list()
for (n_simu in 1:100) {
  print(n_simu)
  # pars = c(0.379, 0.398, 0.157)
  pars = c(0.398, 0.398, 0.157)
  poolday = 30
  # The initial cycle - epidemic outbreak
  seed_vec = matrix(0,2,2)
  seed_vec[,1] = round(34*c(0.4,0.6))
  seed_vec[,1] = rmultinom(1, size = 34, prob = c(0.4,0.6))
  seed_mat_I1 = diag(seed_vec[1,])
  seed_mat_I2 = diag(seed_vec[2,])
  N = rep(32583, 2)
  Onsets_mat_list = list()
  Onsets_mat = simu(seed_mat_I1, seed_mat_I2, N, poolday, pars)
  Onsets_mat_list[[1]] = Onsets_mat
  
  for (n in 1:24) {
    Onsets_mat = Onsets_mat_list[[n]]
    Onset1 = Onsets_mat[poolday + 1:poolday, 1]+1 # avoid 0
    Onset2 = Onsets_mat[poolday + 1:poolday, 2]+1 # avoid 0
    
    mobility = rep(0.01,30) # Mobility: Control force
    Mobility_matrix = diag(mobility)
    
    seed_vec =  (Onset1 + Onset2) %*% Mobility_matrix %>% as.numeric()
    p = Onset1/(Onset1 + Onset2)
    seed_matrix = sapply(1:length(p), function(x){
      rmultinom(1, ceiling(seed_vec[x]), c(p[x],1-p[x]))
    })
    # seed_matrix = rbind(seed_vec * p, seed_vec * (1-p))
    seed_mat_I1 = diag(seed_matrix[1,])
    seed_mat_I2 = diag(seed_matrix[2,])
    fit = fitlist[[n]]
    
    posterior = rstan::extract(fit)
    # contact = mean(posterior$contact)
    contact = posterior$contact[n_simu] #1:5000
    N = seed_vec * contact + 1
    Onsets_mat = simu(seed_mat_I1, seed_mat_I2, N, poolday, pars)
    Onsets_mat_list[[n+1]] = Onsets_mat
  }
  
  # calibration
  for (i in 1:length(Onsets_mat_list)) {
    Onsets_mat_list[[i]] = Onsets_mat_list[[i]]/28
  }
  
  
  dfplot_simu = data.frame()
  for (i in 1:25) {
    Onsets_mat = Onsets_mat_list[[i]]
    n = i-1
    
    dfplot_simu1 = data.frame(x = rep(1:nrow(Onsets_mat)+poolday*n,2),
                              y = c(Onsets_mat[,1],Onsets_mat[,2]),
                              group = rep(paste0(c('A','B'),n), 
                                          each = nrow(Onsets_mat)),
                              color = rep(c('A','B'), 
                                          each = nrow(Onsets_mat)))
    dfplot_simu = rbind(dfplot_simu, dfplot_simu1)
  }
  
  df2 = dfplot_simu %>% group_by(x, color) %>% 
    summarise(y = sum(y)) %>%
    as.data.frame()
  df2_list[[n_simu]] = df2$y
  
}

simu_Onset = data.frame(bind_cols(df2_list))
ci_lower <- apply(simu_Onset, 1, quantile, probs = 0.025, na.rm = T)
ci_upper <- apply(simu_Onset, 1, quantile, probs = 0.975, na.rm = T)



plot_data <- data.frame(
  x = df2$x,
  V = df2$color,
  # Observed = observed_cases,
  Fitted = rowMeans(simu_Onset),
  LowerCI = ci_lower,
  UpperCI = ci_upper
)


plot_data$x = plot_data$x + as.Date('2019-12-31')

values = c(hue_pal()(3)[1], hue_pal()(3)[3])

date_vec = as.Date('2019-12-31')+1:700
df = read.csv('../3_Epidemiological_analysis/Covid19CasesGISAID.csv')
df$Var1 = as.Date(df$Var1)
df = df[df$Var1 %in% date_vec,]
dfobserve = data.frame(
  x = rep(date_vec,2),
  y = c(df[df$Mutations == 'Lineage A','Freq'], 
        df[df$Mutations == 'Lineage B','Freq']),
  group = c(rep('A', 700),rep('B', 700))
)

p = ggplot() +
  geom_point(data = dfobserve, 
             aes(x = x, y = y, 
                 group = group, color = group),
             size = 0.4, shape = 16) +
  geom_ribbon(data = plot_data, 
              aes(x = x, group = V, 
                  ymin = LowerCI, ymax = UpperCI, fill = V)) +  # Confidence interval
  geom_line(data = plot_data, 
            aes(x = x, y = Fitted, 
                group = V, color = V), linewidth = 1) +
  scale_y_continuous(trans='log10', 
                     breaks = c(1,10,100,1000,10000),
                     labels = c(expression(10^0),expression(10^1),
                                expression(10^2), expression(10^3),
                                expression(10^4))) +
  theme_bw() +
  scale_color_manual(name="Variant",
                     labels=c("A", "B"),
                     values = alpha(values, 0.6)) +
  scale_fill_manual(name="Variant",
                    labels=c("A", "B"),
                    values = alpha(values, 0.6)) +
  scale_x_date(breaks = seq(as.Date('2020-01-01'), as.Date('2022-11-01'), by="6 months"),
               minor_breaks = seq(as.Date('2019-12-01'), as.Date('2022-11-01'), by ='1 month'),
               date_labels = "%y-%b") +
  xlab('') + ylab('Cases') + 
  coord_cartesian(xlim = c(as.Date('2019-12-31'), as.Date('2021-10-31')),
                  ylim = c(1,2*10^4)) +
  theme(legend.position = 'none')
p
pdf(paste0("Output/evoSSS_h2.pdf"), width = 2.5, height = 1.8)
print(p)
dev.off()