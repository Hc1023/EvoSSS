rm(list = ls())
library(ggplot2)   
library(dplyr)
library(tidyr)
library(scales)
library(deSolve)   
library(rstan)
library(zoo)

load('evoSIR/F3E_intra_competition.rdata')
df_all = read.csv("../2_Experiment/S5A_vcf_231228.csv")
sampleid = read.csv("../2_Experiment/S4B_sampleid.csv")
sampleid[, 2][sampleid[, 2] == ""] <- NA
sampleid[,2] = na.locf(sampleid[,2])

sampleid1 = sampleid[sampleid$original != '',c(1,2,3)]
sampleid2 = sampleid[,-2]

myfun = function(x){
  v = sampleid[sampleid$cell == x,-3]
  v = v[11:28,-c(3,4)]
  colnames(v) = c('Strain', '0', '24', '48', '72', '1st', '2nd')
  v$Replicate = rep(c('1', '2','3'), times = 6)
  gathered_data <- v %>%
    gather(key = "Time", value = "Value", -Strain, -Replicate)
  gathered_data[,c('C8782T','T28144C')] = NA
  for (i in 1:nrow(gathered_data)) {
    df = df_all[df_all$sample== gathered_data[i,4],]
    df$MUT = paste0(df$REF, df$POS, df$ALT)
    gathered_data[i,5] <- ifelse(sum(df$MUT == 'C8782T') == 0, 0, 
                                 df[df$MUT == 'C8782T', 'percent'])
    gathered_data[i,6] <- ifelse(sum(df$MUT == 'T28144C') == 0, 0, 
                                 df[df$MUT == 'T28144C', 'percent'])
  }
  df <- gathered_data[,-4] %>%
    gather(key = "Mutations", value = "VAF", -Strain, -Replicate, -Time)
  
  df$Time <- factor(df$Time, colnames(v)[2:7])
  df$Strain <- factor(df$Strain, unique(df$Strain)[c(1,3,5,2,4,6)])
  return(df)
}

viral_model <- function(t, state, parameters) {
  V1 <- state[1]
  V2 <- state[2]
  
  r1 <- parameters["r1"]
  r2 <- parameters["r2"]
  K <- parameters["K"]
  alpha12 <- parameters["alpha12"]
  alpha21 <- parameters["alpha21"]
  mu <- parameters["mu"]
  
  dV1dt <- r1 * V1 * (1 - (V1 + alpha12 * V2) / K) - mu * V1
  dV2dt <- r2 * V2 * (1 - (V2 + alpha21 * V1) / K) + mu * V1
  
  return(list(c(dV1dt, dV2dt)))
}

withinhost_fun = function(param_sets, state = c(V1 = 1, V2 = 0)){
  # Initial state and time sequence
  times <- seq(0, 100, by = 1)
  # Run simulations for each parameter set
  results <- lapply(seq(nrow(param_sets)), function(i) {
    params <- unlist(param_sets[i, ])
    out <- ode(y = state, times = times, func = viral_model, parms = params)
    out_df <- as.data.frame(out)
    return(out_df)
  })
  
  # Combine all dataframes into one
  combined_results <- bind_rows(results)
  
  return(combined_results)
  
}

################ Calu ###########

# Calu-3
# > result1$coefficients[2,]
# Estimate   Std. Error      t value     Pr(>|t|) 
# 1.911578e-01 1.103223e-02 1.732722e+01 3.736368e-16 
# > result2$coefficients[2,]
# Estimate   Std. Error      t value     Pr(>|t|) 
# 2.272174e-01 6.407589e-03 3.546066e+01 5.932948e-42 
df = myfun('Calu')
# df =df[df$Time %in% unique(df$Time)[1:4],]
df = df %>%
  group_by(Time, Strain, Replicate) %>%
  summarize(VAF = mean(VAF)) 
median_df <- df %>%
  group_by(Time) %>%
  summarize(median_VAF = median(VAF),
            sd = sd(VAF)) %>%
  as.data.frame()

stan_data <- list(
  V10 = median_df[1,2],                   
  V20 = 1-median_df[1,2],  
  T = 72+1,
  alpha21 = 1,
  mu = 1e-6,
  r1 = 1.911578e-01,
  r2 = 2.272174e-01,
  ratios = median_df$median_VAF[2:4],
  ratios_sd = median_df$sd[2:4]
)

if(F){
  fit <- stan(file = 'F3E_intra_competition.stan', data = stan_data, 
              iter = 3000, chains = 1, 
              warmup = 2000, verbose = TRUE)
  fit1 = fit
  # pairs(fit)
}

fit = fit1
posterior = rstan::extract(fit)
param_sets1 <- expand.grid(r1 = stan_data$r1, 
                           r2 = stan_data$r2, 
                           K = posterior$K, 
                           alpha12 = 1, alpha21 = 1, 
                           mu = stan_data$mu)
param_sets1$alpha12 = 1-posterior$deltaa
combined_results = withinhost_fun(param_sets1, 
                                  state = c(stan_data$V10, stan_data$V20))
colnames(combined_results) = c('time', 'V1', 'V2')
total_v = combined_results$V1 + combined_results$V2
combined_results$V1 = combined_results$V1/total_v 
combined_results$V2 = combined_results$V2/total_v
dfribbon = combined_results %>% group_by(time) %>% 
  summarise(meanV1 = mean(V1),
            minV1 = quantile(V1, probs = 0.025), 
            maxV1 = quantile(V1, probs = 0.975),
            meanV2 = mean(V2),
            minV2 = quantile(V2, probs = 0.025), 
            maxV2 = quantile(V2, probs = 0.975)) %>%
  as.data.frame()
dfribbon = rbind(dfribbon[,c(1:4)], setNames(dfribbon[,c(1,5:7)], names(dfribbon[,c(1:4)])))
dfribbon$group = c(rep('V1', 101),rep('V2', 101))

values = c(hue_pal()(3)[1], hue_pal()(3)[3])

p1 = ggplot() + 
  geom_point(data = median_df[1:4,], 
             aes(x = as.numeric(as.character(Time)), 
                 y = median_VAF), 
             size = 2, color = values[1]) +
  geom_errorbar(data = median_df[1:4,], 
                aes(x = as.numeric(as.character(Time)), 
                    ymin = median_VAF - sd, ymax = median_VAF + sd), 
                width = 2, color = values[1]) +
  geom_line(data = dfribbon, 
            aes(x = time, y = meanV1,
                group = group, color = group)) +
  geom_ribbon(data = dfribbon, 
              aes(x = time, ymin = minV1, ymax = maxV1,
                  group = group, fill = group),
              alpha = 0.3) +
  theme_bw() + 
  scale_color_manual(values = values) +
  scale_fill_manual(values = values) +
  xlab('Time unit') + 
  ylab('') +
  theme(legend.position = "none") +
  scale_x_continuous(minor_breaks = seq(0 , 72, 24),
                     breaks = c(0,24,48,72)) +
  scale_y_continuous(limits = c(0, 1),
                     minor_breaks = seq(0 , 1, 0.25),
                     breaks = c(0,0.5,1),
                     labels = c('0.0','0.5','1.0')) +
  coord_cartesian(xlim = c(0,72))

p1

################ Vero ###########

# Vero
# > result1$coefficients[2,]
# Estimate   Std. Error      t value     Pr(>|t|) 
# 2.611061e-01 2.424562e-02 1.076921e+01 2.821663e-09 
# > result2$coefficients[2,]
# Estimate   Std. Error      t value     Pr(>|t|) 
# 3.397670e-01 2.035067e-02 1.669561e+01 4.344274e-19 

df = myfun('Vero')
# df =df[df$Time %in% unique(df$Time)[1:4],]
df = df %>%
  group_by(Time, Strain, Replicate) %>%
  summarize(VAF = mean(VAF)) 
median_df <- df %>%
  group_by(Time) %>%
  summarize(median_VAF = median(VAF),
            sd = sd(VAF)) %>%
  as.data.frame()

stan_data <- list(
  V10 = median_df[1,2],                   
  V20 = 1-median_df[1,2],  
  T = 72+1,
  alpha21 = 1,
  mu = 1e-6,
  r1 = 2.611061e-01,
  r2 = 3.397670e-01,
  ratios = median_df$median_VAF[2:4],
  ratios_sd = median_df$sd[2:4]
)

if(F){
  fit <- stan(file = 'F3E_intra_competition.stan', data = stan_data, 
              iter = 3000, chains = 1, 
              warmup = 2000, verbose = TRUE)
  pairs(fit)
}

fit = fit2
posterior = rstan::extract(fit)
param_sets1 <- expand.grid(r1 = stan_data$r1, 
                           r2 = stan_data$r2, 
                           K = posterior$K, 
                           alpha12 = 1, alpha21 = 1, 
                           mu = stan_data$mu)
param_sets1$alpha12 = 1-posterior$deltaa
combined_results = withinhost_fun(param_sets1, 
                                  state = c(stan_data$V10, stan_data$V20))
colnames(combined_results) = c('time', 'V1', 'V2')
total_v = combined_results$V1 + combined_results$V2
combined_results$V1 = combined_results$V1/total_v 
combined_results$V2 = combined_results$V2/total_v
dfribbon = combined_results %>% group_by(time) %>% 
  summarise(meanV1 = mean(V1),
            minV1 = quantile(V1, probs = 0.025), 
            maxV1 = quantile(V1, probs = 0.975),
            meanV2 = mean(V2),
            minV2 = quantile(V2, probs = 0.025), 
            maxV2 = quantile(V2, probs = 0.975)) %>%
  as.data.frame()
dfribbon = rbind(dfribbon[,c(1:4)], setNames(dfribbon[,c(1,5:7)], names(dfribbon[,c(1:4)])))
dfribbon$group = c(rep('V1', 101),rep('V2', 101))

values = c(hue_pal()(3)[1], hue_pal()(3)[3])

p2 = ggplot() + 
  geom_point(data = median_df[1:4,], 
             aes(x = as.numeric(as.character(Time)), 
                 y = median_VAF), 
             size = 2, color = values[1]) +
  geom_errorbar(data = median_df[1:4,], 
                aes(x = as.numeric(as.character(Time)), 
                    ymin = median_VAF - sd, ymax = median_VAF + sd), 
                width = 2, color = values[1]) +
  geom_line(data = dfribbon, 
            aes(x = time, y = meanV1,
                group = group, color = group)) +
  geom_ribbon(data = dfribbon, 
              aes(x = time, ymin = minV1, ymax = maxV1,
                  group = group, fill = group),
              alpha = 0.3) +
  theme_bw() + 
  scale_color_manual(values = values) +
  scale_fill_manual(values = values) +
  xlab('Time unit') + 
  ylab('') +
  theme(legend.position = "none") +
  scale_x_continuous(minor_breaks = seq(0 , 72, 24),
                     breaks = c(0,24,48,72)) +
  scale_y_continuous(limits = c(0, 1),
                     minor_breaks = seq(0 , 1, 0.25),
                     breaks = c(0,0.5,1),
                     labels = c('0.0','0.5','1.0')) +
  coord_cartesian(xlim = c(0,72))
# fit2 = fit
p2
if(F){
  save(fit1, fit2, file = 'evoSIR/intra_competition.rdata')
}
############## plot ##########

pdf(file = 'Output/F3E_intra_competition.pdf', width = 2, height = 1.4)
print(p1)
print(p2)
dev.off()

posterior1 = rstan::extract(fit1)
hist(posterior1$K[posterior1$K<1000])
hist(posterior1$deltaa)
posterior2 = rstan::extract(fit2)
hist(posterior2$K[posterior1$K<1000])
hist(posterior2$deltaa)

library(hrbrthemes)
data = data.frame(K1 = posterior1$K, a1 = 1-posterior1$deltaa,
           K2 = posterior2$K, a2 = 1-posterior2$deltaa)

p3 = ggplot(data) +
  geom_density(aes(x = a1, y = ..density..), fill="#69b3a2") +
  geom_density( aes(x = a2, y = -..density..), fill= "#404080") +
  theme_minimal() +
  xlab(expression(alpha[AB])) + ylab('Density') +
  annotate("segment", x = mean(data$a1), xend = mean(data$a1),
           y = 0, yend = Inf, color = "black", linetype = "dashed") +
  annotate("segment", x = mean(data$a2), xend = mean(data$a2),
           y = 0, yend = -Inf, color = "black", linetype = "dashed") +  
  geom_hline(yintercept = 0, color = "black", linetype = "solid")  +
  scale_y_continuous(breaks = c(-10,-5,0,2.5)) +
  theme(panel.grid.minor = element_blank())

# > mean(data$a1)
# [1] 0.4882342
# > mean(data$a2)
# [1] 0.1670126

p4 = ggplot(data) +
   geom_density(aes(x = log10(K1), y = ..density..), fill="#69b3a2" ) +
   geom_density( aes(x = log10(K2), y = -..density..), fill= "#404080") +
   theme_minimal() +
   annotate("segment", x = log10(mean(data$K1)), xend = log10(mean(data$K1)),
            y = 0, yend = Inf, color = "black", linetype = "dashed") +
   annotate("segment", x = log10(mean(data$K2)), xend = log10(mean(data$K2)),
            y = 0, yend = -Inf, color = "black", linetype = "dashed") +  
   geom_hline(yintercept = 0, color = "black", linetype = "solid") +
   scale_x_continuous(breaks = c(2,3,4), 
                      labels = c(expression(10^2), expression(10^3), expression(10^4))) +
   xlab("K") + ylab('Density') +
  theme(panel.grid.minor = element_blank())

p3
p4
# > log10(mean(data$K1))
# [1] 2.838322
# > log10(mean(data$K2))
# [1] 2.069668

# > mean(data$K1)
# [1] 689.1628
# > mean(data$K2)
# [1] 117.4001

pdf(file = 'Output/S7B_intra_competition_par.pdf', width = 2, height = 1.5)
print(p3)
print(p4)
dev.off()


################### transmission ###########
## Calu-3
df = myfun('Calu')
# df =df[df$Time %in% unique(df$Time)[1:4],]
df = df %>%
  group_by(Time, Strain, Replicate) %>%
  summarize(VAF = mean(VAF)) 
median_df <- df %>%
  group_by(Time) %>%
  summarize(median_VAF = median(VAF),
            sd = sd(VAF)) %>%
  as.data.frame()

stan_data <- list(
  V10 = median_df[1,2],                   
  V20 = 1-median_df[1,2],  
  T = 72+1,
  alpha21 = 1,
  mu = 1e-6,
  r1 = 1.911578e-01,
  r2 = 2.272174e-01,
  ratios = median_df$median_VAF[2:4],
  ratios_sd = median_df$sd[2:4]
)

fit = fit1
posterior = rstan::extract(fit)
param_sets1 <- expand.grid(r1 = stan_data$r1, 
                           r2 = stan_data$r2, 
                           K = posterior$K, 
                           alpha12 = 1, alpha21 = 1, 
                           mu = stan_data$mu)
param_sets1$alpha12 = 1-posterior$deltaa
combined_results = withinhost_fun(param_sets1, 
                                  state = c(stan_data$V10, stan_data$V20))
colnames(combined_results) = c('time', 'V1', 'V2')
total_v = combined_results$V1 + combined_results$V2
combined_results$V1 = combined_results$V1/total_v 
combined_results$V2 = combined_results$V2/total_v

dfcycle = combined_results[combined_results$time == 0,]
dfcycle$cycle = 0
dfcycle = rbind(dfcycle, cbind(combined_results[combined_results$time == 24,], cycle = 1))

for (i in 1:9) {
  print(i)
  v10 = mean(combined_results[combined_results$time == 24,'V1'])
  combined_results = withinhost_fun(param_sets1, state = c(v10, 1-v10))
  colnames(combined_results) = c('time', 'V1', 'V2')
  total_v = combined_results$V1 + combined_results$V2
  combined_results$V1 = combined_results$V1/total_v 
  combined_results$V2 = combined_results$V2/total_v
  dfcycle = rbind(dfcycle, cbind(combined_results[combined_results$time == 24,], cycle = i+1))
}

dfribbon = dfcycle %>% group_by(cycle) %>% 
  summarise(meanV1 = mean(V1),
            minV1 = quantile(V1, probs = 0.025), 
            maxV1 = quantile(V1, probs = 0.975),
            meanV2 = mean(V2),
            minV2 = quantile(V2, probs = 0.025), 
            maxV2 = quantile(V2, probs = 0.975)) %>%
  as.data.frame()
dfribbon = rbind(dfribbon[,c(1:4)], setNames(dfribbon[,c(1,5:7)], names(dfribbon[,c(1:4)])))
dfribbon$group = c(rep('V1', nrow(dfribbon)/2),
                   rep('V2', nrow(dfribbon)/2))

values = c(hue_pal()(3)[1], hue_pal()(3)[3])

p5 = ggplot() + 
  geom_point(data = dfribbon, 
             aes(x = cycle, y = meanV1, group = group, color = group), 
             size = 2) +
  geom_errorbar(data = dfribbon, 
                aes(x = cycle, ymin = minV1, ymax = maxV1, group = group, color = group), 
                width = 0.5) +
  geom_line(data = dfribbon, 
            aes(x = cycle, y = meanV1,
                group = group, color = group)) +
  theme_bw() + 
  scale_color_manual(values = values) +
  scale_fill_manual(values = values) +
  xlab('Transmission cycle') + 
  ylab('') +
  theme(legend.position = "none") +
  scale_x_continuous(breaks = c(0,1,2,5,10),
                     minor_breaks = c(0,1,2,5,10)) +
  scale_y_continuous(breaks = c(0,0.5,1), 
                     labels = c('0.0','0.5','1.0'))

p5

param_sets1 <- expand.grid(r1 = stan_data$r1, 
                           r2 = stan_data$r2, 
                           K = posterior$K, 
                           alpha12 = 1, alpha21 = 1, 
                           mu = stan_data$mu)
param_sets1$alpha12 = 1-posterior$deltaa
combined_results = withinhost_fun(param_sets1, 
                                  state = c(stan_data$V10, stan_data$V20))
colnames(combined_results) = c('time', 'V1', 'V2')
total_v = combined_results$V1 + combined_results$V2
combined_results$V1 = combined_results$V1/total_v 
combined_results$V2 = combined_results$V2/total_v

dfcycle = combined_results[combined_results$time == 0,]
dfcycle$cycle = 0
dfcycle = rbind(dfcycle, cbind(combined_results[combined_results$time == 72,], cycle = 1))

for (i in 1:9) {
  print(i)
  v10 = mean(combined_results[combined_results$time == 72,'V1'])
  combined_results = withinhost_fun(param_sets1, state = c(v10, 1-v10))
  colnames(combined_results) = c('time', 'V1', 'V2')
  total_v = combined_results$V1 + combined_results$V2
  combined_results$V1 = combined_results$V1/total_v 
  combined_results$V2 = combined_results$V2/total_v
  dfcycle = rbind(dfcycle, cbind(combined_results[combined_results$time == 72,], cycle = i+1))
}

dfribbon = dfcycle %>% group_by(cycle) %>% 
  summarise(meanV1 = mean(V1),
            minV1 = quantile(V1, probs = 0.025), 
            maxV1 = quantile(V1, probs = 0.975),
            meanV2 = mean(V2),
            minV2 = quantile(V2, probs = 0.025), 
            maxV2 = quantile(V2, probs = 0.975)) %>%
  as.data.frame()
dfribbon = rbind(dfribbon[,c(1:4)], setNames(dfribbon[,c(1,5:7)], names(dfribbon[,c(1:4)])))
dfribbon$group = c(rep('V1', nrow(dfribbon)/2),
                   rep('V2', nrow(dfribbon)/2))

values = c(hue_pal()(3)[1], hue_pal()(3)[3])

p6 = ggplot() + 
  geom_point(data = dfribbon, 
             aes(x = cycle, y = meanV1, group = group, color = group), 
             size = 2) +
  geom_errorbar(data = dfribbon, 
                aes(x = cycle, ymin = minV1, ymax = maxV1, group = group, color = group), 
                width = 0.5) +
  geom_line(data = dfribbon, 
            aes(x = cycle, y = meanV1,
                group = group, color = group)) +
  theme_bw() + 
  scale_color_manual(values = values) +
  scale_fill_manual(values = values) +
  xlab('Transmission cycle') + 
  ylab('') +
  theme(legend.position = "none") +
  scale_x_continuous(breaks = c(0,1,2,5,10),
                     minor_breaks = c(0,1,2,5,10)) +
  scale_y_continuous(breaks = c(0,0.5,1), limits = c(0,1),
                     labels = c('0.0','0.5','1.0'))
p6

pdf(file = 'Output/F3F_intra_competition_cycle.pdf', width = 2, height = 1.4)
print(p5)
print(p6)
dev.off()
