rm(list = ls())
library(data.table)
library(tidyverse)
library(readxl)
library(ggpubr)
library(scales)
library(rstatix)

df = read.csv("calu3_ct.csv")

get_df_long = function(df1){
  colnames(df1) = c('Strain', '0', '8', '24', '48', '72')
  
  df2 = df1[-c(1:10), ] %>%
    group_by(Strain) %>%
    summarise_at(vars(-group_cols()), mean, na.rm = TRUE) %>%
    as.data.frame()
  
  df2 = rbind(df1[1:10, ], df2)
  df2$Strain <- factor(c(rep('A', 6), rep('B', 4), rep('A+B', 6)),
                       levels = c('A', 'B', 'A+B'))
  
  df_long <- gather(df2, key = "Time",
                    value = "Ct",-Strain)
  df_long$Time = factor(df_long$Time, levels = unique(df_long$Time))
  return(df_long)
}


get_signif = function(t) {
  subset_data <- df_long[df_long$Time == t, ]
  anova_result <- aov(Ct ~ Strain, data = subset_data)
  p_value = summary(anova_result)[[1]]$`Pr(>F)`[1]
  signif_code = ifelse(p_value < 0.001, "***",
                       ifelse(p_value < 0.01, "**",
                              ifelse(
                                p_value < 0.05, "*",
                                ifelse(p_value < 0.1, ".", " ")
                              )))
  return(paste0('p=',paste(format(round(p_value, 3), nsmall = 3), 
                           signif_code)))
}
gene = c('ORF1ab', 'N')
# pdf(paste0("Output/replication.pdf"), width = 2.7, height = 2.3)

for (i in 1:1) {
  if (i==1) {
    df1 = df[, 1:6]
  }else{
    df1 = df[, c(1, 9:13)]
  }
  df_long = get_df_long(df1)
  x= sapply(c('0', '8', '24', '48', '72'), get_signif)
  y = df_long %>% group_by(Time) %>%
    summarise_at(vars(Ct), function(c){
      ifelse(min(c) > 25, min(c) - 2, max(c) + 0.2)
    })
  
  annotation_df <- data.frame(Time = names(x), Annotation = x, 
                              y = y$Ct)
  values = c(hue_pal()(3)[1], hue_pal()(3)[3], '#aa85a6')
  p = ggplot(df_long, aes(x = Time, y = Ct, fill = Strain)) +
    geom_boxplot(aes(color = Strain)) +
    geom_point(
      position = position_jitterdodge(),
      show.legend = FALSE,
      aes(color = Strain),
      alpha = 0.5,
      size = 1
    ) +
    scale_color_manual(values = values) +
    scale_fill_manual(values = alpha(values, 0.3)) +
    labs(x = "Time points (h)",
         y = paste(gene[i], "(Ct)")) +
    theme_bw() +
    theme(
      legend.position = c(0.85, 0.7),
      legend.background = element_rect(color = NA, fill = NA),
      legend.key = element_blank()
    ) +
    annotate("text", x = annotation_df$Time, y = annotation_df$y,
             label = x, vjust = -0.5, size = 2.8) +
    coord_cartesian(ylim = c(12, 33.9))
  
  print(p)
  
}
# dev.off()

library(ggplot2)   
library(dplyr)
library(tidyr)
library(scales)
library(deSolve)   
library(rstan)

df = read.csv("calu3_ct.csv")

timepoint = c(0,8,24,48,72)
dat = data.frame()
for (i in 1:length(timepoint)) {
  d = df[,c(1,i+1)]
  colnames(d)[2] = 'Ct' 
  dat = rbind(dat,d)
}

dat$timepoint = rep(timepoint, each = nrow(df))
dat$group = rep(c(rep('A',6),rep('B',4),rep('A+B',18)), length(timepoint))

dat$V = 2^(-dat$Ct + max(dat$Ct))
dat$Vlog = log(dat$V)

rvec = c()
g = c('A','B','A+B')
for (i in 1:3) {
  print(g[i])
  dat1 = dat[dat$timepoint %in% timepoint[1:3] & dat$group == g[i],]
  fit1 = lm(formula = Vlog ~ timepoint, data = dat1)
  result1 = summary(fit1)
  print(result1$coefficients[2,])
  rvec[i] = result1$coefficients[2,1]
  
}

# [1] "A"
# Estimate   Std. Error      t value     Pr(>|t|) 
# 3.452843e-01 1.874275e-02 1.842228e+01 3.383621e-12 
# [1] "B"
# Estimate   Std. Error      t value     Pr(>|t|) 
# 3.994415e-01 2.058320e-02 1.940619e+01 2.880877e-09 
# [1] "A+B"
# Estimate   Std. Error      t value     Pr(>|t|) 
# 3.580951e-01 1.165525e-02 3.072394e+01 5.186158e-35 
viral_model <- function(t, state, parameters) {
  V1 <- state[1]
  r1 <- parameters["r1"]
  K <- parameters["K"]
  dV1dt <- r1 * V1 * (1 - V1 / K)
  return(list(dV1dt))
}
out_df = data.frame()
for (i in 1:3) {
  state <- c(V1 = 1)
  times <- seq(0, 100, by = 1)
  params = c(r1 = rvec[i], K = 2*median(dat[dat$timepoint == 48,'V']))
  out <- ode(y = state, times = times, func = viral_model, parms = params)
  out <- as.data.frame(out)
  out$group = g[i]
  out_df = rbind(out_df, out)
}

out_df$Ct = -log2(out_df$V1) + max(dat$Ct)

# Plotting
pd = 3
values = c(hue_pal()(3)[1], hue_pal()(3)[3],  '#aa85a6')
out_df$group = factor(out_df$group, levels = g)
dat$group = factor(dat$group, levels = c('B','A+B','A'))
annotation_df$Time = c(0.6,8,26,48,72) + 21
annotation_df$y = c(31,28.6,20,13,11.6)
p1 = ggplot(data = out_df, 
            aes(x = time, y = Ct, group = group,
                color = group, fill = group)) +
  geom_line() +
  geom_boxplot(data = dat[dat$group == 'A',], 
               aes(x = timepoint+pd, y = Ct, 
                   group = timepoint),
               color = values[1], fill = alpha(values[1], 0.3),
               outlier.shape = NA, width = 1.9*pd) +
  geom_boxplot(data = dat[dat$group == 'B',], 
               aes(x = timepoint-pd, y = Ct, 
                   group = timepoint), 
               color = values[2], fill = alpha(values[2], 0.3),
               outlier.shape = NA, width = 1.9*pd) +
  geom_boxplot(data = dat[dat$group == 'A+B',], 
               aes(x = timepoint, y = Ct, 
                   group = timepoint), 
               color = values[3], fill = alpha(values[3], 0.3),
               outlier.shape = NA, width = 1.9*pd) +
  geom_jitter(data = dat, 
              aes(x = timepoint, y = Ct, color = group, group = group),
              alpha = 0.5,
              position = position_jitterdodge(jitter.width = 1,
                                              dodge.width = 4*pd),
              show.legend = F) + 
  scale_color_manual(name = '', values = values) +
  scale_y_continuous(trans = 'log10') +
  scale_x_continuous(breaks = c(2,8,24,48,72),
                     minor_breaks = c()) +
  theme_bw() +
  labs(x = 'Time points (h)', y = 'Orf1ab (Ct)') +
  theme(legend.position = c(0.8,0.8),
        legend.key.size = unit(0.4,'cm'),
        legend.background = element_blank(),
        legend.key = element_blank()) +
  annotate("text", x = annotation_df$Time, y = annotation_df$y,
           label = x, vjust = -0.5, size = 2.8)
p1

pdf(paste0("Output/replication_fit.pdf"), width = 2.5, height = 2)
print(p1)
dev.off()

