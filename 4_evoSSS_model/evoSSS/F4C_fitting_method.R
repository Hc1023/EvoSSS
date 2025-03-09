rm(list = ls())
library(rstan)
library(ggplot2)
library(scales)
library(ggnewscale)
library(tidyverse)
library(dplyr)
library(RColorBrewer)

df = read.csv('../3_Epidemiological_analysis/F1D_Covid19CasesGISAID.csv')

full_dates <- as.Date('2019-12-31') + 1:100
df1 = data.frame(x = c(1:100,1:100),
                 y = sqrt(df[as.Date(df$Var1) %in% full_dates & 
                          df$Mutations %in% c('Lineage A', 'Lineage B'),'Freq']),
                 group = c(rep('A',100), rep('B',100)))

ggplot() +
  geom_point(data = df1, 
             aes(x = x, y = y, group = group, color = group))

lo1 = loess(y ~ x, data=df1[df1$group == 'A',], span=0.3)
lo2 = loess(y ~ x, data=df1[df1$group == 'B',], span=0.3)
df1$smooth= c(predict(lo1),predict(lo2))

idx = seq(1,200,5)
ggplot() +
  geom_point(data = df1[idx,], 
             aes(x = x, y = y, group = group, color = group)) +
  geom_line(data = df1, 
             aes(x = x, y = smooth, group = group, color = group))


df1$model = c(dnorm(1:100, 30, 10)*max(df1$smooth[1:30])/max(dnorm(1:100, 30, 10)),
              dnorm(1:100, 30, 10)*max(df1$smooth[1:30+100])/max(dnorm(1:100, 30, 10)))
dfx = df1[idx,]
dfx[dfx$x>60,'smooth'] = 0

values = c(hue_pal()(3)[1], hue_pal()(3)[3])

p = ggplot() +
  geom_point(data = df1[idx,], 
             aes(x = x, y = y, group = group, color = group),
             alpha = 0.8, size = 0.6) +
  geom_point(data = dfx, 
             aes(x = x, y = smooth, group = group, color = group),
             shape = 4, size = 0.8) +
  # geom_line(data = df1, 
  #           aes(x = x, y = smooth, group = group, color = group), alpha = 0.3) +
  geom_line(data = df1, 
            aes(x = x, y = model, group = group, color = group)) +
  geom_ribbon(data = df1[df1$x>30,], 
              aes(x= x, ymin=model, ymax=smooth, 
                  group=group, fill = group), 
              alpha = 0.2) +
  scale_x_continuous(breaks = c(0,30,60,90), expand = c(0,0)) +
  scale_y_continuous(breaks = c(20,40),expand = c(0,1)) +
  coord_cartesian(xlim = c(-2, 100)) +
  theme_classic() + xlab('') + ylab('') +
  scale_color_manual(values = values) +
  scale_fill_manual(values = values) +
  theme(legend.position = 'none',
        axis.text.y = element_blank(),
        plot.margin = margin(0.2,0.5,0,0, "cm"),
        legend.background = element_blank(),
        plot.background = element_blank(),
        panel.background = element_blank(),
        strip.background = element_blank())
pdf(paste0("Output/F4C_fitting1.pdf"), width = 1.6, height = 1.3)
print(p)
dev.off()

df2 = df1[df1$x>30,]
df2$y2 = sapply(df2$y - df2$model, function(x) max(0,x))
df2$ymax = sapply(df2$smooth - df2$model, function(x) max(0,x))

p = ggplot() +
  geom_point(data = df2[df2$x %in% idx,], 
             aes(x = x, y = y2, 
                 group = group, color = group),
             alpha = 0.8, size = 0.6) +
  geom_ribbon(data = df2, 
              aes(x = x, ymin = 0, ymax = ymax, 
                  group = group, fill = group), 
              alpha = 0.2) +
  scale_x_continuous(breaks = c(0,30,60,90), expand = c(0,0)) +
  scale_y_continuous(breaks = c(20,40), expand = c(0,0)) +
  coord_cartesian(xlim = c(-2, 100), ylim = c(-0.5,47)) +
  theme_classic() + xlab('') + ylab('') +
  scale_color_manual(values = values) +
  scale_fill_manual(values = values) +
  theme(legend.position = 'none',
        axis.text.y = element_blank(),
        plot.margin = margin(0.2,0.5,0,0, "cm"),
        legend.background = element_blank(),
        plot.background = element_blank(),
        panel.background = element_blank(),
        strip.background = element_blank())
pdf(paste0("Output/F4C_fitting2.pdf"), width = 1.4, height = 0.85)
print(p)
dev.off()

p = ggplot() +
  geom_line(data = df1, 
            aes(x = x, y = model, group = group, color = group)) +
  scale_x_continuous(breaks = c(0,30,60,90), expand = c(0,0)) +
  scale_y_continuous(breaks = c(2,4), expand = c(0,1)) +
  coord_cartesian(xlim = c(-2, 100),ylim = c(0.8,5)) +
  theme_classic() + xlab('') + ylab('') +
  scale_color_manual(values = values) +
  theme(legend.position = 'none',
        axis.text.y = element_blank(),
        plot.margin = margin(0.2,0.5,0,0, "cm"),
        legend.background = element_blank(),
        plot.background = element_blank(),
        panel.background = element_blank(),
        strip.background = element_blank(),
        axis.ticks.y=element_blank())
pdf(paste0("Output/F4C_fitting3.pdf"), width = 1.4, height = 0.85)
print(p)
dev.off()

df3 = data.frame(x = c(1:100,1:100)+30,
                 y = sqrt(df[as.Date(df$Var1) %in% (full_dates+30) & 
                               df$Mutations %in% c('Lineage A', 'Lineage B'),'Freq']),
                 group = c(rep('A',100), rep('B',100)))
df3$y[c(1:70,1:70+100)] = df2$y2


ggplot() +
  geom_point(data = df3, 
             aes(x = x, y = y, group = group, color = group))

lo1 = loess(y ~ x, data=df3[df3$group == 'A',], span=0.3)
lo2 = loess(y ~ x, data=df3[df3$group == 'B',], span=0.3)
df3$smooth= c(predict(lo1),predict(lo2))

idx = seq(1,200,5)
ggplot() +
  geom_point(data = df3[idx,], 
             aes(x = x, y = y, group = group, color = group)) +
  geom_line(data = df3, 
            aes(x = x, y = smooth, group = group, color = group))


df3$model = c(dnorm(1:100, 30, 10)*max(df3$smooth[1:30])/max(dnorm(1:100, 30, 10)),
              dnorm(1:100, 30, 10)*max(df3$smooth[1:30+100])/max(dnorm(1:100, 30, 10)))
df3x = df3[idx,]
df3x[df3x$x>60+30,'smooth'] = 0

p = ggplot() +
  geom_point(data = df3[idx,], 
             aes(x = x, y = y, group = group, color = group),
             alpha = 0.8, size = 0.6) +
  geom_point(data = df3x, 
             aes(x = x, y = smooth, group = group, color = group),
             shape = 4, size = 0.8) +
  # geom_line(data = df1, 
  #           aes(x = x, y = smooth, group = group, color = group), alpha = 0.3) +
  geom_line(data = df3, 
            aes(x = x, y = model, group = group, color = group)) +
  geom_ribbon(data = df3[df3$x>30+30,], 
              aes(x= x, ymin=model, ymax=smooth, 
                  group=group, fill = group), 
              alpha = 0.2) +
  scale_x_continuous(breaks = c(0,30,60,90)+30, expand = c(0,0)) +
  scale_y_continuous(breaks = c(20,40), expand = c(0,1)) +
  coord_cartesian(xlim = c(-2, 100)+30) +
  theme_classic() + xlab('') + ylab('') +
  scale_color_manual(values = values) +
  scale_fill_manual(values = values) +
  theme(legend.position = 'none',
        axis.text.y = element_blank(),
        plot.margin = margin(0.2,0.5,0,0, "cm"),
        legend.background = element_blank(),
        plot.background = element_blank(),
        panel.background = element_blank(),
        strip.background = element_blank())
pdf(paste0("Output/F4C_fitting4.pdf"), width = 1.6, height = 1.3)
print(p)
dev.off()


df2$y2 = sapply(df2$y - df2$model, function(x) max(0,x))
df2$ymax = sapply(df2$smooth - df2$model, function(x) max(0,x))

df4 = data.frame(x = c(1:100,1:100)+30*2,
                 y = sqrt(df[as.Date(df$Var1) %in% (full_dates+30*2) & 
                               df$Mutations %in% c('Lineage A', 'Lineage B'),'Freq']),
                 group = c(rep('A',100), rep('B',100)))
df4$y[c(1:70,1:70+100)] =  sapply(df3$y - df3$model, function(x) max(0,x))[-c(1:30,1:30+100)]
lo1 = loess(y ~ x, data=df4[df4$group == 'A',], span=0.3)
lo2 = loess(y ~ x, data=df4[df4$group == 'B',], span=0.3)
df4$smooth= c(predict(lo1),predict(lo2))
df4$model = c(dnorm(1:100, 30, 10)*max(df4$smooth[1:30])/max(dnorm(1:100, 30, 10)),
              dnorm(1:100, 30, 10)*max(df4$smooth[1:30+100])/max(dnorm(1:100, 30, 10)))
df5 = data.frame(x = c(1:100,1:100)+30*3,
                 y = sqrt(df[as.Date(df$Var1) %in% (full_dates+30*3) & 
                               df$Mutations %in% c('Lineage A', 'Lineage B'),'Freq']),
                 group = c(rep('A',100), rep('B',100)))
df5$y[c(1:70,1:70+100)] =  sapply(df4$y - df4$model, function(x) max(0,x))[-c(1:30,1:30+100)]
lo1 = loess(y ~ x, data=df5[df5$group == 'A',], span=0.3)
lo2 = loess(y ~ x, data=df5[df5$group == 'B',], span=0.3)
df5$smooth= c(predict(lo1),predict(lo2))
df5$model = c(dnorm(1:100, 30, 10)*max(df5$smooth[1:30])/max(dnorm(1:100, 30, 10)),
              dnorm(1:100, 30, 10)*max(df5$smooth[1:30+100])/max(dnorm(1:100, 30, 10)))

dfmodel = df1
dfmodel$model2[c((30+1):100,(30+1):100+100)] = df3$model[df3$x %in% c((30+1):100,(30+1):100+100)]
dfmodel$model3[c((30*2+1):100,(30*2+1):100+100)] = df4$model[df4$x %in% c((30*2+1):100,(30*2+1):100+100)]
dfmodel$model4[c((30*3+1):100,(30*3+1):100+100)] = df5$model[df5$x %in% c((30*3+1):100,(30*3+1):100+100)]



p = ggplot() +
  geom_point(data = dfmodel[idx,], 
             aes(x = x, y = y, group = group, color = group),
             alpha = 0.8, size = 0.6) +
  geom_line(data = dfmodel, 
            aes(x = x, y = model, group = group, color = group)) +
  geom_line(data = dfmodel, 
            aes(x = x, y = model2, group = group, color = group)) +
  geom_line(data = dfmodel, 
            aes(x = x, y = model3, group = group, color = group)) +
  geom_line(data = dfmodel, 
            aes(x = x, y = model4, group = group, color = group)) +
  scale_x_continuous(breaks = c(0,30,60,90), expand = c(0,0)) +
  scale_y_continuous(breaks = c(20,40), 
                     labels = c(2000,4000), 
                     expand = c(0,1)) +
  coord_cartesian(xlim = c(-2, 100)) +
  theme_classic() + xlab('') + ylab('') +
  scale_color_manual(values = values) +
  scale_fill_manual(values = values) +
  theme(legend.position = 'none',
        plot.margin = margin(0.2,0.5,0,0, "cm"),
        axis.text.y = element_blank(),
        legend.background = element_blank(),
        plot.background = element_blank(),
        panel.background = element_blank(),
        strip.background = element_blank())

pdf(paste0("Output/F4C_fitting5.pdf"), width = 1.6, height = 1.3)
print(p)
dev.off()