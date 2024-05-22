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
load('Onsets_mat_list.rdata')

dfposterior = data.frame()
for (i in 1:length(fitlist)) {
  fit = fitlist[[i]]
  posterior_samples = rstan::extract(fit)
  dfposterior[i,1] = i
  dfposterior$m[i] = mean(posterior_samples$contact)
  dfposterior$q1[i] = quantile(posterior_samples$contact, 0.025)
  dfposterior$q2[i] = quantile(posterior_samples$contact, 0.975)
}

# n = i-1
# Onsets_mat_list[[n]]
# 1:nrow(Onsets_mat)+poolday*n
poolday = 30

for (i in 2:length(Onsets_mat_list)) {
  n = i-1
  Onsets_mat = Onsets_mat_list[[i]]
  y = Onsets_mat[,1] + Onsets_mat[,2]
  x = which.max(y) + poolday*n
  y = max(y)/28
  dfposterior[n,5:6] = c(x,y)
}

dfposterior$size = log2(dfposterior$V6)
dfposterior$size = dfposterior$size - min(dfposterior$size) + 0.3
dfposterior$x = dfposterior$V5 + as.Date('2019-12-31')
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

dfplot_simu$y = dfplot_simu$y/28
df2$y = df2$y/28
dfplot_simu$Date = dfplot_simu$x + as.Date('2019-12-31')
df2$Date = df2$x + as.Date('2019-12-31')
str(df2)
values = c(hue_pal()(3)[1], hue_pal()(3)[3])
colors <- colorRampPalette(c("white", "darkblue"))
color_palette <- colors(5)
breaks_colors <- color_palette
size_labels = c(1,2,4,8,16)*1000
size_breaks = log2(size_labels) - min(log2(dfposterior$V6)) + 0.3

df = read.csv('../3_Epidemiological_analysis/Covid19CasesGISAID.csv')

df$Var1 = as.Date(df$Var1)
# df = df[df$Var1 < as.Date('2021-12-01'),]
dfobserve = df[df$Mutations %in% c('Lineage A', 'Lineage B'),]
dfobserve$group = 'A'
dfobserve$group[dfobserve$Mutations == 'Lineage B'] = 'B'
p = ggplot(dfposterior, aes(x, m)) +
  geom_line(color = 'darkblue', alpha = 0.8) +
  geom_point(aes(size = size, fill = size), 
             shape = 21, alpha = 0.8) + 
  scale_fill_gradient2(low = "white", high = "darkblue",
                       breaks = c(1,3,5), guide = 'none') +
  scale_size_continuous(name = 'Peak', 
                        breaks = size_breaks,
                        label = size_labels) + 
  guides(size = guide_legend(override.aes = list(fill = breaks_colors))) +
  new_scale_color() +
  geom_point(data = dfobserve, 
             aes(x = Var1, y = Freq, 
                 group = group, color = group),
             size = 0.4, shape = 16) +
  geom_line(data = dfplot_simu, 
            aes(x = Date, y = y, 
                group = group, color = color), 
            alpha = 0.2) + 
  geom_line(data = df2[-c(1:2),], aes(x = Date, y = y, 
                            group = color, color = color)) + 
  
  scale_color_manual(name="Variant",
                     labels=c("A", "B"),
                     values = alpha(values, 0.6)) +
  theme_bw() +
  scale_y_continuous(trans='log10', 
                     breaks = c(1,10,100,1000,10000),
                     sec.axis = sec_axis(
                       transform = ~ ., 
                       name = "Spreading capacity", # Susceptible/Seed
                       breaks = c(1,10,100,1000,10000)
                     )) +
  theme(legend.background = element_blank(),
        legend.key = element_blank(),
        legend.key.size = unit(0.4, units = 'cm'),
        legend.spacing = unit(0.01, units = 'cm'),
        axis.text.y.right = element_text(color = alpha('darkblue', 0.8)),
        axis.title.y.right = element_text(color = alpha('darkblue', 0.8))) +
  annotation_logticks(sides = "l", linewidth = 0.1, alpha = 0.5) +
  xlab('Cycle') + ylab('Susceptible ratio') +
  scale_x_date(breaks = seq(as.Date('2020-01-01'), as.Date('2022-11-01'), by="6 months"),
               minor_breaks = seq(as.Date('2019-12-01'), as.Date('2022-11-01'), by ='1 month'),
               date_labels = "%y-%b") +
  xlab('') + ylab('Cases') + 
  coord_cartesian(xlim = c(as.Date('2019-12-31'), as.Date('2021-10-31')),
                  ylim = c(1,2*10^4)) +
  theme(legend.position = 'right')

p
pdf(paste0("Output/evoSSS_spreading_capacity.pdf"), width = 6, height = 2.5)
print(p)
dev.off()
