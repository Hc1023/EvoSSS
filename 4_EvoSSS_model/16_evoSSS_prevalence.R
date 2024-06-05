rm(list = ls())
library(dplyr)
library(ggplot2)
library(scales)

df = read.csv('../3_Epidemiological_analysis/Covid19CasesGISAID.csv')
df$Var1 = as.Date(df$Var1)
dfplot = df
x = format(dfplot$Var1,'%y %b')
dfplot$mon = as.Date(paste(x,01), '%y %b %d')
for (i in 1:nrow(dfplot)) {
  dfplot$mon[i] = seq(dfplot$mon[i], length=2, by="1 month")[2]-1
}
dfsum = dfplot %>% 
  group_by(mon, Mutations) %>% 
  summarise(Freq = sum(Freq))

dfsum2 = dfplot %>% 
  group_by(mon) %>% 
  summarise(Freq = sum(Freq))
dfsum2 = merge(dfsum, dfsum2, by = "mon", all = T)


for (i in 1:nrow(dfsum2)) {
  t = binom.test(dfsum2$Freq.x[i],dfsum2$Freq.y[i], p = 1/3)
  dfsum2$ratio[i] = t$estimate
  dfsum2$lower[i] = t$conf.int[1]
  dfsum2$upper[i] = t$conf.int[2]
}
dfsum2 = dfsum2[-c((nrow(dfsum2)-3):nrow(dfsum2)),]

values = c(hue_pal()(3)[1], hue_pal()(3)[3], hue_pal()(3)[2], hue_pal()(4)[4])

pd <- position_dodge(5) 

# 0 0.005 0.010 0.015 0.020 0.2 0.4 0.6 0.8 1.0
dfsum2$Mutations = factor(dfsum2$Mutations, 
                          levels = unique(dfsum2$Mutations)[c(3,4,2,1)])


load('h1.rdata')
# simu_Onset
simu1 = as.matrix(simu_Onset[plot_data$V == 'A',])
rownames(simu1) = 1:nrow(simu1)
simu2 = as.matrix(simu_Onset[plot_data$V == 'B',])
rownames(simu2) = 1:nrow(simu1)
# str(simu1)

df2 = plot_data[plot_data$V == 'A',]
simup = simu1 / (simu1 + simu2)
ci_lower <- apply(simup, 1, quantile, probs = 0.025, na.rm = T)
ci_upper <- apply(simup, 1, quantile, probs = 0.975, na.rm = T)
df2$Fitted =  rowMeans(simup)
df2$LowerCI =  ci_lower
df2$UpperCI =  ci_upper
df3 = plot_data[plot_data$V == 'B',]
simup = simu2 / (simu1 + simu2)
ci_lower <- apply(simup, 1, quantile, probs = 0.025, na.rm = T)
ci_upper <- apply(simup, 1, quantile, probs = 0.975, na.rm = T)
df3$Fitted =  rowMeans(simup)
df3$LowerCI =  ci_lower
df3$UpperCI =  ci_upper
df2 = rbind(df2,df3)

dfsum2 = dfsum2[dfsum2$Mutations %in% c('Lineage A', 'Lineage B'),]
dfsum2$V  = 'A'
dfsum2$V[dfsum2$Mutations == 'Lineage B'] = 'B'
p = ggplot() + 
  geom_point(data = dfsum2, 
             aes(x = mon, y = ratio, colour = V),
             position = pd) + 
  geom_line(data = dfsum2, 
            aes(x = mon, y = ratio, colour = V),
            position = pd) +
  geom_errorbar(data = dfsum2, 
                aes(x = mon, y = ratio, colour = V,
                    ymin = lower, ymax = upper), 
                width=2, position = pd) +
  scale_x_date(breaks = seq(as.Date('2020-01-01'), as.Date('2022-11-01'), by="6 months"),
               minor_breaks = seq(as.Date('2019-12-01'), as.Date('2022-11-01'), by ='1 month'),
               date_labels = "%y-%b") +
  geom_ribbon(data = df2, 
              aes(x = x, group = V, 
                  ymin = LowerCI, ymax = UpperCI, 
                  fill = V)) +  # Confidence interval
  geom_line(data = df2, 
            aes(x = x, y = Fitted, 
                group = V, color = V), linewidth = 1) +
  xlab('') + ylab('p (%)') + theme_bw() +
  scale_color_manual(values = alpha(values, 0.6)) +
  scale_fill_manual(values = alpha(values, 0.6)) +
  scale_y_continuous(trans='log10',
                     breaks = c(1e-4,1e-2,1),
                     labels = c(expression(10^-2),expression(10^0),
                                expression(10^2))) +
  coord_cartesian(xlim = c(as.Date('2019-12-31'), as.Date('2021-10-31')),
                  ylim = c(1e-5,1)) +
  theme(legend.position = 'none',
        plot.margin = margin(0.5,1,0,0.5, "cm") # trbl
  ) 
p
pdf("Output/evoSSS_prevalence.pdf", width = 2.8, height = 1.3)
print(p)
dev.off()


