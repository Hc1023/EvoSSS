rm(list = ls())
library(dplyr)
library(ggplot2)
library(scales)

df = read.csv('Covid19CasesGISAID.csv')
df$Var1 = as.Date(df$Var1)
dfplot = df#[df$Var1 < as.Date('2021-12-01'),]
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
p = ggplot(data = dfsum2, 
           aes(x = mon, y = ratio, colour = Mutations)) + 
  geom_point(position = pd) + geom_line(position = pd) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width=2, position = pd) +
  scale_x_date(breaks = seq(as.Date('2020-01-01'), as.Date('2022-11-01'), by="6 months"),
               minor_breaks = seq(as.Date('2019-12-01'), as.Date('2022-11-01'), by ='1 month'),
               date_labels = "%y-%b") +
  xlab('') + ylab('Prevalence (%)') + theme_bw() +
  scale_color_manual(values = alpha(values, 0.6)) +
  scale_fill_manual(values = alpha(values, 0.1)) +
  scale_y_continuous(trans='log10',
                     breaks = c(1e-6,1e-5,1e-4,1e-3,1e-2,1e-1,1),
                     labels = c('1e-4','1e-3','1e-2','1e-1','1','10','100')) +
  annotation_logticks(linewidth = 0.1) +
  coord_cartesian(xlim = c(as.Date('2019-12-31'), as.Date('2022-09-01'))) +
  theme(legend.position = 'none',
        plot.margin = margin(0.5,1,0,0.5, "cm") # trbl
  ) 
p
pdf(paste0("Output/prevalence", ".pdf"), width = 4.4, height = 2)
print(p)
dev.off()


