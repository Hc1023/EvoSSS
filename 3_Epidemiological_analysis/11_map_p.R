rm(list = ls())
library(treeio)
library(ggtree)
library(ggplot2)
library(ape)
library(scales)
library("ggridges")
library(tidyverse)

ddf = read.csv('map.csv')
table(ddf$country)
table(ddf$continent)
colnames(ddf)[2:5] = c('Lineage A', 'Lineage B', '8782T','28144C')
ddf$p = ddf$`Lineage A`/(ddf$`Lineage A`+ddf$`Lineage B`)
idx1 = ddf$country=='China'
idx2 = ddf$continent=='Asia' 

region = c(rep('China',table(idx1)[2]),
  rep('Asia',table(idx2)[2]),
  rep('Global',nrow(ddf)))
df = rbind(ddf[idx1,], ddf[idx2,], ddf)
df$region = factor(region, 
                   levels = c('China','Asia','Global'))
df2 = df[df$region != 'China',]
x = df[df$region == 'Asia' & df$country =='China',]
sum(x$`Lineage A`)/(sum(x$`Lineage A`)+sum(x$`Lineage B`))
dftext = df %>% group_by(region) %>% 
  summarise(mean1 = mean(p),
            mean2 = sum(`Lineage B`)/(sum(`Lineage A`)+sum(`Lineage B`))) %>% data.frame()
dfline = data.frame(x = rep(dftext$mean2, each = 2),
                    y = c(1,4.5,2,4.5,3,4.5),
                    group = rep(c('China','Asia','Global'), each = 2))
dfline$region = factor(dfline$group, 
                       levels = c('China','Asia','Global'))

p = ggplot()+
  geom_density_ridges(data=df, aes(x=p, y=region, fill = region),
                      scale = 2, alpha=0.5, bandwidth = 0.1,
                      jittered_points = TRUE, 
                      point_alpha=1, point_shape=21,
                      quantile_lines = T, quantile_fun = mean) +
  theme_ridges() +
  scale_y_discrete(expand = expansion(mult = c(0, 0.55), 
                                      add = c(0.01, 1)))+
  scale_x_continuous(limits = c(-0.02, 1.02)) +
  scale_fill_manual(name = '',
                    values = c('#3172ae','#ea9d5d','#cc706b'),
                    limits = c("Global", "Asia", "China")) +
  scale_color_manual(name = '',
                    values = alpha(c('#3172ae','#ea9d5d','#cc706b'),0.8),
                    limits = c("Global", "Asia", "China")) +
  theme(legend.position = 'none') +
  geom_text(data = data.frame(x = 0.7, y = 4.8, 
                              text = 'By March 1, 2020'), 
            aes(x = x, y = y, label = text), 
            inherit.aes = FALSE) +
  xlab('') + ylab('')
p
pdf(file = 'Output/map_p.pdf', width = 3.6, height = 1.8)
print(p)
dev.off()

df1 = df[df$region == 'Global' & df$p >0.3 & df$country != 'China',]
df1[,c('country','p')]
# country         p
# 22               Spain 0.6428571
# 231           Thailand 0.347826E1
# 241 UnitedArabEmirates 0.3076923
# 25                 USA 0.5165563
df1[df1$country %in% c('Spain','Thailand','USA'),]
# > dftext
# region     mean1     mean2
# 1  China 0.3451470 0.6694796
# 2   Asia 0.2298248 0.7326733
# 3 Global 0.1479226 0.7785016



