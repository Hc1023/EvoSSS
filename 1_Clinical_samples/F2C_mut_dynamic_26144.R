rm(list = ls())
library(treeio)
library(ggtree)
library(ggplot2)
library(ape)
library(scales)
library(RColorBrewer)
library(ggnewscale)
library(tidyverse)
library(ggpubr)
library(scatterpie)

df1 = read.csv('F2C_sample_snpdate.csv')
df2 = read.csv('F1A_clinical_info.csv')

# P30
# P34
s = df2[df2$Names_ID == 'P30','Sample']
s = df2[df2$Names_ID == 'P34','Sample']
# G26144T
df = data.frame(matrix(0, nrow = length(s), ncol = 3))
colnames(df) = c('Date','b1','b2')
for (i in 1:length(s)) {
  tmp = df1[df1$Mut == 'G26144T',
            c(grep(paste0('S',s[i],'_'), colnames(df1)))] %>% unlist()
  df$Date[i] = df2[df2$Sample == s[i],'Date']
  df$b1[i] = tmp[1] - tmp[2] - tmp[3]
  df$b2[i] = tmp[2]
  
}

df$x = seq(0,5,5/(nrow(df)+1))[1:nrow(df)+1]
df$rad = log10(df$b1 + df$b2)/10
df$p = df$b2/(df$b1 + df$b2)
df$y = 1
t = 1.5
df$y2 = df$p + t

d = as.Date(df$Date)
d = c(format(d[1], "%b %d"), paste0('+',d[2:length(d)]-d[1]))

dftext = df[c('x','p','y2')]
v <- ifelse(dftext$p > 0.5, 2, -1)

p1 = ggplot(data = df) + 
  geom_line(aes(x=x,y=y2), color = '#5fbfb1') +
  geom_point(aes(x=x,y=y2),  color = '#5fbfb1') +
  new_scale_color() +
  geom_scatterpie(data = df,
                  aes(x=x, y=y, group= x, r = rad), 
                  cols=c('b2','b1')) + 
  coord_equal() +
  scale_fill_manual(values = alpha(c('#5fbfb1', 'grey'),0.8),
                    guide = 'none') +
  labs(x = '', y = '') + 
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  scale_x_continuous(breaks = df$x, labels =  d,
                     limits = c(0.51, 4+0.55)) +
  scale_y_continuous(breaks = c(1,t,t+0.5,t+1), 
                     labels = c('26144',0.0,0.5,1.0),
                     limits = c(0.6, t+1.1)) +
  geom_text(data=data.frame(dftext),
            aes(label = sprintf("%.2f", p),
                x=x, y=y2), size = 3,
            inherit.aes = F, vjust = v) 
p1
p2

pdf(file = 'Output/F2C_mut_dynamic_26144.pdf', width = 2.3, height = 2.2)
print(p1)
print(p2)
dev.off()




