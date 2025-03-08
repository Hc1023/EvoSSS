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
df3 = read.csv('F2A_clinical_coverage.csv')
df1 = df3[1:149,]
df2 = df3[149+1:149,]
df2 <- df2[order(match(df2$Sample, df1$Sample)), ]
df3 = rbind(df1,df2)

plotfun = function(draw = F){
  ddf = df3[df3$clade %in% c('28144C','8782T'),]
  ddf$region = 1:nrow(ddf)
  
  ddf$x = c(1:10,1:10)
  ddf$y = rep(c(1,2), each = 10)
  ddf$rad = log10(ddf$total_reads)/10
  ddf$rad[ddf$rad == -Inf] = 0.1
  values = c(hue_pal()(3)[1], hue_pal()(3)[3])
  
  if(draw){
    library(scatterpie)
    p = ggplot() + 
      geom_scatterpie(data = ddf,
                      aes(x=x, y=y, group=region, r = rad), 
                      cols=c('b2','b1')) + 
      coord_equal() +
      scale_fill_manual(name="SNP",
                        values = alpha(values,0.8),
                        labels = c('Lineage A', 'Lineage B')) +
      labs(x = '', y = '') +
      theme_bw() +
      geom_scatterpie_legend(ddf$rad,
                             x=10.4, y=1.2, n=3, 
                             labeller=function(x) 10^(10*x)) +
      theme(panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            axis.text.x = element_text(angle = 30, hjust = 1)) +
      scale_x_continuous(breaks = 1:10, labels = ddf$label[1:10],
                         limits = c(0.6, 11.5)) +
      scale_y_continuous(breaks = 1:2, labels = c('C8782T','T28144C'))
    
    pdf(file = 'Output/mut_coverage_2.pdf', width = 8)
    print(p)
    dev.off()
  }
  
  return(ddf)
}

ddf = plotfun()

ddf[10,'p'] = 1
yy = sapply(ddf[1:10,1], function(x){
  y = ddf[ddf$Sample == x,'p']
  if(y[1]>=0.5 && y[2]>=0.5){
    return('Lineage A')
  }
  if(y[1]<0.5 && y[2]<0.5){
    return('Lineage B')
  }
  return(ddf[ddf$Sample == x,'clade'][1])
})

df3[df3$Sample %in% names(yy),'clade'] = rep(yy,2)

df3$rad = log10(df3$total_reads)/10
df3$rad[df3$rad == -Inf] = 0.1

values = c(hue_pal()(3)[1], hue_pal()(3)[3])
set = names(table(df3$Names_ID))[table(df3$Names_ID) > 2]
dflist = list()
for (j in 1:length(set)) {
  
  df = df3[df3$Names_ID == set[j],]
  df = df[order(df$Date),]
  nn = nrow(df)/2
  # between locus
  pp1 = abs(df$p[(1:nn)*2]-df$p[(1:nn)*2-1])
  nn1 = which(!is.na(pp1) & pp1<0.1)
  df = df[sort(c(nn1*2, nn1*2-1)),]
  nn = nrow(df)/2
  if(nn <= 2){
    next
  }
  # between sample at 8782
  pp2 = abs(df$p[(2:nn)*2-1]-df$p[(1:(nn-1))*2-1])
  # between sample at 28144
  pp3 = abs(df$p[(2:nn)*2]-df$p[(1:(nn-1))*2])
  
  if (max(pp2) < 0.1 | max(pp3) < 0.1) {
    next
  }
  
  dflist[[set[j]]] = df
}

pdf(file = 'Output/2ES5C_mut_coverage_dynamic.pdf', width = 3.2, height = 2.2)


for (j in 1:length(dflist)) {
  
  df = dflist[[j]]
  df$region = 1:nrow(df)
  nn = nrow(df)/2
  bb = seq(0,5,5/(nn+1)); bb = bb[2:(length(bb)-1)]
  df$x = rep(bb,each = 2)
  df$y = rep(c(1,2), nn)
  t = 2.8
  df$y2 = df$p + t
  t2 = 0.6
  df$pos = as.factor(df$pos)
  
  dftext = df %>% group_by(x) %>% 
    summarise(y = t+mean(p), p = mean(p))
  v <- ifelse(dftext$p > 0.5, 2, -1)
  
  p = ggplot(data = df) + 
    geom_line(aes(x=x,y=y2, group=pos, color = pos)) +
    geom_point(aes(x=x,y=y2, group=pos, color = pos)) +
    scale_color_manual(name="SNP",
                      values = alpha(c('firebrick3', 'orange'),0.8),
                      labels = c('Lineage A', 'Lineage B'),
                      guide = 'none') +
    new_scale_color() +
    # new_fill_color() +
    geom_scatterpie(data = df,
                    aes(x=x, y=y, group=region, r = rad), 
                    cols=c('b2','b1')) + 
    coord_equal() +
    scale_fill_manual(name="SNP",
                      values = alpha(values,0.8),
                      labels = c('Lineage A', 'Lineage B'),
                      guide = 'none') +
    labs(x = '', y = '') + ggtitle(names(dflist)[j])+
    theme_bw() +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())+
          # axis.text.x = element_text(angle = 30, hjust = 1)) +
    scale_x_continuous(breaks = bb, labels =  format(as.Date(df$Date), "%b %d")[seq(1,nn*2,2)],
                       limits = c(0.6, 4+0.5)) +
    scale_y_continuous(breaks = c(1:2,t,t+t2,t+2*t2), 
                       labels = c('8782','28144',0.0,0.5,1.0),
                       limits = c(0.6, t+2*t2)) +
    geom_text(data=data.frame(dftext),
              aes(label = sprintf("%.2f", p),
              x=x,y=y), size = 3,
              inherit.aes = F,
              vjust = v) 
  
  print(p)

  
  
}

dev.off()



pdf(file = 'Output/2ES5C_mut_coverage_dynamic_legend.pdf', width = 5, height = 2.2)

p = ggplot(data = df) + 
  geom_line(aes(x=x,y=y2, group=pos, color = pos)) +
  geom_point(aes(x=x,y=y2, group=pos, color = pos)) +
  scale_color_manual(name="Mutations",
                     values = alpha(c('firebrick3', 'orange'),0.8),
                     labels = c('Lineage A', 'Lineage B')) +
  new_scale_color() +
  # new_fill_color() +
  geom_scatterpie(data = df,
                  aes(x=x, y=y, group=region, r = rad), 
                  cols=c('b2','b1')) + 
  coord_equal() +
  scale_fill_manual(name="Clades",
                    values = alpha(values,0.8),
                    labels = c('Lineage A', 'Lineage B')) +
  labs(x = '', y = '') + ggtitle(names(dflist)[j])+
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  # axis.text.x = element_text(angle = 30, hjust = 1)) +
  scale_x_continuous(breaks = bb, labels =  format(as.Date(df$Date), "%b %d")[seq(1,nn*2,2)],
                     limits = c(0.6, 4+0.5)) +
  scale_y_continuous(breaks = c(1:2,t,t+t2,t+2*t2), 
                     labels = c('8782','28144',0.0,0.5,1.0),
                     limits = c(0.6, t+2*t2)) +
  geom_text(data=data.frame(dftext),
            aes(label = sprintf("%.2f", p),
                x=x,y=y), size = 3,
            inherit.aes = F,
            vjust = v) +
  geom_scatterpie_legend(ddf$rad,
                         x=2, y=3.6, n=3, 
                         labeller=function(x) 10^(10*x))
print(p)
dev.off()
