rm(list = ls())
library(reshape2) 
library(ggplot2) 
library(grDevices) 
library(RColorBrewer)
library(directlabels) 
library(tidyverse)

load('S8A_bottleneck_stochasticity.Rdata')
rf <- colorRampPalette(rev(brewer.pal(11,'Spectral')))
colormap <- colorRampPalette(rev(brewer.pal(11,'Spectral')))(32)
dfvar_all = data.frame(dfvar_all)
dfvar_all$varlog = log10(dfvar_all$var)
dfvar_all$varlogr = round(dfvar_all$varlog, digits = 1)
dfvar_all$varlogfit = round(dfvar_all$varlog, digits = 0)

dflines = data.frame()
for (i in unique(dfvar_all$varlogfit)) {
  df1 = dfvar_all[dfvar_all$varlogfit == i,]
  df1$h = as.numeric(df1$h)
  fitlm = lm(seed~I(h^(-1)), data=df1)
  dfline = data.frame(h = unique(df1$h), var = i)
  dfline$predy = predict(fitlm, newdata = dfline)
  dflines = rbind(dflines, dfline)
}
dfvar_all$h = as.numeric(dfvar_all$h)
dflines$var = factor(dflines$var, levels = -c(1:5))
p = ggplot() +
  geom_tile(data = dfvar_all,
            aes(x=h, y=seed, fill = varlogr)) +
  scale_fill_gradientn(colours= alpha(colormap,1), 
                       breaks = -c(1:5),
                       labels = paste0('1e-',1:5),
                       name = 'Variance') +
  geom_line(data = dflines[dflines$var != -5,], 
            aes(x = h, y = predy, group = var, color = var)) +
  scale_color_manual(values = c('#a20b1d','#cd5b3e',
                                '#bfc977','#356f76'),
                     guide = 'none') +
  scale_y_continuous(limits = c(0,100),expand=c(0,0)) +
  scale_x_continuous(limits = c(0,100),expand=c(0,0)) +
  theme_classic() +
  xlab('Number of hotspots (h)') + ylab('Seed pop. (s)') +
  theme(plot.margin = margin(t = 1, r = 1, b = 0, l = 0, 
                             unit = "cm"), 
        legend.key.size = unit(0.25,'cm'),
        legend.spacing = unit(0.0,'cm'),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 10,
                                    margin = margin(2, 0, 2, 0)),
        legend.margin = margin(0, 0, 0, 0)
        )
p  

pdf(file = 'Output/S8A_bottleneck_stochasticity.pdf', width = 2.65, height = 1.7)
print(p)
dev.off()
