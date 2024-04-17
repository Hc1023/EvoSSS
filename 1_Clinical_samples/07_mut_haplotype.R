rm(list = ls())
library(scatterpie)
library(RColorBrewer)
library(ggplot2)
library(scales)

df3 = read.csv('clinical_coverage.csv')
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
    
    pdf(file = 'Output/haplotype.pdf', width = 8)
    print(p)
    dev.off()
  }
  
  return(ddf)
}

ddf = plotfun(draw = T)
