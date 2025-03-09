rm(list = ls())
library(ggplot2)
library(scales)
library(ggnewscale)
library(ggpubr)
library(tidyverse)

migration_distribution <- function(qm = 0.5, 
                                   k = c(1/4,1/2,1,2,4)){
  
  phi_q = function(q, qm = 0.5, k){
    # k = m*N
    g = gamma(4*k)/gamma(4*k*qm)/gamma(4*k*(1-qm))
    y = g*q^(4*k*qm-1)*(1-q)^(4*k*(1-qm)-1)
    return(y)
  }
  
  x = 1:99/100
  
  y <- sapply(k, function(k) {
    sapply(x, function(x_val) {
      phi_q(x_val, qm, k)
    })
  })
  
  dat = data.frame(q = x,
                   y)
  
  colnames(dat)[2:ncol(dat)] = as.character(k)
  dat_gather <- gather(dat, k, y, -q)
  dat_gather$k = factor(dat_gather$k, 
                        levels = rev(as.character(k)))
  
  return(dat_gather)
}
qm = 0.5
dat_gather1 = migration_distribution(qm = qm)
p1 = ggplot(dat_gather1, 
            aes(x = q, y = y, group = k, color = k)) +
  geom_line() + theme_bw() +
  geom_vline(xintercept = qm, linetype = "dashed") +
  ylab('Density') + xlab('') +
  theme(legend.position = 'right',
        legend.background = element_rect(color = NA, fill = NA),
        legend.key=element_blank(),
        legend.spacing.y = unit(0.05,'cm'),
        legend.key.size = unit(0.3,'cm')) +
  scale_color_manual(name = '',
                     values = c( '#cb8335','#a3a637',
                                 '#68a588','#579aa7','#8781ba'))
p1
qm = 0.3
dat_gather2 = migration_distribution(qm = qm)
p2 = ggplot(dat_gather2, 
            aes(x = q, y = y, group = k, color = k)) +
  geom_line() + theme_bw() +
  geom_vline(xintercept = qm, linetype = "dashed") +
  ylab('Density') + xlab('') +
  theme(legend.position = 'right',
        legend.background = element_rect(color = NA, fill = NA),
        legend.key=element_blank(),
        legend.spacing.y = unit(0.05,'cm'),
        legend.key.size = unit(0.3,'cm')) +
  scale_color_manual(name = 'k',
                     values = c( '#cb8335','#a3a637',
                                 '#68a588','#579aa7','#8781ba'))


pdf(file = paste0('Output/F4G_bottleneck_phi.pdf'), width = 2.5, height = 1.5)
print(p1)
print(p2)
dev.off() 


