rm(list = ls())
library(rstan)
library(ggplot2)
library(ggpubr)

load(file = 'SIR.rdata')
fit = fitlist[[5]]
posterior_samples <- extract(fit)

c = posterior_samples$beta1/posterior_samples$beta2
dat = data.frame(x = -log(c))

p = ggplot(dat, aes(x=x)) + 
  geom_histogram(binwidth=0.001, color = "#2e6bc5", fill= alpha("skyblue",0.3)) +
  geom_boxplot(width=800, aes(y= 1500), 
               outlier.alpha = 0.1)  +
  theme_bw() +
  xlab(expression(-Delta*r)) + ylab('') +
  scale_x_continuous(limits = c(0.072,0.093), breaks = seq(0.07,0.09,0.01))
p
pdf(file = 'Output/negdeltar.pdf', width = 2, height = 1.2)
print(p)
dev.off()
