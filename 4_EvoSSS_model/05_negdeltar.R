rm(list = ls())
library(rstan)
library(ggplot2)
library(ggpubr)

load(file = 'evoSIR.rdata')
fit = fitlist[[4]]
posterior_samples <- extract(fit)

c = posterior_samples$beta1/posterior_samples$beta2
dat = data.frame(y = -log(c))

p = ggplot(dat, aes(y=y)) + 
  geom_histogram(binwidth=0.001, color = "#2e6bc5", fill= alpha("skyblue",0.3)) +
  geom_boxplot(width=800, aes(x = 1500),
               outlier.alpha = 0.1)  +
  theme_bw() +
  ylab(expression(-Delta*r)) + xlab('') +
  scale_x_continuous(breaks = seq(0,3000, 1000), labels = seq(0,3)) + 
  scale_y_continuous(limits = range(dat$y), n.breaks = 5,
                     sec.axis = sec_axis(
                       transform = ~ exp(.*(-1)), 
                       name = ''
                     )) 
p

pdf(file = 'Output/negdeltar.pdf', width = 2.1, height = 2)
print(p)
dev.off()
