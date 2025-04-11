rm(list = ls())
library(rstan)
library(ggplot2)
library(ggpubr)

load(file = 'evoSIR/F4D_evoSIR.rdata')
fit = fitlist[[4]]
posterior_samples <- rstan::extract(fit)

c = posterior_samples$beta1/posterior_samples$beta2
dat = data.frame(y = -log(c))

p = ggplot(dat, aes(x=y)) + 
  geom_histogram(binwidth=0.001, color = "#2e6bc5", 
                 fill= alpha("skyblue",0.3)) +
  geom_boxplot(width=800, aes(y = 1500),
               outlier.alpha = 0.1)  +
  theme_bw() +
  xlab(expression(r[B]-r[A])) + ylab('Counts (1e+3)') +
  scale_y_continuous(breaks = seq(0,3000, 1000), labels = seq(0,3)) + 
  scale_x_continuous(limits = range(dat$y), n.breaks = 5,
                     sec.axis = sec_axis(
                       transform = ~ exp(.*(-1)), 
                       name = expression(beta[A]/beta[B])
                     )) +
  theme(panel.grid.minor = element_blank())
p

pdf(file = 'Output/F4E_deltar.pdf', width = 2.2, height = 2.1)
print(p)
dev.off()
