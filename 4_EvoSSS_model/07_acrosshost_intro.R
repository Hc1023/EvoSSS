rm(list = ls())
library(rstan)
library(ggplot2)
library(scales)
library(gridExtra)
library(reshape2)

realData_all <- read.csv("Covid19CasesWH.csv", row.names = 1)
CaseNum = realData_all$CaseNum
I0 = sum(CaseNum[22:24])
load(file = 'evoSIR.rdata')

df = data.frame()
for (i in 1:9) {
  fit = fitlist[[i]]
  posterior_samples <- rstan::extract(fit)
  beta1 = posterior_samples$beta1
  beta2 = posterior_samples$beta2
  gamma = posterior_samples$gamma
  R1 = mean(beta1/gamma)
  R2 = mean(beta2/gamma)
  ratio = i/10
  I10 = round(sum(CaseNum[22:24])*ratio)
  I20 = round(sum(CaseNum[22:24])*(1-ratio))
  time_intro1 = log(I10)/(beta1-gamma) # - log(m)
  time_intro2 = log(I20)/(beta2-gamma) # - log(m)
  delta = time_intro1-time_intro2
  ci_lower = quantile(delta, probs = 0.025, na.rm = T)
  ci_upper = quantile(delta, probs = 0.975, na.rm = T)
  
  df = rbind(df, data.frame(A = R1, B = R2, 
                            delta = mean(delta), 
                            group = ratio,
                            ci_lower = ci_lower,
                            ci_upper = ci_upper))
}

df$group = as.factor(df$group)

## barplot
p = ggplot(df, aes(x = group, y = delta, fill = delta)) +
  geom_bar(stat = "identity", position = position_dodge(), width = 0.7) +
  geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper), 
                width = 0.2) +
  theme_minimal() +
  xlab('') + ylab('')+
  theme(legend.key.size = unit(0.2, 'cm')) +
  scale_y_continuous(breaks = c(-10,0,10), labels = c('-10','0','+10')) +
  geom_hline(yintercept = 0) + 
  scale_fill_gradient2(low = "#7c615c", mid = "white", high = "#59788e", 
                       midpoint = 0, breaks = c(-10,0,10), 
                       labels = c('-10','0','+10'),
                       name = expression(d[B]-d[A]))
p
pdf(paste0("Output/R0_delta.pdf"), width = 4.5, height = 1)
print(p)
dev.off()


data_melted <- melt(df, 
                    id.vars = c("group"),
                    measure.vars = c("A", "B"))

## dotplot
p <- ggplot(data_melted, 
            aes(x = group, y = variable, 
                size = value, fill = value)) +
  geom_point(shape = 21) +
  scale_size_continuous(range = c(1, 12), guide = 'none') + 
  theme_minimal() +
  labs(x = "", y = "", size = "R0") +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.position = "right",
        legend.key.size = unit(0.2,'cm'),
        panel.background = element_blank()) +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", 
                       midpoint = 2.5, name = expression(R[0]))

p

pdf(paste0("Output/R0.pdf"), width = 4.5, height = 1.45)
print(p)
dev.off()


