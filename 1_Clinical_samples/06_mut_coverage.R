rm(list = ls())
library(ggpubr)

df = read.csv('clinical_coverage.csv')

p = ggboxplot(df, x = 'pos', y = "total_reads", 
          fill = "pos", color = 'pos', 
          scales = "free_y", ncol = 4,
          outlier.shape = 21, width = 0.7) +
  aes(alpha = 0.5) +
  geom_jitter(size=0.6, alpha=0.8, shape = 16, 
              width = 0.1, height = 0.01, aes(color = pos) 
              ) +
  theme_bw() + theme(legend.position="none") +
  labs(y = "Coverage", x = "") +
  scale_y_continuous(trans=scales::pseudo_log_trans(base=2),
                     breaks = round(c(0, 10^seq(0,5,1)))) #yscale("log10", .format = TRUE)

p
pdf(file = 'Output/mut_coverage.pdf', width = 2.5, height = 2.5)
print(p)
dev.off()
