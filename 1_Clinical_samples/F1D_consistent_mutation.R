rm(list= ls())
library(dplyr)
library(tidyverse)
library(ggplot2)
library(scales)
library(knitr)
library(RColorBrewer)
library(ggsci)

df = read.csv('F1C_clinical_coverage.csv')
# b1: lineage B b2: lineage A
df1 = df[df$pos == 8782,
         c('Sample', 'Names_ID','total_reads',
           'b1','b2','n')]
df2 = df[df$pos == 28144,
         c('Sample', 'Names_ID','total_reads',
           'b1','b2','n')]

df1$b = df1$b1+df1$b2
df1$ratio = df1$b2/df1$b
df1mix = df1[df1$ratio > 0.1 & df1$ratio < 0.9 & !is.na(df1$ratio),
             c('Sample', 'b', 'ratio')]
df2$b = df2$b1+df2$b2
df2$ratio = df2$b2/df2$b
df2mix = df2[df2$ratio > 0.1 & df2$ratio < 0.9 & !is.na(df2$ratio),
             c('Sample', 'b', 'ratio')]

dfmerge = merge(df1mix, df2mix, by = "Sample")
dfmerge[,'group'] = '>= 0.2'
dfmerge[abs(dfmerge$ratio.x-dfmerge$ratio.y) < 0.2,'group'] = '0.1-0.2' 
dfmerge[abs(dfmerge$ratio.x-dfmerge$ratio.y) < 0.1,'group'] = '< 0.1' 
dfmerge$group = factor(dfmerge$group, 
                       levels = c('< 0.1', '0.1-0.2', '>= 0.2'))
dfmerge$shape = 'Other'
dfmerge$shape[dfmerge$b.x>1000 & dfmerge$b.y>1000] = '> 1000'
x_min = 0
x_max = 1
# |y - x| < 0.1
band_0.1 <- data.frame(
  x = c(x_min, x_max, x_max, x_min),
  y = c(x_min + 0.1, x_max + 0.1, x_max - 0.1, x_min - 0.1)
)
# 0.1 < |y - x| < 0.2
band_0.2 <- data.frame(
  x = c(x_min, x_max, x_max, x_min),
  y = c(x_min + 0.2, x_max + 0.2, x_max - 0.2, x_min - 0.2)
)

p = ggplot(dfmerge, aes(x = ratio.x, y = ratio.y, 
                        color = group)) +
  labs(x = "VAF of C8782T", 
       y = "VAF of T28144C", 
       color = "Difference") +
  theme_bw() +
  theme(panel.grid.minor = element_blank()) +
  coord_cartesian(xlim = c(0.2,0.8), ylim = c(0.2,0.8)) +
  geom_polygon(data = band_0.1, aes(x = x, y = y), 
               fill = "grey20", alpha = 0.1, inherit.aes = FALSE) +
  geom_polygon(data = band_0.2, aes(x = x, y = y), 
               fill = "grey60", alpha = 0.1, inherit.aes = FALSE) +
  geom_abline(intercept = c(-0.1, 0.1), slope = 1, linetype = "longdash", color = "grey20") +
  geom_abline(intercept = c(-0.2, 0.2), slope = 1, linetype = "longdash", color = "grey60") +
  geom_point(size = 1, alpha = 0.8) +
  scale_color_manual(values = c("< 0.1" = "red", 
                                "0.1-0.2" = "black", 
                                ">= 0.2" = "grey")) 
  
pdf(file = 'Output/F1D_consistent_mutation.pdf', width = 3, height = 1.8)
print(p)
dev.off()

# > table(dfmerge$group)
# 
# < 0.1 0.1-0.2  >= 0.2 
# 15       9       4 
# > table(dfmerge$shape)
# 
# > 1000  Other 
# 8     20 
