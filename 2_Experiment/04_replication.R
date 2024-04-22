rm(list = ls())
library(data.table)
library(tidyverse)
library(readxl)
library(ggpubr)
library(scales)
library(rstatix)

df = read.csv("calu3_ct.csv")

get_df_long = function(df1){
  colnames(df1) = c('Strain', '0', '8', '24', '48', '72')
  
  df2 = df1[-c(1:10), ] %>%
    group_by(Strain) %>%
    summarise_at(vars(-group_cols()), mean, na.rm = TRUE) %>%
    as.data.frame()
  
  df2 = rbind(df1[1:10, ], df2)
  df2$Strain <- factor(c(rep('A', 6), rep('B', 4), rep('A+B', 6)),
                       levels = c('A', 'B', 'A+B'))
  
  df_long <- gather(df2, key = "Time",
                    value = "Ct",-Strain)
  df_long$Time = factor(df_long$Time, levels = unique(df_long$Time))
  return(df_long)
}


get_signif = function(t) {
  subset_data <- df_long[df_long$Time == t, ]
  anova_result <- aov(Ct ~ Strain, data = subset_data)
  p_value = summary(anova_result)[[1]]$`Pr(>F)`[1]
  signif_code = ifelse(p_value < 0.001, "***",
                       ifelse(p_value < 0.01, "**",
                              ifelse(
                                p_value < 0.05, "*",
                                ifelse(p_value < 0.1, ".", " ")
                              )))
  return(paste0('p=',paste(format(round(p_value, 3), nsmall = 3), 
                           signif_code)))
}
gene = c('ORF1ab', 'N')
pdf(paste0("Output/replication.pdf"), width = 2.7, height = 2.3)

for (i in 1:2) {
  if (i==1) {
    df1 = df[, 1:6]
  }else{
    df1 = df[, c(1, 9:13)]
  }
  df_long = get_df_long(df1)
  x= sapply(c('0', '8', '24', '48', '72'), get_signif)
  y = df_long %>% group_by(Time) %>%
    summarise_at(vars(Ct), function(c){
      ifelse(min(c) > 25, min(c) - 2, max(c) + 0.2)
    })
  
  annotation_df <- data.frame(Time = names(x), Annotation = x, 
                              y = y$Ct)
  values = c(hue_pal()(3)[1], hue_pal()(3)[3], '#aa85a6')
  p = ggplot(df_long, aes(x = Time, y = Ct, fill = Strain)) +
    geom_boxplot(aes(color = Strain)) +
    geom_point(
      position = position_jitterdodge(),
      show.legend = FALSE,
      aes(color = Strain),
      alpha = 0.5,
      size = 1
    ) +
    scale_color_manual(values = values) +
    scale_fill_manual(values = alpha(values, 0.3)) +
    labs(x = "Time points (h)",
         y = paste(gene[i], "(Ct)")) +
    theme_bw() +
    theme(
      legend.position = c(0.85, 0.7),
      legend.background = element_rect(color = NA, fill = NA),
      legend.key = element_blank()
    ) +
    annotate("text", x = annotation_df$Time, y = annotation_df$y,
             label = x, vjust = -0.5, size = 2.8) +
    coord_cartesian(ylim = c(12, 33.9))
  
  print(p)
  
}
dev.off()
