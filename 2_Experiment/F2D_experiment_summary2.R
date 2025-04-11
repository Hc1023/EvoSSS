rm(list = ls())
library(data.table)
library(tidyverse)
library(ggpubr)
library(zoo)
library(ggthemes)
library(scales)
library(ggsci)

df_all = read.csv("S8A_vcf_231228.csv")
sampleid = read.csv("S5C_sampleid.csv")
sampleid[, 2][sampleid[, 2] == ""] <- NA
sampleid[,2] = na.locf(sampleid[,2])

sampleid1 = sampleid[sampleid$original != '',c(1,2,3)]
sampleid2 = sampleid[,-2]

myfun = function(x){
  v = sampleid[sampleid$cell == x,-3]
  v = v[11:28,-c(3,4)]
  colnames(v) = c('Strain', '0', '24', '48', '72', '1st', '2nd')
  v$Replicate = rep(c('1', '2','3'), times = 6)
  gathered_data <- v %>%
    gather(key = "Time", value = "Value", -Strain, -Replicate)
  gathered_data[,c('C8782T','T28144C')] = NA
  for (i in 1:nrow(gathered_data)) {
    df = df_all[df_all$sample== gathered_data[i,4],]
    df$MUT = paste0(df$REF, df$POS, df$ALT)
    gathered_data[i,5] <- ifelse(sum(df$MUT == 'C8782T') == 0, 0, 
                                 df[df$MUT == 'C8782T', 'percent'])
    gathered_data[i,6] <- ifelse(sum(df$MUT == 'T28144C') == 0, 0, 
                                 df[df$MUT == 'T28144C', 'percent'])
  }
  df <- gathered_data[,-4] %>%
    gather(key = "Mutations", value = "VAF", -Strain, -Replicate, -Time)
  
  df$Time <- factor(df$Time, colnames(v)[2:7])
  df$Strain <- factor(df$Strain, unique(df$Strain)[c(1,3,5,2,4,6)])
  return(df)
}

hlist = list()
for (i in c('Vero', 'Calu')) {
  df = myfun(i)
  df = df %>%
    group_by(Time, Strain, Replicate) %>%
    summarize(VAF = mean(VAF)) 
  median_df <- df %>%
    group_by(Time) %>%
    summarize(median_VAF = median(VAF))
  
  
  hlist[[i]] = ggplot(df, aes(x = Time, y = VAF)) +
    geom_boxplot(outlier.shape = NA, 
                 # color = alpha('orange',0.8),
                 # linewidth = 0.3,
                 fill = alpha(hue_pal()(3)[1],0.3), 
                 color = hue_pal()(3)[1]) +
    geom_jitter(aes(color = Strain), width = 0.2, size = 0.6) + 
    scale_color_manual(values = alpha(pal_lancet()(6)[6:1], 0.5)) + 
    geom_line(data = median_df, 
              aes(x = Time, y = median_VAF, group = 1),
              color = hue_pal()(3)[1]) +  # Add median line
    geom_boxplot(df, mapping = aes(x = Time, y = 1-VAF),
                 outlier.shape = NA, 
                 fill = alpha(hue_pal()(3)[3],0.3), 
                 color = hue_pal()(3)[3]) +
    geom_line(data = median_df, 
              aes(x = Time, y = 1-median_VAF, group = 1),
              color = hue_pal()(3)[3]) +  # Add median line
    theme_classic() +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          legend.background = element_blank(),
          plot.background = element_blank(),
          strip.background = element_blank(),
          legend.key.size = unit(0.4, 'cm')
          ) +
    scale_y_continuous(breaks = c(0,0.5,1)) +
    scale_x_discrete(labels = c('P0','24','48','72','P1','P2')) +
    coord_cartesian(ylim = c(0,1))
}

pdf(paste0("Output/F2D_experiment_summary2.pdf"), width = 3.2, height = 1.6)
print(hlist[[1]])
print(hlist[[2]])
dev.off()


