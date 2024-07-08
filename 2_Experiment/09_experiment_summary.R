rm(list = ls())
library(data.table)
library(tidyverse)
library(ggpubr)
library(zoo)
library(ggthemes)
df_all = read.csv("vcf_231228.csv")
sampleid = read.csv("sampleid.csv")
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
myfun2 = function(x){
  v = sampleid[sampleid$cell == x,-3]
  v = v[1:10,-c(3,4)]
  colnames(v) = c('Strain', '0', '24', '48', '72', '1st', '2nd')
  v$Replicate = rep(c('1', '2'), times = 5)
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
  
  return(df)
}


pdf(paste0("Output/experiment_summary.pdf"), width = 1.4, height = 1)

for (i in c('Vero', 'Calu')) {
  df = myfun(i)
  # df =df[df$Time %in% unique(df$Time)[1:4],]
  median_df <- df %>%
    group_by(Time) %>%
    summarize(median_VAF = median(VAF))
  h=ggplot(df, aes(x = Time, y = VAF)) +
    geom_boxplot(outlier.shape = NA, 
                 # color = alpha('orange',0.8),
                 linewidth = 0.1,
                 fill = alpha('orange',0.6)) +
    geom_line(data = median_df, 
              aes(x = Time, y = median_VAF, group = 1)) +  # Add median line
    theme_bw() +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    scale_y_continuous(breaks = c(0,0.5,1)) +
    coord_cartesian(ylim = c(0,1))
  
  print(h)
}

for (i in c('Vero', 'Calu')) {
  df = myfun2(i)
  df = df[grep('A',df$Strain),]
  # df =df[df$Time %in% unique(df$Time)[1:4],]
  median_df <- df %>%
    group_by(Time) %>%
    summarize(median_VAF = median(VAF))
  h=ggplot(df, aes(x = Time, y = VAF)) +
    geom_boxplot(outlier.shape = NA, 
                 linewidth = 0.1,
                 fill = alpha('orange',0.6)) +
    geom_line(data = median_df, 
              aes(x = Time, y = median_VAF, group = 1)) +  # Add median line
    theme_bw() +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    scale_y_continuous(breaks = c(0,0.5,1)) +
    coord_cartesian(ylim = c(0,1))
  
  print(h)
}


dev.off()

library(scales)
show_col('#e2891f')
