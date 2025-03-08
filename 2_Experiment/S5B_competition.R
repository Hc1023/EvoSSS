rm(list = ls())
library(data.table)
library(tidyverse)
library(ggpubr)
library(zoo)

df_all = read.csv("S5A_vcf_231228.csv")
sampleid = read.csv("S4B_sampleid.csv")
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


pdf(paste0("Output/S5B_competition.pdf"), width = 4.3, height = 2.8)

for (i in c('Vero', 'Calu')) {
  df = myfun(i)
  h = ggplot(df, aes(x = Time, y = VAF, 
                     group = interaction(Mutations, Replicate))) +
    geom_line(aes(color = Mutations)) +
    geom_point(aes(shape = Replicate, color = Mutations)) +
    facet_wrap( ~ Strain, scales = "free_y",
                ncol = 3, strip.position = "top", 
                drop = FALSE) +
    labs(title = "",
         x = "Hours / Passaged",
         y = "") +
    scale_color_manual(values = alpha(c('firebrick3', 'orange'), 0.8)) +
    scale_fill_manual(values = alpha(c('firebrick3', 'orange'), 0.8)) +
    theme_minimal() +
    scale_y_continuous(limits = c(0, 1), 
                       minor_breaks = seq(0 , 1, 0.25),n.breaks = 3) +
    theme(#axis.text.x = element_text(angle = 90, hjust = 1),
          axis.text.y = element_blank(),
          strip.background = element_blank(),  # Remove strip background
          strip.text = element_text(face = "bold", size = 8),
          legend.position = "top", 
          legend.margin = margin(t = -10, r = 0, b = -10, l = 0)) +
    guides(color = guide_legend(override.aes = list(shape = NA)),
           shape = 'none')
  print(h)
}

h = ggplot(df, aes(x = Time, y = VAF, 
                   group = interaction(Mutations, Replicate))) +
  geom_line(aes(color = Mutations)) +
  geom_point(aes(shape = Replicate, color = Mutations)) +
  facet_wrap( ~ Strain, scales = "free_y",
              ncol = 3, strip.position = "top", 
              drop = FALSE) +
  labs(title = "",
       x = "Hours / Passaged",
       y = "VAF") +
  scale_color_manual(values = alpha(c('firebrick3', 'orange'), 0.8)) +
  scale_fill_manual(values = alpha(c('firebrick3', 'orange'), 0.8)) +
  theme_minimal() +
  scale_y_continuous(limits = c(0, 1), 
                     minor_breaks = seq(0 , 1, 0.25),n.breaks = 3) +
  theme(strip.background = element_blank(),  # Remove strip background
    strip.text = element_text(face = "bold", size = 8),
    legend.position = "top", 
    legend.margin = margin(t = -10, r = 0, b = -10, l = 0)) +
  guides(color = guide_legend(override.aes = list(shape = NA)),
         shape = 'none')
print(h)
dev.off()
