rm(list = ls())
library(treeio)
library(ggtree)
library(ggplot2)
library(ape)
library(scales)
library(RColorBrewer)
library(ggnewscale)
library(tidyverse)
library(ggpubr)

sample_snpdate = read.csv('clinical_snpdate.csv')
dattmp = read.csv("clinical_info.csv")

dattmp$Date = as.Date(dattmp$Date)
df = dattmp[,-6]
df$sample_id = sapply(dattmp[,5], function(x){
  strsplit(x,'-')[[1]][8]
}) %>% unname()
  
df$Sample = paste0('S', dattmp$Sample)
nextclade <- read.table(file = 'nextclade_qc_early.tsv', 
                        sep = '\t', header = TRUE)
nextclade <- nextclade[nextclade$seqName %in% df$sample_id, c(1,2,3,16)]

for (i in 1:nrow(nextclade)) {
  if(any(grep('C8782T', nextclade[i,4])) & any(grep('T28144C', nextclade[i,4]))){
    nextclade[i,5] <- 'Lineage A'
  }else if(any(grep('C8782T', nextclade[i,4]))){
    nextclade[i,5] <- '8782T'
  }else if(any(grep('T28144C', nextclade[i,4]))){
    nextclade[i,5] <- '28144C'
  }else{
    nextclade[i,5] <- 'Lineage B'
  }
}

df$clade = sapply(df$sample_id,function(x){
  nextclade[nextclade$seqName == x,5]
})

myfun <- function(x){
  df$total_reads = as.numeric(sample_snpdate[sample_snpdate$Mut == x, (1:nrow(df))*3])
  n = as.numeric(sample_snpdate[sample_snpdate$Mut == x, (1:nrow(df))*3+2])
  b2 = as.numeric(sample_snpdate[sample_snpdate$Mut == x, (1:nrow(df))*3+1])
  df$b1 = df$total_reads - n - b2
  df$b2 = b2
  df$n = n
  df1 = df[order(df$total_reads, df$b2, decreasing = T),]
  return(df1)
}

df1 = myfun('C8782T')
df2 = myfun('T28144C')

df1$pos = 8782
df2$pos = 28144
df3 = rbind(df1,df2)
df3$label = paste(df3$Names_ID, 
                  format(df3$Date, "%b %d"))
df3$p = df3$b2/df3$total_reads

write.csv(df3, file = '../data/ZJU_sample/mut_coverage.csv', 
          row.names = F)

# colors = c('#d98276','#bf6d39')
p = ggboxplot(df3, x = 'pos', y = "total_reads", 
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
