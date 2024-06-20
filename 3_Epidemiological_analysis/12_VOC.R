rm(list = ls())
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(circlize)
library(dplyr)
library(ggplot2)
library(scales)
library(RColorBrewer)

if(F){
  dat = readxl::read_excel("Output/metadata/VOC.xlsx")
  table(dat$clade)
  unique(dat$clade)
  dat = dat[!is.na(dat$clade)&!is.na(dat$submission_date),]
  
  df = data.frame(dat[,c('submission_date', 'clade')])
  df$submission_date = as.Date(df$submission_date)
  result <- df %>%
    group_by(submission_date, clade) %>%
    summarise(count = n()) %>% as.data.frame()
  colnames(result)[1] = 'date'
  levels = c("22B (Omicron)","22A (Omicron)","21M (Omicron)","21L (Omicron)",
             "21K (Omicron)","21J (Delta)","21I (Delta)", "21H (Mu)","20E (EU1)","21A (Delta)",
             "20I (Alpha, V1)", "20B","20A" )
  
  VOC = c('Omicron','Delta','Alpha','20A')
  result$V =  NA
  for (i in 1:length(VOC)) {
    v = VOC[i]
    result$V[result$clade %in% levels[grep(v,levels)]] = v
  }
  result = result[!is.na(result$V),]
  result$V[result$V == '20A'] = 'D614G'
  result = result[,-2]
  
  write.csv(result, 'VOC.csv', row.names = F)
}

result = read.csv('VOC.csv')
result$date = as.Date(result$date)
result$x = format(result$date, '%Y-%m')
df <- result %>%
  group_by(x, V) %>%
  summarise(count2 = sum(count)) %>% as.data.frame()
df$x = as.Date(paste0(df$x,'-01'))
VOC = c('Omicron','Delta','Alpha','D614G')
df$V = factor(df$V, levels = VOC)

values = brewer.pal(n=4,name = "Spectral")

p = ggplot() +
  geom_area(data = df, aes(x = x, y = count2, fill = V),
            position = 'fill') +
  scale_fill_manual(values = values) +
  scale_x_date(breaks = seq(as.Date('2020-01-01'), as.Date('2024-11-01'), by="6 months"),
               minor_breaks = seq(as.Date('2019-12-01'), as.Date('2024-11-01'), by ='1 month'),
               date_labels = "%y-%b",
               expand = c(0, 0)) +
  scale_y_continuous(breaks = seq(0,1,0.25), labels = seq(0,100,25),expand = c(0, 0)) +
  xlab('') +  ylab('Prevalence (%)') + theme_bw() +
  theme(legend.background = element_blank(),
        legend.title = element_blank(),
        plot.margin = margin(t = 10, r = 10, b = 10, l = 10))
p
pdf(paste0("Output/VOC.pdf"), width = 4.2, height = 1.8)
print(p)
dev.off()


