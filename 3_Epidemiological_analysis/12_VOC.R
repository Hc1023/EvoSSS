rm(list = ls())
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(circlize)
library(dplyr)
library(ggplot2)
library(scales)
library(RColorBrewer)
dat = readxl::read_excel("Output/metadata/VOC.xlsx")
table(dat$clade)
unique(dat$clade)
dat = dat[!is.na(dat$clade)&!is.na(dat$submission_date),]

df = data.frame(dat[,c('submission_date', 'clade')])
df$submission_date = as.Date(df$submission_date)
df$x = format(df$submission_date, '%Y-%m')

levels = c("22B (Omicron)","22A (Omicron)","21M (Omicron)","21L (Omicron)",
           "21K (Omicron)","21J (Delta)","21I (Delta)", "21H (Mu)","20E (EU1)","21A (Delta)",
           "20I (Alpha, V1)", "20B","20A" )

VOC = c('Omicron','Delta','Alpha','20A')
df$V =  NA
for (i in 1:length(VOC)) {
  v = VOC[i]
  df$V[df$clade %in% levels[grep(v,levels)]] = v
}
df = df[!is.na(df$V) & !is.na(df$x),]

df$V[df$V == '20A'] = 'D614G'
result <- df %>%
  group_by(x, V) %>%
  summarise(count = n()) %>% as.data.frame()
VOC = c('Omicron','Delta','Alpha','D614G')
result$x = as.Date(paste0(result$x,'-01'))
result$V = factor(result$V, levels = VOC) 

values = c('#4148f0','#b2360d','#81e880','#fddd10')
show_col(values)
values = brewer.pal(n=4,name = "Spectral")

p = ggplot() +
  geom_area(data = result, aes(x = x, y = count, fill = V),
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
library(RColorBrewer)

