rm(list = ls())
library(dplyr)
library(ggplot2)
library(scales)
library(RColorBrewer)
if(F){
  df_filter0 <- read.csv('Output/metadata/metadata_part0_1.csv')
  df_filter1 <- read.csv('Output/metadata/metadata_part1_1.csv')
  df_filter2 <- read.csv('Output/metadata/metadata_part2_1.csv')
  df_filter2.2 <- read.csv('Output/metadata/metadata_part2.2_1.csv')
  df_filter3 <- read.csv('Output/metadata/metadata_part3_1.csv')
  df_filter4 <- read.csv('Output/metadata/metadata_part4_1.csv')
  metadata <- rbind(df_filter0, df_filter1, df_filter2,
                    df_filter2.2, df_filter3, df_filter4)
  
  metadata$date = as.Date(metadata$date)
  table(metadata$GISAID_clade)
  
  metadata = metadata[metadata$GISAID_clade %in% c('L','S','G','GRY','GK'),]
  metadata = metadata[!is.na(metadata$date),]
  metadata$x = format(metadata$date, '%Y-%m')
  result <- metadata %>%
    group_by(x, GISAID_clade) %>%
    summarise(count = n()) %>% as.data.frame()
  
  
  write.csv(result, 'VOC_gisaid.csv', row.names = F)
}

if(F){
  result = read.csv('VOC_gisaid.csv')
  VOC = c('D614G', 'A', 'B', 'Delta','Alpha')
  names(VOC) = unique(result$GISAID_clade)
  result$V = NA
  for (i in 1:length(VOC)) {
    result$V[result$GISAID_clade == names(VOC)[i]] = VOC[i]
  }
  result$x = as.Date(paste0(result$x,'-01'))
  VOC = VOC[c(4,5,1,3,2)]
  VOC = c("Delta","Alpha","D614G","B","A")
  result$V = factor(result$V, levels = VOC) 
  result = result[result$x <= as.Date('2021-10-31'),]
  write.csv(result, 'VOC_gisaid.csv', row.names = F)
  
}

result = read.csv('VOC_gisaid.csv')
VOC = c("Delta","Alpha","D614G","B","A")
result$V = factor(result$V, levels = VOC)
result$x = as.Date(result$x)
values = c('#b2360d','#81e880','#fddd10',hue_pal()(3)[c(3,1)])
show_col(brewer.pal(n=4,name = "Spectral"))
brews = brewer.pal(n=4,name = "Spectral")
values = c(brews[2:4], hue_pal()(3)[c(3,1)])

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

pdf(paste0("Output/VOC_gisaid.pdf"), width = 4.2, height = 1.8)
print(p)
dev.off()


