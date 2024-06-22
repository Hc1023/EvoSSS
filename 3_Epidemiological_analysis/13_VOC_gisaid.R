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
  result <- metadata %>%
    group_by(date, GISAID_clade) %>%
    summarise(count = n()) %>% as.data.frame()
  write.csv(result, 'VOC_gisaid.csv', row.names = F)


}

result = read.csv('VOC_gisaid.csv')
result$date = as.Date(result$date)
result$x = format(result$date, '%Y-%m')
df <- result %>%
  group_by(x, GISAID_clade) %>%
  summarise(count2 = sum(count)) %>% as.data.frame()
df$x = as.Date(paste0(df$x,'-01'))
VOC = c("Delta","Alpha","D614G","B 19A","A")
names(VOC) = c("GK","GRY","G","S","L")
df = df[df$GISAID_clade %in% names(VOC),]
df$V = NA
for (i in 1:length(VOC)) {
  df$V[df$GISAID_clade == names(VOC)[i]] = VOC[i]
}

df2 = df[df$x <= as.Date('2021-10-31'),]
df2$V = factor(df2$V, levels = VOC)

brews = brewer.pal(n=4,name = "Spectral")
values = c(brews[2:4], hue_pal()(3)[c(3,1)])

p = ggplot() +
  geom_area(data = df2, aes(x = x, y = count2, fill = V),
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

pdf(paste0("Output/VOC_gisaid.pdf"), width = 3.5, height = 1.8)
print(p)
dev.off()


