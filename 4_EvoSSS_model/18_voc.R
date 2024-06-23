rm(list = ls())
library(rstan)
library(ggplot2)
library(scales)
library(ggnewscale)
library(tidyverse)
library(dplyr)

df = read.csv('../3_Epidemiological_analysis/VOC_gisaid.csv')

VOC = c("Omicron","Delta","Alpha","D614G")
names(VOC) = c("GRA","GK","GRY","G")
df = df[df$GISAID_clade %in% names(VOC),]
df$V = NA
for (i in 1:length(VOC)) {
  df$V[df$GISAID_clade == names(VOC)[i]] = VOC[i]
}
df$V = factor(df$V, levels = VOC)
dfall = df %>% group_by(date) %>%
  summarise(yall = sum(count))

df <- na.omit(df)

full_dates <- as.Date('2020-06-30') + 1:800
# Create a dataframe with all dates
full_df <- expand.grid(date = full_dates, 
                       V = unique(df$V))
df$date = as.Date(df$date)
merged_df <- full_df %>%
  left_join(df, by = c("date", "V")) %>%
  mutate(count = ifelse(is.na(count), 0, count))
voc = unique(merged_df$V)
observed_matrix = data.frame(v1 = merged_df[merged_df$V == 'D614G', 'count'],
                             v2 = merged_df[merged_df$V == 'Alpha', 'count'],
                             v3 = merged_df[merged_df$V == 'Delta', 'count'],
                             v3 = merged_df[merged_df$V == 'Omicron', 'count'])
rownames(observed_matrix) = full_dates

values = brewer.pal(n=4, name = "Spectral")
p1 = ggplot(data = merged_df, 
            aes(x = date, y = count, 
                group = V, color = V)) +
  geom_line() +
  scale_color_manual(name = '',
                     breaks = rev(levels(merged_df$V)),
                     values = rev(values)) +
  theme_bw() +
  scale_y_continuous(trans='log10', 
                     breaks = c(1,10,100,1000,10000),
                     labels = c(expression(10^0),expression(10^1),
                                expression(10^2), expression(10^3),
                                expression(10^4))) +
  scale_x_date(breaks = seq(as.Date('2020-01-01'), as.Date('2022-11-01'), by="6 months"),
               minor_breaks = seq(as.Date('2019-12-01'), as.Date('2022-11-01'), by ='1 month'),
               date_labels = "%y-%b", expand = c(0, 0)) +
  xlab('') + ylab('Cases') + 
  coord_cartesian(xlim = c(as.Date('2020-06-30'), as.Date('2022-05-01')),
                  ylim = c(1,2*10^4)) +
  theme(legend.position = c(0.18,0.9),
        legend.background = element_rect(color = NA, fill = NA),
        legend.key = element_blank(),
        legend.key.size = unit(0.12, units = 'cm'),
        legend.key.height= unit(0.06, 'cm'),
        legend.key.width = unit(0.2, 'cm'),
        legend.text = element_text(size=7))
p1

pdf(paste0("Output/voc_curve.pdf"), width = 2.5, height = 1.6)
print(p1)
dev.off()


data2 = merged_df %>%
  group_by(date) %>%
  mutate(p = count/sum(count)) %>%
  as.data.frame()
data2$group = factor(data2$V, levels = c('Omicron','Delta','Alpha','D614G'))

p2= ggplot() +
  geom_area(data = data2, 
            aes(x = date, y = p, fill = group),
            position = 'fill') +
  scale_fill_manual(name="", breaks = rev(levels(data2$group)),
                    values = alpha(rev(values), 1)) +
  coord_cartesian(xlim = c(as.Date('2020-06-30'), as.Date('2022-05-01')),
                  ylim = c(0,1)) +
  scale_x_date(breaks = seq(as.Date('2020-01-01'), as.Date('2024-11-01'), by="6 months"),
               minor_breaks = seq(as.Date('2019-12-01'), as.Date('2024-11-01'), by ='1 month'),
               date_labels = "%y-%b",
               expand = c(0, 0)) +
  scale_y_continuous(breaks = seq(0,1,0.5), expand = c(0, 0)) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = 'none',
        legend.background = element_rect(color = NA, fill = NA),
        legend.key = element_blank(),
        legend.key.size = unit(0.5, units = 'cm'),
        legend.key.height= unit(0.5, 'cm'),
        legend.text = element_text(size=0.5)) +
  xlab('') + ylab('Prevalence') 


pdf(paste0("Output/voc_prevalence.pdf"), width = 2.5, height = 1.2)
print(p2)
dev.off()

