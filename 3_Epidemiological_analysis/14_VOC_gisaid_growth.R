rm(list = ls())
library(dplyr)
library(ggplot2)
library(scales)
library(RColorBrewer)

result = read.csv('VOC_gisaid.csv')
result$date = as.Date(result$date)
df = result

VOC = c("Delta","Alpha","D614G","B 19A","A")
names(VOC) = c("GK","GRY","G","S","L")
df = df[df$GISAID_clade %in% names(VOC),]
df$V = NA
for (i in 1:length(VOC)) {
  df$V[df$GISAID_clade == names(VOC)[i]] = VOC[i]
}
df$V = factor(df$V, levels = VOC)
dfall = df %>% group_by(date) %>%
  summarise(yall = sum(count))

brews = brewer.pal(n=4,name = "Spectral")
values = c(brews[2:4], hue_pal()(3)[c(3,1)])
p = ggplot() +
  geom_point(data = df, 
            aes(x = date, y = count, 
                group = V, color= V),
            size = 0.5) +
  geom_line(data = dfall,
            aes(x = date, y = yall),
            linewidth = 0.5, alpha = 0.3) +
  scale_y_continuous(trans='log10', 
                     breaks = c(1,10,100,1000,10000),
                     labels = c(expression(10^0),expression(10^1),
                                expression(10^2), expression(10^3),
                                expression(10^4))) +
  scale_fill_manual(values = alpha(values, 0.3)) +
  scale_color_manual(values = alpha(values, 0.7)) +
  labs(x = "", y = "") +
  theme_bw() +
  theme(legend.position = "right",
        plot.title = element_text(hjust = 0.5)) +
  scale_x_date(breaks = seq(as.Date('2020-01-01'), as.Date('2022-11-01'), by="6 months"),
               minor_breaks = seq(as.Date('2019-12-01'), as.Date('2022-11-01'), by ='1 month'),
               date_labels = "%y-%b") +
  xlab('') + ylab('Cases') + 
  coord_cartesian(xlim = c(as.Date('2019-12-31'), as.Date('2021-10-31')),
                  ylim = c(1,2*10^4)) 
pdf(paste0("Output/VOC_gisaid_growth.pdf"), width = 3.5, height = 1.8)
print(p)
dev.off()
