rm(list = ls())
library(dplyr)
library(ggplot2)
library(scales)


df = read.csv('Covid19CasesGISAID.csv')
df$Var1 = as.Date(df$Var1)
df$Mutations = factor(df$Mutations, levels = unique(df$Mutations))

i = 2
for (i in c(1,3,4)) {
  c = cor(df$Freq[df$Mutations == unique(df$Mutations)[2]],
          df$Freq[df$Mutations == unique(df$Mutations)[i]])
  print(c)
}
# [1] -0.08236075
# [1] 0.7396961
# [1] 0.4516835
values = c(hue_pal()(3)[1], hue_pal()(3)[3], hue_pal()(3)[2], hue_pal()(4)[4])

anno = data.frame(x0 = as.Date(c('2020-8-10','2020-9-10','2020-10-29','2020-10-30','2021-9-1')),
                  y0 = c(128,8192,256,16384,2048),
                  text = c('Beta', 'Alpha', 'Gamma', 'Delta', 'Omicron'))
anno_arr = data.frame(x = as.Date(c('2020-8-10','2020-9-20','2020-10-29','2020-10-30','2021-9-1')),
                      y = c(165,6000,320,13000,2550),
                      xend = as.Date(c('2020-8-10','2020-9-20','2020-10-29','2020-10-30','2021-9-1')),
                      yend = c(500,3000,780,5200,9000),
                      text = c('Beta', 'Alpha', 'Gamma', 'Delta', 'Omicron'))


p <- ggplot(data = df, aes(x = Var1, y = Freq, colour = Mutations)) + 
  geom_point(size = 0.5) + 
  geom_smooth(aes(fill = Mutations, color = Mutations), 
              n = 300, size = 0.5) +
  scale_y_continuous(trans=scales::pseudo_log_trans(base=2),
                     breaks = c(1,10,100,1000,10000),
                     labels = c(expression(10^0),expression(10^1),
                                expression(10^2), expression(10^3),
                                expression(10^4))) +
  # annotation_logticks(linewidth = 0.1) +
  scale_x_date(breaks = seq(as.Date('2020-01-01'), as.Date('2022-11-01'), by="6 months"),
               minor_breaks = seq(as.Date('2019-12-01'), as.Date('2022-11-01'), by ='1 month'),
               date_labels = "%y-%b") +
  xlab('Date (2019-2022)') + 
  ylab('Cases') + 
  theme_bw() +
  scale_color_manual(name="",
                     labels=c("A", "B", '8782 C>T', "28144 T>C"),
                     values = alpha(values, 0.6)) +
  scale_fill_manual(name="",
                    labels=c("A", "B", '8782 C>T', "28144 T>C"),
                    values = alpha(values, 0.1)) +
  coord_cartesian(xlim = c(as.Date('2019-12-31'), as.Date('2021-10-31')),
                  ylim = c(1,2*10^4)) +
  geom_text(data = anno, aes(x = x0,  y = y0, label = text),
            inherit.aes = F, size = 2.8, alpha = .7) +
  geom_segment(data = anno_arr, aes(x = x, y = y, xend = xend, yend = yend),
               arrow = arrow(length = unit(0.1, "cm")), alpha = .7,
               inherit.aes = F) +
  theme(legend.key.size = unit(0.3, 'cm'),
        legend.position = 'top',
        legend.margin = margin(t = 0, r = 0, b = -2, l = 0)
  ) 
p


pdf(paste0("Output/gr.pdf"), width = 3, height = 2.5)
print(p)
dev.off()


