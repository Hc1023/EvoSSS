rm(list = ls())
library(dplyr)
library(ggplot2)
library(scales)

dfplot<- read.csv('Covid19CasesGISAID.csv')
dfplot$Var1 = as.Date(dfplot$Var1)
dfplot$Mutations = factor(dfplot$Mutations, levels = c("Lineage A","Lineage B","8782T","28144C"))
values = c(hue_pal()(3)[1], hue_pal()(3)[3], hue_pal()(3)[2], hue_pal()(4)[4])

anno = data.frame(x0 = as.Date(c('2020-8-10','2020-9-10','2020-10-29','2020-10-30','2021-9-1')),
                  y0 = c(128,8192,256,16384,2048),
                  text = c('Beta', 'Alpha', 'Gamma', 'Delta', 'Omicron'))
anno_arr = data.frame(x = as.Date(c('2020-8-10','2020-9-20','2020-10-29','2020-10-30','2021-9-1')),
                      y = c(165,6000,320,13000,2550),
                      xend = as.Date(c('2020-8-10','2020-9-20','2020-10-29','2020-10-30','2021-9-1')),
                      yend = c(500,3000,780,5200,9000),
                      text = c('Beta', 'Alpha', 'Gamma', 'Delta', 'Omicron'))

p <- ggplot(data = dfplot, aes(x = Var1, y = Freq, colour = Mutations)) + 
  geom_point(size = 0.5) + 
  geom_smooth(aes(fill = Mutations, color = Mutations), 
              n = 300, size = 0.5) +
  scale_y_continuous(trans=scales::pseudo_log_trans(base=2),
                     breaks = c(0, 2^seq(0,14,2))) +
  annotation_logticks(linewidth = 0.1, alpha = 0.5) +
  scale_x_date(breaks = c(as.Date('2019-12-01'),
                          seq(as.Date('2019-05-01'), as.Date('2022-11-01'), by="6 months")),
               minor_breaks = seq(as.Date('2019-12-01'), as.Date('2022-11-01'), by ='1 month'),
               date_labels = "%y %b") +
  xlab('Date (2019-2021)') +
  ylab('') + theme_bw() +
  scale_color_manual(name="Clades",
                     labels=c("Lineage A", "Lineage B", '8782 C>T', "28144 T>C"),
                     values = alpha(values, 0.6)) +
  scale_fill_manual(name="Clades",
                    labels=c("Lineage A", "Lineage B", '8782 C>T', "28144 T>C"),
                    values = alpha(values, 0.1)) +
  coord_cartesian(ylim = c(0, max(dfplot$Freq)), 
                  xlim = c(min(dfplot$Var1), as.Date('2021-11-01'))) +
  geom_text(data = anno, aes(x = x0,  y = y0, label = text),
            inherit.aes = F, size = 5, alpha = .7) +
  geom_segment(data = anno_arr, aes(x = x, y = y, xend = xend, yend = yend),
               arrow = arrow(length = unit(0.2, "cm")), alpha = .7,
               inherit.aes = F) +
  theme(axis.title=element_text(size=16),
        axis.text=element_text(size=14),
        legend.key.size = unit(0.8, 'cm'),
        legend.title = element_text(size=14),
        legend.text = element_text(size=12),
        legend.position = 'none',
        plot.margin = margin(0.5,1,0,0.5, "cm") # trbl
  ) 
p


pdf(paste0("Output/gr.pdf"), width = 5.8, height = 4)
print(p)
dev.off()
