rm(list = ls())
library(ggplot2)
library(scales)

load('dfplot_region.RData')
point_size = 0.1

p1 = ggplot(dfplot_A, aes(x = Var1, y = Freq, colour = Region)) + 
  scale_y_continuous(trans=scales::pseudo_log_trans(base = 2),
                     breaks = c(0, 2^seq(0,14,2)),
                     labels = c(0, 1, expression(2^{2}),
                                expression(2^{4}),expression(2^{6}),
                                expression(2^{8}),expression(2^{10}),
                                expression(2^{12}),expression(2^{14}))) +
  scale_x_date(breaks = c(as.Date('2019-12-01'),
                          seq(as.Date('2019-05-01'), as.Date('2022-11-01'), by="6 months")),
               minor_breaks = seq(as.Date('2019-12-01'), as.Date('2022-11-01'), by ='1 month'),
               date_labels = "%y %b") +
  xlab('Date (2019-2022)') + ylab('Lineage A') + theme_bw() +
  geom_point(size = point_size) + geom_smooth(aes(fill = Region), lwd = 0.5) +
  scale_color_manual(name = 'Regions', 
                     values = alpha(hue_pal()(6), 0.6)) +
  scale_fill_manual(name = 'Regions', 
                    values = alpha(hue_pal()(6), 0.1)) +
  coord_cartesian(ylim = c(0, max(dfplot_A$Freq))) +
  theme(axis.title=element_text(size=16),
        axis.text=element_text(size=13),
        #axis.text.x = element_text(angle = 60, hjust = 1),
        legend.key.size = unit(0.8, 'cm'),
        legend.title = element_text(size=18),
        legend.text = element_text(size=16),
        legend.position = 'top',
        plot.margin = margin(0.5,1,0,0.5, "cm") # trbl
  ) + guides(colour = guide_legend(nrow = 1))
p1

p2 = ggplot(dfplot_B, aes(x = Var1, y = Freq, colour = Region)) + 
  scale_y_continuous(trans=scales::pseudo_log_trans(base = 2),
                     breaks = c(0, 2^seq(0,14,2)),
                     labels = c(0, 1, expression(2^{2}),
                                expression(2^{4}),expression(2^{6}),
                                expression(2^{8}),expression(2^{10}),
                                expression(2^{12}),expression(2^{14}))) +
  scale_x_date(breaks = c(as.Date('2019-12-01'),
                          seq(as.Date('2019-05-01'), as.Date('2022-11-01'), by="6 months")),
               minor_breaks = seq(as.Date('2019-12-01'), as.Date('2022-11-01'), by ='1 month'),
               date_labels = "%y %b") +
  xlab('Date (2019-2022)') + ylab('Lineage B') + theme_bw() +
  geom_point(size = point_size) + geom_smooth(aes(fill = Region), lwd = 0.5) +
  scale_color_manual(name = 'Regions', 
                     values = alpha(hue_pal()(6), 0.6)) +
  scale_fill_manual(name = 'Regions', 
                    values = alpha(hue_pal()(6), 0.1)) +
  coord_cartesian(ylim = c(0, max(dfplot_B$Freq))) +
  theme(axis.title=element_text(size=16),
        axis.text=element_text(size=13),
        #axis.text.x = element_text(angle = 60, hjust = 1),
        legend.key.size = unit(0.8, 'cm'),
        legend.title = element_text(size=18),
        legend.text = element_text(size=16),
        legend.position = 'top',
        plot.margin = margin(0.5,1,0,0.5, "cm") # trbl
  ) + guides(colour = guide_legend(nrow = 1))
p2


p3 = ggplot(dfplot_T, aes(x = Var1, y = Freq, colour = Region)) + 
  scale_y_continuous(trans=scales::pseudo_log_trans(base = 2),
                     breaks = c(0, 2^seq(0,14,2)),
                     labels = c(0, 1, expression(2^{2}),
                                expression(2^{4}),expression(2^{6}),
                                expression(2^{8}),expression(2^{10}),
                                expression(2^{12}),expression(2^{14}))) +
  scale_x_date(breaks = c(as.Date('2019-12-01'),
                          seq(as.Date('2019-05-01'), as.Date('2022-11-01'), by="6 months")),
               minor_breaks = seq(as.Date('2019-12-01'), as.Date('2022-11-01'), by ='1 month'),
               date_labels = "%y %b") +
  xlab('Date (2019-2022)') + ylab('8782 C>T') + theme_bw() +
  geom_point(size = point_size) + geom_smooth(aes(fill = Region), lwd = 0.5) +
  scale_color_manual(name = 'Regions', 
                     values = alpha(hue_pal()(6), 0.6)) +
  scale_fill_manual(name = 'Regions', 
                    values = alpha(hue_pal()(6), 0.1)) +
  coord_cartesian(ylim = c(0, max(dfplot_T$Freq))) +
  theme(axis.title=element_text(size=16),
        axis.text=element_text(size=13),
        #axis.text.x = element_text(angle = 60, hjust = 1),
        legend.key.size = unit(0.8, 'cm'),
        legend.title = element_text(size=18),
        legend.text = element_text(size=16),
        legend.position = 'top',
        plot.margin = margin(0.5,1,0,0.5, "cm") # trbl
  ) + guides(colour = guide_legend(nrow = 1))

p3


p4 = ggplot(dfplot_C, aes(x = Var1, y = Freq, colour = Region)) + 
  scale_y_continuous(trans=scales::pseudo_log_trans(base = 2),
                     breaks = c(0, 2^seq(0,14,2)),
                     labels = c(0, 1, expression(2^{2}),
                                expression(2^{4}),expression(2^{6}),
                                expression(2^{8}),expression(2^{10}),
                                expression(2^{12}),expression(2^{14}))) +
  scale_x_date(breaks = c(as.Date('2019-12-01'),
                          seq(as.Date('2019-05-01'), as.Date('2022-11-01'), by="6 months")),
               minor_breaks = seq(as.Date('2019-12-01'), as.Date('2022-11-01'), by ='1 month'),
               date_labels = "%y %b") +
  xlab('Date (2019-2022)') + ylab('28144 T>C') + theme_bw() +
  geom_point(size = point_size) + geom_smooth(aes(fill = Region), lwd = 0.5) +
  scale_color_manual(name = 'Regions', 
                     values = alpha(hue_pal()(6), 0.6)) +
  scale_fill_manual(name = 'Regions', 
                    values = alpha(hue_pal()(6), 0.1)) +
  coord_cartesian(ylim = c(0, max(dfplot_C$Freq))) +
  theme(axis.title=element_text(size=16),
        axis.text=element_text(size=13),
        #axis.text.x = element_text(angle = 60, hjust = 1),
        legend.key.size = unit(0.8, 'cm'),
        legend.title = element_text(size=18),
        legend.text = element_text(size=16),
        legend.position = 'top',
        plot.margin = margin(0.5,1,0,0.5, "cm") # trbl
  ) + guides(colour = guide_legend(nrow = 1))
p4


pdf(paste0("Output/gr_region", ".pdf"), width = 5.5, height = 3.8)
print(p1)
print(p2)
print(p3)
print(p4)
dev.off()


