rm(list = ls())
library(ggplot2)
library(scales)

load('dfplot_region.RData')

point_size = 0.1
values = hue_pal()(6)

getplot = function(dfplot){
  p = ggplot(dfplot, aes(x = Var1, y = Freq, colour = Region)) + 
    scale_y_continuous(trans=scales::pseudo_log_trans(base = 2),
                       breaks = c(1,10,100,1000,10000),
                       labels = c(expression(10^0),expression(10^1),
                                  expression(10^2), expression(10^3),
                                  expression(10^4))) +
    scale_x_date(breaks = seq(as.Date('2020-01-01'), as.Date('2022-11-01'), by="6 months"),
                 minor_breaks = seq(as.Date('2019-12-01'), as.Date('2022-11-01'), by ='1 month'),
                 date_labels = "%y-%b") +
    xlab('') + ylab('') + theme_bw() +
    geom_point(size = point_size) + geom_smooth(aes(fill = Region), lwd = 0.5) +
    scale_color_manual(name = 'Regions', 
                       values = alpha(values, 0.3)) +
    scale_fill_manual(name = 'Regions', 
                      values = alpha(values, 0.1)) +
    coord_cartesian(xlim = c(as.Date('2019-12-31'), as.Date('2021-10-31')),
                    ylim = c(0, max(dfplot$Freq))) +
    theme(legend.position = 'none')
  return(p)
}

p1 = getplot(dfplot_A)
p2 = getplot(dfplot_B)
p3 = getplot(dfplot_C)
p4 = getplot(dfplot_T)

pdf(paste0("Output/growth_region.pdf"), width = 2.5, height = 1.6)
print(p1)
print(p2)
print(p3)
print(p4)
dev.off()

dfplot = dfplot_A
p = ggplot(dfplot, aes(x = Var1, y = Freq, colour = Region)) + 
  scale_y_continuous(trans=scales::pseudo_log_trans(base = 2),
                     breaks = c(1,10,100,1000,10000),
                     labels = c(expression(10^0),expression(10^1),
                                expression(10^2), expression(10^3),
                                expression(10^4))) +
  scale_x_date(breaks = seq(as.Date('2020-01-01'), as.Date('2022-11-01'), by="6 months"),
               minor_breaks = seq(as.Date('2019-12-01'), as.Date('2022-11-01'), by ='1 month'),
               date_labels = "%y-%b") +
  xlab('') + ylab('') + theme_bw() +
  geom_point(size = point_size) +
  geom_smooth(aes(fill = Region), lwd = 0.5) +
  scale_color_manual(name = 'Regions', 
                     values = alpha(values, 0.3)) +
  scale_fill_manual(name = 'Regions', 
                    values = alpha(values, 0.1)) +
  coord_cartesian(xlim = c(as.Date('2019-12-31'), as.Date('2021-10-31')),
                  ylim = c(0, 100)) +
  theme(legend.position = 'none',
        panel.grid.minor = element_blank(),
        plot.background = element_blank(),
        panel.background = element_blank())

p
pdf(paste0("Output/growth_A.pdf"), width = 2.5, height = 1.6)
print(p)
dev.off()

