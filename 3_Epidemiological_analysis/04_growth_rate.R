library(dplyr)
library(ggplot2)
library(scales)

load("GRrate.RData")

## Growth rate curve  (Bootstrapping)

# GAM fitting
# `geom_smooth()` using method = 'gam' and formula 'y ~ s(x, bs = "cs")'

values = c(hue_pal()(3)[1], hue_pal()(3)[3], hue_pal()(3)[2], hue_pal()(4)[4])

df.polygon.all = data.frame()
for(l in unique(df.plot.all$Mutations)){
  df.plot = df.plot.all[df.plot.all$Mutations == l,]
  idx <- !(is.na(df.plot$GR.median))
  df.polygon = data.frame(x = c(df.plot$Var1[idx], rev(df.plot$Var1[idx])), 
                          y = c(df.plot$GR.upper[idx], rev(df.plot$GR.lower[idx])))
  
  df.polygon$Mutations = l
  df.polygon.all = rbind(df.polygon.all, df.polygon)
}
df.polygon.all$Mutations <- factor(df.polygon.all$Mutations,
                                   levels = c('Lineage A', 'Lineage B',
                                              '8782T', '28144C'))

anno = data.frame(x0 = as.Date(c('2020-7-28','2020-8-29','2020-11-05','2020-10-30','2021-9-1')),
                  y0 = c(0.022,0.031,-0.01,0.04,0.025),
                  text = c('Beta', 'Alpha', 'Gamma', 'Delta', 'Omicron'))
anno_arr = data.frame(x = as.Date(c('2020-8-10','2020-9-20','2020-10-29','2020-10-30','2021-9-1')),
                      y = c(0.018,0.027,-0.006,0.036,0.021),
                      xend = as.Date(c('2020-8-10','2020-9-20','2020-10-29','2020-10-30','2021-9-1')),
                      yend = c(0.005,0.012,0.01,0.015,0.008),
                      text = c('Beta', 'Alpha', 'Gamma', 'Delta', 'Omicron'))

p = ggplot(df.plot.all, aes(x = Var1, y = GR.median, colour = Mutations)) + 
  scale_x_date(breaks = seq(as.Date('2020-01-01'), as.Date('2022-11-01'), by="6 months"),
               minor_breaks = seq(as.Date('2019-12-01'), as.Date('2022-11-01'), by ='1 month'),
               date_labels = "%y-%b") +
  xlab('Date (2019-2022)') + 
  ylab('Growth rate') + 
  theme_bw() +
  geom_polygon(data = df.polygon.all, aes(x = x, y = y, fill = Mutations),
               color = NA) +
  scale_color_manual(name="Variant",
                     labels=c("A", "B", '8782 C>T', "28144 T>C"),
                     values = alpha(values, 0.7)) +
  scale_fill_manual(name="Variant",
                    labels=c("A", "B", '8782 C>T', "28144 T>C"),
                    values = alpha(values, 0.3)) +
  geom_line() + 
  geom_hline(yintercept = 0, linetype="dashed", size=0.3) +
  coord_cartesian(xlim = c(as.Date('2019-12-31'), as.Date('2021-10-31'))) +  
  geom_text(data = anno, aes(x = x0,  y = y0, label = text),
            inherit.aes = F, size = 2.8, alpha = .7) +
  geom_segment(data = anno_arr, aes(x = x, y = y, xend = xend, yend = yend),
               arrow = arrow(length = unit(0.1, "cm")), alpha = .7,
               inherit.aes = F) +
  theme(legend.key.size = unit(0.3, 'cm'),
        legend.position = 'right',
        legend.key.spacing.x = unit(0.1,'cm'),
        legend.margin = margin(t = 0, r = 0, b = -2, l = 0)
  ) 
p
pdf(paste0("Output/gr_rate_bootstrap.pdf"), width = 3.7, height = 1.8)
print(p)
dev.off()

