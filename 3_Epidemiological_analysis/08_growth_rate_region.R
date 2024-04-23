rm(list = ls())
library(dplyr)
library(ggplot2)
library(scales)

load(file = 'GRrate_region.RData')
regions = unique(df.plot.all.region$Region)

# GAM fitting
# `geom_smooth()` using method = 'gam' and formula 'y ~ s(x, bs = "cs")'
plist = list()
for (i in 1:length(regions)) {
  values = c(hue_pal()(3)[1], hue_pal()(3)[3], hue_pal()(3)[2], hue_pal()(4)[4])
  df.plot.all = df.plot.all.region[df.plot.all.region$Region == regions[i],]
  
  df.polygon.all = data.frame()
  for(l in unique(df.plot.all$Mutations)){
    df.plot = df.plot.all[df.plot.all$Mutations == l,]
    idx <- !(is.na(df.plot$GR.median))
    df.polygon = data.frame(x = c(df.plot$Var1[idx], rev(df.plot$Var1[idx])), 
                            y = c(df.plot$GR.upper[idx], rev(df.plot$GR.lower[idx])))
    
    df.polygon$Mutations = l
    df.polygon.all = rbind(df.polygon.all, df.polygon)
  }
  df.plot.all$Mutations <- factor(df.plot.all$Mutations,
                                     levels = c('A', 'B',
                                                'T', 'C'))
  df.polygon.all$Mutations <- factor(df.polygon.all$Mutations,
                                     levels = c('A', 'B',
                                                'T', 'C'))
  
  p = ggplot(df.plot.all, aes(x = Var1, y = GR.median, colour = Mutations)) + 
    scale_x_date(breaks = c(as.Date('2019-12-01'),
                            seq(as.Date('2019-05-01'), as.Date('2022-11-01'), by="6 months")),
                 minor_breaks = seq(as.Date('2019-12-01'), as.Date('2022-11-01'), by ='1 month'),
                 date_labels = "%y %b") +
    xlab('') + ylab('Growth rate') + theme_bw() +
    scale_color_manual(values = alpha(values, 0.6)) +
    geom_polygon(data = df.polygon.all, aes(x = x, y = y, fill = Mutations),
                 color = NA) +
    scale_fill_manual(values = alpha(values, 0.2)) +
    geom_line() + 
    geom_hline(yintercept = 0, linetype="dashed", size=0.3) +
    coord_cartesian(xlim = c(min(df.plot.all$Var1), max(df.plot.all$Var1))) +
    theme(axis.title=element_text(size=18),
          axis.text=element_text(size=14),
          #axis.text.x = element_text(angle = 60, hjust = 1),
          legend.key.size = unit(0.8, 'cm'),
          legend.title = element_text(size=18),
          legend.text = element_text(size=16),
          legend.position = 'none',
          axis.text.x = element_text(angle = 60, vjust=1, hjust =1),
          plot.title = element_text(size = 20),
          plot.margin = margin(0.5,1,0,0.5, "cm") # trbl
    ) +
    ggtitle(regions[i])
  
  plist[[i]] = p
}


pdf(paste0("Output/gr_rate_regions.pdf"), width = 5, height = 3.5)
for (i in 1:length(plist)) {
  print(plist[[i]])
}
dev.off()
