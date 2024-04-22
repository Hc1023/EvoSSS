rm(list = ls())
library(dplyr)
library(ggplot2)
library(scales)

get.all.bootstraps <- function(df, R = 10){
  
  sample.size = nrow(df)
  
  # Original samples, ensure complete dates
  r = 1
  q = qplot(Var1, Freq, data = df) + 
    stat_smooth(n=nrow(df)) +
    scale_y_continuous(trans=scales::pseudo_log_trans(base = 2))
  gq = ggplot_build(q)$data[[2]]
  gq$x = as.Date(gq$x, origin = "1970-01-01")
  gq$y = gq$y
  all.bootstraps = data.frame(V1 = gq$y)
  rownames(all.bootstraps) = gq$x
  
  repeat {
    r = r + 1
    # Bootstrapping is any test or metric that uses random sampling with replacement
    df.sample = df[sample(nrow(df), sample.size, replace = T),]
    q = qplot(Var1, Freq, data = df.sample) + 
      stat_smooth(n=nrow(df.sample)) +
      scale_y_continuous(trans=scales::pseudo_log_trans(base = 2))
    gq = ggplot_build(q)$data[[2]]
    gq$x = as.character(as.Date(gq$x, origin = "1970-01-01"))
    all.bootstraps[gq$x, r] = gq$y
    
    if(r == R+1){
      break
    }
  }
  return(all.bootstraps)
}


extract.CI.bootstrap <- function(all.bootstraps){
  CI.all = data.frame(
    Var1 = as.Date(rownames(all.bootstraps)),
    median = apply(all.bootstraps, 1, function(x){2^quantile(x, 0.5, na.rm = T)-1}),
    lower = apply(all.bootstraps, 1, function(x){2^quantile(x, 0.025, na.rm = T)-1}),
    upper = apply(all.bootstraps, 1, function(x){2^quantile(x, 0.975, na.rm = T)-1})
  )
  return(CI.all)
}

getGR <- function(one.bootstrap, x.values){
  idx <- !is.na(one.bootstrap)
  one.bootstrap <- one.bootstrap[idx]
  x.values <- x.values[idx]
  GR.one = data.frame()
  for (i in 2:(length(x.values)-1)) {
    # exponential growth (exp() <-> log)
    # log(2^y) = y*log(2)
    dy = (one.bootstrap[i+1]-one.bootstrap[i-1])*log(2)
    dt = x.values[i+1]-x.values[i-1]
    GR.one = rbind(GR.one, data.frame(x.values = x.values[i], GR = dy/dt))
  }
  return(GR.one)
}


extract.GR.bootstrap <- function(all.bootstraps){
  GR.all <- data.frame(matrix(nrow = nrow(all.bootstraps), ncol = 0))
  rownames(GR.all) <- rownames(all.bootstraps)
  for (i in 1:ncol(all.bootstraps)) {
    GR.one = getGR(unlist(all.bootstraps[,i]), x.values = 1:nrow(all.bootstraps))
    GR.all[GR.one$x.values,i] = GR.one$GR
  }
  GR.CI.all = data.frame(
    Var1 = as.Date(rownames(GR.all)),
    GR.median = apply(GR.all, 1, function(x){2^quantile(x, 0.5, na.rm = T)-1}),
    GR.lower = apply(GR.all, 1, function(x){2^quantile(x, 0.025, na.rm = T)-1}),
    GR.upper = apply(GR.all, 1, function(x){2^quantile(x, 0.975, na.rm = T)-1})
  )
  return(GR.CI.all)
}

if(F){
  df.plot.all = data.frame()
  l = 'Lineage A'
  for(l in unique(dfplot$Mutations)){
    df = dfplot[dfplot$Mutations == l,c(1:2)]
    df.bootstrap <- get.all.bootstraps(df, R = 1000)
    df.CI.all <- extract.CI.bootstrap(df.bootstrap)
    
    print(l)
    
    df.GR.all <- extract.GR.bootstrap(df.bootstrap)
    df.plot = merge(df, df.CI.all, by = "Var1", all = T)
    df.plot = merge(df.plot, df.GR.all, by = 'Var1', all = T)
    df.plot$Mutations = l
    df.plot.all = rbind(df.plot.all, df.plot)
  }
  
  df.plot.all$Mutations <- factor(df.plot.all$Mutations,
                                  levels = c('Lineage A', 'Lineage B',
                                             '8782T', '28144C'))
}

values = c(hue_pal()(3)[1], hue_pal()(3)[3], hue_pal()(3)[2], hue_pal()(4)[4])

## Growth curve (Bootstrapping)
l = unique(df.plot.all$Mutations)
df.plot1 = df.plot.all[df.plot.all$Mutations == l[1],]
df.plot1[2,6:8]
# GR.median   GR.lower   GR.upper
# 2 0.03878277 0.03292318 0.04576787
df.plot2 = df.plot.all[df.plot.all$Mutations == l[2],]
df.plot2[2,6:8]
# GR.median   GR.lower   GR.upper
# 1017 0.05529341 0.04778657 0.06182397
# > 0.03878277/0.05529341
# [1] 0.7013995
# > 0.04576787/0.04778657
# [1] 0.9577559
# > 0.03292318/0.06182397
# [1] 0.532531
df.plot1[2,6:8]
idx = df.plot1[,7]>0 & df.plot1[,1]<as.Date('2020-03-01')
idx[1] = F
df.plot1[idx,6]/df.plot2[idx,6]

df.polygon.all = data.frame()
for(l in unique(df.plot.all$Mutations)){
  df.plot = df.plot.all[df.plot.all$Mutations == l,]
  idx <- !(is.na(df.plot$median))
  df.polygon = data.frame(x = c(df.plot$Var1[idx], rev(df.plot$Var1[idx])), 
                          y = c(df.plot$upper[idx], rev(df.plot$lower[idx])))
  
  df.polygon$Mutations = l
  df.polygon.all = rbind(df.polygon.all, df.polygon)
  
}

df.polygon.all$Mutations <- factor(df.polygon.all$Mutations,
                                   levels = c('Lineage A', 'Lineage B',
                                              '8782T', '28144C'))



p <- ggplot(df.plot.all, aes(x = Var1, y = Freq, colour = Mutations)) + 
  geom_point(size = 0.5) +
  scale_y_continuous(trans=scales::pseudo_log_trans(base = 2),
                     breaks = c(0, 2^seq(0,14,2))) +
  annotation_logticks(size = 0.1) +
  scale_x_date(breaks = c(as.Date('2019-12-01'),
                          seq(as.Date('2019-05-01'), as.Date('2022-11-01'), by="6 months")),
               minor_breaks = seq(as.Date('2019-12-01'), as.Date('2022-11-01'), by ='1 month'),
               date_labels = "%y %b") +
  # scale_x_date(date_breaks = "3 month", minor_breaks = "1 month", date_labels = "%y %b") +
  xlab('Date (2019-2022)') + ylab('Sequences') + 
  theme_bw() +
  scale_color_manual(name="Clades",
                     labels=c("Lineage A", "Lineage B", 
                              'T/T', "C/C"),
                     values = alpha(values, 0.6)) +
  geom_polygon(data = df.polygon.all, 
               aes(x = x, y = y, fill = Mutations),
               color = NA) +
  scale_fill_manual(name="Clades",
                    labels=c("Lineage A", "Lineage B", 
                             'T/T', "C/C"),
                    values = alpha(values, 0.3)) +
  geom_line(aes(x = Var1, y = median), lwd=0.5) +
  coord_cartesian(ylim = c(0, max(df.plot.all$Freq))) +
  theme(axis.title=element_text(size=20),
        axis.text=element_text(size=14),
        #axis.text.x = element_text(angle = 60, hjust = 1),
        legend.key.size = unit(0.8, 'cm'),
        legend.title = element_text(size=18),
        legend.text = element_text(size=16),
        legend.position = 'top',
        plot.margin = margin(0.5,1,0,0.5, "cm") # trbl
  )

p

pdf(paste0("Output/gr_bootstrap", ".pdf"), width = 8, height = 5)
print(p)
dev.off()


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
  scale_x_date(breaks = c(as.Date('2019-12-01'),
                          seq(as.Date('2019-05-01'), as.Date('2022-11-01'), by="6 months")),
               minor_breaks = seq(as.Date('2019-12-01'), as.Date('2022-11-01'), by ='1 month'),
               date_labels = "%y %b") +
  xlab('') + ylab('') + theme_bw() +
  scale_color_manual(values = alpha(values, 0.6)) +
  geom_polygon(data = df.polygon.all, aes(x = x, y = y, fill = Mutations),
               color = NA) +
  scale_fill_manual(values = alpha(values, 0.2)) +
  geom_line() + 
  geom_hline(yintercept = 0, linetype="dashed", size=0.3) +
  # coord_cartesian(xlim = c(min(dfplot$Var1), max(dfplot$Var1))) +
  geom_text(data = anno, aes(x = x0,  y = y0, label = text),
            inherit.aes = F, size = 5, alpha = .7) +
  geom_segment(data = anno_arr, aes(x = x, y = y, xend = xend, yend = yend),
               arrow = arrow(length = unit(0.2, "cm")), alpha = .7,
               inherit.aes = F) +
  theme(axis.title=element_text(size=20),
        axis.text=element_text(size=14),
        #axis.text.x = element_text(angle = 60, hjust = 1),
        legend.key.size = unit(0.8, 'cm'),
        legend.title = element_text(size=18),
        legend.text = element_text(size=16),
        legend.position = 'none',
        plot.margin = margin(0.5,1,0,0.5, "cm") # trbl
  ) 
p
pdf(paste0("Output/gr_rate_bootstrap", ".pdf"), width = 6.3, height = 3.5)
print(p)
dev.off()

