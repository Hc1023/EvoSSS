rm(list = ls())
library(dplyr)
library(ggplot2)
library(scales)
library(mgcv)
library(rstan)

df = read.csv('Covid19CasesGISAID.csv')
df$Var1 = as.Date(df$Var1)
df = df[df$Var1 < as.Date('2021-12-01'),]

getplot = function(df, fit){
  ndays = nrow(df)
  time = 1:ndays
  cases =  df$Freq
  # Fit a GAM with P-splines
  gam_model <- gam(cases ~ s(time, bs = "ps"), 
                   family = poisson(link = "log")
  )
  # Extracting model matrix
  X <- model.matrix(gam_model)
  posterior_samples <- extract(fit)
  simu_Onset <- matrix(NA, nrow = 15000, ncol = ndays)
  for (i in 1:15000) {
    pars = posterior_samples$beta[i,]
    simu_Onset[i, ] <- exp(X %*% pars)  
  }
  ci_lower <- apply(simu_Onset, 2, quantile, probs = 0.025, na.rm = T)
  ci_upper <- apply(simu_Onset, 2, quantile, probs = 0.975, na.rm = T)
  start_date <- as.Date("2019-12-24") 
  date_vector <- seq.Date(start_date, by = "day", length.out = ndays)
  
  plot_data <- data.frame(
    Day = 1:ndays,
    date_vector = date_vector,
    Observed = cases,
    Fitted = colMeans(simu_Onset),
    LowerCI = ci_lower,
    UpperCI = ci_upper
  )
  return(plot_data)
}


load(file = 'grBspline.Rdata')
df1 = df[df$Mutations == 'Lineage A',]
df2 = df[df$Mutations == 'Lineage B',]
plot_data1 = getplot(df = df1, fit = fit1)
plot_data2 = getplot(df = df2, fit = fit2)
plot_data1$Mutations = 'Lineage A'
plot_data2$Mutations = 'Lineage B'
plot_data = rbind(plot_data1, plot_data2)

values = c(hue_pal()(3)[1], hue_pal()(3)[3])

anno = data.frame(x0 = as.Date(c('2020-8-10','2020-9-10','2020-10-29','2020-10-30','2021-9-1')),
                  y0 = c(128,8192,256,16384,2048),
                  text = c('Beta', 'Alpha', 'Gamma', 'Delta', 'Omicron'))
anno_arr = data.frame(x = as.Date(c('2020-8-10','2020-9-20','2020-10-29','2020-10-30','2021-9-1')),
                      y = c(165,6000,320,13000,2550),
                      xend = as.Date(c('2020-8-10','2020-9-20','2020-10-29','2020-10-30','2021-9-1')),
                      yend = c(500,3000,780,5200,9000),
                      text = c('Beta', 'Alpha', 'Gamma', 'Delta', 'Omicron'))

p = ggplot(plot_data, aes(x = date_vector, group = Mutations)) +
  geom_ribbon(aes(ymin = LowerCI, ymax = UpperCI, fill = Mutations),  alpha = 0.6) +  # Confidence interval
  scale_y_continuous(trans=scales::pseudo_log_trans(base=2),
                     breaks = c(0,2^seq(2,14,3)),
                     labels = c(0, expression('10'^'2'),expression('10'^'5'),
                                expression('10'^'8'),expression('10'^'11'),
                                expression('10'^'14'))) +
  annotation_logticks(linewidth = 0.1, alpha = 0.5) +
  geom_point(aes(y = Observed, color = Mutations), size = 1, alpha = 0.5, shape = 16) + 
  geom_line(aes(y = Observed, color = Mutations), size = 0.2, alpha = 0.5) +  # Observed data
  geom_line(aes(y = Fitted, color = Mutations), size = 0.8, alpha = 0.7) +  # Fitted line
  scale_x_date(breaks = seq(as.Date('2020-01-01'), as.Date('2022-11-01'), by="4 months"),
               minor_breaks = seq(as.Date('2019-12-01'), as.Date('2022-11-01'), by ='1 month'),
               date_labels = "%y-%b") +
  xlab('Date (2019-2021)') +
  ylab('') + theme_bw() +
  scale_color_manual(name="",
                     labels=c("Lineage A", "Lineage B", '8782 C>T', "28144 T>C"),
                     values = alpha(values, 0.6)) +
  scale_fill_manual(name="",
                    labels=c("Lineage A", "Lineage B", '8782 C>T', "28144 T>C"),
                    values = alpha(values, 0.1)) +
  coord_cartesian(ylim = c(0, max(plot_data$Observed)), 
                  xlim = c(min(plot_data$date_vector), as.Date('2021-10-05'))) +
  geom_text(data = anno, aes(x = x0,  y = y0, label = text),
            inherit.aes = F, size = 3, alpha = .7) +
  geom_segment(data = anno_arr, aes(x = x, y = y, xend = xend, yend = yend),
               arrow = arrow(length = unit(0.2, "cm")), alpha = .7,
               inherit.aes = F) +
  theme(legend.background = element_blank(),
        # legend.box.background = element_rect(fill = alpha("white", 0.5), 
        #                                      color = alpha('black', 0.5), 
        #                                      linewidth = 0.2),
        legend.position = c(0.16,0.86),
        legend.title = element_blank())


p
pdf(paste0("Output/gr.pdf"), width = 4.5, height = 3)
print(p)
dev.off()
