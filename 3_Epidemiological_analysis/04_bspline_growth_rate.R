rm(list = ls())
library(dplyr)
library(ggplot2)
library(scales)
library(mgcv)
library(rstan)

load(file = 'grBspline.Rdata')
df = read.csv('Covid19CasesGISAID.csv')
df$Var1 = as.Date(df$Var1)
df = df[df$Var1 < as.Date('2021-12-01'),]

df1 = df[df$Mutations == 'Lineage A',]
df2 = df[df$Mutations == 'Lineage B',]
getsimu = function(df, fit){
  ndays = nrow(df)
  time = 1:ndays
  cases =  df$Freq
  # Fit a GAM with P-splines
  gam_model <- gam(cases ~ s(time, bs = "ps"), 
                   family = poisson(link = "log"))
  # Extracting model matrix
  X <- model.matrix(gam_model)
  posterior_samples <- extract(fit)
  simu_Onset <- matrix(NA, nrow = 15000, ncol = ndays)
  for (i in 1:15000) {
    pars = posterior_samples$beta[i,]
    simu_Onset[i, ] <- exp(X %*% pars)  
  }
  return(simu_Onset)
}

# Function to calculate the growth rate
calculate_growth_rate <- function(simu) {
  # Set a threshold below which values are considered zero
  threshold <- 1
  simu[simu < threshold] <- NA
  
  # Initialize the matrix to store growth rates
  growth_rates <- matrix(NA, nrow = nrow(simu), ncol = ncol(simu) - 1)
  
  for (i in 2:ncol(simu)) {
    valid <- !is.na(simu[, i]) & !is.na(simu[, i - 1]) & simu[, i - 1] > 0
    growth_rates[valid, i - 1] <- log(simu[valid, i] / simu[valid, i - 1])
  }
  
  return(growth_rates)
}

get_stats <- function(data) {
  mean_data <- apply(data, 2, function(x) mean(x, na.rm = TRUE))
  ci_lower <- apply(data, 2, function(x) quantile(x, probs = 0.025, na.rm = TRUE))
  ci_upper <- apply(data, 2, function(x) quantile(x, probs = 0.975, na.rm = TRUE))
  return(list(mean = mean_data, lower = ci_lower, upper = ci_upper))
}


get_plot_data = function(df, fit){

  simu = getsimu(df, fit)
  gr <- calculate_growth_rate(simu)
  stats <- get_stats(gr)
  
  # Convert to data frame for ggplot
  days <- 1:(ncol(simu) - 1)
  start_date <- as.Date("2019-12-24") 
  date_vector <- seq.Date(start_date, by = "day", length.out = ncol(simu) - 1)
  
  plot_data <- data.frame(
    Day = days,
    date_vector = date_vector,
    Fitted = stats$mean,
    LowerCI = stats$lower,
    UpperCI = stats$upper
  )
  
  return(plot_data)
}

plot_data1 = get_plot_data(df1, fit1)
plot_data2 = get_plot_data(df2, fit2)

plot_data1$Mutations = 'Lineage A'
plot_data2$Mutations = 'Lineage B'
plot_data = rbind(plot_data1, plot_data2)

values = c(hue_pal()(3)[1], hue_pal()(3)[3])
anno = data.frame(x0 = as.Date(c('2020-7-20','2020-9-1','2020-10-29','2020-10-30','2021-8-10')),
                  y0 = c(0.08,0.12,0.16,-0.06,0.08),
                  text = c('Beta', 'Alpha', 'Gamma', 'Delta', 'Omicron'))
anno_arr = data.frame(x = as.Date(c('2020-8-10','2020-9-20','2020-10-29','2020-10-30','2021-9-1')),
                      y = c(0.06,0.1,0.14,-0.04,0.06),
                      xend = as.Date(c('2020-8-10','2020-9-20','2020-10-29','2020-10-30','2021-9-1')),
                      yend = c(0.02,0.04,0.06,0,0.02),
                      text = c('Beta', 'Alpha', 'Gamma', 'Delta', 'Omicron'))

# Plotting
p = ggplot(plot_data, aes(x = date_vector, y = Fitted, group = Mutations)) +
  geom_line(aes(color = Mutations)) +
  geom_ribbon(aes(ymin = LowerCI, ymax = UpperCI, fill = Mutations), 
             alpha = 0.3) +
  labs(x = "Day", y = "Growth Rate") +
  scale_x_date(breaks = seq(as.Date('2020-01-01'), as.Date('2022-11-01'), by="6 months"),
               minor_breaks = seq(as.Date('2019-12-01'), as.Date('2022-11-01'), by ='1 month'),
               date_labels = "%y-%b") +
  xlab('Date (2019-2021)') +
  ylab('Growth rate') + theme_bw() +
  scale_color_manual(name="",
                     labels=c("A", "B"),
                     values = alpha(values, 0.6)) +
  scale_fill_manual(name="",
                    labels=c("A", "B"),
                    values = alpha(values, 0.1)) +
  coord_cartesian(ylim = c(-0.1, max(plot_data$Fitted)), 
                  xlim = c(min(plot_data$date_vector), as.Date('2021-10-05'))) +
  geom_text(data = anno, aes(x = x0,  y = y0, label = text),
            inherit.aes = F, size = 3, alpha = .7) +
  geom_segment(data = anno_arr, aes(x = x, y = y, xend = xend, yend = yend),
               arrow = arrow(length = unit(0.2, "cm")), alpha = .7,
               inherit.aes = F) +
  theme(legend.background = element_blank(),
        legend.position = c(0.8,0.8),
        legend.title = element_blank())

pdf(paste0("Output/gr_rate.pdf"), width = 2.8, height = 2.3)
print(p)
dev.off()


