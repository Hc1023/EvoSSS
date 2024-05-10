rm(list = ls())
library(deSolve)   
library(ggplot2)   
library(dplyr)
library(tidyr)
library(scales)
# Define the model
viral_model <- function(t, state, parameters) {
  V1 <- state[1]
  V2 <- state[2]
  
  r1 <- parameters["r1"]
  r2 <- parameters["r2"]
  K <- parameters["K"]
  alpha12 <- parameters["alpha12"]
  alpha21 <- parameters["alpha21"]
  mu <- parameters["mu"]
  

  dV1dt <- r1 * V1 * (1 - (V1 + alpha21 * V2) / K) - mu * V1
  dV2dt <- r2 * V2 * (1 - (V2 + alpha12 * V1) / K) - mu * V1
  
  return(list(c(dV1dt, dV2dt)))
}

withinhost_fun = function(param_sets){
  # Initial state and time sequence
  state <- c(V1 = 1, V2 = 1)
  times <- seq(0, 72, by = 1)
  
  # Run simulations for each parameter set
  results <- lapply(seq(nrow(param_sets)), function(i) {
    params <- unlist(param_sets[i, ])
    out <- ode(y = state, times = times, func = viral_model, parms = params)
    out_df <- as.data.frame(out)
    
    return(out_df)
  })
  combined_results = bind_rows(results)
  
  return(combined_results)
  
}



# Parameters
load(file = 'evoSIR.rdata')
fit = fitlist[[4]]
posterior_samples <- extract(fit)

c = posterior_samples$beta1/posterior_samples$beta2
deltar = -log(c)

rbase = log2(exp(1))/3
Kbase = 2^18
param_sets <- expand.grid(r1 = rbase, 
                          r2 = rbase + deltar[1:50], 
                          K = Kbase, 
                          alpha12 = 1, 
                          alpha21 = 0, 
                          mu = 0)

combined_results = withinhost_fun(param_sets)
combined_results$r = combined_results$V1/(combined_results$V1 + combined_results$V2)

dataplot = data.frame()
for (i in 0:72) {
  df = combined_results[combined_results$time == i,]
  y = sapply(2:4, function(x){
    v = df[,x]
    m = mean(v)
    ci_lower = quantile(v, probs = 0.025, na.rm = T)
    ci_upper = quantile(v, probs = 0.975, na.rm = T)
    return(c(m=m, ci_lower=ci_lower, ci_upper = ci_upper))
  })
  
  dfout = data.frame(time = rep(i,3), group = c('A','B','r'))
  dfout = cbind(dfout,t(y))
  dataplot = rbind(dataplot, dfout)
}

colnames(dataplot)[4:5] = c('LowerCI', 'UpperCI')
dataplot[dataplot$group != 'r',c(3:5)] = dataplot[dataplot$group != 'r',c(3:5)] / max(dataplot[dataplot$group != 'r',c(3:5)])
dataplot$group = factor(dataplot$group, levels = c('A','B','r'))

values = c(hue_pal()(3)[1], hue_pal()(3)[3], "#9ab9c6")
p = ggplot(dataplot, aes(x = time, y = m, group = group)) +
  geom_ribbon(data = filter(dataplot, group != "r"),
              aes(ymin = LowerCI, ymax = UpperCI, fill = group),  
              show.legend = F) +
  geom_line(data = filter(dataplot, group != "r"), 
            aes(x = time, y = m, color = group)) + 
  geom_line(data = filter(dataplot, group == "r", time %in% c(0,24,48,72)), 
            aes(x = time, y = m, color = group)) + 
  geom_point(data = filter(dataplot, group == "r", time %in% c(0,24,48,72)), 
            aes(x = time, y = m, color = group)) + 
  theme_bw() + 
  theme(legend.position = 'none',
        legend.background = element_blank(),
        legend.key = element_blank(),
        axis.text.y.right = element_text(color = values[3]),
        axis.title.y.right = element_text(color = values[3])) +
  scale_color_manual(name="",
                     labels=c("A", "B","Ratio"),
                     values = alpha(values, 0.9)) +
  scale_fill_manual(name="",
                    labels=c("A", "B","Ratio"),
                    values = alpha(values, 0.3)) +
  scale_x_continuous(breaks = c(0,24,48,72)) + 
  scale_y_continuous(limits = c(0, 1), 
                     minor_breaks = seq(0 , 1, 0.25),
                     breaks = c(0,0.5,1),
                     labels = c('0','50','100'),
                     name = 'Population \n(% of MAX)',
                     sec.axis = sec_axis(
                       transform = ~ ., 
                       name = "Ratio of A",
                       breaks = seq(0, 1, by = 0.5)
                     )) +
  xlab('Time unit')

pdf(paste0("Output/withinhostWH.pdf"), width = 2.8, height = 1.5)
print(p)
dev.off()
