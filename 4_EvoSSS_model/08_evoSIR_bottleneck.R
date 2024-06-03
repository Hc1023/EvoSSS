rm(list = ls())
library(deSolve)   
library(ggplot2)   
library(dplyr)
library(tidyr)
library(scales)
library(ggnewscale)
library(ggpubr)

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
  # state <- c(V1 = 1, V2 = 1)
  # Run simulations for each parameter set
  results <- lapply(seq(nrow(param_sets)), function(i) {
    params <- unlist(param_sets[i, ])
    out <- ode(y = c(params['V1'],params['V2']), times = times, 
               func = viral_model, parms = params)
    out_df <- as.data.frame(out)
    out_df$group <- params['group']
    return(out_df)
  })
  combined_results = bind_rows(results)
  
  return(combined_results)
  
}

transform_data = function(combined_results){
  long_data <- pivot_longer(combined_results, cols = c("V1", "V2"), 
                            names_to = "Strain", values_to = "Population")
  long_data$Population_label = long_data$Population/max(long_data$Population)
  long_data$group = factor(long_data$group, levels = unique(long_data$group))
  return(long_data) 
}
ratio_fun = function(combined_results){
  data_ratio <- combined_results %>%
    mutate(ratio = V1/(V1+V2), group = group) %>%
    select(time, group, ratio)
  return(data_ratio)
}
# Parameters
load(file = 'evoSIR.rdata')
fit = fitlist[[4]]
posterior_samples <- rstan::extract(fit)

c = posterior_samples$beta1/posterior_samples$beta2
deltar = -log(c)


n = 5000
seed = 20
pA = 0.25

getdf = function(n, seed, pA){
  pA_vec = rep(pA, n)
  seed_vec = rpois(n,seed)
  seed_vec = seed_vec[seed_vec>0]
  group = c(0.02,0.05,0.1,0.2)
  params = c(r1 = 0.2, 
             r2 = 0.2 + group[2], 
             K = seed*2500, 
             alpha12 = 2, 
             alpha21 = 0, 
             mu = 0)
  state_vec = sapply(1:length(seed_vec), function(x){
    y = rmultinom(1, size = seed_vec[x], 
                  prob = c(pA_vec[x], 1-pA_vec[x]))
    return(matrix(y))
  })
  times = seq(0,48,2)
  df = data.frame()
  for (i in 1:length(seed_vec)) {
    out <- ode(y = c(state_vec[1,i], state_vec[2,i]), 
               times = times, func = viral_model, 
               parms = params)
    out = data.frame(out)
    pdf = out[out$time %in% c(2,8,24,48),]
    if(is.na(pdf$X1[1]/(pdf$X1[1]+pdf$X2[1]))){
      return(state_vec[,i])
      break
    }
    df = rbind(df, pdf)
  }
  
  df$p = df$X1/(df$X1 + df$X2)
  df$time = factor(df$time)
  return(df)
}

df1 = getdf(n = 5000, seed = 20, pA = 0.2)
df2 = getdf(n = 5000, seed = 5, pA = 0.2)
df3 = getdf(n = 5000, seed = 20, pA = 0.4)
df4 = getdf(n = 5000, seed = 5, pA = 0.4)
save(df1, df2, df3, df4, file = 'evoSIR_bottleneck.rdata')
getplot = function(df, p0){

  means <- df %>%
    group_by(time) %>%
    summarize(mean_p = mean(p, na.rm = T))
  values = c('#98afc7','#82caff','#b8cbfb', 
             '#737ca1')
  p = ggplot(df, aes(p, fill = time)) + 
    geom_density(alpha = 0.3, bw = 0.1) +
    geom_vline(data = means, aes(xintercept = mean_p, color = time), 
               linetype = "dashed") +
    geom_vline(xintercept = p0, color = 'black', 
               linetype = "dashed") +
    ylab('') + xlab('') +
    scale_y_continuous(breaks = c(0:3),
                       labels = 0:3,
                       limits = c(0,3.8)) +
    scale_x_continuous(breaks = seq(0,1,0.5),
                       limits = c(0,1)) +
    theme_bw() +
    scale_fill_manual(name = '',
                      values = values) +
    scale_color_manual(name = '',
                      values = values) +
    theme(plot.background = element_blank(),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          legend.position = 'right',
          legend.background = element_blank())
  p
  return(p)
}

pdf(paste0("Output/withinhost_bottleneck.pdf"), 
    width = 2.5, height = 1.2)

print(getplot(df1, 0.2))
print(getplot(df2, 0.2))
print(getplot(df3, 0.4))
print(getplot(df4, 0.4))

dev.off()

