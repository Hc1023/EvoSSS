rm(list = ls())
library(deSolve)   
library(ggplot2)   
library(dplyr)
library(tidyr)
library(scales)
library(ggnewscale)
library(ggpubr)

# Define the model
viral_model <- function(state, params) {
  V1 <- state[1]
  V2 <- state[2]
  
  r1 <- params["r1"]
  r2 <- params["r2"]
  K <- params["K"]
  alpha12 <- params["alpha12"]
  alpha21 <- params["alpha21"]
  mu <- params["mu"]

  dV1dt <- r1 * V1 * (1 - (V1 + alpha21 * V2) / K)
  dV2dt <- r2 * V2 * (1 - (V2 + alpha12 * V1) / K)
  dV1dt = max(0, dV1dt)
  dV2dt = max(0, dV2dt)
  total_change_rate <- dV1dt + dV2dt
  if (total_change_rate > 0) {
    pV1 <- dV1dt / total_change_rate
    pV2 <- dV2dt / total_change_rate
    # Multinomial sampling
    size = floor(total_change_rate) + 
      as.numeric(runif(1) < (total_change_rate - floor(total_change_rate)))
    changes <- rmultinom(1, size = size, prob = c(pV1, pV2))
    
    # Update populations based on sampled changes
    V1 <- V1 + changes[1] 
    V2 <- V2 + changes[2]
  }
  
  return(c(V1, V2))
}



getdf = function(seed, params, n, p){
  # Time points
  times <- 1:48
  selected_times = c(2, 8, 24, 48)
  # Initial state
  # n = 100; seed = 20; p = 0.5
  seed_mat = rmultinom(n, seed, c(p,1-p))
  V1_values <- numeric(length(times))
  V2_values <- numeric(length(times))
  p = sapply(1:n, function(j){
    state = seed_mat[,j]
    for (i in times) {
      state <- viral_model(state, params)
      V1_values[i] <- state[1]
      V2_values[i] <- state[2]
    }
    p = V1_values[selected_times]/(V1_values[selected_times] + V2_values[selected_times])
    return(p)
  })

  df = data.frame(time = rep(selected_times, n),
                  p = c(p))

  return(df)
}

seed = 3


if(F){

  seed = 20
  params = c(r1 = 0.2, 
             r2 = 0.24, 
             K = 200*seed, 
             alpha12 = 1, 
             alpha21 = 0, 
             mu = 0)
  df = getdf(seed = seed, params = params, n = 5000, p = 0.3)
  df1 = df
  seed = 6
  params = c(r1 = 0.2, 
             r2 = 0.24, 
             K = 200*seed, 
             alpha12 = 1, 
             alpha21 = 0, 
             mu = 0)
  df = getdf(seed = seed, params = params, n = 5000, p = 0.3)
  df2 = df
  seed = 3
  params = c(r1 = 0.2, 
             r2 = 0.24, 
             K = 200*seed, 
             alpha12 = 1, 
             alpha21 = 0, 
             mu = 0)
  df = getdf(seed = seed, params = params, n = 5000, p = 0.3)
  df2 = df
  save(df1, df2, df3, file = 'evoSIR_bottleneck.rdata')
}
load('evoSIR_bottleneck.rdata')
plotfun = function(df){
  df$time = factor(df$time)
  means <- df %>%
    group_by(time) %>%
    summarize(mean_p = mean(p, na.rm = T))
  values = c('#98afc7','#79baec','#0041c2', 
             '#a37ca1')
  p = ggplot(df, aes(p, fill = time)) + 
    geom_density(alpha = 0.3, bw = 0.1) +
    geom_vline(data = means, aes(xintercept = mean_p, color = time), 
               linetype = "dashed") +
    geom_vline(xintercept = 0.3, color = 'black', 
               linetype = "dashed") +
    ylab('') + xlab('') +
    scale_y_continuous(n.breaks = 3) +
    scale_x_continuous(breaks = seq(0,1,0.25),
                       labels = c('0.00','0.25','0.50','0.75','1.00'),
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
          legend.position = 'none',
          legend.key.size = unit(0.2,'cm'),
          legend.key.width = unit(0.3,'cm'),
          legend.background = element_blank())
  return(p)
}

p1 = plotfun(df1)
p2 = plotfun(df2)
p3 = plotfun(df3)

pdf(paste0("Output/withinhost_bottleneck.pdf"), 
    width = 1.8, height = 1.2)
print(p1)
print(p2)
print(p3)
dev.off()

plotfun2 = function(df, p0){
  df$seed = factor(df$seed, levels = c(50,20,10,6,3))
  means <- df %>%
    group_by(seed) %>%
    summarize(mean_p = mean(p, na.rm = T))
  values = c( '#cb8335','#a3a637',
              '#68a588','#579aa7','#8781ba')
  p = ggplot(df, aes(p, fill = seed, color = seed)) + 
    geom_density(alpha = 0.3, bw = 0.1) +
    geom_vline(data = means,
               aes(xintercept = mean_p, color = seed), 
               linetype = "dashed") +
    geom_vline(xintercept = p0, color = 'black', 
               linetype = "dashed") +
    ylab('') + xlab('') +
    scale_y_continuous(n.breaks = 3) +
    scale_x_continuous(breaks = seq(0,1,0.25),
                       labels = c('0.00','0.25','0.50','0.75','1.00'),
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
          legend.key.size = unit(0.2,'cm'),
          legend.key.width = unit(0.3,'cm'),
          legend.background = element_blank())
  return(p)
}

if(F){
  n = 1000
  dfall = data.frame()
  for (seed in c(3,6,10,20,50)) {
    print(seed)
    params = c(r1 = 0.2, 
               r2 = 0.2, 
               K = 200*seed, 
               alpha12 = 1, 
               alpha21 = 1, 
               mu = 0)
    df = getdf(seed = seed, params = params, n = n, p = 0.3)
    df$seed = seed
    dfall = rbind(dfall, df)
  }
  
  n = 1000
  dfall2 = data.frame()
  for (seed in c(3,6,10,20,50)) {
    print(seed)
    params = c(r1 = 0.2, 
               r2 = 0.2, 
               K = 200*seed, 
               alpha12 = 1, 
               alpha21 = 1, 
               mu = 0)
    df = getdf(seed = seed, params = params, n = n, p = 0.5)
    df$seed = seed
    dfall2 = rbind(dfall2, df)
  }
}

p1 = plotfun2(dfall, p0 = 0.3)
p2 = plotfun2(dfall2, p0 = 0.5)

pdf(paste0("Output/withinhost_bottleneck_seed.pdf"), 
    width = 1.8, height = 1.4)
print(p1)
print(p2)
dev.off()
