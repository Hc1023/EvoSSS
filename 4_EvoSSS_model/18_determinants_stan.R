rm(list = ls())
library(rstan)
library(ggplot2)
library(scales)
library(ggnewscale)
library(tidyverse)
library(dplyr)

df = read.csv('../3_Epidemiological_analysis/VOC_gisaid.csv')

VOC = c("Delta","Alpha","D614G","B 19A","A")
names(VOC) = c("GK","GRY","G","S","L")
df = df[df$GISAID_clade %in% names(VOC),]
df$V = NA
for (i in 1:length(VOC)) {
  df$V[df$GISAID_clade == names(VOC)[i]] = VOC[i]
}
df$V = factor(df$V, levels = VOC)
dfall = df %>% group_by(date) %>%
  summarise(yall = sum(count))

df <- na.omit(df)

full_dates <- as.Date('2020-07-01') + 0:800
# Create a dataframe with all dates
full_df <- expand.grid(date = full_dates, 
                       V = unique(df$V))
df$date = as.Date(df$date)
merged_df <- full_df %>%
  left_join(df, by = c("date", "V")) %>%
  mutate(count = ifelse(is.na(count), 0, count))
voc = unique(merged_df$V)
observed_matrix = data.frame(v1 = merged_df[merged_df$V == 'D614G', 'count'],
                             v2 = merged_df[merged_df$V == 'Alpha', 'count'],
                             v3 = merged_df[merged_df$V == 'Delta', 'count'])
rownames(observed_matrix) = full_dates


fexpect0 = data.frame(x = rep(1:length(full_dates),3), 
                      y = c(observed_matrix[,1],
                            observed_matrix[,2],
                            observed_matrix[,3]),
                      group = factor(rep(c('D614G','Alpha','Delta'), 
                                         each = length(full_dates)),
                                     levels = c('D614G','Alpha','Delta')))
ggplot(data = fexpect0, 
       aes(x = x, y = y, 
           group = group, color = group)) +
  geom_line() +
  scale_y_continuous(trans = 'log10')

update_fun <- function(pars, states_old, N, n) {
  S <- states_old[,1]
  I <- states_old[,2:(n+1)]
  Onset <- states_old[,(n+2):(2*n+1)]
  
  betas = pars
  
  gamma = 0.157
  
  S_new <- S
  for (i in 1:n) {
    S_new <- S_new - betas[i] * S * I[,i] / N
  }
  
  I_new <- matrix(0, nrow(states_old), n)
  I_new[,1] <- I[,1] + betas[1] * S * I[,1] / N - gamma * I[,1] 
  for (i in 2:n) {
    I_new[,i] <- I[,i] + betas[i] * S * I[,i] / N - gamma * I[,i] 
  }
  
  Onset_new <- matrix(0, nrow(states_old), n)
  Onset_new[,1] <- betas[1] * S * I_new[,1] / N
  
  for (i in 2:n) {
    Onset_new[,i] <- betas[i] * S * I_new[,i] / N 
  }
  
  mat <- cbind(S_new, I_new, Onset_new)
  
  return(mat)
}

simu <- function(seed_mats, N, poolday, pars, n) {
  # Initial conditions
  I_old <- matrix(0, ncol=n, nrow=nrow(seed_mats[[1]]))
  S_old <- N - rowSums(I_old)
  
  for (i in 1:n) {
    I_old[,i] <- seed_mats[[i]][1,]
  }
  
  states_old <- cbind(S_old, I_old, matrix(0, ncol=n, nrow=nrow(seed_mats[[1]])))
  
  ndays <- 200
  Onsets_mat <- matrix(0, ndays, n)
  
  for (t in 2:ndays) {
    
    states_old <- update_fun(pars = pars, states_old = states_old, N = N, n = n)
    if (t <= nrow(seed_mats[[1]])) {
      for (i in 1:n) {
        states_old[t,(i+1):(i+1)] <- states_old[t,(i+1):(i+1)] + seed_mats[[i]][t,t]
      }
    }
    Onsets_mat[t,] <- colSums(states_old[,(n+2):(2*n+1)])
  }
  
  return(Onsets_mat)
}

fitlist = list()


determinant_fun = function(cond = T){
  
  poolday = 30
  # The initial cycle - epidemic outbreak
  n= 3
  restrict_coef = 0.7
  
  seed_vec = c(34,rep(0, poolday-1))
  probs <- c(1,1e-3,1e-3)
  
  seed_matrix <- matrix(0, nrow=n, ncol=length(seed_vec))
  for (i in 1:n) {
    seed_matrix[i,] <- seed_vec * probs[[i]]
  }
  seed_mats <- list()
  
  for (i in 1:n) {
    seed_mats[[i]] <- diag(seed_matrix[i,])
  }
  
  nday = 100
  
  expected_matrix = observed_matrix[1:nday,] 
  expected_matrix[expected_matrix < 0] = 0
  expected_matrix = round(expected_matrix)
  expected_matrix[round(restrict_coef*nday+1):nday,] = 0
  
  if(cond){
    stan_data <- list(
      poolday = poolday,
      nday = nday,
      expected_matrix = expected_matrix,
      pars_last = c(200,0.3,0.3,0.3),
      seed_mat_I1 = seed_mats[[1]],
      seed_mat_I2 = seed_mats[[2]],
      seed_mat_I3 = seed_mats[[3]],
      seed_vec = seed_vec,
      gamma = 0.157
    )
    # Fit the model
    fit <- stan(file = 'VOC.stan', data = stan_data, 
                iter = 10000, chains = 1, warmup = 5000,
                verbose = TRUE)
    fitlist[[1]] = fit
  }
  
  
  fexpect = data.frame(x = rep(1:nday,3), 
                       y = c(expected_matrix[1:nday,1],
                             expected_matrix[1:nday,2],
                             expected_matrix[1:nday,3]),
                       group = factor(rep(c('D614G','Alpha','Delta'), 
                                          each = nday),
                                      levels = c('D614G','Alpha','Delta')))
  
  
  Onsets_mat_list = list()
  fit= fitlist[[1]]
  posterior = rstan::extract(fit)
  pars_last = c(mean(posterior$contact), mean(posterior$beta1),
                mean(posterior$beta2), mean(posterior$beta3))
  Onsets_mat = simu(seed_mats, 
                    N = seed_vec * pars_last[1] + 1, 
                    poolday, pars = c(pars_last[2], pars_last[3], pars_last[4]), n)
  fonset = data.frame(x = rep(1:nday,3), 
                      y = c(Onsets_mat[1:nday,1],
                            Onsets_mat[1:nday,2],
                            Onsets_mat[1:nday,3]),
                      group = factor(rep(c('D614G','Alpha','Delta'), 
                                         each = nday), 
                                     levels = c('D614G','Alpha','Delta')))
  
  
  ggplot() +
    geom_point(data = fexpect, 
               aes(x = x, y = y, group = group, color = group)) +
    geom_line(data = fonset,
              aes(x = x, y = y, group = group, color = group))
  
  
  Onsets_mat_list[[1]] = Onsets_mat
  cond = T
  cond = F
  for (j in 19:24) {
    print(j)
    Onsets_mat = Onsets_mat_list[[j]]
    
    # Generalized extraction of Onset columns for all variants
    Onsets <- list()
    for (i in 1:n) {
      Onsets[[i]] <- Onsets_mat[poolday + 1:poolday, i] + 1e-10
    }
    
    
    # Mobility control force
    mobility <- rep(0.1, 30) 
    Mobility_matrix <- diag(mobility)
    
    # Create seed_vec considering all variants
    seed_vec <- rowSums(sapply(Onsets, function(Onset) Onset)) %*% Mobility_matrix %>% as.numeric()
    
    # Calculate probabilities for each variant
    probs <- lapply(Onsets, function(Onset) Onset / rowSums(do.call(cbind, Onsets)))
    seed_matrix <- matrix(0, nrow=n, ncol=length(seed_vec))
    
    for (i in 1:n) {
      seed_matrix[i,] <- seed_vec * probs[[i]]
    }
    
    # Create seed_matrices for each variant
    seed_mats <- list()
    
    for (i in 1:n) {
      seed_mats[[i]] <- diag(seed_matrix[i,])
    }
    
    
    expected_matrix = observed_matrix[poolday*j+1:nday,] - Onsets_mat[poolday+1:nday,]
    expected_matrix[expected_matrix<0] = 0
    expected_matrix = round(expected_matrix)
    
    expected_matrix[round(restrict_coef*nday+1):nday,] = 0
    fexpect = data.frame(x = rep(1:nday,3), 
                         y = c(expected_matrix[1:nday,1],
                               expected_matrix[1:nday,2],
                               expected_matrix[1:nday,3]),
                         group = factor(rep(c('D614G','Alpha','Delta'), 
                                            each = nday),
                                        levels = c('D614G','Alpha','Delta')))
    
    if(cond){
      fit = fitlist[[j]]
      posterior = rstan::extract(fit)
      pars_last = c(mean(posterior$contact), mean(posterior$beta1),
                    mean(posterior$beta2), mean(posterior$beta3))
      stan_data <- list(
        poolday = poolday,
        nday = nday,
        expected_matrix = expected_matrix,
        pars_last = pars_last,
        seed_mat_I1 = seed_mats[[1]],
        seed_mat_I2 = seed_mats[[2]],
        seed_mat_I3 = seed_mats[[3]],
        seed_vec = seed_vec,
        gamma = 0.157
      )
      # Fit the model
      fit <- stan(file = 'VOC.stan', data = stan_data, 
                  iter = 10000, chains = 1, warmup = 5000,
                  verbose = TRUE)
      fitlist[[j+1]] = fit
    }
    
    
    fit = fitlist[[j+1]]
    posterior = rstan::extract(fit)
    pars_last = c(mean(posterior$contact), mean(posterior$beta1),
                  mean(posterior$beta2), mean(posterior$beta3))
    Onsets_mat = simu(seed_mats, 
                      N = seed_vec * pars_last[1] + 1, 
                      poolday, pars = c(pars_last[2], pars_last[3], pars_last[4]), n)
    fonset = data.frame(x = rep(1:nday,3), 
                        y = c(Onsets_mat[1:nday,1],
                              Onsets_mat[1:nday,2],
                              Onsets_mat[1:nday,3]),
                        group = factor(rep(c('D614G','Alpha','Delta'), 
                                           each = nday), 
                                       levels = c('D614G','Alpha','Delta')))
    
    ggplot() +
      geom_point(data = fexpect, 
                 aes(x = x, y = y, group = group, color = group)) +
      geom_line(data = fonset,
                aes(x = x, y = y, group = group, color = group))
    
   
     Onsets_mat_list[[j+1]] = Onsets_mat
  }
  
  save(fitlist, file = 'voc_fitlist.rdata')
  dfplot_simu = data.frame()
  for (i in 1:25) {
    Onsets_mat = Onsets_mat_list[[i]]
    k = i-1
    
    dfplot_simu1 = data.frame(x = rep(1:nrow(Onsets_mat)+poolday*k,n),
                              y = as.vector(Onsets_mat),
                              group = rep(paste0(1:n,'_',k), 
                                          each = nrow(Onsets_mat)),
                              color = rep(1:n, 
                                          each = nrow(Onsets_mat)))
    dfplot_simu = rbind(dfplot_simu, dfplot_simu1)
  }
  
  df2 = dfplot_simu %>% group_by(x, color) %>% 
    summarise(y = sum(y)) %>%
    as.data.frame()
  
  
  # Normalizing proportions for each date
  data <- df2 %>%
    group_by(x) %>%
    mutate(p = y/sum(y)) %>%
    as.data.frame()
  
  levels = rev(unique(data$color))
  
  data$color = factor(data$color,levels = levels)
  
  return(data)
}

data = determinant_fun()

parsmat = matrix(0.3,nrow = 25, ncol = 3)
parsmat[10:25,2] = 0.36
parsmat[15:25,3] = 0.5
data = determinant_fun(parsmat)
ggplot() +
  geom_area(data = data, 
            aes(x = x, y = p, fill = color),
            position = 'fill')

dfall = data %>% group_by(x) %>% 
  summarise(yall = sum(y))

ggplot() +
  geom_line(data = data, 
            aes(x = x, y = y, group = color, color= color)) +
  geom_line(data = dfall, aes(x = x, y = yall)) +
  scale_y_continuous(trans='log10') +
  coord_cartesian(ylim = c(2,max(data$y))) +
  labs(x = "Date", y = "Proportion") +
  theme_minimal() +
  theme(legend.position = "right",
        plot.title = element_text(hjust = 0.5))

