rm(list = ls())
library(rstan)
library(ggplot2)
library(scales)
realData_all <- read.csv("evoSIR/F3D_Covid19CasesWH.csv", row.names = 1)
realData <- realData_all[-c(1:24), ] 
observed_cases = realData$CaseNum
load(file = 'evoSIR/F3D_evoSIR.rdata')
update_fun = function(pars, states_old){
  
  S <- states_old[1]
  I1 <- states_old[2]
  I2 <- states_old[3]
  
  beta1 = pars[1]
  beta2 = pars[2]
  gamma = pars[3]
  
  S_new = S - beta1*S*I1/N - beta2*S*I2/N
  I1_new = I1 + beta1*S*I1/N - gamma*I1
  I2_new = I2 + beta2*S*I2/N - gamma*I2
  Onset1 = beta1*S*I1/N
  Onset2 = beta2*S*I2/N
  
  return(c(S_new, I1_new, I2_new, Onset1, Onset2))
}

update_fun_stochastic = function(pars, states_old){
  
  S <- states_old[1]
  I1 <- states_old[2]
  I2 <- states_old[3]
  
  beta1 = pars[1]
  beta2 = pars[2]
  gamma = pars[3]
  
  pS_vec = c(max(0,beta1*I1/N), max(0,beta2*I2/N), max(0,1-beta1*I1/N-beta2*I2/N))
  sample_S <- rmultinom(1, size = S, prob = pS_vec)
  pI_vec <- c(gamma, 1-gamma)
  sample_I1 <- rmultinom(1, size = I1, prob = pI_vec)
  sample_I2 <- rmultinom(1, size = I2, prob = pI_vec)
  
  ## new values
  S_new <- sample_S[3]
  I1_new <- sample_I1[2] + sample_S[1] 
  I2_new <- sample_I2[2] + sample_S[2] 
  
  Onset1 <- sample_S[1]
  Onset2 <- sample_S[2]
  
  return(c(S_new, I1_new, I2_new, Onset1, Onset2))
}

simu <- function(pars, states_old, ndays, f = update_fun) {
  # Initial conditions
  states_old = c(states_old, NA, NA)
  
  # Simulate the dynamics over ndays
  mycol <- c("time", "S", "I1", "I2", "Onset1", "Onset2")
  states_mat <- matrix(0, ndays, length(mycol))
  states_mat[, 1] <- 1:ndays
  states_mat[1, 2:length(mycol)] <- states_old
  
  colnames(states_mat) <- mycol
  
  for (t in 1:(ndays-1)) {
    states_mat[t+1, -1] <- f(pars = pars, states_old = states_old)
    states_old <- states_mat[t+1, -1]
  }
  
  return(states_mat)
}

generate_plotdata = function(df){
  simu_Onset = df
  # Calculate the 2.5th and 97.5th percentiles for the confidence interval
  ci_lower <- apply(simu_Onset, 2, quantile, probs = 0.025, na.rm = T)
  ci_upper <- apply(simu_Onset, 2, quantile, probs = 0.975, na.rm = T)
  
  plot_data <- data.frame(
    Day = 1:ndays,
    date_vector = date_vector,
    # Observed = observed_cases,
    Fitted = colMeans(simu_Onset),
    LowerCI = ci_lower,
    UpperCI = ci_upper
  )
  
  return(plot_data)
}


if(F){
  id = 1:9
  # summarylist = list()
  # plist = list()
  plot_data_list = list()
  for (j in 1:9) {
    print(j)
    fit = fitlist[[id[j]]]
    posterior_samples <- extract(fit)
    
    # Define data
    {
      N = sum(observed_cases)
      CaseNum = realData_all$CaseNum
      observed_cases = CaseNum[-c(1:24)]
      ratio = id[j]/10
      I10 = round(sum(CaseNum[22:24])*ratio)
      I20 = round(sum(CaseNum[22:24])*(1-ratio))
      S0 <- N - I10 - I20
      ndays <- length(observed_cases)
    }
    
    simu_Onset1 <- matrix(NA, nrow = 5000, ncol = ndays)
    simu_Onset2 <- matrix(NA, nrow = 5000, ncol = ndays)
    # Simulate the epidemic for each set of sampled parameters
    for (i in 1:5000) {
      
      beta1 = posterior_samples$beta1[i]
      beta2 = posterior_samples$beta2[i]
      gamma = posterior_samples$gamma[i]
      pars = c(beta1 = beta1, beta2 = beta2, gamma = gamma)
      states_old = c(S0,I10,I20)
      result <- simu(pars, states_old, ndays,
                     f = update_fun_stochastic)
      simu_Onset1[i, ] <- result[,'Onset1']
      simu_Onset2[i, ] <- result[,'Onset2']
    }
    
    start_date <- as.Date("2020-01-01") 
    date_vector <- seq.Date(start_date, by = "day", length.out = ndays)
    
    
    data_observed <- data.frame(
      Day = 1:ndays,
      date_vector = date_vector,
      Observed = observed_cases
    )
    simu_Onset = simu_Onset1 + simu_Onset2
    df1 = generate_plotdata(df = simu_Onset1)
    df2 = generate_plotdata(df = simu_Onset2)
    df = generate_plotdata(df = simu_Onset)
    df1$V = 'V1'
    df2$V = 'V2'
    df$V = 'All'
    plot_data = rbind(df1, df2, df)
    plot_data_list[[j]] = plot_data
    # r = sapply(1:nrow(simu_Onset), function(i){
    #   sum(simu_Onset1[i,], na.rm = T)/sum(simu_Onset[i,], na.rm = T)
    # })
    
    # sum(df1$Fitted, na.rm = T)
    # sum(df2$Fitted, na.rm = T)
    
    # print(fit)
    # dfsummary = data.frame(fit@.MISC$summary$msd)[1:3,]
    # dfsummary['r',1:2] =  c(mean(r), sd(r))
    # summarylist[[j]] = dfsummary
    
  }
  save(plot_data_list, file = 'evoSIR/evoSIR_interhost_plotdata.rdata')
}
load(file = 'evoSIR/F4D_evoSIR_interhost_plotdata.rdata')
values = c("#aa85a6", hue_pal()(3)[1], hue_pal()(3)[3])
pdf(paste0("Output/F4D_acrosshost.pdf"), width = 2.4, height = 2.2)
for (j in 1:9) {
  
  p = ggplot() +
    geom_ribbon(data = plot_data_list[[j]], 
                aes(x = date_vector, group = V, 
                    ymin = LowerCI, ymax = UpperCI, fill = V)) +  # Confidence interval
    geom_point(data = data_observed, 
               aes(x = date_vector, y = Observed), 
               size = 0.5, alpha = 0.5) +  # Observed data
    geom_line(data = plot_data_list[[j]], 
              aes(x = date_vector, group = V, 
                  y = Fitted, color = V)) +  
    labs(x = "Date (2020)", y = "Cases") +
    scale_x_date(date_labels = "%b %d", date_breaks = "1 month") +  
    theme_bw() +
    scale_color_manual(name="",
                       labels=c("Sum","A", "B"),
                       values = alpha(values, 0.8)) +
    scale_fill_manual(name="",
                      labels=c("Sum","A", "B"),
                      values = alpha(values, 0.3)) +
    theme(legend.background = element_blank(),
          legend.position = c(0.82,0.72))
  print(p + ggtitle(j))
}
dev.off()

dfposterior = data.frame()
for (j in 1:9) {
  fit = fitlist[[id[j]]]
  posterior_samples = extract(fit)

  dfposterior[j,1] = j/10
  ratio = j/10
  
  m = mean(posterior_samples$beta1)
  q1 = quantile(posterior_samples$beta1, 0.025)
  q2 = quantile(posterior_samples$beta1, 0.975)
  dfposterior[j,2] = paste0(sprintf('%.3f',m),' (',
                            sprintf('%.3f',q1),' ~ ',
                            sprintf('%.3f',q2),')')
  
  m = mean(posterior_samples$beta2)
  q1 = quantile(posterior_samples$beta2, 0.025)
  q2 = quantile(posterior_samples$beta2, 0.975)
  dfposterior[j,3] = paste0(sprintf('%.3f',m),' (',
                            sprintf('%.3f',q1),' ~ ',
                            sprintf('%.3f',q2),')')
  
  m = mean(posterior_samples$gamma)
  q1 = quantile(posterior_samples$gamma, 0.025)
  q2 = quantile(posterior_samples$gamma, 0.975)
  dfposterior[j,4] = paste0(sprintf('%.3f',m),' (',
                            sprintf('%.3f',q1),' ~ ',
                            sprintf('%.3f',q2),')')
  
  c = posterior_samples$beta1/posterior_samples$beta2
  m = mean(c)
  q1 = quantile(c, 0.025)
  q2 = quantile(c, 0.975)
  dfposterior[j,5] = paste0(sprintf('%.3f',m),' (',
                            sprintf('%.3f',q1),' ~ ',
                            sprintf('%.3f',q2),')')
  
  c = -log(posterior_samples$beta1/posterior_samples$beta2)
  m = mean(c)
  q1 = quantile(c, 0.025)
  q2 = quantile(c, 0.975)
  dfposterior[j,6] = paste0(sprintf('%.3f',m),' (',
                            sprintf('%.3f',q1),' ~ ',
                            sprintf('%.3f',q2),')')
  
  I10 = round(34*ratio)
  I20 = round(34*(1-ratio))
  time_intro1 = log(I10)/(posterior_samples$beta1-posterior_samples$gamma) # - log(m)
  time_intro2 = log(I20)/(posterior_samples$beta2-posterior_samples$gamma) # - log(m)
  delta = time_intro1-time_intro2
  ci_lower = quantile(delta, probs = 0.025, na.rm = T)
  ci_upper = quantile(delta, probs = 0.975, na.rm = T)
  dfposterior[j,7] = paste0(sprintf('%.3f',mean(delta)),' (',
                            sprintf('%.3f',ci_lower),' ~ ',
                            sprintf('%.3f',ci_upper),')')

}

colnames(dfposterior) = c('Ratio','beta1','beta2','gamma',
                          'beta1/beta2','r2-r1', 'd2-d1')
write.csv(dfposterior, file = 'Output/F4D_acrosshost_parameters.csv',
          row.names = F)
