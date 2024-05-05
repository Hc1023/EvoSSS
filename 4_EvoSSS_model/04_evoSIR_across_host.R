rm(list = ls())
library(rstan)
library(ggplot2)
library(scales)
realData_all <- read.csv("Covid19CasesWH.csv", row.names = 1)
realData <- realData_all[-c(1:24), ] 
observed_cases = realData$CaseNum

# Define data
N = sum(observed_cases)
I0 = sum(realData_all$CaseNum[22:24])
I10 = I0*0.5
I20 = I0*0.5
S0 <- N - I0
ndays <- length(observed_cases)

load(file = 'SIR.rdata')
fit = fitlist[[5]]
posterior_samples <- extract(fit)

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
  
  pS_vec = c(beta1*I1/N, beta2*I2/N, 1-beta1*I1/N-beta2*I2/N)
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


simu_Onset1 <- matrix(NA, nrow = 20000, ncol = ndays)
simu_Onset2 <- matrix(NA, nrow = 20000, ncol = ndays)
# Simulate the epidemic for each set of sampled parameters
for (i in 1:20000) {

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

sum(df1$Fitted, na.rm = T)
sum(df2$Fitted, na.rm = T)
# Plot using ggplot2

values = c("#aa85a6", hue_pal()(3)[1], hue_pal()(3)[3])
show_col(values)
p = ggplot() +
  geom_ribbon(data = plot_data, 
              aes(x = date_vector, group = V, 
                  ymin = LowerCI, ymax = UpperCI, fill = V)) +  # Confidence interval
  geom_point(data = data_observed, 
             aes(x = date_vector, y = Observed), 
             size = 0.5, alpha = 0.5) +  # Observed data
  geom_line(data = plot_data, 
            aes(x = date_vector, group = V, 
                y = Fitted, color = V)) +  
  labs(x = "Date (2020)", y = "Cases") +
  scale_x_date(date_labels = "%b-%d", date_breaks = "1 month") +  
  theme_bw() +
  scale_color_manual(name="",
                     labels=c("A+B","A", "B"),
                     values = alpha(values, 0.8)) +
  scale_fill_manual(name="",
                    labels=c("A+B","A", "B"),
                    values = alpha(values, 0.3)) +
  theme(legend.background = element_blank(),
        legend.position = c(0.83,0.75))
p
pdf(paste0("Output/acrosshostWH.pdf"), width = 2.8, height = 1.8)
print(p)
dev.off()






