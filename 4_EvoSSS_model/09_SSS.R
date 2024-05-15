rm(list = ls())
library(rstan)
library(ggplot2)
library(scales)
library(ggplot2)
library(scales)
library(ggnewscale)

realData_all <- read.csv("Covid19CasesWH.csv", row.names = 1)
realData <- realData_all[-c(1:24), ]
observed_cases = realData$CaseNum

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
  
  return(data.frame(S_new, I1_new, I2_new, Onset1, Onset2))
}

simu <- function(pars, seed_mat_list, dim, f = update_fun) {
  # Initial conditions
  N = rep(32583, dim)
  seed_mat_I1 = seed_mat_list[[1]]
  seed_mat_I2 = seed_mat_list[[2]]
  I1_old = seed_mat_I1[1,]
  I2_old = seed_mat_I2[1,]
  S_old = N - I1_old - I2_old
  states_old = data.frame(S_old, I1_old, I2_old)

  # Simulate the dynamics over ndays
  ndays = dim*2
  mycol <- c("time", "Onset1", "Onset2")
  Onsets_mat <- matrix(0, ndays, length(mycol))
  Onsets_mat[,1] = 1:ndays
  
  colnames(Onsets_mat) <- mycol
  pars = c(0.379, 0.398, 0.157)
  for (t in 1:(ndays-1)) {
    states_mat = f(pars = pars, states_old = states_old)
    states_old = states_mat[1:3]
    states_old[,2:3] = states_old[,2:3] + cbind(seed_mat_I1[t+1,], seed_mat_I2[t+1,])
    
    Onsets_mat[t+1,-1] = c(sum(states_mat[4]), sum(states_mat[5]))
  }
  
  return(Onsets_mat)
}

dim = 100
spacing = round(dim/2)
mobility = rep(0, dim)
mobility[30:31] = 1 
Mobility_matrix = diag(mobility)

I_last = matrix(rep(20, dim), nrow = dim)
seed_vec = Mobility_matrix %*% I_last
p = rep(0.2, dim)

seeding_time = spacing + 1:dim
Onset1 = Onsets_mat[seeding_time,2]
Onset2 = Onsets_mat[seeding_time,3]
p = Onset1/(Onset1 + Onset2)
seed_matrix = rbind(t(seed_vec * p), t(seed_vec * (1-p)))
seed_matrix = sapply(1:dim, function(x){
  rmultinom(1, seed_vec[x], c(p[x], 1-p[x]))
})

# The initial cycle - epidemic outbreak
seed_matrix = rbind(rep(0,dim),rep(0,dim))
seed_matrix[,1] = round(34*c(0.4,0.6))

seed_mat_list = list()
seed_mat_list[[1]] = rbind(diag(seed_matrix[1,]), matrix(0, dim, dim))
seed_mat_list[[2]] = rbind(diag(seed_matrix[2,]), matrix(0, dim, dim))
Onsets_mat = simu(pars, seed_mat_list, dim)
dfplot_simu = data.frame(x = rep(Onsets_mat[,1],2),
                 y = c(Onsets_mat[,2],Onsets_mat[,3]),
                 group = rep(c('A0','B0'), each = nrow(Onsets_mat)),
                 color = rep(c('A','B'), each = nrow(Onsets_mat)))
# n = 1
spacing = 100
seeding_time = spacing*(n-1) + 1:dim
Onset1 = Onsets_mat[seeding_time,2]
Onset2 = Onsets_mat[seeding_time,3]
I_last = Onset1 + Onset2
mobility = rep(0, dim)
mobility[50:51] = 1 
Mobility_matrix = diag(mobility)
seed_vec =  I_last %*% Mobility_matrix
p = Onset1/(Onset1 + Onset2)
seed_matrix = rbind(seed_vec * p, seed_vec * (1-p))
seed_matrix[,1] = 0
seed_mat_list = list()
seed_mat_list[[1]] = rbind(diag(seed_matrix[1,]), matrix(0, dim, dim))
seed_mat_list[[2]] = rbind(diag(seed_matrix[2,]), matrix(0, dim, dim))
Onsets_mat = simu(pars, seed_mat_list, dim)
dfplot_simu1 = data.frame(x = rep(seeding_time[1]-1+Onsets_mat[,1],2),
                         y = c(Onsets_mat[,2],Onsets_mat[,3]),
                         group = rep(paste0(c('A','B'),n), each = nrow(Onsets_mat)),
                         color = rep(c('A','B'), each = nrow(Onsets_mat)))
dfplot_simu = rbind(dfplot_simu, dfplot_simu1)
df1 = dfplot_simu
df1$x = as.Date('2019-12-31') + df1$x

df = read.csv('../3_Epidemiological_analysis/Covid19CasesGISAID.csv')

df$Var1 = as.Date(df$Var1)
df = df[df$Var1 < as.Date('2021-12-01'),]
df = df[df$Mutations %in% c('Lineage A', 'Lineage B'),]

document_scale = max(df1[,2])/max(df$Freq[df$Var1<as.Date('2020-02-01')])
df1[,2] = df1[,2]/document_scale

values = c(hue_pal()(3)[1], hue_pal()(3)[3])



ggplot() + 
  geom_point(data = df, 
             aes(x = Var1, y = Freq, colour = Mutations),
             size = 0.5) +
  scale_color_manual(name="Variant",
                     labels=c("A", "B"),
                     values = alpha(values, 0.6)) +
  scale_y_continuous(trans=scales::pseudo_log_trans(base=2),
                     breaks = c(0, 2^seq(0,14,2))) +
  new_scale_color() +
  new_scale_fill() +
  geom_line(data = df1, 
            aes(x, y, color = color, 
                group = group), linewidth = 1) + 
  scale_color_manual(name="Variant",
                     labels=c("A", "B"),
                     values = alpha(values, 0.6)) +
  scale_fill_manual(name="Variant",
                    labels=c("A", "B"),
                    values = alpha(values, 0.6)) +
  scale_x_date(breaks = seq(as.Date('2020-01-01'), as.Date('2022-11-01'), by="6 months"),
               minor_breaks = seq(as.Date('2019-12-01'), as.Date('2022-11-01'), by ='1 month'),
               date_labels = "%y-%b") +
  xlab('Date (2019-2021)') +
  ylab('') + theme_bw() 
