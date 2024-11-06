rm(list = ls())
library(ggplot2)   
library(dplyr)
library(tidyr)
library(scales)
library(deSolve)   
library(rstan)

read_excel_allsheets <- function(filename, tibble = FALSE) {
  sheets <- readxl::excel_sheets(filename)
  x <- lapply(sheets, function(X) readxl::read_excel(filename, sheet = X))
  if(!tibble) x <- lapply(x, as.data.frame)
  names(x) <- sheets
  x
}

mysheets <- read_excel_allsheets("../2_Experiment/cd.xlsx")
df = mysheets[['replication']]
idset = c('2','8','4','5','3','9')
df = df[df$ID %in% idset,]
idnew = c('A1','A2','B1','B2','B3','B4')
for (i in 1:6) {
  idx = (df$ID==idset[i])
  df$ID[idx] = idnew[i]
  df$group[idx] = strsplit(idnew[i], split = '')[[1]][1]
}

cell_lines = unique(df$cell_line)

df1 = df[df$cell_line == cell_lines[3],]
df1$V = 2^(-df1$ORF1b + max(df1$ORF1b))
df1$Vlog = log(df1$V)
dat1 = df1[df1$timepoint != 72 & df1$group == 'A',]
dat2 = df1[df1$timepoint != 72 & df1$group == 'B',]

fit1 = lm(formula = Vlog ~ timepoint, data = dat1)
result1 = summary(fit1)
r1 = result1$coefficients[2,1]
fit2 = lm(formula = Vlog ~ timepoint, data = dat2)
result2 = summary(fit2)
r2 = result2$coefficients[2,1]
# > result1$coefficients[2,]
# Estimate   Std. Error      t value     Pr(>|t|) 
# 1.911578e-01 1.103223e-02 1.732722e+01 3.736368e-16 
# > result2$coefficients[2,]
# Estimate   Std. Error      t value     Pr(>|t|) 
# 2.272174e-01 6.407589e-03 3.546066e+01 5.932948e-42 

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
  dV2dt <- r2 * V2 * (1 - (V2 + alpha12 * V1) / K) + mu * V1
  
  return(list(c(dV1dt, dV2dt)))
}


state <- c(V1 = 1, V2 = 1)

# Time sequence for the simulation
times <- seq(0, 100, by = 1)
params = c(r1 = r1, 
           r2 = r2, 
           K = 2*median(df1[df1$timepoint == 72,'V']),  # 297235.6
           alpha12 = 0, 
           alpha21 = 0, 
           mu = 0)
out <- ode(y = state, times = times, func = viral_model, parms = params)
out_df <- as.data.frame(out)

# Plotting
pd = 2
values = c(hue_pal()(3)[1], hue_pal()(3)[3])

p1 = ggplot(data = out_df, aes(x = time)) +
  geom_line(aes(y = V1, color = "A")) +
  geom_line(aes(y = V2, color = "B")) +
  geom_boxplot(data = df1[df1$group == 'A',], 
               aes(x = timepoint-pd, y = V, 
                   group = timepoint),
               color = values[1], fill = alpha(values[1], 0.3),
               outlier.shape = NA, width = 1.9*pd) +
  geom_boxplot(data = df1[df1$group == 'B',], 
               aes(x = timepoint+pd, y = V, 
                   group = timepoint), 
               color = values[2], fill = alpha(values[2], 0.3),
               outlier.shape = NA, width = 1.9*pd) +
  geom_jitter(data = df1, 
              aes(x = timepoint, y = V, color = group, group = group),
              alpha = 0.5,
              position = position_jitterdodge(jitter.width = 1,
                                              dodge.width = 4*pd)) + 
  scale_color_manual(name = 'Variant', values = values) +
  scale_y_continuous(trans = 'log10') +
  scale_x_continuous(breaks = c(2,8,24,48,72),
                     minor_breaks = c()) +
  theme_bw() +
  labs(x = 'Time unit', y = 'Population unit') +
  theme(legend.position = 'none',
        legend.background = element_blank(),
        legend.key = element_blank())
p1
###### Vero #######

df1 = df[df$cell_line == cell_lines[4] & df$timepoint %in% c(2,8,24,48),]
df1$V = 2^(-df1$ORF1b + max(df1$ORF1b))
df1$Vlog = log(df1$V)
dat1 = df1[df1$timepoint != 48 & df1$group == 'A',]
dat2 = df1[df1$timepoint != 48 & df1$group == 'B',]
fit1 = lm(formula = Vlog ~ timepoint, data = dat1)
result1 = summary(fit1)
r1 = result1$coefficients[2,1]
fit2 = lm(formula = Vlog ~ timepoint, data = dat2)
result2 = summary(fit2)
r2 = result2$coefficients[2,1]

# > result1$coefficients[2,]
# Estimate   Std. Error      t value     Pr(>|t|) 
# 2.611061e-01 2.424562e-02 1.076921e+01 2.821663e-09 
# > result2$coefficients[2,]
# Estimate   Std. Error      t value     Pr(>|t|) 
# 3.397670e-01 2.035067e-02 1.669561e+01 4.344274e-19 

K = 2*median(df1[df1$timepoint == 48,'V']) # 47248.41
params = c(r1 = r1, 
           r2 = r2, 
           K =  K, 
           alpha12 = 0, 
           alpha21 = 0, 
           mu = 0)
out <- ode(y = state, times = times, func = viral_model, parms = params)
out_df <- as.data.frame(out)

p2 = ggplot(data = out_df, aes(x = time)) +
  geom_line(aes(y = V1, color = "A")) +
  geom_line(aes(y = V2, color = "B")) +
  geom_boxplot(data = df1[df1$group == 'A',], 
               aes(x = timepoint-pd, y = V, 
                   group = timepoint),
               color = values[1], fill = alpha(values[1], 0.3),
               outlier.shape = NA, width = 1.9*pd) +
  geom_boxplot(data = df1[df1$group == 'B',], 
               aes(x = timepoint+pd, y = V, 
                   group = timepoint), 
               color = values[2], fill = alpha(values[2], 0.3),
               outlier.shape = NA, width = 1.9*pd) +
  geom_jitter(data = df1, 
              aes(x = timepoint, y = V, color = group, group = group),
              alpha = 0.5,
              position = position_jitterdodge(jitter.width = 1,
                                              dodge.width = 4*pd)) + 
  scale_color_manual(name = 'Variant', values = values) +
  scale_y_continuous(trans = 'log10', 
                     breaks = c(10,1000,100000),
                     labels = c('1e+01','1e+03','1e+05')) +
  scale_x_continuous(breaks = c(2,8,24,48,72),
                     minor_breaks = c()) +
  theme_bw() +
  labs(x = 'Time unit', y = 'Population unit') +
  theme(legend.position = 'none',
        legend.background = element_blank(),
        legend.key = element_blank()) +
  coord_cartesian(ylim = c(1,100000))

pdf(paste0("Output/withinhost_par.pdf"), width = 2.5, height = 1.5)
print(p1)
print(p2)
dev.off()
