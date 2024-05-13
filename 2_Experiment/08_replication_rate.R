rm(list = ls())
library(ggplot2)   
library(dplyr)
library(tidyr)
library(scales)

read_excel_allsheets <- function(filename, tibble = FALSE) {
  sheets <- readxl::excel_sheets(filename)
  x <- lapply(sheets, function(X) readxl::read_excel(filename, sheet = X))
  if(!tibble) x <- lapply(x, as.data.frame)
  names(x) <- sheets
  x
}

mysheets <- read_excel_allsheets("cd.xlsx")
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
df$Vlog = -df$ORF1b*log(2)
dfsummary = data.frame(expand.grid(group = c('A','B'),
                                   cell = cell_lines[3:4]))
dfplot = data.frame()
for (i in 1:nrow(dfsummary)) {
  df_filter = df[df$cell_line == dfsummary[i,2] & 
                 df$group == dfsummary[i,1], ]
  times = unique(df_filter$timepoint)
  df2 = sapply(2:length(times), function(t){
    df1 = df_filter[df_filter$timepoint %in% c(times[t-1],times[t]),]
    model <- lm(Vlog ~ timepoint, data = df1)
    s = summary(model)
    return(c(times[t],s$coefficients[2,c(1,2)]))
  }) %>% t() %>% as.data.frame()
  
  K = df_filter$Vlog[df_filter$timepoint == max(times)]-df_filter$Vlog[df_filter$timepoint == min(times)]
  df2[length(times),] = c(-1,mean(exp(K)),sd(exp(K)))
  df2$cell_line = dfsummary[i,2]
  df2$group = dfsummary[i,1]
  dfplot = rbind(dfplot, df2)
}
dfplot[dfplot$V1 == -1,-1]
# Capacity
#   Estimate Std. Error cell_line group
#   78156.831 56719.3912     calu3     A
#   44090.609 17868.7685     calu3     B
#   1160.215   340.3424      vero     A
#   6996.790  4928.9688      vero     B
dfplot1 = dfplot[!dfplot$V1 %in% c(-1,2),]


colnames(dfplot1) = c('x','y','sd','cell','group')

dfplot1$x = factor(dfplot1$x, levels = sort(unique(dfplot1$x)))
values = c(hue_pal()(3)[1], hue_pal()(3)[3])

p1 = ggplot(filter(dfplot1, cell == 'calu3'), aes(x, y)) +
  geom_errorbar(
    aes(ymin = y-sd, ymax = y+sd, color = group),
    position = position_dodge(0.6), width = 0.5
  )+
  geom_point(aes(color = group), position = position_dodge(0.6),
             shape = 19, size = 2) +
  scale_color_manual(name = 'Variants', values = alpha(values, 0.7)) +
  theme_bw() +
  xlab('Time') + ylab('Replication rate')

p2 = ggplot(filter(dfplot1, cell == 'vero'), aes(x, y)) +
  geom_errorbar(
    aes(ymin = y-sd, ymax = y+sd, color = group),
    position = position_dodge(0.6), width = 0.5
  )+
  geom_point(aes(color = group), position = position_dodge(0.6),
             shape = 19, size = 2) +
  scale_color_manual(name = 'Variants', values = alpha(values, 0.7)) +
  theme_bw() +
  xlab('Time') + ylab('Replication rate')



df = read.csv("calu3_ct.csv")
get_df_long = function(df1){
  colnames(df1) = c('Strain', '0', '8', '24', '48', '72')
  
  df2 = df1[-c(1:10), ] %>%
    group_by(Strain) %>%
    summarise_at(vars(-group_cols()), mean, na.rm = TRUE) %>%
    as.data.frame()
  
  df2 = rbind(df1[1:10, ], df2)
  df2$Strain <- factor(c(rep('A', 6), rep('B', 4), rep('A+B', 6)),
                       levels = c('A', 'B', 'A+B'))
  
  df_long <- gather(df2, key = "Time",
                    value = "Ct",-Strain)
  df_long$Time = as.numeric(df_long$Time)
  return(df_long)
}
df1 = df[, 1:6]
df = get_df_long(df1)

df$Vlog = -df$Ct*log(2)
v =  levels(df$Strain)

dfplot = data.frame()
for (i in 1:length(v)) {
  df_filter = df[df$Strain == v[i], ]
  times = unique(df$Time)
  df2 = sapply(2:length(times), function(t){
    df1 = df_filter[df_filter$Time %in% c(times[t-1],times[t]),]
    model <- lm(Vlog ~ Time, data = df1)
    s = summary(model)
    return(c(times[t],s$coefficients[2,c(1,2)]))
  }) %>% t() %>% as.data.frame()
  
  K = df_filter$Vlog[df_filter$Time == max(times)]-df_filter$Vlog[df_filter$Time == min(times)]
  df2[length(times),] = c(-1,mean(exp(K)),sd(exp(K)))
  df2$group = v[i]
  dfplot = rbind(dfplot, df2)
}
dfplot[dfplot$V1 == -1,-1]

# Capacity
# Estimate    Std. Error group
# 926766.3    319103.0     A
# 1519375.0   656484.2     B
# 726761.3    141975.5   A+B

dfplot1 = dfplot[!dfplot$V1 %in% c(-1),]
dfplot1$group = factor(dfplot1$group, levels = v)

colnames(dfplot1) = c('x','y','sd','group')

dfplot1$x = factor(dfplot1$x, levels = sort(unique(dfplot1$x)))
values = c(hue_pal()(3)[1], hue_pal()(3)[3], '#aa85a6')

p0 = ggplot(filter(dfplot1), aes(x, y)) +
  geom_errorbar(
    aes(ymin = y-sd, ymax = y+sd, color = group),
    position = position_dodge(0.6), width = 0.5
  )+
  geom_point(aes(color = group), position = position_dodge(0.6),
             shape = 19, size = 2) +
  scale_color_manual(name = 'Variants', values = alpha(values, 0.7)) +
  theme_bw() +
  xlab('Time') + ylab('Replication rate')
p0

pdf(paste0("Output/replication_rate.pdf"), width = 3, height = 1.8)
print(p0)
print(p1)
print(p2)
dev.off()
