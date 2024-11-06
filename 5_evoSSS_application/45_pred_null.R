rm(list = ls())
library(rstan)
library(ggplot2)
library(scales)
library(ggnewscale)
library(tidyverse)
library(dplyr)
library(data.table)
library(lubridate)
library(ISOweek)
df1 = fread('VIW_FNT.csv')
df1[is.na(df1)] = 0

ob1 = df1[,c('ISO_WEEKSTARTDATE','INF_A','INF_B','INF_ALL')]
colnames(ob1)[1] = 'x'
ob1 = ob1 %>% 
  group_by(x) %>%
  summarise(yA = sum(INF_A),
            yB = sum(INF_B),
            yALL = sum(INF_ALL)) %>%
  as.data.frame()


if(F){
  dates = c(as.Date('2017-08-01'), as.Date('2021-08-01'))
  {
    voc = c('A','B')
    observed_matrix = ob1[ob1$x >= dates[1] & ob1$x < dates[2],-4] 
    dfob = data.frame(y = c(as.matrix(observed_matrix[,2:3])),
                      x = rep(as.Date(observed_matrix$x),2),
                      group = rep(voc, each = nrow(observed_matrix)))
    dfob$m = format(dfob$x, "%Y-%m")
    dfob1 = dfob %>% group_by(m, group) %>% 
      summarise(y = sum(y)) %>% as.data.frame()
  }
  
  data = observed_matrix[observed_matrix$x >= as.Date('2017-08-01') + months(3) &
                        observed_matrix$x < as.Date('2017-08-01') + months(39),]
  data$m = format(data$x, '%Y-%m')
  data = data %>% group_by(m) %>%
    summarise(yA = sum(yA), yB = sum(yB))
  data_pred = data[-nrow(data),-1]
  data = cbind(data[-1,], data_pred)
  data = rbind(data[,c(1,2,4)], setNames(data[,c(1,3,5)], names(data[,c(1,2,4)])))
  colnames(data) = c('date','x', 'y')
  data$group = rep(c('A','B'), each = nrow(data)/2)
  
  
  values = c('#98aff7','#a37c91')
  
  fit1 <- lm(log10(y) ~ log10(x), data = data[data$group == 'A',])
  slope1A <- coef(fit1)[2]
  intercept1A <- coef(fit1)[1]
  adj_r2_1A <- summary(fit1)$adj.r.squared
  
  fit1 <- lm(log10(y) ~ log10(x), data = data[data$group == 'B',])
  slope1B <- coef(fit1)[2]
  intercept1B <- coef(fit1)[1]
  adj_r2_1B <- summary(fit1)$adj.r.squared
  
  coordA = log10(c(min(data[data$group == 'A', c('x','y')]),
                   max(data[data$group == 'A', c('x','y')])))
  p1 = ggplot(data = data[data$group == 'A',], 
              aes(x = log10(y), y =log10(x))) +
    geom_point(color = alpha(values[1], 0.6)) +
    geom_abline(slope = slope1A, intercept = intercept1A, color = values[1]) +
    annotate("text", x = 3, y = 2.3,
             label = paste0("y = ", round(slope1A, 2), "x + ", round(intercept1A, 2),
                            "\nR² = ", round(adj_r2_1A, 3)),
             color = values[1], hjust = 0) +
    xlab(expression(log[10](Expected))) + ylab(expression(log[10](Predicted))) +
    theme_bw() +
    coord_cartesian(xlim = coordA, ylim = coordA) +
    scale_x_continuous(breaks = 2:5) + scale_y_continuous(breaks = 2:5) +
    theme(panel.grid.minor = element_blank(),
          legend.position = 'none')
  p1
  # y = 1.08x-0.22  r2 = 0.825
  coordB = log10(c(min(data[data$group == 'B',c('x','y')]),
                   max(data[data$group == 'B',c('x','y')])))
  p2 = ggplot(data = data[data$group == 'B',], 
              aes(x = log10(x), y =log10(y))) +
    geom_point(color = alpha(values[2], 0.6)) +
    geom_abline(slope = slope1B, intercept = intercept1B, color = values[2]) +
    annotate("text", x = 2.8, y = 2.1,
             label = paste0("y = ", round(slope1B, 2), "x + ", round(intercept1B, 2),
                            "\nR² = ", round(adj_r2_1B, 3)),
             color = values[2], hjust = 0) +
    xlab(expression(log[10](Expected))) + ylab(expression(log[10](Predicted))) +
    theme_bw() +
    scale_x_continuous(breaks = 2:5) + scale_y_continuous(breaks = 2:5) +
    coord_cartesian(xlim = coordB, ylim = coordB) +
    theme(panel.grid.minor = element_blank(),
          legend.position = 'none')
  p2
  # y=1.1x-0.25  r2=0.853
  pdf(paste0("Output/pred_month_null.pdf"), width = 2.2, height = 2)
  print(p1)
  print(p2)
  dev.off()
}

if(F){
  dates = c(as.Date('2017-08-01'),as.Date('2021-08-01'))
  w=1; w=2; w=3; w=4
  dataw = data.frame()
  for (w in 1:4) {
    {
      voc = c('A','B')
      observed_matrix = ob1[ob1$x >= dates[1] & ob1$x < dates[2],-4] 
      observed_matrix$adjx = as.Date(observed_matrix$x) - 7*(w-1)
      observed_matrix = observed_matrix[w:nrow(observed_matrix),]
      dfob = data.frame(y = c(as.matrix(observed_matrix[,2:3])),
                        x = rep(as.Date(observed_matrix$adjx),2),
                        group = rep(voc, each = nrow(observed_matrix)))
      dfob = dfob[as.numeric(format(dfob$x, '%d')) <= 7,]
      dfob$m = format(dfob$x, "%Y-%m")
      observed_matrix$x = format(observed_matrix$x, "%Y-%m")
      observed_matrix = observed_matrix %>% group_by(x) %>%
        summarise(yA=sum(yA), yB=sum(yB)) %>% 
        as.data.frame()
      rownames(observed_matrix) = as.Date(paste0(unlist(observed_matrix$x),'-01')) + months(1)-1
      observed_matrix = observed_matrix[,-1]
      
      colnames(observed_matrix) = voc
      fitlist = fitlistw[[w]]
    }
    
    
    df2_list = list()
    
    for (k in 3:39) {
      print(k)
      data = determinant_fun(pred = k)
      data2 = data[data$date >= as.Date('2017-08-01') + months(k) &
                     data$date < as.Date('2017-08-01') + months(k+1),]
      data2$m = format(data2$date, '%Y-%m')
      data2 = data2 %>% group_by(group, m) %>%
        summarise(y = sum(y)*7/30) %>% as.data.frame()
      data3 <- data2 %>%
        inner_join(dfob[,c('m','group','y')], by = c("m", "group"))
      data3$pred = k
      df2_list[[k-2]] = data3
    }
    
    data = data.frame(bind_rows(df2_list))
    colnames(data)[3:4] = c('y','x')
    data$w  = w
    dataw = rbind(dataw, data)
  }
  
  data = dataw
  values = c('#98aff7', '#a37c91')
  
  fit1 = lm(log10(y) ~ log10(x), data = data[data$group == 'A',])
  slope1A = coef(fit1)[2]
  intercept1A = coef(fit1)[1]
  adj_r2_1A = summary(fit1)$adj.r.squared
  
  fit1 = lm(log10(y) ~ log10(x), data = data[data$group == 'B',])
  slope1B = coef(fit1)[2]
  intercept1B = coef(fit1)[1]
  adj_r2_1B = summary(fit1)$adj.r.squared
  
  coordA = log10(c(min(data[data$group == 'A', c('x','y')]),
                   max(data[data$group == 'A', c('x','y')])))
  coordB = log10(c(min(data[data$group == 'B',c('x','y')]),
                   max(data[data$group == 'B',c('x','y')])))
  p1 = ggplot(data = data[data$group == 'A',], 
              aes(x = log10(y), y =log10(x))) +
    geom_point(color = alpha(values[1], 0.6)) +
    geom_abline(slope = slope1A, intercept = intercept1A, color = values[1]) +
    annotate("text", x = 0.9, y = 4.5,
             label = paste0("y = x - 0.23",
                            "\nR² = ", round(adj_r2_1A, 3)),
             color = values[1], hjust = 0) +
    # annotate("text", x = 3, y = 2.3,
    #          label = paste0("y = ", round(slope1A, 2), "x + ", round(intercept1A, 2),
    #                         "\nR² = ", round(adj_r2_1A, 3)),
    #          color = values[1], hjust = 0) +
    xlab(expression(log[10](Expected))) + ylab(expression(log[10](Predicted))) +
    theme_bw() +
    coord_cartesian(xlim = coordA, ylim = coordA) +
    scale_x_continuous(breaks = 2:5) + scale_y_continuous(breaks = 2:5) +
    theme(panel.grid.minor = element_blank(),
          legend.position = 'none')
  p1
  # y = x - 0.23; r2 = 0.877
  
  p2 = ggplot(data = data[data$group == 'B',], 
              aes(x = log10(x), y =log10(y))) +
    geom_point(color = alpha(values[2], 0.6)) +
    geom_abline(slope = slope1B, intercept = intercept1B, color = values[2]) +
    annotate("text", x = 0.8, y = 4.2,
             label = paste0("y = ", round(slope1B, 2), "x - 0.25", 
                            "\nR² = ", round(adj_r2_1B, 3)),
             color = values[2], hjust = 0) +
    # annotate("text", x = 3, y = 2.1,
    #          label = paste0("y = ", round(slope1B, 2), "x + ", round(intercept1B, 2),
    #                         "\nR² = ", round(adj_r2_1B, 3)),
    #          color = values[2], hjust = 0) +
    xlab(expression(log[10](Expected))) + ylab(expression(log[10](Predicted))) +
    theme_bw() +
    scale_x_continuous(breaks = 2:5) + scale_y_continuous(breaks = 2:5) +
    coord_cartesian(xlim = coordB, ylim = coordB) +
    theme(panel.grid.minor = element_blank(),
          legend.position = 'none')
  p2
  # y=1.02x-0.25; r2 = 0.895
  pdf(paste0("Output/pred_multiple_week.pdf"), width = 2.2, height = 2)
  print(p1)
  print(p2)
  dev.off()
}


