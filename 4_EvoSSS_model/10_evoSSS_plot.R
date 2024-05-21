
for (n in 1:24) {
  Onsets_mat = Onsets_mat_list[[n]]
  Onset1 = Onsets_mat[poolday + 1:poolday, 1]
  Onset2 = Onsets_mat[poolday + 1:poolday, 2]
  
  mobility = rep(0.01,30) # Mobility: Control force
  Mobility_matrix = diag(mobility)
  
  seed_vec =  (Onset1 + Onset2) %*% Mobility_matrix %>% as.numeric()
  p = Onset1/(Onset1 + Onset2)
  
  seed_matrix = rbind(seed_vec * p, seed_vec * (1-p))
  seed_mat_I1 = diag(seed_matrix[1,])
  seed_mat_I2 = diag(seed_matrix[2,])
  fit = fitlist[[n]]
  
  posterior = rstan::extract(fit)
  contact = mean(posterior$contact)
  N = seed_vec * contact + 1
  Onsets_mat = simu(seed_mat_I1, seed_mat_I2, N, poolday, pars)
  Onsets_mat_list[[n+1]] = Onsets_mat
}


# contact_vec_init = c(10000,100,100,150,280,260,260,260,260,260,260,
#                      260,260,260,150,160,200,400,300,200,200,260)
# 


dfobserve = data.frame(
  x = rep(1:700,2),
  y = c(observed_matrix[1:700,1], observed_matrix[1:700,2])/28,
  group = c(rep('A', 700),rep('B', 700))
)

{
  dfplot_simu = data.frame()
  for (i in 1:25) {
    Onsets_mat = Onsets_mat_list[[i]]
    n = i-1
    
    dfplot_simu1 = data.frame(x = rep(1:nrow(Onsets_mat)+poolday*n,2),
                              y = c(Onsets_mat[,1],Onsets_mat[,2]),
                              group = rep(paste0(c('A','B'),n), 
                                          each = nrow(Onsets_mat)),
                              color = rep(c('A','B'), 
                                          each = nrow(Onsets_mat)))
    dfplot_simu = rbind(dfplot_simu, dfplot_simu1)
  }
  
  df2 = dfplot_simu %>% group_by(x, color) %>% 
    summarise(y = sum(y)) %>%
    as.data.frame()
  values = c(hue_pal()(3)[1], hue_pal()(3)[3])
  dfplot_simu$y = dfplot_simu$y/28
  df2$y = df2$y/28
  dfobserve$Date = dfobserve$x + as.Date('2019-12-31')
  dfplot_simu$Date = dfplot_simu$x + as.Date('2019-12-31')
  df2$Date = df2$x + as.Date('2019-12-31')
  
  p = ggplot() +
    geom_line(data = dfplot_simu, 
              aes(x = Date, y = y, 
                  group = group, color = color), alpha = 0.2) + 
    geom_point(data = dfobserve, 
               aes(x = Date, y = y, 
                   group = group, color = group),
               size = 0.4, shape = 16) +
    geom_line(data = df2, 
              aes(x = Date, y = y, 
                  group = color, color = color), linewidth = 1) + 
    scale_y_continuous(trans=scales::pseudo_log_trans(base=2),
                       breaks = c(0, 2^seq(2,20,3)),
                       labels = c(0, expression('2'^'2'),expression('2'^'5'),
                                  expression('2'^'8'),expression('2'^'11'),
                                  expression('2'^'14'), expression('2'^'17'),
                                  expression('2'^'20'))) + 
    theme_bw() +
    scale_color_manual(name="Variant",
                       labels=c("A", "B"),
                       values = alpha(values, 0.6)) +
    scale_fill_manual(name="Variant",
                      labels=c("A", "B"),
                      values = alpha(values, 0.6)) +
    scale_x_date(breaks = seq(as.Date('2020-01-01'), as.Date('2022-11-01'), by="6 months"),
                 minor_breaks = seq(as.Date('2019-12-01'), as.Date('2022-11-01'), by ='1 month'),
                 date_labels = "%y-%b") +
    xlab('Date (2019-2021)') + ylab('') + 
    coord_cartesian(xlim = c(as.Date('2019-12-31'), as.Date('2021-10-05'))) +
    theme(legend.position = 'none')
  
  pdf(paste0("Output/evoSSS_stan.pdf"), width = 3, height = 2.3)
  print(p)
  dev.off()
}
