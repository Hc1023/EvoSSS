rm(list = ls())
library(data.table)
library(tidyverse)
library(ComplexHeatmap)
library(circlize)

df_all = read.csv("S5B_vcf_230428.csv")
names = c('ZJU_2', 'ZJU_8', 'ZJU_4', 'ZJU_5')
names2 = c('A1','A2','B1','B2')

pdf(paste0("Output/S5B_vcf_230428.pdf"), width = 5, height = 2.5)

for (i in 1:4) {
  ii = names2[i]
  df = df_all[substr(df_all$sample, 1, 2) == ii,]
  df$MUT = paste0(df$REF, df$POS, df$ALT)
  df1 = df[df$sample == ii, c('MUT','percent')]
  colnames(df1)[2] = ii
  plot_data = df1
  
  for (j in 1:5) {
    jj = paste(names2[i], j, sep = '-')
    df1 = df[df$sample == jj,c('MUT','percent')]
    colnames(df1)[2] = jj
    plot_data <- full_join(plot_data, df1, by = "MUT")
  }
  
  plot_data[is.na(plot_data)] <- 0
  
  mat <- t(plot_data[1:nrow(plot_data), 2:ncol(plot_data)])
  class(mat) <- "numeric"
  colnames(mat) <- plot_data[, 1]
  
  n = sapply(colnames(mat), function(n){
    as.numeric(substr(n, 2, nchar(n)-1))
  })
  mat <- mat[, order(n)]
  
  ht <- Heatmap(mat, cluster_columns = FALSE, cluster_rows = FALSE,
                   rect_gp = gpar(col = "grey", lwd = 0.5),
                   column_names_gp = grid::gpar(fontsize = 10),
                   col=colorRamp2(c(0, 1), c("white", "firebrick3")),
                   heatmap_legend_param = list(
                     title = "", at = c(0, 0.5, 1), 
                     labels = c("0", "0.5", "1")
                   ),
                   cell_fun = function(j, i, x, y, width, height, fill){
                     grid.text(if((colnames(mat)[j] %in% c('C8782T','T28144C') &&
                                  mat[i,j] > 0||mat[i,j] > 0.5)){
                       sprintf("%.3f", mat[i,j])}, 
                       x, y, gp = gpar(fontsize = 10))
                   },show_row_dend = FALSE,
                   show_column_dend = FALSE)
  print(ht)
}
dev.off()
