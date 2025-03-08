rm(list = ls())
library(data.table)
library(tidyverse)
library(ComplexHeatmap)
library(circlize)

df_all = read.csv("S5A_vcf_231228.csv")
sampleid = read.csv("sampleid.csv")
sampleid1 = sampleid[sampleid$original != '',c(1,2,3)]

for (i in 1:11) {
  df = df_all[df_all$sample== sampleid1[i,2],]
  df$MUT = paste0(df$REF, df$POS, df$ALT)
  df1 = df[, c('MUT','percent')]
  colnames(df1)[2] = sampleid1[i,1]
  if (i == 1){
    plot_data = df1
  }else{
    plot_data <- full_join(plot_data, df1, by = "MUT")
  }
}

plot_data[is.na(plot_data)] <- 0

mat <- t(plot_data[1:nrow(plot_data), 2:ncol(plot_data)])
class(mat) <- "numeric"
colnames(mat) <- plot_data[, 1]

n = sapply(colnames(mat), function(n){
  as.numeric(substr(n, 2, nchar(n)-1))
})


mat1 <- mat[, order(n)]

ht1 <- Heatmap(mat1, cluster_columns = FALSE, cluster_rows = FALSE,
              rect_gp = gpar(col = "grey", lwd = 0.5),
              column_names_gp = grid::gpar(fontsize = 10),
              col=colorRamp2(c(0, 1), c("white", "firebrick3")),
              heatmap_legend_param = list(
                title = "", at = c(0, 0.5, 1), 
                labels = c("0", "0.5", "1")
              ),
              cell_fun = function(j, i, x, y, width, height, fill){
                grid.text(if((colnames(mat1)[j] %in% c('C8782T','T28144C') &&
                              mat1[i,j] > 0||mat1[i,j] > 0.5)){
                  sprintf("%.3f", mat1[i,j])}, 
                  x, y, gp = gpar(fontsize = 10))
              },show_row_dend = FALSE,
              show_column_dend = FALSE)



for (i in 12:22) {
  df = df_all[df_all$sample== sampleid1[i,2],]
  df$MUT = paste0(df$REF, df$POS, df$ALT)
  df1 = df[, c('MUT','percent')]
  colnames(df1)[2] = sampleid1[i,1]
  if (i == 12){
    plot_data = df1
  }else{
    plot_data <- full_join(plot_data, df1, by = "MUT")
  }
}

plot_data[is.na(plot_data)] <- 0

mat <- t(plot_data[1:nrow(plot_data), 2:ncol(plot_data)])
class(mat) <- "numeric"
colnames(mat) <- plot_data[, 1]

n = sapply(colnames(mat), function(n){
  as.numeric(substr(n, 2, nchar(n)-1))
})
mat <- mat[, order(n)]

ht2 <- Heatmap(mat, cluster_columns = FALSE, cluster_rows = FALSE,
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

pdf(paste0("Output/S5A_vcf_231228.pdf"), width = 5, height = 2.6)
print(ht1)
print(ht2)
dev.off()

