rm(list = ls())
library(data.table)
library(tidyverse)
library(ComplexHeatmap)
library(circlize)
library(zoo)

getplot = function(x, cell){
  
  names = c('A1-2','A1-4','A1-5','B1-1','B2-1')
  df_all = read.csv("vcf_230428.csv")
  df = df_all[substr(df_all$sample, 1, 1) == x,]
  df$MUT = paste0(df$REF, df$POS, df$ALT)
  df1 = df[df$sample == paste0(x,1), c('MUT','percent')]
  colnames(df1)[2] = paste0(x,1)
  plot_data = df1
  df1 = df[df$sample == paste0(x,2), c('MUT','percent')]
  colnames(df1)[2] = paste0(x,2)
  plot_data <- full_join(plot_data, df1, by = "MUT")
  
  for (j in grep(x, names)) {
    df1 = df[df$sample == names[j], c('MUT','percent')]
    colnames(df1)[2] = names[j]
    plot_data <- full_join(plot_data, df1, by = "MUT")
  }
  plot_data[is.na(plot_data)] <- 0
  plot_data0 = data.frame('MUT' = plot_data[,1])
  plot_data0$Yao = rowMeans(plot_data[,2:3])
  plot_data0$Huang = rowMeans(plot_data[,4:ncol(plot_data)])
  
  df_all = read.csv("vcf_231228.csv")
  sampleid = read.csv("sampleid.csv")
  sampleid = sampleid[-grep('[+]',sampleid$X),]
  samples = sampleid[grep(x,sampleid$X),]
  samples = samples[samples$cell == cell, ]
  
  samples1 = data.frame(cell = cell,
                        id = c(samples[2*1:(nrow(samples)/2)-1,2],
                               samples$passaged1, 
                               samples$passaged2),
                        group = c(rep('Original',nrow(samples)/2),
                                  rep(c('Passage1','Passage2'),
                                      each = nrow(samples))))
  
  df = df_all[df_all$sample %in% samples[,10],]
  df2 = df %>% group_by(POS, REF, ALT) %>% summarise(mean = mean(percent))
  df2 = df1[df1$mean > 0.5,]
  for (i in 1:nrow(samples1)) {
    df = df_all[df_all$sample == samples1[i,2],]
    df$MUT = paste0(df$REF, df$POS, df$ALT)
    df1 = df[, c('MUT','percent')]
    colnames(df1)[2] = samples1[i,3]
    if (i == 1){
      plot_data = df1
    }else{
      plot_data <- full_join(plot_data, df1, by = "MUT")
    }
  }
  
  plot_data[is.na(plot_data)] <- 0
  
  plot_data1 = data.frame(MUT = plot_data[,1])
  plot_data1$Original = rowMeans(plot_data[,1 + 1:(nrow(samples)/2)])
  plot_data1$Passage1 = rowMeans(plot_data[,1 + nrow(samples)/2 + 1:nrow(samples)])
  plot_data1$Passage2 = rowMeans(plot_data[,1 + nrow(samples)*1.5 + 1:nrow(samples)])
  
  plot_data <- full_join(plot_data0, plot_data1, by = "MUT")
  plot_data[is.na(plot_data)] <- 0
  rows_to_drop <- apply(plot_data[, -1], 1, function(x) all(x < 0.05))
  plot_data <- plot_data[!rows_to_drop, ]
  
  colnames(plot_data) =c("MUT", "Yao et al.", "This study", "Passage1","Passage2", "Passage3" )
  
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
  return(ht)
}
p1 = getplot(x = 'A', cell = 'Calu')
p2 = getplot(x = 'A', cell = 'Vero')
p3 = getplot(x = 'B', cell = 'Calu')
p4 = getplot(x = 'B', cell = 'Vero')

pdf(paste0("Output/evolution_passage.pdf"), width = 5, height = 2.2)
print(p1)
print(p2)
print(p3)
print(p4)
dev.off()
