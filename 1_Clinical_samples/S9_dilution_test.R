rm(list = ls())
library(data.table)
library(tidyverse)
library(ComplexHeatmap)
library(circlize)

############ Amplicon Sequencing ########
load(file = 'S9_dilution_amplicon.Rdata')

plot_data <- lst.21.1.8[[1]]
for (i in 2:8) {
  plot_data <- full_join(plot_data, lst.21.1.8[[i]], by = "POS")
}
for (i in 1:8) {
  plot_data <- full_join(plot_data, lst.21.1.16[[i]], by = "POS")
}
plot_data[is.na(plot_data)] <- 0
pos.p <- union(plot_data[plot_data[2]!=0,1], plot_data[plot_data[10]!=0,1])

mat <- t(plot_data[1:nrow(plot_data), 2:ncol(plot_data)])
class(mat) <- "numeric"
colnames(mat) <- plot_data[, 1]
mat[-c(1,9), !(colnames(mat) %in% as.character(pos.p))] <- -mat[-c(1,9), !(colnames(mat) %in% as.character(pos.p))] 
mat <- mat[, order(as.numeric(colnames(mat)))]
rownames(mat) <- c("Amplicon-seq (1) Ct=12","1e3 X (1) Ct=22",
                   "1e4 X (1) Ct=25","1e5 X (1) Ct=29",
                   "1e6 X (1) Ct=32","1e7 X (1) Ct=35",
                   "1e8 X (1) Ct=39","1e9 X (1) Ct=42",
                   "Amplicon-seq (2) Ct=12","1e3 X (2) Ct=22",
                   "1e4 X (2) Ct=25","1e5 X (2) Ct=29",
                   "1e6 X (2) Ct=32","1e7 X (2) Ct=35",
                   "1e8 X (2) Ct=39","1e9 X (2) Ct=42")

mat0 <- mat[c(1,9),]
ht0<-Heatmap(mat0, cluster_columns = FALSE, cluster_rows = FALSE,
             rect_gp = gpar(col = "grey", lwd = 0.5),
             column_names_gp = grid::gpar(fontsize = 10),
             col=colorRamp2(c(-1, 0, 1), c("navy", "white", "firebrick3")),
             heatmap_legend_param = list(
               title = "", at = c(-1, -0.5, 0, 0.5, 1), 
               labels = c(1, 0.5, 0, 0.5, 1)
             ),
             cell_fun = function(j, i, x, y, width, height, fill){
               grid.text(if(abs(mat0[i, j])>0.5){sprintf("*")}, x, y, gp = gpar(fontsize = 20))
             },show_row_dend = FALSE,
             show_column_dend = FALSE)

mat1 <- mat[2:8,]
ht1<-Heatmap(mat1, cluster_columns = FALSE, cluster_rows = FALSE,
             rect_gp = gpar(col = "grey", lwd = 0.5),
             column_names_gp = grid::gpar(fontsize = 10),
             col=colorRamp2(c(-1, 0, 1), c("navy", "white", "firebrick3")),
             show_heatmap_legend = FALSE,
             cell_fun = function(j, i, x, y, width, height, fill){
               grid.text(if(abs(mat1[i, j])>0.5){sprintf("*")}, x, y, gp = gpar(fontsize = 20))
             },show_row_dend = FALSE,
             show_column_dend = FALSE)

mat2 <- mat[c(10:16),]
ht2<-Heatmap(mat2, cluster_columns = FALSE, cluster_rows = FALSE,
             rect_gp = gpar(col = "grey", lwd = 0.5),
             column_names_gp = grid::gpar(fontsize = 10),
             col=colorRamp2(c(-1, 0, 1), c("navy", "white", "firebrick3")),
             show_heatmap_legend = FALSE,
             cell_fun = function(j, i, x, y, width, height, fill){
               grid.text(if(abs(mat2[i, j])>0.5){sprintf("*")}, x, y, gp = gpar(fontsize = 20))
             },show_row_dend = FALSE,
             show_column_dend = FALSE)
ht_list =  ht0 %v% ht1 %v% ht2

pdf("Output/dilution_amplicon_S9A.pdf", width = 10, height = 5)
print(ht_list)
dev.off()


############ Capture and Metagenomic Sequencing ########
load(file = 'S9_dilution_capture.Rdata')

plot_data <- lst.20.9.24[[1]]
plot_data <- full_join(plot_data, lst.20.9.24[[2]], by = "POS")
for (i in 1:4) {
  plot_data <- full_join(plot_data, lst.20.10.30[[i]], by = "POS")
}
for (i in 1:3) {
  plot_data <- full_join(plot_data, lst.21.2.22[[i]], by = "POS")
}

plot_data <- lst.20.9.24[[1]]
plot_data <- full_join(plot_data, lst.20.9.24[[2]], by = "POS")
for (i in 1:4) {
  plot_data <- full_join(plot_data, lst.20.10.30[[i]], by = "POS")
}
for (i in 1:3) {
  plot_data <- full_join(plot_data, lst.21.2.22[[i]], by = "POS")
}


plot_data[is.na(plot_data)] <- 0
pos.p <- plot_data[plot_data[2]!=0,1]

mat <- t(plot_data[1:nrow(plot_data), 2:ncol(plot_data)])
class(mat) <- "numeric"
colnames(mat) <- plot_data[, 1]
mat[-1, !(colnames(mat) %in% as.character(pos.p))] <- -mat[-1, !(colnames(mat) %in% as.character(pos.p))] 
mat <- mat[, order(as.numeric(colnames(mat)))]

rownames(mat) <- c("Capture-seq Ct=12",
                   "Metagenomic-seq Ct=12",
                   "1e3 X (1) Ct=22",
                   "1e4 X (1) Ct=25",
                   "1e5 X (1) Ct=29",
                   "1e6 X (1) Ct=32",
                   "1e4 X (2) Ct=25",
                   "1e5 X (2) Ct=29",
                   "1e6 X (2) Ct=32")

mat0 <- mat[c(1,2),]
ht0<-Heatmap(mat0, cluster_columns = FALSE, cluster_rows = FALSE,
             rect_gp = gpar(col = "grey", lwd = 0.1),
             column_names_gp = grid::gpar(fontsize = 1),
             col=colorRamp2(c(-1, 0, 1), c("navy", "white", "firebrick3")),
             heatmap_legend_param = list(
               title = "", at = c(-1, -0.5, 0, 0.5, 1),
               labels = c(1, 0.5, 0, 0.5, 1)
             ),
             cell_fun = function(j, i, x, y, width, height, fill){
               grid.text(if(abs(mat0[i, j])>0.5){sprintf("*")}, x, y, gp = gpar(fontsize = 20))
             },show_row_dend = FALSE,
             show_column_dend = FALSE)

mat1 <- mat[3:6,]
ht1<-Heatmap(mat1, cluster_columns = FALSE, cluster_rows = FALSE,
             rect_gp = gpar(col = "grey", lwd = 0.1),
             column_names_gp = grid::gpar(fontsize = 1),
             col=colorRamp2(c(-1, 0, 1), c("navy", "white", "firebrick3")),
             show_heatmap_legend = FALSE,
             cell_fun = function(j, i, x, y, width, height, fill){
               grid.text(if(abs(mat1[i, j])>0.5){sprintf("*")}, x, y, gp = gpar(fontsize = 20))
             },show_row_dend = FALSE,
             show_column_dend = FALSE)

mat2 <- mat[c(7:9),]
ht2<-Heatmap(mat2, cluster_columns = FALSE, cluster_rows = FALSE,
             rect_gp = gpar(col = "grey", lwd = 0.1),
             column_names_gp = grid::gpar(fontsize = 1),
             col=colorRamp2(c(-1, 0, 1), c("navy", "white", "firebrick3")),
             show_heatmap_legend = FALSE,
             cell_fun = function(j, i, x, y, width, height, fill){
               grid.text(if(abs(mat2[i, j])>0.5){sprintf("*")}, x, y, gp = gpar(fontsize = 20))
             },show_row_dend = FALSE,
             show_column_dend = FALSE)

ht_list =  ht0 %v% ht1 %v% ht2

pdf("Output/dilution_capture_S9B.pdf", width = 10, height = 3)
print(ht_list)
dev.off()



