rm(list = ls())
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)

dat= read.csv("clinical_info.csv")
load(file = 'clinical_mutation.Rdata')
cd = readxl::read_excel('mut_cd.xlsx')
cd = data.frame(cd[1:46,c('SampleID','POS','REF','ALT','Freq')])

for (i in 1:nrow(dat)) {
  a <- dat[i,]
  b <- strsplit(a$Date, '/')[[1]]
  if(b[2] == 1) b[2] = 'Jan'
  if(b[2] == 2) b[2] = 'Feb'
  if(nchar(b[3]) == 1) b[3] = paste0('0',b[3])
  
  dat$tip[i] <- paste0(dat$Names_ID[i], ' ', b[2], ' ', b[3])
  dat$tip.lable[i] <- strsplit(a$sample_id, '-')[[1]][8]
}


dat_mut$SampleID = sapply(dat_mut$sample,  function(x){
  dat[dat$Sample == x, 'tip'][1]
})

df2 = dat_mut[,c(7,2:5)]
colnames(df2) = colnames(cd)
dat_all = rbind(cd, df2)

g = data.frame(labels = c('P04 Jan 29','ZJU_2', 'P11 Jan 25','ZJU_3',
                          'P15 Jan 29', 'ZJU_8', 'P20 Jan 29','ZJU_6',
                          'P52 Jan 25', 'ZJU_1', 'P54 Feb 03','ZJU_7',
                          'P17 Jan 25', 'ZJU_9', 'P55 Feb 03','ZJU_11',
                          'P56 Feb 03', 'ZJU_10'),
                labels2 = c("P04 Jan 29","ZJU-2 Jan 27", "P11 Jan 25", "ZJU-3 Jan 25",     
                            "P15 Jan 29", "ZJU-8 Jan 26", "P20 Jan 29", "ZJU-6 Feb 02",    
                            "P52 Jan 25", "ZJU-1 Jan 25", "P54 Feb 03", "ZJU-7 Feb 03",
                            "P17 Jan 25", "ZJU-9 Jan 28", "P55 Feb 03", "ZJU-11 Feb 04",    
                            "P56 Feb 03", "ZJU-10 Feb 03"))


df = data.frame()

for (i in 1:nrow(g)) {
  
  n = g$labels[i]
  df1 = dat_all[dat_all$SampleID == n,]
  df1 = data.frame(labels = rep(n,nrow(df1)),
                   Mutations = paste0(df1$REF,df1$POS,df1$ALT),
                   Freq = df1$Freq)

  df = rbind(df,df1)
}

df = spread(df, Mutations, Freq)
df[is.na(df)] = 0
rownames(df) = df[,1]
df = df[,-1]
n = sapply(colnames(df), function(n){
  as.numeric(substr(n, 2, nchar(n)-1))
})
df <- df[, order(n)]
df = df[g$labels,]
rownames(df) = g$labels2
df = as.matrix(df)

ht_opt$message = FALSE

g$Patients = rep(letters[1:9],each =2)
g$Source = c(rep('Sputum', 10), rep('Pharynx',2), rep(c('Sputum','Feces'), 3))

brewer.pal(n = 3, name = "RdBu")
display.brewer.pal(n = 6, name = 'Dark2')

v = alpha(brewer.pal(n = 9, name = 'Set1'), 0.9)
names(v) = letters[1:9]
v2 = alpha(brewer.pal(n = 6, name = 'Dark2'), 0.9)[c(1,3,6)]
names(v2) = c('Sputum', 'Pharynx', 'Feces')
col = list(
  Patients = v,
  Source = v2
)

p = Heatmap(df, cluster_columns = FALSE, cluster_rows = FALSE,
        rect_gp = gpar(col = "grey", lwd = 0.5),
        right_annotation = rowAnnotation(Source = g$Source,
                                         Patients = g$Patients,
                                         col = col),
        col=colorRamp2(c(0, 1), c("white", "firebrick3")),
        heatmap_legend_param = list(
          title = "", at = c(0, 0.5, 1), 
          labels = c("0", "0.5", "1")),
        show_row_dend = FALSE,
        show_column_dend = FALSE)
p
pdf("Output/feces.pdf", width = 12, height = 4)
print(p)
dev.off()
