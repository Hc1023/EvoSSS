rm(list = ls())
library(data.table)
library(tidyverse)
library(readxl)
library(ggpubr)
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

# lineage A: 2, 8
# lineage B: 3, 9

idset = c('2','8','4','5','3','9')
df = df[df$ID %in% idset,]
idnew = c('A1','A2','B1','B2','B3','B4')
for (i in 1:6) {
  idx = (df$ID==idset[i])
  df$ID[idx] = idnew[i]
  df$group[idx] = strsplit(idnew[i], split = '')[[1]][1]
}
df$timepoint = factor(df$timepoint)
cell_lines = unique(df$cell_line)


cline = cell_lines[1]
df1 = df[df$cell_line == cline,]

get_anno <- function(df1){
  stat.test2 <- df1 %>%
    group_by(timepoint) %>%
    t_test(ORF1b ~ group) %>%
    add_significance()

  anno = data.frame(stat.test2)
  
  for (i in 1:nrow(anno)) {
    anno$text[i] = paste0('p=',
                          format(round(anno[i,'p'], 3), 
                                 nsmall = 3),' ',
                          anno[i,'p.signif'])

  }
  y = df1 %>% group_by(timepoint) %>%
    summarise_at(vars(ORF1b), function(c){
      ifelse(min(c) > 40, min(c) - 2, max(c) + 0.2)
    })
  anno$texty = y$ORF1b
  return(anno)
}

get_plot = function(df1){
  anno = get_anno(df1)
  
  values = c(hue_pal()(3)[1], hue_pal()(3)[3])
  p = ggplot(df1, aes(x=timepoint, y=ORF1b, color = group, fill = group)) +
    geom_boxplot(outlier.shape = NA) + 
    geom_jitter(position = position_jitterdodge(jitter.width = 0.1), 
                size=1, alpha=0.5, shape = 16,
                show.legend = F) +
    theme_bw() +
    scale_color_manual(values = values, name = 'Strain') +
    scale_fill_manual(values = alpha(values, 0.3), name = 'Strain') +
    labs(y = "ORF1b (Ct)", x = "Time points (h)") +
    theme(legend.background = element_rect(color = NA, fill = NA),
          legend.key = element_blank()) +
    annotate("text", x = anno$timepoint, y = anno$texty,
             label = anno$text, vjust = -0.5, size = 2) +
    coord_cartesian(ylim = c(min(df1$ORF1b), max(anno$texty)+1))
  
  return(p)
}


df1 = df[df$cell_line == cell_lines[1],]
p1 = get_plot(df1) + ggtitle('Huh-7')

df1 = df[df$cell_line == cell_lines[3],]
p2 = get_plot(df1) + ggtitle('Calu-3')

df1 = df[df$cell_line == cell_lines[4],]
df1 = df1[df1$timepoint != 1,]
p3 = get_plot(df1) + ggtitle('Vero')

pdf(file = paste0('Output/cd_replication.pdf'), width = 3.5, height = 2.2)
print(p1)
print(p2)
print(p3)
dev.off()



