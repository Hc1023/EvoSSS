rm(list = ls())
library(dplyr)
library(tidyverse)
library(treeio)
library(ggtree)
library(ggtreeExtra)
library(ggnewscale)
library(stringr)

ddf = read.csv('sarsrcov.csv')
id = ddf[,c(1,22)]

seq = unlist(strsplit('TTTAGCCAG', split = ""))
loc = 8777:8785
seqloc1 = paste0(seq, loc)
seq = unlist(strsplit('TGTTTACCT', split = ""))
loc = 28140:28148
seqloc2 = paste0(seq, loc)

df1 <- ddf[,c(22,4:12)]
df1 <-  df1 %>% gather(mut,base,-Name)
for (i in 1:nrow(df1)) {
  df1$pos[i] <- which(seqloc1 == df1[i,2])
}

df1$col = 'blue'
df1$col[df1$mut %in% seqloc1[4:6]] <- 'red'

df2 <- ddf[,c(22,13:21)]
df2 <-  df2 %>% gather(mut,base,-Name)
for (i in 1:nrow(df2)) {
  df2$pos[i] <- which(seqloc2 == df2[i,2])
}

df2$col = 'blue'
df2$col[df2$mut %in% seqloc2[4:6]] <- 'red'


Ttree <- read.tree('sarsrcov_run_3_tr.phy.treefile')
# clock rate = 0.0008
tiplabel = Ttree$tip.label
idx <- sapply(tiplabel, function(x){
  strsplit(x, '[hNSAOGBDM]')[[1]][1]
})
idx <- unname(idx) %>% as.numeric()

Ttree$tip.label <- ddf$Name[idx+1]


df2$text_color = 'black'
df2$text_color[df2$base == '-'] = 'red'
sample_dat = data.frame(lable = Ttree$tip.label,
                        group = c(rep('0',3), rep('1',17), rep('0',37)))
p <- ggtree(Ttree) %<+% 
  sample_dat + 
  geom_tiplab(aes(color = group), 
              offset = 0.1, size = 3, geom = "text", align = F) +
  scale_color_manual(values = c('black','#c23a31'),
                     guide = 'none')
p
p1 <- p + geom_text(aes(label=node), size = 2)
p1
rotate(p1,72)
# flip(p1, 54, 60) %>% flip(8,9) %>% flip(61,83)

p1 = rotate(p,72)
p1
hexp1 = 6.2
hexp2 = 10
pw = 1
p2 = p1 + 
  new_scale_fill() + 
    geom_fruit(
      data=df1,
      geom=geom_tile,
      mapping=aes(y=Name, x=pos, fill=col),
      position_identityx(hexpand = hexp1, vexpand = NA),
      color = 'white',
      size = 0.03,
      offset = 0.08,   # default 0.03
      pwidth = pw
    ) +
  scale_fill_manual(
    values=alpha(c("#C2D7EC", "#F1CDB2"), 0.9),
    guide="none") +
  new_scale_fill() +
  geom_fruit(
    data=df1,
    geom=geom_text,
    mapping=aes(y=Name, x=pos, label=base), color = 'black',
    position_identityx(hexpand = hexp1, vexpand = NA),
    size = 3,
    offset = 0.08,   # default 0.03
    pwidth = pw
  ) +
  new_scale_fill() +
  geom_fruit(
    data=df2,
    geom=geom_tile,
    mapping=aes(y=Name, x=pos, fill=col),
    position_identityx(hexpand = hexp2, vexpand = NA),
    color = 'white',
    size = 0.03,
    offset = 0.08,   # default 0.03
    pwidth = pw
  ) +
  scale_fill_manual(
    values=alpha(c("#C2D7EC", "#F1CDB2"), 0.9),
    guide="none") +
  new_scale_fill() +
  new_scale_color() +
  geom_fruit(
    data=df2,
    geom=geom_text,
    mapping=aes(y=Name, x=pos, label=base, color = text_color),# color = 'black',
    position_identityx(hexpand = hexp2, vexpand = NA),
    size = 3,
    offset = 0.08,   # default 0.03
    pwidth = pw
  ) +
  scale_color_manual(
    values=alpha(c("black", "red"), 1),
    guide="none") 

p2



pdf(file = 'Output/SARSrCoV.pdf', width = 8, height = 8)
print(p2)
dev.off()

