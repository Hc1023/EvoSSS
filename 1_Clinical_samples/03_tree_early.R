rm(list = ls())
library(treeio)
library(ggtree)
library(ggplot2)
library(ape)
library(scales)
library(ggnewscale)
library(tidyverse)

Ttree <- read.tree('tree_early.nwk')
# clock rate = 0.0008
Ttree$edge.length <- Ttree$edge.length/0.0008

nextclade <- read.table(file = 'nextclade_qc_early.tsv', 
                        sep = '\t', header = TRUE)
nextclade <- nextclade[nextclade$seqName %in% Ttree$tip.label, c(1,2,3,16)]

for (i in 1:nrow(nextclade)) {
  if(any(grep('G26144T', nextclade[i,4])) & 
     any(grep('G11083T', nextclade[i,4]))){
    nextclade[i,'cls'] <- 'V'
  }else{
    nextclade[i,'cls'] <- strsplit(nextclade[i,3], split = "[.]")[[1]][1]
  }
}


table(nextclade[,'cls']) 

# V clade all in lineage B
# Haplotype: artifacts

Ttree.drop <- drop.tip(phy = Ttree, tip = c('Wuhan/WH01/2019'))
cls <- list(Clade_1=nextclade$seqName[nextclade$cls == 'A'],
            Clade_2=nextclade$seqName[nextclade$cls == 'B'],
            Clade_3=nextclade$seqName[nextclade$cls == 'V'])
Ttree.drop <- groupOTU(Ttree.drop, cls)

## GISAID continents
metadata <- read.table(file = 'custom.metadata.tsv', 
                       sep = '\t', header = TRUE, quote = "")
metadata <- metadata %>% distinct(strain, .keep_all = TRUE)

j = 149
regions = c()
nrow(nextclade)
for (i in (j+1):nrow(nextclade)) {
  regions = c(regions,
              metadata$region[metadata$strain %in% paste('hCoV-19', nextclade$seqName[i], 
                                                         sep = '/')])
}

Names_ID = c(rep('This study', j), regions)
sample_dat <- data.frame(seqName = nextclade$seqName, 
                         Names_ID = factor(Names_ID,
                                           levels = c('Africa', 'Asia',
                                                      'Europe', 'North America', 'Oceania',
                                                      'South America',
                                                      'This study')))



values = c(hue_pal()(3)[1], hue_pal()(3)[3], '#5fbfb1')
values2 = c(hue_pal()(6)[1:3], "#CB1414", hue_pal()(6)[5:6], "#F0F014")

options(ignore.negative.edge=TRUE)

sizetext = 4

p1 <- ggtree(Ttree.drop, mrsd = '2020-03-01', as.Date = T, 
             aes(color = group), size = 0.3) %<+% 
  sample_dat + theme_tree2() +
  scale_color_manual(name="Clades",
                     labels=c("Lineage A", "Lineage B",
                              "Clade V"),
                     values = values) +
  new_scale_color() + 
  new_scale_fill() +
  geom_tippoint(aes(fill = Names_ID, 
                    subset = !(label %in% nextclade$seqName[1:149])), 
                size = 3, shape=21) +
  geom_tippoint(aes(fill = Names_ID, 
                    subset = label %in% nextclade$seqName[1:149]), 
                size = 4, shape=21, show.legend = F) +
  scale_color_manual(name="Samples", 
                     values=alpha(values2, alpha = 1)) + 
  scale_fill_manual(name="Samples", 
                    values=alpha(values2, alpha = 1)) + 
  scale_x_date(date_labels = "%y %b", breaks = "1 month") +
  coord_cartesian(ylim = c(-10, 2100)) +
  theme(legend.key.size = unit(0.4, 'cm'),
        legend.background = element_blank()) +
  guides(shape = guide_legend(order = 2),
         fill = guide_legend(order = 1))
p1



pdf(file = 'Output/tree_early.pdf', width = 4.5, height = 3)
print(p1)
dev.off()


