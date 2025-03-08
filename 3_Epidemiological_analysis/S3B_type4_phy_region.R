rm(list = ls())
library(treeio)
library(ggtree)
library(ggplot2)
library(ape)
library(scales)
library(ggnewscale)

Ttree <- read.tree('S3B_tree.nwk')
# clock rate = 0.0008
Ttree$edge.length <- Ttree$edge.length/0.0008

Ttree.tip.label <- data.frame(Ttree$tip.label)

nextclade <- read.table(file = 'nextclade_qc.tsv', 
                        sep = '\t', header = TRUE, quote = "")

nextclade <- nextclade[nextclade$seqName %in% Ttree$tip.label, c(1,2,3,16)]

for (i in 1:nrow(nextclade)) {
  if(any(grep('C8782T', nextclade[i,4])) & any(grep('T28144C', nextclade[i,4]))){
    nextclade[i,5] <- 'Lineage A'
  }else if(any(grep('C8782T', nextclade[i,4]))){
    nextclade[i,5] <- '8782T'
  }else if(any(grep('T28144C', nextclade[i,4]))){
    nextclade[i,5] <- '28144C'
  }else{
    nextclade[i,5] <- 'Lineage B'
  }
}

# > table(nextclade$V5)
# 
# 28144C  8782T Lineage A Lineage B 
# 75       237       578       948 

cls <- list(Clade_1=nextclade$seqName[nextclade$V5 == 'Lineage A'],
            Clade_2=nextclade$seqName[nextclade$V5 == 'Lineage B'],
            Clade_3=nextclade$seqName[nextclade$V5 == '8782T'],
            Clade_4=nextclade$seqName[nextclade$V5 == '28144C'])

Ttree <- groupOTU(Ttree, cls)

values = c(hue_pal()(3)[1], hue_pal()(3)[3], hue_pal()(3)[2], hue_pal()(4)[4])

## GISAID continents

metadata <- read.table(file = 'S3B_strains4.tsv', 
                       sep = '\t', header = TRUE, quote = "")
metadata2 <- read.table(file = 'S3B_metadata_ref.tsv', 
                        sep = '\t', header = TRUE, quote = "")
for (i in 1:nrow(metadata)) {
  metadata[i,1] = paste(strsplit(metadata[i,1], '/')[[1]][-1], collapse = "/")
}

sample_dat = rbind(metadata[,c('strain', 'region')], metadata2[,c('strain', 'region')])
sample_dat = sample_dat[sample_dat$strain %in% Ttree$tip.label,]
table(sample_dat$region)

for (i in 1:nrow(sample_dat)) {
  sample_dat$Mutations[i] = nextclade$V5[nextclade$seqName %in% sample_dat$strain[i]]
}

sample_dat$region = factor(sample_dat$region,
                           levels = c('Africa', 'Asia',
                                      'Europe', 'North America', 
                                      'Oceania', 'South America'))
sample_dat$Mutations = factor(sample_dat$Mutations,
                              levels = c('Lineage A', 'Lineage B',
                                         '8782T', '28144C'))

values = c(hue_pal()(3)[1], hue_pal()(3)[3], hue_pal()(3)[2], hue_pal()(4)[4])

values2 = c(hue_pal()(6)[1:3], "#CB1414", hue_pal()(6)[5:6])

options(ignore.negative.edge=TRUE)
p <- ggtree(Ttree, mrsd = max(c(metadata$date,metadata2$date)), as.Date = T, 
            aes(color = group), size = 0.4) %<+% 
  sample_dat + theme_tree2() +
  scale_color_manual(name="Variant",
                     labels=c("Lineage A", "Lineage B", '8782 C>T', "28144 T>C"),
                     values = values) +
  new_scale_color() + 
  new_scale_fill() +
  geom_tippoint(mapping = aes(fill = region), size = 3, shape=21) +
  scale_fill_manual(name="Region", 
                    values=alpha(values2, alpha = 0.9)) + 
  xlab('Date (2019-2022)') +
  scale_x_date(date_labels = "%y %b", breaks = "3 month", minor_breaks = "1 month") +
  theme(axis.title=element_text(size=22),
        axis.text=element_text(size=18),
        axis.text.x = element_text(angle = 60, hjust = 1),
        legend.key.size = unit(1.2, 'cm'),
        legend.title = element_text(size=24),
        legend.text = element_text(size=22)
  ) +
  coord_cartesian(ylim = c(0, 1850)) +
  geom_text(data = data.frame(x = as.Date(c('2019-12-1', '2020-1-1', '2020-2-1',
                                            '2020-3-1', '2020-8-15',
                                            '2020-11-1', '2020-9-20',
                                            '2020-8-15', '2020-12-1', 
                                            '2020-3-10', '2019-12-1')), 
                              y = c(150, 240, 765, 
                                    745, 855,
                                    767, 740,
                                    1088, 1125,
                                    660, 20), 
                              text = c('19A', '19B', '20A', 
                                       '20B', '20I (Alpha, V1)',
                                       '21M (Omicron)', '20J (Gamma, V3)',
                                       '21A (Delta)', '21J (Delta)', 
                                       '20C', 'Hu-1')), 
            aes(x = x, y = y, label = text), 
            size = 5, alpha=.8,
            inherit.aes = FALSE) +
  guides(shape = guide_legend(order = 1), 
         fill = guide_legend(order = 2, override.aes = list(size=5)))

pdf(file = 'Output/S3B_phy_tree.pdf', width = 9, height = 13)
print(p)
dev.off()



