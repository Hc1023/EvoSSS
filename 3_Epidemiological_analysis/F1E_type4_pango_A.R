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

nextclade <- read.table(file = 'S3B_nextclade_qc.tsv', 
                        sep = '\t', header = TRUE, quote = "")

nextclade <- nextclade[nextclade$seqName %in% Ttree$tip.label, c(1,2,3,16)]

for (i in 1:nrow(nextclade)) {
  if(any(grep('C8782T', nextclade[i,4])) | any(grep('T28144C', nextclade[i,4]))){
    nextclade[i,5] <- 'A'
  }else{
    nextclade[i,5] <- 'B'
  }
}

# > table(nextclade$V5)
# 
# 28144C  8782T Lineage A Lineage B 
# 75       237       578       948 

cls <- list(Clade_1=nextclade$seqName[nextclade$V5 == 'A'],
            Clade_2=nextclade$seqName[nextclade$V5 == 'B'])

Ttree <- groupOTU(Ttree, cls)

values = c(hue_pal()(3)[1], hue_pal()(3)[3])

load('F1E_pango_A.Rdata')
values2 = values2[names(values2) %in% unique(sample_dat$lineage)]

sample_dat$lineage = factor(sample_dat$lineage,
                            levels = names(values2))

metadata <- read.table(file = 'S3B_strains4.tsv', 
                       sep = '\t', header = TRUE, quote = "")
metadata2 <- read.table(file = 'S3B_metadata_ref.tsv', 
                        sep = '\t', header = TRUE, quote = "")

options(ignore.negative.edge=TRUE)

x = ggtree(Ttree) %>% collapse(node=2306)  
x
y = ggtree(Ttree)
y
y$data$angle = x$data$angle 
y
z = ggtree(Ttree, as.Date = T)
Ttree.drop <- drop.tip(phy = Ttree,
                       tip = Ttree$tip[!Ttree$tip %in% x$data$label])

table(Ttree$tip %in% x[["data"]][["label"]])

maxdate = as.Date(max(c(metadata$date,metadata2$date)))

p <- ggtree(Ttree, #mrsd = maxdate, 
            as.Date = T, 
            aes(color = group), size = 0.5) %>% 
  collapse(node=2306) %<+%
  sample_dat + #theme_tree2() +
  scale_color_manual(name="Clades",
                     labels=c("Lineage A", "Lineage B"),
                     values = values,
                     guide = 'none') +
  geom_polygon(data = data.frame(x=c(0.028, 0.2, 0.2), y= c(611, 591, 631)),
               aes(x, y), fill = values2[1], color = 'black',
               inherit.aes = F) +
  new_scale_color() + 
  new_scale_fill() +
  geom_tippoint(mapping = aes(fill = lineage), 
                size = 5, shape=21) +
  scale_fill_manual(name="",
                    values=alpha(values2, alpha = 0.9)) + 
  # guides(shape = guide_legend(order = 2), fill = guide_legend(order = 1)) +
  theme(axis.title=element_text(size=14),
        axis.text=element_text(size=12),
        legend.key.size = unit(0.2, 'cm'),
        # legend.text = element_text(size = 10),
        legend.spacing = unit(0.0001, 'cm'),
        legend.position = c(0.7,0.55),
        plot.background = element_blank(),
        strip.background = element_blank(),
        panel.background = element_blank(),
        # legend.title = element_text(size=14),
        # legend.text = element_text(size=12),
        legend.background = element_blank(),
        legend.box.background = element_blank(),
        legend.key = element_blank()
  ) +
  guides(fill = guide_legend(override.aes = list(size=2),
                             ncol = 4)) +
  coord_cartesian(ylim = c(-5, 650)) +
  geom_text(data = data.frame(x = c(0.16,
                                    0.015,0.07,0.14,
                                    0.145,0.168,0.3,
                                    0.45,0.63,0.75,
                                    0.74,0.8,0.84,
                                    0.91,1.1,1.02), 
                              y = c(611,
                                    90,300,29,
                                    441,255,80,
                                    561,568,110,
                                    517,460,494,
                                    177,169,544), 
                              text = c('B',
                                       'A', 'A.1', 'A.2',
                                       'A.5','A.3','A.2.4',
                                       'A.23','A.23.1','A.2.5',
                                       'A.28','A.21','A.27',
                                       'A.2.5.2','A.2.5.1','A.29')), 
            aes(x = x, y = y, label = text), 
            size = 3.8, alpha=.8, fontface='bold',
            inherit.aes = FALSE) +
  scale_x_continuous() 

p


pdf(file = 'Output/F1E_tree_A.pdf', width = 4.4, height = 3)
print(p)
dev.off()

