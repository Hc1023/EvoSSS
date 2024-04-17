rm(list = ls())
library(treeio)
library(ggtree)
library(ggplot2)
library(ape)
library(scales)
library(ggnewscale)


tr <- rtree(10)
ggtree(tr) + geom_tiplab()
ggtree(tr) + geom_tiplab(align = TRUE)

Ttree <- read.tree('tree_early.nwk')
# clock rate = 0.0008
Ttree$edge.length <- Ttree$edge.length/0.0008

dat = read.csv("clinical_info.csv")
for (i in 1:nrow(dat)) {
  a <- dat[i,]
  b <- strsplit(a$Date, '/')[[1]]
  if(b[2] == 1) b[2] = 'Jan'
  if(b[2] == 2) b[2] = 'Feb'
  if(nchar(b[3]) == 1) b[3] = paste0('0',b[3])
  
  dat$tip[i] <- paste0(dat$Names_ID[i], ' ', b[2], ' ', b[3])
  dat$tip.lable[i] <- strsplit(a$sample_id, '-')[[1]][8]
}

sample_dat = dat[,c(7,6,2)]
tmp <- table(sample_dat$Names_ID)
sample_dat$Names_ID[!sample_dat$Names_ID %in% names(tmp)[tmp > 1]] <- 'Others'
sample_dat[nrow(sample_dat)+1,] <- c("Wuhan-Hu-1/2019", "Hu-1", 'Others')

Ttree.drop <- drop.tip(phy = Ttree,
                       tip = Ttree$tip[!Ttree$tip %in% sample_dat$tip.lable])

nextclade <- read.table(file = 'nextclade_qc_early.tsv', 
                        sep = '\t', header = TRUE)
nextclade <- nextclade[nextclade$seqName %in% Ttree.drop$tip.label, c(1,2,3,16)]

for (i in 1:nrow(nextclade)) {
  if(any(grep('G26144T', nextclade[i,4])) & 
     any(grep('G11083T', nextclade[i,4]))){
    nextclade[i,'cls'] <- 'V'
  }else{
    nextclade[i,'cls'] <- nextclade[i,3]
  }
}

cls <- list(Clade_1=nextclade$seqName[nextclade$cls == 'A'],
            Clade_2=nextclade$seqName[nextclade$cls == 'B'],
            Clade_3=nextclade$seqName[nextclade$cls == 'V'])
Ttree.drop <- groupOTU(Ttree.drop, cls)

for (i in 1:length(Ttree.drop$tip.label)) {
  Ttree.drop$tip.label[i] <- sample_dat$tip[sample_dat$tip.lable == Ttree.drop$tip.label[i]]
}

n = table(sample_dat$Names_ID)
# Others    P01    P02    P03    P04    P06    P08    P09    P10    P11    P12 
# 29      3      2      2      6      2      5      3      3      4      2 
# P15    P16    P17    P19    P20    P21    P22    P24    P25    P27    P30 
# 4      3      4      2      4      4      2      2      3      3      3 
# P31    P32    P33    P34    P35    P39    P40    P46    P47    P48    P51 
# 4      3      4      5      3      2      3      3      3      2      3 
# P52    P53    P54    P55    P56    P58 
# 7      2      4      3      2      2 

# > length(n) - 1
# [1] 38

sample_dat1 <- sample_dat
sample_dat <- sample_dat[,-1]
namesid <- table(sample_dat$Names_ID)
sample_dat$Names_ID = factor(sample_dat$Names_ID)

values = c(hue_pal()(3)[1], hue_pal()(3)[3], '#5fbfb1')
options(ignore.negative.edge=TRUE)

sizetext = 4

c39 = alpha(c('grey', rainbow(38)), 0.6)
p = ggtree(Ttree.drop, mrsd = '2020-03-01', as.Date = T, show.legend = F, 
       layout = "circular", aes(color = group), size = 0.3) %<+% 
  sample_dat +
  scale_color_manual(name="Clades",
                     labels=c("Lineage A", "Lineage B",
                              "Clade V"),
                     values = values) +
  new_scale_color() + 
  new_scale_fill() +
  geom_tiplab(aes(color = Names_ID), offset = 10, 
              size = 1.5, geom = "text", align = T, 
              show.legend = FALSE) +
  scale_color_manual(name='Patients', values = c39) +
  geom_tippoint(
    mapping = aes(color = Names_ID),         
    size = 2, shape = 19, show.legend = F) +
  scale_fill_manual(name='Patients', 
                    values=c39)

pdf(file = 'Output/tree_clinical_only.pdf', width = 4, height = 4)
print(p)
dev.off()
