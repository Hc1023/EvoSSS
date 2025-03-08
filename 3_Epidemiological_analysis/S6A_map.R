rm(list = ls())
library(treeio)
library(ggtree)
library(ggplot2)
library(ape)
library(scales)
library(scatterpie)

if(F){
  Ttree <- read.tree('../data/ZJU_sample/tree_early.nwk')
  scale_length <- Ttree$edge.length/0.0008
  Ttree$edge.length <- scale_length
  Ttree.drop <- drop.tip(phy = Ttree, tip = 'Wuhan/WH01/2019')
  Ttree.tip.label <- data.frame(Ttree.drop$tip.label)
  nextclade <- read.table(file = '../data/ZJU_sample/nextclade_qc_early.tsv', 
                          sep = '\t', header = TRUE)
  nextclade <- nextclade[nextclade$seqName %in% Ttree.drop$tip.label, c(1,2,3,16)]
  nextclade <- nextclade[-c(1:149),]
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
  dfregion = data.frame()
  
  for (l in unique(nextclade$V5)) {
    df = nextclade[nextclade$V5 == l,]
    clregion = c()
    for (i in 1:nrow(df)) {
      if(any(grep('/', df[i,1]))){
        clregion <-c(clregion, strsplit(df[i,1], '/')[[1]][1])
      }
    }
    df0 <- data.frame(table(clregion))
    df0$group <- l
    colnames(df0)[1] <- 'region'
    dfregion = rbind(dfregion, df0)
    
  }
  df = dfregion
  rm(dfregion)
  
  ddf <-data.frame(region = unique(df$region))
  for (l in unique(df$group)) {
    ddf[,l] = 0
  }
  
  for (i in 1:nrow(ddf)) {
    x <- df[df$region %in% ddf[i,1],]
    for (l in unique(df$group)) {
      if(any(x[x$group == l,2])){
        ddf[i,l] = x[x$group == l,2]
      }else{
        ddf[i,l] = 0
      }
      
    }
  }
  
  library(RJSONIO)
  library(ggmap) 
  ddf$lat <- ''
  ddf$lon <- ''
  for (i in 1:nrow(ddf)) {
    url <- paste0('http://nominatim.openstreetmap.org/search?q=', ddf[i,1], '&limit=9&format=json')
    loc <- fromJSON(url)
    if(length(loc) > 0){
      loc <- loc[[1]]
      ddf$lat[i] <- loc$lat
      ddf$lon[i] <- loc$lon
    }
  }
  
  ddf[ddf$region=='CanaryIslands',c('lat','lon')] = c(28.291565,-16.629129)
  ddf[ddf$region=='CzechRepublic',c('lat','lon')] = c(50.073658,14.418540)
  ddf[ddf$region=='UnitedArabEmirates',c('lat','lon')] = c(37.1532,46.4964)
  ddf[ddf$region=='NewZealand',c('lat','lon')] = c(-40.9006,174.8860)
  ddf[ddf$region=='DominicanRepublic',c('lat','lon')] = c(18.7357,-70.1627)
  ddf[ddf$region=='SaintBarthelemy',c('lat','lon')] = c(17.9139,-62.8339)
  ddf <- ddf[-which(ddf$region == 'env'),]
  # ddf[ddf$region == 'Wuhan',2] <- ddf[ddf$region == 'Wuhan',2]+3
  # ddf[ddf$region == 'USA',3] <- ddf[ddf$region == 'USA',3]+1
  ddf <- ddf[-which(ddf$region == 'Wuhan-Hu-1'),]
  ddf[ddf$region == 'Wuhan',3] <- ddf[ddf$region == 'Wuhan',3] + 1
  
  ddf$lat <- as.numeric(ddf$lat)
  ddf$lon <- as.numeric(ddf$lon)
  
  write.csv(ddf,'tree_map_ori.csv')
  
}

ddf = read.csv('S6A_map.csv')
colnames(ddf)[2:5] = c('Lineage A', 'Lineage B', '8782T','28144C')
ddf$radius <- log(apply(ddf[,2:5], 1, sum))/log(5) + 1

# Pacific, Atlantic, Indian, ArcticSouthern Ocean which is off the coast of Antarctica.

anns = data.frame(CONTINENT = c("Africa", "Asia", "Europe", "North America",
                                "Oceania", "South America"),
                  continent_long = c(15, 80, 20, -100, 140, -60),
                  continent_lat  = c(15, 35, 50, 40, -22, -15),
                  stringsAsFactors = F)
anns2 = data.frame(OCEAN = c("North\nPacific\nOcean", 
                             "South\nPacific\nOcean",
                             "North\nAtlantic\nOcean", 
                             "South\nAtlantic\nOcean",
                             "Indian\nOcean", 
                             "Arctic Ocean"),
                   ocean_long = c(-149, -139, -40, -19, 82, 90),
                   ocean_lat  = c(40, -38, 32, -34, -33, 85),
                   stringsAsFactors = F)

values = c(hue_pal()(3)[1], hue_pal()(3)[3], hue_pal()(3)[2], hue_pal()(4)[4])
ddf$radius = ddf$radius*3

world <- map_data('world')
p <- ggplot(world, aes(long, lat)) +
  geom_map(map=world, aes(map_id=region), 
           fill="#F7F6F4", color="#CCD2D3", size = 0.3) +
  geom_text(data = anns, color = "#777777", alpha=.8,
            aes(x = continent_long, y = continent_lat, label = CONTINENT)) +
  geom_text(data = anns2, color = "#777777", fontface = 'italic', alpha=.5,
            aes(x = ocean_long, y = ocean_lat, label = OCEAN)) +
  geom_scatterpie(data=ddf, 
                  aes(x=lon, y=lat, group=region, r=radius), 
                  cols=c('Lineage A', 'Lineage B'), 
                  color = alpha('black', 0.9), size = 0.1) + 
  scale_fill_manual(name="Variant",
                    values = alpha(values,.8),
                    labels = c('A', 'B')) +
  coord_equal() + 
  geom_scatterpie_legend(ddf$radius, x=-160, y=-55, n=2, breaks = c(6,12),
                         labeller=function(x) 5^(x/3-1)) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = alpha("#CCD2D3", 0.5)),
        panel.ontop = FALSE,
        axis.title=element_text(size=16),
        axis.text=element_text(size=14),
        legend.title = element_text(size=16),
        legend.text = element_text(size=14),
        legend.position = 'top') +
  xlab('Longitude') + ylab('Lattitude') +
  geom_text(data = data.frame(x = 155,
                              y = -85, 
                              text = 'By March 1, 2020'), 
            aes(x = x, y = y, label = text), 
            size = 4, alpha=.8,
            inherit.aes = FALSE)

# ocean #777777
pdf(file = 'Output/S6A_map.pdf', height = 3.5)
print(p)
dev.off()






