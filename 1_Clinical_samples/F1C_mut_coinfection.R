rm(list= ls())
library(dplyr)
library(tidyverse)
library(ggplot2)
library(scales)
library(knitr)
library(RColorBrewer)
library(ggsci)

df = read.csv('F1C_clinical_coverage.csv')

df1 = df[df$pos == 8782,
         c('Sample', 'Names_ID','total_reads',
           'b1','b2','n')]
df2 = df[df$pos == 28144,
         c('Sample', 'Names_ID','total_reads',
           'b1','b2','n')]

mygather = function(df1, y=100){
  a = gather(df1, key = "Mutations", value = "Counts", -colnames(df1)[1:3])
  a$Sample = factor(a$Sample, levels = df1$Sample)
  tmp1 = df1[df1$b1 > y & df1$b2 > y,]
  tmp2 = df1[(df1$b1 > 20 & df1$b2 > 20)&(!(df1$b1 > 100 & df1$b2 > 100)),]
  return(list(a, tmp1, tmp2))
}

mylist = mygather(df1)
dfa = mylist[[1]]
dftmp1 = mylist[[2]]
dftmp2 = mylist[[3]]

p = ggplot(data = dfa, 
           mapping = aes(x = Sample, y = Counts, 
                         fill = Mutations)) +
  geom_col() + theme_classic() +
  scale_fill_manual(values = c(alpha(c('#45cbba', 'firebrick3', 'black'), 0.8)),
                    labels = c('C8782','C8782T','Others'),
                    name = '') +
  labs(x = 'Samples', y = 'Number of reads') +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        legend.position = 'top') +
  scale_y_continuous(expand = c(0.01, 0), 
                     limits = c(0, max(dfa$total_reads)+300)) +
  geom_text(data = dftmp1, aes(x = Sample, y = total_reads),
            label = '*', size = 2.5, color = 'red', inherit.aes = F) +
  geom_text(data = dftmp2, aes(x = Sample, y = total_reads),
            label = '.', inherit.aes = F,
            size = 3, vjust = -0.2)
p
pdf(file = 'Output/F1C_num_reads_8782.pdf', width = 3.6, height = 1.8)
print(p)
dev.off()

mylist = mygather(df1[df1$total_reads < 200,], 20)

dfa_sub = mylist[[1]]
dftmp1_sub = mylist[[2]]
p = ggplot(data = dfa_sub, 
           mapping = aes(x = Sample, y = Counts, 
                         fill = Mutations)) +
  geom_col() +
  theme_classic()+
  scale_fill_manual(values = c(alpha(c('#45cbba', 'firebrick3', 'black'), 0.8)),
                    labels = c('C8782','C8782T','Others'),
                    guide  = 'none') +
  labs(x = '', y = '') +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) +
  scale_y_continuous(expand = c(0.01, 0), 
                     limits = c(0, max(dfa_sub$total_reads)+5)) +
  geom_text(data = dftmp1_sub, aes(x = Sample, y = total_reads),
            label = '.', inherit.aes = F, size = 3, vjust = -0.2)
p
pdf(file = 'Output/F1C_num_reads_8782_sub.pdf', width = 1.7, height = 1)
print(p)
dev.off()


mylist = mygather(df2)
dfb = mylist[[1]]
dfbtmp1 = mylist[[2]]
dfbtmp2 = mylist[[3]]
p = ggplot(data = dfb, mapping = aes(x = Sample, y = Counts, fill = Mutations)) +
  geom_col() +
  theme_classic()+
  scale_fill_manual(values = alpha(c('#206fc8', 'orange', 'black'), 0.8),
                    labels = c('T28144','T28144C','Others'),
                    name = '') +
  labs(x = 'Samples', y = 'Number of reads') +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        legend.position = 'top') +
  scale_y_continuous(expand = c(0.01, 0), 
                     limits = c(0, max(dfb$total_reads)+300)) +
  geom_text(data = dfbtmp1, aes(x = Sample, y = total_reads),
            label = '*', size = 2.5, color = 'red', inherit.aes = F) +
  geom_text(data = dfbtmp2, aes(x = Sample, y = total_reads),
            label = '.', inherit.aes = F,size = 3, vjust = -0.2)
p
pdf(file = 'Output/F1C_num_reads_28144.pdf', width = 3.6, height = 1.8)
print(p)
dev.off()

mylist = mygather(df2[df2$total_reads < 200,], 20)

dfb_sub = mylist[[1]]
dfbtmp1_sub = mylist[[2]]
p = ggplot(data = dfb_sub, mapping = aes(x = Sample, y = Counts, fill = Mutations)) +
  geom_col() +
  theme_classic()+
  scale_fill_manual(values = alpha(c('#206fc8', 'orange', 'black'), 0.8),
                    labels = c('T28144','T28144C','Others'),
                    guide = 'none') +
  labs(x = '', y = '') +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) +
  scale_y_continuous(expand = c(0.01, 0), 
                     limits = c(0, max(dfb_sub$total_reads)+8)) +
  geom_text(data = dfbtmp1_sub, aes(x = Sample, y = total_reads),
            label = '.', inherit.aes = F, size = 3, vjust = -0.2)
p
pdf(file = 'Output/F1C_num_reads_28144_sub.pdf', width = 1.2, height = 1.2)
print(p)
dev.off()
