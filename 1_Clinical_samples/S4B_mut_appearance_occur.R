rm(list = ls())
library(data.table)
library(tidyverse)
library(ggrepel)

dattmp = read.csv("F1A_clinical_info.csv")
snpdate_all = read.csv('S4B_snpdate_seq3gr10.tsv', sep = '\t')

load(file = 'S4A_clinical_mutation.Rdata')
## count in patients
for (i in 1:nrow(dat_mut)) {
  tmp <-  dat_mut[i,]
  dat_mut[i,'mut'] <- paste0(tmp[3],tmp[2],tmp[4])
}
for (i in 1:nrow(dat_mut)) {
  dat_mut$patient[i] <- dattmp$Names_ID[dattmp$Sample == dat_mut$sample[i]]
}

### 0.05
dat_mut_table <- unique(dat_mut[c('mut', "patient")])
dat_mut_table <- sort(table(dat_mut_table$mut), decreasing = T) %>% 
  as.data.frame()
dat_mut_table <- dat_mut_table[dat_mut_table$Freq>5,]

colnames(dat_mut_table)[1] = 'Mut'
snpdate = snpdate_all[snpdate_all$Mut %in% dat_mut_table$Mut,]

df = merge(dat_mut_table, snpdate, by = 'Mut', all = T)

df = df[df$Date != '',]
df$Date = as.Date(df$Date)

p = ggplot(df, aes(Date, Freq)) + geom_point(alpha = 0.6) + 
  theme_classic() +
  labs(y = 'Number of patients (>5)\n (VAF>0.05)', 
       x='Emergence date in GISAID (Accessed by October 8, 2022)') +
  scale_x_date(date_labels = "%y %b", breaks = "3 month") +
  geom_text_repel(aes(label=Mut)) 
p
pdf(file = 'Output/S4B_mut_appearance_occur.pdf', width = 6, height = 2.7)
print(p)
dev.off()
