rm(list = ls())
library(tidyverse)
library(ggplot2)
library(ggpubr)

dattmp = read.csv("F1A_clinical_info.csv")
load(file = 'S2A_clinical_mutation.Rdata')

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

p1 = dat_mut_table  %>%
  mutate(name = fct_reorder(Var1, Freq)) %>%
  ggplot(aes(x=name, y=Freq)) +
  geom_bar(stat="identity", fill="#f68060", alpha=.6, width=.6) +
  scale_y_continuous(limits = c(0, 46)) +
  geom_text(aes(label = Freq),
            position = position_dodge(width = 0.4), 
            hjust = -0.2) +
  coord_flip() +
  theme_bw() +
  labs(x = '', y = "VAF > 0.05")

# Number of patients (>5)
### 0.5
dat_mut_table <- unique(dat_mut[dat_mut$percent > 0.5, c('mut', "patient")])
dat_mut_table <- sort(table(dat_mut_table$mut), decreasing = T) %>% 
  as.data.frame()
dat_mut_table <- dat_mut_table[dat_mut_table$Freq>5,]

p2 = dat_mut_table  %>%
  mutate(name = fct_reorder(Var1, Freq)) %>%
  ggplot(aes(x=name, y=Freq)) +
  geom_bar(stat="identity", fill="#f68060", alpha=.6, width=.4) +
  scale_y_continuous(limits = c(0, 30)) +
  geom_text(aes(label = Freq),
            position = position_dodge(width = 0.4), 
            hjust = -0.2) +
  coord_flip() +
  theme_bw() +
  labs(x = '', y = "VAF > 0.5")

### 0.8
dat_mut_table <- unique(dat_mut[dat_mut$percent > 0.8, c('mut', "patient")])
dat_mut_table <- sort(table(dat_mut_table$mut), decreasing = T) %>% 
  as.data.frame()
dat_mut_table <- dat_mut_table[dat_mut_table$Freq>5,]

p3 = dat_mut_table  %>%
  mutate(name = fct_reorder(Var1, Freq)) %>%
  ggplot(aes(x=name, y=Freq)) +
  geom_bar(stat="identity", fill="#f68060", alpha=.6, width=.4) +
  scale_y_continuous(limits = c(0, 30)) +
  geom_text(aes(label = Freq),
            position = position_dodge(width = 0.4), 
            hjust = -0.2) +
  coord_flip() +
  theme_bw() +
  labs(x = '', y = "VAF > 0.8")


### Output

pdf(file = paste0('Output/S2A_mutation_appearance.pdf'), 
    width = 3.2, height = 3.4)
print(p1)
dev.off()

pdf(file = paste0('Output/S2A_mutation_appearance_2.pdf'), 
    width = 3.2, height = 1.7)
print(p2)
print(p3)
dev.off()
