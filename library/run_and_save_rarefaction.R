#library(vegan)
library(dplyr)
source("R_scripts/library.R")

#read the data
asv_tab = read.csv("data/asv_tab.csv", row.names = 1)
ncol(asv_tab)
nrow(asv_tab)
colnames(asv_tab)
rownames(asv_tab)

rowSums(asv_tab) %>% sort
#how much greater is max sample size than min
max(rowSums(asv_tab))/min(rowSums(asv_tab)) 
#297
#1000 rarefactions should be fine to get a represenative sample
# note, we ran 10K and the resulting list was 3.8 Gb which would be a little annoying to deal with
# plus bitbucket recommended max repo size is 2 G

min_seqs = min(rowSums(asv_tab)) 
#perform rarefaction on vol and wash tables
rarefactions_list = multiple_subsamples(x = asv_tab, depth = min_seqs, iterations = 1000) 

saveRDS(rarefactions_list, "data/rarefactions.rds")


