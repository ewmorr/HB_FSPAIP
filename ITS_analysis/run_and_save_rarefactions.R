library(dplyr)
source("library/library.R")

asv_tab = read.table("data/ITS/ASVs_counts.tsv", header = T)
ncol(asv_tab)
nrow(asv_tab)
colnames(asv_tab)
rownames(asv_tab)

metadata = read.csv("data/sample_metadata/metadata.ITS.csv")
head(metadata)

senescent_ids = metadata %>% filter(., soil_layer == "senescent") %>% pull(., seqLabel)
oi_ids = metadata %>% filter(., soil_layer == "Oi") %>% pull(., seqLabel)
oHor_ids = metadata %>% filter(., soil_layer == "Oe-Oa") %>% pull(., seqLabel)
mineral_ids = metadata %>% filter(., soil_layer == "mineral") %>% pull(., seqLabel)


asv_tab.sen = asv_tab[,colnames(asv_tab) %in% senescent_ids]
asv_tab.sen = asv_tab.sen[rowSums(asv_tab.sen) > 1, colSums(asv_tab.sen) > 999] %>% t() #remove 0 count asvs, filter to min 1k seqs per sample, and transpose
ncol(asv_tab.sen)
nrow(asv_tab.sen)
asv_tab.oi = asv_tab[,colnames(asv_tab) %in% oi_ids]
asv_tab.oi = asv_tab.oi[rowSums(asv_tab.oi) > 1, colSums(asv_tab.oi) > 999] %>% t()
ncol(asv_tab.oi)
nrow(asv_tab.oi)
asv_tab.ohor = asv_tab[,colnames(asv_tab) %in% oHor_ids]
asv_tab.ohor = asv_tab.ohor[rowSums(asv_tab.ohor) > 1, colSums(asv_tab.ohor) > 400] %>% t() #ohor is low sow go with 500 min
ncol(asv_tab.ohor)
nrow(asv_tab.ohor)
asv_tab.min = asv_tab[,colnames(asv_tab) %in% mineral_ids]
asv_tab.min = asv_tab.min[rowSums(asv_tab.min) > 1, colSums(asv_tab.min) > 999] %>% t()
ncol(asv_tab.min)
nrow(asv_tab.min)

rowSums(asv_tab.sen) %>% sort
rowSums(asv_tab.oi) %>% sort
rowSums(asv_tab.ohor) %>% sort
rowSums(asv_tab.min) %>% sort

max(rowSums(asv_tab.sen))/min(rowSums(asv_tab.sen)) 
max(rowSums(asv_tab.oi))/min(rowSums(asv_tab.oi)) 
max(rowSums(asv_tab.ohor))/min(rowSums(asv_tab.ohor)) 
max(rowSums(asv_tab.min))/min(rowSums(asv_tab.min)) 
#1000 rarefactions should be fine to get a represenative sample all ratios <255


###############
#perform rarefactions
#

min_seqs.sen = min(rowSums(asv_tab.sen)) 
rarefactions_list.sen = multiple_subsamples(x = asv_tab.sen, depth = min_seqs.sen, iterations = 1000) 
saveRDS(rarefactions_list.sen, "data/ITS/rarefactions.sen.rds")

min_seqs.oi = min(rowSums(asv_tab.oi)) 
rarefactions_list.oi = multiple_subsamples(x = asv_tab.oi, depth = min_seqs.oi, iterations = 1000) 
saveRDS(rarefactions_list.oi, "data/ITS/rarefactions.oi.rds")

min_seqs.ohor = min(rowSums(asv_tab.ohor)) 
rarefactions_list.ohor = multiple_subsamples(x = asv_tab.ohor, depth = min_seqs.ohor, iterations = 1000) 
saveRDS(rarefactions_list.ohor, "data/ITS/rarefactions.ohor.rds")

min_seqs.min = min(rowSums(asv_tab.min)) 
rarefactions_list.min = multiple_subsamples(x = asv_tab.min, depth = min_seqs.min, iterations = 1000) 
saveRDS(rarefactions_list.min, "data/ITS/rarefactions.min.rds")

