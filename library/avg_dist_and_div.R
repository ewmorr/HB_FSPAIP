library(vegan)
library(dplyr)
source("R_scripts/library.R")

rarefactions_list = readRDS("data/rarefactions.rds")

##############################################
#Run distance calcs
#saving square matrix for downstream averaging
bray_list = lapply(rarefactions_list, vegdist, method = "bray", binary = F, diag = T, upper = T)
bray_binary_list = lapply(rarefactions_list, vegdist, method = "bray", binary = T, diag = T, upper = T)
bray_logCts_list = lapply(rarefactions_list, log_dist, method = "bray")

#############################################
#Run diversity calcs
shannon_list = lapply(rarefactions_list, diversity, index = "shannon")
simpson_list = lapply(rarefactions_list, diversity, index = "simpson")
richness_list = lapply(rarefactions_list, richness_calc)

#############################################
#Take avgs (for dists convert to dist object)
bray_avg = avg_matrix_list(bray_list) %>% as.dist()
bray_binary_avg = avg_matrix_list(bray_binary_list) %>% as.dist()
bray_logCts_avg = avg_matrix_list(bray_logCts_list) %>% as.dist()

div_avg = data.frame(
    sample = rownames(avg_matrix_list(shannon_list)),
    shannon = avg_matrix_list(shannon_list),
    simpson = avg_matrix_list(simpson_list),
    richness = avg_matrix_list(richness_list),
    stringsAsFactors = F
)

###########################################
#Also average the counts table
avg_counts = avg_matrix_list(rarefactions_list)

##########################################
# write files
saveRDS(bray_avg, "data/avg_dist/bray.rds")
saveRDS(bray_binary_avg, "data/avg_dist/bray-binary.rds")
saveRDS(bray_logCts_avg, "data/avg_dist/bray-logCounts.rds")

write.csv(div_avg, "data/avg_div/diversity.csv", row.names = F)

write.csv(avg_counts, "data/asv_tab.rarefaction_avg.csv")
