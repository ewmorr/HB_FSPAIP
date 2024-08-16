library(vegan)
library(dplyr)
source("library/library.R")

######################
# senescent litter
rarefactions_list = readRDS("data/ITS/rarefactions.sen.rds")

#calc distance and alpha-div
bray_binary_list = lapply(rarefactions_list, vegdist, method = "bray", binary = T, diag = T, upper = T)
bray_logCts_list = lapply(rarefactions_list, log_dist, method = "bray")
shannon_list = lapply(rarefactions_list, diversity, index = "shannon")
simpson_list = lapply(rarefactions_list, diversity, index = "simpson")
richness_list = lapply(rarefactions_list, richness_calc)

#############################################
#Take avgs (for dists convert to dist object)
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
saveRDS(bray_binary_avg, "data/ITS/avg_dist/bray-binary.sen.rds")
saveRDS(bray_logCts_avg, "data/ITS/avg_dist/bray-logCounts.sen.rds")
write.csv(div_avg, "data/ITS/avg_div/diversity.sen.csv", row.names = F)

write.csv(avg_counts, "data/ITS/asv_tab.rarefaction_avg.sen.csv")


######################
# Oi litter
rarefactions_list = readRDS("data/ITS/rarefactions.oi.rds")

#calc distance and alpha-div
bray_binary_list = lapply(rarefactions_list, vegdist, method = "bray", binary = T, diag = T, upper = T)
bray_logCts_list = lapply(rarefactions_list, log_dist, method = "bray")
shannon_list = lapply(rarefactions_list, diversity, index = "shannon")
simpson_list = lapply(rarefactions_list, diversity, index = "simpson")
richness_list = lapply(rarefactions_list, richness_calc)

#############################################
#Take avgs (for dists convert to dist object)
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
saveRDS(bray_binary_avg, "data/ITS/avg_dist/bray-binary.oi.rds")
saveRDS(bray_logCts_avg, "data/ITS/avg_dist/bray-logCounts.oi.rds")
write.csv(div_avg, "data/ITS/avg_div/diversity.oi.csv", row.names = F)

write.csv(avg_counts, "data/ITS/asv_tab.rarefaction_avg.oi.csv")

######################
# ohor
rarefactions_list = readRDS("data/ITS/rarefactions.ohor.rds")

#calc distance and alpha-div
bray_binary_list = lapply(rarefactions_list, vegdist, method = "bray", binary = T, diag = T, upper = T)
bray_logCts_list = lapply(rarefactions_list, log_dist, method = "bray")
shannon_list = lapply(rarefactions_list, diversity, index = "shannon")
simpson_list = lapply(rarefactions_list, diversity, index = "simpson")
richness_list = lapply(rarefactions_list, richness_calc)

#############################################
#Take avgs (for dists convert to dist object)
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
saveRDS(bray_binary_avg, "data/ITS/avg_dist/bray-binary.ohor.rds")
saveRDS(bray_logCts_avg, "data/ITS/avg_dist/bray-logCounts.ohor.rds")
write.csv(div_avg, "data/ITS/avg_div/diversity.ohor.csv", row.names = F)

write.csv(avg_counts, "data/ITS/asv_tab.rarefaction_avg.ohor.csv")

######################
# mineral litter
rarefactions_list = readRDS("data/ITS/rarefactions.min.rds")

#calc distance and alpha-div
bray_binary_list = lapply(rarefactions_list, vegdist, method = "bray", binary = T, diag = T, upper = T)
bray_logCts_list = lapply(rarefactions_list, log_dist, method = "bray")
shannon_list = lapply(rarefactions_list, diversity, index = "shannon")
simpson_list = lapply(rarefactions_list, diversity, index = "simpson")
richness_list = lapply(rarefactions_list, richness_calc)

#############################################
#Take avgs (for dists convert to dist object)
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
saveRDS(bray_binary_avg, "data/ITS/avg_dist/bray-binary.min.rds")
saveRDS(bray_logCts_avg, "data/ITS/avg_dist/bray-logCounts.min.rds")
write.csv(div_avg, "data/ITS/avg_div/diversity.min.csv", row.names = F)

write.csv(avg_counts, "data/ITS/asv_tab.rarefaction_avg.min.csv")
