library(dplyr)

lib_read_counts = read.csv("data/sample_metadata/all_libs.reads_per_sample.csv", header = F)
head(lib_read_counts)
colnames(lib_read_counts) = c("seqLabel", "seq_count", "lib_type")

plot2gel = read.csv("data/sample_metadata/plot-gel.csv")
head(plot2gel)
sample2gel.ITS = read.csv("data/sample_metadata/sample-gel.ITS.csv")
head(sample2gel.ITS)
sample2gel.16S = read.csv("data/sample_metadata/sample-gel.16S.csv")
sample2gel.18S = read.csv("data/sample_metadata/sample-gel.18S.csv")

lib_read_counts$SampleID = sub("_S\\d+", "", lib_read_counts$seqLabel, perl = T)
head(lib_read_counts)

metadata.ITS = left_join(lib_read_counts %>% 
    filter(lib_type == "FUN"), sample2gel.ITS, by = "SampleID") %>%
    left_join(., plot2gel, by = "gelKey")

metadata.16S = left_join(lib_read_counts %>% 
    filter(lib_type == "BAC"), sample2gel.16S, by = "SampleID") %>%
    left_join(., plot2gel, by = "gelKey")

metadata.18S = left_join(lib_read_counts %>% 
    filter(lib_type == "EUK"), sample2gel.18S, by = "SampleID") %>%
    left_join(., plot2gel, by = "gelKey")

write.csv(metadata.ITS, "data/sample_metadata/metadata.ITS.csv", row.names = F, quote = F)
write.csv(metadata.ITS, "data/sample_metadata/metadata.16S.csv", row.names = F, quote = F)
write.csv(metadata.ITS, "data/sample_metadata/metadata.18S.csv", row.names = F, quote = F)
