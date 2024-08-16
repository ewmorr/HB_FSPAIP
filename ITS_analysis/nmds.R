library(dplyr)
library(vegan)
library(ggplot2)
source("library/ggplot_theme.txt")

bray.sen = readRDS("data/ITS/avg_dist/bray-binary.sen.rds")
bray.oi = readRDS("data/ITS/avg_dist/bray-binary.oi.rds")
bray.ohor = readRDS("data/ITS/avg_dist/bray-binary.ohor.rds")
bray.min = readRDS("data/ITS/avg_dist/bray-binary.min.rds")

bray.log.sen = readRDS("data/ITS/avg_dist/bray-logCounts.sen.rds")
bray.log.oi = readRDS("data/ITS/avg_dist/bray-logCounts.oi.rds")
bray.log.ohor = readRDS("data/ITS/avg_dist/bray-logCounts.ohor.rds")
bray.log.min = readRDS("data/ITS/avg_dist/bray-logCounts.min.rds")

nmds.sen = metaMDS(bray.log.sen)
nmds.oi = metaMDS(bray.log.oi) 
nmds.ohor = metaMDS(bray.log.ohor)
nmds.min = metaMDS(bray.log.min)

saveRDS(nmds.sen, file = "data/ITS/nmds.sen.rds")
saveRDS(nmds.oi, file = "data/ITS/nmds.oi.rds")
saveRDS(nmds.ohor, file = "data/ITS/nmds.ohor.rds")
saveRDS(nmds.min, file = "data/ITS/nmds.min.rds")

plot(nmds.sen)
plot(nmds.oi)
plot(nmds.ohor)
plot(nmds.min)

pcoa.sen = capscale(bray.log.sen ~ 1)
pcoa.oi = capscale(bray.log.oi ~ 1)
pcoa.ohor = capscale(bray.log.ohor ~ 1)
pcoa.min = capscale(bray.log.min ~ 1)

pcoa.sen
pcoa.oi
pcoa.ohor
pcoa.min

pcoa.sen.eig = pcoa.sen$CA$eig/sum(pcoa.sen$CA$eig)
pcoa.oi.eig = pcoa.oi$CA$eig/sum(pcoa.oi$CA$eig)
pcoa.ohor.eig = pcoa.ohor$CA$eig/sum(pcoa.ohor$CA$eig)
pcoa.min.eig = pcoa.min$CA$eig/sum(pcoa.min$CA$eig)

pcoa.sen.eig
pcoa.oi.eig
pcoa.ohor.eig
pcoa.min.eig

plot(pcoa.sen)
plot(pcoa.min)

################
#construct dfs
metadata = read.csv("data/sample_metadata/metadata.ITS.csv")


all.sites_scores = left_join(
        data.frame(
            seqLabel = c(
            rownames(scores(pcoa.sen)$sites),
            rownames(scores(pcoa.oi)$sites),
            rownames(scores(pcoa.ohor)$sites),
            rownames(scores(pcoa.min)$sites)
        ),
        rbind(
            scores(pcoa.sen)$sites,
            scores(pcoa.oi)$sites,
            scores(pcoa.ohor)$sites,
            scores(pcoa.min)$sites
        )
    ),
    metadata,
    by = "seqLabel"
)

layer_labels = c("Ash leaves", "Litter (Oi)", "O horizon (Oe/Oa)", "Mineral (0-10 cm)")
names(layer_labels) = c("senescent", "Oi", "Oe-Oa", "mineral")

p1 = ggplot(all.sites_scores,
       aes(x = MDS1, y = MDS2, color = treatment)
) +
    geom_point(size = 3) +
    facet_wrap(
        ~factor(
            soil_layer, 
            levels = c("senescent", "Oi", "Oe-Oa", "mineral"),
            labels = c("Senescent ash leaves", "Litter (Oi)", "O horizon (Oe/Oa)", "Mineral (0-10 cm)")
        ), 
        labeller = labeller(soil_layer = layer_labels)
    ) +
    scale_color_manual(values = c("#D95F02", "#1B9E77"), labels = c("Untreated", "Treated")) +
    labs(color = "Emamectin benzoate", x = "PCo1 (12-14% variance)", y = "PCo2 (10-13% variance)") +
    my_gg_theme +
    theme(
        legend.position = "bottom",
        legend.title = element_text(size = 18),
        legend.text = element_text(size = 18)
    )
p1
pdf("figures/ITS/HB_fungi_pcoa_leaves_litter_Ohor_min.pdf", width = 7, height = 5)
p1
dev.off()

#################
#adonis
#

adonis2(bray.log.sen ~ treatment, data = all.sites_scores %>% filter(., soil_layer == "senescent"))
# P = 0.797, R2 = 0.05
adonis2(bray.log.oi ~ treatment, data = all.sites_scores %>% filter(., soil_layer == "Oi"))
# P = 0.974, R2 = 0.04
adonis2(bray.log.ohor ~ treatment, data = all.sites_scores %>% filter(., soil_layer == "Oe-Oa"))
# P = 0.92, R2 = 0.07



