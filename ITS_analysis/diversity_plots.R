library(dplyr)
library(ggplot2)
source("library/ggplot_theme.txt")

div.sen = read.csv("data/ITS/avg_div/diversity.sen.csv")
div.oi = read.csv("data/ITS/avg_div/diversity.oi.csv")
div.ohor = read.csv("data/ITS/avg_div/diversity.ohor.csv")
div.min = read.csv("data/ITS/avg_div/diversity.min.csv")

metadata = read.csv("data/sample_metadata/metadata.ITS.csv")

colnames(metadata)[1] = "sample"

div.sen.meta = left_join(div.sen, metadata)
div.oi.meta = left_join(div.oi, metadata)
div.ohor.meta = left_join(div.ohor, metadata)
div.min.meta = left_join(div.min, metadata)

#ANOVA
##richness
aov.sen = aov(richness ~ treatment, data = div.sen.meta)
qqnorm(residuals(aov.sen))
plot(residuals(aov.sen))
summary(aov.sen)

aov.oi = aov(richness ~ treatment, data = div.oi.meta)
qqnorm(residuals(aov.oi))
plot(residuals(aov.oi))
summary(aov.oi)

aov.ohor = aov(richness ~ treatment, data = div.ohor.meta)
qqnorm(residuals(aov.ohor))
plot(residuals(aov.ohor))
summary(aov.ohor)

aov.min = aov(richness ~ treatment, data = div.min.meta)
qqnorm(residuals(aov.min))
plot(residuals(aov.min))
summary(aov.min)

#shannon
aov.sen = aov(shannon ~ treatment, data = div.sen.meta)
qqnorm(residuals(aov.sen))
plot(residuals(aov.sen))
summary(aov.sen)

aov.oi = aov(shannon ~ treatment, data = div.oi.meta)
qqnorm(residuals(aov.oi))
plot(residuals(aov.oi))
summary(aov.oi)

aov.ohor = aov(shannon ~ treatment, data = div.ohor.meta)
qqnorm(residuals(aov.ohor))
plot(residuals(aov.ohor))
summary(aov.ohor)

aov.min = aov(shannon ~ treatment, data = div.min.meta)
qqnorm(residuals(aov.min))
plot(residuals(aov.min))
summary(aov.min)
#marginal p = 0.08

#simpson
aov.sen = aov(simpson ~ treatment, data = div.sen.meta)
qqnorm(residuals(aov.sen))
plot(residuals(aov.sen))
summary(aov.sen)

aov.oi = aov(simpson ~ treatment, data = div.oi.meta)
qqnorm(residuals(aov.oi))
plot(residuals(aov.oi))
summary(aov.oi)

aov.ohor = aov(simpson ~ treatment, data = div.ohor.meta)
qqnorm(residuals(aov.ohor))
plot(residuals(aov.ohor))
summary(aov.ohor)

aov.min = aov(simpson ~ treatment, data = div.min.meta)
qqnorm(residuals(aov.min))
plot(residuals(aov.min))
summary(aov.min)

#plot
all_div = rbind(div.sen.meta, div.oi.meta, div.ohor.meta, div.min.meta)

all_div_summary = all_div %>% 
    group_by(., soil_layer, treatment) %>%
    summarize(., mean_rich = mean(richness), se_rich = sd(richness)/sqrt(n()))

p1 = ggplot(
    all_div_summary,
    aes(
        x = factor(
        soil_layer, 
        levels = c("senescent", "Oi", "Oe-Oa", "mineral"),
        labels = c("Senescent ash leaves", "Litter (Oi)", "O horizon (Oe/Oa)", "Mineral (0-10 cm)")
        ),
    y = mean_rich, 
    fill = treatment
    )
) +
    geom_bar(stat = "identity", position = position_dodge(width = 0.9), width = 0.8) +
    geom_errorbar(
        aes(ymin = mean_rich-se_rich, ymax = mean_rich+se_rich), 
        position = position_dodge(width = 0.9),
        width = 0.1,
        color = "black"
    ) +
    scale_fill_manual(values = c("#D95F02", "#1B9E77"), labels = c("Untreated", "Treated")) +
    labs(y = "Fungal ASV richness", fill = "Emamectin benzoate") +
    my_gg_theme +
    theme(
        axis.title.x = element_blank(),
        legend.position = "bottom",
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 16)
    )
p1
pdf("figures/ITS/HB_fungi_richness_leaves_litter_Ohor_min.pdf", width = 10, height = 4)
p1
dev.off()

p2 = ggplot(all_div, 
       aes(
           x = factor(
               soil_layer, 
               levels = c("senescent", "Oi", "Oe-Oa", "mineral"),
               labels = c("Senescent ash leaves", "Litter (Oi)", "O horizon (Oe/Oa)", "Mineral (0-10 cm)")
           ),
           y = richness,
           fill = treatment
       )
) +
    geom_boxplot(width = 0.7) +
    scale_fill_manual(values = c("#D95F02", "#1B9E77"), labels = c("Untreated", "Treated")) +
    labs(y = "Fungal ASV richness", fill = "Emamectin benzoate") +
    my_gg_theme +
    theme(
        axis.title.x = element_blank(),
        legend.position = "bottom",
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 16)
    )
p2

p3 = ggplot(all_div, 
            aes(
                x = factor(
                    soil_layer, 
                    levels = c("senescent", "Oi", "Oe-Oa", "mineral"),
                    labels = c("Senescent ash leaves", "Litter (Oi)", "O horizon (Oe/Oa)", "Mineral (0-10 cm)")
                ),
                y = shannon,
                fill = treatment
            )
) +
    geom_boxplot(width = 0.7) +
    scale_fill_manual(values = c("#D95F02", "#1B9E77"), labels = c("Untreated", "Treated")) +
    labs(y = "Shannon diversity (H')", fill = "Emamectin benzoate") +
    my_gg_theme +
    theme(
        axis.title.x = element_blank(),
        legend.position = "bottom",
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 16)
    )
p3

p4 = ggplot(all_div, 
            aes(
                x = factor(
                    soil_layer, 
                    levels = c("senescent", "Oi", "Oe-Oa", "mineral"),
                    labels = c("Senescent ash leaves", "Litter (Oi)", "O horizon (Oe/Oa)", "Mineral (0-10 cm)")
                ),
                y = simpson,
                fill = treatment
            )
) +
    geom_boxplot(width = 0.7) +
    scale_fill_manual(values = c("#D95F02", "#1B9E77"), labels = c("Untreated", "Treated")) +
    labs(y = "Simpson dominance (D)", fill = "Emamectin benzoate") +
    my_gg_theme +
    theme(
        axis.title.x = element_blank(),
        legend.position = "bottom",
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 16)
    )
p4

pdf("figures/ITS/HB_fungi_diversity_leaves_litter_Ohor_min.pdf", width = 9, height = 4)
p2
p3
p4
dev.off()
