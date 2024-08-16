

seq_counts = read.csv("data/seq_count_comparison.csv")
seq_counts

pdf("figures/seq_count_comp.pdf")
plot(log10(BAC) ~ log10(FUN), data = seq_counts)
dev.off()

#plot((BAC) ~ (FUN), data = seq_counts)
cor(seq_counts$BAC, seq_counts$FUN, method = "pearson")

cor(log10(seq_counts$BAC), log10(seq_counts$FUN), method = "pearson")
cor.test(log10(seq_counts$BAC), log10(seq_counts$FUN), method = "pearson")
