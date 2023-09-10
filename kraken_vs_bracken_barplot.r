# Generation of the comparative bar plot for Kraken2 and Bracken results

# Load required libraries
library(ggplot2)
library(reshape2)
library(dplyr)

# Import data
kraken <- read.csv('analysis/classification/otu_table_S.csv', header = T, row.names = 1)
bracken <- read.csv('analysis/plots/bracken_count_matrix.csv', header = T, row.names = 1)

# Create dataframe with sum of counts
kraken <- kraken[, !grepl('SRR', colnames(kraken))]
colnames(kraken) <- sub('X', '', colnames(kraken))
kraken$Kraken <- rowSums(kraken)
kraken.t <- as.data.frame(kraken[, grepl('Kraken', colnames(kraken))])
names.k <- rownames(kraken)
rownames(kraken.t) <- names.k
names(kraken.t)[1] <- 'Kraken'
kraken.t$rowname <- rownames(kraken.t)

bracken <- bracken[, !grepl('SRR', colnames(bracken))]
colnames(bracken) <- sub('X', '', colnames(bracken))
bracken$Bracken <- rowSums(bracken)
bracken.t <- as.data.frame(bracken[, grepl('Bracken', colnames(bracken))])
names.b <- rownames(bracken)
rownames(bracken.t) <- names.b
names(bracken.t)[1] <- 'Bracken'
bracken.t$rowname <- rownames(bracken.t)

# merge both dataframes
merged <- merge(kraken.t, bracken.t, by = 'rowname', all = T)
orgnames <- merged$rowname
merged <- subset(merged, select = -rowname)
rownames(merged) <- orgnames

# Remove NA and do some filtering for visualization
merged[is.na(merged)] <- 0
merged <- merged[rowSums(merged) >= 5000, ]
tp <- head(merged, 15)

# Convert dataframe to long format
tp$OTU <- rownames(tp)
tp.ordered <- arrange(tp, desc(Kraken))
tp_long <- melt(tp.ordered, id.vars = c('OTU'), variable.name = 'method', value.name = 'value')
head(tp_long, 3)

# Plot
gbp <- ggplot(tp_long, aes(x=reorder(OTU, value, decreasing = TRUE), y = value, fill = method)) +
        geom_bar(stat = "identity", position = "dodge", width = 0.7) +
        scale_fill_manual(values = c("Kraken" = "grey", "Bracken" = "black")) +
        scale_x_discrete(limits = tp_long$OTU) +
        xlab("Organisms") +
        ylab("Total") +
        ggtitle("Comparison of Kraken and Bracken totals by organisms") + 
        theme_minimal() +
        theme(axis.text.x = element_text(angle = 45, face = "bold.italic", hjust = 1)) +
        theme(panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank())

ggsave('Kraken_vs_Bracken_bp.png', gbp, width = 11, height = 8.5, bg = 'white', dpi = 300)