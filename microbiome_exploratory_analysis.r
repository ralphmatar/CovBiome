# Load required pacckaged
library(compositions)
library(zCompositions)
library(CoDaSeq)
library(ggplot2)
library(PCAtools)
library(ggpubr)

# Import data and adjust col names
d.pro.0 <- read.csv('analysis/plots/bracken_count_matrix.csv', header = T, row.names = 1)
colnames(d.pro.0) <- sub('X', '', colnames(d.pro.0))


# filter usinf codaseq to resduce the dataset retaining taxa that have at least 00 reads, 
# 0.1% abundant in any sample and are in 10% in all samples
d.pro.abund.unordered <- codaSeq.filter(d.pro.0, min.reads = 5000, min.prop = 0.001, min.occurrence = 0.1, samples.by.row = F)
d.pro.abund.unordered <- d.pro.abund.unordered[, colSums(d.pro.abund.unordered != 0) > 0] # remove samples with 0 counts for all taxa

# Replace zero values with pseudo-count
d.pro.abund.unordered <- cmultRepl(d.pro.abund.unordered, label = '0', method = 'CZM', output = 'p-counts')

# add in the names again and sort by abundance
d.names <- rownames(d.pro.abund.unordered)[
order(apply(d.pro.abund.unordered, 1, sum), decreasing=T) ]
# get the taxa in the reduced dataset by name
d.pro.abund_unordered <- d.pro.abund.unordered[d.names,]
d.pro.abund <- as.data.frame(d.pro.abund_unordered)
d.clr.abund <- codaSeq.clr(d.pro.abund, samples.by.row = F)

### top 10 organisms ###
org.cov.10 <- names(head(sort(rowSums(d.clr.abund), decreasing = TRUE), 10))
# create top 10 matrices
cov.10 <- as.data.frame(d.clr.abund[rownames(d.clr.abund) %in% org.cov.10, ])
# plot
bp.cov <- ggplot(cov.10, aes(x=reorder(rownames(cov.10), rowSums(cov.10), decreasing = TRUE), y=rowSums(cov.10))) +
  geom_bar(stat="identity", fill='black', width = 0.5) +
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 45, face = 'bold.italic', size = 12, vjust = 1.25, hjust = 1))+
  labs(x='Taxa', y='CLR Abundace')+
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank())
ggsave('cov_bp.png', bp.cov, width = 11, height = 8.5, dpi = 300, bg = 'white')

# Principal Componant Analysis and biplot ###
pca <- PCAtools::pca(d.clr.abund)
bip <- biplot(pca, showLoadings = TRUE, labSize = 5, pointSize = 3, sizeLoadingsNames = 3, title = 'Biplot COVID-19', axisLabSize = 15)
scree <- screeplot(pca, xlab = NULL, components = getComponents(pca, 1:7), titleLabSize = 15, title = 'Scree plot', axisLabSize = 15)
cov.biscree <- ggarrange(bip, cov, ncol = 2, widths = c(2/3, 1/3))
ggsave('covid_biscree.png', cov.biscree, width = 11, height = 8.5, dpi = 300, bg = 'white')

### Unsupervided clustering ###
# create a color table for the taxa
tax.cols <- data.frame(matrix(ncol = 2, nrow = length(rownames(d.clr.abund))))
colnames(tax.cols) <- c('taxa', 'color')
tax.cols$taxa <- rownames(d.clr.abund)
tax.cols$color <- sample(colors(), length(rownames(tax.cols)))
# generate the distance matrix and cluster as dendrogram
dend <- t(d.clr.abund) %>%
    dist(method = 'euclidian') %>%
    hclust(method = 'ward.D2') %>%
    as.dendrogram()

# now re-order the data to plot the barplot in the same order
d.order <- d.pro.abund[rownames(d.clr.abund) ,hc$order]
d.order.acomp <- acomp(t(d.order))

### add metadata ###
metadata <- read.csv('data/metadata.csv', header = T, sep = '')
metadata$sequencing_successful <- NULL

colnames(metadata)[colnames(metadata) == 'ID'] <- 'Date'
colnames(metadata)[colnames(metadata) == 'collection_date'] <- 'Location'
colnames(metadata)[colnames(metadata) == 'division'] <- 'Variant'
colnames(metadata)[colnames(metadata) == 'lineage'] <- 'Sequenced'

for (i in seq_along(metadata[, 'Sample'])) {
  metadata[i, 'Sample'] <- sub('sample_', '', metadata[i, 'Sample'])
    }

row <- which(metadata$Sample == '518B')
metadata$Sample[row] <- '518b'

metadata.clr <- metadata[metadata$Sample %in% colnames(d.clr.abund), ]

snames <- metadata.clr$Sample
metadata.clr <- metadata.clr[, 2:ncol(metadata.clr)]
colnames(metadata.clr) <- c('date', 'location', 'variant', 'seq')
rownames(metadata.clr) <- snames

### set colors for dendro annotation ###
location.cols <- data.frame(matrix(ncol = 2, nrow = 2))
colnames(location.cols) <- c('location', 'color')
location.cols$location <- c('Addis', 'Jimma')
location.cols$color <- c('yellow', 'blue')

seq.cols <- data.frame(matrix(ncol = 2, nrow = 2))
colnames(seq.cols) <- c('seq', 'color')
seq.cols$seq <- c('Yes', 'No')
seq.cols$color <- c('red', 'black')

metadata.clr <- t(metadata.clr)
metadata.clr <- metadata.clr[, order(colnames(metadata.clr))]
metadata.order <- metadata.clr[rownames(metadata.clr) ,hc$order]
sample.loc <- metadata.order[2, ]
location.colors <- location.cols[match(sample.loc, location.cols[, 1]), 2]
sample.loc.1 <- metadata.order[4, ]
seq.colors <- seq.cols[match(sample.loc.1, seq.cols[, 1]), 2] 

# create legends including metadata and taxa
lnames.artic <- rownames(d.order)
lnames.artic <- c('ARTIC+', 'ARTIC-', lnames.artic)

lnames.location <- rownames(d.order)
lnames.location <- c('Addis', 'Jimma', lnames.location)

df.location <- data.frame(taxa = c('Addis', 'Jimma'), color = c('yellow', 'blue'))
tax.location <- tax.cols
tax.location <- rbind(df.location, tax.location)

df.artic <- data.frame(taxa = c('ARTIC+', 'ARTIC-'), color = c('red', 'black'))
tax.artic <- tax.cols
tax.artic <- rbind(df.artic, tax.artic)

# grouped plots
# w/ location
pdf('analysis/plots/dend_location.pdf', width = 11, height = 7.5, bg = 'white')
layout(matrix(c(1,3,2,3),2,2, byrow=T), widths=c(6,2), height=c(4,4))
par(mar=c(2,1,1,1)+0.1)
dend %>% set('leaves_pch', 15) %>%
    set('labels', NULL) %>%
    set('leaves_cex', 1) %>%
    set('leaves_col', location.colors) %>%
    plot()
barplot(d.order.acomp, legend.text=F, col=tax.cols$color, axisnames=F, border=NA, xpd=T) # barplot
par(mar=c(0,1,1,1)+0.1)
plot(1,2, pch = 1, lty = 1, ylim=c(-20,20), type = "n", axes = FALSE, ann = FALSE) # legend
legend(x="right", legend=lnames.location, col=tax.location$color, lwd=5, cex=.9, bty = 'n' )
dev.off()
# w/ artic
pdf('analysis/plots/dend_artic.pdf', width = 11, height = 7.5, bg = 'white')
layout(matrix(c(1,3,2,3),2,2, byrow=T), widths=c(75,25), height=c(50,50))
par(mar=c(2,1,1,1)+0.1)
dend %>% set('leaves_pch', 15) %>%
    set('labels', NULL) %>%
    set('leaves_cex', 1) %>%
    set('leaves_col', seq.colors) %>%
    plot()
barplot(d.order.acomp, legend.text=F, col=tax.cols$color, axisnames=F, border=NA, xpd=T) # barplot
par(mar=c(0,1,1,1)+0.1)
plot(1,2, pch = 1, lty = 1, ylim=c(-20,20), type = "n", axes = FALSE, ann = FALSE) # legend
legend(x="right", legend=lnames.artic, col=tax.artic$color, lwd=5, cex=.9, bty = 'n' )
dev.off()

### Univariate differences analysis ###
# get clusters and label samples by cluster
clust <- as.matrix(cutree(hc.50, k=2))
clust[clust == 1] <- 'A'
clust[clust == 2] <- 'B'
clust_rownames <- rownames(clust)

for (i in 1:ncol(d.pro.0)) {
  colname <- colnames(d.pro.0)[i]
  if (colname %in% clust_rownames) {
    index <- which(clust_rownames == colname)
    value <- clust[index, 1]
    new_colname <- paste(colname, value, sep = "_")
    colnames(d.pro.0)[i] <- new_colname
  }
}

# generate the dataset by making a data frame of
d.B <- colnames(d.pro.0)[grep("B", colnames(d.pro.0))] # group B
d.A <- colnames(d.pro.0)[grep("A", colnames(d.pro.0))] # group A
d.aldex <- data.frame(d.pro.0[,d.B], d.pro.0[,d.A]) # make a data frame

# make the vector of set membership in the same order as
conds.aldex <- c(rep("A", 166), rep("B", 72))

# generate 128 Dirichlet Monte-Carlo replicates
x <- aldex.clr(d.aldex, conds.aldex, mc.samples=128, verbose=FALSE)
# calculate p values for each replicate and report the mean
x.t <- aldex.ttest(x)
x.e <- aldex.effect(x, verbose = F)
x.all <- data.frame(x.e, x.t)
# get significant taxa
sig <- which(x.all$wi.eBH <= 0.05)

# generate table
data <- x.all[sig, c(4:7, 10,11)]
xtable(data)

# volcano, effect, MAplots
pdf('univariate.pdf', width = 11, height = 8, bg = 'white')
layout(matrix(c(1,2,3,1,2,3),2,3, byrow=T), widths=c(5,2,2), height=c(4,4))
par(mar=c(5,4,4,1)+0.1)
aldex.plot(x.all, test="wilcox", cutoff=0.05, all.cex=0.7, called.cex=1.1,
rare.col="grey", called.col="black")
plot(x.all$effect, x.all$wi.eBH, log="y", pch=19, main="Effect",
cex=0.5, xlab="Effect size", ylab="Expected Benjamini-Hochberg P")
abline(h=0.05, lty=2)
plot(x.all$diff.btw, x.all$wi.eBH, log="y", pch=19, main="Volcano",
cex=0.5, xlab="Difference", ylab="Expected Benjamini-Hochberg P")
abline(h=0.05, lty=2)
dev.off()

# Generate Phyla and Genera
total <- rowSums(d.pro.abund)
totalcount <- sum(total)
percentages <- round(total/totalcount*100, 2)
d.percent <- data.frame(taxa = rownames(d.pro.abund), percentages)
d.percent <- d.percent[, 1:ncol(d.percent)]
rownames(d.percent) <- NULL
acido <- which(d.percent$taxa == '[Acidovorax] ebreus')
d.percent$taxa[acido] <- 'Acidovorax ebreus'
d.percent$Genus <- sapply(strsplit(as.character(d.percent$taxa), " "), "[", 1)
colnames(d.percent)[1] <- 'Taxon'
d.Genus <- aggregate(percentages ~ Genus, d.percent, sum)
d.Genus <- d.Genus[order(-d.Genus$percentages),]
rownames(d.Genus) <- NULL
d.Genus$Phylum <- c("Proteobacteria", "Bacteroidetes", "Proteobacteria", "Proteobacteria", "Actinobacteria", "Proteobacteria", "Actinobacteria", "Proteobacteria", "Proteobacteria", "Proteobacteria", "Proteobacteria", "Proteobacteria", "Proteobacteria", "Proteobacteria", "Actinobacteria", "Firmicutes","Basidiomycota",  "Proteobacteria","Bacteroidetes","Proteobacteria","Proteobacteria")
d.Phylum <- aggregate(percentages ~ Phylum, d.Genus, sum)
d.Phylum <- d.Phylum[order(-d.Phylum$percentages), ]

