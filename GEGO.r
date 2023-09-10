### Load required packages ###
library(DESeq2)
library(clusterProfiler)
library(org.Hs.eg.db)
library(ggplot2)
library(tidyverse)
library(pheatmap)

### Import featureCounts data and filter ###
countdata <- read.table('analysis/fcounts.txt', header = T, row.names = 1)
countdata <- countdata[, 7:ncol(countdata)]
colnames(countdata) <- sub('analysis.mapped.', '', colnames(countdata))
colnames(countdata) <- sub('.sorted.bam', '', colnames(countdata))
countdata <- countdata[, !grepl('SRR', colnames(countdata))] # filter unwanted samples
countdata <- countdata[rowSums(countdata) >= 10, ] # keep genes that have a count more than 10

### comparison between ARTIC+ and ARTIC- total gene count ###
metadata <- read.csv('data/metadata.csv', header = T, sep = '') # import data adjust names
metadata$sequencing_successful <- NULL
colnames(metadata)[colnames(metadata) == 'ID'] <- 'Date'
colnames(metadata)[colnames(metadata) == 'collection_date'] <- 'Location'
colnames(metadata)[colnames(metadata) == 'division'] <- 'Variant'
colnames(metadata)[colnames(metadata) == 'lineage'] <- 'Sequenced'
for (i in seq_along(metadata[, 'Sample'])) {
  metadata[i, 'Sample'] <- sub('sample_', '', metadata[i, 'Sample'])
    }
metadata$Date <- NULL # keep required data 
metadata$Variant <- NULL
countdata <- countdata[, order(colnames(countdata))]
metadata$Sample <- ifelse(substr(metadata$Sample, 1, 1) == "N", metadata$Sample, tolower(metadata$Sample))
meta.dds <- metadata[metadata$Sample %in% colnames(countdata), ]
meta.dds <- meta.dds[order(meta.dds$Sample), ]
meta.dds <- meta.dds[-2, ] # remove duplicates
snames <- meta.dds$Sample
meta.dds$Sample <- NULL
rownames(meta.dds) <- snames


gene_count <- countdata %>% # prepare plot data 
  rownames_to_column(var = "gene") %>%
  gather(key = "sample", value = "count", -gene) %>%
  filter(count > 0) %>% # Only include rows where the count is greater than 0
  group_by(sample) %>%
  summarize(gene_count = n())
plot_data <- gene_count %>%
  left_join(meta.dds, by = c("sample" = "Location"))
gene.counts <- gene_count$gene_count
meta.dds$gene_count <- gene.counts
# plot the data
bp_genes <- ggplot(meta.dds, aes(x = rownames(meta.dds), y = gene_count, fill = Sequenced))+
                geom_bar(stat = 'identity')+
                scale_fill_manual(values = c("Yes" = "red", "No" = "black")) +
                xlab('Samples') +
                ylab('Number of identified genes with featureCounts') +
                scale_x_discrete(breaks = NULL) +
                labs(fill = 'ARTIC enriched') +
                facet_wrap(~Sequenced) +
                ylim(0,4500) +
                theme_minimal()

ggsave('gene_counts.png', bp_genes, width = 11, height = 7.5, dpi = 300, bg = 'white')


### Create metada for DESeq2 ###
cond <- c(rep(1, length(colnames(countdata))))
coldata <- data.frame(row.names = colnames(countdata), cond)
countdata <- countdata + 1 # add a pseudo-count to handle 0 in every gene

### Normalize counts and run DESeq2 ###
dds <- DESeqDataSetFromMatrix(
    countData = countdata, 
    colData = coldata, 
    design = ~ 1,
)
dds <- estimateSizeFactors(dds)
normalized_counts <- counts(dds, normalized = T)
dds <- DESeq(dds, fitType = 'local')

### Get DESeq2 results and filter to keep genes with base mean > or = to 100 ###
res <- results(dds)
sig <- na.omit(res)
sig <- sig[(sig$baseMean >= 100), ]
sig <- sig[order(-sig$baseMean), ]
genes <- rownames(sig)
mapdata <- normalized_counts - 1 # remove the pseudocount
mapdata <- mapdata[genes, ]
mapdata.sclaed <- scale(mapdata)
hmap <- pheatmap(mapdata.scaled, cluster_rows = F, cluster_cols = T, show_colnames = F, cellheight = 10, cellwidth = 3, fontsize = 15)
save_pheatmap_png <- function(x, filename, width=8500, height=3500, res = 300) {
  png(filename, width = width, height = height, res = res)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}
save_pheatmap_png(hmap, 'heatmap.png')

### Gene Ontology ###
GO_results <- enrichGO(gene = genes, OrgDb = 'org.Hs.eg.db', keyType = 'SYMBOL', ont = 'ALL', readable = T)
df_GO <- as.data.frame(GO_results)
write_csv(df_GO, 'go.csv') # downloaded and calculated Enrichment ration in excel 
go <- read.csv('data/go_enra.csv', sep = ';')
go$p.adjust <- df_GO$p.adjust # readjust the p.adjusted values accoriding to original data

b1 <- ggplot(df_GO, aes(x=Description, y=EnrichmentRatio, fill=p.adjust)) +
  ylab('Enrichment Ratio') +
  geom_bar(position = position_dodge() ,stat = 'identity') +
  coord_flip() +
  ggtitle('Gene ontology') +
  scale_y_continuous(limits = c(0, 20)) +
  theme_classic() +
  theme(axis.text = element_text(face="bold", size = 11))

ggsave('GO.png', b1, width = 11, height = 6, dpi = 300)