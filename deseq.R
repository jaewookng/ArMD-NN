library(DESeq2)

# Load the raw counts matrix
raw_counts <- read.table("~/Downloads/GSE115828_raw_counts_GRCh38.p13_NCBI.tsv", 
                         header = TRUE, 
                         row.names = 1, 
                         sep = "\t")
View(raw_counts)

#systematically match conditions from geo query
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("GEOquery", force = T)
library(GEOquery)
# Fetch geo query for a specific GEO series
gds = getGEO("GSE115828", GSEMatrix = TRUE)



# Extract MGS levels (minnesota grading system, not microbiome)
characteristics <- pData(phenoData(gds[[1]]))
mgs_levels <- characteristics$characteristics_ch1
characteristics$mgs_level <- as.numeric(gsub("mgs_level: ", "", 
                                             sapply(strsplit(as.character(mgs_levels), ";"), 
                                                    function(x) grep("mgs_level", x, value = TRUE))))

#clean out IDs
colnames(raw_counts) <- gsub("^X", "", colnames(raw_counts))
colnames(raw_counts) <- gsub("\\..*", "", colnames(raw_counts))
rownames(characteristics) <- gsub("^X", "", rownames(characteristics))
rownames(characteristics) <- gsub("\\..*", "", rownames(characteristics))
rownames(characteristics) <- toupper(rownames(characteristics))
colnames(raw_counts) <- toupper(colnames(raw_counts))

# Filter for valid MGS levels (1, 2, 3, 4)
valid_samples <- rownames(characteristics)[characteristics$mgs_level %in% c(1, 2, 3, 4)]

# Align valid_samples with raw_counts
valid_samples <- intersect(valid_samples, colnames(raw_counts))
valid_indices <- match(valid_samples, colnames(raw_counts))

# Subset raw_counts and metadata
filtered_raw_counts <- raw_counts[, valid_indices, drop = FALSE]
filtered_metadata <- characteristics[valid_samples, , drop = FALSE]


# Create metadata for DESeq2
metadata <- data.frame(
  row.names = valid_samples,
  condition = filtered_metadata$mgs_level
)

#DESeq2 pipeline
dds <- DESeqDataSetFromMatrix(
  countData = filtered_raw_counts,
  colData = metadata,
  design = ~ condition
)
dds <- estimateSizeFactors(dds)
dds <- DESeq(dds)
results <- results(dds, alpha = 0.05)
summary(results)

top_genes <- results[order(abs(results$log2FoldChange), decreasing = TRUE), ]
head(top_genes)

library(ggplot2)

# Create a data frame for ggplot
results_df <- as.data.frame(results)
results_df$Significant <- results_df$padj < 0.05

ggplot(results_df, aes(x = log2FoldChange, y = -log10(pvalue), color = Significant)) +
  geom_point() +
  labs(title = "Volcano Plot", x = "Log2 Fold Change", y = "-Log10(p-value)") +
  theme_minimal()

plotMA(results, ylim = c(-2, 2))

ranked_genes <- results[order(results$log2FoldChange, decreasing = TRUE), "log2FoldChange"]

# Load necessary library
library(dplyr)

# Read in the TSV file
annotation_data <- read.csv("~/Downloads/Human.GRCh38.p13.annot.tsv", sep = "\t", header = TRUE)

# Filter the annotation data to include only rows that match the GeneIDs in the list
matched_genes <- annotation_data %>%
  filter(GeneID %in% top_genes@rownames) %>%
  select(GeneID, Symbol)  # You can add other columns like 'Description' if needed

ordered_matched_genes <- matched_genes %>%
  arrange(match(GeneID, top_genes@rownames))
# View the results
print(ordered_matched_genes)
