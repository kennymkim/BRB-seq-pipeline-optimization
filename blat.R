#Myscript.R
library(data.table)
library(Matrix)


barcodes <- read.table(file = "BRBseq_Robinson_24_samples.txt", header= TRUE)

dmel_brbseq <- read.table(file = "result.txt", header= TRUE, row.names = 1)
bger_rnaseq <- read.csv(file = "bgermanica_counts.csv", header= TRUE, row.names = 1)


matrix_dir <- "/Users/kkim4/Documents/arthro_molt/rnaseq/brbseq/bger/lib2Solo.out/Gene/raw/"
f <- file(paste0(matrix_dir, "umiDedup-NoDedup.mtx"), "r")
mat <- as.data.frame(as.matrix(readMM(f)))
close(f)
feature.names = fread(paste0(matrix_dir, "features.tsv"), header = FALSE, stringsAsFactors = FALSE, data.table = F)
barcode.names = fread(paste0(matrix_dir, "barcodes.tsv"), header = FALSE, stringsAsFactors = FALSE, data.table = F)
colnames(mat) <- barcode.names$V1
rownames(mat) <- feature.names$V1
fwrite(mat, file = "lib1.counts.txt", sep = "\t", quote = F, row.names = T, col.names = T)

# Define the mapping of sequences to Blat names based on the new dataframe
mapping <- data.frame(
  Name = c("Blat1", "Blat2", "Blat67", "Blat73", "Blat74", "Blat75", "Blat76", "Blat77", "Blat78", "Blat79", "Blat80", "Blat81", "Blat82", "Blat83", "Blat84", "Blat85", "Blat86", "Blat87", "Blat88", "Blat89", "Blat90", "Blat166", "Blat171", "Blat176"),
  B1 = c("ACTACTCGATCTAG", "AACAACCAAGCGGC", "TTGACCTTGCTGGA", "CCGGATGGTGGTAC", "ATGAACTTGGCACA", "ATACCTTAATCTCC", "TACCGGTATGGTCC", "CTCAACCACGTTCC", "TTCCGAAGCCTAAG", "CAACAAGACCTAAG", "CACTTGGACCATCC", "ACTTCCAGGCCTGA", "ACAGCGTACGCGCA", "TACTATGCACGACG", "TTATAACTGGCAAC", "AAGAAGAGTGGACC", "AACACTTGGCCACC", "TAACCGGCCGCAAC", "ACCTTGTCCAATAC", "CCTACGGTGTTCGG", "CTTAAGCTTAAGAG", "TAACCGCAATAGGA", "CACGTATAGATCGC", "CTGGTCTGTCTTCC")
)

# Create the new dataframe that needs column renaming
df_new <- data.frame(
  AACAACCAAGCGGC = 9, AACACTTGGCCACC = 3, AAGAAGAGTGGACC = 2, ACAGCGTACGCGCA = 9, ACCTTGTCCAATAC = 2,
  ACTACTCGATCTAG = 8, ACTTCCAGGCCTGA = 7, ATACCTTAATCTCC = 0, ATGAACTTGGCACA = 2, CAACAAGACCTAAG = 4,
  CACGTATAGATCGC = 14, CACTTGGACCATCC = 0, CCGGATGGTGGTAC = 2, CCTACGGTGTTCGG = 5, CTCAACCACGTTCC = 3,
  CTGGTCTGTCTTCC = 11, CTTAAGCTTAAGAG = 4, TAACCGCAATAGGA = 11, TAACCGGCCGCAAC = 2, TACCGGTATGGTCC = 5,
  TACTATGCACGACG = 3, TTATAACTGGCAAC = 1, TTCCGAAGCCTAAG = 2, TTGACCTTGCTGGA = 8
)

# Create a named vector for renaming
name_map <- setNames(mapping$Name, mapping$B1)

# Rename the columns of df_new using the mapping
colnames(mat) <- sapply(colnames(mat), function(x) if (x %in% names(name_map)) name_map[[x]] else x)
blatella_lib2 <- mat

# Print the updated dataframe to verify changes
print(df_new)

rownames(blat_metadata) <- blat_metadata$Samples
blat_metadata$Samples <- NULL
blatnew_tpm <- log2(wide_data +1)


blatella_allTPM_x  <- blatella_allTPM %>% select(-Blat2_lib1,-Blat2_lib2,-Blat1_lib1,-Blat1_lib2  )

blatella_metadata <- blat_metadata  %>% filter(Sex != "nd")


p_blat <- pca(blatnew_tpm,metadata = blat_metadata2, removeVar = 0.05)
p_blat <- pca(blatella_allTPM_x,metadata = blatella_metadata , removeVar = 0.05)

biplot(p_blat,
       x = 'PC1', y = 'PC2',
       lab = NULL,
       hline = 0, vline = 0,
       legendPosition = 'right',
       showLoadings = FALSE)

biplot(p_blat,
       x = 'PC1', y = 'PC2',
       lab = NULL,
       colby = 'Stage',
       shape = 'Sex',
       hline = 0, vline = 0,
       legendPosition = 'right',
       showLoadings = FALSE)



all_rldcor <- cor(blatnew_tpm)
pheatmap(all_rldcor)
  
  
  
  
# Add the row means as a new column to the data frame (optional)
mat_tpm$row_means <- rowMeans(mat_tpm)

# View the data frame with row means
print("Data Frame with Row Means:")
print(df)

# Filter out rows where the mean is below 1
filtered_df <- blatella_allTPM %>% filter(row_means == 0)



filtered_df <- blatella_allTPM %>% filter(row_means != 0)



fwrite(filtered_df2, file = "blat_nonzero_1MM_all.txt", sep = "\t", quote = F, row.names = T, col.names = T)



# Drop the row_means column if it was added
filtered_df <- filtered_df %>% select(-row_means)

fwrite(filtered_df, file = "blat_nonzero_counts.txt", sep = "\t", quote = F, row.names = T, col.names = T)

2837/25916
12482/28676
4230/27445

13062/28676

blatella_allcounts <- merge(blatella_lib1, blatella_lib2, by = "row.names")


rownames(blatella_allTPM) <- blatella_allTPM$Row.names
blatella_allTPM <- blatella_allTPM[ , -1] # Remove the Row.names column

blatella_allTPM <- subset(blatella_allTPM, select = -"row_means.y")
blatella_allTPM  <- blatella_mat [, !names(blatella_allTPM ) %in% c("row_means.x")]

all_rldcor <- cor(blatella_allTPM)
pheatmap(all_rldcor)

colnames(blatella_allcounts) <- gsub("\\.x$", "_lib1", colnames(blatella_allcounts))
colnames(blatella_allcounts) <- gsub("\\.y$", "_lib2", colnames(blatella_allcounts))



fwrite(blat_metadata, file = "blat_metadata.txt", sep = "\t", quote = F, row.names = F, col.names = T)
write.table(as.data.frame(blat_metadata),file = paste0("blat_metadata.txt"), quote = FALSE,sep = '\t', row.names = FALSE)



blat_metadata2 <- read.table(file = "blat_metadata.txt", header= TRUE, row.names = 1)



# Load necessary library
library(ggplot2)
# Load necessary libraries
library(ggplot2)
library(tidyr)

# Create data frame
data <- data.frame(
  species = c("Lepisosteus_oculatus", "Esox_lucius", "Danio_rerio", "Blatella_germanica", "Drosophila melanogaster"),
  genes_without_counts = c(8431, 9102, 5606, 16194,4032),
  genes_with_counts = c(14884, 22336, 26914, 12482,10837),
  total_number_of_genes = c(23315, 31438, 32520, 28676,14869),
  percentage_detected = c(63.83873043, 71.04777658, 82.76137761, 43.52768866,72.883)
)

# Calculate the percentage of genes without counts
data <- data %>%
  mutate(percentage_without_counts = 100 - percentage_detected)

# Reshape the data for stacked bar plot
data_long <- data %>%
  select(species, percentage_detected, percentage_without_counts) %>%
  pivot_longer(cols = c(percentage_detected, percentage_without_counts),
               names_to = "counts_type", values_to = "percentage")

# Create a stacked bar plot with percentages
ggplot(data_long, aes(x = species, y = percentage, fill = counts_type)) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  labs(title = "Percentage of Genes with and without Counts per Species",
       x = "Species",
       y = "Percentage of Genes") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_manual(values = c("percentage_detected" = "gray", "percentage_without_counts" = "coral2"),
                    name = "Counts Type",
                    labels = c("Total number of genes", "Genes without Counts")) +
  ylim(0, 100) + coord_flip()


# Load necessary libraries
library(ggplot2)
library(tidyr)

# Create data frame
data <- data.frame(
  species = c("Lepisosteus_oculatus", "Esox_lucius", "Danio_rerio", "Blatella_germanica", "Drosophila_melanogaster"),
  Single = c(82.64, 91.52, 78.93, 82.7, 86.65),
  Duplicated = c(1.48, 2.6, 18.03, 11.38, 9.85),
  Missing = c(15.88, 5.88, 3.04, 5.93, 3.5)
)

# Reshape the data for stacked bar plot
data_long <- data %>%
  select(species, Single, Duplicated, Missing) %>%
  pivot_longer(cols = c(Single, Duplicated, Missing),
               names_to = "metric_type", values_to = "percentage")

# Create a stacked bar plot with percentages
ggplot(data_long, aes(x = species, y = percentage, fill = metric_type)) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  labs(title = "Single, Duplicated, and Missing Genes per Species",
       x = "Species",
       y = "Percentage") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_manual(values = c("Single" = "green", "Duplicated" = "orange", "Missing" = "red"),
                    name = "Metric Type",
                    labels = c("Single", "Duplicated", "Missing")) +
  ylim(0, 100)


# Load required library
library(ggplot2)

# Define the data
data <- data.frame(
  species = c("Lepisosteus_oculatus", "Esox_lucius", "Danio_rerio", "Blatella_germanica", "Drosophila_melanogaster"),
  Single = c(82.64, 91.52, 78.93, 82.7, 86.65),
  Duplicated = c(1.48, 2.6, 18.03, 11.38, 9.85),
  Missing = c(15.88, 5.88, 3.04, 5.93, 3.5)
)

# Reshape the data into long format
library(tidyr)
data_long <- pivot_longer(data, cols = c(Single, Duplicated, Missing),
                          names_to = "Status", values_to = "Percentage")

# Plot using ggplot2
ggplot(data_long, aes(x = species, y = Percentage, fill = Status)) +
  geom_bar(stat = "identity") +
  labs(title = "Genome completeness",
       x = "Species", y = "Percentage") +
  theme_minimal() +
  coord_flip() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))



# Load necessary libraries
library(ggplot2)
library(tidyr)

# Create data frame
data <- data.frame(
  species = c("Lepisosteus_oculatus", "Esox_lucius", "Danio_rerio", "Blatella_germanica", "Drosophila_melanogaster"),
  consistent = c(98, 98, 96.62, 49.92, 97.21),
  inconsistent = c(0.95, 0.38, 0.7, 3.93, 0.67),
  contamination = c(0,0,0,0,0),
  unknown = c(0.99,1.62,2.69,46.15,2.12)
)

# Reshape the data for stacked bar plot

data_long <- data %>%
  select(species, consistent, inconsistent, contamination, unknown) %>%
  pivot_longer(cols = c(consistent, inconsistent, contamination, unknown),
               names_to = "metric_type", values_to = "percentage")

# Create a stacked bar plot with percentages
ggplot(data_long, aes(x = species, y = percentage, fill = metric_type)) +
  geom_bar(stat = "identity") +
  labs(title = "Proteome assessment per species",
       x = "Species", y = "Percentage") +
  theme_minimal() +
  scale_fill_manual(values = c("consistent" = "lightyellow2", "contamination" = "darkseagreen1","inconsistent" = "darkslategray4", "unknown" = "darkslategrey"),
                    name = "Metric Type",
                    labels = c("Consistent Lineage Placement", "Contamination", "Inconsistent Lineage Placement","Total unknown"))+
  coord_flip() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))




if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("NOISeq")

]












library(GenomicFeatures)

## make TxDb from GTF file 
txdb <- makeTxDbFromGFF('GCA_003018175.1_Bger_1.1_genomic.gff')

## get intron information
all.introns <- intronicParts(txdb)
109005

summary(intronslen$width)
intronslen <- as.data.frame(all.introns@ranges)
intronslensorted <- sort(intronslen$width)

introndensity <- density(intronslensorted[1:100000])

`plot(introndensity, main="Density Plot of Intron Lengths", xlab="Length", ylab="Density", col="grey38")
abline(v = c(22, 151), col = c("red", "blue"),  lty = 2, lwd = 2)





sorted_width_desc <- sort(data$width, decreasing=TRUE)

top_100_width <- intronslensorted[1:1000]

breaks <- c(0, 1000, 5000, 10000, 50000, Inf)
labels <- c("1-1000", "1001-5000", "5001-10000", "10001-50000", "above 50000")


quantile((intronslen$width), probs = seq(0.5,1,0.01))

33799












