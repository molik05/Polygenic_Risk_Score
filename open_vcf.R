# Install Bioconductor manager if you don't have it
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

# Install VariantAnnotation
BiocManager::install("VariantAnnotation")

# Load library
library(VariantAnnotation)

# Path to your compressed VCF
vcf_file <- "S20.deepvariant.vcf.gz"
vcf_file <- "S21.deepvariant.vcf.gz"

# Read the file
vcf <- readVcf(vcf_file, genome = "hg38")  # or hg38 depending on your data

# Inspect
vcf

df <- as.data.frame(rowRanges(vcf))
head(df)
