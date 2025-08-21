if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("VariantAnnotation")
library(VariantAnnotation)

# Path to your compressed VCF
vcf_file <- "S20.deepvariant.vcf.gz"
vcf_file <- "S21.deepvariant.vcf.gz"

# Read the file
vcf <- readVcf(vcf_file, genome = "hg38") 

df <- as.data.frame(rowRanges(vcf))
head(df)
