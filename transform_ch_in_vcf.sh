#!/bin/bash
# Base input folder
INPUT_DIR="./Novaseq_data_PRS/test_gvcf"
# Loop through each subdirectory in input_variants
for folder in "$INPUT_DIR"/*; do
    if [ -d "$folder" ]; then
        # Find the .vcf.gz file inside this folder (assumes only one matching file)
        vcf_file=$(find "$folder" -type f -name "*g.vcf.gz" | head -n 1)

        if [ -n "$vcf_file" ]; then
            echo "Processing $vcf_file..."
            zcat "$vcf_file" | awk 'BEGIN{OFS="\t"} /^#/ {print; next} {sub(/^chr/, "", $1); print}' | bgzip > "${vcf_file%.g.vcf.gz}_non_ch.g.vcf.gz" 

        else
            echo "No VCF file found in $folder, skipping..."
        fi
    fi
done

echo "Done processing all folders."
