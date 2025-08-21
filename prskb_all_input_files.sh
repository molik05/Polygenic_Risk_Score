#!/bin/bash
# Base input folder
exec > >(tee -a prs_pbs_250526.log) 2>&1
#INPUT_DIR="./input_variants_250526/"
INPUT_DIRS=(
  "./input_variants_250526/"
)

for INPUT_DIR in "${INPUT_DIRS[@]}"; do
    echo "Processing input directory: $INPUT_DIR"

    has_output=0
    missing_output=0
    missing_folders=()

    # Loop through each subdirectory in input_variants
    for folder in "$INPUT_DIR"/*; do
        if [ -d "$folder" ]; then
            # Check if output.tsv already exists
            if [ -f "$folder/output.tsv" ]; then
                echo "Output already exists in $folder, skipping..."
                continue
            fi

            # Find the .vcf.gz file inside this folder (assumes only one matching file)
            vcf_file=$(find "$folder" -type f -name "*_non_ch.vcf.gz" | head -n 1)

            if [ -n "$vcf_file" ]; then
                SECONDS=0
                echo "Processing $vcf_file..."
                ./PrskbCLI_unzip/runPrsCLI.sh \
                    -f "$vcf_file" \
                    -o "$folder/output.tsv" \
                    -c 0.05 \
                    -r hg38 \
                    -p EUR \
                    -h 0.50 \
                    -q EUR \
                    -v \
                    -t "Alzheimer Disease" -t "Alzheimer's Disease" \
                    -t "Adult Asthma" -t "Adult Onset Asthma" \
                    -t "Allergic Disease (Asthma, Hay Fever or Eczema)" \
                    -t "Breast Cancer" -t "Breast Carcinoma" \
                    -t "Colorectal Cancer" \
                    -t "Chronic Obstructive Pulmonary Disease" -t "Chronic Obstructive Pulmonary Disease (Moderate To Severe)" \
                    -t "Chronic Obstructive Pulmonary Disease (Severe)" \
                    -t "Obesity" \
                    -t "Type 1 Diabetes" -t "Type 1 Diabetes Mellitus" \
                    -t "Type 2 Diabetes" -t "Type 2 Diabetes Mellitus" \
                    -t "Psoriasis" \
                    -t "Glaucoma" \
                    -t "Coronary Atherosclerosis" -t "Coronary Artery Disease" \
                    -t "Congestive hearth failure" \
                    -t "Cardiovascular Disease" \
                    -t "Chronic Kidnes Disease" \
                    -t "Hypertesion" \
                    -t "Atrial Fibrillation" -t "Atrial Fibrillation/Atrial Flutter" -t "Atrial Fibrillation and Flutter (Phecode 427.2)" \
                    -t "Major Depressive Disorder" \
                    -t "Acute Myocardial Infarction" -t "Myocardial Infarction" \
                    -t "Osteoporosis" \
                    -t "Pancreatic Carcinoma" -t "Pancreatic Cancer" \
                    -t "Prostate Carcer" -t "Prostate Carcinoma" -t "Cancer Of Prostate (Phecode 185)" \
                    -t "Rheumatoid Arthritis"
                duration=$SECONDS
                mins=$((duration / 60))
                secs=$((duration % 60))
                echo "Time taken for $vcf_file: ${mins} minute(s) and ${secs} second(s)"
            else
                echo "No VCF file found in $folder, skipping..."
            fi
        fi
    done

    for folder in "$INPUT_DIR"/*; do
        if [ -d "$folder" ]; then
            if [ -f "$folder/output.tsv" ]; then
                ((has_output++))
            else
                ((missing_output++))
                missing_folders+=("$(basename "$folder")")
            fi
        fi
    done

    echo "Done processing all folders."
    echo "Folders with output.tsv: $has_output"
    echo "Folders missing output.tsv: $missing_output"
    echo

    if [ ${#missing_folders[@]} -gt 0 ]; then
        echo "Folders missing output.tsv:"
        for f in "${missing_folders[@]}"; do
            echo "  - $f"
        done
        echo
    fi
done