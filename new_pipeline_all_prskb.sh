#!/bin/bash
# Base input folder
exec > >(tee -a prs_test.log) 2>&1
#INPUT_DIR="./input_variants_250526/"
INPUT_DIRS=(
  "./testing_new_PRS_settings_my/"
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
                    -q EUR \
                    -v \
                    -i GCST012182 -i GCST007319 -i GCST007511 -i GCST002245 -i GCST001026 \
                    -i GCST004295 -i GCST004296 -i GCST004297 -i GCST004298 -i GCST004299 -i GCST001499 -i GCST006414 -i GCST006061\
                    -i GCST005212 -i GCST006862 -i GCST007799 -i GCST007798 -i GCST000804 -i GCST007993 -i GCST007994 \
                    -i GCST000035 -i GCST001937 -i GCST004950 -i GCST007236 -i GCST004988 -i GCST000678 \
                    -i GCST007856 -i GCST012877 -i GCST012878 -i GCST012879 -i GCST012880 -i GCST007992 -i GCST006131 \
                    -i GCST004184 -i GCST004185 -i GCST007429 -i GCST007692 -i GCST90018807 -i GCST90018807 \
                    -i GCST003116 -i GCST000998 -i GCST005194 -i GCST005195 -i GCST005196 \
                    -i GCST000651 -i GCST000397 -i GCST000649 -i GCST008064 -i GCST008065 \
                    -i GCST007342 -i GCST001877 -i GCST005839 -i GCST010416 \
                    -i GCST000392 -i GCST90014023 -i GCST005536 -i GCST90013445 -i GCST90013446 -i 	GCST000038 \
                    -i GCST009379 -i GCST009380 -i GCST005047 -i GCST000028 -i GCST010553 -i GCST010554 -i GCST010555 -i GCST010556 -i GCST010557 -i GCST006867 -i GCST006868 \
                    -i GCST009722 -i GCST009726 -i GCST011438 -i GCST011439 -i 	GCST011441 -i GCST90011766 -i GCST90011767 -i GCST90011768 -i GCST90011769 -i GCST90011770 -i GCST006395 -i GCST003342 -i GCST003344 \
                    -i GCST011364 -i GCST011365 -i GCST000340 -i GCST001347 \
                    -i GCST000698 \
                    -i GCST000456 -i GCST002991 -i GCST005434 -i GCST000574 -i GCST002553 \
                    -i GCST005567 -i GCST005568 -i GCST005569 -i GCST002318 -i GCST000679 
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