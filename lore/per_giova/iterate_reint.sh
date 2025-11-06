#!/bin/bash

# Define the list of lambda values to iterate over
LAMBDA_VALUES=(75)

# Base directory for output
BASE_OUTPUT_DIR="lore/per_giova/models/reint_activity"

# --- Loop Start ---
echo "Starting adabmDCA reintegration loop..."

for lambda in "${LAMBDA_VALUES[@]}"; do
    
    echo "--------------------------------------------------"
    echo "Running for --lambda_ $lambda"
    echo "--------------------------------------------------"
    
    # Define a unique output path for this specific lambda value
    OUTPUT_PATH="${BASE_OUTPUT_DIR}/reint_${lambda}"
    
    # Create the output directory if it doesn't exist (optional, but good practice)
    mkdir -p "$OUTPUT_PATH"
    
    # Run the command
    # Using backslashes (\) to make the long command readable
    adabmDCA reintegrate \
        -o "$OUTPUT_PATH" \
        -l eaDCA \
        -d lore/per_giova/fasta_azeo.txt \
        --model eaDCA \
        --adj lore/per_giova/activity_adjustment_vector_all_nodup.txt \
        --reint lore/per_giova/activity_all_sequences_to_be_reintegrated_nodup.fa \
        --lambda_ "$lambda" \
        --alphabet ACGU-
        
    echo "Completed run for lambda = $lambda. Output in $OUTPUT_PATH"

done

echo "--------------------------------------------------"
echo "All runs finished."