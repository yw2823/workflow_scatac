#!/usr/bin/env bash

output="${snakemake_output[0]}"

# Header for the new column
header=$(head -n 1 "${snakemake_inpupt[0]}")
echo -e "sample_id\t$header" > "$output"

for file in "${snakemake_input[@]}"; do
    sample_id=$(basename "$file")
    sample_id=${sample_id%.*}
    awk -v sample_id="$sample_id" 'NR > 1 {print sample_id "\t" $0}' "$file" >> "$output"
done

