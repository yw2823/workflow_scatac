#!/bin/bash

# Use Snakemake's input and output variables
input_files=("${snakemake_input[@]}")
output_file="${snakemake_output[0]}"

# Check if at least one input file is provided
if [ "${#input_files[@]}" -lt 1 ]; then
    echo "Error: No input files provided."
    exit 1
fi

# Write header = sample_id + attributes
first_file="${input_files[0]}"
header="sample_id"
while IFS=":" read -r attribute _; do
    header+="\t$attribute"
done < "$first_file"
echo -e "$header" > "$output_file"

# Process each input file and append rows to the output file
for file in "${input_files[@]}"; do
    row=$(basename "$file")  # Get the file basename
    row=${row%.*}
    while IFS=":" read -r _ value; do
        row+="\t$(echo "$value" | xargs)"  # Append the value, trimming whitespace
    done < "$file"
    echo -e "$row" >> "$output_file"
done

