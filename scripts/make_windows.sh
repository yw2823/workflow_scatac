#!/usr/bin/env bash

chrom_sizes="${snakemake_input[chrom_sizes]}"
# blacklist="${snakemake_input[blacklist]}"
centrotelo="${snakemake_input[centrotelo]}"
width="${snakemake_params[width]}"
step="${snakemake_params[step]}"
output="${snakemake_output[0]}"
log="${snakemake_log[0]}"

(sort -k1,1 "$chrom_sizes" | \
    bedtools makewindows -g - -w "$width" -s "$step" | \
    bedtools intersect -a stdin -b "$centrotelo" -wa -v -sorted) > "$output" 2> "$log"
