#!/usr/bin/env bash
set -e
BAM="${snakemake_input[bam]}"
CHROM_SIZE="${snakemake_input[chrom_sizes]}"
SLOP="${snakemake_params[slop]}"
BLACKLIST="${snakemake_input[exclude]}"
OUT="${snakemake_output[0]}"
# BED=${BED%.bed}_unfiltered_unsorted.bed
LOG="${snakemake_log[0]}"
# BAM="$1"
# CHROM_SIZE="$2"
# SLOP="$3"
# BLACKLIST="$4"
# OUT="$5"

(bedtools bamtobed -i "$BAM" | \
    awk 'BEGIN {OFS = "\t"} ; {if ($6 == "+") print $1, $2 + 4, $2 + 5, $4, $5, $6; else print $1, $3 - 5, $3 - 4, $4, $5, $6}' | \
    bedtools slop  -i - -g "$CHROM_SIZE" -b "$SLOP" | \
    sort -k1,1 -k2,2n | \
    bedtools intersect -a - -b "$BLACKLIST" -v -sorted ) > "$OUT" 2> "$LOG"
# bedtools sort -i "$BED" | bedtools intersect -a - -b "$BLACKLIST" -v -sorted > "$OUT"
# rm "$BED"
