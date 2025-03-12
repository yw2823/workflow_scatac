#!/bin/bash
set -e

if [[ -z "$snakemake_input" ]]; then
    BAM="$1"
    Tn5="$2" # Tn5 bed
    FACTOR="$3" # Coverage normalization factor
    CHROM_SIZES="$4"
    OUT_BW="$5"
    LOG=/dev/stderr
else
    BAM="${snakemake_input[bam]}"
    Tn5="${snakemake_input[bed]}" # Tn5 bed
    FACTOR="${snakemake_params[rpm]}" # Coverage normalization factor
    CHROM_SIZES="${snakemake_input[chrom_sizes]}"
    OUT_BW="${snakemake_output[bw]}"
    LOG="${snakemake_log[0]}"
fi
OUT_BG="${OUT_BW%.bw}".bedgraph
OUT_SF="${OUT_BW%.bw}"_scale-factor.txt

NREADS=$(samtools view -c -F 260 "$BAM")
SF=$(bc -l <<< "${FACTOR} / ${NREADS}")
echo "$SF" > "$OUT_SF"

# unstranded
bedtools genomecov -i "$Tn5" -g "$CHROM_SIZES" -bg -scale "${SF}" > "$OUT_BG" 2> "$LOG"
(bedGraphToBigWig "$OUT_BG" "$CHROM_SIZES" "$OUT_BW" && rm "$OUT_BG") 2>> "$LOG"

# # forward strand
# bedtools genomecov -i "$Tn5" -g "$CHROM_SIZES" -bg -scale "${SF}" -strand + > "$OUT_BG"
# bedGraphToBigWig "$OUT_BG" "$CHROM_SIZES" "$OUT_BWf" && rm "$OUT_BG"

# # reverse strand
# bedtools genomecov -i "$Tn5" -g "$CHROM_SIZES" -bg -scale "${SF}" -strand - > "$OUT_BG"
# bedGraphToBigWig "$OUT_BG" "$CHROM_SIZES" "$OUT_BWr" && rm "$OUT_BG"
