rule get_2bit_sequence:
    output:
        maxatac_genome / f"{genome}.2bit"
    log:
        "logs/wget_2bit.log"
    params:
        genome=genome,
        outdir=maxatac_genome
    shell:
        "wget -c https://hgdownload.soe.ucsc.edu/goldenPath/{params.genome}/bigZips/latest/{params.genome}.2bit -P {params.outdir} &> {log}"

rule get_maxatac_data:
    output:
        multiext(str(maxatac_genome/genome),
                 ".chrom.sizes",
                 "_maxatac_blacklist.bed",
                 "_maxatac_blacklist.bw"),
        ctcf_model
    log:
        "logs/wget/maxatac_data.log"
    params:
        basedir=maxatac_data
    shell:
        "mkdir -p {params.basedir} && git clone git@github.com:MiraldiLab/maxATAC_data.git {params.basedir}"

rule sort_maxatac_blacklist:
    input:
        maxatac_genome / (genome + "_maxatac_blacklist.bed")
    output:
        blacklist_bed
    log:
        "logs/preprocess/sort_blacklist.log"
    shell:
        "cut -f1-3 {input} | sort -k1,1 -k2,2n > {output} 2> {log}"

def get_bioreplicates(wildcards):
    sid = wildcards.sid
    RID = meta.loc[[sid], 'RID']
    return [f"results/maxatac/{rid}_atac_coverage.bw" for rid in RID]

rule average_bioreplicates:
    input:
        get_bioreplicates
    output:
        temp("results/maxatac/{sid}_atac.bw")
    log:
        "logs/maxatac/{sid}_average.log"
    params:
        Ninputs=lambda wildcards, input: len(input),
        chroms=' '.join(chromosomes),
        chrom_sizes=chrom_sizes,
    resources:
        mem_mb=lambda wc, attempt: attempt * 2**14,
        runtime=lambda wc, attempt: attempt * 60,
        cpus_per_task=4,
    conda:
        "../envs/maxatac.yaml"
    shell:
        """
        if [[ {params.Ninputs} -gt 1 ]]; then
            maxatac average -i {input} -n {wildcards.sid}_atac -o results/maxatac -cs {params.chrom_sizes} -c {params.chroms}
        elif [[ {input[0]} != {output} ]]; then
            mv {input[0]} {output}
        fi
        """

rule normalize_signal:
    input:
        bw="results/maxatac/{sid}_atac.bw",
        blacklist=blacklist_bw,
        chrom_sizes=chrom_sizes
    output:
        bw="results/maxatac/{sid}_atac_norm.bw",
        chrom_ranges=temp("results/maxatac/{sid}_atac_norm_chromosome_min_max.txt"),
        genome_stats=temp("results/maxatac/{sid}_atac_norm_genome_stats.txt")
    log:
        "logs/maxatac/{sid}_normalize.log"
    params:
        chroms=' '.join(chromosomes)
    resources:
        mem_mb=lambda wc, attempt: attempt * 2**15,
        runtime=lambda wc, attempt: attempt * 60,
        cpus_per_task=4,
    conda:
        "../envs/maxatac.yaml"
    shell:
        """
        maxatac normalize \
            --signal {input.bw} \
            --output results/maxatac \
            --prefix {wildcards.sid}_atac_norm \
            --blacklist_bw {input.blacklist} \
            --chrom_sizes {input.chrom_sizes} \
            --chroms {params.chroms} \
            --method min-max \
            --max_percentile 99 &> {log}
        """

rule merge_chrom_stats:
    input:
        expand("results/maxatac/{sid}_atac_norm_chromosome_min_max.txt", sid=SID)
    output:
        "results/maxatac/atac_norm_chromosome_min_max.txt"
    script:
        "../scripts/merge_chromosomes_stats.sh"

rule merge_genome_stats:
    input:
        expand("results/maxatac/{sid}_atac_norm_genome_stats.txt", sid=SID)
    output:
        "results/maxatac/atac_norm_genome_stats.txt"
    script:
        "../scripts/merge_genome_stats.sh"

rule maxatac_predict:
    input:
        bw="results/maxatac/{sid}_atac_norm.bw",
        model=ctcf_model,
        seq=maxatac_genome / f"{genome}.2bit",
        blacklist=blacklist_bed,
        chrom_sizes=chrom_sizes
    output:
        bw="results/maxatac/{sid}_ctcf_prediction.bw",
    log:
        "logs/maxatac/{sid}_predict.log"
    conda:
        "../envs/maxatac.yaml"
    params:
        chroms=' '.join(chromosomes)
    resources:
        mem_mb=lambda wc, attempt: 45000 + attempt*25000,  # 70Gb, 90G, 120G
        runtime=lambda wc, attempt: attempt * 150,
        cpus_per_task=8,
    shell:
        """
        maxatac predict -m {input.model} \
            --signal {input.bw} \
            -o results/maxatac \
            --prefix {wildcards.sid}_ctcf_prediction \
            --blacklist {input.blacklist} --sequence {input.seq} \
            --chromosome_sizes {input.chrom_sizes} \
            --chromosomes {params.chroms}
        """

