# Post process
rule generate_Tn5_sites_bed:
    input:
        bam="results/bowtie2/dedup_{rid}.bam",
        bai="results/bowtie2/dedup_{rid}.bam.bai",
        chrom_sizes=chrom_sizes,
        exclude=blacklist_bed
    output:
        temp("results/maxatac/{rid}_tn5_sites.bed")
    log:
        "logs/process_atac/{rid}_get_tn5_sites.log"
    conda:
        "../envs/bedtools.yaml"
    resources:
        cpus_per_task=4,
        mem_mb=lambda wc, attempt: attempt * 2**13,
        runtime=lambda wc, attempt: attempt * 2*60,
    params:
        slop=config['slop']
    script:
        "../scripts/get_tn5_sites.sh"

rule Tn5_coverage:
    input:
        bam="results/bowtie2/dedup_{rid}.bam",
        bed="results/maxatac/{rid}_tn5_sites.bed",
        chrom_sizes=chrom_sizes
    output:
        bw=temp("results/maxatac/{rid}_atac_coverage.bw"),
        sf="results/maxatac/{rid}_atac_coverage_scale-factor.txt"
    log:
        "logs/process_atac/{rid}_tn5_coverage.log"
    conda:
        "../envs/coverage.yaml"
    resources:
        cpus_per_task=4,
        mem_mb=lambda wc, attempt: attempt * 2**12,
        runtime=lambda wc, attempt: attempt * 2*60,
    params:
        rpm=config['rpm_sf']
    script:
        "../scripts/tn5_coverage.sh"

rule macs2_call_peaks:
    input:
        treatment="results/maxatac/{rid}_tn5_sites.bed"
    output:
        xls="results/macs2/{rid}_atac_peaks.xls",
        narrow="results/macs2/{rid}_atac_peaks.narrowPeak",
        summits="results/macs2/{rid}_atac_summits.bed"
    log:
        "logs/macs2/{rid}_callpeak.log"
    conda:
        "../envs/macs2.yaml"
    resources:
        mem_mb=lambda wc, attempt: attempt * 2**12,
        runtime=lambda wc, attempt: attempt * 60,
    params:
        gsize='hs' if genome.startswith('hg') else 'mm'
    shell:
        """
        macs2 callpeak -t {input} --name {wildcards.rid}_atac --outdir results/macs2 \
            -g {params.gsize} --nomodel --keep-dup all --shift 0 --extsize 40 --verbose 2 --call-summits &> {log}
        """

# QC
rule calculate_frip:
    input:
        tn5="results/maxatac/{rid}_tn5_sites.bed",
        peaks="results/macs2/{rid}_atac_peaks.narrowPeak"
    output:
        temp("results/macs2/frip_{rid}.tsv")
    log:
        "logs/process_atac/{rid}_calculate_frip.log"
    conda:
        "../envs/bedtools.yaml"
    resources:
        mem_mb=lambda wc, attempt: attempt * 2**12,
        runtime=lambda wc, attempt: attempt * 60,
    script:
        "../scripts/compute_frip.sh"
    # shell: "../scripts/compute_frip.sh {input} {output} 2> {log}"

rule combine_frip:
    input:
        expand("results/macs2/frip_{rid}.tsv", rid=meta['RID'])
    output:
        "results/macs2/atac_frip.tsv"
    shell:
        "echo -e sample_id\\tpeak_counts\\tpeak_length\\tread_counts\\treads_in_peaks\\tfrip > {output} && cat {input} >> {output}"

rule atacqv:
    input:
        bam="results/bowtie2/dedup_{rid}.bam",
        peaks="results/macs2/{rid}_atac_peaks.narrowPeak",
        exclude=blacklist_bed
    output:
        out="results/ataqv/{rid}.ataqv.out",
        metrics="results/ataqv/{rid}.ataqv.json.gz"
    log:
        "logs/ataqv/{rid}_ataqv.log"
    conda:
        "../envs/atacqv.yaml"
    resources:
        mem_mb=lambda wc, attempt: attempt * 2**12,
        runtime=lambda wc, attempt: attempt * 90
    params:
        tss="--tss-file {}".format(config) if 'tss.bed' in config else '',
        org="human" if genome.startswith("hg") else "mouse"
    shell:
        "ataqv --peak-file {input.peaks} --excluded-region-file {input.exclude} {params.tss} "
        "--name {wildcards.rid} --metrics-file {output.metrics} "
        "{params.org} {input.bam} > {output.out} 2> {log}"

rule mkarv:
    input:
        expand("results/ataqv/{rid}.ataqv.json.gz", rid=meta['RID'])
    output:
        "results/ataqv/index.html"
    log:
        "logs/ataqv/mkarv.log"
    conda:
        "../envs/atacqv.yaml"
    threads: 8
    resources:
        cpus_per_task=8,
        mem_mb=lambda wc, attempt: attempt * 2**14,
        runtime=8*60
    shell:
        "mkarv --concurrency {threads} --force results/_summary/ataqv/ {input} 2> {log}"

