# Alignment
rule bowtie2_align_PE:
    input:
        fq1="results/trim_galore/{rid}_1.fq.gz",
        fq2="results/trim_galore/{rid}_2.fq.gz",
    output:
        sam=temp("results/bowtie2/all_{rid}.sam")
    log:
        "logs/bowtie2/{rid}_bowtie2.log"
    conda: "../envs/align.yaml"
    threads: 16
    resources:
        cpus_per_task=16,
        mem_mb=lambda wc, attempt: 2**13 * attempt,
        runtime=lambda wc, attempt: 24 * 60 * attempt,
        slurm_partition="cpu_medium"
    params:
        index=igenome_base / genome / "Sequence" / "Bowtie2Index" / "genome",
        extra=config["bowtie2_params"]
    shell:
        "bowtie2 {params.extra} -p {threads} -x {params.index} -1 {input.fq1} -2 {input.fq2} -S {output.sam} 2> {log}"

rule process_bam:
    input:
        "results/bowtie2/all_{rid}.sam"
    output:
        bam=temp("results/bowtie2/marked_{rid}.bam"),
        bai=temp("results/bowtie2/marked_{rid}.bam.bai")
    log:
        "logs/samtools/mark_{rid}.log"
    threads: 8
    conda: "../envs/align.yaml"
    resources:
        cpus_per_task=8,
        mem_mb=lambda wc, input, attempt: min(2 * input.size_mb, 2**(13 + attempt)),
        runtime=lambda wc, attempt: 60*4 * attempt
    shell:
        """
        (samtools view -@ {threads} -u {input} |\
         samtools collate -@ {threads} -u -O - |\
         samtools fixmate -@ {threads} -m /dev/stdin /dev/stdout |\
         samtools sort -@ {threads} -l 0 |\
         samtools markdup -@ {threads} - {output.bam} && samtools index -@ {threads} {output.bam}) 2> {log}
        """

rule filter_alignments:
    input:
        bam="results/bowtie2/marked_{rid}.bam",
        bai="results/bowtie2/marked_{rid}.bam.bai"
    output:
        bam="results/bowtie2/dedup_{rid}.bam",
        bai="results/bowtie2/dedup_{rid}.bam.bai"
    log: "logs/samtools/{rid}_filter.log"
    conda: "../envs/align.yaml"
    threads: 8
    resources:
        cpus_per_task=8,
        mem_mb=lambda wc, attempt: 4096 * attempt,
        runtime=lambda wc, attempt: 60*2 * attempt,
    params:
        extra=config['samtools_filter'],
        chroms=' '.join(chromosomes)
    shell:
        "samtools view -@ {threads} {params.extra} -b -o {output.bam} {input.bam} {params.chroms} && samtools index -@ {threads} {output.bam}"

rule alignment_qc:
    input:
        bam="results/bowtie2/marked_{rid}.bam",
        ref=genome_fasta,
        bai="results/bowtie2/marked_{rid}.bam.bai"
    output:
        multiext("results/picard/{rid}",
                 ".alignment_summary_metrics",
                 ".insert_size_metrics",
                 ".insert_size_histogram.pdf",
                 ".quality_distribution_metrics",
                 ".quality_distribution.pdf",
                 ".quality_by_cycle_metrics",
                 ".quality_by_cycle.pdf",
                 ".base_distribution_by_cycle_metrics",
                 ".base_distribution_by_cycle.pdf",
                 ".gc_bias.detail_metrics",
                 ".gc_bias.summary_metrics",
                 ".gc_bias.pdf",
                 ".quality_yield_metrics")
    log: "logs/picard/{rid}_bam_metrics.log"
    resources:
        mem_mb=lambda wc, attempt: 4096 * attempt,
        runtime=lambda wc, attempt: 60*2 * attempt,
    params:
        extra="--VALIDATION_STRINGENCY LENIENT --METRIC_ACCUMULATION_LEVEL ALL_READS"
    wrapper:
        "v4.5.0/bio/picard/collectmultiplemetrics"

