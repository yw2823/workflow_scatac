if 'Run' in meta.columns:
    # If "Run" column is present download data from SRA
    rule get_fastq_pe_gz:
        output:
            "resources/fastq/{rid}_1.fastq.gz",
            "resources/fastq/{rid}_2.fastq.gz",
        log:
            "logs/fastq/{rid}_pe.gz.log"
        params:
            extra="--skip-technical"
        threads: 6  # defaults to 6
        wrapper:
            "v4.5.0/bio/sra-tools/fasterq-dump"
else:
    # Rename fastq to have standardized names
    rule link_fastq_files:
        input:
            unpack(get_fastqs)
        output:
            fq1="resources/fastq/{rid}_1.fastq.gz",
            fq2="resources/fastq/{rid}_2.fastq.gz",
        shell:
            """
            ln -s "$(realpath {input.fq1})" {output.fq1}
            ln -s "$(realpath {input.fq2})" {output.fq2}
            """

# Preprocess
rule fastqc_pre_trim:
    input:
        "resources/fastq/{rid}_{x}.fastq.gz"
    output:
        html="results/fastqc/{rid}_{x}.html",
        zip="results/fastqc/{rid}_{x}_fastqc.zip"
    log: "logs/fastqc/{rid}_{x}.log"
    threads: 2
    resources:
        mem_mb=4096
    params:
        extra="--quiet"
    wrapper:
        "v4.5.0/bio/fastqc"

rule trim_galore_pe:
    input:
        ["resources/fastq/{rid}_1.fastq.gz", "resources/fastq/{rid}_2.fastq.gz"],
    output:
        fasta_fwd=temp("results/trim_galore/{rid}_1.fq.gz"),
        report_fwd="results/trim_galore/{rid}_1_trimming_report.txt",
        fasta_rev=temp("results/trim_galore/{rid}_2.fq.gz"),
        report_rev="results/trim_galore/{rid}_2_trimming_report.txt"
    threads: 6
    resources:
        cpus_per_task=6,
        mem_mb=lambda wc, input: max(2.5 * input.size_mb, 2**13),
        runtime=lambda wc, attempt: 12 * 60 * attempt,
        slurm_partition=lambda wc, attempt: "cpu_medium" if attempt > 1 else "cpu_short"
    params:
        extra="--illumina -q {}".format(config["trim_q"]),
    log:
        "logs/trim_galore/{rid}.log",
    wrapper:
        "v4.5.0/bio/trim_galore/pe"

rule fastqc_post_trim:
    input:
        "results/trim_galore/{rid}_{x}.fq.gz"
    output:
        html="results/fastqc/{rid}_{x}_trimmed.html",
        zip="results/fastqc/{rid}_{x}_trimmed_fastqc.zip"
    log: "logs/fastqc/{rid}_{x}_trimmed.log"
    threads: 2
    resources:
        mem_mb=4096
    params:
        extra="--quiet"
    wrapper:
        "v4.5.0/bio/fastqc"

