configfile: "config.yaml"
from pathlib import Path

# Workflow params
rule_dir = Path('rules')
genome = config['genome']
igenome_base = Path(config['igenome'])

# maxatac params
genome_fasta   = igenome_base / genome / "Sequence" / "WholeGenomeFasta" / "genome.fa"
bowtie2_index  = igenome_base / genome / "Sequence" / "Bowtie2Index"
maxatac_data   = Path(config["maxatac_data"])
ctcf_model     = maxatac_data / "models" / "CTCF" / config["CTCF_model"]
ctcf_calib     = maxatac_data / "models" / "CTCF" / config["CTCF_calibration"]
maxatac_genome = maxatac_data / genome
chrom_sizes    = maxatac_genome / f'{genome}.chrom.sizes'
blacklist_bed  = maxatac_genome / f'{genome}_maxatac_blacklist_sorted.bed'
blacklist_bw   = maxatac_genome / f'{genome}_maxatac_blacklist.bw'
chromosomes    = config['chromosomes']

# corigami params
corigami_base   = Path(config["corigami_base"])
corigami_model  = corigami_base / "model_weights" / config.get("corigami_model", "corigami_base.ckpt")
corigami_genome = corigami_base / "data" / genome
corigami_seq    = corigami_genome / "dna_sequence"

def get_num_clusters(sample):
    #import subprocess
    #cmd = f"find ./results/{sample}/barcode/ -type f |wc -l"
    #result = subprocess.check_output(cmd,shell=True)
    #return range(int(result.decode().strip()))
    """Get number of cluster files using glob"""
    import glob
    num_files = len(glob.glob(f"./results/{sample}/barcode/*"))
    return range(num_files)

rule all:
    input:
        #rule1
        expand("results/{sample}/clusters_done.txt", sample=config["samples"]),
        #rule2
        expand("results/{sample}/fragments/cluster_{cluster}_fragments.tsv",
                 sample=config["samples"],
                 cluster=get_num_clusters(config["samples"][0])),
        #rule3.1
         expand("maxatac_install.txt",sample=config["samples"]),
        # #rule3
         expand("results/{sample}/maxatac_prepare/cluster{cluster}_IS_slop20_RP20M_minmax01_genome_stats.txt",
                 sample=config["samples"],
                 cluster=get_num_clusters(config["samples"][0])),
        #rule corigami
        expand("results/{sample}/corigami/cluster0_hic_prediction.cool",
                sample=config["samples"],
                cluster=get_num_clusters(config["samples"][0]))

rule process_h5:
    input:
        h5_matrix = "data/{sample}.h5",
        fragment = "data/{sample}_fragment.tsv.gz"
        #ensdb = "data/ensdb_hs_v109.rds"
    output:
        done = "results/{sample}/clusters_done.txt"
    params:
        output_dir = "results/{sample}",
        log_file = "results/{sample}/analysis_log.txt"
    shell:
        """

        module load r/4.4.1


        Rscript scripts/r/process_h5.R \
        --input {input.h5_matrix} \
        --fragment {input.fragment} \
        --output {params.output_dir} \
        --summary {output.done} \
        2> {params.log_file}
        """



rule subset_fragment:
    input:
        log = "results/{sample}/clusters_done.txt",
        fragment = "data/{sample}_fragment.tsv.gz",
        barcode= "results/{sample}/barcode/cluster_{cluster}.csv"
    output:
        "results/{sample}/fragments/cluster_{cluster}_fragments.tsv"
    params:
        out_dir="results/{sample}/fragments",
        log_file = "results/{sample}/subset_log.txt"
    conda:
        "./envs/python.yaml"
    shell:
        """
        mkdir -p {params.out_dir}

        python scripts/python/subset_tsv.py \
        --input {input.fragment} \
        --barcode {input.barcode} \
        --output {output} 2> {params.log_file}
        """

rule setup_maxatac:
    output:
        touch("maxatac_install.txt")
    conda:
        "/gpfs/home/yw2823/.conda/envs/maxatac"
    shell:
        """
        if [!-d"/gpfs/home/yw2823/opt/maxatac"]; then
            echo "MaxATAC data not found. Downloading required data..."
            maxatac data
        else
            echo "MaxATAC data found in /gpfs/home/yw2823/opt/maxatac_data"
        fi
        """
    


#make prediction windows unchanged
#corigami parameters
#chrom_sizes = maxatac_genome / f'{genome}.chrom.sizes'
#corigami_genome = corigami_base / "data" / genome

#corigami
#make prediction windows
rule maxatac_prepare:
    input:
        check_file = "maxatac_install.txt",
        fragments= "results/{sample}/fragments/cluster_{cluster}_fragments.tsv"
    output:
        bed_gz = "results/{sample}/maxatac_prepare/cluster{cluster}_IS_slop20.bed.gz",
        raw_bw= "results/{sample}/maxatac_prepare/cluster{cluster}_slop20_RP20M.bw",
        normal_bw = "results/{sample}/maxatac_prepare/cluster{cluster}_IS_slop20_RP20M_minmax01.bw",
        chrom_stat = "results/{sample}/maxatac_prepare/cluster{cluster}_IS_slop20_RP20M_minmax01_chromosome_min_max.txt",
        genome_stat = "results/{sample}/maxatac_prepare/cluster{cluster}_IS_slop20_RP20M_minmax01_genome_stats.txt"
    params:
        out_dir="results/{sample}/maxatac_prepare",
        log_file = "results/{sample}/maxatac_prepare_log.txt"
    conda:
        "/gpfs/home/yw2823/.conda/envs/maxatac"
    shell:
        """
        mkdir -p {params.out_dir}

        maxatac_prepare.bash
        """



rule corigami_predict:
    input:
        windows="/gpfs/home/yw2823/pipeline-atac2hic/resources/corigami_data/data/hg38_tiles.bed",  
        model='/gpfs/home/rodrij92/PROJECTS/SHARE/epoch=53-step=64260.ckpt',
        seq='/gpfs/home/yw2823/pipeline-atac2hic/resources/corigami_data/data/hg38/dna_sequence',
        ctcf="results/patient3/maxatac_predict/output/cluster0_combined.bigwig",
        atac="results/patient3/maxatac_prepare/cluster0_IS_slop20_RP20M_minmax01.bw"
    output:
        "results/patient3/corigami/cluster0_hic_prediction.cool"
    log:
        "logs/corigami/patient3_cluster0_predict.log"
    conda:
        "/gpfs/home/yw2823/.conda/envs/corigami"
    params:
        resolution=config["corigami_resolution"],
        window=config["corigami_window"],
        step=config["corigami_step"],
        genome='hg38',
        chrom_sizes='/gpfs/home/yw2823/opt/maxatac/data/hg38/hg38.chrom.sizes',
        chromosomes='chr19'
    resources:
        cpus_per_task=1,
        mem_mb=lambda wc, attempt: 2**12 + attempt * 2**15,  # 32G, 64G, 98G
        runtime=lambda wc, attempt: attempt * 420,  # 7h, 14h, 21h
        slurm_extra="--gres=gpu:1",
        slurm_partition="gpu4_medium"
    script:
        "./corigami_predict.py"

include: rule_dir / "fastq_process.smk"
include: rule_dir / "align_reads.smk"
include: rule_dir / "atac_signal.smk"
include: rule_dir / "maxatac_predict.smk"
include: rule_dir / "corigami.smk
    
