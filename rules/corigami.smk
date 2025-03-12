rule make_prediction_windows:
    input:
        chrom_sizes=chrom_sizes,
        centrotelo=corigami_genome / "centrotelo.bed"
    output:
        corigami_base/ "data" / "hg38_tiles.bed"  # should be corigami_genome / genome / "tiles.bed"
    log:
        "logs/corigami/tiles.log"
    conda:
        "../envs/bedtools.yaml"
    params:
        width=config["corigami_window"],
        step=config["corigami_step"],
    script:
        "../scripts/make_windows.sh"

rule corigami_predict:
    input:
        windows=corigami_base/ "data" / "hg38_tiles.bed",  # should be corigami_genome / genome / "tiles.bed"
        model=corigami_model,
        seq=corigami_seq,
        ctcf="results/maxatac/{sid}_ctcf_prediction.bw",
        atac="results/maxatac/{sid}_atac_norm.bw",
    output:
        "results/corigami/{sid}_hic_prediction.cool"
    log:
        "logs/corigami/{sid}_predict.log"
    conda:
        "../envs/corigami.yaml"
    params:
        resolution=config["corigami_resolution"],
        window=config["corigami_window"],
        step=config["corigami_step"],
        genome=genome,
        chrom_sizes=chrom_sizes,
        chromosomes=chromosomes,
    resources:
        cpus_per_task=1,
        mem_mb=lambda wc, attempt: 2**12 + attempt * 2**15,  # 32G, 64G, 98G
        runtime=lambda wc, attempt: attempt * 420,  # 7h, 14h, 21h
        slurm_extra="--gres=gpu:1",
        slurm_partition="gpu4_medium"
    script:
        "../scripts/corigami_predict.py"

