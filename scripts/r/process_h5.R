# Load required libraries
library(Seurat)
library(Signac)
library(GenomicRanges)
library(ggplot2)
library(patchwork)
library(AnnotationHub)
library(optparse)
#library(tidyverse)

#1 set up
# Get input/output paths from Snakemake
option_list = list(
  make_option("--input", type="character", help="Input H5 file"),
  make_option("--fragment", type="character", help="Fragment file"),
  make_option("--output", type="character", help="Output directory"),
  make_option("--summary", type="character", help="Analysis summary")
)
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

# Get input/output paths
input_h5 <- opt$input
input_fragment <- opt$fragment
output_dir <- opt$output
done_file <- opt$summary


# Create output directory if it doesn't exist
out_dir_plot <- file.path(output_dir,"plot")
out_dir_barcode <- file.path(output_dir,"barcode")
dir.create(out_dir_plot, showWarnings = FALSE, recursive = TRUE)
dir.create(out_dir_barcode, showWarnings = FALSE, recursive = TRUE)

# Print paths for debugging
cat("Input H5:", input_h5, "\n")
cat("Input fragment:", input_fragment, "\n")
cat("Output directory:", output_dir, "\n")
cat("Summary file:", done_file, "\n")
cat("Output plot:", out_dir_plot, "\n")
cat("Cluster-specific barcode:", out_dir_barcode, "\n")


#2 create object
# Read the H5 file

# Step 1: Read data and create initial objects
tryCatch({
    cat("Reading H5 file...\n")
    counts <- Read10X_h5(filename = opt$input)
    cat("Creating chromatin assay...\n")
    chrom_assay <- CreateChromatinAssay(
        counts = counts,
        min.cells = 10,
        min.features = 200,
        sep = c(":","-"),
        fragments = opt$fragment
    )
    cat("Creating Seurat object...\n")
    patient <- CreateSeuratObject(
        assay = "peaks",
        counts = chrom_assay
    )
    rm(counts, chrom_assay)
    gc()
    cat("Initial objects created successfully\n")
    cat("Shape of seurat object: ", nrow(patient), " features by ", ncol(patient), " cells\n")
}, error = function(e) {
    cat("Error in data reading/object creation:", conditionMessage(e), "\n")
    quit(status = 1)
})


# Add gene annotations


tryCatch({
  cat("Loading ensdb from:", opt$ensdb, "\n")
  #ensdb_v109 <- readRDS(opt$ensdb)
  
  ah <- AnnotationHub()
  ensdb_v109 <- ah[["AH109606"]]
  cat("Creating annotations...\n")
  annotations <- GetGRangesFromEnsDb(ensdb = ensdb_v109)
  rm(ensdb_v109)
  gc()
  seqlevels(annotations) <- paste0('chr', seqlevels(annotations))
  genome(annotations) <- "hg38"
  Annotation(patient) <- annotations
  rm(annotations)
  gc()
  cat("Annotations added successfully\n")
}, error = function(e) {
  cat("Error in annotation process:", conditionMessage(e), "\n")
  quit(status = 1)
})

#3. Quality control
patient <- NucleosomeSignal(object = patient)
patient <- TSSEnrichment(object = patient)

patient$blacklist_ratio <- FractionCountsInRegion(
  object = patient, 
  assay = 'peaks',
  regions = blacklist_hg38_unified
)

patient$nucleosome_group <- ifelse(patient$nucleosome_signal > 4, 'NS > 4', 'NS < 4')


## adding a column of total numbers of fragments in each cell to the meta data of the seurat object
total_fragments <- CountFragments(input_fragment)
#head(total_fragments)
rownames(total_fragments) <- total_fragments$CB
#head(total_fragments)
patient$fragment <- total_fragments[colnames(patient),"frequency_count"]

#FRiP
peak_matrix <- FeatureMatrix(fragments = Fragments(patient),
                             features = granges(patient))

patient$FRiP <- colSums(peak_matrix)/patient$fragment



# Subsetting based on QC
patient.1 <- subset(
  x = patient,
  subset = nCount_peaks > quantile(patient$nCount_peaks,0.05) &
    nCount_peaks < quantile(patient$nCount_peaks,0.95) &
    FRiP > quantile(patient$FRiP,0.1) &
    blacklist_ratio < quantile(patient$blacklist_ratio,0.9) &
    nucleosome_signal < 4 &
    TSS.enrichment > quantile(patient$TSS.enrichment,0.1)
)

#Plotting
# Plot TSS enrichment vs peak counts
pdf(file.path(out_dir_plot, "tss_vs_peaks.pdf"), width=8, height=6)
print(DensityScatter(patient, x = 'nCount_peaks', y = 'TSS.enrichment', log_x = TRUE, quantile=TRUE))
dev.off()

# Plot nucleosome signal
pdf(file.path(out_dir_plot, "nucleosome_signal.pdf"), width=8, height=6)
print(FragmentHistogram(object = patient, group.by = 'nucleosome_group'))
dev.off()


# Plot QC metrics before filtering
pdf(file.path(out_dir_plot, "qc_metrics_before_filtering.pdf"), width=15, height=10)
VlnPlot(
  object = patient,
  features = c('nCount_peaks', 'TSS.enrichment','blacklist_ratio', 'nucleosome_signal', 'FRiP'),
  pt.size = 0.1,
  ncol = 5
)
dev.off()

pdf(file.path(out_dir_plot, "qc_metrics_after_filtering.pdf"), width=15, height=10)
VlnPlot(
  object = patient.1,
  features = c('nCount_peaks', 'TSS.enrichment', 'blacklist_ratio','nucleosome_signal', 'FRiP'),
  pt.size = 0.1,
  ncol = 5
)
dev.off()


#Set seed
set.seed(42)


# Normalization and clustering
patient.1 <- RunTFIDF(patient.1)
patient.1 <- FindTopFeatures(patient.1, min.cutoff = 'q0')
patient.1 <- RunSVD(patient.1)
patient.1 <- RunUMAP(object = patient.1, reduction = 'lsi', dims = c(2,4:30),seed.use =42)
patient.1 <- FindNeighbors(object = patient.1, reduction = 'lsi', dims = c(2,4:30))
patient.1 <- FindClusters(object = patient.1, verbose = FALSE, algorithm = 3, random.seed = 42)

saveRDS(patient.1, file = file.path(output_dir, "seurat_object.rds"))

# Plot UMAP
pdf(file.path(out_dir_plot, "umap_clusters.pdf"), width=10, height=8)
print(DimPlot(object = patient.1, label = TRUE) + 
        ggtitle("UMAP visualization of clusters"))
dev.off()


# Write cluster assignments to separate files
clusters <- unique(patient.1$seurat_clusters)
for(cluster in clusters) {
  cells <- colnames(patient.1)[patient.1$seurat_clusters == cluster]
  write.csv(
    cells,
    file = file.path(out_dir_barcode, paste0("cluster_", cluster, ".csv")),
    row.names = FALSE
  )
}

writeLines(
  c(paste("Analysis completed with", length(clusters), "clusters"), 
  paste("Total cells:", ncol(patient)),
  paste("Cells passing QC:", ncol(patient.1)),
  paste("Median peaks per cell:", median(patient.1$nCount_peaks))
  ),done_file)