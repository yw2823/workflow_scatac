# test_script.R
library(optparse)
library(Seurat)
library(Signac)

# Parse command line arguments
option_list = list(
  make_option("--input", type="character", help="Input H5 file"),
  make_option("--fragment", type="character", help="Fragment file"),
  make_option("--output", type="character", help="Output directory"),
  make_option("--ensdb", type="character", help="Path to ensdb RDS file")
)
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

# Print all arguments
cat("Input file:", opt$input, "\n")
cat("Fragment file:", opt$fragment, "\n")
cat("Output directory:", opt$output, "\n")
cat("Ensdb file:", opt$ensdb, "\n")

# Try to read the H5 file
tryCatch({
  cat("Reading H5 file...\n")
  counts <- Read10X_h5(filename = opt$input)
  cat("H5 file read successfully. Dimensions:", dim(counts), "\n")
  
  cat("Creating chromatin assay...\n")
  chrom_assay <- CreateChromatinAssay(
    counts = counts,
    min.cells = 10,
    min.features = 200,
    sep = c(":","-")
  )
  cat("Chromatin assay created successfully\n")
  
  cat("Loading ensdb...\n")
  ensdb <- readRDS(opt$ensdb)
  cat("Ensdb loaded successfully\n")
  
}, error = function(e) {
  cat("Error occurred:", conditionMessage(e), "\n")
  quit(status = 1)
})