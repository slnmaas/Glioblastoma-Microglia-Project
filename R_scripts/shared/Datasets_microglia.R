library(limma)

# Load and merge all the count tables. This file includes the clean_up function used below
source("R_scripts/shared/Merge_batches_and_cleanup_function.R")

# For this selection include the microglia samples coded P5 and P6
datasets <- cntNiekAll[,c(which(grepl("P5",colnames(cntNiekAll))), which(grepl("P6",colnames(cntNiekAll))))]

rm(cntNiekAll)

cntNiek_microglia = removeLowQualitySets(datasets)

rm(datasets, removeLowQualitySets)
