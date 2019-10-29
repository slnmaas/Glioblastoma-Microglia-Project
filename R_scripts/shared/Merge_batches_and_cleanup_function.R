## PART 1: Load function for cleanup

removeLowQualitySets = function(datasets = NA) {
  
  if(!is.data.frame(datasets)) {
    stop("No datasets data.frame supplied, please supply for function to run.")
  } 
  
  htseq_drop_rows <- c("no_feature", "ambiguous", "too_low_aQual", "not_aligned", "alignment_not_unique","ERCC-")
  drop <- rowSums(sapply(htseq_drop_rows, grepl, rownames(datasets)))>0
  datasets <- datasets[!drop,]
  
  # Print dimensions of data.frame
  dim(datasets)
  
  # Keep genes with more than 1 read in at least 2 samples
  keep <- rowSums(datasets>=1)>=2
  table(keep)
  datasets <- datasets[keep,]
  dim(datasets)
  
  # Drop low read count (less than 5 reads) samples
  hist(colSums(datasets>5), breaks=20)
  
  # Print colSums output to spot the samples with less counts
  colSums(datasets>5)
  
  # Drop samples where less than 6000 genes were detected with more than 15 reads
  drop <- colSums(datasets>5)<6000
  table(drop)
  datasets <- datasets[,!drop]
  
  rm(drop, keep, htseq_drop_rows)
  
  return(datasets)
}


## PART 2: Load and merge the datasets

## Prepare batch 1 samples
cntNiek_batch1 <- read.table(file="input/star_genes_erc.counts_microglia.txt")

colnames(cntNiek_batch1) <- strsplit2(colnames(cntNiek_batch1),"_")[,1]

cntNiekAll <- cntNiek_batch1

rm(cntNiek_batch1)
