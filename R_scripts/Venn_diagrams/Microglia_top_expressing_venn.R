##
##    Microglia_top_expressing_venn.R
##
##

## PART 1: Load required packages

# Load required packages
library(DESeq2)
library(gplots)
library(limma)
library(VennDiagram)

## PART 2: source the file generating the microglia dataframe. 
source("R_scripts/shared/Datasets_microglia.R")

## PART 3: Setup DEseq2 analysis for all sample groups. Comparing to microglia from controls (CP6)

### DESeq2 Normalization and DE of P5 v P6
colDataD56 <- data.frame(condition=rep(NA,length(colnames(cntNiek_microglia))), row.names=colnames(cntNiek_microglia))
colDataD56$condition[grep("T.P5",rownames(colDataD56))] <- "TP5"
colDataD56$condition[grep("T..P5",rownames(colDataD56))] <- "TP5"
colDataD56$condition[grep("C.P6",rownames(colDataD56))] <- "CP6"
colDataD56$condition[grep("C..P6",rownames(colDataD56))] <- "CP6"
colDataD56$condition[grep("T.P6",rownames(colDataD56))] <- "TP6"
colDataD56$condition[grep("T..P6",rownames(colDataD56))] <- "TP6"
dds <- DESeqDataSetFromMatrix(countData = cntNiek_microglia,colData = colDataD56,design = ~condition)
dds$condition <- factor(dds$condition,levels=c("CP6",
                                               "TP6",
                                               "TP5"))

dds <- DESeq(dds)

rm(colDataD56, cntNiek_microglia)

## PART 4: 
# For comparative and descriptive analysis we use the rlog (see manual for descriptions). blind=FALSE to keep sample group info for calculating normalizations. 
rld56 <- rlog(dds, blind = FALSE)

# Extract the rld transformed values
rld56_values <- assay(rld56)

# Reorder the columns for sake of readability
rld56_values <- rld56_values[,c("MC.P6", "FC1.P6", "FC2.P6", "MT.P6", "FT.P6", "FT1.P6", "MT.P5", "FT.P5", "FT1.P5")]

# Subset per group for easy calculations of the means per group
rld56_values_CP6 <- rld56_values[, c("MC.P6", "FC1.P6", "FC2.P6")]
rld56_values_TP6 <- rld56_values[, c("MT.P6", "FT.P6", "FT1.P6")]
rld56_values_TP5 <- rld56_values[, c("MT.P5", "FT.P5", "FT1.P5")]

# Order on the highest expressed genes (mean per group) for CP6 group
avgCP6 <- rep(NA, nrow(rld56_values_CP6))
for (i in 1:nrow(rld56_values_CP6)) {avgCP6[i] <- mean(rld56_values_CP6[i,])}
rld56_CP6_ordered <- rld56_values[order(-avgCP6),]

# Order on the highest expressed genes (mean per group) for TP6 group
avgTP6 <- rep(NA, nrow(rld56_values_TP6))
for (i in 1:nrow(rld56_values_TP6)) {avgTP6[i] <- mean(rld56_values_TP6[i,])}
rld56_TP6_ordered <- rld56_values[order(-avgTP6),]

# Order on the highest expressed genes (mean per group) for TP5 group
avgTP5 <- rep(NA, nrow(rld56_values_TP5))
for (i in 1:nrow(rld56_values_TP5)) {avgTP5[i] <- mean(rld56_values_TP5[i,])}
rld56_TP5_ordered <- rld56_values[order(-avgTP5),]

# order on average rlog value per gene
CP6_ordered <- as.data.frame(rld56_values_CP6[order(-avgCP6),])
TP6_ordered <- as.data.frame(rld56_values_TP6[order(-avgTP6),])
TP5_ordered <- as.data.frame(rld56_values_TP5[order(-avgTP5),])

rm(avgCP6, avgTP6, avgTP5, i, rld56, rld56_CP6_ordered, rld56_TP6_ordered, rld56_TP5_ordered,
   rld56_values)

# Add column with average rlog value
CP6_ordered$mean_rld = rowMeans(CP6_ordered[,c("MC.P6", "FC1.P6", "FC2.P6")], na.rm=FALSE)
TP6_ordered$mean_rld = rowMeans(TP6_ordered[,c("MT.P6", "FT.P6", "FT1.P6")], na.rm=FALSE)
TP5_ordered$mean_rld = rowMeans(TP5_ordered[,c("MT.P5", "FT.P5", "FT1.P5")], na.rm=FALSE)

# Only keep the column with calculated gene average
CP6_ordered = CP6_ordered[c("mean_rld")]
TP6_ordered = TP6_ordered[c("mean_rld")]
TP5_ordered = TP5_ordered[c("mean_rld")]

rm(rld56_values_CP6, rld56_values_TP6, rld56_values_TP5, dds)


## PART 5: Create the lists of 450 top genes per set
CP6_ordered_genes = row.names(CP6_ordered)
CP6_ordered_genes_top450 = CP6_ordered_genes[0:450]

TP6_ordered_genes = row.names(TP6_ordered)
TP6_ordered_genes_top450 = TP6_ordered_genes[0:450]

TP5_ordered_genes = row.names(TP5_ordered)
TP5_ordered_genes_top450 = TP5_ordered_genes[0:450]

## PART 6: Calculate overlapping regions etc.
set1 = CP6_ordered_genes_top450
set2 = TP6_ordered_genes_top450
set3 = TP5_ordered_genes_top450

area1 = length(set1)
area2 = length(set2)
area3 = length(set3)
n12 = length(intersect(set1, set2))
n23 = length(intersect(set2, set3))
n13 = length(intersect(set1, set3))
n123 = length(Reduce(intersect, list(set1, set2, set3)))

## PART 7: Draw Venn and save file
venn.plot <- draw.triple.venn(
  area1 = area1,
  area2 = area2,
  area3 = area3,
  n12 = n12,
  n23 = n23,
  n13 = n13,
  n123 = n123,
  lwd = rep(1,3),
  rotation = 2,
  category = c("Control", "GFPneg", "GFPpos"),
  fill = c("dodgerblue4", "firebrick2", "forestgreen"),
  lty = "blank",
  alpha = rep(0.5,3),
  cex = rep(1.3,7),
  cat.cex = 1.2,
  cat.col = c("dodgerblue4", "firebrick2", "forestgreen"),
  cat.dist = rep(0.075, 3),
  fontface = as.vector(rep("plain", 7)),
  fontfamily = as.vector(rep("sans", 7)),
  cat.fontface = as.vector(rep("plain", 3)),
  cat.fontfamily = as.vector(rep("sans", 3))
)

pdf(file="output/pdf/Venn/Venn_GFPpos_set_170201_top450.pdf", height = 6, width = 8)
grid.draw(venn.plot)
dev.off()

## PART 8: Extract info on subsets and write tables
only_set1 = setdiff(set1, c(set2, set3))
only_set2 = setdiff(set2, c(set1, set3))
only_set3 = setdiff(set3, c(set1, set2))

write.table(only_set1,file="output/tables/Venn_diagrams/Only_CP6_top450_170201.txt",col.names=T,quote=F,sep="\t")
write.table(only_set2,file="output/tables/Venn_diagrams/Only_TP6_top450_170201.txt",col.names=T,quote=F,sep="\t")
write.table(only_set3,file="output/tables/Venn_diagrams/Only_TP5_top450_170201.txt",col.names=T,quote=F,sep="\t")
