##
##    DE_analysis_microglia_V2.R
##

## PART 1: Load required packages

# Load required packages
library(DESeq2)
library(gplots)
library(limma)

## PART 2: source the file generating the dataframe. We use the minus get microglia datasets file
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

## PART 4: Do differential expression analysis for groups TP6 vs CP6, TP5 vs CP6, TP5 vs TP6

# Do differential expression analysis

# Group 1: compare Tumor microglia GFP- (TP6), to Control microglia GFP- (CP6)
dds_resT6C6 <- results(dds, c("condition","TP6","CP6"))

# Group 2: compare Tumor microglia GFP+ (TP5), to Control microglia GFP- (CP6)
dds_resT5C6 <- results(dds, c("condition","TP5","CP6"))

# Group 3: compare Tumor microglia GFP+ (TP5), to Tumor microglia GFP- (TP6)
dds_resT5T6 <- results(dds, c("condition","TP5","TP6"))

rm(dds, colDataD56, cntNiek_microglia)


## PART 5: Extract data for specific gene groups from the DE data obtained in part 4

# load groups of genes as currated in the Groups_for_analysis.R file
source("R_scripts/shared/Groups_for_analysis.R")

# First order the datasets on row names so that output is always alphabetically
dds_resT6C6Ord = dds_resT6C6[order(row.names(dds_resT6C6)),]
dds_resT5C6Ord = dds_resT5C6[order(row.names(dds_resT5C6)),]
dds_resT5T6Ord = dds_resT5T6[order(row.names(dds_resT5T6)),]


# purinergic_P2rx genes
purinergic_P2rx_T6C6Ord = subset(dds_resT6C6Ord, rownames(dds_resT6C6Ord) %in% purinergic_P2rx)
write.table(purinergic_P2rx_T6C6Ord,file="output/tables/DE_microglia/subsets/purinergic_P2rx/DESeq_TP6-CP6_purinergic_P2rx_160112.txt",col.names=T,quote=F,sep="\t")

purinergic_P2rx_T5C6Ord = subset(dds_resT5C6Ord, rownames(dds_resT5C6Ord) %in% purinergic_P2rx)
write.table(purinergic_P2rx_T5C6Ord,file="output/tables/DE_microglia/subsets/purinergic_P2rx/DESeq_TP5-CP6_purinergic_P2rx_160112.txt",col.names=T,quote=F,sep="\t")

purinergic_P2rx_T5T6Ord = subset(dds_resT5T6Ord, rownames(dds_resT5T6Ord) %in% purinergic_P2rx)
write.table(purinergic_P2rx_T5T6Ord,file="output/tables/DE_microglia/subsets/purinergic_P2rx/DESeq_TP5-TP6_purinergic_P2rx_160112.txt",col.names=T,quote=F,sep="\t")

# purinergic_P2ry genes
purinergic_P2ry_T6C6Ord = subset(dds_resT6C6Ord, rownames(dds_resT6C6Ord) %in% purinergic_P2ry)
write.table(purinergic_P2ry_T6C6Ord,file="output/tables/DE_microglia/subsets/purinergic_P2ry/DESeq_TP6-CP6_purinergic_P2ry_160112.txt",col.names=T,quote=F,sep="\t")

purinergic_P2ry_T5C6Ord = subset(dds_resT5C6Ord, rownames(dds_resT5C6Ord) %in% purinergic_P2ry)
write.table(purinergic_P2ry_T5C6Ord,file="output/tables/DE_microglia/subsets/purinergic_P2ry/DESeq_TP5-CP6_purinergic_P2ry_160112.txt",col.names=T,quote=F,sep="\t")

purinergic_P2ry_T5T6Ord = subset(dds_resT5T6Ord, rownames(dds_resT5T6Ord) %in% purinergic_P2ry)
write.table(purinergic_P2ry_T5T6Ord,file="output/tables/DE_microglia/subsets/purinergic_P2ry/DESeq_TP5-TP6_purinergic_P2ry_160112.txt",col.names=T,quote=F,sep="\t")

# ccr_family genes
ccr_family_T6C6Ord = subset(dds_resT6C6Ord, rownames(dds_resT6C6Ord) %in% ccr_family)
write.table(ccr_family_T6C6Ord,file="output/tables/DE_microglia/subsets/ccr_family/DESeq_TP6-CP6_ccr_family_160112.txt",col.names=T,quote=F,sep="\t")

ccr_family_T5C6Ord = subset(dds_resT5C6Ord, rownames(dds_resT5C6Ord) %in% ccr_family)
write.table(ccr_family_T5C6Ord,file="output/tables/DE_microglia/subsets/ccr_family/DESeq_TP5-CP6_ccr_family_160112.txt",col.names=T,quote=F,sep="\t")

ccr_family_T5T6Ord = subset(dds_resT5T6Ord, rownames(dds_resT5T6Ord) %in% ccr_family)
write.table(ccr_family_T5T6Ord,file="output/tables/DE_microglia/subsets/ccr_family/DESeq_TP5-TP6_ccr_family_160112.txt",col.names=T,quote=F,sep="\t")

# cxc_receptors genes
cxc_receptors_T6C6Ord = subset(dds_resT6C6Ord, rownames(dds_resT6C6Ord) %in% cxc_receptors)
write.table(cxc_receptors_T6C6Ord,file="output/tables/DE_microglia/subsets/cxc_receptors/DESeq_TP6-CP6_cxc_receptors_160112.txt",col.names=T,quote=F,sep="\t")

cxc_receptors_T5C6Ord = subset(dds_resT5C6Ord, rownames(dds_resT5C6Ord) %in% cxc_receptors)
write.table(cxc_receptors_T5C6Ord,file="output/tables/DE_microglia/subsets/cxc_receptors/DESeq_TP5-CP6_cxc_receptors_160112.txt",col.names=T,quote=F,sep="\t")

cxc_receptors_T5T6Ord = subset(dds_resT5T6Ord, rownames(dds_resT5T6Ord) %in% cxc_receptors)
write.table(cxc_receptors_T5T6Ord,file="output/tables/DE_microglia/subsets/cxc_receptors/DESeq_TP5-TP6_cxc_receptors_160112.txt",col.names=T,quote=F,sep="\t")

# fc_receptors genes
fc_receptors_T6C6Ord = subset(dds_resT6C6Ord, rownames(dds_resT6C6Ord) %in% fc_receptors)
write.table(fc_receptors_T6C6Ord,file="output/tables/DE_microglia/subsets/fc_receptors/DESeq_TP6-CP6_fc_receptors_160112.txt",col.names=T,quote=F,sep="\t")

fc_receptors_T5C6Ord = subset(dds_resT5C6Ord, rownames(dds_resT5C6Ord) %in% fc_receptors)
write.table(fc_receptors_T5C6Ord,file="output/tables/DE_microglia/subsets/fc_receptors/DESeq_TP5-CP6_fc_receptors_160112.txt",col.names=T,quote=F,sep="\t")

fc_receptors_T5T6Ord = subset(dds_resT5T6Ord, rownames(dds_resT5T6Ord) %in% fc_receptors)
write.table(fc_receptors_T5T6Ord,file="output/tables/DE_microglia/subsets/fc_receptors/DESeq_TP5-TP6_fc_receptors_160112.txt",col.names=T,quote=F,sep="\t")

# ifitms genes
ifitms_T6C6Ord = subset(dds_resT6C6Ord, rownames(dds_resT6C6Ord) %in% ifitms)
write.table(ifitms_T6C6Ord,file="output/tables/DE_microglia/subsets/ifitms/DESeq_TP6-CP6_ifitms_160112.txt",col.names=T,quote=F,sep="\t")

ifitms_T5C6Ord = subset(dds_resT5C6Ord, rownames(dds_resT5C6Ord) %in% ifitms)
write.table(ifitms_T5C6Ord,file="output/tables/DE_microglia/subsets/ifitms/DESeq_TP5-CP6_ifitms_160112.txt",col.names=T,quote=F,sep="\t")

ifitms_T5T6Ord = subset(dds_resT5T6Ord, rownames(dds_resT5T6Ord) %in% ifitms)
write.table(ifitms_T5T6Ord,file="output/tables/DE_microglia/subsets/ifitms/DESeq_TP5-TP6_ifitms_160112.txt",col.names=T,quote=F,sep="\t")

# TLRs genes
TLRs_T6C6Ord = subset(dds_resT6C6Ord, rownames(dds_resT6C6Ord) %in% TLRs)
write.table(TLRs_T6C6Ord,file="output/tables/DE_microglia/subsets/TLRs/DESeq_TP6-CP6_TLRs_160112.txt",col.names=T,quote=F,sep="\t")

TLRs_T5C6Ord = subset(dds_resT5C6Ord, rownames(dds_resT5C6Ord) %in% TLRs)
write.table(TLRs_T5C6Ord,file="output/tables/DE_microglia/subsets/TLRs/DESeq_TP5-CP6_TLRs_160112.txt",col.names=T,quote=F,sep="\t")

TLRs_T5T6Ord = subset(dds_resT5T6Ord, rownames(dds_resT5T6Ord) %in% TLRs)
write.table(TLRs_T5T6Ord,file="output/tables/DE_microglia/subsets/TLRs/DESeq_TP5-TP6_TLRs_160112.txt",col.names=T,quote=F,sep="\t")

# siglecs genes
siglecs_T6C6Ord = subset(dds_resT6C6Ord, rownames(dds_resT6C6Ord) %in% siglecs)
write.table(siglecs_T6C6Ord,file="output/tables/DE_microglia/subsets/siglecs/DESeq_TP6-CP6_siglecs_160112.txt",col.names=T,quote=F,sep="\t")

siglecs_T5C6Ord = subset(dds_resT5C6Ord, rownames(dds_resT5C6Ord) %in% siglecs)
write.table(siglecs_T5C6Ord,file="output/tables/DE_microglia/subsets/siglecs/DESeq_TP5-CP6_siglecs_160112.txt",col.names=T,quote=F,sep="\t")

siglecs_T5T6Ord = subset(dds_resT5T6Ord, rownames(dds_resT5T6Ord) %in% siglecs)
write.table(siglecs_T5T6Ord,file="output/tables/DE_microglia/subsets/siglecs/DESeq_TP5-TP6_siglecs_160112.txt",col.names=T,quote=F,sep="\t")

# MMPs genes
MMPs_T6C6Ord = subset(dds_resT6C6Ord, rownames(dds_resT6C6Ord) %in% MMPs)
write.table(MMPs_T6C6Ord,file="output/tables/DE_microglia/subsets/MMPs/DESeq_TP6-CP6_MMPs_160112.txt",col.names=T,quote=F,sep="\t")

MMPs_T5C6Ord = subset(dds_resT5C6Ord, rownames(dds_resT5C6Ord) %in% MMPs)
write.table(MMPs_T5C6Ord,file="output/tables/DE_microglia/subsets/MMPs/DESeq_TP5-CP6_MMPs_160112.txt",col.names=T,quote=F,sep="\t")

MMPs_T5T6Ord = subset(dds_resT5T6Ord, rownames(dds_resT5T6Ord) %in% MMPs)
write.table(MMPs_T5T6Ord,file="output/tables/DE_microglia/subsets/MMPs/DESeq_TP5-TP6_MMPs_160112.txt",col.names=T,quote=F,sep="\t")

# IHC_genes
IHC_genes_T6C6Ord = subset(dds_resT6C6Ord, rownames(dds_resT6C6Ord) %in% IHC_genes)
write.table(IHC_genes_T6C6Ord,file="output/tables/DE_microglia/subsets/IHC_genes/DESeq_TP6-CP6_IHC_genes_160112.txt",col.names=T,quote=F,sep="\t")

IHC_genes_T5C6Ord = subset(dds_resT5C6Ord, rownames(dds_resT5C6Ord) %in% IHC_genes)
write.table(IHC_genes_T5C6Ord,file="output/tables/DE_microglia/subsets/IHC_genes/DESeq_TP5-CP6_IHC_genes_160112.txt",col.names=T,quote=F,sep="\t")

IHC_genes_T5T6Ord = subset(dds_resT5T6Ord, rownames(dds_resT5T6Ord) %in% IHC_genes)
write.table(IHC_genes_T5T6Ord,file="output/tables/DE_microglia/subsets/IHC_genes/DESeq_TP5-TP6_IHC_genes_160112.txt",col.names=T,quote=F,sep="\t")


# MHC_co_receptors
mhc_co_receptors_T6C6Ord = subset(dds_resT6C6Ord, rownames(dds_resT6C6Ord) %in% mhc_co_receptors)
write.table(mhc_co_receptors_T6C6Ord,file="output/tables/DE_microglia/subsets/mhc_co_receptors/DESeq_TP6-CP6_mhc_co_receptors_160112.txt",col.names=T,quote=F,sep="\t")

mhc_co_receptors_T5C6Ord = subset(dds_resT5C6Ord, rownames(dds_resT5C6Ord) %in% mhc_co_receptors)
write.table(mhc_co_receptors_T5C6Ord,file="output/tables/DE_microglia/subsets/mhc_co_receptors/DESeq_TP5-CP6_mhc_co_receptors_160112.txt",col.names=T,quote=F,sep="\t")

mhc_co_receptors_T5T6Ord = subset(dds_resT5T6Ord, rownames(dds_resT5T6Ord) %in% mhc_co_receptors)
write.table(mhc_co_receptors_T5T6Ord,file="output/tables/DE_microglia/subsets/mhc_co_receptors/DESeq_TP5-TP6_mhc_co_receptors_160112.txt",col.names=T,quote=F,sep="\t")


# mouse_sensome
microglia_sensome_updated_T6C6Ord = subset(dds_resT6C6Ord, rownames(dds_resT6C6Ord) %in% microglia_sensome_updated)
write.table(microglia_sensome_updated_T6C6Ord,file="output/tables/DE_microglia/subsets/mouse_sensome/DESeq_TP6-CP6_microglia_sensome_updated_160505.txt",col.names=T,quote=F,sep="\t")

microglia_sensome_updated_T5C6Ord = subset(dds_resT5C6Ord, rownames(dds_resT5C6Ord) %in% microglia_sensome_updated)
write.table(microglia_sensome_updated_T5C6Ord,file="output/tables/DE_microglia/subsets/mouse_sensome/DESeq_TP5-CP6_microglia_sensome_updated_160505.txt",col.names=T,quote=F,sep="\t")

microglia_sensome_updated_T5T6Ord = subset(dds_resT5T6Ord, rownames(dds_resT5T6Ord) %in% microglia_sensome_updated)
write.table(microglia_sensome_updated_T5T6Ord,file="output/tables/DE_microglia/subsets/mouse_sensome/DESeq_TP5-TP6_microglia_sensome_updated_160505.txt",col.names=T,quote=F,sep="\t")

# mouse_TGFbeta1
GSEA_mouse_TGFbeta1_T6C6Ord = subset(dds_resT6C6Ord, rownames(dds_resT6C6Ord) %in% GSEA_mouse_TGFbeta1)
write.table(GSEA_mouse_TGFbeta1_T6C6Ord,file="output/tables/DE_microglia/subsets/mouse_TGFbeta1/DESeq_TP6-CP6_GSEA_mouse_TGFbeta1_160112.txt",col.names=T,quote=F,sep="\t")

GSEA_mouse_TGFbeta1_T5C6Ord = subset(dds_resT5C6Ord, rownames(dds_resT5C6Ord) %in% GSEA_mouse_TGFbeta1)
write.table(GSEA_mouse_TGFbeta1_T5C6Ord,file="output/tables/DE_microglia/subsets/mouse_TGFbeta1/DESeq_TP5-CP6_GSEA_mouse_TGFbeta1_160112.txt",col.names=T,quote=F,sep="\t")

GSEA_mouse_TGFbeta1_T5T6Ord = subset(dds_resT5T6Ord, rownames(dds_resT5T6Ord) %in% GSEA_mouse_TGFbeta1)
write.table(GSEA_mouse_TGFbeta1_T5T6Ord,file="output/tables/DE_microglia/subsets/mouse_TGFbeta1/DESeq_TP5-TP6_GSEA_mouse_TGFbeta1_160112.txt",col.names=T,quote=F,sep="\t")

# MHC_class_II_genes
MHC_class_II_genes_T6C6Ord = subset(dds_resT6C6Ord, rownames(dds_resT6C6Ord) %in% MHC_class_II_genes)
write.table(MHC_class_II_genes_T6C6Ord,file="output/tables/DE_microglia/subsets/MHC_class_II_genes/DESeq_TP6-CP6_MHC_class_II_genes_170124.txt",col.names=T,quote=F,sep="\t")

MHC_class_II_genes_T5C6Ord = subset(dds_resT5C6Ord, rownames(dds_resT5C6Ord) %in% MHC_class_II_genes)
write.table(MHC_class_II_genes_T5C6Ord,file="output/tables/DE_microglia/subsets/MHC_class_II_genes/DESeq_TP5-CP6_MHC_class_II_genes_170124.txt",col.names=T,quote=F,sep="\t")

MHC_class_II_genes_T5T6Ord = subset(dds_resT5T6Ord, rownames(dds_resT5T6Ord) %in% MHC_class_II_genes)
write.table(MHC_class_II_genes_T5T6Ord,file="output/tables/DE_microglia/subsets/MHC_class_II_genes/DESeq_TP5-TP6_MHC_class_II_genes_170124.txt",col.names=T,quote=F,sep="\t")

# sensome_metabolic
sensome_metabolic_T6C6Ord = subset(dds_resT6C6Ord, rownames(dds_resT6C6Ord) %in% sensome_metabolic)
write.table(sensome_metabolic_T6C6Ord,file="output/tables/DE_microglia/subsets/sensome_metabolic/DESeq_TP6-CP6_sensome_metabolic_170628.txt",col.names=T,quote=F,sep="\t")

sensome_metabolic_T5C6Ord = subset(dds_resT5C6Ord, rownames(dds_resT5C6Ord) %in% sensome_metabolic)
write.table(sensome_metabolic_T5C6Ord,file="output/tables/DE_microglia/subsets/sensome_metabolic/DESeq_TP5-CP6_sensome_metabolic_170628.txt",col.names=T,quote=F,sep="\t")

sensome_metabolic_T5T6Ord = subset(dds_resT5T6Ord, rownames(dds_resT5T6Ord) %in% sensome_metabolic)
write.table(sensome_metabolic_T5T6Ord,file="output/tables/DE_microglia/subsets/sensome_metabolic/DESeq_TP5-TP6_sensome_metabolic_170628.txt",col.names=T,quote=F,sep="\t")

# sensome_misc
sensome_misc_T6C6Ord = subset(dds_resT6C6Ord, rownames(dds_resT6C6Ord) %in% sensome_misc)
write.table(sensome_misc_T6C6Ord,file="output/tables/DE_microglia/subsets/sensome_misc/DESeq_TP6-CP6_sensome_misc_170628.txt",col.names=T,quote=F,sep="\t")

sensome_misc_T5C6Ord = subset(dds_resT5C6Ord, rownames(dds_resT5C6Ord) %in% sensome_misc)
write.table(sensome_misc_T5C6Ord,file="output/tables/DE_microglia/subsets/sensome_misc/DESeq_TP5-CP6_sensome_misc_170628.txt",col.names=T,quote=F,sep="\t")

sensome_misc_T5T6Ord = subset(dds_resT5T6Ord, rownames(dds_resT5T6Ord) %in% sensome_misc)
write.table(sensome_misc_T5T6Ord,file="output/tables/DE_microglia/subsets/sensome_misc/DESeq_TP5-TP6_sensome_misc_170628.txt",col.names=T,quote=F,sep="\t")

# microglia_debris
microglia_debris_T6C6Ord = subset(dds_resT6C6Ord, rownames(dds_resT6C6Ord) %in% microglia_debris)
write.table(microglia_debris_T6C6Ord,file="output/tables/DE_microglia/subsets/microglia_debris/DESeq_TP6-CP6_microglia_debris_170628.txt",col.names=T,quote=F,sep="\t")

microglia_debris_T5C6Ord = subset(dds_resT5C6Ord, rownames(dds_resT5C6Ord) %in% microglia_debris)
write.table(microglia_debris_T5C6Ord,file="output/tables/DE_microglia/subsets/microglia_debris/DESeq_TP5-CP6_microglia_debris_170628.txt",col.names=T,quote=F,sep="\t")

microglia_debris_T5T6Ord = subset(dds_resT5T6Ord, rownames(dds_resT5T6Ord) %in% microglia_debris)
write.table(microglia_debris_T5T6Ord,file="output/tables/DE_microglia/subsets/microglia_debris/DESeq_TP5-TP6_microglia_debris_170628.txt",col.names=T,quote=F,sep="\t")

