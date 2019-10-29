##
##    GSEA_sets_human_to_mouse.R
##
##

## PART 1: 

MGI_human_mouse_homologs <- read.csv(file="input/HMD_HumanPhenotype_160218.csv",sep = ",", header = F)
# Remove columns we don't need
MGI_human_mouse_homologs$V2 = NULL
MGI_human_mouse_homologs$V3 = NULL
MGI_human_mouse_homologs$V5 = NULL
MGI_human_mouse_homologs$V6 = NULL

# Rename columns
names(MGI_human_mouse_homologs) = c("Human_gene", "Mouse_gene")


## PART 2: Do the merging
# TGFbeta1
TGFbeta1 <- read.table(file="input/GSEA_sets/TGFb1_signalling.txt")
TGFbeta1 = merge(x = TGFbeta1, y = MGI_human_mouse_homologs, by.x = c("V1"), by.y =c("Human_gene"), all.x = T)

write.table(TGFbeta1,file="output/tables/Conversions/GSEA_sets/GSEA_TGFbeta1_mouse.txt",row.names=F,col.names=T,quote=F,sep="\t")

