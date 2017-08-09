# Inigo Martincorena - 07.08.2017
# Calculating AIC for each cancer type under the 192 and 12 substitution models

library("seqinr")
library("Biostrings")
library("MASS")
library("GenomicRanges")
library("dndscv")

tissues = c("ACC","BLCA","BRCA","CESC","COREAD","ESCA","GBM","HNSC","KIRC","KICH","KIRP","LAML","LGG","LIHC","LUAD","LUSC","MESO","OV","PAAD","PCPG","PRAD","SARC","SKCM","STAD","TGCT","THCA","THYM","UCEC","UCS")
setwd("/nfs/team78pc18/im3/Projects/Selection_paper/Referee_response")

AICvals = array(NA, dim=c(length(tissues),2))
rownames(AICvals) = tissues
colnames(AICvals) = c("12r_3w","192r_3w")

for (j in 1:nrow(AICvals)) {
    for (h in 1:ncol(AICvals)) {
        mutations = read.table(sprintf("/nfs/team78pc20/im3/Projects/Selection_paper/Driver_discovery/DATASET_20160727/%s/dNdScv/Normal_panel/%s.swfiltered.5cols", rownames(AICvals)[j], rownames(AICvals)[j]), header=1, sep="\t", stringsAsFactors=F)
        dndsout = dndscv(mutations, sm=colnames(AICvals)[h])
        AICvals[j,h] = AIC(dndsout$poissmodel)
    }
    print(j)
}

write.table(AICvals, file="AIC_values_cancer_types.txt", col.names=T, row.names=T, sep="\t", quote=F)
