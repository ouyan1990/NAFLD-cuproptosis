

install.packages('e1071')

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("preprocessCore")


inputFile="merge.txt"     

source("geoCRG11.CIBERSORT.R")      


outTab=CIBERSORT("ref.txt", inputFile, perm=1000, QN=T)


outTab=outTab[outTab[,"P-value"]<0.05,]
outTab=as.matrix(outTab[,1:(ncol(outTab)-3)])
outTab=rbind(id=colnames(outTab),outTab)
write.table(outTab, file="CIBERSORT-Results.txt", sep="\t", quote=F, col.names=F)


