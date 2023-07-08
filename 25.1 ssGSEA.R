if (!requireNamespace("BiocManager", quietly = TRUE))
   install.packages("BiocManager")
BiocManager::install("limma")
BiocManager::install("GSEABase")
BiocManager::install("GSVA")

install.packages("ggpubr")

library(reshape2)
library(ggpubr)
library(limma)
library(GSEABase)
library(GSVA)



#Read the expression input file and organize the input file
rt=read.table("merge.txt", header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data)

#Read the gene set file
geneSets=getGmt("immune.gmt", geneIdType=SymbolIdentifier())

#ssGSEA analysis
ssgseaScore=gsva(data, geneSets, method='ssgsea', kcdf='Gaussian', abs.ranking=TRUE)
#The ssGSEA score was corrected
normalize=function(x){
  return((x-min(x))/(max(x)-min(x)))}
ssgseaScore=normalize(ssgseaScore)

#The ssGSEA scoring results were output
ssgseaOut=rbind(id=colnames(ssgseaScore), ssgseaScore)
write.table(ssgseaOut,file="ssGSEA.result.txt",sep="\t",quote=F,col.names=F)


