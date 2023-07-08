
install.packages("colorspace")
install.packages("stringi")
install.packages("ggplot2")
install.packages("circlize")
install.packages("RColorBrewer")

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("org.Hs.eg.db")
BiocManager::install("DOSE")
BiocManager::install("clusterProfiler")
BiocManager::install("enrichplot")
BiocManager::install("ComplexHeatmap")



#引用包
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(ggplot2)
library(circlize)
library(RColorBrewer)
library(dplyr)
library(ComplexHeatmap)

pvalueFilter=0.05     #p值过滤条件
qvalueFilter=1        #矫正后的p值过滤条件

#定义颜色
colorSel="qvalue"
if(qvalueFilter>0.05){
	colorSel="pvalue"
}

setwd("E:\\bioinformation\\cuproptosis\\NAFLD1\\010.GO")      
rt=read.table("diff.txt", header=T, sep="\t", check.names=F)    

#Extract the gene name and convert the gene name to gene id
genes=unique(as.vector(rt[,1]))
entrezIDs=mget(genes, org.Hs.egSYMBOL2EG, ifnotfound=NA)
entrezIDs=as.character(entrezIDs)
gene=entrezIDs[entrezIDs!="NA"]       #Delete the gene whose gene id is NA


#GO enrichment analysis
kk=enrichGO(gene=gene, OrgDb=org.Hs.eg.db, pvalueCutoff=1, qvalueCutoff=1, ont="all", readable=T)
GO=as.data.frame(kk)
#According to the set filtering conditions, significant enrichment results are obtained
GO=GO[(GO$pvalue<pvalueFilter & GO$qvalue<qvalueFilter),]
#output the Go analysis results
write.table(GO, file="GO.txt", sep="\t", quote=F, row.names = F)
		
#Bubble diagram
pdf(file="bubble.pdf", width=9, height=7)
bub=dotplot(kk, showCategory=10, orderBy="GeneRatio",font.size = 9, label_format=60, split="ONTOLOGY", color=colorSel) + facet_grid(ONTOLOGY~., scale='free')
print(bub)
dev.off()


