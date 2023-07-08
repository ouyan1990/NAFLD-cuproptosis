
#引用包
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(ggplot2)
library(circlize)
library(RColorBrewer)
library(dplyr)
library(ComplexHeatmap)

pvalueFilter=0.05    
qvalueFilter=1    

#定义颜色
colorSel="qvalue"
if(qvalueFilter>0.05){
	colorSel="pvalue"
}

setwd("E:\\bioinformation\\cuproptosis\\NAFLD1\\11.KEGG")      
rt=read.table("diff.txt", header=T, sep="\t", check.names=F)     

#Extract the gene name and convert the gene name to gene id
genes=unique(as.vector(rt[,1]))
entrezIDs=mget(genes, org.Hs.egSYMBOL2EG, ifnotfound=NA)
entrezIDs=as.character(entrezIDs)
rt=data.frame(genes, entrezID=entrezIDs)
gene=entrezIDs[entrezIDs!="NA"]        #Delete the gene whose gene id is NA


#kegg enrichment analysis
kk <- enrichKEGG(gene=gene, organism="hsa", pvalueCutoff=1, qvalueCutoff=1)
KEGG=as.data.frame(kk)
KEGG$geneID=as.character(sapply(KEGG$geneID,function(x)paste(rt$genes[match(strsplit(x,"/")[[1]],as.character(rt$entrezID))],collapse="/")))
#According to the set filtering conditions, remarkable enrichment results were obtained
KEGG=KEGG[(KEGG$pvalue<pvalueFilter & KEGG$qvalue<qvalueFilter),]
write.table(KEGG, file="KEGG.txt", sep="\t", quote=F, row.names = F)

#Defines the number of display paths
showNum=30
if(nrow(KEGG)<showNum){
	showNum=nrow(KEGG)
}

#bubble diagram
pdf(file="bubble.pdf", width = 9, height = 7)
dotplot(kk, showCategory=showNum, orderBy="GeneRatio", label_format=130, color=colorSel)
dev.off()




