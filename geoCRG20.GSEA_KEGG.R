if (!requireNamespace("BiocManager", quietly = TRUE))
   install.packages("BiocManager")
BiocManager::install("limma")
BiocManager::install("org.Hs.eg.db")
BiocManager::install("DOSE")
BiocManager::install("clusterProfiler")
BiocManager::install("enrichplot")


library(limma)
library(org.Hs.eg.db)
library(clusterProfiler)
library(enrichplot)
         
setwd("E:\\bioinformation\\cuproptosis\\NAFLD1\\15.GSEA")     

##### 1.GSEA analysis of NFE2L2 using KEGG database
#Read the file and organize the input file
rt=read.table("merge.txt", header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data)
data=data[rowMeans(data)>0,]

#The control sample was removed
group=gsub("(.*)\\_(.*)", "\\2", colnames(data))
data=data[,group=="Treat",drop=F]

#The samples were grouped according to the target gene expression, and the logFC of the high-low expression group was obtained
dataL=data[,data["NFE2L2",]<median(data["NFE2L2",]),drop=F]     #Data of low expression group
dataH=data[,data["NFE2L2",]>=median(data["NFE2L2",]),drop=F]    #Data from the high expression group
meanL=rowMeans(dataL)
meanH=rowMeans(dataH)
meanL[meanL<0.00001]=0.00001
meanH[meanH<0.00001]=0.00001
logFC=log2(meanH)-log2(meanL)
#Sequencing genes based on logfc
logFC=sort(logFC, decreasing=T)
genes=names(logFC)

#Read the gene set file
gmt=read.gmt("c2.cp.kegg.symbols.gmt"  )

#GSEA enrichment analysis
kk=GSEA(logFC, TERM2GENE=gmt, pvalueCutoff = 1)
kkTab=as.data.frame(kk)
kkTab=kkTab[kkTab$pvalue<0.05,]
#kkTab=kkTab[kkTab$p.adjust<0.05,]
write.table(kkTab,file="GSEA.result_NFE2L2_KEGG.txt",sep="\t",quote=F,row.names = F)
	
#Plot the GSEA enrichment
termNum=6     #Display the number of paths
if(nrow(kkTab)>=termNum){
	showTerm=row.names(kkTab)[1:termNum]
	gseaplot=gseaplot2(kk, showTerm, base_size=8, title="NFE2L2")
	pdf(file="GSEA_NFE2L2_KEGG.pdf", width=7.5, height=5.5)
	print(gseaplot)
	dev.off()
}



##### 2.GSEA analysis of DLD using KEGG database
#Read the file and organize the input file
rt=read.table("merge.txt", header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data)
data=data[rowMeans(data)>0,]

#The control sample was removed
group=gsub("(.*)\\_(.*)", "\\2", colnames(data))
data=data[,group=="Treat",drop=F]

#The samples were grouped according to the target gene expression, and the logFC of the high-low expression group was obtained
dataL=data[,data["DLD",]<median(data["DLD",]),drop=F]     #Data of low expression group
dataH=data[,data["DLD",]>=median(data["DLD",]),drop=F]    #Data from the high expression group
meanL=rowMeans(dataL)
meanH=rowMeans(dataH)
meanL[meanL<0.00001]=0.00001
meanH[meanH<0.00001]=0.00001
logFC=log2(meanH)-log2(meanL)
#Sequencing genes based on logfc
logFC=sort(logFC, decreasing=T)
genes=names(logFC)

#Read the gene set file
gmt=read.gmt("c2.cp.kegg.symbols.gmt"  )

#GSEA enrichment analysis
kk=GSEA(logFC, TERM2GENE=gmt, pvalueCutoff = 1)
kkTab=as.data.frame(kk)
kkTab=kkTab[kkTab$pvalue<0.05,]
#kkTab=kkTab[kkTab$p.adjust<0.05,]
write.table(kkTab,file="GSEA.result_DLD_KEGG.txt",sep="\t",quote=F,row.names = F)

#Plot the GSEA enrichment
termNum=6     #Display the number of paths
if(nrow(kkTab)>=termNum){
  showTerm=row.names(kkTab)[1:termNum]
  gseaplot=gseaplot2(kk, showTerm, base_size=8, title="DLD")
  pdf(file="GSEA_DLD_KEGG.pdf", width=7.5, height=5.5)
  print(gseaplot)
  dev.off()
}





##### 3.GSEA analysis of POLD1 using KEGG database
#Read the file and organize the input file
rt=read.table("merge.txt", header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data)
data=data[rowMeans(data)>0,]

#The control sample was removed
group=gsub("(.*)\\_(.*)", "\\2", colnames(data))
data=data[,group=="Treat",drop=F]

#The samples were grouped according to the target gene expression, and the logFC of the high-low expression group was obtained
dataL=data[,data["POLD1",]<median(data["POLD1",]),drop=F]     #Data of low expression group
dataH=data[,data["POLD1",]>=median(data["POLD1",]),drop=F]    #Data from the high expression group
meanL=rowMeans(dataL)
meanH=rowMeans(dataH)
meanL[meanL<0.00001]=0.00001
meanH[meanH<0.00001]=0.00001
logFC=log2(meanH)-log2(meanL)
#Sequencing genes based on logfc
logFC=sort(logFC, decreasing=T)
genes=names(logFC)

#Read the gene set file
gmt=read.gmt("c2.cp.kegg.symbols.gmt"  )

#GSEA enrichment analysis
kk=GSEA(logFC, TERM2GENE=gmt, pvalueCutoff = 1)
kkTab=as.data.frame(kk)
kkTab=kkTab[kkTab$pvalue<0.05,]
#kkTab=kkTab[kkTab$p.adjust<0.05,]
write.table(kkTab,file="GSEA.result_POLD1_KEGG.txt",sep="\t",quote=F,row.names = F)

#Plot the GSEA enrichment
termNum=6     #Display the number of paths
if(nrow(kkTab)>=termNum){
  showTerm=row.names(kkTab)[1:termNum]
  gseaplot=gseaplot2(kk, showTerm, base_size=8, title="POLD1")
  pdf(file="GSEA_POLD1_KEGG.pdf", width=7.5, height=5.5)
  print(gseaplot)
  dev.off()
}





