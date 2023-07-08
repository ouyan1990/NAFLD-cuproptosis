if (!requireNamespace("BiocManager", quietly = TRUE))
   install.packages("BiocManager")
BiocManager::install("limma")

install.packages("ggpubr")

library(limma)
library(pheatmap)
library(reshape2)
library(ggpubr)

setwd("E:\\bioinformation\\cuproptosis\\NAFLD1\\18.testDiff")    

#Read the presentation data file
rt=read.table("normalize_GSE135251.txt", header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp), colnames(exp))
data=matrix(as.numeric(as.matrix(exp)), nrow=nrow(exp), dimnames=dimnames)
rt=avereps(data)


#If the data is not log2, the data is automatically log2
qx=as.numeric(quantile(rt, c(0, 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC=( (qx[5]>100) || ( (qx[6]-qx[1])>50 && qx[2]>0) )
if(LogC){
  rt[rt<0]=0
  rt=log2(rt+1)}
data=normalizeBetweenArrays(rt)



#Read the list of genes and extract the expression of the characteristics of the disease
geneRT=read.table("interGenes.txt", header=F, sep="\t", check.names=F)
data<- as.data.frame(data)
data=data[as.vector(geneRT[,1]),,drop=F]
exp=data

#Sample grouping information was extracted
Type=gsub("(.*)\\_(.*)", "\\2", colnames(data))
conNum=length(Type[Type=="Control"])
treatNum=length(Type[Type=="Treat"])


#Convert the presentation data to a ggplot2 input file
exp=as.data.frame(t(exp))
exp=cbind(exp, Type=Type)
data=melt(exp, id.vars=c("Type"))
colnames(data)=c("Type", "Gene", "Expression")

#Box plotting
p=ggboxplot(data, x="Gene", y="Expression", color = "Type", 
            xlab="",
            ylab="Gene expression",
            legend.title="Type",
            palette = c("#377EB8", "#E41A1C"),
            add="None",
            width=0.8)
p=p+rotate_x_text(60)
p1=p+stat_compare_means(aes(group=Type),
                        method="wilcox.test",
                        symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", " ")),
                        label = "p.signif")

#Output box diagram
pdf(file="boxplot_GSE135251.pdf", width=8, height=5)
print(p1)
dev.off()
