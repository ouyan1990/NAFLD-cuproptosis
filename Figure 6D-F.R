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

#Removing the control sample 
group=gsub("(.*)\\_(.*)", "\\2", colnames(data))
data=data[,group=="Treat",drop=F]

#Read the gene set file
geneSets=getGmt("c2.cp.kegg.symbols.gmt", geneIdType=SymbolIdentifier())

#GSVA analysis
ssgseaScore=gsva(data, geneSets, method='ssgsea', kcdf='Gaussian', abs.ranking=TRUE)
#The scoring was corrected
normalize=function(x){
  return((x-min(x))/(max(x)-min(x)))}
ssgseaScore=normalize(ssgseaScore)
ssgseaScore=ssgseaScore[order(apply(ssgseaScore,1,sd),decreasing=T),]
ssgseaScore=ssgseaScore[1:50,]

#The samples were grouped according to the amount of target gene expression
lowName=colnames(data)[data["DLD",]<median(data["DLD",])]       #Samples from the low expression group
highName=colnames(data)[data["DLD",]>=median(data["DLD",])]     #Samples from the high expression group
lowScore=ssgseaScore[,lowName]
highScore=ssgseaScore[,highName]
data=cbind(lowScore, highScore)
conNum=ncol(lowScore)
treatNum=ncol(highScore)
Type=c(rep("Control",conNum), rep("Treat",treatNum))

#Gene difference analysis
outTab=data.frame()
for(i in row.names(data)){
	test=t.test(data[i,] ~ Type)
	pvalue=test$p.value
	t=test$statistic
	Sig=ifelse(pvalue>0.05, "Not", ifelse(t>0,"Up","Down"))
	outTab=rbind(outTab, cbind(Pathway=i, t, pvalue, Sig))
}

#Draw a bar chart
pdf(file="barplot_DLD_KEGG.pdf", width=12, height=9)
outTab$t=as.numeric(outTab$t)
outTab$Sig=factor(outTab$Sig, levels=c("Down", "Not", "Up"))
gg1=ggbarplot(outTab, x="Pathway", y="t", fill = "Sig", color = "white",
		palette=c("#377EB8","grey","#E41A1C"), sort.val = "asc", sort.by.groups = T,
		rotate=TRUE, legend="right", title="DLD",
		xlab="Term", ylab="t value of GSVA score",  legend.title="Group", x.text.angle=60)
print(gg1)
dev.off()



