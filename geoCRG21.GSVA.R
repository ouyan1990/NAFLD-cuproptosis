

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

setwd("E:\\bioinformation\\cuproptosis\\NAFLD1\\21.GSVA")    

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


























#=========================================
######Video source: https://ke.biowolf.cn

expFile="merge.txt"              #表达数据文件
gmtFile="c5.go.symbols.gmt"     #基因集文件
setwd("E:\\bioinformation\\cuproptosis\\NAFLD1\\19.GSVA")     #设置工作目录

#读取表达输入文件,并对输入文件整理
rt=read.table("merge.txt", header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data)

#去除对照组的样品
group=gsub("(.*)\\_(.*)", "\\2", colnames(data))
data=data[,group=="Treat",drop=F]

#读取基因集文件
geneSets=getGmt("c2.cp.kegg.symbols.gmt", geneIdType=SymbolIdentifier())

#GSVA分析
ssgseaScore=gsva(data, geneSets, method='ssgsea', kcdf='Gaussian', abs.ranking=TRUE)
#对打分进行矫正
normalize=function(x){
  return((x-min(x))/(max(x)-min(x)))}
ssgseaScore=normalize(ssgseaScore)
ssgseaScore=ssgseaScore[order(apply(ssgseaScore,1,sd),decreasing=T),]
ssgseaScore=ssgseaScore[1:30,]

#根据目标基因的表达量对样品进行分组
lowName=colnames(data)[data[gene,]<median(data[gene,])]       #低表达组的样品
highName=colnames(data)[data[gene,]>=median(data[gene,])]     #高表达组的样品
lowScore=ssgseaScore[,lowName]
highScore=ssgseaScore[,highName]
data=cbind(lowScore, highScore)
conNum=ncol(lowScore)
treatNum=ncol(highScore)
Type=c(rep("Control",conNum), rep("Treat",treatNum))

#基因差异分析
outTab=data.frame()
for(i in row.names(data)){
  test=t.test(data[i,] ~ Type)
  pvalue=test$p.value
  t=test$statistic
  Sig=ifelse(pvalue>0.05, "Not", ifelse(t>0,"Up","Down"))
  outTab=rbind(outTab, cbind(Pathway=i, t, pvalue, Sig))
}

#绘制柱状图
pdf(file="barplot_GO_DLD.pdf", width=12, height=9)
outTab$t=as.numeric(outTab$t)
outTab$Sig=factor(outTab$Sig, levels=c("Down", "Not", "Up"))
gg1=ggbarplot(outTab, x="Pathway", y="t", fill = "Sig", color = "white",
              palette=c("#377EB8","grey","#E41A1C"), sort.val = "asc", sort.by.groups = T,
              rotate=TRUE, legend="right", title=gene,
              xlab="Term", ylab="t value of GSVA score",  legend.title="Group", x.text.angle=60)
print(gg1)
dev.off()
