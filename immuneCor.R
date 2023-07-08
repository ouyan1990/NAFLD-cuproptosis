if (!requireNamespace("BiocManager", quietly = TRUE))
   install.packages("BiocManager")
BiocManager::install("limma")

install.packages("tidyverse")
install.packages("ggplot2")
install.packages("ggExtra")


library(limma)
library(reshape2)
library(tidyverse)
library(ggplot2)

setwd("E:\\bioinformation\\cuproptosis\\NAFLD1\\25. ssGSEA\\25.3 ssGSEA immuneCor")     

#Read the immune cell score file
immune=read.table("ssGSEA.result.txt", header=T, sep="\t", check.names=F, row.names=1)
immune=t(immune)

#Read the presentation data file
rt=read.table("merge.txt", header=T, sep="\t", check.names=F, row.names=1)

#Read the gene list file and extract the expression amount of disease characteristic genes
geneRT=read.table("interGenes.txt", header=F, sep="\t", check.names=F)
exp=t(rt[as.vector(geneRT[,1]),])

#Data intersection
sameSample=intersect(row.names(immune), row.names(exp))
treat=grepl("_treat", sameSample, ignore.case=T)
sameSample=sameSample[treat]
immune=immune[sameSample,,drop=F]
exp=exp[sameSample,,drop=F]

#Correlation analysis
outTab=data.frame()
for(cell in colnames(immune)){
	for(gene in colnames(exp)){
		x=as.numeric(immune[,cell])
		y=as.numeric(exp[,gene])
		corT=cor.test(x,y,method="spearman")
		cor=corT$estimate
		pvalue=corT$p.value
		text=ifelse(pvalue<0.001,"***",ifelse(pvalue<0.01,"**",ifelse(pvalue<0.05,"*","")))
		outTab=rbind(outTab,cbind(Gene=gene, Immune=cell, cor, text, pvalue))
	}
}

#Draw a correlation heat map
outTab$cor=as.numeric(outTab$cor)
pdf(file="cor.pdf", width=7, height=6)
ggplot(outTab, aes(Gene, Immune)) + 
	geom_tile(aes(fill = cor), colour = "grey", size = 1)+
	scale_fill_gradient2(low = "#02C9AF", mid = "white", high = "#E41A1C") + 
	geom_text(aes(label=text),col ="black",size = 3) +
	theme_minimal() +    #È¥µô±³¾°
	theme(axis.title.x=element_blank(), axis.ticks.x=element_blank(), axis.title.y=element_blank(),
	      axis.text.x = element_text(angle = 45, hjust = 1, size = 10, face = "bold"),  
	      axis.text.y = element_text(size = 10, face = "bold")) +       
	labs(fill =paste0("***  p<0.001","\n", "**  p<0.01","\n", " *  p<0.05","\n", "\n","Correlation")) +   
	scale_x_discrete(position = "bottom")     
dev.off()

write.table(outTab, file = "CorResults.txt",quote = F,sep = "\t")


