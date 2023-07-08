

install.packages("corrplot")
install.packages("circlize")


library(corrplot)
library(circlize)

setwd("E:\\bioinformation\\cuproptosis\\NAFLD1\\09.corrplot")     

#Read input file
data=read.table("diffGeneExp.txt"   , header=T, sep="\t", check.names=F, row.names=1)

#removing control sample
group=gsub("(.*)\\_(.*)", "\\2", colnames(data))
data=data[,group=="Treat",drop=F]
rt=t(data)

#correlation analysis
cor1=cor(rt)

#Set graphic color
col = c(rgb(1,0,0,seq(1,0,length=32)),rgb(0,1,0,seq(0,1,length=32)))
cor1[cor1==1]=0
c1 = ifelse(c(cor1)>=0,rgb(1,0,0,abs(cor1)),rgb(0,1,0,abs(cor1)))
col1 = matrix(c1,nc=ncol(rt))

#Drawing diagram
pdf(file="circos.pdf", width=7, height=7)
par(mar=c(2,2,2,4))
circos.par(gap.degree=c(3,rep(2, nrow(cor1)-1)), start.degree = 180)
chordDiagram(cor1, grid.col=rainbow(ncol(rt)),  transparency = 0.5, symmetric = T)
par(xpd=T)
colorlegend(col, vertical = T,labels=c(1,0,-1),xlim=c(1.1,1.3),ylim=c(-0.4,0.4))
dev.off()
circos.clear()




