install.packages("reshape2")
install.packages("ggpubr")

library(reshape2)
library(ggpubr)


rt=read.table("CIBERSORT-Results.txt", header=T, sep="\t", check.names=F, row.names=1)

con=grepl("_Control", rownames(rt), ignore.case=T)
treat=grepl("_Treat", rownames(rt), ignore.case=T)
conData=rt[con,]
treatData=rt[treat,]
conNum=nrow(conData)
treatNum=nrow(treatData)
data=t(rbind(conData,treatData))

pdf(file="barplot.pdf", width=14.5, height=8.5)
col=rainbow(nrow(data), s=0.7, v=0.7)
par(las=1,mar=c(8,5,4,16),mgp=c(3,0.1,0),cex.axis=1.5)
a1=barplot(data,col=col,xaxt="n",yaxt="n",ylab="Relative Percent",cex.lab=1.8)
a2=axis(2,tick=F,labels=F)
axis(2,a2,paste0(a2*100,"%"))
par(srt=0,xpd=T)
rect(xleft = a1[1]-0.5, ybottom = -0.01, xright = a1[conNum]+0.5, ytop = -0.06,col="green")
text(a1[conNum]/2,-0.035,"Control",cex=1.8)
rect(xleft = a1[conNum]+0.5, ybottom = -0.01, xright =a1[length(a1)]+0.5, ytop = -0.06,col="red")
text((a1[length(a1)]+a1[conNum])/2,-0.035,"Treat",cex=1.8)
ytick2 = cumsum(data[,ncol(data)])
ytick1 = c(0,ytick2[-length(ytick2)])
legend(par('usr')[2]*0.98,par('usr')[4],legend=rownames(data),col=col,pch=15,bty="n",cex=1.2)
dev.off()


Type=gsub("(.*)\\_(.*)", "\\2", rownames(rt))
data=cbind(as.data.frame(t(data)), Type)
data=melt(data, id.vars=c("Type"))
colnames(data)=c("Type", "Immune", "Expression")

group=levels(factor(data$Type))
bioCol=c( "#377EB8" ,"#E41A1C","#4DAF4A" ,"#984EA3" ,"#FF7F00", "#FFFF33" ,"#A65628", "#F781BF", "#999999")
bioCol=bioCol[1:length(group)]
boxplot=ggboxplot(data, x="Immune", y="Expression", color="Type",
				  xlab="",
				  ylab="Fraction",
				  legend.title="Type",
				  add="none",
				  width=0.8,
				  palette=bioCol)+
				  rotate_x_text(50)+
	stat_compare_means(aes(group=Type),symnum.args=list(cutpoints=c(0, 0.001, 0.01, 0.05, 1), symbols=c("***", "**", "*", "")), label="p.signif")

pdf(file="immune.diff.pdf", width=8, height=6)
print(boxplot)
dev.off()



