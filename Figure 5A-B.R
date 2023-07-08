install.packages("glmnet")
install.packages("pROC")


library(glmnet)
library(pROC)

setwd("E:\\bioinformation\\cuproptosis\\NAFLD1\\16.ROC")   

#Read input file
rt=read.table("diffGeneExp.txt", header=T, sep="\t", check.names=F, row.names=1)

#Obtain sample grouping information
y=gsub("(.*)\\_(.*)", "\\2", colnames(rt))
y=ifelse(y=="Control", 0, 1)

#Read the list file of genes
geneRT=read.table("interGenes.txt", header=F, sep="\t", check.names=F)

#Crossover genes were cycled and ROC curves were drawn
bioCol=rainbow(nrow(geneRT), s=0.9, v=0.9)   
aucText=c()
k=0
for(x in as.vector(geneRT[,1])){
	k=k+1
	#Roc curve
	roc1=roc(y, as.numeric(rt[x,]))     #The parameters of the ROC curve are obtained
	if(k==1){
		pdf(file="ROC.genes.pdf", width=5, height=4.75)
		plot(roc1, print.auc=F, col=bioCol[k], legacy.axes=T, main="")
		aucText=c(aucText, paste0(x,", AUC=",sprintf("%.3f",roc1$auc[1])))
	}else{
		plot(roc1, print.auc=F, col=bioCol[k], legacy.axes=T, main="", add=TRUE)
		aucText=c(aucText, paste0(x,", AUC=",sprintf("%.3f",roc1$auc[1])))
	}
}
#Draw a legend to get the area under the ROC curve
legend("bottomright", aucText, lwd=2, bty="n", col=bioCol[1:(ncol(rt)-1)])
dev.off()

#Build a logical model
rt=rt[as.vector(geneRT[,1]),]
rt=as.data.frame(t(rt))
logit=glm(y ~ ., family=binomial(link='logit'), data=rt)
pred=predict(logit, newx=rt)     #Get a score on the model

#Draw the ROC curve of the model
roc1=roc(y, as.numeric(pred))      #The parameters of the model ROC curve are obtained
ci1=ci.auc(roc1, method="bootstrap")     #The fluctuation range of the area under the ROC curve is obtained
ciVec=as.numeric(ci1)
pdf(file="ROC.model.pdf", width=5, height=4.75)
plot(roc1, print.auc=TRUE, col="red", legacy.axes=T, main="Model")
text(0.39, 0.43, paste0("95% CI: ",sprintf("%.03f",ciVec[1]),"-",sprintf("%.03f",ciVec[3])), col="red")
dev.off()



