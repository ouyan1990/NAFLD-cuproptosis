
install.packages("e1071")

set.seed(12345)
library(e1071)

setwd("E:\\bioinformation\\cuproptosis\\NAFLD1\\14.SVM")     
source("geoCRG12.msvmRFE.R")

#Read input file
data=read.table("diffGeneExp.txt", header=T, sep="\t", check.names=F, row.names=1)
data=as.data.frame(t(data))
#Obtain sample grouping information
group=gsub("(.*)\\_(.*)", "\\2", row.names(data))
data=cbind(group, data)
data$group=factor(data$group, levels=c("Control","Treat"))

#Construct Machine Learning-Support Vector Machine Recursive Feature Elimination Algorithm (SVM-RFE)
svmRFE(data, k=10, halve.above=50)
nfold=10
geneNum=nrow(data)
folds=rep(1:nfold, len=geneNum)[sample(geneNum)]
folds=lapply(1:nfold, function(x) which(folds == x))
results=lapply(folds, svmRFE.wrap, data, k=10, halve.above=50)

#Rank the importance of the characteristic genes
top.features=WriteFeatures(results, data, save=F)
#Output the sorted result
write.table(top.features, file="feature_svm.txt", sep="\t", quote=F,row.names=F)

#Cross validation
featsweep=lapply(1:9, FeatSweep.wrap, results, data)

#Get the error of cross validation
no.info=min(prop.table(table(data[,1])))
errors=sapply(featsweep, function(x) ifelse(is.null(x), NA, x$error))

#Plot cross-validation errors
pdf(file="errors.pdf", width=5, height=5)
PlotErrors(errors, no.info=no.info)
dev.off()

#Draw graphs that cross-verify accuracy
pdf(file="accuracy.pdf", width=5, height=5)
Plotaccuracy(1-errors, no.info=no.info)
dev.off()

#The characteristic genes of SVM were output
featureGenes=top.features[1:which.min(errors),1,drop=F]
write.table(file="SVM-RFE.gene.txt", featureGenes, sep="\t", quote=F, row.names=F, col.names=F)

