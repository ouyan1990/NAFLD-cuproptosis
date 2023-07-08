install.packages("glmnet")


set.seed(12345)
library(glmnet)     

setwd("E:\\bioinformation\\cuproptosis\\NAFLD1\\12. lasso")      

#Read input file
rt=read.table("diffGeneExp.txt" , header=T, sep="\t", check.names=F, row.names=1)
rt=t(rt)

#Model construct
x=as.matrix(rt)
y=gsub("(.*)\\_(.*)", "\\2", row.names(rt))
fit=glmnet(x, y, family = "binomial", alpha=1)
cvfit=cv.glmnet(x, y, family="binomial", alpha=1,type.measure='deviance',nfolds = 10)

#Draw a graph of Lasso regression
pdf(file="lasso.pdf", width=6, height=5.5)
plot(fit)
dev.off()
#Draw cross-validated graphs
pdf(file="cvfit.pdf",width=6,height=5.5)
plot(cvfit)
dev.off()

#Output the screened feature genes
coef=coef(fit, s=cvfit$lambda.min)
index=which(coef != 0)
lassoGene=row.names(coef)[index]
lassoGene=lassoGene[-1]
write.table(lassoGene, file="LASSO.gene.txt", sep="\t", quote=F, row.names=F, col.names=F)


