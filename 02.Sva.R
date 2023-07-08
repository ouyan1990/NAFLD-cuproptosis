
#### Merging the enrolled four GSE datasets
library(limma)
library(sva)
outFile="merge.txt"       

files=dir()
files=grep("txt$", files, value=T)
geneList=list()

#Read genetic information in all txt files and save to geneList
for(file in files){
  if(file==outFile){next}
  rt=read.table(file, header=T, sep="\t", check.names=F)      #Read input file
  geneNames=as.vector(rt[,1])      #Extraction of gene name
  uniqGene=unique(geneNames)      
  header=unlist(strsplit(file, "\\.|\\-"))
  geneList[[header[1]]]=uniqGene
}

#Acquisition of intersection gene
interGenes=Reduce(intersect, geneList)

#Data Merging 
allTab=data.frame()
batchType=c()
for(i in 1:length(files)){
  inputFile=files[i]
  header=unlist(strsplit(inputFile, "\\.|\\-"))
  rt=read.table(inputFile, header=T, sep="\t", check.names=F)
  rt=as.matrix(rt)
  rownames(rt)=rt[,1]
  exp=rt[,2:ncol(rt)]
  dimnames=list(rownames(exp),colnames(exp))
  data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
  rt=avereps(data)
  rt=normalizeBetweenArrays(rt)
  if(i==1){
    allTab=rt[interGenes,]
  }else{
    allTab=cbind(allTab, rt[interGenes,])
  }
  batchType=c(batchType, rep(i,ncol(rt)))
}

boxplot(allTab)
#The data is corrected and the corrected results are output
outTab=ComBat(allTab, batchType, par.prior=TRUE)
boxplot(outTab)
outTab=rbind(geneNames=colnames(outTab), outTab)

write.table(outTab, file="normalize.txt", sep="\t", quote=F, col.names=F)



