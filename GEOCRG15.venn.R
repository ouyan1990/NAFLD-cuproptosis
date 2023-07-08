install.packages("VennDiagram")


library(VennDiagram)        
setwd("E:\\bioinformation\\cuproptosis\\NAFLD1\\15.venn")    
files=dir()     
files=grep("txt",files,value=T)     
geneList=list()

#Read the genetic information in all tthe files and be saved to genelist
for(i in 1:length(files)){
  inputFile=files[i]
  if(inputFile=="intersect.txt"){next}
  rt=read.table(inputFile,header=F)
  header=unlist(strsplit(inputFile,"\\.|\\-"))
  geneList[[header[1]]]=as.vector(rt[,1])
  uniqLength=length(unique(as.vector(rt[,1])))
  print(paste(header[1],uniqLength,sep=" "))
}

#Venn diagram
venn.plot=venn.diagram(geneList,filename=NULL,fill=c("#a50026", "#e0f3f8",  "#f7f7f7"),
                      cat.cex = 1.2)
pdf(file="venn.pdf",width=6,height=6)
grid.draw(venn.plot)
dev.off()

#Intersection gene preservation
intersectGenes=Reduce(intersect,geneList)
write.table(file="intersect.txt",intersectGenes,sep="\t",quote=F,col.names=F,row.names=F)


