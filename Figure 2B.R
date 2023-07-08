
install.packages("VennDiagram")

library(VennDiagram)      
geneList=list()

#read file
rt=read.table("diffGeneExp.txt" , header=F, sep="\t", check.names=F)
geneNames=as.vector(rt[,1])              
geneNames=geneNames[-1]
geneNames=gsub("^ | $","",geneNames)    
uniqGene=unique(geneNames)               
geneList[["DEGs"]]=uniqGene             

rt=read.table("gene.txt" , header=F, sep="\t", check.names=F)
geneNames=as.vector(rt[,1])            
geneNames=gsub("^ | $","",geneNames)     
uniqGene=unique(geneNames)              
geneList[["cuproptosis-related gene"]]=uniqGene          

#Draw venn diagram 
venn.plot=venn.diagram(geneList,filename=NULL,fill=c("#E41A1C", "#377EB8"),scaled=FALSE,cat.pos=c(-1,1),cat.col = c("cornflowerblue", "darkorchid1"),cat.cex = 1.2)
pdf(file="venn.pdf", width=5, height=5)
grid.draw(venn.plot)
dev.off()

#Intersection genes
interGenes=Reduce(intersect,geneList)
write.table(file="interGenes.txt", interGenes, sep="\t", quote=F, col.names=F, row.names=F)



