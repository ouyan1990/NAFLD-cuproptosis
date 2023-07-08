

install.packages("RCircos")

library("RCircos")      
setwd("E:\\bioinformation\\cuproptosis\\NAFLD1\\07 circle\\07.2 Rcircos")    

#Initializes the circle graph
cytoBandIdeogram=read.table("refer.txt", header=T, sep="\t", check.names=F)
chr.exclude <- NULL
cyto.info <- cytoBandIdeogram
tracks.inside <- 4
tracks.outside <- 0
RCircos.Set.Core.Components(cyto.info, chr.exclude, tracks.inside, tracks.outside)

#Set the parameters of the circle graph
rcircos.params <- RCircos.Get.Plot.Parameters()
rcircos.params$text.size=0.7
rcircos.params$point.size=5
RCircos.Reset.Plot.Parameters(rcircos.params)

#File output
pdf(file="RCircos.pdf", width=7, height=7)
RCircos.Set.Plot.Area()
RCircos.Chromosome.Ideogram.Plot()

#Read the gene annotation file and label the name of the gene

RCircos.Gene.Label.Data=read.table("Rcircos.geneLabel.txt", header=T, sep="\t", check.names=F)
name.col <- 4
side <- "in"
track.num <- 1
RCircos.Gene.Connector.Plot(RCircos.Gene.Label.Data, track.num, side)
track.num <- 2
RCircos.Gene.Name.Plot(RCircos.Gene.Label.Data, name.col, track.num, side)
dev.off()



