library(ggpubr)
library(reshape2)
setwd("E:\\bioinformation\\cuproptosis\\NAFLD1\\25. ssGSEA\\25.2 ssGSEA_vioplot")

data11=read.table("ssGSEA.result.txt",header = T,sep = "\t",check.names =F)
data11=t(data11)
colnames(data11)=data11[1,]
data11=data11[-1,]

risk=gsub("(.*)\\_(.*)", "\\2", rownames(data11))

#Convert the presentation data to a ggplot2 input file
exp=cbind(risk,data11)
exp<-exp[ ,c(1,2)]
rownames(exp)=gsub("_Treat", "", rownames(exp))   
rownames(exp)=gsub("_Control", "", rownames(exp))

exp<- cbind(rownames(exp),exp)
exp<-exp[,-3]
colnames(exp)<-c("id","mygroup")


write.table(exp,"clinical.txt",row.names = F,col.names = T,sep = "\t",quote = F)


rownames(data11)=gsub("_Treat", "", rownames(data11))   
rownames(data11)=gsub("_Control", "", rownames(data11))
data11<-cbind(rownames(data11),data11)
colnames(data11)[1]<-"id"

write.table(data11,"data1.txt",row.names = F,col.names = T,sep = "\t",quote = F)

data1=read.table("data1.txt",header = T,sep = "\t",check.names =F)
data2=read.table("clinical.txt",header = T,sep = "\t",check.names = F)
data3=merge(data2,data1,by="id")

### Box plot of immune cells
imc=c("aDCs","B_cells","CD8+_T_cells","DCs","iDCs","Macrophages",
      "Mast_cells","Neutrophils","NK_cells","pDCs","T_helper_cells",
      "Tfh","Th1_cells","Th2_cells","TIL","Treg")
imcd=data3[,c(imc,"mygroup")]
newdata=melt(imcd,id.vars=c("mygroup"))
colnames(newdata)=c("Group","Type","Score")
a=ggboxplot(newdata, x="Type", y="Score", color = "Group",
            ylab="Score",add = "none",xlab="",palette = c("#377EB8","#E41A1C") )
a=a+rotate_x_text(51)
pdf("imc.pdf",8,8)            
a+stat_compare_means(aes(group=newdata$Group),symnum.args=list(cutpoints = c(0,0.001,0.01,0.05, 1), symbols = c("***", "**","*", " ")),label = "p.signif")
dev.off()




######================#Box plot of immune-related functions
imf=c("APC_co_inhibition","APC_co_stimulation","CCR",
      "Check-point","Cytolytic_activity","HLA","Inflammation-promoting",
      "MHC_class_I","Parainflammation","T_cell_co-inhibition",
      "T_cell_co-stimulation","Type_I_IFN_Reponse","Type_II_IFN_Reponse")
imfd=data3[,c(imf,"mygroup")]
newdata=melt(imfd,id.vars=c("mygroup"))
colnames(newdata)=c("Group","Type","Score")
a=ggboxplot(newdata, x="Type", y="Score", color = "Group",
            ylab="Score",add = "none",xlab="",palette = c("#377EB8","#E41A1C") )
a=a+rotate_x_text(51)
pdf("imf.pdf",8,8)            
a+stat_compare_means(aes(group=newdata$Group),symnum.args=list(cutpoints = c(0,0.001,0.01,0.05, 1), symbols = c("***", "**","*", " ")),label = "p.signif")
dev.off()


