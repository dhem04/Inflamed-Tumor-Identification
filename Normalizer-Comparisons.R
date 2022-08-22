library(limma)
library(edgeR)
library(Biobase)
list.files("", full.names=TRUE, pattern="*.txt")
files0<-list.files("Desktop Directory", full.names=TRUE, pattern="*.txt")
matrix_HTBE<-read.table(files0[1], sep="\t", header=TRUE)
matrix_B<-matrix_HTBE
matrix0<-read.table("Desktop Directory", sep="\t", header=TRUE, row.names=1)
matrix0<-matrix0[,-1]
matrixA<-matrix0[,3]
matrixA<-cbind(matrixA, matrix0[5])
colnames(matrixA)<-c("Genes","Cell Line")
samples0<-list.files("Desktop Directory", full.names=FALSE, pattern="*.txt")
samples0<-gsub("\\.txt", "", samples0)
colnames(matrix_HTBE)<-samples0[1]
for (i in 2:6){
  n1<-read.table(files0[i], sep="\t", header=TRUE)
  colnames(n1)<-samples0[i]
  matrix_HTBE<-cbind(matrix_HTBE, n1)
}
data.frame(ID=samples0, Treat=c(rep("Control", 3), rep("Treatment Type", 3)))
design1<-data.frame(ID=samples1, Treat=c(rep("Control", 3), rep("Treatment Type", 3)))
pd1<-new("AnnotatedDataFrame", data=design1)
model.matrix(~Treat, data=design1)
model1<-model.matrix(~Treat, data=design1)
dge1<-DGEList(counts=matrix0)
dim(matrix1)
keep1<-filterByExpr(dge1, model1)
dge2<-dge1[keep1,,keep.lib.sizes=FALSE]
matrix2<-getCounts(dge2)
v1<-voom(matrix2, model1, plot=FALSE)
str(design1)
fit1<-lmFit(v1, model1)
efit1<-eBayes(fit1)
res1<-topTable(efit1, coef=2, n=50000)
pl1<-read.table("Desktop Directory", sep=",", header=TRUE)
res1$Gene<-pl1$Gene.name[match(rownames(res1), pl1$Gene.stable.ID)]
res1$Gene<-NULL
res1<-cbind(Gene=pl1$Gene.name[match(rownames(res1), pl1$Gene.stable.ID)], res1)
write.table(res1, file="File Name", sep="\t", row.names=FALSE)
