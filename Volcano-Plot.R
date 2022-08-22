library(ggrepel)
save_volcano_plot<-function(res1, path1, title1){
  library(ggplot2)
  gg1<-ggplot(data=res1, aes(x=logFC, y=-log10(P.Value), col=diffexpressed, label=Glabel)) +geom_hline(yintercept=-log10(0.05), col="black")+geom_text_repel()+geom_point(alpha=0.7, size=0.75) +xlim(c(-9, 9)) + ylim(c(0, 50)) +xlab("log2 fold change") + ylab("-log10 p-value")+scale_color_manual(values=c("blue", "black", "red")) +ggtitle(title1)
  ggsave(gg1, filename=path1)
}
res11<-read.table("File Directory", sep="\t", row.names=NULL)
matrix15<-res11[1]
res11$Gene<- "N"
res11[8]<-matrix15
for (i in 1:14400){
  w<-""
  u<-matrix15[i,]
  u_split<-strsplit(u, "")[[1]]
  for(number in u_split){
    if(number=="."){
      break
    }
    else{
      w<-paste(w,number, sep="")
    }
  }
  matrix15[i,]<-w
}
for (i in 1:50000){
  l<-FALSE
  w<-""
  u<-matrix15[i,]
  u_split<-strsplit(u, "")[[1]]
  for(number in u_split){
    if(number=="." && l)
      break
    if(l){
      w<-paste(w,number, sep="")
    }
    if(number=="|"){
      l<-TRUE
    }
  }
  matrix15[i,]<-w
}
res11$Gene<- NA
res11[8]<-matrix14
pl1<-read.table("File Directory", sep=",", header=TRUE)
matrix15<-pl1$Gene.name[match(matrix16[,], pl1$Gene.stable.ID)]
res11<-read.table("File Directory", sep="\t", header=TRUE, row.names=NULL)
res11$diffexpressed <- "Not Significant"
res11$diffexpressed[res11$logFC > 0 & res11$P.Value < 0.05] <- "UP"
res11$diffexpressed[res11$logFC < 0 & res11$P.Value < 0.05] <- "DOWN"
res11$check<- "Wrong"
res11$check[res11$logFC > 4 & res11$P.Value <5.39e-9] <- "UP"
res11$check[res11$logFC < -3 & res11$P.Value <3.06e-10] <- "DOWN"
res11$Glabel <- NA
res11$Glabel[res11$check != "Wrong"] <- res11$Gene[res11$check != "Wrong"]
save_volcano_plot(res11, "File Name", "Myeloid_Cells TGF-b Treated 12 Hours")
save_volcano_plot<-function(res1, path1, title1){
  library(ggplot2)
  gg1<-ggplot(data=res1, aes(x=logFC, y=-log10(P.Value), colour=FDR_0.1_pval_0.05_FC_2)) +geom_point(alpha=0.7, size=0.75) +xlim(c(-9, 9)) + ylim(c(0, 50)) +xlab("log2 fold change") + ylab("-log10 p-value")+ggtitle(paste0("RNAseq ",title1))
  ggsave(gg1, filename=path1)}
