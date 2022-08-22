library(biomaRt)
library(cluster)
library(ggplot2)
?read.table
library(org.Hs.eg.db)
hs <- org.Hs.eg.db
my.symbols <- c(down1)
g<-AnnotationDbi::select(hs, keys = my.symbols, columns = c("ENTREZID", "SYMBOL"),keytype = "SYMBOL")
DOWNN<-g[,2]
biomart_converter<-ensids_to_entrez<-function(genelists){
  library(biomaRt)
  ensembl <- useEnsembl(biomart = "ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl")
  #ensembl <- useEnsembl(biomart = "ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl", host = "http://www.ensembl.org", mirror = "uswest")
  res_list<-lapply(genelists, function(x) getBM(attributes = c('ensembl_gene_id', 'entrezgene_id'),filters = 'ensembl_gene_id', values = x, mart = ensembl))
  res_list<-lapply(res_list, function(x) x[complete.cases(x),])
  res_list<-lapply(res_list, function(x) x[match(unique(x[,1]), x[,1]),2])
  res_list}
list.files("Desktop Directory", full.names=TRUE)
for (i in 1:4686){
  w<-""
  u<-up1[i]
  u_split<-strsplit(u, "")[[1]]
  for(number in u_split){
    if(number=="."){
      break
    }
    else{
      w<-paste(w,number, sep="")
    }
  }
  up1[i]<-w
}
w<-0
for(i in 1:611){
  if(is.na(UPP[i])){
    UPP<-UPP[-i]
    i<-(i-1)
  }
}
for(i in 1:554){
  if(is.na(DOWNN[i])){
    DOWNN<-DOWNN[-i]
    i<-i-1
  }
}
up1[-210]
matrix15<-res1[,1]
matrix3<-read.table(up1, sep="|", header=TRUE)
pl1<-read.table("Desktop Directory", sep=",", header=TRUE)

UP123<-pl1$Gene.name[match(up1, pl1$Gene.stable.ID)]
DOWN123<-pl1$Gene.name[match(down1, pl1$Gene.stable.ID)]

res1<-read.table("Desktop Directory", sep="\t", header=TRUE, row.names=1)
entids_micro_up<-list("Desktop Directory"=entid_micro1[[1]], "Desktop Directory"=entid_micro2[[1]], "Desktop Directory"=entid_micro3[[1]],"Desktop Directory"=entid_micro4[[1]], "Desktop Directory"=entid_micro5[[1]], "Desktop Directory"=entid_micro6[[1]], "Desktop Directory"=entid_micro7[[1]], "Desktop Directory"=entid_micro8[[1]], "Desktop Directory"=entid_micro9[[1]], "Desktop Directory"=entid_micro10[[1]])
up1<-rownames(res1)[res1$logFC>1&res1$P.Value<0.05&res1$adj.P.Val<0.1]
up11<-res1[,1][res1$logFC>1]
down1<-rownames(res1)[res1$logFC<(-1)&res1$P.Value<0.05&res1$adj.P.Val<0.1]
write.table(up1, file="File Name", sep="\t")
down1<-rownames(res1[,1])[res1$logFC<(-1)]
library(org.Hs.eg.db)
hs <- org.Hs.eg.db
my.symbols <- c("ANKRD62P1-PARP4P3")
select(hs, keys = my.symbols, columns = c("ENTREZID", "SYMBOL"),keytype = "SYMBOL")
entidnew<-biomart_human_hngc_to_entrez(list(File Name=up1, File Name=down1))
alias2Symbol(up1, species="Hs")
entid22<-biomart_human_hngc_to_entrez(list(File Name=up1, File Name=down1))
entid23<-biomart_human_ENST_to_ENSG(list(File Name=up1, File Name=down1))

entid35<-(list(File Name=UPP, File Name=DOWNN))
entid1<-list(GSE96885_Myeloid_Cells_TGF_b_Treated_36_Hours_UP=UPP)
entid_up1<-list(GSE156295=entid1[[1]], GSE156295=entid2[[1]], GSE157052=entid3[[1]], GSE161733=entid4[[1]], GSE110021=entid5[[1]], GSE148767=entid6[[1]], GSE148767=entid7[[1]], GSE148767=entid8[[1]], GSE148767=entid9[[1]], GSE173953=entid10[[1]], GSE173953=entid11[[1]], GSE152409=entid12[[1]], GSE152409=entid13[[1]], GSE152409=entid14[[1]], GSE152409=entid23[[1]], GSE109182=entid16[[1]], GSE132704=entid17[[1]], GSE132704=entid18[[1]], GSE132704=entid19[[1]], GSE163067=entid20[[1]], GSE161732=entid21[[1]], GSE178640=entid22[[1]])
entid_up<-list("Desktop Directory"=entid1[[1]], "Desktop Directory"=entid2[[1]], "Desktop Directory"=entid3[[1]], "Desktop Directory"=entid4[[1]], "Desktop Directory"=entid5[[1]], "Desktop Directory"=entid6[[1]], "Desktop Directory"=entid7[[1]], "Desktop Directory"=entid8[[1]], "Desktop Directory"=entid9[[1]], "Desktop Directory"=entid10[[1]], "Desktop Directory"=entid11[[1]], "Desktop Directory"=entid12[[1]], "Desktop Directory"=entid13[[1]], "Desktop Directory"=entid14[[1]], "Desktop Directory"=entid23[[1]], "Desktop Directory"=entid16[[1]], "Desktop Directory"=entid17[[1]], "Desktop Directory"=entid18[[1]], "Desktop Directory"=entid19[[1]], "Desktop Directory"=entid20[[1]], "Desktop Directory"=entid21[[1]], "Desktop Directory"=entid22[[1]])
entid_down<-list("Desktop Directory"=entid1[[2]], "Desktop Directory"=entid2[[2]], "Desktop Directory"=entid3[[2]], "Desktop Directory"=entid4[[2]], "Desktop Directory"=entid5[[2]], "Desktop Directory"=entid6[[2]], "Desktop Directory"=entid7[[2]], "Desktop Directory"=entid8[[2]], "Desktop Directory"=entid9[[2]], "Desktop Directory"=entid10[[2]], "Desktop Directory"=entid11[[2]], "Desktop Directory"=entid12[[2]], "Desktop Directory"=entid13[[2]], "Desktop Directory"=entid14[[2]], "Desktop Directory"=entid23[[2]], "Desktop Directory"=entid16[[2]], "Desktop Directory"=entid17[[2]], "Desktop Directory"=entid18[[2]], "Desktop Directory"=entid19[[2]], "Desktop Directory"=entid20[[2]], "Desktop Directory"=entid21[[2]], "Desktop Directory"=entid22[[2]])
GSE110021_D1_TGFb_Treated=entid15[[1]]
ensembl<- useEnsembl(biomart = "ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl", host = http://www.ensembl.org, mirror = "uswest")
Genes<-biomart_human_hngc_to_entrez(entid35)
biomart_human_hngc_to_entrez<-function(genelists){
  library(biomaRt)
  ensembl <- useEnsembl(biomart = "ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl")
  res_list<-lapply(genelists, function(x) getBM(attributes = c('entrezgene_id', 'hgnc_symbol'),filters = 'entrezgene_id', values = x, mart = ensembl))
  res_list<-lapply(res_list, function(x) x[complete.cases(x),])
  res_list<-lapply(res_list, function(x) x[match(unique(x[,1]), x[,1]),2])
  res_list}
biomart_human_ENST_to_ENSG<-function(genelists){
  library(biomaRt)
  ensembl <- useEnsembl(biomart = "ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl")
  res_list<-lapply(genelists, function(x) getBM(attributes = c('ensembl_transcript_id', 'entrezgene_id'),filters = 'ensembl_transcript_id', values = x, mart = ensembl))
  res_list<-lapply(res_list, function(x) x[complete.cases(x),])
  res_list<-lapply(res_list, function(x) x[match(unique(x[,1]), x[,1]),2])
  res_list}
biomart_human_hngc_to_entrez<-function(genelists){
  library(biomaRt)
  ensembl <- useEnsembl(biomart = "ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl")
  res_list<-lapply(genelists, function(x) getBM(attributes = c('entrezgene_id', 'hgnc_symbol'),filters = 'entrezgene_id', values = x, mart = ensembl))
  res_list<-lapply(res_list, function(x) x[complete.cases(x),])
  res_list<-lapply(res_list, function(x) x[match(unique(x[,1]), x[,1]),2])
  res_list}
biomart_human_ENST_to_ENSG<-function(genelists){
  library(biomaRt)
  ensembl <- useEnsembl(biomart = "ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl")
  res_list<-lapply(genelists, function(x) getBM(attributes = c('ensembl_transcript_id', 'entrezgene_id'),filters = 'ensembl_transcript_id', values = x, mart = ensembl))
  res_list<-lapply(res_list, function(x) x[complete.cases(x),])
  res_list<-lapply(res_list, function(x) x[match(unique(x[,1]), x[,1]),2])
  res_list}
library(org.Hs.eg.db)
library(clusterProfiler)
library(enrichplot)
pl1<-read.table("Desktop Directory", sep=",", header=TRUE)
UP1<-pl1$Gene.name[match(up1, pl1$Gene.stable.ID)]
DOWN1<-pl1$Gene.name[match(down1, pl1$Gene.stable.ID)]
DP1<-compareCluster(entid_down, fun = "enrichGO", OrgDb="org.Hs.eg.db",  ont="BP", pvalueCutoff=0.1)
BP1<-compareCluster(entid20, fun = "enrichGO", OrgDb="org.Hs.eg.db",  ont="BP", pvalueCutoff=0.1)
MF1<-compareCluster(entids, fun = "enrichGO", OrgDb="org.Hs.eg.db",  ont="MF", pvalueCutoff=0.01)
CC1<-compareCluster(entids, fun = "enrichGO", OrgDb="org.Hs.eg.db",  ont="CC", pvalueCutoff=0.01)
BP1<-compareCluster(entid1, fun = "enrichGO", OrgDb="org.Hs.eg.db",  ont="BP", pvalueCutoff=0.1)
ggsave(dotplot(BP1, showCategory=12)+ggtitle(paste0("GO Analysis")), filename="", width=15, height=9)
