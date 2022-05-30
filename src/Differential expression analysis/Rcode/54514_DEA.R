if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install()
a
BiocManager::install(c("GEOquery","limma","impute" ),ask = F,update = F)
BiocManager::install(c("org.Hs.eg.db","illuminaHumanv3.db" ),ask = F,update = F)
library(GEOquery)
eSet <- getGEO("GSE54514", 
               destdir = '.',  
               getGPL = F) 
raw_exp <- exprs(eSet[[1]])
exp <- exprs(eSet[[1]])
library(illuminaHumanv3.db)
ids=toTable(illuminaHumanv3SYMBOL)
exp = exp[rownames(exp) %in% ids$probe_id,]
ids=ids[match(rownames(exp),ids$probe_id),]
tmp = by(exp,
         ids$symbol,
         function(x) rownames(x)[which.max(rowMeans(x))])
probes = as.character(tmp)
exp = exp[rownames(exp) %in% probes,] 
rownames(exp)=ids[match(rownames(exp),ids$probe_id),2]
pd <- pData(eSet[[1]])   
library(stringr)
exp2=exp[,c(19:49,50:145)]

group_list=c(rep("non",times=31),rep("sur",times=96))


library(limma)
design <- model.matrix(~0+factor(group_list))
colnames(design)=levels(factor(group_list))
rownames(design)=colnames(exp2)
contrast.matrix<-makeContrasts(paste0(c("non","sur"),collapse = "-"),levels = design)
contrast.matrix
fit <- lmFit(exp2,design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)
tempOutput = topTable(fit2, coef=1, n=Inf)
nrDEG = na.omit(tempOutput) 


BiocManager::install('clusterProfiler')
n
library(clusterProfiler)

options(stringsAsFactors = F)
load(file = "DEGoutput.Rdata")
colnames(nrDEG)

DEG=nrDEG

logFC_cutoff =0.7
DEG$change = as.factor(ifelse(DEG$P.Value < 0.05 & abs(DEG$logFC) > logFC_cutoff,
                              ifelse(DEG$logFC > logFC_cutoff ,'UP','DOWN'),'NOT')
)

nrDEG_lfc<-DEG[which(DEG$P.Value<=0.05),]

write.csv(nrDEG_lfc,file="D:/54514_DEG.csv")


