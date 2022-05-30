if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install()
n
BiocManager::install(c("GEOquery","limma","impute" ),ask = F,update = F)
library(GEOquery)

#2.Get dataset
eSet2 <- getGEO("GSE46955",
                destdir = '.',
                getGPL = F)
exp2 <- exprs(eSet2[[1]])
exp_raw <- exprs(eSet2[[1]])
exp2[1:4,1:4]
boxplot(exp2,las=2)
pd2 <- pData(eSet2[[1]])
gpl2 <- eSet2[[1]]@annotation

#4.gene id convert
gpl2
gpl<- getGEO('GPL6104', destdir=".")
colnames(Table(gpl)) #查一下列名
head(Table(gpl)[,c(1,12)])
ids2=Table(gpl)[,c(1,12)]
save(ids2,file='ids2.Rdata')
length(unique(ids2$Symbol))
tail(sort(table(ids2$Symbol)))
table(sort(table(ids2$Symbol)))
table(rownames(exp2) %in% ids2$ID)
exp2 = exp2[rownames(exp2) %in% ids2$ID,]
ids2=ids2[match(rownames(exp2),ids2$ID),]
tmp = by(exp2,
         ids2$Symbol,
         function(x) rownames(x)[which.max(rowMeans(x))])
probes = as.character(tmp)
exp2 = exp2[rownames(exp2) %in% probes,]
rownames(exp2)=ids2[match(rownames(exp2),ids2$ID),2]

exp4=exp2[,c(3,9,15,19,21,23,5,11,17,27,31,35,39,43)]

group_list <- c(rep('healthy_basal',6),rep('sepsis_basal',8))

exp <- exp4


library(limma)
design <- model.matrix(~0+factor(group_list))
colnames(design) <- levels(factor(group_list))
rownames(design) <-colnames(exp)
contrast.matrix<-makeContrasts(paste0(c("sepsis_basal","healthy_basal"),collapse = "-"),levels = design)
contrast.matrix
fit <- lmFit(exp,design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)
tempOutput = topTable(fit2, coef=1, n=Inf)

nrDEG = na.omit(tempOutput)


BiocManager::install('clusterProfiler')
a
library(clusterProfiler)
load(file = "DEGoutput.Rdata")
colnames(nrDEG)
DEG=nrDEG
logFC_cutoff =0.7
DEG$change = as.factor(ifelse(DEG$P.Value < 0.05 & abs(DEG$logFC) > logFC_cutoff,
                              ifelse(DEG$logFC > logFC_cutoff ,'UP','DOWN'),'NOT')
)

nrDEG_lfc<-DEG[which(DEG$P.Value<=0.05),]

write.csv(nrDEG_lfc,file="D:/46955_DEG.csv")











