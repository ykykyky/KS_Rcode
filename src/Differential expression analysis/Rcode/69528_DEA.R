if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install()
n
BiocManager::install(c("GEOquery","limma","impute" ),ask = F,update = F)
library(GEOquery)

#Get dataset
eSet1 <- getGEO("GSE69528",
                destdir = '.',
                getGPL = F)

exp1 <- exprs(eSet1[[1]])
exp_raw <- exprs(eSet1[[1]])
exp1[1:4,1:4]
library(limma)  
exp1=normalizeBetweenArrays(exp1)  
pd1 <- pData(eSet1[[1]])
gpl1 <- eSet1[[1]]@annotation


#gene id conversion
gpl1
BiocManager::install(c("org.Hs.eg.db","illuminaHumanv4.db" ),ask = F,update = F)
library(illuminaHumanv4.db)
ids=toTable(illuminaHumanv4SYMBOL)
exp1 = exp1[rownames(exp1) %in% ids$probe_id,]
ids=ids[match(rownames(exp1),ids$probe_id),]
tmp = by(exp1,
         ids$symbol,
         function(x) rownames(x)[which.max(rowMeans(x))])
probes = as.character(tmp)
exp1 = exp1[rownames(exp1) %in% probes,]
rownames(exp1)=ids[match(rownames(exp1),ids$probe_id),2]

exp3=exp1[,c(1,9:17,19,21:27,33:36,42:47,56:66,69,73,75,88,91,92,104,107,111,114,118,120,127,129,133,134,136,137)]
group_list <- c(rep('control',28),rep('case',29))

design <- model.matrix(~0+factor(group_list))
colnames(design) <- levels(factor(group_list))
rownames(design) <-colnames(exp())
contrast.matrix<-makeContrasts(paste0(c("case","control"),collapse = "-"),levels = design)
contrast.matrix
fit <- lmFit(exp,design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)
tempOutput = topTable(fit2, coef=1, n=Inf)

nrDEG = na.omit(tempOutput)

BiocManager::install('clusterProfiler')
a
library(clusterProfiler)
colnames(nrDEG)
DEG=nrDEG
logFC_cutoff = 0.7
DEG$change = as.factor(ifelse(DEG$P.Value < 0.05 & abs(DEG$logFC) > logFC_cutoff,
                              ifelse(DEG$logFC > logFC_cutoff ,'UP','DOWN'),'NOT')
)

nrDEG_lfc<-DEG[which(DEG$P.Value<=0.05),]
write.csv(nrDEG_lfc,file="D:/69528_DEG.csv")




