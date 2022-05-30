rm(list = ls())

#gene id convert to entrez id
library('org.Hs.eg.db')
idlist=toTable(org.Hs.egSYMBOL)

#Read the gene list of the current pathway
deg = read.csv("C:\\Users\\yk\\Desktop\\TNF signaling pathway.csv", header = TRUE)


###########add LogFC to genes in pathway

#Read the gene expression matrix of the current dataset
fc = read.csv("C:\\Users\\yk\\Desktop\\\\46955_GEG.csv", header = TRUE)

sandf = fc[match(deg$Column1,fc$X),c(1,2)]
sandf$ENTREZID= idlist[match(sandf$X,idlist$symbol),1]
eandf = sandf[c(3,2)]
eandf=na.omit(eandf)
write.csv(eandf,file="C:\\Users\\yk\\Desktop\\46955\\pathway_logfc\\TNF signaling pathway.csv")


