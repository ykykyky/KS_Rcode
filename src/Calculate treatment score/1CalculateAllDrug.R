#   calculate score and p-value for all drug(size = 978 genes)                   
#       disease : TCGA diff_sig   drug: LINCS(diff idoses)          
#计算所有药物的分数和p值（978个基因）
#疾病：TCGA diff_sig   药物：LINCS

# CMap drug signature follows 
cmap_score <- function(sig_up, sig_down, drug_signature) {
  # calculate interaction score between two sets of top-ranked genes and a
  # drug-induced expression signature.
  # Inputs:
  #   sig_up - numeric of GeneIDs from top end of ranked list (lo expression)
  #   sig_down - numeric GeneIDs from bottom end of ranked list (hi expression)
  #   drug_signature - data.frame GeneIDs and rank
  # Returns:
  #   connectivity score
  
  #  sig_up = ms_up
  #  sig_down = ms_down
  
  #计算相互作用分数：两个top排序基因集合和一个药物介导的表达标签之间的相互作用分数
  #输入：
  #      sig_up:需要进行KS检验的所有上调基因的id（低表达）
  #      sig_down:需要进行KS检验的所有下调基因的id（高表达）
  #      drug_signature：数据框GeneID和顺序
  #返回值：
  #      关联度分数
  
  num_genes <- nrow(drug_signature)
  #num_genes:药物标签的行数
  ks_up <- 0
  ks_down <- 0
  connectivity_score <- 0
  #数据初始化
  
  # re-rank because the GeneID mapping changed the original rank range
  drug_signature[, "rank"] <- rank(drug_signature[, "value"])           #对表达值从小到大排序
  #药物标签中新加了一个属性“名次”，这个名次就是按表达值从小到大排序
  
  # Merge the drug signature with the disease signature by GeneID. This becomes
  # the V(j) from the algorithm description.
  #合并药物标签和疾病标签（根据GeneID）这变成了算法描述中的V(j)
  up_tags_rank <- merge(drug_signature, sig_up, by.x="GeneID", by.y=1)   #up_tags_rank的排序与之前两个都不一样
  #  print("up 的交集")
  #  print(up_tags_rank)
  down_tags_rank <- merge(drug_signature, sig_down, by.x="GeneID", by.y=1)
  #  print("bottom 的交集")
  #  print(down_tags_rank)
  
  up_tags_position <- sort(up_tags_rank$rank)
  #获取公共基因数据在药物数据当中的排名
  #  print("up position")
  #  print(up_tags_position)
  down_tags_position <- sort(down_tags_rank$rank)
  #获取公共基因数据在药物数据当中的排名
  #  print("down position")
  #  print(down_tags_position)
  
  num_tags_up <- length(up_tags_position)
  #获取公共基因的个数
  num_tags_down <- length(down_tags_position)
  #获取公共基因的个数
  
  if(num_tags_up > 1 && num_tags_down > 1) {
    #如果上调和下调的公共基因数量都大于1
    a_up <- 0
    b_up <- 0
    
    # small speed up by changing from sapply to inline math (~0.5 sec)
    a_up <- max(((1:num_tags_up) / num_tags_up) - (up_tags_position / num_genes))
    b_up <- max((up_tags_position / num_genes) - (((1:num_tags_up)-1) / num_tags_up))
    #up_tags_position是一个列表，
    
    if(a_up > b_up) {
      ks_up <- a_up
    } else {
      ks_up <- -b_up
    }
    
    a_down <- 0
    b_down <- 0
    
    # small speed up by changing from sapply to inline math (~0.5 sec)
    a_down <- max(((1:num_tags_down) / num_tags_down) - (down_tags_position / num_genes))
    b_down <- max((down_tags_position / num_genes) - (((1:num_tags_down)-1) / num_tags_down))
    
    if(a_down > b_down) {
      ks_down <- a_down
    } else {
      ks_down <- -b_down
    }
    
    if (sum(sign(c(ks_down, ks_up))) == 0) {   #一正一负
      connectivity_score <- ks_up - ks_down # different signs
    }
  }
  return(connectivity_score)
  #返回关联度
}

######## read data ##############################################################
#读取数据
######## read disease data ########
#读取疾病数据
#疾病相关通路包含的基因集的ENTREZID和logFC组成的数据文件
MS <- read.csv("D:\\kegg-TCGA\\TCGA-ESCA\\ks-KIRP\\Aldosterone-regulated sodium reabsorption.csv", header = TRUE)

### use for calculate rand_scores  
#用于计算随机分数？
## not rank
#未排序
geneid <- MS$ENTREZID
#geneid:MS中的ENTREZID
zscore <- MS$logFC
#zscore:MS中的logFC
ZData <- data.frame(geneid,zscore)
#ZData:基因id和差异倍数组成的数据框
colnames(ZData) <- c("GeneID","value")
#设置列名
### delete NA``
ZData = na.omit(ZData)
#删除空缺值

### use for calculate score
#用来计算分数
### sort(large->small)
#从大到小排序
d <- MS[order(MS$logFC),]
#按基因差异倍数从大到小排序
geneid <- d$ENTREZID
#排序之后的geneid
zscore <- d$logFC
#排序之后的差异倍数
ZData1 <- data.frame(geneid,zscore)
#排序之后的geneid和差异倍数
colnames(ZData1) <- c("GeneID","value")
设置列名
### delete NA
ZData1 = na.omit(ZData1)
#删除空缺值

#### 读取药物文件---> 计算ks_score，pvalue
setwd("D:\\LINCS\\fenjie")
#设置当前工作文件夹
temp <- list.files(pattern = '*.csv')
#将工作目录下所有csv文件的文件名存到temp中

# all diff genes used
#用到所有差异基因
 top_size = c(4)#-个数
 #选择哪些上下调基因参与KS检验
 bottom_size = c(15)#+个数
 #选择哪些上下调基因参与KS检验
# top_size = c(100,200,300)
# bottom_size = c(32,32,32)

for(size in 1:length(top_size)){
  #遍历需要进行KS检验的上/下调基因
  # size = 1
  #从第一个基因开始
  for(drug_size in 1:length(temp)){
    #遍历所有文件名（文件名是某一种药物名称）
    # drug_size = 1
    CMapData <- read.csv(temp[drug_size], header = FALSE)
    #读取第drug_size个文件（CMapData是某种药物不同处理的文件数据？）
    CMapData <- CMapData[,-c(2)]  #去掉第2列数据（DMSO）
    head <- CMapData[c(1:7),]   #只取前7行的数据
    CMapData <- CMapData[-c(1:7),]#去掉前7行的数据
    colnames(CMapData) = c("GeneID")#将去除前7行数据的剩下数据加上列名：GeneID
    ### pull out matching_genes
    #选出匹配基因？
    matching_genes = unique(merge(CMapData, ZData, by.x = "GeneID")[,1])
    #将基因数据和差异倍数信息合并
    print(length(matching_genes)) # 1286
    #  print(length(which(ZData1$value > 0)))
    #打印出匹配的基因数
  
    N = ncol(CMapData) - 1
    score=numeric(N)
    score_pvalue=numeric(N)
    #score和score的p_value都设置为长度为N（当前药物有N中处理）的列表
    
    ### for each RA-drug with different idose and different cell line
    #对每一个不同剂量和不同细胞系的药物
    for(i in 1:N){
      # i=1
      drug_signature = CMapData[,c(1, i+1)]
      #药物标签：取第1列（基因ID）和第i列数据（基因表达值）
      colnames(drug_signature) = c("GeneID","value")
      #药物标签的列名设置为基因id和value值
      drug_signature = as.data.frame(drug_signature)
      #转化为数据框
      
      ### select top and bottom of rank and collect enrichment score
      #选择名次顶端和底端，收集富集分数
      ms_up = ZData1$GeneID[1:top_size[size]]
      #需要进行KS检验的所有上调基因的id
      gsel = (nrow(ZData1) - bottom_size[size] + 1):nrow(ZData1)
      ms_down = ZData1$GeneID[gsel]
      #需要进行KS检验的所有下调基因的id
      
      score[i] = cmap_score(ms_up, ms_down, drug_signature)
      #当前药物第i个处理的分数
      
    } # end for
    ### write score and p-value
    #写分数和p_value
    file1 = strsplit(temp[drug_size],'.csv')
    #取出当前药物文件名中的药物名（无后缀）
    ########score.csv and pvalue.csv both have a col names = "x"########
    #分数文件和pvalue文件都有同一个列名
    ofile = paste("D:\\RstudioWorkspace\\pathway\\mid\\D_",top_size[size],"\\", file1[[1]], "_score.csv", sep = "")
    #创建保存当前药物score的结果文件
    write.csv(score, file = ofile, quote = F, col.names = F, row.names = F)
    #将score列表（当前药物所有处理的分数）写入对应的文件中
    score = t(score)  # 转置分数列表（一列数据变成一行数据） 
    rownames(score) = c("score")   #给分数加上行名
    name = as.character(head[,1])  #name取head的第一列，即行名
    name[8] = "score"     #第8个字段设为score
    head <- head[,-c(1)]  #去除第一列（去除行名）
    head = as.matrix(head)#转化为矩阵
    score_fdr = as.matrix(score)#得到分数的fdr
    data <- rbind(head,score_fdr)#将head和score按行合并
    data <- cbind(name,data)     #将name和data(head和score)按列合并
    
    ofile = paste("D:\\RstudioWorkspace\\pathway\\mid\\D_",top_size[size],"\\", file1[[1]], "_head_score.csv", sep = "")
    #创建保存当前药物score（添加了头信息）的结果文件
    write.csv(data, ofile,  col.names = F, row.names = F)
    #写入数据
  }
}
