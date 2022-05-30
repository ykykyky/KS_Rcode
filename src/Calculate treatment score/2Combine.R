 ######## collect all sizes` score for all drugs ############
 #收集所有药物的所有分数
################ 先合并（所有结果合并到一个文件中） ##############
#### 无显著性
 rm(list = ls())
# ra就是之前的top_size
 ra =c(12)
for(topsize in ra){
  # topsize = 990
  pathname = paste("D:\\RstudioWorkspace\\pathway\\ks-GBM\\mid\\D_",topsize, sep = "")
  setwd(pathname)
  filename <- list.files(pattern = '*_head_score.csv')
  #获取所有文件名称
  merge.data = read.csv(file = filename[1], header = FALSE, row.names = 1)
  #读取第一个文件
  merge.data = merge.data[-c(1),]
  #去除文件中的第一行
  for (i in 2:length(filename)){
    #从第2个文件开始读取
    new.data = read.csv(file = filename[i], header = FALSE, row.names = 1)
    new.data = new.data[-c(1),]
    merge.data = cbind(merge.data, new.data)
    #和前面的文件合并
  }
  score_fdr = merge.data
  #分数fdr
  score_fdr = as.matrix(score_fdr)
  #转化为矩阵
  ### 归一化
  head = score_fdr[c(1:7),]
  #取分数fdr的前7行（头信息）
  score = as.numeric(as.character(score_fdr[8,]))
  #将分数fdr的第8行转化为数值
  pos_p = which(score > 0)
  pos_n = which(score < 0)
  max = max(score)
  min = min(score)
  score[pos_p] = score[pos_p] / max
  score[pos_n] = -(score[pos_n] / min)
  head = rbind(head, score)
  #将头信息和分数合并
  # rownames
  name = as.character(row.names(score_fdr))
  row.names(head) = name
  score_fdr = head
  #设置score_fdr的行名
  # output
  #输出结果
  f = paste("all_score_", topsize, ".csv", sep = "")
  write.csv(score_fdr, f)
  # 下面是另一个目录（都在mid下）
  f = paste("D:\\RstudioWorkspace\\pathway\\ks-GBM\\mid\\all_score_", topsize, ".csv", sep = "")
  write.csv(score_fdr,f)
}

 