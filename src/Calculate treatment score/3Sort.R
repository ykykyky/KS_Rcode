#####################################
## 对结果进行两种排序
## ①对结果的|ks|从大到小排序【这部分包含了治疗和加剧疾病的药物】
## ②只用ks<0的那部分数据（对其按照|ks|从大到小排序）【这部分是治疗的药物】
# size就是之前的top_size
rm(list = ls())

 size = c(8,11,25,56,62,63,67,68,71,86)

for(s in size){
  # s = 200
  f = paste("D:\\RstudioWorkspace\\pathway\\ks-GBM\\mid\\all_score_", s, ".csv", sep = "")
  data = read.csv(f, header = F)
  #读取数据
  data = data[-c(1),-c(1)]
  #去掉第一行和第一列
#  data = data[,-c(1)]
  data = as.matrix(data)
  #数据转化为矩阵
  len = length(data[8,])
  #len是第8行的长度（第8行是ks值）
  score = c(0)
  #分数列表
  for(i in 1:len){
    #遍历第8行
    score[i] = as.numeric(as.character(data[8,i]))
    #每一列得到一个分数
  }
  ########################################## 所有|ks|排序
  #data1 = data[,order(-abs(score))]
  #对数据按照|ks|排序
  #f = paste("D:\\RstudioWorkspace\\pathway\\ks-CHOL\\mid\\all_score_sort\\all_normal_abs_sort_", s,".csv", sep = "")
  #write.csv(data1, f,row.names = F, col.names = F)
  #写数据
  #########################################  不需要对绝对值排序
  
  # 只用ks<0的排序
  data1 = data[,order(score)]
  #按照分数对数据排序
  data1 = data1[,which(data1[8,] < 0)]
  #只取第8行<0的那些数据
  print(length(data1[3,]))
  f = paste("D:\\RstudioWorkspace\\pathway\\ks-GBM\\mid\\all_score_sort\\all_normal_less0_sort_", s,".csv", sep = "")
  write.csv(data1, f,row.names = F, col.names = F)
  #写数据
}
 ########### 对排序结果提取药物最优排序 #############
## 得到的药物排序就是最终的排名，接下来可以验证看top 10,20,30里有多少与疾病真正有关
# size就是top_size
 size = c(70)
 ## 需要分成两次（要修改文件名）
for(s in size){
  # s = 800
  f = paste("E:\\RstudioWorkspace\\pathway1\\mid\\all_score_sort\\all_normal_less0_sort_", s,".csv", sep = "")
#   f = paste("E:\\work\\chol\\result\\nodelete\\mid\\all_score_sort\\all_normal_abs_sort_", s,".csv", sep = "")
  data = read.csv(f, header = F)
  #读取文件
  data = data[-c(1),]
  #去掉第一行
  data = as.matrix(data)
  #数据转化为矩阵
  
  ## 一共有多少药
  drugs = as.character(data[4,])
  #取第4行的数据（第4行是药物名称）
  drugs = unique(drugs)
  #对数据去重
  print(length(drugs))
  #打印药物数量
  
  ## 按顺序提出药物最优排序
  result = data
  drugs = c("")
  for(i in 1:length(result[4,])){
    #遍历第4行所有药物
    drugs[i] = as.character(result[4,i])
    #将药物存放到列表里
  }
  drugs = unique(drugs)   #药物列表去重
  drugs = tolower(drugs)  #药物列表转化为小写
  print(length(drugs))    #打印药物列表长度
  of = paste("E:\\RstudioWorkspace\\pathway1\\mid\\all_score_", s, "_less0_sort_drugs.csv", sep = "")
#  of = paste("E:\\work\\chol\\result\\nodelete\\mid\\sort_drugs\\all_score_", s, "_abs_sort_drugs.csv", sep = "")
  write.csv(drugs, of, col.names = F, row.names = F)
  #写数据
}

