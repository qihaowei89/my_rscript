##################################################################################################
# 学习DESeq2包
# Data 2018-02-05
# Author Howard MENG
# E-mail menghaowei1992@126.com
##################################################################################################
rm(list=ls())

library(DESeq2)
###################################################################################
# 1.计算gene的raw reads count
###################################################################################
#### 方法1. 使用MENG Howard的代码，或者根据代码修改自己的需求；
  #优点：在R中完成操作，支持多线程；
  #缺点：内存占用比较大；

#### 方法2. 使用htseq-count程序
  #优点：节省内存
  #缺点：不够优雅，需要对读取进来的内容再做处理

#### 方法1,2本质没有差别，我们以方法1为例子进行演示


###################################################################################
# 2.读取计算结果
###################################################################################
#### 读取上一步的计算结果，如果用MENG的代码，直接就是生成SummarizedExperiment.obj
load(file="./ReadsCount.RData")
SummarizedExperiment.obj
colData(SummarizedExperiment.obj)

#### 构建DESeqDataSet
DES_set = DESeqDataSet(SummarizedExperiment.obj,design = ~case_type)

#### 进行差异表达分析
DEG_set.run = DESeq(DES_set)

#### 输出结果
DEG_set.result = results(DEG_set.run)
DEG_set.result

#### 对结果进行统计
summary(DEG_set.result)

#### 将结果保存成data.frame
DEG_set.result.df = as.data.frame(DEG_set.result)

library(org.Hs.eg.db)
columns(org.Hs.eg.db)

gene_symbol = mapIds(x=org.Hs.eg.db,
                     keys = rownames(DEG_set.result.df),
                     keytype = "ENTREZID",
                     column = "SYMBOL")

DEG_set.result.df.fix = data.frame(gene_symbol,DEG_set.result.df)

