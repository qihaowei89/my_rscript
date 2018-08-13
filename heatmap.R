library(rJava)
library(xlsxjars)
library(xlsx)
library(gplots)

data <- read.xlsx("甲基化检测结果-LHY0066.xlsx",1,encoding="UTF-8")
row.names(data) <- data[,1] #将第一列设为行标题
data <- data[,-1] 
df <- as.matrix(scale(data)) #归一化，矩阵化

heatmap.2(df,scale = "none",col = bluered(100), trace = "none",density.info = "none" )