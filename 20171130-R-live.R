##################################################################
# 第0部分 安装R 并安装R Studio
##################################################################
# A good start ！
print("Hello World!")


##################################################################
# 第1部分 R的基础操作
##################################################################

## 常用的赋值操作
a <- 1
2 -> a
a <- 1e5
a <- 1e-10
a <- log2(10)
a <- log10(10)
a <- log(10)
a <- pi
a <- exp(1) 

## 常用的方法生成向量
b <- c(1:10)
a <- 1 
c <- a + b

## 创建序列 
d <- seq(from=1,to=100,by=5)
e <- rep(c(1:10),5)
f <- rep(c(1:10),each=5)

## 用c函数可以合并2个向量
a = c(1:10)
b = rep(c(1:10),each=5)
c = rep(c(1:5),2)
d = c(a,b,c)

## 创建矩阵
g <- matrix(c(1:100), nrow = 10)
h <- matrix(c(1:15),5,5)

g_up <- upper.tri(g) #上三角
g_low <- lower.tri(g) #下三角
det(h) # 求方阵行列式
t(g) # 转置
eigen(g)  #特征值

## 创建data.frame
gene_id = c(1:100)
gene_fpkm = rnorm(100,10,5)
gene_fpkm = abs(gene_fpkm)

gene_table = data.frame(gene_id,
                        gene_fpkm)
plot(x=gene_table$gene_id,y=gene_table$gene_fpkm,type = "o")


# 查看变量空间中都有哪些变量
ls()

# 清除所有环境空间内的数据
rm(list=ls()) 


## R中的if判断
## R中的for循环
d <- seq(from=1,to=100,by=5)
for(i in d){
  if(i >= 50){
    print(i)
  }
}

d[d>50 & d<80] #同样的

# 对于数据框，也同样可以做筛选
gene_table
gene_table.select = gene_table[gene_table$gene_fpkm >= 10,]

## R中写自己的函数
my_ABS <- function(x){
  if(x <= 0){
    -x
  }else{
    x
  }
}

my_ABS(10)
my_ABS(-10)

##################################################################
# 第2部分 R基础绘图系统初探
##################################################################
# R与统计分析
case_1 <- rnorm(50,mean = 20,sd=5) #生成正态分布样本1
plot(case_1,main="case_1")
case_2 <- rnorm(40,mean = 25,sd=5) #生成正态分布样本2
plot(case_2,main="case_2")

hist(case_1,breaks = c(1:40)) #分布图
hist(case_2,breaks = c(1:40))

case_T.test <- t.test(case_1,case_2) #独立样本T检验
case_T.test$p.value

# R画图
## 使用R画boxplot
case_1 <- rnorm(10,5,1)
case_2 <- rnorm(10,6,1)
case_3 <- rnorm(10,7,1)
case_4 <- rnorm(10,8,1)
case_all <- cbind(case_1,case_2,case_3,case_4)
boxplot(case_all)

col_1 <- "red"
col_2 <- rainbow(10)[5]
col_3 <- rgb(1,0,0,alpha = 0.5)
col_4 <- rgb(0,0,1,alpha = 0.5)
col_list <- c(col_1,col_2,col_3,col_4)
boxplot(case_all,col = col_list)

col_list <- c(col_3,col_3,col_4,col_4)
boxplot(case_all,col = col_list)

##################################################################
# 第3部分 绘制RNA-Seq基因表达量的散点图
##################################################################
rm(list=ls())
cuffnorm_result = read.csv(file="~/Desktop/live_R_data/cuffnorm_genes_fpkm.csv")
head(cuffnorm_result)

x.vector = cuffnorm_result$Empty_KD_0
y.vector = cuffnorm_result$Empty_KD_1

plot(x=x.vector,y=y.vector,col="#0000FF11",pch=16,cex=0.5,xlim=c(0,10),ylim=c(0,10),xlab = "Repeat_1 Log2(FPKM)",ylab = "Repeat_2 Log2(FPKM)")
text(2,8,round(cor(x.vector.filter,y.vector.filter,method = "spearman"),4))
hist(x.vector.filter[x.vector.filter>=0],breaks = seq(0,15,0.1),border = F,col = "blue",xlab = "FPKM",ylab = "Frequency",main = "")
hist(y.vector.filter[y.vector.filter>=0],breaks = seq(0,15,0.1),border = F,col = "blue",xlab = "FPKM",ylab = "Frequency",main = "")

##################################################################
# 第4部分 绘制多维度的测序仪性能图
##################################################################



##################################################################
# 第5部分 修改了绘图系统的heatmap绘图
##################################################################
library(pheatmap)
test = matrix(rnorm(200), 20, 10)
test[1:10, seq(1, 10, 2)] = test[1:10, seq(1, 10, 2)] + 3 
test[11:20, seq(2, 10, 2)] = test[11:20, seq(2, 10, 2)] + 2 
test[15:20, seq(2, 10, 2)] = test[15:20, seq(2, 10, 2)] + 4 
colnames(test) = paste("Test", 1:10, sep = "")
rownames(test) = paste("Gene", 1:20, sep = "")

# 设置每一列的注释
annotation_col = data.frame(CellType = factor(rep(c("CT1", "CT2"), 5)),Time = 1:5)
rownames(annotation_col) = paste("Test", 1:10, sep = "") # 设置每一行的注释

annotation_row = data.frame(GeneClass = factor(rep(c("Path1", "Path2", "Path3"), c(10, 4, 6))))
rownames(annotation_row) = paste("Gene", 1:20, sep = "") # 设置注释的颜色

ann_colors = list(Time = c("white", "firebrick"),
                  CellType = c(CT1 = "#1B9E77", CT2 = "#D95F02"),
                  GeneClass = c(Path1 = "#7570B3", Path2 = "#E7298A", Path3 = "#66A61E")
                  ) 

pheatmap(test,
         annotation_col = annotation_col,
         annotation_row = annotation_row,
         annotation_colors = ann_colors)






