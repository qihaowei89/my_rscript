source("http://www.bioconductor.org/biocLite.R")
biocLite("limma")
biocLite("GEOquery")
biocLite("gcrma")
biocLite("ArrayExpress")
library(limma)
library(GEOquery)
library(affy)
library(gcrma)
 
# 下载基因芯片数据，destdir参数指定下载到本地的地址
gse382 <- getGEO('GSE49382', destdir = ".") ##根据GSE号来下载数据，下载_series_matrix.txt.gz
gpl515 <- getGEO('GPL17515', destdir = ".")##根据GPL号下载的是芯片设计的信息, soft文件

###已经下载好的数据从此开始##########
# 打开已下载的本地数据
gse382 <- getGEO(filename = 'GSE49382_series_matrix.txt.gz')
gpl515 <- getGEO(filename = 'GPL17515.soft')

###################################################################################################################
# 【数据分析】用limma进行差异表达分析
# 自己做好三个数据矩阵（表达矩阵，分组矩阵，差异比较矩阵），然后limma的三个步骤（lmFit,eBayes,topTable）就可以啦
###################################################################################################################
# 查看列名

colnames(Table(gpl515))
head(Table(gpl515))
write.csv(Table(gpl515)[,c(1,4)],"GPL515.csv", row.names = F)

# gse382中的行名ID与gene name的对应关系
genename = read.csv("GPL515.csv")

# 构建表达矩阵
exprSet <- as.data.frame(exprs(gse382)) # 得到表达矩阵，行名为ID，需要转换

# 转换ID为gene name
exprSet$ID = rownames(exprSet)
express = merge(x=exprSet, y=genename, by="ID", all.x = T)
express$ID = NULL
rowMeans=apply(express,1,function(x) mean(as.numeric(x),na.rm=T))

express = express[!duplicated(express[,28]),] #express第28列为gene name

rownames(express) = express[,28]
express=express[,-28]

# 构建分组矩阵
pdata = pData(gse382) # 每个sample所对应的信息，包括处理条件等
group_list = subset(pdata, select=title) # Sample的分组信息

group_list$condition = rep(c("c0","h0","r0","c1","h1","r1","c6","h6","r6"), each=3)

design = model.matrix(~0+factor(group_list$condition))

colnames(design) = levels(factor(group_list$condition))
rownames(design) = colnames(express)
# 至此，分组矩阵（design）已构建好


# 构建差异比较矩阵
contrast.matrix = makeContrasts(c0-c1,c0-c6,h0-h1,h0-h6,r0-r1,r0-r6,levels = design)
# 至此，差异表达矩阵已构建好

fit = lmFit(express,design)

fit2 = contrasts.fit(fit, contrast.matrix)

fit2 = eBayes(fit2)


# 得到两两差异表达的结果
# c0 vs. c1
x = topTable(fit2, coef = 1, n=Inf, adjust.method = "BH", sort.by = "P")
sum(x$adj.P.Val<0.05,na.rm = T)
re = x[which(x$adj.P.Val  < 0.05& (x$logFC > 1 | x$logFC < 1)),]   #选取adj.p.value<0.05且|logFC|>1的基因
write.csv(re, "c0-c1_DEG_limma.re.csv",quote = F)
# coef可是column number，也可以是column name，这样就可以指定你所感兴趣的两两比较的结果
# 在此例中coef =1 就是c0-c1的差异表达比较结果

# 查看差异表达结果分组情况
results = decideTests(fit2,p.value = 0.05)








