library(ggplot2)
library(xlsx)

dat  <- read.xlsx( "数据源2.xlsx",2, as.data.frame = TRUE, header = TRUE )
aa <- dat[,-1]
row.names(aa) <- aa[,1]
data <- aa[,-1]
t.data <- t(data)   #转置数据

data.pca <- prcomp(t.data, scale.= T)   #主成分分析
pca.sum=summary(data.pca)
Eigenvalue <- round(data.pca$sdev^2,2) #获得特征值
#Eigenvector <- data.pca$rotation[,i]   #获得特征向量
Var <- pca.sum$importance[2,]*100
rbind(Eigenvalue,Var) -> PCA_analysis
row.names(PCA_analysis) <- c("Eigenvalue","% Var")
##其他计算Eigenvalues方法
#library(factoextar)
#eig.val <- get_eigenvalue(data.pca)
## eigenvalues  可视化
#fviz_eig(data.pca)
write.table(PCA_analysis,"PCA_analysis.txt",sep = "\t",quote = F)
data.pca$x -> a
as.data.frame(a) -> b

#Group <- c(rep(c("Yi"),5),rep(c("Ya"),5),rep(c("P"),4))
Group <- c(rep(c("Ya"),5),rep(c("P"),4))
cbind(b,Group ) -> b_g
pc1 <-paste("PC1","(",round(pca.sum$importance[2,]["PC1"],4)*100,"%",")",sep = "")
pc2 <-paste("PC2","(",round(pca.sum$importance[2,]["PC2"],4)*100,"%",")",sep = "")

pdf("PCA_阴虚-平和-聚类-63.pdf")


ggplot(b_g)+
    aes(x=PC1,y=PC2,colour=Group,label = row.names(b_g))+
    geom_point(aes(x=PC1,y=PC2,colour=Group),size = 5)+
    geom_text(vjust = 0,nudge_x = 1,size=3)+
    labs(x = pc1,y = pc2)

dev.off()

png("PCA_阴虚-平和-聚类-63.png")
ggplot(b_g)+
    aes(x=PC1,y=PC2,colour=Group,label = row.names(b_g))+
    geom_point(aes(x=PC1,y=PC2,colour=Group),size = 5)+
    geom_text(vjust = 0,nudge_x = 1,size=3)+
    labs(x = pc1,y = pc2)

dev.off()









