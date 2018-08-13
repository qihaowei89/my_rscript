library(xlsx)

dat  <- read.xlsx( "สýพÝิด.xlsx", 1, as.data.frame = TRUE, header = TRUE )
aa <- dat[,-1]
row.names(aa) <- aa[,1]
data <- aa[,-1]
t.data <- t(data)

data.pca <- prcomp(t.data, scale.= T)
pca.sum=summary(data.pca)
#pdf(file="$outdir/$name.pcaPlot.pdf")
plot(data.pca$x,
     main= 'PCA',
     xlab= "PC1", 
     ylab= "PC2",
     type= "p",pch=16,
     xlim=c(1.2*min(data.pca$x),1.2*max(data.pca$x)),
     ylim=c(1.2*min(data.pca$x),1.2*max(data.pca$x)),
     col=c(rep("#FF0000",5),rep("#0000FF",5),rep("#8B0A50",5)),
     cex=2)
text(x= data.pca$x[,1], y= data.pca$x[,2], labels= rownames(data.pca$x), cex= 1, col=c(rep("#FF0000",5),rep("#0000FF",5),rep("#8B0A50",5)),pos=4)




