library(gplots)
library(xlsx)

dat  <- read.xlsx("数据源2.xlsx",2,as.data.frame=TRUE,header=TRUE)

dat1 <- dat[,-1]

row.names(dat1) <- dat1[,1]

dat2 <- as.matrix(dat1[,-1])
#dat2 <- scale(dat2)

#heatmap.2(dat2,col=greenred(75),density.info = "none",xlab ="Samples",lhei = c(2,8),na.rm = TRUE,symbreaks = min(dat2, na.rm=TRUE),na.color="black",cexRow =0.9, cexCol =0.9, trace="none",main ="Heatmap of whole CpG methylation sites",Rowv=F,Colv=F,margins = c(10, 10))
pdf("heatmap_阴虚-平和-聚类-63.pdf")

heatmap.2(
          dat2, 
          col=greenred(75),
          #col=bluered,
          symkey=FALSE, 
          scale = c("row"),
          density.info="none", 
          trace="none",
          cexRow = 1,
          keysize = 1,
          #key.title=NA,
          key.xlab="Score",
          #ColSideColors = c(rep("#3B3B3B",5),rep("#B22222",5),rep("#ADFF2F",4)),
          ColSideColors = c(rep("#3B3B3B",5),rep("#ADFF2F",4)),
          srtRow=45, adjRow=c(0, 1), srtCol=45, adjCol=c(1,1),  ## Show effect of row and column label rotation
          offsetRow=10, offsetCol=1,    #离坐标轴距离
          #reorderfun=function(d,w) reorder(d, w, agglo.FUN=mean), 
          #distfun=function(x) as.dist(1-cor(t(x))), 
          #hclustfun=function(x) hclust(x, method="complete")
          reorderfun=function(d,w) reorder(d, w, agglo.FUN=mean), 
          distfun=function(x) dist(x, method="euclidean"), 
          hclustfun=function(x) hclust(x, method="ward.D2")
          )
dev.off()



