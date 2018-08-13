library(rJava)
library(xlsxjars)
library(xlsx)
library(gplots)

data <- read.xlsx("WY20170823#19.xlsx",1,encoding="UTF-8")
row.names(data) <- data[,1] #将第一列设为行标题
data <- data[,-1] 
df <- as.matrix(scale(data)) #归一化，矩阵化
#df <- as.matrix(data)
#heatmap.2(df,scale = "none",col = bluered(100), trace = "none",density.info = "none" )
#heatmap(df, scale = "none", col=bluered(100),ColSideColors = c(rep("purple", 15), rep("orange", 15)))
heatmap.2(df,col = topo.colors(100),trace ="none",xlab ="Samples",lhei = c(2,8),density.info = "none",scale = c("none"),na.rm = TRUE,symbreaks = min(df, na.rm=TRUE),na.color="black",cexRow =0.9, cexCol =0.9, main ="Heatmap of whole CpG methylation sites",ColSideColors=c(rep("purple", 14), rep("orange", 14)),margins = c(10, 10))
heatmap.2(df,col = topo.colors(100),density.info = "none",xlab ="Samples",lhei = c(2,8),scale = c("none"),na.rm = TRUE,symbreaks = min(df, na.rm=TRUE),na.color="black",cexRow =0.9, cexCol =0.9, trace="none",main ="Heatmap of whole CpG methylation sites",ColSideColors=c(rep("purple", 14), rep("orange", 14)),Rowv=FALSE,Colv=FALSE,margins = c(10, 10))




col=heat.colors(100)
col=rev(heat.colors(16))
col=terrain.colors(256)
col= cm.colors(255)
col=rainbow(30)
col=redgreen(100)
col=redblue(100),
col=topo.colors(100)
col=colorpanel(100,low="white",high="steelblue")
col=colorRampPalette(c("yellow","darkblue"))  #提供渐变色