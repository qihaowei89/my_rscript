#########################
## mice 包学习笔记
## Date : 2018-2-28
#########################
data <- airquality
data0 <- airquality
data[4:10,3] <- rep(NA,7)
data[1:5,4] <- NA
data0 <- data0[-c(5,6)]
data <- data[-c(5,6)]

pMiss <- function(x) {sum(is.na(x))/length(x)*100}
apply(data,2,pMiss)
apply(data,1,pMiss)


library(VIM)
aggr_plot <- aggr(data, col = c('navyblue', 'red'), numbers=TRUE, sortVars=TRUE, 
                  labels=names(data), cex.axis=.7, gap=3, 
                  ylab=c("Histogram of missing data", "Pattern"))

tempData <- mice(data,m=5,maxit=50,meth='pmm',seed=500)
completedData <- complete(tempData,1)

#可视化
densityplot(tempData,~Ozone|.imp)
densityplot(tempData,~Solar.R|.imp)
densityplot(tempData,~Wind|.imp)
densityplot(tempData,~Temp|.imp)
stripplot(tempData, pch = 20, cex = 1.2)

modelFit1 <- with(tempData,lm(Temp~ Ozone+Solar.R+Wind))
modelFit2 <- with(data0,lm(Temp~ Ozone+Solar.R+Wind))

summary(pool(modelFit1))
summary(modelFit1)

