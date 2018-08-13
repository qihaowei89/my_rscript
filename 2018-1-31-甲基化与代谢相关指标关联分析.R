##############################
#Data:2018-1-31
#qihao
#甲基化与相关代谢指标关联分析
#############################

library(xlsx)
library(psych)

################定义数据过滤函数和相关分析函数####
rm_na <- function(raw){
    raw <- data.frame(apply(raw,2,function(n) as.numeric(n)))
    temp <- raw[,-1]
    ID <- raw[,1]
    row_call <- apply(temp,1,function(x)  {(1-length(which(is.na(x)))/length(colnames(temp)))})
    temp <- temp[!row_call==0,]
    ID <- ID[!row_call==0]
    col_call <- apply(temp,2,function(x)  {(1-length(which(is.na(x)))/length(rownames(temp)))})
    data <- data.frame(编号=ID[row_call[!row_call==0]>0.8],temp[row_call[!row_call==0]>0.8, col_call > 0.7])
}
cor_test <- function(file_A,file_B,A2B){
    file_A <- rm_na(file_A )
    file_B <- rm_na(file_B )
    use <- intersect(file_A[,1],file_B[,1])
    result <- corr.test(file_A[file_A[,1]%in%use,-1],file_B[file_B[,1]%in%use,-1],method = "spearman")
    frame <- data.frame(result$r,Null=rep("|",nrow(result$r)),result$p)
    A2B <- paste(A2B,".csv",sep = "")
    write.csv(frame,A2B,quote = F)
}    

################主体############
file_name <- "C:/Users/weiqi/Desktop/甲基化检测结果_XJJ0039/XJJ0039(汇总).xlsx"
file <- read.xlsx(file_name,1,header = T,encoding = "UTF-8",check.names=F)

intake_VB    <- file[,c(1,19:24)]
metabolism   <- file[,c(1,12:18)]
OST          <- file[,c(1,7 :11)]
cg12166806_3 <- file[,c(1,25:41)]
cg01596986_7 <- file[,c(1,42:66)]
cg25073435_3 <- file[,c(1,67:74)]

###############维生素摄入与甲基化的关联分分析####

intake_VB_cg12166806_3 <- cor_test(intake_VB,cg12166806_3,"维生素摄入与cg12166806#3")
intake_VB_cg01596986_7 <- cor_test(intake_VB,cg01596986_7,"维生素摄入与cg01596986#7")
intake_VB_cg25073435_3 <- cor_test(intake_VB,cg25073435_3,"维生素摄入与cg25073435#3")


###############  代谢指标与甲基化的关联分析  ####

metabolism_cg12166806_3 <- cor_test(metabolism ,cg12166806_3,"代谢指标与cg12166806#3")
metabolism_cg01596986_7 <- cor_test(metabolism ,cg01596986_7,"代谢指标与cg01596986#7")
metabolism_cg25073435_3 <- cor_test(metabolism ,cg25073435_3,"代谢指标与cg25073435#3")

###############  氧化应激与甲基化的关联分析  ####

OST_cg12166806_3 <- cor_test(OST ,cg12166806_3,"氧化应激与cg12166806#3")
OST_cg01596986_7 <- cor_test(OST ,cg01596986_7,"氧化应激与cg01596986#7")
OST_cg25073435_3 <- cor_test(OST ,cg25073435_3,"氧化应激与cg25073435#3")




