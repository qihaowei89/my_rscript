#读取文件
library(reshape)
setwd("C:/Users/weiqi/Desktop/新建文件夹 (5)/")
#file <- read.table("年龄-性别.txt",header = T,encoding = "GB2312",stringsAsFactors = F)
#file <- read.table("病例-对照.txt",header = T,encoding = "GB2312",stringsAsFactors = F)
#file <- read.table("高白-低白-正常.txt",header = T,encoding = "GB2312",stringsAsFactors = F)
file <- read.table("T-B-正常.txt",header = T,encoding = "GB2312",stringsAsFactors = F)
file$积分 <- as.factor(file$积分)
file$总体甲基化水平 <- as.factor(file$总体甲基化水平)

# file$组别[file$组别=="正常"] <- 0
# file$组别[file$组别=="ALL低白"] <- 1
# file$组别[file$组别=="ALL高白"] <- 2
# file$组别 <- factor(file$组别,levels = c(0,1,2),labels = c("正常","ALL低白","ALL高白"))
file$组别[file$组别=="正常"] <- 0
file$组别[file$组别=="T-ALL"] <- 1
file$组别[file$组别=="B-ALL高危" | file$组别=="B-ALL中危" | file$组别=="B-ALL低危"] <- 2
file$组别 <- factor(file$组别,levels = c(0,1,2),labels = c("正常","T-ALL","B-ALL"))

file_chi <- file[,c(1,3,5)]
table1 <- with(file_chi,table(组别,积分))
table2 <- with(file_chi,table(组别,总体甲基化水平))
# file$组别[file$组别=="B-ALL低危"] <- 0
# file$组别[file$组别=="B-ALL中危"] <- 1
# file$组别[file$组别=="B-ALL高危"] <- 2
# file$组别 <- factor(file$组别,levels = c(0,1,2),labels = c("B-ALL低危","B-ALL中危","B-ALL高危"))
# file <- file[which(file$组别=="B-ALL低危"|file$组别=="B-ALL中危"|file$组别=="B-ALL高危"),]
file1 <- file[,-c(2,3,4,5)]
file2 <- melt(file1)
result <- list()
site <- levels(file2$variable)
for (i in site) {
    fit <- aov(value~组别,data=subset(file2,variable==i)[,-2])
    result_a <- summary(fit)
    Pr_F  <- unlist(summary(fit))["Pr(>F)1"]
    result_a1 <- TukeyHSD(fit)
    a1 <- result_a1$组别[,"p adj"]
    row_name <- c("方差分析 p_value",paste(rownames(result_a1$组别)," p_adj",""))
    temp <- cbind(row_name,p_value=c(Pr_F , a1))
    result <-cbind(result,temp[,2]) 
    result <- apply(result,2,function(n) as.numeric(n) )
    row.names(result) <- row_name
}
colnames(result) <- colnames(file1[,-1])



# sink("高白_vs_低白_vs_正常_方差分析.txt")
# print(result)
# sink()

sink("B-ALL_高_中_低_方差分析.txt")
print(result)
sink()

# file$分组[file$分组=="ALL缓解"] <- 0
# file$分组[file$分组=="B-ALL难治复发"] <- 1
# file$分组 <- factor(file$分组,levels = c(0,1),labels = c("ALL缓解","B-ALL难治复发")) 
# # file$性别分组[file$性别分组=="男"] <- 0
# file$性别分组[file$性别分组=="女"] <- 1
# file$性别分组 <- factor(file$性别分组,levels = c(0,1),labels = c("男","女"))
# file$年龄分组[file$年龄分组=="＜8"] <- 1
# file$年龄分组[file$年龄分组=="8至15"] <- 2 
# file$年龄分组 <- factor(file$年龄分组,levels = c(1,2),labels = c("＜8","8至15"))   
# file$分组[file$分组=="正常"] <- 0
# #file$疾病分组[file$疾病分组=="病例"] <- 1 
# file$分组[file$分组=="ALL病例"] <- 1 
# file$分组 <- factor(file$分组,levels = c(0,1),labels = c("正常","ALL病例"))  
##
B  <- file[file$分组=="B-ALL难治复发",]
re <- file[file$分组=="ALL缓解",]
##
# case <- file[file$分组 == "ALL病例",]
# control <- file[file$分组 == "正常",]
##
# all_male       <- file[file$性别分组== "男",]
# all_female     <- file[file$性别分组== "女",]
# all_male.control <- file[file$性别分组== "男"&file$疾病分组=="正常",]
# all_female.control <- file[file$性别分组== "女"&file$疾病分组=="正常",]
# age_less_8 <- file[file$年龄分组 == "＜8",]
# age_big_8 <- file[file$年龄分组 == "8至15",]
# age_less_8.control <- file[file$年龄分组 == "＜8"&file$疾病分组=="正常",]
# age_big_8.control <- file[file$年龄分组 == "8至15"&file$疾病分组=="正常",]

###各甲基化位点差异

# male_and_female.diff <- test(all_male,all_female,"男ALL_vs_女ALL")                                  #分组信息：男ALL   女ALL
# control_male_and_female.diff <- test(all_male.control,all_female.control,"男正常_vs_女正常")        #分组信息：男正常  女正常
# age_8.diff <- test(age_less_8, age_big_8,"＜8_vs_8至15")                                            #分组信息：＜8     8至15
# control_age_8.diff <- test(age_less_8.control, age_big_8.control,"＜8正常_vs_8至15正常")            #分组信息：＜8正常 8至15正常
# case_and_control <- test(case,control,"ALL病例_VS_正常")
B_and_re <- test(B,re,"B-ALL难治复发_VS_ALL缓解")
###t检验函数
test <- function(raw_groupA,raw_groupB,group){
    raw_groupA <- raw_groupA[,-c(1:7)]
    raw_groupB <- raw_groupB[,-c(1:7)]
    T_test_P_value <-c()
    wilcox_test_P_value <-c()
    pnorm_A <- c()
    pnorm_B <- c()
    nrom.test <- c()
    n.col <- ncol(raw_groupA)
    group.list <- rep(group,n.col)
    for (i in 1:n.col){
        temp_A <- raw_groupA[,i][!is.na(raw_groupA[,i])]
        temp_B <- raw_groupB[,i][!is.na(raw_groupB[,i])]
        pnorm_A[i] <- ks.test(temp_A ,"pnorm")$p.value 
        pnorm_B[i] <- ks.test(temp_B ,"pnorm")$p.value
        if (pnorm_A[i] >= 0.05 & pnorm_A[i]>= 0.05) nrom.test[i] <- "yes"
        if (pnorm_A[i] < 0.05  | pnorm_A[i] < 0.05) nrom.test[i] <- "no"
        tt <- var.test(temp_A ,temp_B )$p.value  >= 0.05 #方差齐性检验
        T_test_P_value[i] <- t.test(temp_A , temp_B ,paired = F,var.equal = tt)$p.value
        wilcox_test_P_value[i] <- wilcox.test(raw_groupA[,i],raw_groupB[,i],paired = F)$p.value
    }
    result <- data.frame(Contrast = rep(group,n.col), Site=colnames(raw_groupA), T_test_P_value, wilcox_test_P_value,nrom.test)
    sink(sprintf("%s_T_检验.txt",group))
    print(result)
    sink()

}


