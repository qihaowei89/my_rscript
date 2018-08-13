#读取文件
library(reshape)
setwd("C:/Users/weiqi/Desktop/新建文件夹 (5)/")
source("npmc.R")
do_nonpara_test <- function(data,group){
    data <- data[c(data$组别 %in% group),]
    data.met <- melt(data,id="组别")
    data.met <- data.met[!is.na(data.met$value),]
    site <- levels(data.met$variable)
    result <- NULL
    for (i in site) {
        sub_set  <- subset(data.met,variable==i)[,-2]
        sub_set$value <- as.numeric(sub_set$value)
        colnames(sub_set) <- c("class","var")
        if(length(group) > 2){
            kru_test <- kruskal.test(var~class,data = sub_set)
            if(kru_test$p.value  > 0){
                kru.p  <- kru_test$p.value
                npmc_test <- summary(npmc(sub_set),type = "BF")
                temp <- data.frame(ID = seq(1:3),cbind(kruskal_pValue =  c(kru.p,NA,NA),site = rep(i,3),npmc_test$`Results of the multiple Behrens-Fisher-Test`,note = rep("|",3),npmc_test$`Data-structure`),row.names = 1)
            }
        }else{
            wilcox_test <- wilcox.test(var~class,data = sub_set)
            temp <- data.frame(site=i,wilcox_pValue = wilcox_test$p.value)
        }
        out_put <- rbind(result,temp)
        result <-out_put
    }
    
    filename <- paste(group,collapse = "_")
    file_output <- sprintf("%s差异分析.txt",filename)
    sink(file_output)
    write.table(out_put,quote = F,sep = "\t",row.names = F)
    sink()
}

#file <- read.table("T-B-正常.txt",header = T,encoding = "GB2312",stringsAsFactors = F)
#file <- read.table("高白-低白-正常.txt",header = T,encoding = "GB2312",stringsAsFactors = F)
#file <- read.table("复发-缓解.txt",header = T,encoding = "GB2312",stringsAsFactors = F)
# file <- read.table("病例-对照.txt",header = T,encoding = "GB2312",stringsAsFactors = F)
# file <- read.table("年龄-性别.txt",header = T,encoding = "GB2312",stringsAsFactors = F)
file$积分 <- as.factor(file$积分)
file$总体甲基化水平 <- as.factor(file$总体甲基化水平)
file.bak  <- file

#############"T-B-正常.txt"##################################
# #比较T-ALL与B-ALL与正常对照三组各甲基化位点有无差异，甲基化位点数，积分，平均甲基化率，以及总体甲基化水平有无差异											
# file$组别[file$组别=="正常"]  <- 0
# file$组别[file$组别=="T-ALL"] <- 1
# file$组别[file$组别=="B-ALL高危" | file$组别=="B-ALL中危" | file$组别=="B-ALL低危"] <- 2
# file$组别 <- factor(file$组别,levels = c(0,1,2),labels = c("正常","T-ALL","B-ALL"))
# do_nonpara_test(file,group = c("正常","T-ALL","B-ALL"))

#比较B-ALL高、中、低危三组各甲基化位点有无差异，甲基化位点数，积分，平均甲基化率，以及总体甲基化水平有无差异											
# file <- file.bak
# file$组别[file$组别=="B-ALL低危"]  <- 0
# file$组别[file$组别=="B-ALL中危"]  <- 1 
# file$组别[file$组别=="B-ALL高危"]  <- 2
# file$组别 <- factor(file$组别,levels = c(0,1,2),labels = c("B-ALL低危","B-ALL中危","B-ALL高危"))
# do_nonpara_test(file,group = c("B-ALL低危","B-ALL中危","B-ALL高危"))

#############"高白-低白-正常.txt"####################
#比较高白ALL与低白ALL与正常对照三组各甲基化位点有无差异，甲基化位点数，积分，平均甲基化率，以及总体甲基化水平有无差异
# file$组别[file$组别=="正常"]    <- 0
# file$组别[file$组别=="ALL低白"] <- 1
# file$组别[file$组别=="ALL高白"] <- 2
# file$组别 <- factor(file$组别,levels = c(0,1,2),labels = c("正常","ALL低白","ALL高白"))
# do_nonpara_test(file,group = c("正常","ALL低白","ALL高白"))

#############"复发-缓解.txt"####################
#比较复发难治ALL与缓解ALL两组各甲基化位点有无差异，甲基化位点数，积分，平均甲基化率，以及总体甲基化水平有无差异										
# file$组别[file$组别=="ALL缓解"]    <- 0
# file$组别[file$组别=="B-ALL难治复发"] <- 1
# file$组别 <- factor(file$组别,levels = c(0,1),labels = c("ALL缓解","B-ALL难治复发"))
# do_nonpara_test(file,group = c("ALL缓解","B-ALL难治复发"))


#############"病例-对照.txt"########################
#比较病例和对照组各甲基化位点有无差异，甲基化位点数，积分，平均甲基化率，以及总体甲基化水平有无差异									
# file$组别[file$组别=="正常"]    <- 0
# file$组别[file$组别=="ALL病例"] <- 1
# file$组别 <- factor(file$组别,levels = c(0,1),labels = c("正常","ALL病例"))
# do_nonpara_test(file,group = c("正常","ALL病例"))

#############"年龄-性别.txt"###############
#比较男性ALL与女性ALL两组各甲基化位点有无差异，甲基化位点数，积分，平均甲基化率，以及总体甲基化水平有无差异										
# file$性别分组[file$性别分组=="男"]  <- 0
# file$性别分组[file$性别分组=="女"]  <- 1
# file$性别分组 <- factor(file$性别分组,levels = c(0,1),labels = c("男","女"))
# file <- data.frame(cbind(组别=file$性别分组,file[,c(-1,-2,-3)]))
# do_nonpara_test(file,group = c("男","女"))

#比较男性正常与女性正常两组各甲基化位点有无差异，甲基化位点数，积分，平均甲基化率，以及总体甲基化水平有无差异										
# file <- file.bak
# file$组别[file$性别分组=="男"&file$疾病分组=="正常"]  <- 0
# file$组别[file$性别分组=="女"&file$疾病分组=="正常"]  <- 1
# file$组别<- factor(file$组别,levels = c(0,1),labels = c("男性正常","女性正常"))
# file <- data.frame(cbind(组别=file$组别,file[,c(-1,-2,-3,-28)]))
# do_nonpara_test(file,group = c("男性正常","女性正常"))

#比较8岁以下ALL与8-15岁ALL两组各甲基化位点有无差异，甲基化位点数，积分，平均甲基化率，以及总体甲基化水平有无差异										
# file <- file.bak
# file$组别[file$年龄分组=="＜8" ]  <- 0
# file$组别[file$年龄分组=="8至15"]  <- 1
# file$组别<- factor(file$组别,levels = c(0,1),labels = c("＜8" ,"8至15"))
# file <- data.frame(cbind(组别=file$组别,file[,c(-1,-2,-3,-28)]))
# do_nonpara_test(file,group = c("＜8" ,"8至15"))

#比较8岁以下正常与8-15岁正常两组各甲基化位点有无差异，甲基化位点数，积分，平均甲基化率，以及总体甲基化水平有无差异										
# file <- file.bak
# file$组别[file$年龄分组=="＜8"&file$疾病分组=="正常" ]  <- 0
# file$组别[file$年龄分组=="8至15"&file$疾病分组=="正常"]  <- 1
# file$组别<- factor(file$组别,levels = c(0,1),labels = c("＜8正常" ,"8至15正常"))
# file <- data.frame(cbind(组别=file$组别,file[,c(-1,-2,-3,-28)]))
# do_nonpara_test(file,group = c("＜8正常" ,"8至15正常"))











###########################方差分析####
# file2 <- melt(file[,-c(2,3,4,5)])
# result <- list()
# site <- levels(file2$variable)
# for (i in site) {
#     fit <- aov(value~组别,data=subset(file2,variable==i)[,-2])
#     result_a <- summary(fit)
#     Pr_F  <- unlist(summary(fit))["Pr(>F)1"]
#     result_a1 <- TukeyHSD(fit)
#     a1 <- result_a1$组别[,"p adj"]
#     row_name <- c("方差分析 p_value",paste(rownames(result_a1$组别)," p_adj",""))
#     temp <- cbind(row_name,p_value=c(Pr_F , a1))
#     result <-cbind(result,temp[,2]) 
#     result <- apply(result,2,function(n) as.numeric(n) )
#     row.names(result) <- row_name
# }
# colnames(result) <- colnames(file1[,-1])
# 

