#��ȡ�ļ�
library(reshape)
setwd("C:/Users/weiqi/Desktop/�½��ļ��� (5)/")
#file <- read.table("����-�Ա�.txt",header = T,encoding = "GB2312",stringsAsFactors = F)
#file <- read.table("����-����.txt",header = T,encoding = "GB2312",stringsAsFactors = F)
#file <- read.table("�߰�-�Ͱ�-����.txt",header = T,encoding = "GB2312",stringsAsFactors = F)
file <- read.table("T-B-����.txt",header = T,encoding = "GB2312",stringsAsFactors = F)
file$���� <- as.factor(file$����)
file$����׻���ˮƽ <- as.factor(file$����׻���ˮƽ)

# file$���[file$���=="����"] <- 0
# file$���[file$���=="ALL�Ͱ�"] <- 1
# file$���[file$���=="ALL�߰�"] <- 2
# file$��� <- factor(file$���,levels = c(0,1,2),labels = c("����","ALL�Ͱ�","ALL�߰�"))
file$���[file$���=="����"] <- 0
file$���[file$���=="T-ALL"] <- 1
file$���[file$���=="B-ALL��Σ" | file$���=="B-ALL��Σ" | file$���=="B-ALL��Σ"] <- 2
file$��� <- factor(file$���,levels = c(0,1,2),labels = c("����","T-ALL","B-ALL"))

file_chi <- file[,c(1,3,5)]
table1 <- with(file_chi,table(���,����))
table2 <- with(file_chi,table(���,����׻���ˮƽ))
# file$���[file$���=="B-ALL��Σ"] <- 0
# file$���[file$���=="B-ALL��Σ"] <- 1
# file$���[file$���=="B-ALL��Σ"] <- 2
# file$��� <- factor(file$���,levels = c(0,1,2),labels = c("B-ALL��Σ","B-ALL��Σ","B-ALL��Σ"))
# file <- file[which(file$���=="B-ALL��Σ"|file$���=="B-ALL��Σ"|file$���=="B-ALL��Σ"),]
file1 <- file[,-c(2,3,4,5)]
file2 <- melt(file1)
result <- list()
site <- levels(file2$variable)
for (i in site) {
    fit <- aov(value~���,data=subset(file2,variable==i)[,-2])
    result_a <- summary(fit)
    Pr_F  <- unlist(summary(fit))["Pr(>F)1"]
    result_a1 <- TukeyHSD(fit)
    a1 <- result_a1$���[,"p adj"]
    row_name <- c("������� p_value",paste(rownames(result_a1$���)," p_adj",""))
    temp <- cbind(row_name,p_value=c(Pr_F , a1))
    result <-cbind(result,temp[,2]) 
    result <- apply(result,2,function(n) as.numeric(n) )
    row.names(result) <- row_name
}
colnames(result) <- colnames(file1[,-1])



# sink("�߰�_vs_�Ͱ�_vs_����_�������.txt")
# print(result)
# sink()

sink("B-ALL_��_��_��_�������.txt")
print(result)
sink()

# file$����[file$����=="ALL����"] <- 0
# file$����[file$����=="B-ALL���θ���"] <- 1
# file$���� <- factor(file$����,levels = c(0,1),labels = c("ALL����","B-ALL���θ���")) 
# # file$�Ա����[file$�Ա����=="��"] <- 0
# file$�Ա����[file$�Ա����=="Ů"] <- 1
# file$�Ա���� <- factor(file$�Ա����,levels = c(0,1),labels = c("��","Ů"))
# file$�������[file$�������=="��8"] <- 1
# file$�������[file$�������=="8��15"] <- 2 
# file$������� <- factor(file$�������,levels = c(1,2),labels = c("��8","8��15"))   
# file$����[file$����=="����"] <- 0
# #file$��������[file$��������=="����"] <- 1 
# file$����[file$����=="ALL����"] <- 1 
# file$���� <- factor(file$����,levels = c(0,1),labels = c("����","ALL����"))  
##
B  <- file[file$����=="B-ALL���θ���",]
re <- file[file$����=="ALL����",]
##
# case <- file[file$���� == "ALL����",]
# control <- file[file$���� == "����",]
##
# all_male       <- file[file$�Ա����== "��",]
# all_female     <- file[file$�Ա����== "Ů",]
# all_male.control <- file[file$�Ա����== "��"&file$��������=="����",]
# all_female.control <- file[file$�Ա����== "Ů"&file$��������=="����",]
# age_less_8 <- file[file$������� == "��8",]
# age_big_8 <- file[file$������� == "8��15",]
# age_less_8.control <- file[file$������� == "��8"&file$��������=="����",]
# age_big_8.control <- file[file$������� == "8��15"&file$��������=="����",]

###���׻���λ�����

# male_and_female.diff <- test(all_male,all_female,"��ALL_vs_ŮALL")                                  #������Ϣ����ALL   ŮALL
# control_male_and_female.diff <- test(all_male.control,all_female.control,"������_vs_Ů����")        #������Ϣ��������  Ů����
# age_8.diff <- test(age_less_8, age_big_8,"��8_vs_8��15")                                            #������Ϣ����8     8��15
# control_age_8.diff <- test(age_less_8.control, age_big_8.control,"��8����_vs_8��15����")            #������Ϣ����8���� 8��15����
# case_and_control <- test(case,control,"ALL����_VS_����")
B_and_re <- test(B,re,"B-ALL���θ���_VS_ALL����")
###t���麯��
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
        tt <- var.test(temp_A ,temp_B )$p.value  >= 0.05 #�������Լ���
        T_test_P_value[i] <- t.test(temp_A , temp_B ,paired = F,var.equal = tt)$p.value
        wilcox_test_P_value[i] <- wilcox.test(raw_groupA[,i],raw_groupB[,i],paired = F)$p.value
    }
    result <- data.frame(Contrast = rep(group,n.col), Site=colnames(raw_groupA), T_test_P_value, wilcox_test_P_value,nrom.test)
    sink(sprintf("%s_T_����.txt",group))
    print(result)
    sink()

}

