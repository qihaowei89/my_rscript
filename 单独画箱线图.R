library(reshape)
library(ggplot2)
setwd("C:/Users/weiqi/Desktop/SNP_csv/")
sah <- read.table("sah.txt",header = T,stringsAsFactors = F,sep = "\t")
sam <- read.table("sam.txt",header = T,stringsAsFactors = F,sep = "\t")

a <- sah[,-1]
a <- melt(a)
result_a <- aov(value~variable,data=a)
summary(result_a) -> sum_sah
#result_b1 <- LSD.test(result_a,"variable",p.adj="bonferroni")
result_a1 <- TukeyHSD(result_a)
b <- sam[,-1]
b <- melt(b)
result_b <- aov(value~variable,data=b)
summary(result_b) ->sum_sam
#result_b1 <- LSD.test(result_b,"variable",p.adj="bonferroni") 
result_b1 <- TukeyHSD(result_b)
############################################################
SAH_T_test_P_value <- c()
SAH_T_test_P_value[1] <- t.test(sah[,2][!is.na(sah[,2])], sah[,3][!is.na(sah[,3])],paired = F)$p.value
SAH_T_test_P_value[2] <- t.test(sah[,2][!is.na(sah[,2])], sah[,4][!is.na(sah[,4])],paired = F)$p.value
SAH_T_test_P_value[3] <- t.test(sah[,3][!is.na(sah[,3])], sah[,4][!is.na(sah[,4])],paired = F)$p.value
t_test.sah <- data.frame(T_test = c("AC-AH","AC-AHL","AH-AHL"),p_value = SAH_T_test_P_value)
###

SAM_T_test_P_value <- c()
SAM_T_test_P_value[1] <- t.test(sam[,2][!is.na(sam[,2])], sam[,3][!is.na(sam[,3])],paired = F)$p.value
SAM_T_test_P_value[2] <- t.test(sam[,2][!is.na(sam[,2])], sam[,4][!is.na(sam[,4])],paired = F)$p.value
SAM_T_test_P_value[3] <- t.test(sam[,3][!is.na(sam[,3])], sam[,4][!is.na(sam[,4])],paired = F)$p.value
t_test.sam <- data.frame(T_test = c("AC-AH","AC-AHL","AH-AHL"),p_value = SAM_T_test_P_value)
######################################################################################

#作图  
rbind(sah,sam) -> sss
ss <- melt(sss,id="group")
names(ss) <- c("metabolite","group","value")
ss$group <- factor(ss$group)

pdf("boxplot.pdf")
p <- ggplot(ss, aes(x=group, y=value,fill=metabolite)) + geom_boxplot()+ geom_dotplot(binaxis='y', stackdir='center',position=position_dodge(1),dotsize=0.35)+
    theme(axis.text.x = element_text(angle=40, hjust=1, vjust=1))
p
dev.off()

pdf("boxplot2.pdf")
q <- p + facet_wrap(~metabolite,scales = "free")
q
dev.off()

###############
# output file #
sink("方差分析.txt")
print(sum_sah)
print(result_a1)
print(sum_sam)
print(result_b1)
sink()

sink("T_检验.txt")
print(t_test.sah)
print(t_test.sam)
sink()