library(xlsx)
library(stringr)
source("Script/2018-2-5-SNP转换程序.R")

snp <- read.xlsx("QTL分析---段晓冉.xlsx",1,header = T)
snp1 <- snp[snp$group==1,]
snp2 <- snp[snp$group==2,]

aa <- read.xlsx("QTL分析---段晓冉.xlsx",2,header = T,encoding="UTF-8")

plot(aa$ex,aa$age)
aa <- read.table("g1.covar",header = T)
aa <- read.table("g2.covar",header = T)
corr <- corr.test(aa[3:8],method = "spearman")
sink("g2_总暴露剂量与年龄相关性.txt")
print(corr$r)
print(corr$p)
sink()

aa1 <- aa[aa$分组==1,][,-1]
aa2 <- aa[aa$分组==2,][,-1]
write.table(aa1,"g1.covar",quote = F, row.names = F,sep = "\t")
write.table(aa2,"g2.covar",quote = F, row.names = F,sep = "\t")

bb <- read.xlsx("QTL分析---段晓冉.xlsx",3,header = T,encoding="UTF-8")
bb1 <- bb[bb$分组==1,][,-1]
bb2 <- bb[bb$分组==2,][,-1]
write.table(bb1,"g1.pheno",quote = F, row.names = F,sep = "\t")
write.table(bb2,"g2.pheno",quote = F, row.names = F,sep = "\t")


map <- read.csv("result.csv",header = F,stringsAsFactors = F)

# snp$FID  <- as.character(snp$FID)
# snp$group[grepl("(TP.)|(DQ.)",snp$FID)] <- 1
# snp$group[!grepl("(TP.)|(DQ.)",snp$FID)] <- 2

snp2ped(snp1,outputfile = "g1")
snp2ped(snp2,outputfile = "g2")



result2map(map)



