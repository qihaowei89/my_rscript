library(xlsx)
library(readxl)
#file <- read.xlsx("PlateData-XJJ0048-W1+W2-384_tr.xlsx",2,encoding = "UTF-8",colClasses = "character")

# file1 <- read_excel("PlateData-XJJ0048-W1+W2(汇总)tr--徐佳佳用(bak).xlsx",2)
file <- read.csv("XJJ0048_file.csv",stringsAsFactors = F,na.strings = "")
file1 <- file[,-64]

sheep <- file1[file1$species=="羊",][,-c(2,3)]
sheep <- as.data.frame(sheep)

sheep <- apply(X = sheep,MARGIN = 2, FUN = function(n) as.factor(n))

#other <- file1[!file1$species=="羊",][,-c(2,3)]
#other <- apply(X = other,MARGIN = 2, FUN = function(n) as.factor(n))
file$species

chicken <- apply(file1[file1$species=="鸡",][,-c(2,3)],MARGIN = 2, FUN = function(n) as.factor(n))
horse <- apply(file1[file1$species=="马",][,-c(2,3)] ,MARGIN = 2, FUN = function(n) as.factor(n)) 
cow <- apply(file1[file1$species=="牛",][,-c(2,3)],MARGIN = 2, FUN = function(n) as.factor(n))
pig<- apply(file1[file1$species=="猪",][,-c(2,3)],MARGIN = 2, FUN = function(n) as.factor(n))

#summary(other) -> S

sink("统计某个物种中检测不出来的位点.csv")
write.csv(summary(chicken))
write.csv(summary(horse ))
write.csv(summary(cow))
write.csv(summary(pig))
sink()

sink("bbb.csv")
write.csv(S)
sink()

