####
file1_name <- "PlateData-BML0046-BT3-8-180_tr.csv"
file2_name <- "PlateData-BML0046 _sum_tr.csv"


file1 <- read.csv(file1_name,header = T,sep = ",",stringsAsFactors = FALSE)
file2 <- read.csv(file2_name,header = T,sep = ",",stringsAsFactors = FALSE)

site1 <- as.character(file1[,2])
site <- file2[,2]


sink("result.csv")

for (i in site){
  if(i %in% site1){
  }else{
    temp <- subset(file2, X96位置== i)
    if(file.info("result.csv")$size == 0){
      write.table(temp,"result.csv",col.names = T,row.names = F,quote = F,sep = ",",append = T)
    }else{
      write.table(temp,"result.csv",col.names = F,row.names = F,quote = F,sep = ",",append = T)
    }
    sink()
  }
}