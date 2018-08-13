#Date: 2017-10-13
#Qihao WEI
##########################################################
#
#library(xlsx)
#library(xlsxjars)


file1_name <- "PlateData-BML0025-W2-B-13-264_tr.csv"
file2_name <- "PlateData-BML0025-W2_tr.csv"
take_ID <- "样品ID"
file1 <- read.csv(file1_name,header = T,sep = ",")
file2 <- read.csv(file2_name,header = T,sep = ",")

site1 <- file1[,4] #目标
site2 <- file2[,4] #源文件

#############################################
######   提取需要更换的列  ##################
sink("result.csv")

for (i in site1){
  if(i %in% site2){
    temp <- subset(file2, 样品ID==i)
    if(file.info("result.csv")$size == 0){
      write.table(temp,"result.csv",col.names = T,row.names = F,quote = F,sep = ",",append = T)
    }else{
      write.table(temp,"result.csv",col.names = F,row.names = F,quote = F,sep = ",",append = T)
    }
    sink()
  }
  
}

#############################################
######   提取需要除更换的列以外  ############
sink("result_EX.csv")

for (i in site2){
    if(i %in% site1){
    }else{
        temp <- subset(file2, 样品ID== i)
        if(file.info("result_EX.csv")$size == 0){
            write.table(temp,"result_EX.csv",col.names = T,row.names = F,quote = F,sep = ",",append = T)
        }else{
            write.table(temp,"result_EX.csv",col.names = F,row.names = F,quote = F,sep = ",",append = T)
        }
        sink()
    }
}


