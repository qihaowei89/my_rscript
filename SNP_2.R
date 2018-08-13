library(stringr)
setwd("C:/Users/weiqi/Desktop/SNP_csv/work_dir/SNP")

#接收手动输入参数
argv <- commandArgs(TRUE)
#argv <- "BML0081"
#group_1 <- "XJJ0022-W1"


#获得文件列表
csv_list <- list.files(getwd(),full.names = FALSE,pattern = "_tr.csv")
 
#根据输入参数的名字创建汇总文件
      
grouplist <- as.character(argv[1])
file_name <- paste(grouplist,"_sum_tr.csv",sep = "")


for (i in 1:length(csv_list)) { 
    if (i == 1){
        file <- read.csv(csv_list[i],header = FALSE,sep = ",")
    }else{
        file_temp <- read.csv(csv_list[i],header = FALSE,sep = ",")
        file <- file_temp[-1,]
    }
    write.table(file,file = file_name,append = TRUE,quote = FALSE,row.names = FALSE,col.names = FALSE,sep = ",")
}
# }else{for (i in 1:length(csv_list)) { 
#         file <- read.table(csv_list[i],header = FALSE,sep = ",")
#         if (i == 1){
#             file <- read.csv(csv_list[i],header = FALSE,sep = ",")
#         }else{
#             file_temp <- read.csv(csv_list[i],header = FALSE,sep = ",")
#             file <- file_temp[-1,]
#         }
#         write.table(file,file = file_namelist,append = TRUE,quote = FALSE,row.names = FALSE,col.names = FALSE,sep = ",")
#         group <- str_extract(csv_list[i],"[A-Z]{3}[0-9]{4}-[A-Z]+[0-9]+")
#       for (i in 1:length(argv)) {
#           if (group == grouplist[i]) {
#              write.table(file,file = file_namelist[i],append = TRUE,quote = FALSE,row.names = FALSE,col.names = FALSE,sep = ",")
#           }
#         }
#       }
#  }



