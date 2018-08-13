##2017-11-26
##Qihao

library(stringr)

file_name1 <- list.files(pattern = "*_1.csv") #file_1 已校读文件，包含需替换的位点
#file_name2 <- list.files(pattern = "*_1.csv") #file_1 已校读文件，包含需替换的位点
file_name2 <- list.files(pattern = "*_2.csv") #file_2 需要替换位点的文件

#site <- c("rs1805097","rs1800624","rs2070424")
site <- read.table("SNP.txt",stringsAsFactors = F,sep = "\t",header = F)

replace_SNP_site <- function(file_name1, file_name2, site){
    file1 <- read.csv(file_name1, skip = 2, header = T, stringsAsFactors = FALSE)
    file1_1 <- read.csv(file_name1, skip = 2, header = F, stringsAsFactors = FALSE)
    temp1_0 <- file1_1[1,]
    file2 <- read.csv(file_name2, skip = 2, header = T, stringsAsFactors = FALSE)
    title_of_file <- read.csv(file_name2, header = FALSE, stringsAsFactors = FALSE)
    temp0 <- title_of_file[,1][1]                   #获得第一行
    
    list_snp <- file2["Assay.Id"]
    list <- list_snp[,1]
    as.factor(list) -> ff  #转换成因子，获得唯一值
    list <- levels(ff)
    #output_name <- paste("Replace_",site,".csv",sep = "") #输出文件名字
    
    output_name1 <- "temp_file.csv"
    sink(output_name1)  #输出文件到excel表格中
    write.table(temp1_0,row.names = F,col.names = F,na = "",append = T,sep = ",")
    for (i in list){
        if (i %in% site) {
            temp1 <- subset( file1, file1$Assay.Id == i )
            write.table(temp1,row.names = F,col.names = F,na = "",append = T,sep = ",")
        }else{
            temp2 <- subset( file2, file2$Assay.Id == i )  
            write.table(temp2,row.names = F,col.names = F,na = "",append = T,sep = ",")
        }
    }
    sink()

    order_file_name <- list.files(pattern = "temp_file*")
    file <- read.csv(order_file_name,header = T,stringsAsFactors = FALSE)
    order_file <- file[order(file$Sample.Id,file$Assay.Id),]             #排序
    
    file2_name <- paste(str_sub(file_name2,1,-7),"_re.csv",sep = "") 
    sink(file2_name)                               #输出文件到excel表格中
    write.table(temp0,row.names = F,col.names = F,na = "",append = T,sep = ",")     #添加第一列
    write.table(order_file,row.names = F,col.names = T,na = "",append = T,sep = ",")
    sink()
}

replace_SNP_site(file_name1, file_name2, site)

