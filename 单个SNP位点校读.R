args <- commandArgs(TRUE)
file_name1 <- list.files(pattern = "*_1.csv") #file_1 已校读文件，包含需替换的位点
file_name2 <- list.files(pattern = "*_2.csv") #file_2 需要替换位点的文件
site <- as.character(args[1])

replace_SNP_site <- function(file_name1=file_name1, file_name2=file_name1, site=site){
    file1 <- read.csv(file_name1, skip = 2, header = T, stringsAsFactors = FALSE)
    file2 <- read.csv(file_name2, skip = 2, header = T, stringsAsFactors = FALSE)
    title_of_file <- read.csv(file_name2, header = FALSE, stringsAsFactors = FALSE)
    temp0 <- title_of_file[,1][1]
    temp1 <- subset(file1, file1$Assay.Id == site) #获得文件1替换位点
    temp2 <- subset(file2, file2$Assay.Id !=site)  #获得文件2中排除需校读位点
    list_snp <- file2["Assay.Id"]
    list <- list_snp[,1]
    output_name <- paste("Replace_",site,".csv",sep = "") #输出文件名字
    ########
    sink(output_name)                               #输出文件到excel表格中
    write.table(temp1,row.names = F,col.names = T,na = "",append = T,sep = ",")
    write.table(temp2,row.names = F,col.names = F,na = "",append = T,sep = ",")
    sink()
    ########
    order_file_name <- list.files(pattern = "Replace_*")
    file <- read.csv(order_file_name,header = T,stringsAsFactors = FALSE)
    order_file <- file[order(file[1]),]             #按384板位置排序
    ########
    sink(output_name)                               #输出文件到excel表格中
    write.table(temp0,row.names = F,col.names = F,na = "",append = T,sep = ",")     #添加第一列
    write.table(order_file,row.names = F,col.names = T,na = "",append = T,sep = ",")
    sink()
    ########
}
replace_SNP_site(file_name1=file_name1, file_name2=file_name1, site=site)
