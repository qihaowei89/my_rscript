#coding:GB2312
#Data 2017-10
#Qihao Wei
###贴转板明细

library(xlsx)
library(stringr)

setwd("C:/Users/weiqi/Desktop/SNP_csv/work_dir/###Methy报告整理###")
file_csv_list <- list.files(pattern = "*.csv$")
file_xlsx_list <- list.files(pattern = "*.xlsx$")
#file_xlsx_list <- list.files(pattern = "*.xls")
# filepath <- file_csv_list[1]
# file_note <- file_xlsx_list[1]
############################################

get_col_from_file <- function(filepath,file_note){
    
    file <- read.csv(filepath,header = TRUE,skip = 1)
    file_colName <- file[,c(1,2,3)]  #文件前三列
    note1 <- read.xlsx(file_note,1,header = F,encoding = "UTF-8")
    note1 <- note1[apply(note1,1,function(x) !all(is.na(x))),apply(note1,2,function(x) !all(is.na(x)))]
    note1 <- note1[order(note1[,1]),]
    t(note1)[2,] -> a
    t(note1)[1,] -> b
    title <- rbind(a,b)
    as.data.frame(title) -> ooo
    
    
    gene_name <- c("基因名称：","SampleID")
    aaa <- c("","CPG Position")
    cbind(gene_name,aaa,ooo) -> uuu
    
    note <- read.xlsx(file_note,1,header = T,encoding = "UTF-8")
    note <-note[apply(note,1,function(x) !all(is.na(x))),apply(note,2,function(x) !all(is.na(x)))]
    note_order <- note[order(note[,1]),]
    get_col_ID <- as.character(note_order[[1]])   ##需要取得列
    get_col    <- subset(file,select = get_col_ID)
    
    get_table <- cbind(file_colName,get_col)
    csv_head <- str_sub(filepath,1,-5)
    result_file_name <- paste("甲基化检测结果_",csv_head,".csv",sep = "")
    sink(result_file_name)
    write.table(uuu,result_file_name,col.names = F,row.names = F,quote = F,sep = ",",append = T)
    write.table(get_table,result_file_name,col.names = T,row.names = F,quote = F,sep = ",",append = T)
    
    sink()
}
#################################################################################################################
#get_col_from_file(filepath ,file_note)

for(i in 1:length(file_csv_list)){
    
    get_col_from_file(file_csv_list[i] ,file_xlsx_list[i])
    
}
