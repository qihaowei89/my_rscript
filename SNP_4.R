######################################
#Data:2018-2-11 
#

Sys.setenv(JAVA_HOME='C:/Program Files (x86)/Java/jre1.8.0_144') 

library(stringr)
library(reshape)
library(xlsx)
library(rJava)

path <-"C:/Users/weiqi/Desktop/SNP_csv/work_dir/SNP/"
setwd(path)

csv_list <- list.files(path,full.names = FALSE,pattern = "*.csv")

# file_csv <- csv_list[1]

data_re <- function(file_csv) {
  tamp_csv <- read.table(file_csv,header = TRUE,check.names = T,sep = ",",skip = 2) 
  csv_name <- str_sub(file_csv,1,-5) 
  csv <- cbind(tamp_csv["Sample.Id"],tamp_csv["Assay.Id"],tamp_csv["Call"])
  csv_order <- csv[order(csv[,1],csv[,2]),]
  #file <- str_c(csv_name,".txt",collapse = "") 
  #write.table(data_order,file=file,row.names = FALSE,quote = FALSE,sep = "\t")
  md_csv <- cast(csv_order,Sample.Id~Assay.Id)
  #names(md_csv)[1] <- name_xlsx
  table1 <- md_csv
  file <- str_c(csv_name,"_tr",".csv",collapse = "")
  write.csv(table1,file=file ,row.names = FALSE)
}



lapply(csv_list, function(n) data_re(n))


