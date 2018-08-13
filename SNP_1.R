Sys.setenv(JAVA_HOME='C:/Program Files (x86)/Java/jre1.8.0_144') 
library(stringr)
library(reshape)
library(xlsx)
library(rJava)

path <-"C:/Users/weiqi/Desktop/SNP_csv/work_dir/SNP"

setwd(path)

csv_list <- list.files(path,full.names = FALSE,pattern = "*.csv$")
xlsx_list <- list.files(path,full.names = FALSE,pattern = "*.xlsx$")
#xlsx_list <- list.files(path,full.names = FALSE,pattern = "*.xls")
pro_ID <- str_extract(csv_list[1],"[A-Z]{3}[0-9]{4}")
dir_name <- paste(pro_ID,"_SNP¼ì²â½á¹û",sep="")
dir.create(dir_name)
# dir_path <- paste(path,dir_name,collapse = "/")
#  i=1 
# file_csv <- csv_list[i]
# file_xlsx  <-  xlsx_list[i]
data_re <- function(file_csv,file_xlsx) {
  tamp_csv <- read.table(file_csv,header = TRUE,sep = ",",skip = 2) 
  tamp_xlsx <- read.xlsx(file_xlsx,1,encoding = "UTF-8")
  tamp_xlsx <- tamp_xlsx[apply(tamp_xlsx,1,function(x) !all(is.na(x))),apply(tamp_xlsx,2,function(x) !all(is.na(x)))]
  csv_name <- str_sub(file_csv,1,-5) 
  name_xlsx <- names(tamp_xlsx)[1]
  csv <- cbind(tamp_csv["Sample.Id"],tamp_csv["Assay.Id"],tamp_csv["Call"])
  csv_order <- csv[order(csv[,1],csv[,2]),]
  #file <- str_c(csv_name,".txt",collapse = "") 
  #write.table(data_order,file=file,row.names = FALSE,quote = FALSE,sep = "\t")
  md_csv <- cast(csv_order,Sample.Id~Assay.Id)
  names(md_csv)[1] <- name_xlsx
  table1 <- md_csv
  table2 <- tamp_xlsx[order(tamp_xlsx[,1]),]
  table <- merge(table2,table1,by = name_xlsx)
  file <- str_c(csv_name,"_tr",".csv",collapse = "")
  write.csv(table,file=file ,row.names = FALSE)
}

for (i in 1:length(csv_list)){

  data_re(csv_list[i],xlsx_list[i])

}

# lapply(csv_list,xlsx_list,FUN =function(x,y) data_re(x, y ) )


