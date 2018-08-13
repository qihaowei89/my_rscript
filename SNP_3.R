#################################
#Data:2018-2-11
#Qihao
##################################

library(stringr)
library(xlsx)
# argv <- commandArgs(TRUE)
# deletNum <- as.numeric(argv[1])
work_dir <- "C:/Users/weiqi/Desktop/SNP_csv/work_dir/SNP"
fileList <- list.files(work_dir,pattern = ".*sum_tr.csv$")
# file <- fileList[1]
countRate <- function(file,i=1){
  inputFile <- read.csv(file,na.strings = c("NA",""),stringsAsFactors = FALSE)
  inputFile <- inputFile[,c(-1:-i)]
  nrow(inputFile)->rowNum
  NAME<-names(inputFile)
  countNa <- c()
  for (i in NAME) {
    which(is.na(inputFile[[i]])) -> naIndex
    inputFile[[i]][naIndex] <- 0
    length(which(inputFile[[i]]==0)) -> countNa[i]
  }
  (rowNum-countNa)/rowNum -> rate 
  rate <- rate[order(rate,decreasing = TRUE)]
  outfile <- paste(str_sub(file,1,-8),"检出率.csv",sep = "_")
  write.csv(rate,outfile)
  
}
# countRate <- function(file,i){
#   inputFile <- read.csv(file,na.strings = c("NA",""),stringsAsFactors = FALSE,encoding = "utf-8")
#   inputFile <- inputFile[,c(-1:-i)]
#   nrow(inputFile)->rowNum
#   NAME<-names(inputFile)
#   countNa <- c()
#   for (i in NAME) {
#     which(is.na(inputFile[[i]])) -> naIndex
#     inputFile[[i]][naIndex] <- 0
#     length(which(inputFile[[i]]==0)) -> countNa[i]
#   }
#   (rowNum-countNa)/rowNum -> rate 
#   rate <- rate[order(rate,decreasing = TRUE)]
#   outfile <- paste(str_sub(file,1,-8),"检出率.csv",sep = "_")
#   write.csv(rate,outfile)
#   
# }
###########
# for (file in fileList){
#   countRate(file,deletNum)
# }
lapply(fileList,function(n) countRate(n))
###########
######################################
# inputFile <- read.csv("WAL0014 _sum_tr.csv",na.strings = c("NA",""),stringsAsFactors = FALSE,encoding = "utf-8")
# 
# inputFile <- inputFile[,c(-1:-3)]
# 
# ncol(inputFile)->colNum
# nrow(inputFile)->rowNum
# 
# NAME<-names(inputFile)
# countNa <- c()
# 
# for (i in NAME) {
#   which(is.na(inputFile[[i]])) -> naIndex
#   inputFile[[i]][naIndex] <- 0
#   length(which(inputFile[[i]]==0)) -> countNa[i]
# }
# 
# (rowNum-countNa)/rowNum -> rate
# rate <- rate[order(rate,decreasing = TRUE)]
# outfile <- paste("PlateData-WAL0040-21-200_tr.csv","检出率.csv",sep = "_")
# write.csv(rate,outfile)






