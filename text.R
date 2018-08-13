#Date: 2017-10-13
#Qihao WEI
##########################################################

#library(xlsx)
#library(xlsxjars)

argv <- commandArgs(TRUE)

file1_name <- argv[1]
file2_name <- argv[2]
file1_name <- "vitB12_analysis_diff.txt"
file2_name <- "vitB6_analysis_diff.txt"
name1 <- strsplit(file1_name,split = "_",fixed = FALSE)[[1]][1]
name2 <- strsplit(file2_name,split = "_",fixed = FALSE)[[1]][1]

output_name <- paste(name1,name2,sep = "_VS_")
output_name <- paste(output_name,".csv",sep = "")

file1 <- read.csv(file1_name,header = T,sep = "\t")
file2 <- read.csv(file2_name,header = T,sep = "\t")

site1 <- row.names(file1)
site1 <- site1[1:20]
site <- row.names(file2)

file3 <- cbind(site,file2)


sink(output_name)

for (i in site1){
  if(i %in% site){
    
    temp <- subset(file3,site2==i)
    if(file.info(output_name)$size == 0){
      write.table(temp,output_name,col.names = T,row.names = F,quote = F,sep = ",",append = T)
    }else{
      write.table(temp,output_name,col.names = F,row.names = F,quote = F,sep = ",",append = T)
    }
    sink()
  }

}



