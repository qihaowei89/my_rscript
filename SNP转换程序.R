########################################################
#Data: 2018-1-9,2018-1-30
#Qihao Wei
#Translate: "AA" to "A A", "NA" to "0 0" ,et.
########################################################

###############函数主体#################
library(xlsx)
library(stringr)
translat2snp <- function(input.file,start.col=6) {
    outputfile <- "SNP格式转换后.csv"
    snp.site <- input.file[,-(1:start.col)]
    colnames(snp.site) <-str_extract(colnames(snp.site),pattern = "rs\\d+")
    n <- apply(snp.site, 2, function(n) as.character(n))
    translat.snp <- function(n){
        for (i in 1:length(n)){
            if ( is.na(n[i]) ) n[i] <- "0 0"
            if ( n[i] == "A" ) n[i] <- "A A"
            if ( n[i] == "T" ) n[i] <- "T T"
            if ( n[i] == "C" ) n[i] <- "C C"
            if ( n[i] == "G" ) n[i] <- "G G" 
            if ( n[i] == "AT"| n[i] == "TA" ) n[i] <- "A T"
            if ( n[i] == "AC"| n[i] == "CA" ) n[i] <- "A C"
            if ( n[i] == "AG"| n[i] == "GA" ) n[i] <- "A G"
            if ( n[i] == "TC"| n[i] == "CT" ) n[i] <- "C T"
            if ( n[i] == "TG"| n[i] == "GT" ) n[i] <- "G T"
            if ( n[i] == "CG"| n[i] == "GC" ) n[i] <- "C G"}
        n
    }
    snp.site_translated <- cbind(input.file[,(1:start.col)], apply(snp.site, 2,function(n) translat.snp(n)))
    write.table(snp.site_translated, outputfile,sep = ",", quote = F, row.names = F)
}









# start.col <- 2 
# #read file from xlsx 
# 
# input.file <- read.xlsx("QTL分析---段晓冉 .xlsx",1, encoding = "UTF-8", header = T)
# outputfile <- "单倍体数据---段晓冉.csv"
# snp.site <- input.file[,-(1:start.col)]
# colnames(snp.site) <-str_extract(colnames(snp.site),pattern = "rs\\d+")
# snp.site <- apply(snp.site, 2, function(n) as.character(n))
# 
# translat.snp <- function(n) {
#     for (i in 1:length(n)){
#         if ( is.na(n[i]) ) n[i] <- "0 0"
#         if ( n[i] == "A" ) n[i] <- "A A"
#         if ( n[i] == "T" ) n[i] <- "T T"
#         if ( n[i] == "C" ) n[i] <- "C C"
#         if ( n[i] == "G" ) n[i] <- "G G" 
#         if ( n[i] == "AT"| n[i] == "TA" ) n[i] <- "A T"
#         if ( n[i] == "AC"| n[i] == "CA" ) n[i] <- "A C"
#         if ( n[i] == "AG"| n[i] == "GA" ) n[i] <- "A G"
#         if ( n[i] == "TC"| n[i] == "CT" ) n[i] <- "C T"
#         if ( n[i] == "TG"| n[i] == "GT" ) n[i] <- "G T"
#         if ( n[i] == "CG"| n[i] == "GC" ) n[i] <- "C G"
#     }
#     n
# }
# 
# 
# snp.site_translated <- cbind(input.file[,(1:start.col)], apply(snp.site, 2,function(n) translat.snp(n)))
# 
# write.table(snp.site_translated, outputfile,sep = ",", quote = F, row.names = F)


# file <- read.csv("单倍体数据---段晓冉.csv")[,c(-1,-2)]
# map  <- read.table("1.map")
# 
# rs1 <- file[,c(1,2,3,4,5,6,7)]
# rs2 <- file[,c(1,2,3,4,5,8,9,10,11,12)]
# rs3 <- file[,c(1,2,3,4,5,13,14,15,16,17,18)]
# rs4 <- file[,c(1,2,3,4,5,19,20,21)]
# 
# rs2map <- function(map,rs,filename){
#      rs.map <- map[map$V2%in%colnames(rs)[7:ncol(rs)], ]
#      map.name <- paste(filename, ".map", sep = "")
#      rs.name  <- paste(filename, ".ped", sep = "")
#      write.table(rs.map, map.name, col.names = F, row.names = F,quote = F)
#      write.table(rs, rs.name,col.names = F, row.names = F,quote = F)
# }
# 
# 
# rs2map(map, rs1, "RNA.TERF2")
# rs2map(map, rs2, "RNA.POT1")
# rs2map(map, rs3, "RNA.TPP1")
# rs2map(map, rs4, "RNA.TERF1")


