########################################################
#Data: 2018-1-9,2018-1-30，2018-2-5
#Qihao Wei
#Translate: "AA" to "A A", "NA" to "0 0" ,et.
########################################################
setwd("C:/Users/weiqi/Desktop/SNP_csv/work_dir/")
###############函数主体#################
library(xlsx)
library(stringr)
snp2ped <- function(input.file,start.col=6, outputfile=1, .col.names=F) {
    outputfile <- paste(outputfile,".ped",sep = "")
    snp.site <- input.file[,-(1:start.col)]
    colnames(snp.site) <-str_extract(colnames(snp.site),pattern = "rs\\d+")
    snp.site <- apply(snp.site, 2, function(n) as.character(n))
    snp.site <- snp.site[,colnames(snp.site)[order(colnames(snp.site))]]
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
    snp.site_translated <- cbind(input.file[,(1:6)], apply(snp.site, 2,function(n) translat.snp(n)))
    write.table(snp.site_translated, outputfile,sep = "\t", quote = F, row.names = F,col.names = .col.names)
}
result2map <- function(file){
    chr <- list()
    pos <- list()
    file <- file[order(file[,1]),] #reordered SNP sites 
    chr_pos <- file[,2]
    snp <- file[,1]
    for (i in 1:length(chr_pos)){
        temp <- strsplit(chr_pos[i],":")[[1]]
        chr[i] <- temp[1]
        pos[i] <- temp[2]
    }
    results <- cbind(unlist(chr),snp,GD=rep(0,length(chr_pos)),unlist(pos))
    colnames(results) <- c("chr","snp","GD","pos")
    results <- as.data.frame(results)
    results <- apply(results,2,function(n) as.character(n))
    write.table(results, "1.map", sep = "\t",col.names = F,row.names = F,quote = F)
}
#result2covar <- function(file){} 

# res <- read.csv(file = "work_dir/result.csv",header = F,stringsAsFactors = F)
# result2map(res) 

snp2ped(input.file1,start.col = 7,outputfile = "group1")
snp2ped(input.file2,start.col = 7,outputfile = "group2")
snp2ped(input.file3,start.col = 7,outputfile = "group3")

#######
input.file1 <-  read.xlsx("Group1.xlsx",1)
input.file2 <-  read.xlsx("Group2.xlsx",1)
input.file3 <-  read.xlsx("Group3.xlsx",1)

covar1 <- input.file1[,c(1,2,5,7)]
write.table(covar1, "covar1",sep = "\t", quote = F, row.names = F)
covar2 <- input.file2[,c(1,2,5,7)]
write.table(covar2, "covar2",sep = "\t", quote = F, row.names = F)
covar3 <- input.file3[,c(1,2,5,7)]
write.table(covar3, "covar3",sep = "\t", quote = F, row.names = F)


#############################
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


