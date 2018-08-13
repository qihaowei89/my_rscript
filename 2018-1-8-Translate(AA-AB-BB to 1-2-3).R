#################################################
#qihao Wei
#Data :2018-1-8
#Translate [A,AB,BB] to [1,2,3]
################################################
inputfile <- "需转换.txt"
outputfile <- "已转换.txt"
table <- read.table(inputfile,sep = "\t",stringsAsFactors = F,header = T,row.names = 1)

for(i in colnames(table)) {
    tmp <- table[,i]
    unit <- unique(tmp)
    unit <- unit[order(unit)]
    unit <- unit[unit != ""]
    
    if (length(unit)==3){
        for (j in 1:length(row.names(table))){
            if ( tmp[j] == unit[3]) {tmp[j] <- 3} 
            if ( tmp[j] == unit[2]) {tmp[j] <- 2}
            if ( tmp[j] == unit[1]) {tmp[j] <- 1}
            #if ( tmp[j] == "") {tmp[j] <- 0}
            }
    tmp -> table[,i]
    }
    if (length(unit)==2){
        for (j in 1:length(row.names(table))){
            if ( tmp[j] == unit[2]) {tmp[j] <- 2}
            if ( tmp[j] == unit[1]) {tmp[j] <- 1}
            #if ( tmp[j] == "") {tmp[j] <- 0}
            }
    tmp -> table[,i]
    }
    if (length(unit)==1){
        for (j in 1:length(row.names(table))){
            if (tmp[j] == unit[1]) {tmp[j] <- 1}
            #if ( tmp[j] == "") {tmp[j] <- 0}
            }
    tmp -> table[,i]
    }
}
table <- data.frame(ID = rownames(table),apply(table,2,function(n) as.factor(n)))
#table <- data.frame(apply(table,2,function(n) as.factor(n)))
write.table(table,outputfile,row.names = T,quote = F,sep = ",",na = "")


