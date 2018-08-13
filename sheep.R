###########################################
#Data: 20187-1-9
#Qihao Wei
#compare sheep and other animal by SNPs
###########################################
library(xlsx)


table <- read.xlsx("data.xlsx",1,header = T)

table.sheep <- table[table$CLASS == "sheep",][,c(-1,-2)]
apply(table.sheep,MARGIN = 2,function(n) as.character(n)) -> table.sheep
apply(table.sheep,MARGIN = 2, function(n) unique(n)) -> sheep

table.other <- table[table$CLASS == "others",][,c(-1,-2)]
apply(table.other,MARGIN = 2,function(n) as.character(n)) -> table.other
apply(table.other,MARGIN = 2, function(n) unique(n)) -> other

diff.site <- list()
for(i in colnames(table.sheep)){
    intersect(sheep[i][[1]],other[i][[1]]) -> inter
    inter <- inter[which(inter != "NA")]
    if (length(inter) == 0) diff.site[i] <- 1
    if (length(inter) == 1) diff.site[i] <- 1-length(which(table.other[,i] != "NA" & table.other[,i] == inter))/length(which(table.other[,i] != "NA"))
    if (length(inter) == 2) diff.site[i] <- 0
    if (length(inter) == 3) diff.site[i] <- 0
    if (length(inter) == 4) diff.site[i] <- 0
}

write.table(diff.site,"diff_site.csv",sep = ",",row.names = F)
