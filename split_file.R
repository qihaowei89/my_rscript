###############################################
#Data: 2017-10
#Qihao Wei
#将甲基化原始文件分割成每个 订单#方案 一个文件
#读取原始文件，并生成各个单号的散文件

filepath <- list.files(pattern = "*.csv$")
file <- read.csv(filepath,skip = 1,header = TRUE)

ID_GROUP <- as.character(file[[1]])
group = NULL
############################################################################
##################
for (i in 1:length(ID_GROUP)){
  temp <- ID_GROUP[i]
  group[i] <- strsplit(temp,"_")[[1]][1]
}
group_list <- levels(factor(group)) #获得唯一值
########################################################################
cbind(group,file) -> f #添加分组信息

for (i in group_list){
  save_name <- paste(i,".csv",sep = "")
  temp <- subset(f,group==i)[,-1]
  write.csv(temp,save_name,row.names = F)
}

#####################################################################





