file <- read.table("����ƽ��.hwe",header = T)

file.unaff <- file[file$TEST == "UNAFF",]

write.table(file.unaff,"����ƽ��_UNAFF.hwe",sep = "\t",quote = F, row.names = F,col.names = T)
