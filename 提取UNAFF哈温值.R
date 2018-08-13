file <- read.table("¹þÎÂÆ½ºâ.hwe",header = T)

file.unaff <- file[file$TEST == "UNAFF",]

write.table(file.unaff,"¹þÎÂÆ½ºâ_UNAFF.hwe",sep = "\t",quote = F, row.names = F,col.names = T)
