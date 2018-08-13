library(Hmisc)
library(stringr)
list.files(pattern = "*.sav") -> filelist

sav_to_csv <- function(i) {
    name <- str_sub(i,1,-5)
    name <- paste(name,".csv",sep = "")
    
    file <- spss.get(i,use.value.labels = TRUE)
    write.csv(file,name,row.names = FALSE,na = "")
}

for (i in filelist){
    sav_to_csv(i)
}
