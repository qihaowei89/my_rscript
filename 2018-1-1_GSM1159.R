# Version info: R 3.2.3, Biobase 2.30.0, GEOquery 2.40.0, limma 3.26.8
# R scripts generated  Sat Jan 13 01:11:44 EST 2018

################################################################
diffexper <- function(gset, groupA,groupB,output_name){
    gset <- gset[ ,c(groupA,groupB)]
    # log2 transform
    ex <- exprs(gset)
    qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
    LogC <- (qx[5] > 100) ||
        (qx[6]-qx[1] > 50 && qx[2] > 0) ||
        (qx[2] > 0 && qx[2] < 1 && qx[4] > 1 && qx[4] < 2)
    if (LogC) { ex[which(ex <= 0)] <- NaN
    exprs(gset) <- log2(ex) }
    #log2.ex <- log2(ex)
    # set up the data and proceed with analysis
    sml <- c(rep("1",length(groupA)), rep("0",length(groupB)))
    sml <- paste("G", sml, sep="")    # set group names
    fl <- as.factor(sml)
    gset$description <- fl
    design <- model.matrix(~ description + 0, gset)
    colnames(design) <- levels(fl)
    fit <- lmFit(gset, design)
    cont.matrix <- makeContrasts(G1-G0, levels=design)
    fit2 <- contrasts.fit(fit, cont.matrix)
    fit2 <- eBayes(fit2, 0.01)
    #tT <- topTable(fit2, adjust="fdr", sort.by="B",number = Inf)
    tT <- topTable(fit2, adjust="BH", sort.by="P",number = Inf)
    #tT <- subset(tT, select=c("ID","Gene.title","Gene.symbol","logFC","P.Value"))
    #diff.tT <- subset(tT,P.Value < 0.05 & (logFC >1 | logFC < -1), select=c("ID","Gene.title","Gene.symbol","logFC","P.Value"))
    diff.ID <- which(tT$P.Value < 0.05 & (tT$logFC >1 | tT$logFC < -1))
    ex.ID <- rownames(ex)
    ex <- cbind(ID=ex.ID,ex)
    diff.tT <- tT[diff.ID,]
    #diff.ex <- ex[diff.ID,]
    diff.output <- merge(diff.tT,ex)
    file.name <- sprintf("%s_diff.csv",output_name)
    write.table(diff.output, file=file.name, row.names=F, sep=",")
    
    # # set parameters and draw the plot
    # labels <- c("NB","Fav")
    # palette(c("#c7ff9d","#f4dff4", "#AABBCC"))
    # dev.new(width=4+dim(gset)[[2]]/5, height=6)
    # par(mar=c(2+round(max(nchar(sampleNames(gset)))/2),4,2,1))
    # title <- paste ("GSE1159", '/', annotation(gset), " selected samples", sep ='')
    # boxplot(ex, boxwex=0.6, notch=T, main=title, outline=FALSE, las=2, col=fl)
    # legend("topleft", labels, fill=palette(), bty="n")
    
    
    ###volcano plot
    dif.result <- subset(tT,select=c("logFC","P.Value"))
    dif.result$log10_P.Value <- -(log10(dif.result$P.Value))
    col.vector <- rep(rgb(0,0,1,0.2),length(row.names(dif.result)))
    col.vector[dif.result$log10_P.Value >=3 &(dif.result$logFC >3 | dif.result$logFC < -3) ] <- rgb(1,0,0)
    
    # plot(x =  dif.result$logFC, y= dif.result$log10_P.Value,
    #      col = col.vector ,
    #      abline(h = 3,v = c(-3,3),lwd=2,lty=3,col= "#4C5B61"),
    #      xlab = "logFC",
    #      ylab = "-log10_P.Value")
    #      
}

#   Differential expression analysis with limma
library(Biobase)
library(GEOquery)
library(limma)

# load series and platform data from GEO

gset <- getGEO("GSE1159", GSEMatrix =TRUE, AnnotGPL=TRUE)
if (length(gset) > 1) idx <- grep("GPL96", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]

# make proper column names to match toptable 
fvarLabels(gset) <- make.names(fvarLabels(gset))

#readfile

#groupfile <- read.table("group.txt",header = T,stringsAsFactors = F)
# group names for all samples
sample_name <- sampleNames(gset)

Fav <- paste("GSM", c("20857","20727","20745","20852","20860","20866","20867"),sep = "")
Adv <- paste("GSM", c("20791","20795","20861","20897"),sep = "")
NBM <- paste("GSM", c("20970","20971","20972","20973","20974"),sep = "")
# Fav <- groupfile$SampleID[1:5]
# Rel <- groupfile$SampleID[6:9]
# Ref <- groupfile$SampleID[10:30]
# NBM <- groupfile$SampleID[31:35]

Fav_list <- which(sample_name %in% Fav)
Adv_list <- which(sample_name %in% Adv)
NBM_list <- which(sample_name %in% NBM)
# Fav_list <- which(sample_name %in% Fav)
# Rel_list <- which(sample_name %in% Rel)
# Ref_list <- which(sample_name %in% Ref)
# NBM_list <- which(sample_name %in% NBM)
# 
###分组1
diffexper(gset, Adv_list,Fav_list,"Adv_Fav")
diffexper(gset, Adv_list,NBM_list,"Adv_NBM")
diffexper(gset, Fav_list,NBM_list,"Fav_NBM")

###分组2
# diffexper(gset,Rel_list,Fav_list,"Rel_Fav")
# diffexper(gset,Ref_list,Fav_list,"Ref_Fav")
# diffexper(gset,Fav_list,NBM_list,"Fav_NBM")
# diffexper(gset,Rel_list,NBM_list,"Rel_NBM")
# diffexper(gset,Ref_list,NBM_list,"Ref_NBM")
# groupA <- NBM_list
# groupB <- Fav_list
# output_name <- "Rel_Fav"


################################################################
# #   Boxplot for selected GEO samples
# library(Biobase)
# library(GEOquery)
# 
# # load series and platform data from GEO
# 
# gset <- getGEO("GSE1159", GSEMatrix =TRUE, getGPL=FALSE)
# if (length(gset) > 1) idx <- grep("GPL96", attr(gset, "names")) else idx <- 1
# gset <- gset[[idx]]
# 
# # # group names for all samples in a series
# # gsms <- paste0("XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX",
# #                "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX",
# #                "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX",
# #                "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX",
# #                "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX",
# #                "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX11100000")
# # sml <- c()
# # for (i in 1:nchar(gsms)) { sml[i] <- substr(gsms,i,i) }
# # sml <- paste("G", sml, sep="")  #set group names
# 
# # eliminate samples marked as "X"
# sel <- which(sml != "X")
# sml <- sml[sel]
# gset <- gset[ ,sel]
# 
# # order samples by group
# ex <- exprs(gset)[ , order(sml)]
# sml <- sml[order(sml)]
# fl <- as.factor(sml)
# labels <- c("NB","CD")
# 
# # set parameters and draw the plot
# palette(c("#c7ff9d","#f4dff4", "#AABBCC"))
# dev.new(width=4+dim(gset)[[2]]/5, height=6)
# par(mar=c(2+round(max(nchar(sampleNames(gset)))/2),4,2,1))
# title <- paste ("GSE1159", '/', annotation(gset), " selected samples", sep ='')
# boxplot(ex, boxwex=0.6, notch=T, main=title, outline=FALSE, las=2, col=fl)
# legend("topleft", labels, fill=palette(), bty="n")




