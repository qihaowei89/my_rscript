#####################################################################################################################
# read and check parameter
#####################################################################################################################
# --input : input meth file
# --group : group file, the same length as --input file 
# --output : output dir
# --ctr : control labs
# --exp : experiment labs

if (!require(gplots)){
    install.packages("gplots")
    library(gplots)
    }
if (!require(ggplot2)){
    install.packages("ggplot2")
    library(ggplot2)
    }
if (!require(reshape2)){
    install.packages("reshape2")
    library(reshape2)
    }
if (!require(plyr)){
    install.packages("plyr")
    library(plyr)
    }
if (!require(stringr)){
    install.packages("stringr")
    library(stringr)
    }
if (!require(mice)){
    install.packages("mice") 
    library(mice)
    }
if (!require(VIM)){
    install.packages("VIM") 
    library(VIM)
}
# if (!require(magrittr)){
#     install.packages("magrittr") 
#     library(magrittr)
# }

# help function 
my_help <- function(){
    cat("meth.R Usage like:\n")
    cat("----------------------------------------------------------------------------------------------------\n")
    cat("Rscript meth.R --input XXX --group XXX --ctr XXX --exp XXX --output XXX\n")
    cat("----------------------------------------------------------------------------------------------------\n")
    cat("\t--input : input meth file\n")
    cat("\t--group : group file, the same length as --input file\n")
    cat("\t--ctr : control labs\n")
    cat("\t--exp : experiment labs\n")
    cat("\t--output : output dir\n")
}

# args <- commandArgs(T)

args <-  c("--input", "BML0065_2.txt", "--group", "group-BML0065_2.txt", "--ctr" ,"para-carcinoma","--exp", "carcinoma", "--output", "BML0065_2")

args.length = length(args)
# 1. check args length 
if(args.length != (5 * 2)){
    my_help()
    stop("Input parameters error!")
}

# 2. make parameter list
args.list = list()
for (index in seq(1,length(args), by = 2)){args.list[[args[index]]]  = args[index+1]}

# 3.check if all parameter in the list 
parameter_list = c("--input","--group","--ctr","--exp","--output")
if(!all(parameter_list %in% names(args.list))){
    my_help()
    stop("Input parameters keys error!")
}

meth_file  <- unlist(args.list["--input"])

group_file <- unlist(args.list["--group"])

group_labs <- unlist(c(args.list["--ctr"], args.list["--exp"]))

out_dir    <- unlist(args.list["--output"])


work.dir   <- sprintf("C:/Users/weiqi/Desktop/SNP_csv/work_dir/meth")
meth_file  <- sprintf("%s/%s", work.dir, meth_file)

group_file <- sprintf("%s/%s", work.dir, group_file)

dir.create(sprintf("%s/%s", work.dir, out_dir))


##########
#读取文件#
##########
raw <- read.table(meth_file, header = T, row.names = 1, sep = "\t",check.names = F)

############################
### 去除NA多的样本和位??? ###
############################
temp <- raw
ID <- colnames(temp)
row_call <- apply(temp,1,function(x)  {(1-length(which(is.na(x)))/length(colnames(temp)))})
temp <- temp[row_call != 0,]
col_call <- apply(temp,2,function(x)  {(1-length(which(is.na(x)))/length(rownames(temp)))})
#rm_ID <- ID[!col_call >= 1]

raw <- temp
# col_call.80 <- apply(raw,2,function(x)  {(1-length(which(is.na(x)))/length(rownames(raw)))})
# row_call.80 <- apply(raw,1,function(x)  {(1-length(which(is.na(x)))/length(colnames(raw)))})
# raw <- raw[row_call.80 > 0.8, col_call.80 > 0.8]


# row_out_list <- which(apply(raw,1,function(x)  all(is.na(x)))) 
# col_out_list <- which(apply(raw,2,function(x)  all(is.na(x)))) 
# col_list     <- which(apply(raw,2,function(x) !all(is.na(x)))) 
# row_list     <- which(apply(raw,1,function(x) !all(is.na(x)))) 
# raw <- raw[row_list, col_list]   

#row.names(raw) <- str_extract(row.names(raw),pattern = "\\d_CpG.*") 

################
# 读取分组文件 #
################
group.file <- read.table(group_file, header = T, stringsAsFactors = F)
#group.file<- group.file[col_call >= 1 ,]
raw_t <- t(raw)
if (identical(as.character(group.file$ID), rownames(raw_t))){
    raw_3 <- data.frame(group.file,raw_t,row.names = 1)
}

#############################
## 定义函数summarySE计算SE ##
#############################
summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=T, conf.interval=.95, .drop=TRUE) {
    length2 <- function (x, na.rm=FALSE) {
        if (na.rm) sum(!is.na(x))
        else length(x)
    }
    datac <- ddply(data, groupvars, .drop=.drop,
                   .fun = function(xx, col) {
                       c(N    = length2(xx[[col]], na.rm=na.rm),
                         mean = mean   (xx[[col]], na.rm=na.rm),
                         sd   = sd     (xx[[col]], na.rm=na.rm)
                       )
                   },
                   measurevar
    )
    datac <- rename(datac, c("mean" = measurevar))
    datac$se <- datac$sd / sqrt(datac$N)
    ciMult <- qt(conf.interval/2 + .5, datac$N-1)
    datac$ci <- datac$se * ciMult
    return(datac)
}

##################
## 简单汇总数??? ##
##################
md <- melt(raw_3,id = "group")
names(md) <- c("group","site","methy")
dd <- md[order(md[,2],md[,1]),]  

################
## 绘制箱式??? ##
################
dd$site <- factor(dd$site)
pdf(sprintf("%s/%s/%s_boxplot.pdf",work.dir,out_dir,out_dir))
p <- ggplot(dd, aes(x=site, y=methy,fill=group)) + geom_boxplot()
p +  theme(axis.text.x = element_text(angle=40, hjust=1, vjust=1))
dev.off()
# da <- dd
# da$methy <- log10(dd$methy)
# p <- ggplot(da, aes(x=site, y=methy,fill=group)) + geom_boxplot()
# p + theme(axis.text.x = element_text(angle=40, hjust=1, vjust=1))
# p

png(sprintf("%s/%s/%s_boxplot.png",work.dir,out_dir,out_dir))
dd$site <- factor(dd$site)
p <- ggplot(dd, aes(x=site, y=methy,fill=group)) + geom_boxplot()
#p + geom_dotplot(binaxis='y', stackdir='center',position=position_dodge(1),dotsize=0.2)+
p + theme(axis.text.x = element_text(angle=40, hjust=1, vjust=1))
dev.off()


############
## 计算SE ##
############
tgc2 <- summarySE(dd, measurevar="methy", groupvars=c("group","site"));
tgc2$site <- factor(tgc2$site)

################
## 绘制条形??? ##
################
pdf(sprintf("%s/%s/%s_barplot.pdf",work.dir,out_dir,out_dir))
ggplot(tgc2, aes(x=site, y=methy, fill=group)) + 
    geom_bar(position=position_dodge(), stat="identity") +
    geom_errorbar(aes(ymin=methy-se, ymax=methy+se),width=.2,position=position_dodge(1))+
    theme(axis.text.x = element_text(angle=40, hjust=1, vjust=1))
dev.off()

png(sprintf("%s/%s/%s_barplot.png",work.dir,out_dir,out_dir))
ggplot(tgc2, aes(x=site, y=methy, fill=group)) + 
    geom_bar(position=position_dodge(), stat="identity") +
    geom_errorbar(aes(ymin=methy-se, ymax=methy+se),width=.2,position=position_dodge(1))+
    theme(axis.text.x = element_text(angle=40, hjust=1, vjust=1))
dev.off()


################
## 绘制密度??? ##
################
# Group <- group_labs
pdf(sprintf("%s/%s/%s_meth_level.pdf",work.dir,out_dir,out_dir))
p <- ggplot(tgc2,aes(methy))
p + geom_histogram(position = 'identity',alpha=0.5, aes(y = ..density.., fill = group)) +
    stat_density(geom = 'line', position = 'identity', aes(colour = group))+xlim(0,1)

p <- ggplot(tgc2,aes(methy))
p + geom_histogram(position = 'identity',alpha=0.5, aes(y = ..density..,fill = group)) +
    stat_density(geom = 'line', position = 'identity', aes(colour = group))+facet_grid(.~group)+xlim(0,1)
dev.off()

png(sprintf("%s/%s/%s_meth_level.png",work.dir,out_dir,out_dir))
p <- ggplot(tgc2,aes(methy))
p + geom_histogram(position = 'identity', alpha=0.5, aes(y = ..density.., fill = group)) +
    stat_density(geom = 'line', position = 'identity', aes(colour = group))+xlim(0,1)

p <- ggplot(tgc2,aes(methy))
p + geom_histogram(position = 'identity', alpha=0.5, aes(y = ..density.., fill = group)) +
    stat_density(geom = 'line', position = 'identity', aes(colour = group))+facet_grid(.~group)+xlim(0,1)
dev.off()

###############
## 绘制PCA??? ##
###############

raw_3_1 <- raw_3[,c(-1)]

aggr_plot <- aggr(raw_3_1, col = c('navyblue', 'red'), numbers=TRUE, sortVars=TRUE, 
                  labels=names(raw_3_1), cex.axis=.7, gap=3, 
                  ylab=c("Histogram of missing data", "Pattern"))


imp <- mice(raw_3_1, seed = 1234,m=5,maxit=50,printFlag=F)


# densityplot(imp,~CpG_2|.imp)
# densityplot(imp,~CpG_7.8|.imp)
# densityplot(imp,~CpG_5.6|.imp)
# densityplot(imp,~CpG_3|.imp)
# densityplot(imp,~CpG_4|.imp)
# densityplot(imp,~CpG_9.10|.imp)
stripplot(imp, pch = 20, cex = 1.2)

raw_3_1 <- complete(imp, action = 1)
data.pca <- prcomp(raw_3_1, scale.= F)
data.pca$x -> a
as.data.frame(a) -> b
groups <- raw_3$group

pdf(sprintf("%s/%s/%s_PCA.pdf",work.dir,out_dir,out_dir),width = 8 ,height = 8 )
ggplot(b)+
  aes(x=PC1,y=PC2,colour=groups,label = row.names(b))+
  geom_point(aes(x=PC1,y=PC2,colour=groups),size = 3)+
  geom_text(vjust = 3,size=3)
dev.off()

png(sprintf("%s/%s/%s_PCA.png",work.dir,out_dir,out_dir),width = 700,height = 700,units = "px",pointsize = 10)
ggplot(b)+
  aes(x=PC1,y=PC2,colour=groups,label = row.names(b))+
  geom_point(aes(x=PC1,y=PC2,colour=groups),size = 4)+
  geom_text(vjust = 3,size=3)
dev.off()

# pdf(sprintf("%s/%s/%s_pca.barpolt.#pdf",work.dir,out_dir,out_dir))
# pca.sum <- summary(data.pca)
# percent <- pca.sum$importance[2,]*100
# percent_pc <- as.data.frame(percent)
# percent_pc <- percent_pc[which(percent_pc$percent!=0),]
# ggplot(data=percent_pc, aes(x=row.names(percent_pc), y=percent)) +
#     geom_bar(stat="identity",fill= "skyblue")
# dev.off()


# barplot(percent, col = "skyblue",
#         xlab = "PC", ylab = "Percent(%)",
#         ylim =c(1.2*min(pca.sum$importance[2,]*100),1.2*max(pca.sum$importance[2,]*100)))


#################### 
####  绘制热图  ####
####################
data.matrix <- as.matrix(raw)
# group1 <- grep(pattern = "*-N", x = colnames(data.matrix))
# group2 <- grep(pattern = "*-C", x = colnames(data.matrix))
# data.matrixt <- cbind(data.matrix[,c(group1)],data.matrix[,c(group2)])

##设定分组颜色
side_cols <- rep(1,length(group.file$group))
side_cols[group.file$group == levels(factor(group.file$group))[1]] <- "purple"
side_cols[group.file$group == levels(factor(group.file$group))[2]] <- "orange"
pdf(sprintf("%s/%s/%s_heatmap.pdf", work.dir,out_dir,out_dir),width = 7,height = 8)
heatmap <- heatmap.2(
    #t(scale(data.matrix[,-82],scale = F,center =F)),
    scale(data.matrix,scale = T,center =F),
    col=greenred,
    #col=heat.colors,
    symkey=FALSE, 
    #scale = c("none"),
    #scale = c("row"),
    #scale = c("column"),
    density.info="none", 
    trace="none",
    cexRow = 1,
    keysize = 1,
    #key.title=NA,
    key.xlab="value",
    #RowSideColors = side_cols,
    ColSideColors = side_cols,
    srtRow=0, adjRow=c(0, 1), srtCol=90, adjCol=c(1,1),  ## Show effect of row and column label rotation
    offsetRow=0.5, offsetCol= 0.5,    #离坐标轴距离
    reorderfun=function(d,w) reorder(d, w, agglo.FUN=mean),
    # distfun=function(x) as.dist(1-cor(t(x))),
    # hclustfun=function(x) hclust(x, method="complete"),
    distfun=function(x) dist(x, method="euclidean"),
    hclustfun=function(x) hclust(x, method="ward.D2"),
    main ="Heatmap of CpG methylation sites",
    margin=c(5,6),
    #na.color="black"
    na.color="gray"
)

heatmap.2(
    scale(data.matrix,scale =T,center = F), 
    col=greenred,
    #col=bluered,
    dendrogram = "none",
    symkey=FALSE, 
    #scale = c("none"),
    #scale = c("row"),
    #scale = c("column"),
    density.info="none", 
    trace="none",
    cexRow = 1,
    keysize = 1,
    #key.title=NA,
    key.xlab="value",
    ColSideColors = side_cols,
    srtRow=0, adjRow=c(0, 1), srtCol=90, adjCol=c(1,1),  ## Show effect of row and column label rotation
    offsetRow=1, offsetCol= 0.5,    #离坐标轴距离
    #reorderfun=function(d,w) reorder(d, w, agglo.FUN=mean), 
    #distfun=function(x) as.dist(1-cor(t(x))), 
    #hclustfun=function(x) hclust(x, method="complete")
    reorderfun=function(d,w) reorder(d, w, agglo.FUN=mean), 
    distfun=function(x) dist(x, method="euclidean"), 
    hclustfun=function(x) hclust(x, method="ward.D2"),
    main ="Heatmap of CpG methylation sites",
    margin=c(5,8),
    #na.color="black",
    na.color="gray",
    Rowv=F,Colv=FALSE
)
dev.off()


png(sprintf("%s/%s/%s_heatmap.png",work.dir,out_dir,out_dir))
heatmap <- heatmap.2(
    #t(scale(data.matrix[,-82],scale = F,center =F)),
    scale(data.matrix,scale = T,center =F),
    col=greenred,
    #col=heat.colors,
    symkey=FALSE, 
    #scale = c("none"),
    #scale = c("row"),
    #scale = c("column"),
    density.info="none", 
    trace="none",
    cexRow = 1,
    keysize = 1,
    #key.title=NA,
    key.xlab="value",
    #RowSideColors = side_cols,
    ColSideColors = side_cols,
    srtRow=0, adjRow=c(0, 1), srtCol=90, adjCol=c(1,1),  ## Show effect of row and column label rotation
    offsetRow=0.5, offsetCol= 0.5,    #离坐标轴距离
    reorderfun=function(d,w) reorder(d, w, agglo.FUN=mean),
    # distfun=function(x) as.dist(1-cor(t(x))),
    # hclustfun=function(x) hclust(x, method="complete"),
    # reorderfun=function(d,w) reorder(d, w, agglo.FUN=mean),
    distfun=function(x) dist(x, method="euclidean"),
    hclustfun=function(x) hclust(x, method="ward.D2"),
    main ="Heatmap of CpG methylation sites",
    margin=c(5,12),
    #na.color="black"
    na.color="gray"
)

heatmap.2(
    scale(data.matrix,scale = T,center = F), 
    col=greenred,
    #col=bluered,
    symkey=FALSE, 
    #scale = c("none"),
    #scale = c("row"),
    #scale = c("column"),
    density.info="none", 
    trace="none",
    cexRow = 1,
    keysize = 1,
    #key.title=NA,
    key.xlab="value",
    ColSideColors = side_cols,
    srtRow=0, adjRow=c(0, 1), srtCol=90, adjCol=c(1,1),  ## Show effect of row and column label rotation
    offsetRow=1, offsetCol= 0.5,    #离坐标轴距离
    #reorderfun=function(d,w) reorder(d, w, agglo.FUN=mean), 
    #distfun=function(x) as.dist(1-cor(t(x))), 
    #hclustfun=function(x) hclust(x, method="complete")
    reorderfun=function(d,w) reorder(d, w, agglo.FUN=mean), 
    distfun=function(x) dist(x, method="euclidean"), 
    hclustfun=function(x) hclust(x, method="ward.D2"),
    main ="Heatmap of CpG methylation sites",
    margin=c(5,12),
    #na.color="black",
    na.color="gray",
    Rowv=FALSE,Colv=FALSE
)
dev.off()

##################
## 简单差异分??? ##
##################
options(digits=3)
raw_groupA <- as.matrix(raw[,which(group.file$group == group_labs[1])])
apply(raw_groupA,1,mean,na.rm = T) -> Mean_methy_level_Control
raw_groupB <- as.matrix(raw[,which(group.file$group == group_labs[2])]) 
apply(raw_groupB,1,mean,na.rm = T) -> Mean_methy_level_Resistance 
fold <- Mean_methy_level_Resistance / Mean_methy_level_Control
diff <- Mean_methy_level_Resistance - Mean_methy_level_Control

#################
## 计算四分位数 #
#################
# t(apply(raw_groupA,1,fivenum,na.rm = T)) -> five_num_Control
# colnames(five_num_Control) <- c("0%", "25%", "50%","75%","100%")
# t(apply(raw_groupB,1,fivenum,na.rm = T) )-> five_num_Resistance
# colnames(five_num_Resistance) <- c("0%", "25%", "50%","75%","100%")
# five_num <- cbind(group=tgc2$group,site=as.character(tgc2$site),rbind(five_num_Control, five_num_Resistance),tgc2[,c(-1,-2,-3)])

#####################################
## 样本正态性检验、t检验、秩和检??? ##
#####################################
T_test_P_value <-c()
wilcox_test_P_value <-c()
pnorm_A <- c()
pnorm_B <- c()
for (i in 1:length(raw[,1])){
    pnorm_A[i] <- ks.test(raw_groupA[i,],"pnorm")$p.value 
    pnorm_B[i] <- ks.test(raw_groupB[i,],"pnorm")$p.value
    T_test_P_value[i] <- t.test(raw_groupA[i,],raw_groupB[i,],paired = F)$p.value
    wilcox_test_P_value[i] <- wilcox.test(raw_groupA[i,],raw_groupB[i,],paired = F)$p.value
}


control_col_mean <- sprintf("Mean_methy_level_%s",group_labs[1])
experiment_col_mean <- sprintf("Mean_methy_level_%s",group_labs[2])
pnorm_col <- sprintf("pnorm_%s",group_labs[1])
pnorm_exp <- sprintf("pnorm_%s",group_labs[2])

fram <- data.frame(site= row.names(raw),Mean_methy_level_Control,Mean_methy_level_Resistance,fold,diff, T_test_P_value, wilcox_test_P_value, pnorm_A, pnorm_B) 
fram <- rename(fram,replace = c("Mean_methy_level_Control" = control_col_mean ,"Mean_methy_level_Resistance" = experiment_col_mean,"pnorm_A" =pnorm_col, "pnorm_B" = pnorm_exp))

# if (wilcox_test_P_value >= 3){
#     
#     
#     pdf(sprintf("%s/%s/%s_heatmap.pdf", work.dir,out_dir,out_dir),width = 12,height = 9)
#     heatmap <- heatmap.2(
#         scale(data.matrix,scale = T,center =T),
#         col=greenred,
#         #col=bluered,
#         symkey=FALSE, 
#         #scale = c("none"),
#         #scale = c("row"),
#         #scale = c("column"),
#         density.info="none", 
#         trace="none",
#         cexRow = 1,
#         keysize = 1,
#         #key.title=NA,
#         key.xlab="value",
#         ColSideColors = side_cols,
#         srtRow=0, adjRow=c(0, 1), srtCol=45, adjCol=c(1,1),  ## Show effect of row and column label rotation
#         offsetRow=1, offsetCol= 0.5,    #离坐标轴距离
#         distfun=function(x) as.dist(1-cor(t(x))),
#         hclustfun=function(x) hclust(x, method="complete"),
#         reorderfun=function(d,w) reorder(d, w, agglo.FUN=mean),
#         # distfun=function(x) dist(x, method="euclidean"),
#         # hclustfun=function(x) hclust(x, method="ward.D2"),
#         # reorderfun=function(d,w) reorder(d, w, agglo.FUN=mean),
#         main ="Heatmap of CpG methylation sites",
#         margin=c(5,10),
#         #na.color="black"
#         na.color="white"
#     )
#     
#     heatmap.2(
#         scale(data.matrix,scale = T,center = T), 
#         col=greenred,
#         #col=bluered,
#         symkey=FALSE, 
#         #scale = c("none"),
#         #scale = c("row"),
#         #scale = c("column"),
#         density.info="none", 
#         trace="none",
#         cexRow = 1,
#         keysize = 1,
#         #key.title=NA,
#         key.xlab="value",
#         ColSideColors = side_cols,
#         srtRow=0, adjRow=c(0, 1), srtCol=45, adjCol=c(1,1),  ## Show effect of row and column label rotation
#         offsetRow=1, offsetCol= 0.5,    #离坐标轴距离
#         #reorderfun=function(d,w) reorder(d, w, agglo.FUN=mean), 
#         #distfun=function(x) as.dist(1-cor(t(x))), 
#         #hclustfun=function(x) hclust(x, method="complete")
#         reorderfun=function(d,w) reorder(d, w, agglo.FUN=mean), 
#         distfun=function(x) dist(x, method="euclidean"), 
#         hclustfun=function(x) hclust(x, method="ward.D2"),
#         main ="Heatmap of CpG methylation sites",
#         margin=c(5,10),
#         #na.color="black",
#         na.color="white",
#         Rowv=FALSE,Colv=FALSE
#     )
#     dev.off()
# }
diff_file <- sprintf("%s/%s/%s_test.csv",work.dir,out_dir,out_dir)
# five_file <- sprintf("%s/%s/%s_统计.csv",out_dir,out_dir1)

write.table(fram,diff_file, sep = ",", row.names = F,quote = F)
# write.table(five_num,five_file, sep = ",", row.names = F,quote = F)









