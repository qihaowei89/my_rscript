#QTL·ÖÎö
library(reshape)
file <- read.table("TERF2.txt",header = T,sep = "\t",na.strings = "")
file <- na.omit(file)
file$TERF2rs153045 <-as.factor(file$TERF2rs153045) 
#file$RNA.TERF2 <- log2(file$RNA.TERF2)
file.melt <- melt(file,id="TERF2rs153045")
file.melt <- data.frame(ID=1:length(file.melt$TERF2rs153045),file.melt)[,-3]
colnames(file.melt) <- c("ID","SNP","Value")
file.aa <- cast(file.melt,~SNP,fun.aggregate = mean)
plot(x = file.melt$SNP,y = file.melt$Value)

hist(log2(file.melt$Value[file.melt$SNP=="C"]))
hist(log2(file.melt$Value[file.melt$SNP=="CT"]))
hist(log2(file.melt$Value[file.melt$SNP=="T"]))


###################################################
### code chunk number 1: setup
###################################################
source("https://bioconductor.org/biocLite.R")
biocLite("DOQTL")
biocLite("MUGAExampleData")

### R code from vignette source 'QTL_Mapping_DO_Mice.Rnw'

###################################################
### code chunk number 1: setup
###################################################
library(DOQTL)
library(MUGAExampleData)
data(pheno)
data(model.probs)
#browseVignettes("DOQTL")
"""QTL mapping requires phenotype and genotype data"""
""" phenotypes called pheno"""
"""a 3D array of founder haplotype contributions 
   (num.samples x 8 founders x num.markers) called model.probs."""
"""The sample IDs must be in rownames(pheno) and dimnames(model.probs)[[1]]"""
"""hemoglobin distribution width at time point 2 (HDW2)"""
###################################################
### code chunk number 2: kinship
###################################################
#First, we need to create a kinship matrix using the founder contributions.
K = kinship.probs(model.probs)


###################################################
### code chunk number 3: covar
###################################################
#Second, we need to create a matrix of additive covariates to run in the model.
covar = data.frame(sex = as.numeric(pheno$Sex == "M"), diet = as.numeric(pheno$Diet == "hf"))
rownames(covar) = rownames(pheno)


###################################################
### code chunk number 4: snps
###################################################
#Third, we need to get the marker locations on the array.
load(url("ftp://ftp.jax.org/MUGA/muga_snps.Rdata"))


###################################################
### code chunk number 5: scanone
###################################################
#Fourth, we map the phenotype using scanone
qtl = scanone(pheno = pheno, pheno.col = "HDW2", probs = model.probs, K = K,
              addcovar = covar, snps = muga_snps)


###################################################
### code chunk number 6: perms
###################################################
#Fifth, we run permutations to determine significane thresholds. We recommend running at least 1,000
#permutations. In this demo, we run 100 permutations to save time.
perms = scanone.perm(pheno = pheno, pheno.col = "HDW2", probs = model.probs, 
                     addcovar = covar, snps = muga_snps, nperm = 100)
thr = quantile(perms, probs = 0.95)


###################################################
### code chunk number 7: qtlplot
###################################################
#We then plot the LOD curve for the QTL.

plot(qtl, sig.thr = thr, main = "HDW2")


###################################################
### code chunk number 9: coefplot
###################################################
coefplot(qtl, chr = 9)


###################################################
### code chunk number 11: QTL_Mapping_DO_Mice.Rnw:177-179
###################################################
interval = bayesint(qtl, chr = 9)
interval   #The QTL support interval is 4.7 Mb wide.


###################################################
### code chunk number 12: mergeplot
###################################################
ma = assoc.map(pheno = pheno, pheno.col = "HDW2", probs = model.probs, K = K,
               addcovar = covar, snps = muga_snps, chr = interval[1,2], 
               start = interval[1,3], end = interval[3,3])
 
tmp = assoc.plot(ma, thr = 4)
unique(tmp$sdps)


###################################################
### code chunk number 14: get_genes
###################################################
mgi = get.mgi.features(chr = interval[1,2], start = interval[1,3],
                       end = interval[3,3], type = "gene", source = "MGI")
nrow(mgi)
head(mgi)


###################################################
### code chunk number 15: sessionInfo
###################################################
sessionInfo()


