st <- c(2000000, 2070000, 2100000, 2160000)
ed <- c(2050000, 2130000, 2150000, 2170000)
str <- c("-", "+", "-", "-")
gr <- c("Group1","Group2","Group1", "Group3")
annTrack <- AnnotationTrack(start=st, end=ed, strand=str, chromosome=7,
                            genome="hg19", feature="test", group=gr,
                            id=paste("annTrack item", 1:4),
                            name="annotation track foo",
                            stacking="squish")

ax <- GenomeAxisTrack()

dt <- DataTrack(start=seq(min(st), max(ed), len=10), width=18000,
                data=matrix(runif(40), nrow=4), genome="hg19", chromosome=7,
                type="histogram", name="data track bar")




## Now plot the tracks
res <- plotTracks(list(ax, annTrack, dt))

## Plot only a subrange
res <- plotTracks(list(ax, annTrack, dt), from=2080000, to=2156000)

## Extend plotting ranges
res <- plotTracks(list(ax, annTrack, dt), extend.left=200000, extend.right=200000)

## Add a header
res <- plotTracks(list(ax, annTrack, dt), main="A GenomGraphs plot",
                  col.main="darkgray")

## Change vertical size and title width
res <- plotTracks(list(ax, annTrack, dt), sizes=c(1,1,5))

names(annTrack) <- "foo"
res <- plotTracks(list(ax, annTrack), title.width=0.6)

## Adding and lattice like plots
library(grid)
grid.newpage()
pushViewport(viewport(height=0.5, y=1, just="top"))
grid.rect()
plotTracks(annTrack, add=TRUE)
popViewport(1)
pushViewport(viewport(height=0.5, y=0, just="bottom"))
grid.rect()
plotTracks(dt, add=TRUE)
popViewport(1)

## Not run: 
library(lattice)
myPanel <- function(x, ...) plotTracks(annTrack, panel.only=TRUE, from=min(x), to=max(x), shape="box")
a <- seq(1900000, 2250000, len=40)
xyplot(b~a|c, data.frame(a=a, b=1, c=cut(a, 4)), panel=myPanel,
       scales=list(x="free"))

## End(Not run)





library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
library(org.Hs.eg.db)
library(limma)
ann <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)

# Randomly select 1000 CpGs to be significantly differentially methylated
sigcpgs <- sample(rownames(ann450k),1000,replace=FALSE)

# All CpG sites tested
allcpgs <- rownames(ann450k)

# Use org.Hs.eg.db to extract a GO term
GOtoID <- toTable(org.Hs.egGO2EG)
setname1 <- GOtoID$go_id[1]
setname1
keep.set1 <- GOtoID$go_id %in% setname1
set1 <- GOtoID$gene_id[keep.set1]
setname2 <- GOtoID$go_id[2]
setname2
keep.set2 <- GOtoID$go_id %in% setname2
set2 <- GOtoID$gene_id[keep.set2]

# Make the gene sets into a list
sets <- list(set1, set2)
names(sets) <- c(setname1,setname2)

# Testing with prior probabilities taken into account
# Plot of bias due to differing numbers of CpG sites per gene
gst <- gsameth(sig.cpg = sigcpgs, all.cpg = allcpgs, collection = sets, plot.bias = TRUE, prior.prob = TRUE)
topGSA(gst)

# Testing ignoring bias
gst.bias <- gsameth(sig.cpg = sigcpgs, all.cpg = allcpgs, collection = sets, prior.prob = FALSE)
topGSA(gst.bias)




########################


myDMP <- champ.DMP(arraytype="EPIC")
save(myDMP,file="myDMP.rda")
DMP.GUI()
myDMR <- champ.DMR(arraytype = "EPIC",method="DMRcate",cores=1)
save(myDMR,file="myDMR.rda")
DMR.GUI(arraytype="EPIC")


getwd()

files  <- list.files(pattern = "*.xls$")

library(xlsx)
for (i in files){
    temp <- read.xlsx(i,1,encoding = "UTF-8")
    write.table(temp,"conbind.csv",append = TRUE,quote = FALSE,row.names = FALSE,col.names = T,sep = ",")
}


