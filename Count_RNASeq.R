##################################################################################################################################
# 使用DESeq2 进行RNA-Seq差异表达分析
# Data 2017-12-14
# Author Howard MENG
# E-mail meng_howard@126.com
##################################################################################################################################
rm(list=ls())

#####################################################################################################################
# read and check parameter
#####################################################################################################################
# --input : input bam file list, format like bam1,bam2,bam3...
# --label : the same length as --input, format likt label1,label2,label3...
# --type : the sample type ,the same length as --input, format likt ctrl,ctrl,ko,ko...
# --threads : the CPU limit 
# --output : name of RData that contain SummarizedExperiment.obj

# help function 
my_help <- function(){
  cat("Count_RNASeq.R Usage like:\n")
  cat("----------------------------------------------------------------------------------------------------\n")
  cat("Rscript Count_RNASeq.R --input XXX --label XXX --type XXX --threads XXX --output XXX\n")
  cat("----------------------------------------------------------------------------------------------------\n")
  cat("\t--input : input bam file list, format like bam1,bam2,bam3...\n")
  cat("\t--label : the same length as --input, format likt label1,label2,label3...\n")
  cat("\t--type : the sample type ,the same length as --input, format likt ctrl,ctrl,ko,ko...\n")
  cat("\t--threads : the CPU limit \n")
  cat("\t--output : name of RData that contain SummarizedExperiment.obj\n")
}

# required package
require(magrittr)

# 1. check args length 
args = commandArgs()
args.length = length(args)

if(args.length != (5 + 5 * 2)){
  my_help()
  stop("Input parameters error!")
}

# 2. make parameter list
args.list = list()
for (index in seq(6,length(args),by = 2) ){
  args.list[[args[index]]]  = args[index+1]
}

# 3.check if all parameter in the list 
parameter_list = c("--input","--label","--type","--threads","--output")
if(! all(parameter_list %in% names(args.list))){
  my_help()
  stop("Input parameters keys error!")
}

# args.list[["--input"]] = "/home/menghw/RNA_project/songjh_project/tmp.data/tophat_result/293T-KO1_WT_rep1/accepted_hits.bam,/home/menghw/RNA_project/songjh_project/tmp.data/tophat_result/293T-KO1_WT_rep2/accepted_hits.bam,/home/menghw/RNA_project/songjh_project/tmp.data/tophat_result/293T-KO1_Pus1_rep1/accepted_hits.bam,/home/menghw/RNA_project/songjh_project/tmp.data/tophat_result/293T-KO1_Pus1_rep2/accepted_hits.bam"
# args.list[["--label"]] = "KO1.WT.rep1,KO1.WT.rep2,KO1.PUS1.rep1,KO1.PUS1.rep2"
# args.list[["--type"]] = "WT,WT,KO.PUS1,KO.PUS1"
# args.list[["--threads"]] = "4"
# args.list[["--output"]] = "~/test.RData"

# 4. reads parameters
INPUT_BAM = strsplit(args.list[["--input"]],split=",") %>% unlist()
INPUT_LABEL = strsplit(args.list[["--label"]],split=",") %>% unlist()
INPUT_TYPE = strsplit(args.list[["--type"]],split=",") %>% unlist()
CPU_THREADS = strsplit(args.list[["--threads"]],split=",") %>% unlist() %>% as.integer()
OUTPUT_FILE = strsplit(args.list[["--output"]],split=",") %>% unlist() %>% path.expand()

OUTPUT_DIR = dirname(OUTPUT_FILE)
OUTPUT_BASENAME = basename(OUTPUT_FILE)

# 5. check if BAM file exist 
if(all(file.exists(INPUT_BAM))){
  print("All file checked exsit!")
}else{
  stop(sprintf("File: '%s' Not Exist!",INPUT_BAM[!file.exists(INPUT_BAM)]))
}

######################################################################################
# required packages 
######################################################################################
require(BiocParallel)
require(GenomicFeatures)
require(GenomicAlignments)
require(Rsamtools)
require(TxDb.Hsapiens.UCSC.hg19.knownGene)

######################################################################################
# read BAM file
######################################################################################
input_bamfiles = Rsamtools::BamFileList(INPUT_BAM)
names(input_bamfiles) = INPUT_LABEL

######################################################################################
# read GTF file and make TxDb object
######################################################################################
# gtffile = file.path(GTF_FILE)
# txdb = makeTxDbFromGFF(gtffile,format = "gtf",circ_seqs = character())
# load(file="~/menghw_HD/reference/gene_gtf.TxDb/RefSeq.gene.hg19.txdb.RData")

######################################################################################
# load GTF file
######################################################################################
print("Loading GTF file...")
gene_range = GenomicFeatures::exonsBy(TxDb.Hsapiens.UCSC.hg19.knownGene,by="gene")
print("Done!")


######################################################################################
# counting gene reads
######################################################################################
print("counting...")

# BiocParallel::register(SerialParam())
BiocParallel::register(MulticoreParam(workers=CPU_THREADS))

SummarizedExperiment.obj <- GenomicAlignments::summarizeOverlaps(features=gene_range, 
                                                                 reads=input_bamfiles,
                                                                 mode="Union",
                                                                 singleEnd=FALSE,
                                                                 ignore.strand=TRUE,
                                                                 fragments=TRUE )
print("counting...Done!")

######################################################################################
# save data
######################################################################################
# make sample table 
sample_table = DataFrame(case_label = INPUT_LABEL,
                         case_type = INPUT_TYPE,
                         bam_path = INPUT_BAM)

colData(SummarizedExperiment.obj) = sample_table

# save SummarizedExperiment.obj
save(file=OUTPUT_FILE,list="SummarizedExperiment.obj")
print("Save R image... Done!")

rm(list=ls())
print("Clear R env... Done!")



