###########################################################################
# 知乎Live  R入门与基础绘图系统
# https://www.zhihu.com/lives/913412166131023872
# 
# 上机操作，第2部分
#
# content
# 1. volcano plot
# 2. heatmap plot
# 3. cytogram plot
# 4. plot with layout
###########################################################################
rm(list=ls())

###########################################################################
# part I
# volcano plot
###########################################################################
# data input is from cuffdiff output
cuffdiff_result = read.table(file="~/RNA_project/live_project/20171203-Live-R_partII.R/data_file/gene_exp.diff",header = TRUE)

wt_FPKM = cuffdiff_result$value_2
treat_FPKM = cuffdiff_result$value_1 
log2_foldchange = log2(treat_FPKM / wt_FPKM)
log2_foldchange[wt_FPKM == 0 ] = 0 

log2_foldchange[treat_FPKM == 0 ] = 0 

log10_p_value = log10(cuffdiff_result$p_value) * -1

# 1st edition
plot(x=log2_foldchange,y=log10_p_value,xlim=c(-4,4),ylim=c(0.001,4))

# 2nd edition
log10_p_value.filter = log10_p_value[log10_p_value >= 0.001]
log2_foldchange.filter = log2_foldchange[log10_p_value >= 0.001]
plot(x=log2_foldchange.filter,y=log10_p_value.filter,xlim=c(-4,4),ylim=c(0,4))

# 3rd edition 
plot(x=log2_foldchange.filter,y=log10_p_value.filter,
     xlim=c(-4,4),ylim=c(0,4),
     col=rgb(0,0,1,0.1),pch=16
     )

# 4th eidtion
length(log2_foldchange.filter)
col_vector = rep(rgb(0,0,1,0.1),length(log2_foldchange.filter))
col_vector[log10_p_value.filter >= -1* log10(0.001)] = rgb(1,0,0)
plot(x=log2_foldchange.filter,y=log10_p_value.filter,
     xlim=c(-4,4),ylim=c(0,4),
     col=col_vector,pch=16
     )

# 5th edition
# 1. p-value <= 0.05, 2. wt,treat FPKM >0 3. fold change > 2 OR 0.5
select_sign_vector = (cuffdiff_result$value_1 > 0 ) & (cuffdiff_result$value_2 >0) & (abs(log2_foldchange) >= 1) & (cuffdiff_result$value_1 >= 1 | cuffdiff_result$value_2 >= 1) & (cuffdiff_result$p_value <= 0.05)

log10_p_value.filter = log10_p_value[log10_p_value >= 0.001]
log2_foldchange.filter = log2_foldchange[log10_p_value >= 0.001]

select_sign_vector.filter = select_sign_vector[log10_p_value >= 0.001]

col_vector = rep(rgb(0,0,1,0.1),length(log2_foldchange.filter))
col_vector[select_sign_vector.filter] = rgb(1,0,0)

plot(x=log2_foldchange.filter,y=log10_p_value.filter,
     xlim=c(-4,4),ylim=c(0,4),
     col=col_vector,pch=16
)

abline(h=-1*log10(0.05),lwd=3,lty=3,col="#4C5B61")

###########################################################################
# part II
# heatmap plot
###########################################################################
# what is heatmap ??
# heatmap: marix -> rect + color  
# step1
plot(x=c(1:10),y=c(1:10),type="n")

# plot (0,0) (1,0) (1,1) (0,1)
rect(xleft = 0,ybottom = 0,xright = 1,ytop = 1,col=rgb(1,0,0,0.5))
rect(xleft = 5,ybottom = 5,xright = 6,ytop = 6,col="red")

# input matrix
# xleft,ybottom,xright,ytop
# color
# output image

input_matrix = matrix(c(1:36),6,6)

# set image size 
x_size = dim(input_matrix)[1]
y_size = dim(input_matrix)[2]

# 如果，我们从0,0 向 6,6方向找规律
# 先画，第1列，再画第2列

# 1st row
# xleft 0,0,0,0,0,0
# ybottom 5,4,3,2,1,0
# xright 1,1,1,1,1,1
# ytop 6,5,4,3,2,1

my_xleft = rep(c(0:(x_size-1)),each = x_size)
my_xright = my_xleft + 1
my_ybottom = rep(c((y_size-1):0),y_size)
my_ytop = my_ybottom + 1

# 下一步，我们确定颜色的问题
# 画最简单的颜色，最小值是白色，最大值是红色，中间线性变化

mat.max = max(input_matrix)
input_matrix.rate = input_matrix / mat.max
col.mat = rgb(1,0,0,input_matrix.rate)

plot(x=c(0:x_size),y=c(0:y_size),type="n",frame.plot = F,xaxt="n",yaxt="n",xlab="",ylab="")
rect(xleft = my_xleft,ybottom = my_ybottom,xright = my_xright,ytop = my_ytop,col=t(col.mat),border = F)


# 真实的Hi-C数据！！！
hic_mat.raw = read.table(file="~/RNA_project/live_project/20171203-Live-R_partII.R/data_file/chr_16_100000_MAPQ20.txt",sep = ",",header = F)


input_matrix = as.matrix(hic_mat.raw)
input_matrix = matrix(abs(rnorm(10000,1000,500)),100,100)

# 1st edition
x_size = dim(input_matrix)[1]
y_size = dim(input_matrix)[2]

my_xleft = rep(c(0:(x_size-1)),each = x_size)
my_xright = my_xleft + 1
my_ybottom = rep(c((y_size-1):0),y_size)
my_ytop = my_ybottom + 1

mat.max = quantile(input_matrix,prob=0.9)
input_matrix.rate = input_matrix / mat.max
input_matrix.rate[input_matrix.rate>1] = 1
col.mat = rgb(1,0,0,as.vector(as.matrix(input_matrix.rate)))

# plot(x=c(0:x_size),y=c(0:y_size),type="n",frame.plot = F,xaxt="n",yaxt="n",xlab="",ylab="")
# rect(xleft = my_xleft,ybottom = my_ybottom,xright = my_xright,ytop = my_ytop,col=t(col.mat),border = F)

# 为了快速画图，我们尝试直接把图片生成
# tiff(file="~/test_hic.png",width = 2000,height = 2000)
plot(x=c(0:x_size),y=c(0:y_size),type="n",frame.plot = F,xaxt="n",yaxt="n",xlab="",ylab="")
rect(xleft = my_xleft,ybottom = my_ybottom,xright = my_xright,ytop = my_ytop,col=col.mat,border = NA)
# dev.off()

# 我选择出 其中的第20列 到 第50列之间的区域
# xleft = 20
# xright = 50
# ybottom = 80
# ytop = 50

segments(x0 = 20,x1 = 20,y0 = 50,y1 = 80,lwd=5)
segments(x0 = 20,x1 = 50,y0 = 80,y1 = 80,lwd=5)
segments(x0 = 20,x1 = 50,y0 = 50,y1 = 50,lwd=5)
segments(x0 = 50,x1 = 50,y0 = 50,y1 = 80,lwd=5)
segments(x0 = 20,x1 = 50,y0 = 50,y1 = 80,lwd=5)



plot.matrix <- function(mat,bound.min=0,bound.max=1,mat.lim=NULL,color_type=1,col.min = "red",col.max = "red",col.boundary = NULL,n_block_color="#FFFFFF"){
  # mat = hic matrix
  # min_bound min quantile to miss data
  # max_bound max quantile to miss data
  
  mat <- as.matrix(mat)
  
  #matrix info calculate
  row_num <- dim(mat)[1]
  col_num <- dim(mat)[2]
  
  #matrix rects' coordinate
  x1 <- rep(c(0:(col_num-1)),each=row_num)
  x2 <- x1 + 1
  y1 <- rep(c(0:(-row_num+1)),col_num)
  y2 <- y1 -1
  
  #matrix colour vector
  if(color_type == 1){
    ## 生成矩阵需要的颜色
    ###确定 matrix的上下界
    mat.quantile = quantile(as.vector(mat),prob=c(bound.min,bound.max))
    
    if(is.null(mat.lim)){
      mat.lim = mat.quantile
    }else{
      mat.lim = c(min(c(mat.quantile,mat.lim)),max(c(mat.quantile,mat.lim)))
    }
    
    if(is.null(col.boundary)){
      col.boundary = mat.lim[1]
    }
    
    if(col.boundary>mat.lim[2]){
      print("Error! col.boundary have to smaller than matrix mat.lim[2] value!")
      return(NULL)
    }else if(col.boundary<mat.lim[1]){
      print("Error! col.boundary have to larger than matrix mat.lim[1] value!")
      print(mat.lim)
      print(col.boundary)
      return(NULL)
    }
    
    ## fix value,too large or too small
    mat[mat < mat.lim[1]] = mat.lim[1]
    mat[mat > mat.lim[2]] = mat.lim[2]
    
    ## create color vector
    mat.col_alpha = rep(0,length(as.vector(mat)))
    mat.col_alpha[mat >= col.boundary] = ceiling((mat[mat >= col.boundary] - col.boundary) / (mat.lim[2] - col.boundary) * 255)
    mat.col_alpha[mat < col.boundary] = ceiling((col.boundary - mat[mat < col.boundary]) / (col.boundary - mat.lim[1]) * 255)
    
    mat.color = rep("#FFFFFF",length(as.vector(mat)))
    mat.color[mat>=col.boundary] = rgb(t(col2rgb(col.max)),alpha = mat.col_alpha[mat>=col.boundary],maxColorValue = 255)
    mat.color[mat< col.boundary] = rgb(t(col2rgb(col.min)),alpha = mat.col_alpha[mat< col.boundary],maxColorValue = 255)
    
    mat.color = matrix(mat.color,nrow = row_num,ncol = col_num)
    mat.color[,colSums(mat)==0] = rgb(t(col2rgb(n_block_color)),maxColorValue = 255)
    mat.color[rowSums(mat)==0,] = rgb(t(col2rgb(n_block_color)),maxColorValue = 255)
    
  }else if(color_type == 2){
    # input a col_list 
  }
  
  # plot the final matrix
  plot(x=c(0,col_num),y=c(-row_num,0),type="n",frame.plot = F,xaxt="n",yaxt="n",cex.main = 2,xlab="",ylab="")
  rect(x1,y1,x2,y2,col = mat.color,border = NA)
}

###########################################################################
# part III
# cytogram plot
###########################################################################
chrom.name <- function(chrom_index){
  if(chrom_index <=22){
    chrom_name = paste0("chr",chrom_index)
  }else if(chrom_index == 23){
    chrom_name = "chrX"
  }else if(chrom_index == 24){
    chrom_name = "chrY"
  }else if(chrom_index == 25){
    chrom_name = "rDNA"
  }
}

chrom.name.factor <- function(chrom_index_list){
  chrom_index_list = as.vector(chrom_index_list)
  chrom_name_list = rep("chrN",length(chrom_index_list))
  for(index in c(1:length(chrom_index_list))){
    chrom_index = chrom_index_list[index]
    if(chrom_index < 23){
      chrom_name = paste0("chr",chrom_index)
    }else if(chrom_index == 23){
      chrom_name = "chrX"
    }else if(chrom_index == 24){
      chrom_name = "chrY"
    }else if(chrom_index == 25){
      chrom_name = "rDNA"
    }else{
      chrom_name = NA
    }
    chrom_name_list[index] = chrom_name
  }
  return(chrom_name_list)
}

chrom.length <- function(chrom_index,GenomeVersion="hg19"){
  # 根据chr_index返回chr_len
  HG19_LEN=c(249250621,243199373,198022430,191154276,180915260,171115067,159138663,146364022,141213431,135534747,
             135006516,133851895,115169878,107349540,102531392,90354753,81195210,78077248,59128983,63025520,48129895,
             51304566,155270560,59373566,42999)
  
  HG38_LEN=c(248956422,242193529,198295559,190214555,181538259,170805979,159345973,145138636,138394717,133797422,
             135086622,133275309,114364328,107043718,101991189,90338345,83257441,80373285,58617616,64444167,46709983,
             50818468,156040895,57227415,42999)
  if(GenomeVersion=="hg19"){
    return(HG19_LEN[chrom_index])
  }else if(GenomeVersion=="hg38"){
    return(HG38_LEN[chrom_index])
  }
}

hg19_g_band = read.table("~/menghw_HD/reference/hg19_G-band.txt",header = F,sep = "\t")
plot.chromosome <- function(chrom_index,chrom_start=0,chrom_end=250e6,track_height=8,border.lwd = 1,genome_version="hg19"){
  # initialize
  # chrom_index = 1
  # chrom_start = 0
  # chrom_end = chrom.length(chrom_index)
  df.peak = hg19_g_band
  
  # 将chrom_index生成bed文件中用的chrom_name
  if(chrom_index <=22){
    chrom_name = paste0("chr",chrom_index)
  }else if(chrom_index == 23){
    chrom_name = "chrX"
  }else if(chrom_index == 24){
    chrom_name = "chrY"
  }
  
  # 只画某1条染色体，筛选数据
  bed_table = df.peak[df.peak[,1]==chrom_name,]
  filter_vector = (bed_table[,2] >= chrom_start & bed_table[,2] <= chrom_end) & (bed_table[,3] >= chrom_start & bed_table[,3] <= chrom_end)
  bed_table = bed_table[filter_vector,]
  
  # 只使用有用的数据信息
  track_start = bed_table[,2]
  track_end = bed_table[,3]
  track_value = bed_table[,5]
  track_info = bed_table[,6]
  
  ## use 255 gray pattern
  alpha_vector = floor(track_value / 600 * 255) # 这里设置600纯是为了好看
  alpha_vector[alpha_vector>255] = 255
  color_vector = rgb(t(col2rgb("black")),alpha = alpha_vector,maxColorValue = 255)
  
  # rect 坐标
  track_x1 <- as.vector(track_start)
  track_x2 <- as.vector(track_end)
  track_y1 <- rep((-1)*(track_height %/% 2),nrow(bed_table))
  track_y2 <- rep(track_height %/% 2,nrow(bed_table))
  
  track_y1[track_info=="acen"] <- (-1)*(track_height %/% 2)/2
  track_y2[track_info=="acen"] <- (track_height %/% 2)/2
  
  plot(x=c(chrom_start,chrom_end),y=c((-1)*(track_height %/% 2),(track_height %/% 2)),type="n",frame.plot = F,cex.axis=2,cex.lab=2,ylab = "",xlab = "",xaxt = "n",yaxt="n")
  rect(track_x1,track_y1,track_x2,track_y2,col = color_vector,border = T,lwd = border.lwd)
}