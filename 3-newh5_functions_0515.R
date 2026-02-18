library(Matrix)
library(clusterProfiler)
library(org.Hs.eg.db)
library(tidyverse)
library(patchwork)
library(rhdf5) 
library(tools)
matrix_to_csr <- function(mat) {
  mat <- as(mat, "CsparseMatrix")  # 转换为稀疏矩阵
  list(
    data = as.vector(mat@x),
    indices = as.integer(mat@i),
    indptr = as.integer(mat@p)
  )
}

# input=原始.h5所在文件夹
# h5_file=原始.h5文件名
# output = 新文件夹
# newh5_file = 新的.h5文件名
# df = read.csv('pre.csv or gt.csv',row.names = 1,check.names =F) 

newh5<- function(input,h5_file,output,newh5_file,df){
  input_h5<-file.path(input,h5_file)
  barcode <- h5read(input_h5,"/matrix/barcodes")# 读取原始数据(filtered_feature_bc_matrix.h5)的spot barcode
  df = df[,barcode] ####################### 保证了新的df与原始数据中的spot顺序
  csr <- as.matrix(apply(df, 2, as.numeric)) %>% matrix_to_csr()
  # 拼接目标文件路径
  target_file <- file.path(output, basename(input_h5))
  file.copy(input_h5, output)
  message(paste("文件",input_h5,"已复制到",output))
  # 拼接目标文件路径
  spatial_file <- file.path(output, 'spatial')
  file.copy(file.path(input, 'spatial'), output,recursive = TRUE)
  message(paste("文件",input_h5,"已复制到",output))
  # file.copy(input_h5,output)
  # file.copy(paste0(input,'spatial'),output,recursive=T)
  output_newh5 <- file.path(output,newh5_file)
  h5createFile(output_newh5)
  h5createGroup(output_newh5, "/matrix")
  h5createGroup(output_newh5, "/matrix/features")
  # 含0稀疏矩阵csr的
  h5write(csr$data, output_newh5, "/matrix/data")
  h5write(csr$indices, output_newh5, "/matrix/indices")
  h5write(csr$indptr, output_newh5, "/matrix/indptr")
  genenames <- rownames(df)
  shap=dim(df)
  all_tag_keys <-h5read(input_h5,'/matrix/features/_all_tag_keys')
  feature_type<-h5read(input_h5,'/matrix/features/feature_type')[seq(shap[1])]
  genome<-h5read(input_h5,'/matrix/features/genome')[seq(shap[1])]
  symbol2ENSG <- bitr(rownames(df),fromType = 'SYMBOL',toType = 'ENSEMBL',OrgDb = 'org.Hs.eg.db',drop=F)
  id <- subset(symbol2ENSG,!duplicated(SYMBOL));dim(id)
  ENSG <- id$ENSEMBL#;ENSG
  nas <- sum(is.na(ENSG))
  ENSG[is.na(ENSG)]<-paste0('map-failed_',seq(nas))# ;ENSG;length(ENSG)
  h5write(genenames, output_newh5, "/matrix/features/name")
  h5write(shap, output_newh5, "/matrix/shape")
  h5write(barcode, output_newh5, "/matrix/barcodes")
  h5write(all_tag_keys, output_newh5, "/matrix/features/_all_tag_keys")
  h5write(feature_type, output_newh5, "/matrix/features/feature_type")
  h5write(genome, output_newh5, "/matrix/features/genome")
  h5write(ENSG, output_newh5, "/matrix/features/id")
}

# map_newh5 <- function(input,h5_file,output,newh5_file,dfs){
#   map(dfs,~newh5(input,h5_file,output,newh5_file,.x))
# }
