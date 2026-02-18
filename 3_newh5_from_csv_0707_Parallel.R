rm(list = ls())
source('3-newh5_functions_0515.R')

########################## 一、csv to newh5 ###############################
# input=原始.h5所在文件夹
# h5_file=原始.h5文件名
# output = 新文件夹
# newh5_file = 新的.h5文件名
# df = read.csv('pre.csv or gt.csv',row.names = 1,check.names =F) 
csv_num <- list.files(path = 'newh5/34CRC/',pattern = '\\.csv') %>% length()
(output <-list.dirs("newh5",recursive = F) %>% rep(each =csv_num));names(output) <-basename(output)
(input <-  paste0('/Users/a1-6/R/HE2ST_20251001_DXJ/ST_CMS_2023_Roche/',
                  basename(output)))
(h5_file <- list.files(path = input, pattern = "filtered_feature_bc_matrix\\.h5$", full.names = TRUE) %>% basename()) #

(csv_names <- map(unique(output),function(x){
  files <- list.files(path=x,pattern = "\\.csv$", full.names = TRUE)
  file_names <- file_path_sans_ext(basename(files))
  names(files) <- file_names
  return(files)
}))

newh5_file <- csv_names %>% unlist() %>% str_split("/", simplify = TRUE) %>% as.data.frame() %>% 
  mutate(nh5=paste0(V2,'_',sub("\\.csv$", "", V3),'.h5')) %>% .$nh5

library(furrr)
library(tidyverse)
# availableCores() # 可检查本机拥有的线程数
plan(multisession, workers = 12)  # 开启12线程
df <- furrr::future_map(unlist(csv_names),\(x){read.csv(x,row.names = 1,check.names =F)},.progress=T)   # 在6线程下对iris计算每一列的长度
furrr::future_pmap(list(input=input,h5_file = h5_file,output = output,newh5_file = newh5_file ,
          df = df),newh5,.progress=T)
plan(sequential)                 # 回到单线程


########################################### 在newh5/样本名/下创建outs子文件夹，
#################### 并将filtered_feature_bc_matrix.h5结尾的文件和spatial文件夹剪切至newh5/样本名/outs子文件夹下

# 创建 outs 子文件夹
output_outs <- paste0(output, "/outs")
# walk(output_outs, function(out) {
#   if (dir.exists(out)) {
#     unlink(out, recursive = TRUE)
#   }
# })
walk(output_outs, dir.create, showWarnings = FALSE)

# 将文件和文件夹剪切到 outs 子文件夹中，并统一重命名为 filtered_feature_bc_matrix.h5
# (inp <- input[1])
# (outp <- output[1])
# (out <- output_outs[1])
walk2(output,output_outs,function(outp,out) {
  # 获取需要移动的文件和文件夹的路径
  h5_file_to_move <- list.files(path = outp, pattern = ".*filtered_feature_bc_matrix\\.h5$", full.names = TRUE)
  spatial_folder_to_move <- list.dirs(path = outp, full.names = TRUE, recursive = FALSE) %>%
    keep(~ grepl("spatial$", .))
  
  # 重命名并剪切文件
  if (length(h5_file_to_move) > 0) {
    new_h5_file <- file.path(out, "filtered_feature_bc_matrix.h5")
    file.rename(h5_file_to_move, new_h5_file)
  }
  
  # 剪切 spatial 文件夹
  if (length(spatial_folder_to_move) > 0) {
    file.rename(spatial_folder_to_move, file.path(out, basename(spatial_folder_to_move)))
  }
})

