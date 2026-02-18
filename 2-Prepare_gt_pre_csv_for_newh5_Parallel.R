rm(list = ls())
library(tidyverse)
load('RData/gt_pre_TrainValTest_ColName0515.RData')

########################## 将gt,pre等csv 按照样本,模型,gt/pre的名称统一整理 #################

# 1. 14CRC，6折交叉验证，每个样本都可以作为一次内验集

# imap函数可以“炸出” outputs的元素(.x)和元素名称(.y),.x和.y可以进一步作为下游函数的全局变量使用
# ~map(.x$data_frames,function(x))将.x$data_frames作为输入；逐个处理.x$data_frames的每个元素,如test_0
# function(x)中的x指代test_0,通过$gt或$pre,可取出其中的数据框，加以处理
# 在function(x)函数中，.y作为全局变量可直接使用

library(furrr)
plan(multisession, workers = 12)

future_imap(outputs,~future_map(.x$data_frames,function(x){
  ############ x$gt ###################
  split_colnames = colnames(x$gt) %>% str_split("_", simplify = TRUE)
  samples = paste0(split_colnames[,1],'_',split_colnames[,2],'_',split_colnames[,3] )
  unique_sample = unique(samples)
  if (!dir.exists(paste0('newh5/',unique_sample[1]))) {
    # 如果不存在，则创建文件夹
    dir.create(paste0('newh5/',unique_sample[1]),recursive = T)
  }
  if (!dir.exists(paste0('newh5/',unique_sample[2]))) {
    # 如果不存在，则创建文件夹
    dir.create(paste0('newh5/',unique_sample[2]),recursive = T)
  }
  colnames(x$gt) =  split_colnames[,4]
  x$gt[,samples==unique_sample[1]] %>% write.csv(file = paste0('newh5/',unique_sample[1],'/',.y,'_gt.csv'))
  x$gt[,samples==unique_sample[2]] %>% write.csv(file = paste0('newh5/',unique_sample[2],'/',.y,'_gt.csv'))
  ############ x$pre ################ 
  split_colnames = colnames(x$pre) %>% str_split("_", simplify = TRUE)
  samples = paste0(split_colnames[,1],'_',split_colnames[,2],'_',split_colnames[,3] )
  unique_sample = unique(samples)
  colnames(x$pre) =  split_colnames[,4]
  x$pre[,samples==unique_sample[1]] %>% write.csv(file = paste0('newh5/',unique_sample[1],'/',.y,'_pre.csv'))
  x$pre[,samples==unique_sample[2]] %>% write.csv(file = paste0('newh5/',unique_sample[2],'/',.y,'_pre.csv'))
}),.progress=T)

# 2. 外验集 outputs_external 整理
furrr::future_imap(outputs_external,function(x,y){
  sample <- sub("_([^_]*)$", "", y)
  model <- sub(".*_", "", y)
  if (!dir.exists(paste0('newh5/',sample))) {
    # 如果不存在，则创建文件夹
    dir.create(paste0('newh5/',sample),recursive = T)
  }
  x$gt%>% write.csv(file = paste0('newh5/',sample,'/',model,'_gt.csv'))
  x$pre%>% write.csv(file = paste0('newh5/',sample,'/',model,'_pre.csv'))
},.progress=T)

plan(sequential)

