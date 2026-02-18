library(tidyverse)
library(ggpubr)
library(ggsci)
library(ggstatsplot)
library(ggrepel)
library(ggExtra)
library(palmerpenguins)


Cor_14CRC <- function(path) {
  data_list <- list()
  results <- data.frame(file_gt = character(), file_pre = character(), same_columns = logical())
  
  # 遍历文件
  for (i in seq(0, 5)) {
    file_gt <- file.path(path, paste0("test_", i, "_gt_with_gene.csv"))
    file_pre <- file.path(path, paste0("test_", i, "_pre_with_gene.csv"))
    
    if (file.exists(file_gt) && file.exists(file_pre)) {
      gt <- read.csv(file_gt, check.names = FALSE)
      pre <- read.csv(file_pre, check.names = FALSE)
      
      # 检查列名是否相等
      columns_equal <- all(colnames(gt) == colnames(pre))
      results <- rbind(results, data.frame(file_gt = file_gt, file_pre = file_pre, same_columns = columns_equal))
      
      # 数据预处理
      gt <- gt[, -1]  # 去除第一列
      pre <- pre[, -1]  # 去除第一列
      rownames(gt) <- gt$gene  # 第二列设为行名
      rownames(pre) <- pre$gene
      gt <- gt[, -1]  # 去除"gene"列
      pre <- pre[, -1]
      
      
      # 存储数据框
      data_list[[paste0("test_", i)]] <- list(gt = gt, pre = pre)
    } else {
      warning(paste("Files for test", i, "not found."))
    }
  }
  
  correlation_results <- vector("list", length = 6)
  
  for (i in seq(0, 5)) {
    correlation_result <- list()
    if (!is.null(data_list[[paste0("test_", i)]])) {
      df_gt <- data_list[[paste0("test_", i)]][["gt"]]
      df_pre <- data_list[[paste0("test_", i)]][["pre"]]
      smpl <- unlist(strsplit(colnames(df_gt), "_"))[2] %>% unique()
      for (gene in rownames(df_gt)) {
        # 计算每个样本的相关系数
        cor_values <- cor(unlist(df_gt[gene, ]), unlist(df_pre[gene,]))
        correlation_result[[gene]] = c(Sample_Name = smpl,
        Gene_Name = gene,Correlation = cor_values) # 如果某向量标准差为0(如全0向量),则cor为NA
      }
      # print(str(correlation_result))
      correlation_results[[i+1]] = purrr::map_dfr(correlation_result,~.x)
    }
  }
  final_results <- purrr::map_dfr(correlation_results,~.x)
  colnames(final_results) <- c('Sample_Name', 'Gene_Name','Correlation')
  final_results <- mutate(final_results,Correlation=as.numeric(Correlation))
  return(list(results = results, data_frames = data_list, correlation_results =final_results))
  # return(list(results = results, data_frames = data_list, correlation_results = correlation_results))
}
calculate_correlation <- function(model_data, model_name) {
  map_df(names(model_data$data_frames), function(fold) {
    gt_data <- model_data$data_frames[[fold]]$gt
    pre_data <- model_data$data_frames[[fold]]$pre
    
    rep1_col <- colnames(gt_data) %>% 
      stringr::str_subset(pattern = '.*_Rep1_.*')
    rep2_col <- colnames(gt_data) %>% 
      stringr::str_subset(pattern = '.*_Rep2_.*')
    rep1_sample_names <- stringr::str_extract(rep1_col, ".*_Rep1") %>% unique()
    rep2_sample_names <- stringr::str_extract(rep2_col, ".*_Rep2")%>% unique()
    # 为每个rep分别计算相关性
    rep1_cor <- gt_data[,rep1_col] %>% as.matrix() %>% t() %>% 
      cor(t(as.matrix(pre_data[,rep1_col]))) %>% diag() 
    rep2_cor <- gt_data[,rep2_col] %>% as.matrix() %>% t() %>% 
      cor(t(as.matrix(pre_data[,rep2_col]))) %>% diag() 
    
    # 组合rep1和rep2的数据
    data_frame(
      Gene_Name = c(names(rep1_cor), names(rep2_cor)),
      Correlation = c(rep1_cor, rep2_cor),
      Sample_Name = c(rep(rep1_sample_names , each = length(rep1_cor)), 
                      rep(rep2_sample_names , each = length(rep2_cor))),
      model = rep(model_name, times = length(rep1_cor) + length(rep2_cor))
    )
  })
}

calculate_ssim <- function(model_data, model_name) {
  map_df(names(model_data$data_frames), function(fold) {
    gt_data <- model_data$data_frames[[fold]]$gt
    pre_data <- model_data$data_frames[[fold]]$pre
    
    rep1_col <- colnames(gt_data) %>% 
      stringr::str_subset(pattern = '.*_Rep1_.*')
    rep2_col <- colnames(gt_data) %>% 
      stringr::str_subset(pattern = '.*_Rep2_.*')
    rep1_sample_names <- stringr::str_extract(rep1_col, ".*_Rep1") %>% unique()
    rep2_sample_names <- stringr::str_extract(rep2_col, ".*_Rep2")%>% unique()
      
    safe_ssim <- possibly(SPUTNIK::SSIM, otherwise = NA_real_) #双精度类型的NA值
    rep1_gt_t <- gt_data[,rep1_col] %>% as.matrix() %>% t()  %>% as.data.frame()
    rep1_pre_t <- pre_data[,rep1_col]  %>% as.matrix() %>% t() %>% as.data.frame()
    rep1_ssim <- purrr::map2_dbl(rep1_gt_t,rep1_pre_t,safe_ssim)
    
    rep2_gt_t <- gt_data[,rep2_col] %>% as.matrix() %>% t() %>% as.data.frame()
    rep2_pre_t <- pre_data[,rep2_col]  %>% as.matrix() %>% t() %>% as.data.frame()
    rep2_ssim <- purrr::map2_dbl(rep2_gt_t,rep2_pre_t,safe_ssim)

    # 组合rep1和rep2的数据
    data_frame(
      Gene_Name = c(names(rep1_ssim), names(rep2_ssim)),
      SSIM = c(rep1_ssim, rep2_ssim),
      Sample_Name = c(rep(rep1_sample_names , each = length(rep1_ssim)), 
                      rep(rep2_sample_names , each = length(rep2_ssim))),
      model = rep(model_name, times = length(rep1_ssim) + length(rep2_ssim))
    )
  })
}

Cor_exteral <- function(gt, pre,smpl_model) {
  last_underscore <- regexpr(".*(?=_)", smpl_model, perl = TRUE)
  smpl <- regmatches(smpl_model, last_underscore)
  model <- sub(".*_", "", smpl_model)
  cor_data <- map2_dfr(as.data.frame(t(gt)), as.data.frame(t(pre)), ~cor(.x, .y)) %>%
    pivot_longer(everything(), names_to = "Gene_Name", values_to = "Correlation") %>% 
    mutate(Sample_Name = smpl,model = model)
  return(cor_data)
}

SSIM_14CRC <- function(path) {
  data_list <- list()
  results <- data.frame(file_gt = character(), file_pre = character(), same_columns = logical())
  
  # 遍历文件
  for (i in seq(0, 5)) {
    file_gt <- file.path(path, paste0("test_", i, "_gt_with_gene.csv"))
    file_pre <- file.path(path, paste0("test_", i, "_pre_with_gene.csv"))
    
    if (file.exists(file_gt) && file.exists(file_pre)) {
      gt <- read.csv(file_gt, check.names = FALSE)
      pre <- read.csv(file_pre, check.names = FALSE)
      
      # 检查列名是否相等
      columns_equal <- all(colnames(gt) == colnames(pre))
      results <- rbind(results, data.frame(file_gt = file_gt, file_pre = file_pre, same_columns = columns_equal))
      
      # 数据预处理
      gt <- gt[, -1]  # 去除第一列
      pre <- pre[, -1]  # 去除第一列
      rownames(gt) <- gt$gene  # 第二列设为行名
      rownames(pre) <- pre$gene
      gt <- gt[, -1]  # 去除"gene"列
      pre <- pre[, -1]
      
      
      # 存储数据框
      data_list[[paste0("test_", i)]] <- list(gt = gt, pre = pre)
    } else {
      warning(paste("Files for test", i, "not found."))
    }
  }
  
  correlation_results <- vector("list", length = 6)
  
  for (i in seq(0, 5)) {
    correlation_result <- list()
    if (!is.null(data_list[[paste0("test_", i)]])) {
      df_gt <- data_list[[paste0("test_", i)]][["gt"]]
      df_pre <- data_list[[paste0("test_", i)]][["pre"]]
      smpl <- unlist(strsplit(colnames(df_gt), "_"))[2] %>% unique()
      for (gene in rownames(df_gt)) {
        # 计算每个样本的相关系数
        safe_ssim <- possibly(SPUTNIK::SSIM, otherwise = NA_real_)
        ssim_values <- safe_ssim(unlist(df_gt[gene, ]), unlist(df_pre[gene,]))
        correlation_result[[gene]] = c(Sample_Name = smpl,
                                       Gene_Name = gene,SSIM = ssim_values) # 如果某向量标准差为0(如全0向量),则cor为NA
      }
      # print(str(correlation_result))
      correlation_results[[i+1]] = purrr::map_dfr(correlation_result,~.x)
    }
  }
  final_results <- purrr::map_dfr(correlation_results,~.x)
  colnames(final_results) <- c('Sample_Name', 'Gene_Name','SSIM')
  final_results <- mutate(final_results,SSIM=as.numeric(SSIM))
  return(list(results = results, data_frames = data_list, correlation_results =final_results))
  # return(list(results = results, data_frames = data_list, correlation_results = correlation_results))
}
calculate_correlation <- function(model_data, model_name) {
  map_df(names(model_data$data_frames), function(fold) {
    gt_data <- model_data$data_frames[[fold]]$gt
    pre_data <- model_data$data_frames[[fold]]$pre
    
    rep1_col <- colnames(gt_data) %>% 
      stringr::str_subset(pattern = '.*_Rep1_.*')
    rep2_col <- colnames(gt_data) %>% 
      stringr::str_subset(pattern = '.*_Rep2_.*')
    rep1_sample_names <- stringr::str_extract(rep1_col, ".*_Rep1") %>% unique()
    rep2_sample_names <- stringr::str_extract(rep2_col, ".*_Rep2")%>% unique()
    # 为每个rep分别计算相关性
    rep1_cor <- gt_data[,rep1_col] %>% as.matrix() %>% t() %>% 
      cor(t(as.matrix(pre_data[,rep1_col]))) %>% diag() 
    rep2_cor <- gt_data[,rep2_col] %>% as.matrix() %>% t() %>% 
      cor(t(as.matrix(pre_data[,rep2_col]))) %>% diag() 
    
    # 组合rep1和rep2的数据
    data_frame(
      Gene_Name = c(names(rep1_cor), names(rep2_cor)),
      Correlation = c(rep1_cor, rep2_cor),
      Sample_Name = c(rep(rep1_sample_names , each = length(rep1_cor)), 
                      rep(rep2_sample_names , each = length(rep2_cor))),
      model = rep(model_name, times = length(rep1_cor) + length(rep2_cor))
    )
  })
}

calculate_ssim <- function(model_data, model_name) {
  map_df(names(model_data$data_frames), function(fold) {
    gt_data <- model_data$data_frames[[fold]]$gt
    pre_data <- model_data$data_frames[[fold]]$pre
    
    rep1_col <- colnames(gt_data) %>% 
      stringr::str_subset(pattern = '.*_Rep1_.*')
    rep2_col <- colnames(gt_data) %>% 
      stringr::str_subset(pattern = '.*_Rep2_.*')
    rep1_sample_names <- stringr::str_extract(rep1_col, ".*_Rep1") %>% unique()
    rep2_sample_names <- stringr::str_extract(rep2_col, ".*_Rep2")%>% unique()
    
    safe_ssim <- possibly(SPUTNIK::SSIM, otherwise = NA_real_) #双精度类型的NA值
    rep1_gt_t <- gt_data[,rep1_col] %>% as.matrix() %>% t()  %>% as.data.frame()
    rep1_pre_t <- pre_data[,rep1_col]  %>% as.matrix() %>% t() %>% as.data.frame()
    rep1_ssim <- purrr::map2_dbl(rep1_gt_t,rep1_pre_t,safe_ssim)
    
    rep2_gt_t <- gt_data[,rep2_col] %>% as.matrix() %>% t() %>% as.data.frame()
    rep2_pre_t <- pre_data[,rep2_col]  %>% as.matrix() %>% t() %>% as.data.frame()
    rep2_ssim <- purrr::map2_dbl(rep2_gt_t,rep2_pre_t,safe_ssim)
    
    # 组合rep1和rep2的数据
    data_frame(
      Gene_Name = c(names(rep1_ssim), names(rep2_ssim)),
      SSIM = c(rep1_ssim, rep2_ssim),
      Sample_Name = c(rep(rep1_sample_names , each = length(rep1_ssim)), 
                      rep(rep2_sample_names , each = length(rep2_ssim))),
      model = rep(model_name, times = length(rep1_ssim) + length(rep2_ssim))
    )
  })
}


SSIM_exteral <- function(gt, pre,smpl_model) {
  last_underscore <- regexpr(".*(?=_)", smpl_model, perl = TRUE)
  smpl <- regmatches(smpl_model, last_underscore)
  model <- sub(".*_", "", smpl_model)
  ssim_data <- map2_dfr(as.data.frame(t(gt)), as.data.frame(t(pre)), ~SPUTNIK::SSIM(.x, .y)) %>%
    pivot_longer(everything(), names_to = "Gene_Name", values_to = "SSIM") %>% 
    mutate(Sample_Name = smpl,model = model)
  return(ssim_data)
}

# 绘制Correlations对比散点图
Cor_Scatter <- function(results,model_pairs){
  sample_list<- results %>% pivot_wider(names_from = model, values_from = Correlation) 
  sample_list <- split(sample_list,sample_list$Sample_Name)
  sample_names <- names(sample_list)
  plots_list <- list()
  for (sample in sample_names){
    df = sample_list[[sample]]
    plots_list[[sample]] <- map(model_pairs, ~{
      model_y <- .x[1]
      model_x <- .x[2]
      # 报错最主要的原因是sample_data中某些模型的Correlation值是NA值
      # complete.cases是为了筛选同时不为NA的两列用于绘图
      df[complete.cases(df[,model_x],df[,model_y]),] %>% 
        ggplot(aes_string(x = model_x, y = model_y, label = "Gene_Name"))  +
        geom_point(alpha = 0.6) + geom_text_repel(aes(label=Gene_Name), size=3)+ 
        geom_abline(intercept = 0, slope = 1, color = "red")+
        theme_minimal() + 
        labs(title=paste0("Comparison of Correlations between ",model_y," and ",model_x,' in ',sample))+theme(legend.position="none")
    })
    
  } 
  return(plots_list)
}

