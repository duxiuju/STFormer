# rm(list = ls())
library(SpatialExperiment)
library(tidyverse)
library(jsonlite)
library(tools)
library(furrr)
library(patchwork)
library(RColorBrewer)

library(ggsci)
library(ggstatsplot)
library(ggpubr)
# gm <- expand.grid(models = models,genes = genes418,stringsAsFactors = F)
# plan(multisession, workers = 12)  # 开启12线程
# future_map(paste0('newh5/',samples),\(current_dir){
#   new_outs <- paste0(current_dir,"/outs")
#   dir.create(new_outs, showWarnings = FALSE)
#   old_h5 <- list.files(path = current_dir, pattern = ".*filtered_feature_bc_matrix.h5$", full.names = TRUE)
#   newh5 <- paste0(new_outs,"/filtered_feature_bc_matrix.h5")
#   file.rename(old_h5, newh5)
#   old_spatial <- paste0(current_dir,'/spatial')
#   new_spatial <- paste0(new_outs,'/spatial')
#   file.rename(old_spatial, new_spatial)
# },.progress=T)
# plan(sequential) 
# rowwise(sm) %>%  mutate(spe = list(future_map2(samples,models,safe_construct_spe,.progress = T)))
if (!dir.exists('Figures/SSIM')) {
  # 如果文件夹不存在，则创建文件夹
  dir.create('Figures/SSIM')
  message("Folder created successfully.")
}
if (!dir.exists('Figures/Correlations')) {
  # 如果文件夹不存在，则创建文件夹
  dir.create('Figures/Correlations')
  message("Folder created successfully.")
}

if (!dir.exists('Figures/Heatmap_per_model')) {
  # 如果文件夹不存在，则创建文件夹
  dir.create('Figures/Heatmap_per_model')
  message("Folder created successfully.")
}

if (!dir.exists('Figures/Heatmap_Combined_Models')) {
  # 如果文件夹不存在，则创建文件夹
  dir.create('Figures/Heatmap_Combined_Models')
  message("Folder created successfully.")
}
# sample = '34CRC';model = 'STFormer'
# spe1 <- construct_spe(sample,model)
construct_spe <- function(sample,model){
  current_dir <- paste0('newh5/',sample)
  current_spatial <- paste0(current_dir,'/outs/spatial')
  spatial_pattern <- "tissue_positions.*.csv$"
  spatial_data <- list.files(current_spatial, pattern = spatial_pattern, full.names = T) |> 
    read.csv(header = F) 
  if(spatial_data[1,] |> unlist() |> is.character()){
    spatial_data <- spatial_data[-1,]
  }
  colnames(spatial_data) <- c('Barcode','in_tissue','array_row','array_col',
                              'pxl_row_in_fullres','pxl_col_in_fullres')
  spatial_data <- spatial_data %>% mutate_at(vars(-Barcode),as.numeric)
  col_data <- subset(spatial_data,in_tissue==1,select = -c(pxl_row_in_fullres,pxl_col_in_fullres)) |> 
    mutate(sample_id = model)
  spatial_coords <- subset(spatial_data,in_tissue==1,select = c(pxl_col_in_fullres,pxl_row_in_fullres)) |> 
    as.matrix()
  if(model=='GroundTruth'){ # 原 STFormer_gt 替换为 GroundTruth
    count_name <- paste0(current_dir,'/STFormer_gt.csv')
    count <- read.csv(count_name,check.names = F)
    rownames(count) <- count[,1] %>% gsub(pattern = '-',replacement = '_',.)
    counts <- count[,col_data$Barcode] |> as.matrix()
    row_data <- DataFrame(
      gene_name = rownames(count)
    ) 
    spe <- SpatialExperiment(
      assays = list(counts=counts),
      colData = col_data,
      rowData = row_data,
      spatialCoords = spatial_coords
    )
    url <- paste0(current_spatial,'/tissue_lowres_image.png')
    json_data <- fromJSON(paste0(current_spatial,'/scalefactors_json.json'))
    scalefactor <- json_data[['tissue_lowres_scalef']]
    spe <- addImg(spe, sample_id = model, scaleFactor = scalefactor,image_id = "lowres", 
                  imageSource = url, load = FALSE)
  }else if(model == 'original'){
    spe <- read10xVisium(paste0(current_dir,'/outs/'),type = 'HDF5',
                         data = 'filtered',images = 'lowres')
    rownames(spe) <- rownames(spe) %>% gsub(pattern = '-',replacement = '_',.)
  }else{
    count_name <- paste0(current_dir,'/',model,'_pre.csv')
    count <- read.csv(count_name,check.names = F)
    rownames(count) <- count[,1] %>% gsub(pattern = '-',replacement = '_',.)
    counts <- count[,col_data$Barcode] |> as.matrix()
    row_data <- DataFrame(
      gene_name = rownames(count)
    ) 
    spe <- SpatialExperiment(
      assays = list(counts=counts),
      colData = col_data,
      rowData = row_data,
      spatialCoords = spatial_coords
    )
    url <- paste0(current_spatial,'/tissue_lowres_image.png')
    json_data <- fromJSON(paste0(current_spatial,'/scalefactors_json.json'))
    scalefactor <- json_data[['tissue_lowres_scalef']]
    spe <- addImg(spe, sample_id = model, scaleFactor = scalefactor,image_id = "lowres", 
                  imageSource = url, load = FALSE)
  }
  return(spe)
}
safe_construct_spe <- purrr::possibly(construct_spe,otherwise = 'function error in this loop')

corr_func <- function(count,genes,gt_count){
  c_func <- function(x){
    cor(unlist(as.data.frame(gt_count)[x,]),unlist(as.data.frame(count)[x,colnames(gt_count)]))
  }
  safe_c_func <- possibly(c_func,otherwise = NA)
  names(genes) <- genes
  corrs <- map_dbl(genes,safe_c_func)
  return(corrs)
}

ssim_func <- function(count,genes,gt_count){
  # ssims <- vector(mode = 'list',length = length(genes))
  # ssims <- map(genes,~(ifelse(.x %in% rownames(gt_count),
  #                 SPUTNIK::SSIM(unlist(as.data.frame(gt_count)[.x,]),unlist(as.data.frame(count)[.x,])),
  #                 .x)))
  sfunc <- function(x){
    SPUTNIK::SSIM(unlist(as.data.frame(gt_count)[x,]),unlist(as.data.frame(count)[x,colnames(gt_count)]))
  }
  safe_sfunc <- possibly(sfunc,otherwise = NA)
  names(genes) <- genes
  ssims <- map_dbl(genes,safe_sfunc)
  return(ssims)
  }



# sample = "Control_Rep1"
# model = "STFormer"
# spe <- construct_spe(sample,model)
# library(httpgd)
# hgd()
# ggspavis::plotVisium(spe,annotate = 'EPCAM',point_size = 1.5)

# heatmap_func(spe,'ENO1')

heatmap_func <- function(spe,gene){
  # spe = x$spes[[1]]
  p <- ggspavis::plotVisium(spe,annotate = gene,point_size = 0.6,zoom=F)
  p$layers[[2]]$aes_params$stroke <- 0 #边框宽度设为0
  fill_limits <- range(p$data[[gene]], na.rm = TRUE) # 获取gene值的范围
  p <- p + scale_fill_gradientn(
    colours = rev(colorRampPalette(brewer.pal(11,'RdYlBu'))(256)),  # color Pallete的11个颜色变为256色后，反转(rev)作为legend 色条
    limits = fill_limits, # 色条的数值范围与gene值保持一致
    labels = scales::number_format(accuracy = 0.01)
  )+ theme(
    legend.key.size = unit(0.4, "cm"),  # 调整色条高度
    legend.key.width = unit(0.3, "cm"),  # 调整色条宽度
    legend.title = element_text(size = 6), # 调整色条标题文字大小
    legend.text = element_text(size = 4),  # 调整色条刻度的文字大小
    legend.spacing.x = unit(-0.1, 'cm')  # 减少色条之间的水平间距
  )
  return(p)
}
safe_heatmap <- possibly(heatmap_func,otherwise = NULL)


#################  ！！！！！！！！！！！！！！！！！！！！plot Pathology Annotation with spe ！！！！！！！！！！！！！！！！！！ ##############################
# extracted_data <- fromJSON("Pathology_Annotation_JSON/SN84_A120838_Rep2/SN84_A120838_Rep2_extracted_data.json")
# glimpse(extracted_data)
# str(extracted_data)

# library(httpgd)
# hgd()
# extracted_data$scaled_coordinates[[1]] %>% matrix(nrow = dim(.)[2], ncol = dim(.)[3]) %>%
#   as.data.frame()  %>% rename(X = V1,Y = V2) %>%
#   ggplot(aes(x = X, y = Y)) +geom_path() + theme_minimal() +labs(title = "Scaled Coordinates Path Plot", x = "X", y = "Y")

# anno <- mutate(extracted_data,scaled_coordinates = map(scaled_coordinates,\(x){
#   matrix(x,nrow = dim(x)[2], ncol = dim(x)[3]) %>% as.data.frame()%>% rename(X = V1,Y = V2)}),
#   unique_id = paste(name, row_number(), sep = "_")) %>% unnest(scaled_coordinates)  %>%
#   ggplot(aes(x = X, y = Y,color = unique_id, group = unique_id))  +
#   geom_path() + theme_minimal() +labs(title = "Scaled Coordinates Path Plot", x = "X", y = "Y")
# spe <- subset(spes_df,samples=="SN048_A121573_Rep1",select=spes)[1,1][[1]][[1]];spe
# gene='EPCAM'
# heat <- heatmap_func(spe,gene)

# combined_plot <- heat +
#   geom_path(data = mutate(anno$data, X = X, Y = -Y + 600),
#             aes(x = X, y = Y, group = unique_id, color = factor(name)), inherit.aes = FALSE,linewidth = 2) +
#   scale_color_manual(
#     values = c("AggregateImmuneCell" = "cyan",
#                "Tumor" = "red",
#                "NormalEpithelium" = "green"),
#     drop = TRUE  # 确保去掉不在指定范围内的颜色，如transparent
#   ) +
#   guides(
#     fill = guide_colorbar(order = 1),
#     color = guide_legend(order = 2, title = "Region")
#   ) +
#   labs(color = "Region")  # 更改图例标题为 "Region"

# combined_plot
#####################################################################################################################

spe_pipe <- function(spe){
  # spe <- spes_df$spes[[2]];spe
  set.seed(123)
  genes <- sample(rownames(spe),50)
  spe <- scater::runPCA(spe,subset_row = genes,exprs_values = "counts")
  spe <- scater::runUMAP(spe, dimred = "PCA")
  colnames(reducedDim(spe, "UMAP")) <- paste0("UMAP", 1:2)
  set.seed(123)
  k <- 10
  g <- scran::buildSNNGraph(spe, k = k, use.dimred = "PCA")
  g_walk <- igraph::cluster_walktrap(g)
  clus <- g_walk$membership # clustering is finished in this line
  colLabels(spe) <- factor(clus)
#   ggspavis::plotDimRed(spe, plot_type = "UMAP", 
#            annotate = "label")
# # plot clusters in spatial x-y coordinates
#   ggspavis::plotSpots(spe, annotate = "label", point_size = 1.5)
  return(spe)
  ################ 六、PathoAnnotation, then ARI between Clusters and PathoAnno
}
safe_spe_pipe <- possibly(spe_pipe,otherwise = NULL)

# custom_theme <- function(rl = 1.1) {
#   stopifnot(rl > 0)
#   ggplot2::theme_minimal() + ggplot2::theme(
#     panel.border = element_rect(colour = "black", fill = NA),
#     panel.grid = element_blank(),
#     axis.title = element_blank(),
#     axis.text = element_text(size = rel(rl)),
#     plot.title = element_text(size = rel(rl) * 1.2),
#     strip.background = element_rect(fill = NA, colour = "black"),
#     strip.text = element_text(size = rel(rl)),
#     legend.text = element_text(size = rel(rl)),
#     legend.title = element_text(size = rel(rl), face = "italic"),
#     legend.position = "bottom"
#   )
# }
# cluster_plot <- function(spe){
#   spatialCoords(spe) |>
#   as.data.frame() |>
#   cbind("Region" = spe$label) |>
#   ggplot(aes(pxl_col_in_fullres, pxl_row_in_fullres, colour = Region)) + geom_point() +
#   scale_colour_brewer(palette = "Paired", 
#   guide = guide_legend(override.aes = list(shape = 15, size = 2))) +
#   labs(title = "Region") # +custom_theme()
# }
# # map(~cluster_plot(.x)) 

# spe = spes_df$spes[[31]];smpl = spes_df$samples[[31]]
spe_Patho <- function(spe,smpl){
  patho_anno <- readRDS(paste0('PathologyAnno/14CRC_Pathology_Annotations/Pathologist_Annotations_',smpl,'.rds'))
  rownames(patho_anno) <- patho_anno$spot_id
  patho_anno <- patho_anno[rownames(colData(spe)),]
  colData(spe)$Pathologist_Annotations<- as.factor(patho_anno$Pathologist_Annotations)
  return(spe)
  # plot ground truth labels in spatial coordinates
  # library(httpgd)
  # hgd()
  # ggspavis::plotSpots(spe, annotate = "patho", point_size = 1.5,
  #           pal = "libd_layer_colors")
}

######################## Pathology Percent Functions ###############################
# 包含多分类函数
classify_column <- function(column) {
  # 找到最大值对应的行名
  max_row <- rownames(pathoAnno)[which.max(column)]
  
  # 找到其他大于0.3的行名
  additional_rows <- rownames(pathoAnno)[which(column > 0.3 & rownames(pathoAnno) != max_row & rownames(pathoAnno) != 'None')]
  
  # 组合分类标签
  if (max_row == "None" && length(additional_rows) > 0) {
    label <- paste(additional_rows, collapse = "_")
  } else {
    label <- max_row
    if (length(additional_rows) > 0) {
      label <- paste(label, paste(additional_rows, collapse = "_"), sep = "_")
    }
  }
  
  return(label)
}

# 仅计算max类函数
max_classify_column <- function(column) {
  # 找到最大值对应的行名
  max_row <- rownames(pathoAnno)[which.max(column)]
  
  # 找到其他大于0.1的行名
  additional_rows <- rownames(pathoAnno)[which(column > 0.3 & rownames(pathoAnno) != max_row & rownames(pathoAnno) != 'None')]
  
  # 组合分类标签
  label <- max_row
  # if (length(additional_rows) > 0) {
  #   label <- paste(label, paste(additional_rows, collapse = "_"), sep = "_")
  # }
  
  return(label)
}

################################ Unit Test ##############################

# sample = "Control_Rep1"
# model = "HIPT"
# spe <- construct_spe(sample,model)
# ggspavis::plotVisium(spe,annotate = 'EPCAM',point_size = 1.5)
# 
# 
# spes <- list()
# spes[[sample]][[model]] <- safe_construct_spe(sample,model)
# ggspavis::plotVisium(spes[[sample]][[model]],annotate = 'EPCAM',point_size = 1.5)
# 
# 
