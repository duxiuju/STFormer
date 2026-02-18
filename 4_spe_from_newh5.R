rm(list = ls())
source('4-spe_functions_0707.R') 
options(future.globals.maxSize = 10000*1024^2)  # 将最大大小设置为102000 MiB
########################################## 1. Prepare dataframe for calculation and plot ################################## 
samples <- c('34CRC','Control_Rep1','Control_Rep2','Human',
             'SN048_A121573_Rep1','SN048_A416371_Rep1','SN123_A595688_Rep1',
             'SN048_A121573_Rep2','SN048_A416371_Rep2','SN124_A595688_Rep2',
             'SN123_A798015_Rep1','SN123_A938797_Rep1','SN84_A120838_Rep1',
             'SN124_A798015_Rep2','SN124_A938797_Rep2','SN84_A120838_Rep2')
samples_14CRC <- samples[-seq(4)]
Patho_path <- 'Figures/Pathology_Annotations'
models <- c('Hist2ST','HisToGene','STFormer','STNet','GroundTruth') # STFormer_gt替换为Ground_Truth 加上'original' 会导致ssim计算报错：original  <SptlExpr[,3138]> <TENxMtrx[,3138]>
genes418 <- read.csv('OtherData/418Genes.csv',check.names = F)  %>% .$GeneSymbol %>% 
  gsub(pattern = '-',replacement = '_',.)
plan(multisession, workers = 12)
# load('RData/spes_df.RData') # 加载完数据后可以直接跳至80行开始绘图
#                             group_by(spes_df,samples) %>% group_map(\(x,y)
for (sample in samples) {
  source_path <- file.path("ST_CMS_2023_Roche", sample, "outs")
  destination_path <- file.path("newh5", sample)
  
  # 创建目标文件夹如果它不存在
  if (!dir.exists(destination_path)) {
    dir.create(destination_path, recursive = TRUE)
  }
  
  # 复制文件夹内容
  file.copy(source_path, destination_path, recursive = TRUE)
}

# 'HLA_A'%in% genes418

sm <- expand.grid(models = models,samples = samples,stringsAsFactors = F) |> 
  dplyr::select(samples,models)
spes_df <- expand.grid(models = models,samples = samples,stringsAsFactors = F) %>% as_tibble() 
spes_df <- mutate(spes_df,spes=furrr::future_map2(samples,models,construct_spe,.progress = T),
                  spe_counts=furrr::future_map(spes,~assay(.x,i = 'counts'),.progress = T)) 
spes_df <- mutate(spes_df,genes = rep(list(genes418),nrow(spes_df)))                 
spes_df <- spes_df %>% group_by(samples)  %>% 
  mutate(gt_counts = spe_counts[models=='GroundTruth']) # STFormer_gt替换为GroundTruth

spes_df # View the DataFrame for next calculation and plot 
# plan(sequential)


############################################ 2. generating SSIM,Correlations, and Heatmap ########################################

# 2.1 SSIM 
# plan(multisession, workers = 124)
# group_by(spes_df,samples) %>% group_map(\(x,y){
#   furrr::future_pmap_
# dfc(list(x$spe_counts,x$gt_counts,x$genes),ssim_func) %>% 
#   `colnames<-`(x$models)  %>% `rownames<-`(x$genes[[1]]) %>% 
#   write.csv(file=paste0('SSIM/',unique(y$samples),'_ssim.csv'))
# })
# plan(sequential)


# plan(multisession, workers = 120)
spes_df <- spes_df %>% group_by(samples) %>% 
  mutate(ssim  = furrr::future_map2(spe_counts,genes,
                                    ssim_func,spe_counts[[5]],.progress=T))
# plan(sequential)

# 2.2 Correlations 
# plan(multisession, workers = 124)
# group_by(spes_df,samples) %>% group_map(\(x,y){
#   furrr::future_pmap_dfc(list(x$spe_counts,x$gt_counts,x$genes),corr_func) %>% 
#   `colnames<-`(x$models)  %>% `rownames<-`(x$genes[[1]]) %>% 
#   write.csv(file=paste0('Correlations/',unique(y$samples),'_corr.csv'))
# })
# plan(sequential)
# plan(multisession, workers = 120)
spes_df <- mutate(spes_df,corr = furrr::future_map2(spe_counts,genes,corr_func,
                                                    spe_counts[[5]],.progress=T))
# plan(sequential)

spes_df 
utils::object.size(spes_df) |> print(units='auto',standard = 'SI')

######################################## 2.3 Heatmap ######################################
############## unit test ######################
# library(httpgd)
# hgd()
# ggspavis::plotVisium(spes_df$spes[[5]],annotate = 'ENO1',point_size = 1.5)

############## unit test #########################
# plan(multisession, workers = 120)
# map_chr(genes418,~paste0('Heatmap_Combined_Models/',.x))
# gene <- 'HLA-A'
# # x = spes_df[1:6,]
# numbers <- 1:114
# sequences_list <- split(numbers, ceiling(seq_along(numbers)/6))
# x = map(sequences_list,~spes_df[.x,])
# x = subset(spes_df,samples=='SN123_A798015_Rep1')
# y = x$samples
# unique(y)
# furrr::future_map(x,~(.x$spes %>% map(~safe_heatmap(.x,gene))%>% purrr::compact()  %>% 
#     reduce(`+`)))
  # group_nest(spes_df,samples) %>% .$data %>%
# save(spes_df,file = 'RData/spes_df.RData')

# rownames(spes_df$spes[[1]])
# inter_02 <- c('ELF3','EPCAM','ATP5F1E','CD24','CD9','CHCHD2','CLDN3','CLDN4','COX4I1','COX7A2','COX7C',
#                'HNRNPA2B1','KLF5','KRT8','LGALS3','LGALS4','MUC13','NME2',
#                'PERP','S100A10','SRSF2','SRSF3','TUFM','UBA52','UQCRB')
# write.csv(as.data.frame(inter_02),file='25genes.csv')


spes_df$models <- factor(spes_df$models,levels = c("GroundTruth","STFormer","STNet","Hist2ST","HisToGene"))
spes_df %>% group_by(samples)  %>% arrange(models) %>% group_map(\(x,y){
  map(genes418,\(gene){ # inter_02
    save_path = paste0('Figures/Heatmap_Combined_Models/',gene)
    print(save_path)
    if (!dir.exists(save_path)){dir.create(save_path,recursive = T)}
    if(gene %in% rownames(x$spes[[1]])){
      print(paste0(gene,'_',y$samples))
      # map(x$spes,~safe_heatmap(.x,gene)) %>% purrr::compact()  %>% 
      map2(x$spes,x$corr,\(w,z){
        # print(w)
        # print(z)
        # safe_heatmap(w,gene)+labs(caption = paste0('Correlations = ', round(z[[gene]], 3))) +
        #                      theme(plot.margin = unit(c(1, 1, 1.5, 1), "lines"),
        #                      plot.caption = element_text(hjust = 0.5, size = 5))
        safe_heatmap(w, gene) + labs(caption = paste0('Correlations = ', round(z[[gene]], 3))) +
                                theme(plot.caption = element_text(size = 5,hjust = 0.5,vjust =2),
                                      plot.title = element_text(size = 7))
                             
      }) %>% purrr::compact()  %>% 
        purrr::reduce(`+`) %>% ggsave(file = paste0(save_path,'/',y$samples,'.tiff'))
  }})
})

############################### 2.4 spe pipeline and then plot UMAP with cluster,then PathologyAnno
save_pa <- 'Figures/UmapCluster'
if (!dir.exists(save_pa)){
  dir.create(save_pa,recursive = T)
}

######################## 2.4.1 spe_pipe then plot UMAP and cluster 
# plan(multisession, workers = 120)
spes_df <- spes_df %>% mutate(spes = furrr::future_map(spes,safe_spe_pipe))
# 
# 
# spes_df %>% group_by(samples) %>% 
#   group_map(\(x,y){
#     max_colors <- x$spes %>% map_dbl(\(x){nlevels(x$label)}) %>%max()
#     customizes_pal <- colorRampPalette(brewer.pal(10,'Paired'))(max_colors)
#     x$spes %>% purrr::compact() %>% map(~ggspavis::plotDimRed(.x,
#     plot_type = "UMAP", annotate = "label")+
#     scale_colour_manual(values = customizes_pal)) %>% 
#     purrr::reduce(`+`) %>% ggsave(file = paste0(save_pa,'/',y$samples,'_UmapCluster.tiff'))
#   }) 
# 
# spes_df %>% group_by(samples) %>% 
#   group_map(\(x,y){
#     max_colors <- x$spes %>% map_dbl(\(x){nlevels(x$label)}) %>%max()
#     customizes_pal <- colorRampPalette(brewer.pal(12,'Paired'))(max_colors)
#     x$spes %>% purrr::compact() %>% map(~ggspavis::plotSpots(.x, sample_id = 'sample_id',annotate = "label", 
#     point_size = 0.15)+scale_colour_manual(values = customizes_pal))  %>%
#     # "Paired"3
#     purrr::reduce(`+`) %>% ggsave(file = paste0(save_pa,'/',y$samples,'_SpatialCluster.tiff'))
#   }) 

############################# 2.4.2 14CRC(Charting) Original Pathology Annotation Visualization 
# plan(multisession, workers = 120)

if (!dir.exists(Patho_path)){
  dir.create(Patho_path,recursive = T)
}

# subset(spes_df,samples%in%samples_14CRC) %>% group_by(samples) %>% 
#   group_map(\(x,y){
    
#     map(x$spes,~spe_Patho(.x,y$samples)) %>% map(~ggspavis::plotSpots(.x, annotate = "patho", 
#     point_size = 0.5)+scale_colour_manual(values = customizes_pal))  %>%
#     # "Paired"
#     reduce(`+`) %>% ggsave(file = paste0(Patho_path,'/',y$samples,'_PathoAnno.tiff'))
#   })
# Unit Test 
# x = subset(spes_df,samples == 'SN048_A121573_Rep1')
# y = data.frame(samples='SN048_A121573_Rep1');y

subset(spes_df, samples %in% samples_14CRC) %>%
  group_by(samples) %>%
  group_map(\(x, y){
    
    # 1) 生成 SPE 列表
    new_spes <- x$spes %>%
      purrr::compact() %>%
      map(~ spe_Patho(.x, y$samples))
    
    # === 新增：如果有模型名，给 new_spes 命名，便于按指定顺序重排 ===
    if ("models" %in% names(x)) {
      names(new_spes) <- as.character(x$models)
      desired_models <- c("GroundTruth","STFormer","STNet","HisToGene","Hist2ST")
      # 仅按你给的顺序保留存在的模型
      keep <- intersect(desired_models, names(new_spes))
      new_spes <- new_spes[keep]
    }
    
    # 2) 统计两种注释列的最大水平数（避免某图颜色不够）
    nlev <- function(spe, col){
      if (col %in% names(colData(spe))) {
        nlevels(as.factor(colData(spe)[[col]]))
      } else 0
    }
    max_levels <- max(
      map_int(new_spes, ~ nlev(.x, "label")),
      map_int(new_spes, ~ nlev(.x, "Pathologist_Annotations"))
    )
    seed <- brewer.pal(min(12, max(3, max_levels)), "Paired")
    pal  <- grDevices::colorRampPalette(seed)(max_levels)
    
    # 3) 画图
    cluster_plots <- map(
      new_spes,
      ~ ggspavis::plotVisium(.x,
                             legend_position = "right",
                             annotate = "label",
                             point_size = 0.4,
                             pal = pal)
    )
    # ---- 病理学家标注图（修改处1：显式标题）----
    new_spes[[1]]$sample_id = 'Pathology Annotation' 
    pathoAnno_plot <- ggspavis::plotVisium(
      new_spes[[1]],
      annotate = "Pathologist_Annotations",
      legend_position = "right",
      point_size = 0.4,
      pal = pal
    ) # ← 新增标题
    
    # ---- 修改处2：按指定顺序把病理学家图放到最前面 ----
    all_plots <- c(list(pathoAnno_plot), cluster_plots)
    
    p <- (purrr::reduce(all_plots, `+`) +
            plot_layout(ncol = 3, nrow = 2))
    
    ggsave(
      filename = file.path(Patho_path, paste0(y$samples, "_PathoAnno.tiff")),
      plot = p, width = 12, height = 12, dpi = 300
    )
  })

plan(sequential)

#################  ！！！！！！！！！！！！！！！！！！！！plot Pathology Annotation with spe ！！！！！！！！！！！！！！！！！！ ##############################

# 样本列表（排除指定样本）
sample_names <- list.dirs('PathologyAnno/External_Pathology_Annotation_JSON', full.names = FALSE) %>%
  basename() %>% .[-1]
sample_names <- sample_names[sample_names != "SN123_A551763_Rep1"]

# 基因列表
genes <- c('KRT8','EPCAM','CD9','VIM','KLF5','HSP90AA1')

# 输出目录
PathologyCombinedPlots <- 'Figures/PathologyCombinedPlots'
PathologyPercentPlots  <- 'Figures/PathologyPercentPlots'
dir.create(PathologyCombinedPlots, showWarnings = FALSE, recursive = TRUE)
dir.create(PathologyPercentPlots,  showWarnings = FALSE, recursive = TRUE)

# 线宽 & 分辨率
linewidth <- 0.5
DPI <- 120

# 统计比较
compare_list <- list(c("NNM", "Stroma"))  # 你要求统一使用 compare_list

# 统一的“病理类别”顺序 & 颜色（确保 stformer_box 与 gt_box 顺序与颜色一致）
patho_levels <- c('Tumor','Stroma','AIC','Muscle','Necrosis','RBC','NNM','NE','None')
my_colors_vec <- setNames(colorRampPalette(brewer.pal(12,'Set3'))(length(patho_levels)), patho_levels)

# up_Y 与 X 偏移映射（严格按你的说明）
upY_map <- c(
  "34CRC"                 = 80,
  "Control_Rep1"          = 20,
  "Human"                 = 20,
  "SN048_A121573_Rep1"    = 170,
  "SN048_A416371_Rep1"    = 120,
  "SN123_A798015_Rep1"    = 115,
  "SN124_A938797_Rep2"    = 45,
  "SN84_A120838_Rep2"     = 125,
  "SN123_A595688_Rep1"    = 210  # 同时需要 X 偏移 +10
)
xShift_map <- c("SN123_A595688_Rep1" = 10)  # 仅此样本需要 X+10，其余默认为 0

# 辅助：根据 sample_nm 推断使用哪个 key（你给的映射有数据集名，也有具体样本名）
resolve_upY <- function(sample_nm){
  # 完全匹配优先；否则尝试包含关键字
  if (sample_nm %in% names(upY_map)) return(upY_map[[sample_nm]])
  keys <- names(upY_map)
  hit <- keys[str_detect(sample_nm, fixed(keys))]
  if (length(hit) > 0) return(upY_map[[hit[1]]])
  45  # 默认值（你在注释中给到的 SN124_A938797_Rep2~45）
}
resolve_xShift <- function(sample_nm){
  if (sample_nm %in% names(xShift_map)) return(xShift_map[[sample_nm]])
  keys <- names(xShift_map)
  hit <- keys[str_detect(sample_nm, fixed(keys))]
  if (length(hit) > 0) return(xShift_map[[hit[1]]])
  0
}

## ================== 主循环：sample_nm × gene ==================
for (sample_nm in sample_names) {
  
  # 读取外部病理注释 JSON → scaled coordinates
  json_fp <- file.path("PathologyAnno/External_Pathology_Annotation_JSON", sample_nm,
                       paste0(sample_nm, "_extracted_data.json"))
  if (!file.exists(json_fp)) {
    warning("缺少 JSON：", json_fp, "，跳过该样本。")
    next
  }
  extracted_data <- fromJSON(json_fp)
  
  scl_coo <- extracted_data$scaled_coordinates
  names(scl_coo) <- extracted_data$name
  scl_coo_list <- scl_coo %>%
    list_flatten(name_spec = "{outer}") %>%
    list_flatten(name_spec = "{outer}") %>%
    map_if(.p = ~ length(dim(.x)) == 3,
           .f = \(x) matrix(x, nrow = dim(x)[2], ncol = dim(x)[3]))
  
  scl_coo_df <- map2_dfr(
    scl_coo_list, seq_along(scl_coo_list),
    \(x, y) { as.data.frame(x) %>% mutate(region = str_c('region_', y)) },
    .id = "Pathology_Annotation"
  ) %>% dplyr::rename(X = V1, Y = V2)
  
  # annotation gg 对象（稍后只取其 data）
  anno <- scl_coo_df %>%
    ggplot(aes(x = X, y = Y, color = Pathology_Annotation, group = region)) +
    geom_path() + theme_minimal() +
    labs(title = "Scaled Coordinates Path Plot", x = "X", y = "Y")
  
  spe_stformer <- subset(spes_df, samples == sample_nm, select = spes)[3,1][[1]][[1]]  # STFormer
  spe_gt       <- subset(spes_df, samples == sample_nm, select = spes)[5,1][[1]][[1]]  # GroundTruth
  spe_he       <- spe_stformer                                                           # 仅用于 HE 叠加
  # 仅这份更名为“Pathology Annotation”
  spe_he$sample_id <- 'Pathology Annotation'
  
  up_Y   <- resolve_upY(sample_nm)
  xShift <- resolve_xShift(sample_nm)
  
  for (gene in genes) {
    
    # ---------- CombinedPlots ----------
    heat_st  <- heatmap_func(spe_stformer, gene)      # ← 用 spe_stformer
    add_Y    <- max(round(heat_st$data$pxl_row_in_fullres), na.rm = TRUE) + up_Y
    heat_gt  <- heatmap_func(spe_gt, gene)            # ← 用 spe_gt
    heat_he  <- heatmap_func(spe_he, gene)            # ← 仅 HE 用 spe_he
    
    # STFormer 叠加注释（第三幅应当是 STFormer）
    combined_plot <- heat_st +
      geom_path(
        data = dplyr::mutate(anno$data, X = X + xShift, Y = -Y + add_Y),
        aes(x = X, y = Y, group = region, color = factor(Pathology_Annotation)),
        inherit.aes = FALSE, linewidth = linewidth
      ) +
      scale_color_manual(
        values = c("AggregatedImmuneCell"="#00FFFF","Tumor"="#FF0000","Stroma"="#00FF00",
                   "Necrosis"="#D8DAD9","Non-neoplastic_Mucosa"="#FBDE28","Normal_Epithelium"="#FBDE28",
                   "RBC"="#ED740C","Muscle"="#0000FF","Not Defined"="#262626"),
        drop = TRUE
      ) +
      guides(fill = guide_colorbar(order = 1),
             color = guide_legend(order = 2, title = "Region")) +
      labs(color = "Region")
    
    # HE 覆盖（第一幅）
    HE <- heat_he +
      geom_path(
        data = dplyr::mutate(anno$data, X = X + xShift, Y = -Y + add_Y),
        aes(x = X, y = Y, group = region, color = factor(Pathology_Annotation)),
        inherit.aes = FALSE, linewidth = linewidth
      ) +
      scale_color_manual(
        values = c("AggregatedImmuneCell"="#00FFFF","Tumor"="#FF0000","Stroma"="#00FF00",
                   "Necrosis"="#D8DAD9","Non-neoplastic_Mucosa"="#FBDE28","Normal_Epithelium"="#FBDE28",
                   "RBC"="#ED740C","Muscle"="#0000FF","Not Defined"="#262626"),
        drop = TRUE
      ) +
      guides(fill = guide_colorbar(order = 1),
             color = guide_legend(order = 2, title = "Region")) +
      labs(color = "Region")
    HE$layers <- HE$layers[sapply(HE$layers, function(layer) !inherits(layer$geom, "GeomPoint"))]
    
    # GroundTruth 叠加注释（第二幅）
    combined_plot_gt <- heat_gt +
      geom_path(
        data = dplyr::mutate(anno$data, X = X + xShift, Y = -Y + add_Y),
        aes(x = X, y = Y, group = region, color = factor(Pathology_Annotation)),
        inherit.aes = FALSE, linewidth = linewidth
      ) +
      scale_color_manual(
        values = c("AggregatedImmuneCell"="#00FFFF","Tumor"="#FF0000","Stroma"="#00FF00",
                   "Necrosis"="#D8DAD9","Non-neoplastic_Mucosa"="#FBDE28","Normal_Epithelium"="#FBDE28",
                   "RBC"="#ED740C","Muscle"="#0000FF","Not Defined"="#262626"),
        drop = TRUE
      ) +
      guides(fill = guide_colorbar(order = 1),
             color = guide_legend(order = 2, title = "Region")) +
      labs(color = "Region")
    
    # 保存三联图（顺序：Pathology Annotation / GroundTruth / STFormer）
    tiff(filename = file.path(PathologyCombinedPlots, paste0(sample_nm, "_", gene, ".tif")),
         width = 1920, height = 1017, units = "px", res = DPI)
    print(HE / combined_plot_gt / combined_plot)
    dev.off()
    
    # ---------- Percent + BoxPlot ----------
    pathoAnno_path <- file.path('PathologyAnno/SpotClass_from_Annotation',
                                paste0(sample_nm, '_percentage.csv'))
    if (!file.exists(pathoAnno_path)) next
    pathoAnno <- read.csv(pathoAnno_path, check.names = FALSE)
    pathoAnno <- pathoAnno[, -1, drop = FALSE]
    rownames(pathoAnno) <- c('Tumor','Stroma','AIC','Muscle','Necrosis','RBC','NNM','NE','None')
    
    # **严格一致性校验**
    if (!all(colnames(pathoAnno) == rownames(colData(spe_stformer)))) next
    
    max_classify_column <- function(v) { rownames(pathoAnno)[which.max(v)] }
    classification_labels <- apply(pathoAnno, 2, max_classify_column)
    
    # 仅写回 STFormer 的 spe（不改 sample_id）
    colData(spe_stformer)$patho <- factor(classification_labels,
                                          levels = c('Tumor','Stroma','AIC','Muscle','Necrosis','RBC','NNM','NE','None'))
    
    df_st <- data.frame(
      gene_value = assay(spe_stformer)[gene, ],
      patho      = colData(spe_stformer)$patho
    )
    df_gt <- data.frame(
      gene_value = assay(spe_gt)[gene, ],
      patho      = colData(spe_stformer)$patho  # 和 ST 一致
    )
    
    # 统一 palette 与顺序
    patho_levels <- levels(colData(spe_stformer)$patho)
    my_colors_vec <- setNames(colorRampPalette(RColorBrewer::brewer.pal(12,'Set3'))(length(patho_levels)),
                              patho_levels)
    
    stformer_box <- df_st %>%
      mutate(patho = factor(patho, levels = patho_levels)) %>%
      ggpubr::ggboxplot(
        x = "patho", y = "gene_value",
        color = "patho", palette = my_colors_vec,
        add = "jitter", outlier.shape = NA
      ) +
      labs(y = gene) +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      ggpubr::stat_compare_means(comparisons = compare_list)
    
    gt_box <- df_gt %>%
      mutate(patho = factor(patho, levels = patho_levels)) %>%
      ggpubr::ggboxplot(
        x = "patho", y = "gene_value",
        color = "patho", palette = my_colors_vec,
        add = "jitter", outlier.shape = NA
      ) +
      labs(y = gene) +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      ggpubr::stat_compare_means(comparisons = compare_list)
    
    # **Percent 图：左 GT（地图） + 右 STFormer（地图）**；第二行是对应的 box
    layout <- "
    12
    34
    "
    tiff(filename = file.path(PathologyPercentPlots, paste0(sample_nm, "_", gene, ".tif")),
         width = 1920, height = 1017, units = "px", res = DPI)
    print(combined_plot_gt + stformer_box + combined_plot + gt_box + plot_layout(design = layout))
    # └─ 上行的 1=GT地图, 2=ST box, 3=ST地图, 4=GT box（若你只想两幅地图，改成 print(combined_plot_gt | combined_plot)）
    dev.off()
  } # end for gene
} # end for sample_nm

################## soft dice ###########
# 生成 One-Hot 编码的病理注释
data <- df_gt
patho_onehot <- model.matrix(~ 0 +patho, data)
# 使用正则表达式去掉前缀 "patho"
colnames(patho_onehot) <- sub("^patho", "", colnames(patho_onehot))

# 计算 Soft-Dice Coefficient
soft_dice_coefficient <- function(pred, true_onehot) {
  dice_scores <- apply(true_onehot, 2, function(true_class) {
    2 * sum(pred * true_class) / (sum(pred^2) + sum(true_class^2) + 1e-6) # 避免分母为0
  })
  return(dice_scores)
}

# 计算 Soft-Dice (对每个病理类别计算)
soft_dice <- soft_dice_coefficient(data$gene_value, patho_onehot)
soft_dice

#############################

# library(pROC)
# roc_obj <- roc(df$patho, df$gene_value,levels = c('Stroma','Tumor'))  # 创建 ROC 对象
# plot(roc_obj, main="ROC Curve", col="blue", print.auc=T,xlim=c(1, 0),
#      print.thres=T)
# 
# roc_obj2 <- roc(df$patho, df$gene_value,smooth=T) 
# cutOffPoint <- coords(roc_obj2, "best")
# cutOffPointText <- paste0("(",round(cutOffPoint[[1]],2),',',round(cutOffPoint[[2]],2),")")
# ggroc(roc_obj2,
#       color="red",
#       size=1,
#       legacy.axes = F # F 横坐标1~0 为specificity；
#                       # T 横坐标0~1 为1-specificity
# )+
#   theme_classic()+
#   geom_segment(aes(x = 1, y = 0, xend = 0, yend = 1),        # 绘制对角线
#                colour='grey', 
#                linetype = 'dotdash'
#   ) +
#   geom_point(aes(x = cutOffPoint[[1]],y = cutOffPoint[[2]]))+ # 绘制临界点/阈值
#   geom_text(aes(x = cutOffPoint[[1]],y = cutOffPoint[[2]],label=cutOffPointText),vjust=-1)+
#   annotate("text",x=0.75,y=0.25,label=paste("AUC = ", round(roc_obj$auc,3)))# 添加临界点/阈值文字标签
# 
# 
# 
