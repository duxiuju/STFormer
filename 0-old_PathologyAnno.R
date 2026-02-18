######### PathologyCombinedPlots Unit Test ############

sample_nm <- 'SN124_A938797_Rep2' # 'SN84_A120838_Rep2'
json <- paste0("PathologyAnno/External_Pathology_Annotation_JSON/",sample_nm,'/',
               sample_nm,'_extracted_data.json')
gene='EPCAM'  #  KRT8 EPCAM CD9 VIM KLF5 HSP90AA1
extracted_data <- fromJSON(json) 

# A121573_extracted_data.json
# glimpse(extracted_data)
# glimpse(extracted_data$scaled_coordinates)
# str(extracted_data$scaled_coordinates)
# list_flatten(extracted_data $scaled_coordinates) %>% str()
# class(extracted_data)

scl_coo <- extracted_data$scaled_coordinates;class(scl_coo);length(scl_coo)
str(scl_coo)
names(scl_coo) <- extracted_data$name

# scl_coo |> list_flatten(name_spec = "{outer}") |> 
#   list_flatten(name_spec = "{outer}")  %>%  
#   map_if(.p = ~length(dim(.x))==3,.f = \(x) matrix(x,nrow=dim(x)[2],ncol=dim(x)[3])) %>% 
#   map_dfr(~as.data.frame(.x),.id = "Pathology_Annotation")%>% dplyr::rename(X = V1,Y = V2)  %>% 
#   mutate(unique_id = paste(Pathology_Annotation, row_number(), sep = "_")) %>% dim()
# library(httpgd)
# hgd()

extracted_data$scaled_coordinates %>% .[[1]] %>% matrix(nrow = dim(.)[2], ncol = dim(.)[3]) %>% 
  as.data.frame()  %>% dplyr::rename(X = V1,Y = V2) %>% 
  ggplot(aes(x = X, y = Y)) +geom_path() + theme_minimal() +labs(title = "Scaled Coordinates Path Plot", x = "X", y = "Y")
scl_coo_list <- scl_coo %>% list_flatten(name_spec = "{outer}")  %>% 
  list_flatten(name_spec = "{outer}")  %>%  
  map_if(.p = ~length(dim(.x))==3,.f = \(x) matrix(x,nrow=dim(x)[2],ncol=dim(x)[3]))
length(scl_coo_list)
scl_coo_df <- map2_dfr(scl_coo_list,seq(length(scl_coo_list)),\(x,y){
  as.data.frame(x) %>% mutate(region = str_c('region_',y))
},.id="Pathology_Annotation")  %>% 
  dplyr::rename(X = V1,Y = V2)

anno <-  scl_coo_df %>%
  ggplot(aes(x = X, y = Y,color = Pathology_Annotation, group = region))  +
  geom_path() + theme_minimal() +labs(title = "Scaled Coordinates Path Plot", x = "X", y = "Y")
anno

# spes_df
spes_df
spe <- subset(spes_df,samples==sample_nm,select=spes)[3,1][[1]][[1]];spe # samples== "SN048_A121573_Rep1"
spe_gt <- subset(spes_df,samples==sample_nm,select=spes)[5,1][[1]][[1]]


heat <- heatmap_func(spe,gene);heat
heat_gt <- heatmap_func(spe_gt,gene);heat_gt

combined_plot <- heat + 
  geom_path(data = mutate(anno$data, X = X, Y = -Y + 600), 
            aes(x = X, y = Y, group = region, color = factor(Pathology_Annotation)), 
            inherit.aes = FALSE,linewidth = 0.5) + 
  scale_color_manual(
    values = c("AggregatedImmuneCell" = "#00FFFF", 
               "Tumor" = "#FF0000", 
               "Stroma" = "#00FF00",
               "Necrosis" = '#D8DAD9',
               'Non-neoplastic_Mucosa' = '#FBDE28',
               'Normal_Epithelium' = '#FBDE28',
               'RBC' = '#ED740C','Muscle' = '#0000FF',
               'Not Defined' = '#262626'
    ),
    drop = TRUE  # 确保去掉不在指定范围内的颜色，如transparent
  ) +  
  guides(
    fill = guide_colorbar(order = 1), 
    color = guide_legend(order = 2, title = "Region")
  ) +
  labs(color = "Region")  # 更改图例标题为 "Region"
combined_plot 

spe$sample_id='Pathology Annotation'
HE <- heatmap_func(spe,gene)+
  geom_path(data = mutate(anno$data, X = X, Y = -Y + 600), 
            aes(x = X, y = Y, group = region, color = factor(Pathology_Annotation)), inherit.aes = FALSE,linewidth = 0.7) + 
  scale_color_manual(
    values = c("AggregatedImmuneCell" = "#00FFFF", 
               "Tumor" = "#FF0000", 
               "Stroma" = "#00FF00",
               "Necrosis" = '#D8DAD9',
               'Non-neoplastic_Mucosa' = '#FBDE28',
               'Normal_Epithelium' = '#FBDE28',
               'RBC' = '#ED740C','Muscle' = '#0000FF',
               'Not Defined' = '#262626'
    ),
    drop = TRUE  # 确保去掉不在指定范围内的颜色，如transparent
  ) +  
  guides(
    fill = guide_colorbar(order = 1), 
    color = guide_legend(order = 2, title = "Region")
  ) +
  labs(color = "Region")  # 更改图例标题为 "Region"
HE$layers <- HE$layers[sapply(HE$layers, function(layer) !inherits(layer$geom, "GeomPoint"))]

combined_plot_gt <- heat_gt + geom_path(data = mutate(anno$data, X = X, Y = -Y + 600), 
                                        aes(x = X, y = Y, group = region, color = factor(Pathology_Annotation)), inherit.aes = FALSE,linewidth = 0.7) + 
  scale_color_manual(
    values = c("AggregatedImmuneCell" = "#00FFFF", 
               "Tumor" = "#FF0000", 
               "Stroma" = "#00FF00",
               "Necrosis" = '#D8DAD9',
               'Non-neoplastic_Mucosa' = '#FBDE28',
               'Normal_Epithelium' = '#FBDE28',
               'RBC' = '#ED740C','Muscle' = '#0000FF',
               'Not Defined' = '#262626'
    ),
    drop = TRUE  # 确保去掉不在指定范围内的颜色，如transparent
  ) +  
  guides(
    fill = guide_colorbar(order = 1), 
    color = guide_legend(order = 2, title = "Region")
  ) +
  labs(color = "Region")  # 更改图例标题为 "Region"

HE/combined_plot_gt/combined_plot

######### PathologyCombinedPlots Cycle ###########
(sample_names <- list.dirs('PathologyAnno/External_Pathology_Annotation_JSON') %>% 
   basename() %>% .[-1] )
(sample_names <- sample_names[sample_names != "SN123_A551763_Rep1"])
genes <- c('KRT8','EPCAM','CD9','VIM','KLF5','HSP90AA1')
PathologyCombinedPlots <- 'Figures/PathologyCombinedPlots'
dir.create(PathologyCombinedPlots)
linewidth = 0.5
DPI = 120

(sample_nm = sample_names[9])
up_Y = 45  # 34CRC~80，Control_Rep1~20，Human~20，SN048_A121573_Rep1~170，SN048_A416371_Rep1~120
# SN123_A595688_Rep1~210(up_Y),X也需要一并调整：data = mutate(anno$data, X = X+10
# SN123_A798015_Rep1~115，SN124_A938797_Rep2~45，SN84_A120838_Rep2~125

for(gene in genes){
  extracted_data <- paste0("PathologyAnno/External_Pathology_Annotation_JSON/",sample_nm,'/',
                           sample_nm,'_extracted_data.json')  %>% fromJSON()
  scl_coo <- extracted_data$scaled_coordinates
  names(scl_coo) <- extracted_data$name
  # extracted_data$scaled_coordinates %>% .[[1]] %>% matrix(nrow = dim(.)[2], ncol = dim(.)[3]) %>%
  #   as.data.frame()  %>% dplyr::rename(X = V1,Y = V2) %>%
  #   ggplot(aes(x = X, y = Y)) +
  # geom_path() + theme_minimal() +labs(title = "Scaled Coordinates Path Plot", x = "X", y = "Y")
  scl_coo_list <- scl_coo %>% list_flatten(name_spec = "{outer}")  %>%
    list_flatten(name_spec = "{outer}")  %>%  
    map_if(.p = ~length(dim(.x))==3,.f = \(x) matrix(x,nrow=dim(x)[2],ncol=dim(x)[3]))
  scl_coo_df <- map2_dfr(scl_coo_list,seq(length(scl_coo_list)),\(x,y){
    as.data.frame(x) %>% mutate(region = str_c('region_',y))
  },.id="Pathology_Annotation")  %>% 
    dplyr::rename(X = V1,Y = V2)
  
  anno <-  scl_coo_df %>%
    ggplot(aes(x = X, y = Y,color = Pathology_Annotation, group = region))  +
    geom_path() + theme_minimal() +labs(title = "Scaled Coordinates Path Plot", x = "X", y = "Y")
  spe <- subset(spes_df,samples==sample_nm,select=spes)[3,1][[1]][[1]] # samples== "SN048_A121573_Rep1"
  spe_gt <- subset(spes_df,samples==sample_nm,select=spes)[5,1][[1]][[1]]
  heat <- heatmap_func(spe,gene)
  add_Y <- heat$data$pxl_row_in_fullres %>% max() %>% round() + up_Y
  
  heat_gt <- heatmap_func(spe_gt,gene)
  combined_plot <- heat + 
    geom_path(data = mutate(anno$data, X = X, Y = -Y + add_Y), 
              aes(x = X, y = Y, group = region, color = factor(Pathology_Annotation)), 
              inherit.aes = FALSE,linewidth = linewidth) + 
    scale_color_manual(
      values = c("AggregatedImmuneCell" = "#00FFFF", 
                 "Tumor" = "#FF0000", 
                 "Stroma" = "#00FF00",
                 "Necrosis" = '#D8DAD9',
                 'Non-neoplastic_Mucosa' = '#FBDE28',
                 'Normal_Epithelium' = '#FBDE28',
                 'RBC' = '#ED740C','Muscle' = '#0000FF',
                 'Not Defined' = '#262626'
      ),
      drop = TRUE  # 确保去掉不在指定范围内的颜色，如transparent
    ) +  
    guides(
      fill = guide_colorbar(order = 1), 
      color = guide_legend(order = 2, title = "Region")
    ) +
    labs(color = "Region")  # 更改图例标题为 "Region"
  spe$sample_id='Pathology Annotation'
  HE <- heatmap_func(spe,gene)+
    geom_path(data = mutate(anno$data, X = X, Y = -Y + add_Y), 
              aes(x = X, y = Y, group = region, color = factor(Pathology_Annotation)), 
              inherit.aes = FALSE,linewidth = linewidth) + 
    scale_color_manual(
      values = c("AggregatedImmuneCell" = "#00FFFF", 
                 "Tumor" = "#FF0000", 
                 "Stroma" = "#00FF00",
                 "Necrosis" = '#D8DAD9',
                 'Non-neoplastic_Mucosa' = '#FBDE28',
                 'Normal_Epithelium' = '#FBDE28',
                 'RBC' = '#ED740C','Muscle' = '#0000FF',
                 'Not Defined' = '#262626'
      ),
      drop = TRUE  # 确保去掉不在指定范围内的颜色，如transparent
    ) +  
    guides(
      fill = guide_colorbar(order = 1), 
      color = guide_legend(order = 2, title = "Region")
    ) +
    labs(color = "Region")  # 更改图例标题为 "Region"
  HE$layers <- HE$layers[sapply(HE$layers, function(layer) !inherits(layer$geom, "GeomPoint"))]
  
  combined_plot_gt <- heat_gt + geom_path(data = mutate(anno$data, X = X, Y = -Y + add_Y), 
                                          aes(x = X, y = Y, group = region, color = factor(Pathology_Annotation)), 
                                          inherit.aes = FALSE,linewidth = linewidth) + 
    scale_color_manual(
      values = c("AggregatedImmuneCell" = "#00FFFF", 
                 "Tumor" = "#FF0000", 
                 "Stroma" = "#00FF00",
                 "Necrosis" = '#D8DAD9',
                 'Non-neoplastic_Mucosa' = '#FBDE28',
                 'Normal_Epithelium' = '#FBDE28',
                 'RBC' = '#ED740C','Muscle' = '#0000FF',
                 'Not Defined' = '#262626'
      ),
      drop = TRUE  # 确保去掉不在指定范围内的颜色，如transparent
    ) +  
    guides(
      fill = guide_colorbar(order = 1), 
      color = guide_legend(order = 2, title = "Region")
    ) +
    labs(color = "Region")  # 更改图例标题为 "Region"
  
  tiff(filename = paste0(PathologyCombinedPlots,'/',sample_nm,'_',gene,'.tif'),
       width = 1920, 
       height = 1017,
       units = "px",         # 直接指定像素宽高，避免单位转换
       res = DPI) 
  print(HE/combined_plot_gt/combined_plot)
  dev.off()
}

########### Pathology Percent Unit Test ############

pathoAnno_path <- paste0('PathologyAnno/SpotClass_from_Annotation/',
                         sample_nm,'_percentage.csv')
pathoAnno <- read.csv(pathoAnno_path,
                      check.names = F)

pathoAnno <- pathoAnno[,-1]
rownames(pathoAnno) <- c('Tumor','Stroma','AIC','Muscle','Necrosis','RBC',
                         'NNM','NE','None')
# pathoAnno
# dim(pathoAnno);View(pathoAnno[,1:5])
table(colnames(pathoAnno) == rownames(colData(spe)))


# 对每一列应用分类函数
classification_labels <- apply(pathoAnno, 2, max_classify_column)
# table(classification_labels)
# table(names(classification_labels) == rownames(colData(spe)))


##### BoxPlot between different PathoAnnotations and ROC 

my_colors <- colorRampPalette(brewer.pal(12,'Set3'))(20)  # 获取 12 种颜色

colData(spe)$patho <- classification_labels

df <- data.frame(
  gene_value = assay(spe)[gene, ], 
  patho = colData(spe)$patho
)

# ggboxplot 绘图
stformer_box <- df %>%
  mutate(patho = reorder(patho, gene_value, FUN = median)) %>%
  ggboxplot(
    x = "patho",
    y = "gene_value",  # 动态列名
    color = "patho",
    palette = my_colors,
    add = "jitter",
    outlier.shape = NA
  ) +
  labs(y = gene) + # 设置动态 y 轴标签
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  stat_compare_means(comparisons = list(c("Tumor", "NNM")))

# ggplot2 绘图
# df %>%
#   ggplot(aes(x = reorder(patho, gene_value), y = gene_value)) +
#   geom_boxplot() +
#   labs(y = gene) + # 设置动态 y 轴标签
#   theme(axis.text.x = element_text(angle = 45, hjust = 1))


df_gt <- data.frame(
  gene_value = assay(spe_gt)[gene, ], 
  patho = colData(spe)$patho
)
# ggboxplot 绘图
gt_box <- df_gt %>%
  mutate(patho = reorder(patho, gene_value, FUN = median)) %>%
  ggboxplot(
    x = "patho",
    y = "gene_value",  # 动态列名
    color = "patho",
    palette = my_colors,
    add = "jitter",
    outlier.shape = NA
  ) +
  labs(y = gene) + # 设置动态 y 轴标签
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  stat_compare_means(comparisons = list(c("Tumor", "NNM")))


# 定义布局：1#
#           23
#           45 的排布代表图像会排成3行2列的形式；其中第1幅图占据1行1列,1行2列为空白；
#                                                    第2幅图占据2行1列，第3幅图占据2行2列

layout <- "
12
34
"

combined_plot_gt+gt_box+combined_plot+stformer_box+plot_layout(design = layout)

################## Pathology Percent Cycle ######################
(sample_names <- list.dirs('PathologyAnno/External_Pathology_Annotation_JSON/') %>% basename() %>% .[-1] )
(sample_names <- sample_names[sample_names != "SN123_A551763_Rep1"])
genes <- c('KRT8','EPCAM','CD9','VIM','KLF5','HSP90AA1')
PathologyPercentPlots <- 'Figures/PathologyPercentPlots'
dir.create(PathologyPercentPlots)
compare_list = list(c("NNM", "Stroma"))# list(c("Tumor", "Stroma"))
linewidth = 0.5
DPI = 120
my_colors <- colorRampPalette(brewer.pal(12,'Set3'))(20)  # 获取 12 种颜色
(sample_nm = sample_names[4])
spe <- subset(spes_df,samples==sample_nm,select=spes)[3,1][[1]][[1]] # samples== "SN048_A121573_Rep1"
spe_gt <- subset(spes_df,samples==sample_nm,select=spes)[5,1][[1]][[1]]

for(gene in genes){
  pathoAnno_path <- paste0('PathologyAnno/SpotClass_from_Annotation/',
                           sample_nm,'_percentage.csv')
  pathoAnno <- read.csv(pathoAnno_path,
                        check.names = F)
  
  pathoAnno <- pathoAnno[,-1]
  rownames(pathoAnno) <- c('Tumor','Stroma','AIC','Muscle','Necrosis','RBC',
                           'NNM','NE','None')
  classification_labels <- apply(pathoAnno, 2, max_classify_column)
  colData(spe)$patho <- classification_labels
  
  df <- data.frame(
    gene_value = assay(spe)[gene, ], 
    patho = colData(spe)$patho
  )
  
  # ggboxplot 绘图
  stformer_box <- df %>%
    mutate(patho = reorder(patho, gene_value, FUN = median)) %>%
    ggboxplot(
      x = "patho",
      y = "gene_value",  # 动态列名
      color = "patho",
      palette = my_colors,
      add = "jitter",
      outlier.shape = NA
    ) +
    labs(y = gene) + # 设置动态 y 轴标签
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    stat_compare_means(comparisons = compare_list)
  
  # ggplot2 绘图
  # df %>%
  #   ggplot(aes(x = reorder(patho, gene_value), y = gene_value)) +
  #   geom_boxplot() +
  #   labs(y = gene) + # 设置动态 y 轴标签
  #   theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  
  df_gt <- data.frame(
    gene_value = assay(spe_gt)[gene, ], 
    patho = colData(spe)$patho
  )
  # ggboxplot 绘图
  gt_box <- df_gt %>%
    mutate(patho = reorder(patho, gene_value, FUN = median)) %>%
    ggboxplot(
      x = "patho",
      y = "gene_value",  # 动态列名
      color = "patho",
      palette = my_colors,
      add = "jitter",
      outlier.shape = NA
    ) +
    labs(y = gene) + # 设置动态 y 轴标签
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    stat_compare_means(comparisons = compare_list)
  layout <- "
  12
  34
  "
  tiff(filename = paste0(PathologyPercentPlots,'/',sample_nm,'_',gene,'.tif'),
       width = 1920, 
       height = 1017,
       units = "px",         # 直接指定像素宽高，避免单位转换
       res = DPI) 
  print(combined_plot_gt+gt_box+combined_plot+stformer_box+plot_layout(design = layout))
  dev.off()
}
dev.off()
