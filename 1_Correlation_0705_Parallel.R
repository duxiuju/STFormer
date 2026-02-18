rm(list = ls()) 
source('1-Correlation_0705_functions_Parallel.R')
################################## Training and Validation #########################################
###################################                        #########################################
outputs <- list()
(models <- list.files('./PredictedST/14CRC/')) # [2:4]
for(model in models){
  print(stringr::str_glue('model:{model} of 14CRC is processing'))
  model
  path <- paste0("./PredictedST/14CRC/",model,"/test/analysis")
  outputs[[model]] <- Cor_14CRC(path)
}

unique(outputs$HisToGene$correlation_results$Sample_Name)
# outputs$HisToGene$correlation_results$Sample_Name <- outputs$STNet$correlation_results$Sample_Name
# outputs$Hist2ST$correlation_results$Sample_Name <- outputs$STNet$correlation_results$Sample_Name

############################## 14CRC Combined Sample Correlation ######################
all_correlation_results <- map_df(names(outputs), ~ {
  outputs[[.x]]$correlation_results %>%
    mutate(model = .x)  # 在数据框中添加新列
})
all_correlation_results$model <- factor(all_correlation_results$model, 
                                        levels = c("Hist2ST","HisToGene", "STNet", "STFormer"))
all_correlation_results %>% filter(!is.na(Correlation)) %>% 
  filter(!is.infinite(Correlation)) %>% 
  mutate(Sample_Name = factor(Sample_Name,
                              levels = c('A121573','A798015','A120838',
                                         'A938797','A416371','A595688'))) %>% 
  ggboxplot(x = "model", y = "Correlation", 
            color = "model", palette = "OrRd",
            add = "jitter", 
            outlier.shape = NA) +
  theme_minimal()+ xlab("Sample") + 
  stat_compare_means(comparisons = list(c("STFormer", "STNet"), 
                                        c("STFormer", "HisToGene"), 
                                        c("STFormer", "Hist2ST"),
                                        c("STNet", "HisToGene"),
                                        c("STNet", "Hist2ST"),
                                        c("HisToGene", "Hist2ST")),
                     method = "wilcox.test",size = 1.5)+
  stat_summary(fun = mean, geom = "point", shape = 18, size = 1, color = "blue")+
  scale_x_discrete(labels = NULL) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1,size = 2),
        strip.placement = "outside", 
        strip.background = element_blank(),
        strip.text.x = element_text(face = "bold")) +
  facet_wrap(~Sample_Name, nrow = 1, strip.position = "bottom")


# 小提琴图：ggbetweenstats(all_correlation_results,x= model,y=Correlation,color = model,type = 'nonparametric')

#######################      External Test ###################################
##########################                 ##################################

# 罗列当前文件夹下所有的文件夹
(datasets<- c("34CRC","Control_Rep1","Control_Rep2","Human"))

models_list <- c("Hist2ST", "HisToGene", "STFormer", "STNet")

# 使用expand.grid生成所有可能的数据集和模型组合
(combinations <- expand.grid(Datasets = datasets, Models = models_list,
                             stringsAsFactors = FALSE))
library(furrr)
plan(multisession, workers = 12)
outputs_external <- furrr::future_map2(
  combinations$Dataset,
  combinations$Model,
  ~{
    gt_path <- file.path('PredictedST',.x, .y, "gt_with_gene.csv")
    pre_path <- file.path('PredictedST',.x, .y, "pre_with_gene.csv")
    gt <- read_csv(gt_path) %>% select(-1) %>% column_to_rownames(var = "gene")
    pre <- read_csv(pre_path) %>% select(-1) %>% column_to_rownames(var = "gene") %>% .[rownames(gt),]
    return(list(gt = gt,pre=pre))},.progress = T
) %>% set_names(paste(combinations$Dataset, combinations$Model, sep = "_"))


results <- imap_dfr(outputs_external, ~Cor_exteral(.x$gt,.x$pre,.y))
results$model <- factor(results$model, levels = c("Hist2ST","HisToGene", "STNet", "STFormer"))
plan(sequential)

# library(httpgd)
# hgd()
ggboxplot(results,x = "model", y = "Correlation", 
          color = "model", palette = "OrRd",
          add = "jitter", 
          outlier.shape = NA) +
  theme_minimal()+ xlab("Sample") + 
  stat_compare_means(comparisons = list(c("STFormer", "STNet"), 
                                        c("STFormer", "HisToGene"), 
                                        c("STFormer", "Hist2ST"),
                                        c("STNet", "HisToGene"),
                                        c("STNet", "Hist2ST"),
                                        c("HisToGene", "Hist2ST")
  ),,size = 1.5)+
  scale_x_discrete(labels = NULL) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        strip.placement = "outside", 
        strip.background = element_blank(),
        strip.text.x = element_text(face = "bold")) +
  facet_wrap(~Sample_Name, nrow = 1, strip.position = "bottom")
# # 小提琴图：ggbetweenstats(results,x= model,y=Correlation,color = model,type = 'nonparametric')

################### deal with the spotNames for HistoGene and Hist2ST, then save the data ################################
# 已source('Hist2ST_toGene_ColNames.R') # 
# unique(outputs_external$`34CRC_Hist2ST`$correlation_results$Sample_Name)
save(outputs,outputs_external,file = 'RData/gt_pre_TrainValTest_ColName0515.RData')

############################## 14CRC separate Sample Correlation ######################
separate_14CRC <- map_df(names(outputs), function(model_name) {
  calculate_correlation(outputs[[model_name]], model_name)
})
separate_14CRC$model <- factor(separate_14CRC$model, 
                               levels = c("Hist2ST","HisToGene", "STNet", "STFormer"))
separate_14CRC$Sample_Name <- factor(separate_14CRC$Sample_Name, 
                                     levels = c('SN048_A121573_Rep1','SN123_A798015_Rep1','SN84_A120838_Rep1',
                                                'SN123_A938797_Rep1','SN048_A416371_Rep1','SN123_A595688_Rep1',
                                                'SN048_A121573_Rep2','SN124_A798015_Rep2','SN84_A120838_Rep2',
                                                'SN124_A938797_Rep2','SN048_A416371_Rep2','SN124_A595688_Rep2'
                                     ))
# separate_14CRC
ggboxplot(separate_14CRC,x = "model", y = "Correlation", 
          color = "model", palette = "OrRd",
          add = "jitter", 
          outlier.shape = NA) +
  theme_minimal()+ xlab("Sample_Name") + 
  stat_compare_means(comparisons = list(c("STFormer", "STNet"), 
                                        c("STFormer", "HisToGene"), 
                                        c("STFormer", "Hist2ST"),
                                        c("STNet", "HisToGene"),
                                        c("STNet", "Hist2ST"),
                                        c("HisToGene", "Hist2ST")
  ),size=1.5)+
  scale_x_discrete(labels = NULL) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        strip.placement = "outside", 
        strip.background = element_blank(),
        strip.text.x = element_text(face = "bold")) +
  facet_wrap(~Sample_Name, nrow = 2, strip.position = "bottom")

save(all_correlation_results,separate_14CRC,results,file = 'RData/Correlations.RData')

############################ 5. gene correlations perform well across all the samples ############
############################ 5.1 combined 14CRC 2sections with external 
# load('RData/Correlations.RData')
cutoff = 0.2
model_ = 'STFormer'
CRC14_genes <- all_correlation_results %>% filter(Sample_Name!='A551763',model ==model_,Correlation > cutoff) %>% 
  group_by(Sample_Name) %>% 
  summarise(Gene_List = list(Gene_Name)) %>% .$Gene_List %>% reduce(intersect)
# all_correlation_results %>% subset(model ==model_,Sample_Name) %>% group_by(Sample_Name) %>% slice_max(order_by=Correlation,n=20) %>% .$Gene_Name %>% unique()
# external_genes <- results %>% subset(model ==model_) %>% group_by(Sample_Name)  %>% slice_max(order_by=Correlation,n=20) %>% .$Gene_Name %>% unique()
CRC14_genes %>% sort()
external_genes <- results %>% filter(model ==model_,Correlation > cutoff) %>% group_by(Sample_Name) %>% 
  summarise(Gene_List = list(Gene_Name)) %>% .$Gene_List %>% reduce(intersect)
external_genes  %>% sort()

(inter_02 <- intersect(external_genes,CRC14_genes) %>% sort()) #取交集
save(inter_02,file = 'RData/inter_02Genes.RData')
# (inter_02 <- union(external_genes,CRC14_genes) %>% sort()) #取并集
table(external_genes %in% CRC14_genes)

all_correlation_results %>% filter(model ==model_,Gene_Name %in% CRC14_genes) %>% 
  pivot_wider(id_cols = Sample_Name,names_from = Gene_Name,values_from = Correlation) %>% 
  column_to_rownames(var = 'Sample_Name') %>% 
  pheatmap::pheatmap(cluster_rows = F,  # 是否对行进行聚类
                     cluster_cols = TRUE,  # 是否对列进行聚类
                     # scale = "row",        # 按行进行归一化
                     # color = colorRampPalette(c("blue", "white", "red"))(50),  # 自定义颜色
                     show_rownames = TRUE,  # 显示行名
                     show_colnames = TRUE,  # 显示列名
                     angle_col = 45,        # x轴字体倾斜45度
                     fontsize_col = 6,
                     main = "Sample Heatmap") # 热图标题


library(RColorBrewer)
library(paletteer)
# View(palettes_c_names)
# display.brewer.all(type= "seq")
# 使用 inferno 调色板创建 50 个颜色渐变点
my_colors <-  colorRampPalette(brewer.pal(6, "YlOrRd"))(50)# paletteer_c( "grDevices::Teal", n = 40) # scico(n = 50, palette = "lajolla")# 
rbind(all_correlation_results,results)  %>% 
  filter(Sample_Name!='A551763',model ==model_,Gene_Name %in% inter_02) %>% 
  mutate(Sample_Name = factor(Sample_Name,levels = c("A121573","A798015","A120838","A938797","A416371","A595688",
                                                     "Control_Rep1","Control_Rep2","34CRC","Human"))) %>% 
  pivot_wider(id_cols = Sample_Name,names_from = Gene_Name,values_from = Correlation) %>% 
  column_to_rownames(var = 'Sample_Name') %>% 
  select(sort(colnames(.)))%>% 
  pheatmap::pheatmap(cluster_rows = T,  # 是否对行进行聚类
                     cluster_cols = F,  # 是否对列进行聚类
                     # scale = "row",        # 按行进行归一化
                     color = my_colors,  # 自定义颜色
                     show_rownames = TRUE,  # 显示行名
                     show_colnames = TRUE,  # 显示列名
                     angle_col = 45,        # x轴字体倾斜45度
                     fontsize_col = 6,
                     main = "Sample Heatmap") # 热图标题

############################ 5.2 separate 14CRC 2sections with external ################
cutoff = 0.2
model_ = 'STFormer'

(CRC14Separate_genes <- separate_14CRC %>% filter(model ==model_,Correlation > cutoff) %>% 
  group_by(Sample_Name) %>% 
  summarise(Gene_List = list(Gene_Name)) %>% .$Gene_List %>% reduce(intersect)  %>% sort() )

(external_genes <- results %>% filter(model ==model_,Correlation > cutoff) %>% group_by(Sample_Name) %>% 
  summarise(Gene_List = list(Gene_Name)) %>% .$Gene_List %>% reduce(intersect) %>% sort() )

(inter_sep_02 <- intersect(CRC14Separate_genes,external_genes) %>% sort())
save(inter_sep_02,file = 'RData/inter_sep_02.RData')
rbind(separate_14CRC,results)  %>% 
  filter(model ==model_,Gene_Name %in% inter_sep_02) %>% 
  mutate(Sample_Name = factor(Sample_Name,levels = c('SN048_A121573_Rep1','SN048_A121573_Rep2','SN123_A798015_Rep1','SN124_A798015_Rep2',
                                                     'SN84_A120838_Rep1','SN84_A120838_Rep2','SN123_A938797_Rep1','SN124_A938797_Rep2',
                                                     'SN048_A416371_Rep1','SN048_A416371_Rep2','SN123_A595688_Rep1','SN124_A595688_Rep2',
                                                     "Control_Rep1","Control_Rep2","34CRC","Human"))) %>% 
  pivot_wider(id_cols = Sample_Name,names_from = Gene_Name,values_from = Correlation) %>% 
  column_to_rownames(var = 'Sample_Name') %>% 
  select(sort(colnames(.)))%>% 
  pheatmap::pheatmap(cluster_rows = T,  # 是否对行进行聚类
                     cluster_cols = F,  # 是否对列进行聚类
                     # scale = "row",        # 按行进行归一化
                     color = my_colors,  # 自定义颜色
                     show_rownames = TRUE,  # 显示行名
                     show_colnames = TRUE,  # 显示列名
                     angle_col = 45,        # x轴字体倾斜45度
                     fontsize_col = 6,
                     main = "Sample Heatmap") # 热图标题

