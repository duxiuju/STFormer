library(tidyverse)
library(flextable)
load('RData/spes_df.RData')
samples_order <- c('SN048_A121573_Rep1','SN048_A121573_Rep2','SN123_A798015_Rep1','SN124_A798015_Rep2',
                   'SN84_A120838_Rep1','SN84_A120838_Rep2','SN123_A938797_Rep1','SN124_A938797_Rep2',
                   'SN048_A416371_Rep1','SN048_A416371_Rep2','SN123_A595688_Rep1','SN124_A595688_Rep2',
                   # 'SN123_A551763_Rep1','SN124_A551763_Rep2',
                   '34CRC','Control_Rep1','Control_Rep2','Human')
Patient_order <- str_replace_all(samples_order, "^SN\\d+_(A\\d+)_Rep\\d+$", "\\1") %>% unique()
models_order <- c('Hist2ST','HisToGene','STNet','STFormer','GroundTruth')

# spes_df 
subset(spes_df,models=='STFormer')  %>% 
  mutate(intergenes = map(corr,~.x[.x>0.2] %>% na.omit())) %>% .$samples

ex <- subset(spes_df,models=='STFormer')  %>% 
  mutate(intergenes = map(corr,~.x[.x>0.2] %>% na.omit())) %>% 
  .$intergenes %>% .[seq(4)] %>% map(names) %>% Reduce(intersect,.) %>% sort()
# save(ex,file = 'ex39Genes.rds')
intern <- subset(spes_df,models=='STFormer')  %>% 
  mutate(intergenes = map(corr,~.x[.x>0.2] %>% na.omit())) %>%  
  .$intergenes %>% .[-seq(4)] %>% map(names) %>% Reduce(intersect,.) %>% sort()
table(ex %in% intern)
intersect(ex,intern);setdiff(ex,intern)
############################# 一、 14CRC summary ################################
i = seq(nrow(subset(spes_df,samples%in%samples_order[seq(12)])))%%(length(models_order))==0;i
ft <- subset(spes_df,samples%in%samples_order[seq(12)],select=c(samples,models,corr))%>% 
  group_by(samples,models) %>%  mutate(models = factor(models, levels = models_order),
                                        samples = factor(samples, levels = samples_order)) %>% 
  arrange(samples, models) %>% 
  mutate('Target Gene Number' = 418-sum(is.na(corr[[1]])),
         'Median correlation' = median(corr[[1]],na.rm = T),
         'Mean correlation' = mean(corr[[1]],na.rm = T),
         'Ratio of correlation ≥ 0.20' = mean(corr[[1]]>=0.2,na.rm = T),
         'Ratio of correlation ≥ 0.30' = mean(corr[[1]]>=0.3,na.rm = T),
         'Ratio of correlation ≥ 0.40' = mean(corr[[1]]>=0.4,na.rm = T),
         'Ratio of correlation ≥ 0.50' = mean(corr[[1]]>=0.5,na.rm = T)) %>% 
  mutate(across(where(is.numeric),~round(.x,3))) %>% 
  dplyr::rename(`Samples ID`='samples') %>% 
  ungroup() %>% select(-corr) %>% distinct() %>% 
  flextable()  %>% merge_v(j='Samples ID') %>% 
  flextable::bold( i = ~ models == "STFormer", bold = TRUE) %>% 
  hline(i = i) %>% fontsize(size = 8, part = "all") %>% 
  border_outer(border = fp_border_default(width = 1, color = "black"))

if (!dir.exists('ThreeLineTable/')) {
  # 如果文件夹不存在，则创建文件夹
  dir.create('ThreeLineTable/')
  message("Folder created successfully.")
}
save_as_docx(ft, path = "ThreeLineTable/14CRCSummary.docx")
save_as_image(ft, path = "ThreeLineTable/14CRCSummary.png")
# save_as_pptx(ft, path = "ThreeLineTable/14CRCSummary.pptx")
############################# 二、 external summary ################################
i = seq(nrow(subset(spes_df,samples%in%samples_order[-seq(12)])))%%(length(models_order))==0;i
ft <- subset(spes_df,samples%in%samples_order[-seq(12)],select=c(samples,models,corr))%>% 
  group_by(samples,models) %>%  mutate(models = factor(models, levels = models_order),
                                       samples = factor(samples, levels = samples_order)) %>% 
  arrange(samples, models) %>% 
  mutate('Target Gene Number' = 418-sum(is.na(corr[[1]])),
         'Median correlation' = median(corr[[1]],na.rm = T),
         'Mean correlation' = mean(corr[[1]],na.rm = T),
         'Ratio of correlation ≥ 0.20' = mean(corr[[1]]>=0.2,na.rm = T),
         'Ratio of correlation ≥ 0.30' = mean(corr[[1]]>=0.3,na.rm = T),
         'Ratio of correlation ≥ 0.40' = mean(corr[[1]]>=0.4,na.rm = T),
         'Ratio of correlation ≥ 0.50' = mean(corr[[1]]>=0.5,na.rm = T)) %>% 
  mutate(across(where(is.numeric),~round(.x,3))) %>% 
  dplyr::rename(`Samples ID`='samples') %>% 
  ungroup() %>% select(-corr) %>% distinct() %>% 
  flextable()  %>% merge_v(j='Samples ID') %>% 
  flextable::bold( i = ~ models == "STFormer", bold = TRUE) %>% 
  hline(i = i) %>% fontsize(size = 8, part = "all") %>% 
  border_outer(border = fp_border_default(width = 1, color = "black"))

save_as_docx(ft, path = "ThreeLineTable/external_Summary.docx")
save_as_image(ft, path = "ThreeLineTable/external_Summary.png")
############################ 三、clinical characteristics ######################### 
num_model <- 5
crc14 <- 12*num_model
ft <- ungroup(spes_df) %>% select(samples) %>% 
  mutate(samples = factor(samples, levels = samples_order)) %>% 
  dplyr::arrange(samples) %>% mutate(
             `Data Sets` = factor(c(rep('Leave-one-patient-out validation',crc14),
                                    rep('External 1',num_model),
                                    rep('External 2',num_model),rep('External 2',num_model),rep('External 3',num_model))),
             `Localization` = factor(c(rep('Rectum',num_model*2),rep('Sigma/Rectum',num_model*2),
                                       rep('Colon(Sigma)',num_model*2),rep('Rectum',num_model*2),
                                       rep('Colon(right)',num_model*2),rep('Colon(right)',num_model*2),
                                       
                                       # rep('Cecum',num_model*2),
                                       rep('Large Intestine',num_model),rep('Colon',num_model),rep('Colon',num_model),rep('Large Intestine',num_model)
                                       )),
             `Patient ID` = factor(c(rep('A121573',num_model*2),rep('A798015',num_model*2),
                                     rep('A120838',num_model*2),rep('A938797',num_model*2),
                                     rep('A416371',num_model*2),rep('A595688',num_model*2),
                                     
                                     # rep('A551763',num_model*2),
                                     rep('34CRC',num_model),rep('Control_Rep1',num_model),rep('Control_Rep2',num_model),rep('Human',num_model)
                                     )),
             `Spots Under Tissue`= c(rep(2203,num_model),rep(2385,num_model),rep(1685,num_model),rep(1656,num_model),
                                     rep(328,num_model),rep(1048,num_model),rep(2128,num_model),rep(1691,num_model),
                                     rep(2317,num_model),rep(1803,num_model),rep(1192,num_model),rep(387,num_model),
                                     # rep(691,num_model),rep(1219,num_model),
                                     rep(2660,num_model),rep(6487,num_model),rep(6414,num_model),rep(9080,num_model)
                                     ),
             'Median Genes per Spot' = c(rep(4264,num_model),rep(3809,num_model),rep(2343,num_model),rep(2692,num_model),
                                         rep(3958,num_model),rep(3348,num_model),rep(3084,num_model),rep(5457,num_model),
                                         rep(4116,num_model),rep(4588,num_model),rep(4388,num_model),rep(4407,num_model),
                                         # rep(4643,num_model),rep(1233,num_model),
                                         rep(7438,num_model),rep(3018,num_model),rep(2404,num_model),rep(9560,num_model)
                                         )) %>% 
  dplyr::rename(`Samples ID`='samples')  %>% 
  select(c('Data Sets','Patient ID','Localization','Samples ID','Spots Under Tissue','Median Genes per Spot'))  %>% 
  distinct() %>% flextable() %>% merge_v(j=c('Data Sets','Patient ID','Localization','Spots Under Tissue')) %>% 
  # align(j = c('Localization','Spots Under Tissue'), align = 'center', part = 'body') %>% 
  fontsize(size = 8, part = "all") %>% 
  hline(i = c(2,4,6,8,10,12,13,15),
        border = fp_border_default(width = 1, color = "black"))%>% 
  border_outer(border = fp_border_default(width = 1, color = "black"))

# for (sample in samples_order2) {
#   last_row <- which(spes_df$samples== sample) %>% tail(1)
#   ft <- ft %>% hline(i = last_row, border = fp_border_default(width = 1, color = "black"))
# }
ft
save_as_image(ft, path = "ThreeLineTable/clinicalCharacteristics.png") 
save_as_docx(ft, path = "ThreeLineTable/clinicalCharacteristics.docx")                                    
# save_as_pptx(ft, path = "ThreeLineTable/clinicalCharacteristics.pptx")


############################ 四、Genes Per Sample ######################### 
new_df <- ungroup(spes_df) %>% select(c(samples,models,corr))%>% 
  mutate(samples = factor(samples, levels = samples_order),
         models = factor(models, levels = models_order)) %>% 
  dplyr::arrange(samples,models) %>% mutate(
    `Patient ID` = factor(c(rep('A121573',num_model*2),rep('A798015',num_model*2),
                            rep('A120838',num_model*2),rep('A938797',num_model*2),
                            rep('A416371',num_model*2),rep('A595688',num_model*2),
                            # rep('A551763',num_model*2),
                            rep('34CRC',num_model),rep('Control_Rep1',num_model),
                            rep('Control_Rep2',num_model),rep('Human',num_model)),
                          levels = Patient_order)) %>% 
  group_by(samples,models) %>% 
  mutate(top10gene = list(round(sort(corr[[1]],decreasing =T)[1:10],3)),
         top10geneSymbol = list(names(sort(corr[[1]],decreasing =T)[1:10])),
         combined = list(map2(top10geneSymbol,top10gene,~paste(.x,.y,sep='\n')))) %>% 
  # .[['combined']]
  select(`Patient ID`,models,samples,combined) %>% unnest(combined) %>% 
  unnest_wider(combined,names_sep = '') %>% 
  dplyr::rename(`Samples ID` = 'samples') %>% 
  dplyr::rename_with(~str_replace(., "combined", "top "), starts_with("combined"))

ii = seq(nrow(subset(new_df,`Samples ID`%in%samples_order[seq(12)])))%%
  (2*length(models_order))==0
subset(new_df,`Samples ID`%in%samples_order[seq(12)])%>% 
  group_by(`Patient ID`) %>% arrange(`Patient ID`,models) %>% 
  flextable() %>% merge_v(j=c('Patient ID','models')) %>% 
  # align(j = c('Localization','Spots Under Tissue'), align = 'center', part = 'body') %>% 
  fontsize(size = 8, part = "all") %>% 
  hline(i = ii) %>% vline(j='Patient ID')%>% 
  flextable::bold( i = ~ models == "STFormer", bold = TRUE) %>% 
  border_outer(border = fp_border_default(width = 1, color = "black")) %>% 
  #save_as_image(path = "ThreeLineTable/Top10Genes_14CRC.png")
  save_as_docx(path = "ThreeLineTable/Top10Genes_14CRC.docx")

ii_ = seq(nrow(subset(new_df,`Samples ID`%in%samples_order[-seq(12)])))%%
  (length(models_order))==0

subset(new_df,`Samples ID`%in%samples_order[-seq(12)])%>% 
  group_by(`Patient ID`) %>% arrange(`Patient ID`,models) %>% 
  flextable() %>% merge_v(j=c('Patient ID','models')) %>% 
  # align(j = c('Localization','Spots Under Tissue'), align = 'center', part = 'body') %>% 
  fontsize(size = 8, part = "all") %>% 
  hline(i = ii_) %>% vline(j='Patient ID') %>% flextable::bold( i = ~ models == "STFormer", bold = TRUE) %>% 
  border_outer(border = fp_border_default(width = 1, color = "black")) %>% 
  save_as_image(path = "ThreeLineTable/Top10Genes_External.png")
  #save_as_docx(path = "ThreeLineTable/Top10Genes_External.docx")



