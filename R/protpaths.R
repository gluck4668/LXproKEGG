#--------------kegg all pathways------------------------#

protpaths <- function(protein_data){

#-----------读取数据，并筛选出差异表达蛋白-------------------------
  pro_df <- read.xlsx(trimws(protein_data))

  if(ncol(pro_df)>2){
    colnames(pro_df) <- c("Protein_ID","FC","pvalue")
    pro_df <- dplyr::filter(pro_df,pvalue<0.05) #筛选出差异表达蛋白
    is_inf <- grep("inf", pro_df$FC,ignore.case = T) #判断是否有无穷大值Inf
    if(length(is_inf)>0)
      pro_df <- pro_df[-c(is_inf),] #去掉无穷大的行
  }
  pro_df <- na.omit(pro_df) #去掉NA
  colnames(pro_df)[1] <- c("Protein_ID")
  pro_df <- distinct(pro_df,Protein_ID,.keep_all = T)

#----------- 蛋白enrichKEGG通路富集分析---------------------------
kegg_pro <- enrichKEGG(pro_df[,1], organism ='hsa',
                       #universe,
                       keyType = 'uniprot',  # 蛋白选择"uniprot" （如果是基因，则选择"kegg"）
                       pvalueCutoff = 0.05,
                       pAdjustMethod = 'BH',
                       qvalueCutoff = 0.2,
                       minGSSize = 3,
                       maxGSSize = 3500,
                       use_internal_data = F)
if(!is.null(kegg_pro))
  species <- "human" #判断物种

if(is.null(kegg_pro)){
  {kegg_pro <- enrichKEGG(pro_df[,1], organism ='rno',
                          #universe,
                          keyType = 'uniprot',  # 蛋白选择"uniprot" （如果是基因，则选择"kegg"）
                          pvalueCutoff = 0.05,
                          pAdjustMethod = 'BH',
                          qvalueCutoff = 0.2,
                          minGSSize = 3,
                          maxGSSize = 3500,
                          use_internal_data = F)
  species <- "rat"} #判断物种

  if(is.null(kegg_pro)){
    {kegg_pro <- enrichKEGG(pro_df[,1], organism ='mmu',
                            #universe,
                            keyType = 'uniprot',  # 蛋白选择"uniprot" （如果是基因，则选择"kegg"）
                            pvalueCutoff = 0.05,
                            pAdjustMethod = 'BH',
                            qvalueCutoff = 0.2,
                            minGSSize = 3,
                            maxGSSize = 3500,
                            use_internal_data = F)
    species <- "mouse" } #判断物种
  }

}


#---------------以物种为名字建立相应的文件夹--------------------
dir.file <- dplyr::case_when ( species== "human" ~ "analysis results (human)",
                               species== "mouse" ~ "analysis results (mouse)",
                               species== "rat" ~ "analysis results (rat)",
                               TRUE ~ "analysis results")

if (dir.exists(dir.file)==FALSE)
  dir.create(dir.file)

#--------------蛋白KEGG通路可视化-----------------------------------
kegg_all_pathways <- kegg_pro@result %>% na.omit()
names(kegg_all_pathways) <- gsub("Gene","Protein",names(kegg_all_pathways),ignore.case = T)
#-----把GeneRatio分数字符串变小数-------
#Ratio <- str_split(kegg_all_pathways$GeneRatio,"/",simplify = T) %>% apply(2,as.numeric)
#Ratio <- Ratio[,1]/Ratio[,2]

str_num <- function(i){eval(str2expression(i))} # 把字符串当做命令运行 eval(str2expression(i))
Ratio <- sapply(kegg_all_pathways$ProteinRatio,str_num,simplify = T)
kegg_all_pathways$ProteinRatio <- round(Ratio,4) #在pathways_kegg中添加一列ProteinRatio

#-----拆分Description-------------------
kegg_all_pathways$Description <- str_extract(kegg_all_pathways$Description,".*(?= -)") %>% trimws()

file_path_name <- paste(species,"Protein KEGG pathways (all)_data.xlsx")
file_path_name <-paste0(dir.file,"/",file_path_name)
write.xlsx(kegg_all_pathways,file_path_name)

row_n <- nrow(kegg_all_pathways)

if(row_n>=30)
  title_text <- c("Top 30 Protein KEGG Pathways") else
    title_text <- c("Protein KEGG Pathways")

title_size <- case_when(row_n>30 ~12,
                        row_n>20 ~12,
                        TRUE ~11)

xy_size <- case_when(row_n>30 ~9,
                     row_n>20 ~10,
                     TRUE ~10)

kegg_mytheme<-theme_bw()+
  theme(text=element_text(family = "sans",colour ="black",face="bold",size =title_size),
        panel.grid = element_blank(),
        panel.border = element_rect(linewidth = 0.8,colour = "gray30"),
        axis.line = element_blank(),
        axis.ticks = element_line(linewidth = 0.6,colour = "gray30"),
        axis.ticks.length = unit(1.5,units = "mm"))+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"))+
  theme(panel.grid =element_line(colour="#dcdcdc",linewidth=0.2,linetype = "dashed"))

kegg_xytheme <-theme(axis.text.x = element_text(face="bold",color="black",size=xy_size,angle =0,hjust=1))+
  theme(axis.text.y = element_text(face="bold",color="black",size=xy_size))+
  theme(legend.text=element_text(face="bold",color="black",size=xy_size))

#--------------kegg all pathways------------------------#
if(row_n<30)
  path_n <- row_n else
    path_n <- 30

kegg_df <- kegg_all_pathways[1:path_n,]

kegg_pathways <- ggplot(kegg_df)+
  geom_point(aes(x=ProteinRatio,
                 y=fct_reorder(Description,ProteinRatio),
                 color=-log10(pvalue),size=Count))+
  scale_color_gradient2(midpoint = 1,low = "blue",mid = "#ffcc00",high ="#dc143c" )+
  labs(x = 'ProteinRatio', y = '',title=title_text)+
  kegg_mytheme+kegg_xytheme

kegg_pathways

path_name <- paste(species,"Protein KEGG pathways analysis (all).png")
path_name <-paste0(dir.file,"/",path_name)
ggsave(path_name, kegg_pathways,width=1200, height =1000, dpi=150,units = "px")

protpaths_df <- list(pro_df=pro_df,
                     species=species,
                     kegg_all_pathways=kegg_all_pathways,
                     dir.file=dir.file,
                     kegg_mytheme=kegg_mytheme,
                     kegg_xytheme=kegg_xytheme
                     )

return(protpaths_df)

}

