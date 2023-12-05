
LXproKEGG <- function(protein_data){

#----------安装及引用相关普通R包-----------------------------------------
all_packages <- data.frame(installed.packages()) # 查看已安装的R包

pack <- c("devtools","BiocManager","ggnewscale","R.utils", "ggtext",#需要安装的R包
          "roxygen2","xfun", "ggsci","openxlsx","dplyr","psych",
          "ggplot2","ggrepel","RColorBrewer", "ggthemes","rticles",
          "grid","patchwork","Hmisc","pak")

is_pack <- pack[!pack %in% all_packages$Package] # 筛选出未安装的R包

fun_install <- function(x) {install.packages(x,update = F,ask = F)} #安装R包函数
sapply(is_pack,fun_install,simplify = T)  #用sapply批量安装R包

# 批量library
fun_library <- function(x){library(x, character.only = T)}
sapply(pack,fun_library,simplify = T)

#-----tidyverse比较难安装，需要pak函数来装------------------------------
if(!"tidyverse" %in% all_packages$Package)
  pak::pak("tidyverse/tidyverse")
library(tidyverse)

#----安装BiocManager相关R包----------------------------------------------
Biopack <- c("DOSE","clusterProfiler","do","enrichplot",
             "pathview","BiocParallel","GO.db","KEGGREST",
             "org.Hs.eg.db","org.Mm.eg.db","org.Rn.eg.db","purrr")
# human: "org.Hs.eg.db"
# mouse: "org.Mm.eg.db"
# rat: "org.Rn.eg.db"
# purrr： map函数
library(BiocManager)

is_Biopack <- Biopack[!Biopack %in% all_packages$Package] # 筛选出未安装的R包

fun_Bioinstall <- function(x) {BiocManager::install(x,update = F,ask = F)} #安装R包函数
sapply(is_Biopack,fun_Bioinstall,simplify = T)  #用sapply批量安装R包

# 批量library
fun_library <- function(x){library(x, character.only = T)}
sapply(Biopack,fun_library,simplify = T)

#-----------读取数据，并筛选出差异表达蛋白-------------------------
pro_df <- read.xlsx(protein_data)

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
pathways_kegg <- kegg_pro@result


#-----把GeneRatio分数字符串变小数-------
#Ratio <- str_split(pathways_kegg$GeneRatio,"/",simplify = T) %>% apply(2,as.numeric)

#Ratio <- Ratio[,1]/Ratio[,2]

str_num <- function(i){eval(str2expression(i))} # 把字符串当做命令运行 eval(str2expression(i))
Ratio <- sapply(pathways_kegg$GeneRatio,str_num,simplify = T)

pathways_kegg <- mutate(pathways_kegg,ProteinRatio=round(Ratio,4)) #在pathways_kegg中添加一列ProteinRatio


#-----拆分Description-------------------
if(species== "human")
  kegg_all_pathways <- tidyr::separate(pathways_kegg,Description,"Description"," - Homo",remove = T)

if(species== "mouse")
  kegg_all_pathways <- tidyr::separate(pathways_kegg,Description,"Description"," - Mus",remove = T)

if(species== "rat")
  kegg_all_pathways <- tidyr::separate(pathways_kegg,Description,"Description"," - Rat",remove = T)

colnames(kegg_all_pathways)[c(2,3,8)] <- c("Pathways","Protein/Ratio","ProteinID")

file_path_name <- paste(species,"KEGG pathways analysis (all)_data.xlsx")
file_path_name <-paste0(dir.file,"/",file_path_name)
write.xlsx(kegg_all_pathways,file_path_name)

row_n <- nrow(kegg_all_pathways)

if(row_n>=30)
 title_text <- c("Top 30 KEGG Enriched Pathways") else
   title_text <- c("KEGG Enriched Pathways")

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

if(row_n<30)
  path_n <- row_n else
    path_n <- 30

kegg_df <- kegg_all_pathways[1:path_n,]

#--------------kegg all pathways------------------------#

kegg_pathways <- ggplot(kegg_df)+
  geom_point(aes(x=ProteinRatio,
                 y=fct_reorder(Pathways,ProteinRatio),
                 color=-log10(pvalue),size=Count))+
  scale_color_gradient2(midpoint = 1,low = "blue",mid = "#ffcc00",high ="#dc143c" )+
  labs(x = 'ProteinRatio', y = '',title=title_text)+
  kegg_mytheme+kegg_xytheme

kegg_pathways

path_name <- paste(species,"KEGG pathways analysis (all).png")
path_name <-paste0(dir.file,"/",path_name)

ggsave(path_name, kegg_pathways,width=1200, height =1000, dpi=150,units = "px")


#------------kegg metabolism pathways----------------#

#筛选出含有"metabol"字段的行
kegg_meta <- dplyr::filter(kegg_all_pathways,grepl("metabol",kegg_all_pathways$Pathways,ignore.case = T))
meta_path <- grep("Metabolic pathways",kegg_all_pathways$Pathways,ignore.case = T)
kegg_meta <- kegg_meta[-c(meta_path),] #去掉总的Metabolic pathways

if(nrow(kegg_meta)>0)
{kegg_meta$Pathways <- factor(kegg_meta$Pathways,levels=kegg_meta$Pathways)

title_meta_text <-c("KEGG Metabolism Pathways")

kegg_meta_pathways <- ggplot(kegg_meta)+
  geom_point(aes(x= ProteinRatio,y=reorder(Pathways,ProteinRatio),
                 color=-log10(pvalue),size=Count))+
  scale_color_gradient2(midpoint = 1,low = "blue",mid = "#ffcc00",high ="#dc143c")+
  labs(x = 'ProteinRatio', y = '',title=title_meta_text)+
  kegg_mytheme+kegg_xytheme

kegg_meta_pathways

meta_path_name <- paste(species,"metabolism pathways analysis.png")
meta_path_name <-paste0(dir.file,"/",meta_path_name)

ggsave(meta_path_name, kegg_meta_pathways,width=1200, height =1000, dpi=150,units = "px")


kegg_meta_file <- dplyr::arrange(kegg_meta,desc(ProteinRatio))

meta_name <- paste(species,"metabolism pathways analysis_data.xlsx")
meta_name <-paste0(dir.file,"/",meta_name)

write.xlsx(kegg_meta_file,meta_name)} else
print("There is no metabolism pathway!")


#------------kegg signaling pathways----------------#

#筛选出含有"signal"字段的行
kegg_signal <- dplyr::filter(kegg_all_pathways,grepl("signal",kegg_all_pathways$Pathways,ignore.case = T))

if(nrow(kegg_signal)>0)
{kegg_signal$Pathways <- factor(kegg_signal$Pathways,levels=kegg_signal$Pathways)

title_signal_text <- c("KEGG Signaling Pathways")

kegg_signal_pathways <- ggplot(kegg_signal)+
  geom_point(aes(x=ProteinRatio,y=reorder(Pathways,ProteinRatio),
                 color=-log10(pvalue),size=Count))+
  scale_color_gradient2(midpoint = 1,low = "blue",mid = "#ffcc00",high ="#dc143c")+
  labs(x = 'ProteinRatio', y = '',title=title_signal_text)+
  kegg_mytheme+kegg_xytheme

kegg_signal_pathways

sign_path_name <- paste(species,"signaling pathways analysis.png")
sign_path_name <-paste0(dir.file,"/",sign_path_name)

ggsave(sign_path_name, kegg_signal_pathways,width=1200, height =1000, dpi=150,units = "px")


kegg_signal_file <- dplyr::arrange(kegg_signal,desc(ProteinRatio))

sign_name <- paste(species,"signaling pathways analysis_data.xlsx")
sign_name <-paste0(dir.file,"/",sign_name)

write.xlsx(kegg_signal_file,sign_name)} else
  print("There is no signaling pathway!")


#--------------up-regulated pathways------------------------#
if(ncol(pro_df)>2){

pro_up <- dplyr::filter(pro_df,FC>1)

organism <- dplyr::case_when (species== "human" ~ "hsa",
                            species== "mouse" ~ "mmu",
                            species== "rat" ~ "rno"
                          )

kegg_up <- enrichKEGG(pro_up[,1], organism =organism,
                       #universe,
                       keyType = 'uniprot',  # 蛋白选择"uniprot" （如果是基因，则选择"kegg"）
                       pvalueCutoff = 0.05,
                       pAdjustMethod = 'BH',
                       qvalueCutoff = 0.2,
                       minGSSize = 3,
                       maxGSSize = 3500,
                       use_internal_data = F)

kegg_up <- kegg_up@result

#-----把GeneRatio分数字符串变小数-------
Ratio <- str_split(kegg_up$GeneRatio,"/",simplify = T) %>%
  apply(2,as.numeric)

Ratio <- Ratio[,1]/Ratio[,2]

kegg_up <- mutate(kegg_up,ProteinRatio=round(Ratio,4)) #在kegg_up中添加一列ProteinRatio

#-----拆分Description-------------------
if(species== "human")
  kegg_up <- tidyr::separate(kegg_up,Description,"Description"," - Homo",remove = T)

if(species== "mouse")
  kegg_up <- tidyr::separate(kegg_up,Description,"Description"," - Mus",remove = T)

if(species== "rat")
  kegg_up <- tidyr::separate(kegg_up,Description,"Description"," - Rat",remove = T)

colnames(kegg_up)[c(2,3,8)] <- c("Pathways","Protein/Ratio","ProteinID")

file_path_name <- paste(species,"Up-regulated KEGG pathways analysis_data.xlsx")
file_path_name <-paste0(dir.file,"/",file_path_name)
write.xlsx(kegg_up,file_path_name)

row_n <- nrow(kegg_up)

if(row_n>=30)
  title_text <- c("Top 30 Up-regulated KEGG Enriched Pathways") else
    title_text <- c("Up-regulated KEGG Enriched Pathways")

if(row_n<30)
  path_n <- row_n else
    path_n <- 30

kegg_up <- kegg_up[1:path_n,]

kegg_up_pathways <- ggplot(kegg_up)+
  geom_point(aes(x=ProteinRatio,
                 y=fct_reorder(Pathways,ProteinRatio),
                 color=-log10(pvalue),size=Count))+
  scale_color_gradient2(midpoint = 1,low = "blue",mid = "#ffcc00",high ="#dc143c" )+
  labs(x = 'ProteinRatio', y = '',title=title_text)+
  kegg_mytheme+kegg_xytheme

kegg_up_pathways

path_name <- paste(species,"Up-regulated KEGG pathways analysis.png")
path_name <-paste0(dir.file,"/",path_name)

ggsave(path_name, kegg_up_pathways,width=1200, height =1000, dpi=150,units = "px")

#--------------down-regulated pathways------------------------#
pro_down <- dplyr::filter(pro_df,FC<=1)

kegg_down <- enrichKEGG(pro_down[,1], organism =organism,
                      #universe,
                      keyType = 'uniprot',  # 蛋白选择"uniprot" （如果是基因，则选择"kegg"）
                      pvalueCutoff = 0.05,
                      pAdjustMethod = 'BH',
                      qvalueCutoff = 0.2,
                      minGSSize = 3,
                      maxGSSize = 3500,
                      use_internal_data = F)

kegg_down <- kegg_down@result

#-----把GeneRatio分数字符串变小数-------
Ratio <- str_split(kegg_down$GeneRatio,"/",simplify = T) %>%
  apply(2,as.numeric)

Ratio <- Ratio[,1]/Ratio[,2]

kegg_down <- mutate(kegg_down,ProteinRatio=round(Ratio,4)) #在kegg_down中添加一列ProteinRatio

#-----拆分Description-------------------
if(species== "human")
  kegg_down <- tidyr::separate(kegg_down,Description,"Description"," - Homo",remove = T)

if(species== "mouse")
  kegg_down <- tidyr::separate(kegg_down,Description,"Description"," - Mus",remove = T)

if(species== "rat")
  kegg_down <- tidyr::separate(kegg_down,Description,"Description"," - Rat",remove = T)

colnames(kegg_down)[c(2,3,8)] <- c("Pathways","Protein/Ratio","ProteinID")

file_path_name <- paste(species,"Down-regulated KEGG pathways analysis_data.xlsx")
file_path_name <-paste0(dir.file,"/",file_path_name)
write.xlsx(kegg_down,file_path_name)

row_n <- nrow(kegg_down)

if(row_n>=30)
  title_text <- c("Top 30 Down-regulated KEGG Enriched Pathways") else
    title_text <- c("Down-regulated KEGG Enriched Pathways")

if(row_n<30)
  path_n <- row_n else
    path_n <- 30

kegg_down <- kegg_down[1:path_n,]

kegg_down_pathways <- ggplot(kegg_down)+
  geom_point(aes(x=ProteinRatio,
                 y=fct_reorder(Pathways,ProteinRatio),
                 color=-log10(pvalue),size=Count))+
  scale_color_gradient2(midpoint = 1,low = "blue",mid = "#ffcc00",high ="#dc143c" )+
  labs(x = 'ProteinRatio', y = '',title=title_text)+
  kegg_mytheme+kegg_xytheme

kegg_down_pathways

path_name <- paste(species,"down-regulated KEGG pathways analysis.png")
path_name <-paste0(dir.file,"/",path_name)

ggsave(path_name, kegg_down_pathways,width=1200, height =1000, dpi=150,units = "px")


#-----------------up and down-regulation--------------

if(length(kegg_up)>0 & length(kegg_down)>0){
  meta_up_n <- grep("Metabolic pathways",kegg_up$Pathways,ignore.case = T) %>% as.numeric()
  if(length(meta_up_n)>0){
    kegg_up_10 <- kegg_up[-meta_up_n,]}

  kegg_up_10 <- dplyr::arrange(kegg_up_10,pvalue)
  kegg_up_10 <- mutate(kegg_up_10,"log10pvalue"=-log10(pvalue))

  meta_down_n <- grep("Metabolic pathways",kegg_down$Pathways,ignore.case = T) %>% as.numeric()

  if(length(meta_down_n)>0){
    kegg_down_10 <- kegg_down[-meta_down_n,]}

  kegg_down_10 <- dplyr::arrange(kegg_down_10,pvalue)
  kegg_down_10 <- mutate(kegg_down_10,"log10pvalue"=-log10(pvalue))
  kegg_down_10$log10pvalue <- 0-kegg_down_10$log10pvalue

  if(nrow(kegg_up_10)>10)
    k_up <- kegg_up_10[1:10,] else
      k_up <- kegg_up_10

  if(nrow(kegg_down_10)>10)
    k_down <- kegg_down_10[1:10,] else
      k_down <- kegg_down_10

  kegg_10 <- rbind(k_up,k_down)
  kegg_10$Pathways <- trimws(kegg_10$Pathways)

  Type <- case_when(kegg_10$log10pvalue>=0 ~"up",
                    kegg_10$log10pvalue<0 ~"down")

  kegg_10 <- mutate(kegg_10,Type=Type)

  kegg_10$log10p_up <- purrr::map(c(1:nrow(kegg_10)),~{if(kegg_10[.x,ncol(kegg_10)]=="up")
    kegg_10[.x,ncol(kegg_10)+1]=round(kegg_10[.x,ncol(kegg_10)-1],2) else
      kegg_10[.x,ncol(kegg_10)+1]=""
  } )

  kegg_10$log10p_down <- purrr::map(c(1:nrow(kegg_10)),~{if(kegg_10[.x,ncol(kegg_10)-1]=="down")
    kegg_10[.x,ncol(kegg_10)+1]=0-round(kegg_10[.x,ncol(kegg_10)-2],2) else
      kegg_10[.x,ncol(kegg_10)+1]=""
  } )

  ymax <- round(max(kegg_up_10$log10pvalue)+5,0)
  ymin <- round(min(kegg_down_10$log10pvalue)-5,0)



  p0 <- ggplot(kegg_10,aes(x =reorder(Pathways,log10pvalue),y =log10pvalue,fill=Type))+
    geom_col()+
    theme_classic()+
    ylim(ymin,ymax)+
    coord_flip()+
    scale_fill_manual(values = c("#2f73bb","#ae4531"))+ # 指定颜色
    geom_segment(aes(y=0,yend=0,x=0,xend=nrow(kegg_10)+0.8)) # 加一条坚线

  p0

  p0+geom_text(aes(label=format(round(log10pvalue, digits = 2), nsmall = 2)), nudge_y= -1, color="black", size =3)

  # 调整主题
  my_theme <- theme(
    legend.position = "none", #去除图例
    plot.title=element_text(hjust = 0.5), #标题居中
    axis.line.y = element_blank(), # 隐去y轴线
    axis.title.y = element_blank(), # 隐去y轴"reorder(pathway,value)"
    axis.ticks.y= element_blank(), # 隐去y轴短线
    axis.text.y = element_blank(), # 隐去y轴"pathway1....."
    axis.text = element_text(colour="black",size=14)
  )+
    theme(plot.title = element_text(size=16,hjust = 0.5))+
    theme(axis.title.x = element_text(colour="black",size=16))+
    theme(plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"))

  #添加箭头
  arrow_down <- geom_segment(aes(y=-0.2,yend=ymin, x=nrow(kegg_10)+1,xend=nrow(kegg_10)+1),
                             arrow=arrow(length = unit(0.2,"cm"),type="closed"),
                             linewidth=0.5)


  arrow_up <- geom_segment(aes(y=0.2,yend=ymax, x=nrow(kegg_10)+1,xend=nrow(kegg_10)+1),
                           arrow=arrow(length = unit(0.2,"cm"),type="closed"),
                           linewidth=0.5)

  #在柱上显示数值
  col_value_up <- geom_text(aes(label=log10p_up),
                            nudge_y=-1, # 和柱顶的距离
                            color="black", size =3)

  col_value_down <- geom_text(aes(label=log10p_down),
                              nudge_y=1, # 和柱顶的距离
                              color="black", size =3)


  # 添加labels
  label_up <- geom_text(data = kegg_10[which(kegg_10$log10pvalue>=0),], #筛选出value>=0，用which函数
                        aes(x=Pathways,y=0,label=trimws(Pathways)), # 添加
                        hjust=1, nudge_y =-0.2,  # hjust=1 是右对齐，nudge_y微调距离
                        size=4)

  label_down <- geom_text(data = kegg_10[which(kegg_10$log10pvalue<0),], #筛选出value>=0，用which函数
                          aes(x=Pathways,y=0,label=trimws(Pathways)), # 添加
                          hjust=0, nudge_y =0.2,  # hjust=0 是左对齐，nudge_y微调距离
                          size=4 )

  p01 <- p0+my_theme+
    label_up+label_down+
    arrow_down +
    arrow_up+
    annotate("text",x=nrow(kegg_10)+1.6,y=ymin+1,label="Down",)+
    annotate("text",x=nrow(kegg_10)+1.6,y=ymax-1,label="Up")+
    ylab("-log10(pvalue)")+
    ggtitle("KEGG pathways enrichment\n (up and down-regulation)")+
    scale_x_discrete(expand=expansion(add=c(0.5,2))) # 边距

  p01

  p02 <- p01+
    col_value_up+
    col_value_down

  p02

  path_10_01 <- paste(species,"KEGG enriched pathways (up and down) 01.png")
  path_10_01 <-paste0(dir.file,"/",path_10_01)
  ggsave(path_10_01,p01,width = 1500,height =1000,units = "px",dpi = 150)

  path_10_02 <- paste(species,"KEGG enriched pathways (up and down) 02.png")
  path_10_02 <-paste0(dir.file,"/",path_10_02)
  ggsave(path_10_02,p02,width = 1500,height =1000,units = "px",dpi = 150)

  file_path_name <- paste(species,"KEGG enriched pathways (up and down)_data.xlsx")
  file_path_name <-paste0(dir.file,"/",file_path_name)
  write.xlsx(kegg_10,file_path_name)

}


}

#-------------蛋白GO分析--------------------------------------
OrgDb <- dplyr::case_when ( species== "human" ~ "org.Hs.eg.db",
                            species== "mouse" ~ "org.Mm.eg.db",
                            species== "rat" ~ "org.Rn.eg.db"
                          )

# keytypes(org.Mm.eg.db) 查看org.Mm.eg.db包含的关键条目

go_enrich <- enrichGO(pro_df[,1],
                      OrgDb = OrgDb,
                      keyType = 'UNIPROT',
                      ont='ALL',
                      pAdjustMethod = 'BH',
                      pvalueCutoff = 0.05,
                      qvalueCutoff = 0.2)


go_result <- go_enrich@result

go_result$log10pvalue  <- round(-log10(go_result$pvalue),2)

go_result <- go_result[order(go_result$ONTOLOGY,-go_result$log10pvalue),]

GO_file <- paste(species,"GO enrichment (all)_data.xlsx")
GO_file <-paste0(dir.file,"/",GO_file)
write.xlsx(go_result,GO_file)

go_BP_all <- go_result[go_result$ONTOLOGY=="BP",]
go_BP_p <- go_BP_all[go_BP_all$pvalue<0.05 & go_BP_all$p.adjust<0.05,]

if(nrow(go_BP_p)>20)
  n_BP=20 else
    n_BP= nrow(go_BP_p)

go_CC_all <- go_result[go_result$ONTOLOGY=="CC",]
go_CC_p <- go_CC_all[go_CC_all$pvalue<0.05 & go_CC_all$p.adjust<0.05,]

if(nrow(go_CC_p)>20)
  n_CC=20 else
    n_CC= nrow(go_CC_p)

go_MF_all <- go_result[go_result$ONTOLOGY=="MF",]
go_MF_p <- go_MF_all[go_MF_all$pvalue<0.05 & go_MF_all$p.adjust<0.05,]

if(nrow(go_MF_p)>20)
  n_MF=20 else
    n_MF= nrow(go_MF_p)

go_BP <- go_BP_all[c(1:n_BP),]
go_CC <- go_CC_all[c(1:n_CC),]
go_MF <- go_MF_all[c(1:n_MF),]

go_df <- rbind(go_BP,go_CC,go_MF)
go_df <- na.omit(go_df)

go_df <- go_df[order(go_df$ONTOLOGY,-go_df$log10pvalue),]

title_size_go <- case_when(nrow(go_df)>50 ~12,
                           TRUE ~11)

xy_siz_go <- case_when(nrow(go_df)>50 ~9,
                       nrow(go_df)>30 ~10,
                       TRUE ~11)

go_mytheme<-theme_bw()+
  theme(text=element_text(family = "sans",face="bold",colour ="black",size =title_size_go),
        panel.grid = element_blank(),
        panel.border = element_rect(linewidth = 0.8,colour = "gray30"),
        axis.line = element_blank(),
        axis.ticks = element_line(linewidth = 0.6,colour = "gray30"),
        axis.ticks.length = unit(1.5,units = "mm"))+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"))


go_xytheme <-theme(axis.text.x = element_text(face="bold",color="black",size=xy_siz_go,angle =80,hjust=1))+
  theme(axis.text.y = element_text(face="bold",color="black",size=xy_siz_go))+
  theme(legend.text=element_text(face="bold",color="black",size=xy_siz_go))

color_text <- c(rep("#20b2aa",nrow(go_BP)),rep("#d2691e",nrow(go_CC)),rep("#6666cc",nrow(go_MF)))
textcolor_theme <- theme(axis.text.x = element_text(colour = color_text))

legend_theme01 <- theme(
  legend.title = element_blank(),
  legend.text = element_text(size = 10, face = "bold"),
  legend.direction = "horizontal",
  legend.position = c(0.5,0.9),
  legend.background = element_blank()
)

legend_theme02 <- theme(
  legend.title = element_blank(),
  legend.text = element_text(size = 10, face = "bold"),
  legend.direction = "horizontal",
  legend.position = "none",
  legend.background = element_blank()
)


height_y01 <- max(go_df$log10pvalue)+5

nn_row <- nrow(go_df) %>% as.numeric()

go_df$Term <- NA

for(i in 1:nn_row){
  if(nchar(go_df[i,3])>40)
    go_df[i,12] <- paste(str_sub(go_df[i,3],1,40),"...") else
      go_df[i,12] <- go_df[i,3]
}

go_df <- distinct (go_df,Term,.keep_all = T)

go_df$Term =factor(go_df$Term,levels=go_df$Term)

legend_text <- c("Biological process(BP)","Cellular component (CC)","Molecular function (MF)")

GO_plot <- ggplot(go_df, aes(x = Term,y =log10pvalue, fill = ONTOLOGY))

GO_enrich01 <- GO_plot+ geom_bar(position = "dodge",stat = "identity",width = 0.8)+
  scale_y_continuous(expand = c(0, 0),limits = c(0, height_y01))+
  scale_fill_manual(values = c("#20b2aa", "#d2691e","#6666cc"),label=legend_text)+
  labs(x = '', y = '-log10(pvalue)',title="GO enrichment")+
  go_mytheme+go_xytheme +textcolor_theme+legend_theme01
# + scale_x_discrete(labels=function(go_enrich) str_wrap(go_enrich, width=50)) #限定x轴字段宽度

GO_enrich01

GO_pic <- paste(species,"GO enrichment.png")
GO_pic <-paste0(dir.file,"/",GO_pic)

ggsave(GO_pic, GO_enrich01,width=1500, height =1000, dpi=150,units = "px")

GO_f01 <- paste(species,"GO enrichment_data.xlsx")
GO_f01 <-paste0(dir.file,"/",GO_f01)

write.xlsx(go_df,GO_f01)

#--------------------------------------
print("The results can be seen in the folder of <analysis results>")

kegg_pathways+GO_enrich01

}




