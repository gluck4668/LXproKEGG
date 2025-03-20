
#--------------up-regulated pathways------------------------#

prot_up_path <- function(prot_path_all){

#----基本数据------------------------
  pro_df <- prot_path_all$pro_df
  species <- prot_path_all$species
  kegg_all_pathways <- prot_path_all$kegg_all_pathways
  dir.file <- prot_path_all$dir.file
  kegg_mytheme <- prot_path_all$kegg_mytheme
  kegg_xytheme <- prot_path_all$kegg_xytheme

#-------------------------------------
if(any(pro_df$FC>1)){

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

kegg_up <- kegg_up@result %>% na.omit()
names(kegg_up) <- gsub("Gene","Protein",names(kegg_up),ignore.case = T)
#-----把GeneRatio分数字符串变小数-------
    #用eval(str2expression(i))函数
str_eval <- function(i){eval(str2expression(i))}
kegg_up$ProteinRatio <- sapply(kegg_up$ProteinRatio,str_eval)

#-----拆分Description-------------------
kegg_up$Description <- str_extract(kegg_up$Description,".*(?= -)") %>% trimws()

#----导出数据----------------------------
file_path_name <- paste(species,"Up-regulated Protein KEGG pathways_data.xlsx")
file_path_name <-paste0(dir.file,"/",file_path_name)
write.xlsx(kegg_up,file_path_name)

row_n <- nrow(kegg_up)

if(row_n>=30){n_up=30
      title_text <- c("Top 30 Up-regulated Protein KEGG Pathways") }else
      { n_up=row_n
      title_text <- c("Up-regulated Protein KEGG Pathways")}

kegg_up_pathways <- ggplot(kegg_up[1:n_up,])+
        geom_point(aes(x=ProteinRatio,
                       y=fct_reorder(Description,ProteinRatio),
                       color=-log10(pvalue),size=Count))+
        scale_color_gradient2(midpoint = 1,low = "blue",mid = "#ffcc00",high ="#dc143c" )+
        labs(x = 'ProteinRatio', y = '',title=title_text)+
        kegg_mytheme+kegg_xytheme

kegg_up_pathways

path_name <- paste(species,"Up-regulated Protein KEGG pathways.png")
path_name <-paste0(dir.file,"/",path_name)
ggsave(path_name, kegg_up_pathways,width=1200, height =1000, dpi=150,units = "px")

prot_up_df <- list(kegg_up=kegg_up)

return(prot_up_df)


   }

 }else
    print("There are no up-regulated proteins.")


}
