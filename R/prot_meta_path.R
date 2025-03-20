
#------------kegg metabolism pathways----------------#

prot_meta_path <- function(prot_path_all){

#----基本数据------------------------
  species <- prot_path_all$species
  kegg_all_pathways <- prot_path_all$kegg_all_pathways
  dir.file <- prot_path_all$dir.file
  kegg_mytheme <- prot_path_all$kegg_mytheme
  kegg_xytheme <- prot_path_all$kegg_xytheme

#---筛选出含有"metabol"字段的行-----
  if(any(grepl("category",colnames(kegg_all_pathways),ignore.case = T)==TRUE))
    kegg_meta <- dplyr::filter(kegg_all_pathways,grepl("metabol",kegg_all_pathways$category,ignore.case = T)) else
    kegg_meta <- dplyr::filter(kegg_all_pathways,grepl("metabol",kegg_all_pathways$Description,ignore.case = T))

  if(nrow(kegg_meta)>0) {
    kegg_meta_file <- dplyr::arrange(kegg_meta,desc(RichFactor))
    meta_name <- paste(species,"protein metabolism pathways_data.xlsx")
    meta_name <-paste0(dir.file,"/",meta_name)
    write.xlsx(kegg_meta_file,meta_name)

    if(any(kegg_all_pathways$Description=="Metabolic pathways")){
      meta_path <- grep("Metabolic pathways",kegg_all_pathways$Description,ignore.case = T)
      kegg_meta <- kegg_meta[-c(meta_path),] #去掉总的Metabolic pathways
    }

    kegg_meta <- distinct(kegg_meta,Description,.keep_all = T)

    kegg_meta$Description <- factor(kegg_meta$Description,levels=kegg_meta$Description)

    if(nrow(kegg_meta)>=30)
    {meta_pathways <- kegg_meta[1:30,]
    title_meta_text <-c("Top 30 Protein KEGG Metabolism Pathways") } else
    {meta_pathways <- kegg_meta
    title_meta_text <-c("Protein KEGG Metabolism Pathways") }

    kegg_meta_pathways <- ggplot(meta_pathways)+
      geom_point(aes(x= RichFactor,y=fct_reorder(Description,-log10(pvalue)),
                     color=-log10(pvalue),size=Count))+
      scale_color_gradient2(midpoint = 1,low = "blue",mid = "#ffcc00",high ="#dc143c")+
      labs(x = 'Rich Factor', y = '',title=title_meta_text)+
      kegg_mytheme+kegg_xytheme

    kegg_meta_pathways

    meta_path_name <- paste(species,"protein metabolism pathways.png")
    meta_path_name <-paste0(dir.file,"/",meta_path_name)
    ggsave(meta_path_name, kegg_meta_pathways,width=1200, height =1000, dpi=150,units = "px")

  } else
    print("There is no metabolism pathway!")

}
