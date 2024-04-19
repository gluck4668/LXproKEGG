
#------------kegg signaling pathways----------------#

prot_signal_path <- function(prot_path_all){

#----基本数据------------------------
  species <- prot_path_all$species
  kegg_all_pathways <- prot_path_all$kegg_all_pathways
  dir.file <- prot_path_all$dir.file
  kegg_mytheme <- prot_path_all$kegg_mytheme
  kegg_xytheme <- prot_path_all$kegg_xytheme

#筛选出含有"signal"字段的行
kegg_signal <- dplyr::filter(kegg_all_pathways,grepl("signal",kegg_all_pathways$Description,ignore.case = T))

if(nrow(kegg_signal)>0)
{kegg_signal_file <- dplyr::arrange(kegg_signal,desc(ProteinRatio))
sign_name <- paste(species,"protein signaling pathways_data.xlsx")
sign_name <-paste0(dir.file,"/",sign_name)
write.xlsx(kegg_signal_file,sign_name)

kegg_signal$Description <- factor(kegg_signal$Description,levels=kegg_signal$Description)

if(nrow(kegg_signal)>=30)
{ nn=30
title_signal_text <- c("Top 30 Protein KEGG Signaling Pathways")
} else
{ nn=nrow(kegg_signal)
title_signal_text <- c("Protein KEGG Signaling Pathways")
}

kegg_signal_pathways <- ggplot(kegg_signal[1:nn,])+
  geom_point(aes(x=ProteinRatio,y=reorder(Description,ProteinRatio),
                 color=-log10(pvalue),size=Count))+
  scale_color_gradient2(midpoint = 1,low = "blue",mid = "#ffcc00",high ="#dc143c")+
  labs(x = 'ProteinRatio', y = '',title=title_signal_text)+
  kegg_mytheme+kegg_xytheme

kegg_signal_pathways

sign_path_name <- paste(species,"protein signaling pathways analysis.png")
sign_path_name <-paste0(dir.file,"/",sign_path_name)
ggsave(sign_path_name, kegg_signal_pathways,width=1200, height =1000, dpi=150,units = "px")

} else
  print("There is no signaling pathway!")

}
