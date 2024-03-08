
prot_GO <- function(prot_path_all){

#----基本数据------------------------
  pro_df <- prot_path_all$pro_df
  species <- prot_path_all$species
  kegg_all_pathways <- prot_path_all$kegg_all_pathways
  dir.file <- prot_path_all$dir.file
  kegg_mytheme <- prot_path_all$kegg_mytheme
  kegg_xytheme <- prot_path_all$kegg_xytheme

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


}
