
up_down_path <- function(prot_path_all,up_path,down_path){

#----基本数据------------------------
  pro_df <- prot_path_all$pro_df
  species <- prot_path_all$species
  dir.file <- prot_path_all$dir.file
  kegg_mytheme <- prot_path_all$kegg_mytheme
  kegg_xytheme <- prot_path_all$kegg_xytheme
  kegg_up <- up_path$kegg_up
  kegg_down <- down_path$kegg_down

#-----------------up and down-regulation------------#
if(any(pro_df$FC>1) & any(pro_df$FC<1)){  # 判断有没上调和下调protins

    meta_up_n <- grep("Metabolic pathways",kegg_up$Description,ignore.case = T) %>% as.numeric()
    if(length(meta_up_n)>0){
      kegg_up_10 <- kegg_up[-meta_up_n,]}

    kegg_up_10 <- dplyr::arrange(kegg_up_10,pvalue)
    kegg_up_10 <- mutate(kegg_up_10,"log10pvalue"=-log10(pvalue))

    meta_down_n <- grep("Metabolic pathways",kegg_down$Description,ignore.case = T) %>% as.numeric()
    if(length(meta_down_n)>0){
      kegg_down_10 <- kegg_down[-meta_down_n,]}

    kegg_down_10 <- dplyr::arrange(kegg_down_10,pvalue)
    kegg_down_10 <- mutate(kegg_down_10,"log10pvalue"=-log10(pvalue))
    kegg_down_10$log10pvalue <- 0-kegg_down_10$log10pvalue

    if(nrow(kegg_up_10)>=20)
      k_up <- kegg_up_10[1:20,] else
        k_up <- kegg_up_10

    if(nrow(kegg_down_10)>=20)
      k_down <- kegg_down_10[1:20,] else
        k_down <- kegg_down_10

    kegg_20 <- rbind(k_up,k_down)
    kegg_20$Description <- trimws(kegg_20$Description)

    Type <- case_when(kegg_20$log10pvalue>=0 ~"up",
                      kegg_20$log10pvalue<0 ~"down")

    kegg_20 <- mutate(kegg_20,Type=Type)

    kegg_20$log10p_up <- purrr::map(c(1:nrow(kegg_20)),~{if(kegg_20[.x,ncol(kegg_20)]=="up")
      kegg_20[.x,ncol(kegg_20)+1]=round(kegg_20[.x,ncol(kegg_20)-1],2) else
        kegg_20[.x,ncol(kegg_20)+1]=""
    } )

    kegg_20$log10p_down <- purrr::map(c(1:nrow(kegg_20)),~{if(kegg_20[.x,ncol(kegg_20)-1]=="down")
      kegg_20[.x,ncol(kegg_20)+1]=0-round(kegg_20[.x,ncol(kegg_20)-2],2) else
        kegg_20[.x,ncol(kegg_20)+1]=""
    } )

    ymax <- round(max(kegg_up_10$log10pvalue),0)
    ymin <- round(min(kegg_down_10$log10pvalue),0)

   hight=table(duplicated(kegg_20$Description))["FALSE"] # 坚线高度

    p0 <- ggplot(kegg_20,aes(x =reorder(Description,log10pvalue),y =log10pvalue,fill=Type))+
      geom_col()+
      theme_classic()+
      ylim(ymin,ymax)+
      coord_flip()+
      scale_fill_manual(values = c("#2f73bb","#ae4531"))+ # 指定颜色
      geom_segment(aes(y=0,yend=0,x=0,xend=(hight+0.8))) # 加一条坚线

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
    arrow_down <- geom_segment(aes(y=-0.2,yend=ymin, x=hight+1,xend=hight+1),
                               arrow=arrow(length = unit(0.2,"cm"),type="closed"),
                               linewidth=0.5)


    arrow_up <- geom_segment(aes(y=0.2,yend=ymax, x=hight+1,xend=hight+1),
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
    label_up <- geom_text(data = kegg_20[which(kegg_20$log10pvalue>=0),], #筛选出value>=0，用which函数
                          aes(x=Description,y=0,label=trimws(Description)), # 添加
                          hjust=1, nudge_y =-0.2,  # hjust=1 是右对齐，nudge_y微调距离
                          size=4)

    label_down <- geom_text(data = kegg_20[which(kegg_20$log10pvalue<0),], #筛选出value>=0，用which函数
                            aes(x=Description,y=0,label=trimws(Description)), # 添加
                            hjust=0, nudge_y =0.2,  # hjust=0 是左对齐，nudge_y微调距离
                            size=4 )

    p01 <- p0+my_theme+
      label_up+label_down+
      arrow_down +
      arrow_up+
      annotate("text",x=hight+1.6,y=ymin+1,label="Down",)+
      annotate("text",x=hight+1.6,y=ymax-1,label="Up")+
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
    ggsave(path_10_01,p0,width = 1500,height =1000,units = "px",dpi = 150)

    path_10_02 <- paste(species,"KEGG enriched pathways (up and down) 02.png")
    path_10_02 <-paste0(dir.file,"/",path_10_02)
    ggsave(path_10_02,p01,width = 1500,height =1000,units = "px",dpi = 150)

    path_10_03 <- paste(species,"KEGG enriched pathways (up and down) 03.png")
    path_10_03 <-paste0(dir.file,"/",path_10_03)
    ggsave(path_10_03,p02,width = 1500,height =1000,units = "px",dpi = 150)

    file_path_name <- paste(species,"KEGG enriched pathways (up and down)_data.xlsx")
    file_path_name <-paste0(dir.file,"/",file_path_name)
    write.xlsx(kegg_20,file_path_name)


    }


}
