
LXproKEGG <- function(protein_data){

#----------安装及引用相关普通R包-----------------------------------------
pack_install <- function(){
all_packages <- data.frame(installed.packages()) # 查看已安装的R包

pack <- c("devtools","BiocManager","ggnewscale","R.utils", "ggtext",#需要安装的R包
          "roxygen2","xfun", "ggsci","openxlsx","dplyr","psych",
          "ggplot2","ggrepel","RColorBrewer", "ggthemes","rticles",
          "grid","patchwork","Hmisc","pak")

is_pack <- pack[!pack %in% all_packages$Package] # 筛选出未安装的R包

fun_install <- function(x) {install.packages(x,update = F,ask = F)} #安装R包函数
sapply(is_pack,fun_install)  #用sapply批量安装R包

# 批量library
fun_library <- function(x){library(x, character.only = T)}
sapply(pack,fun_library)

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
sapply(is_Biopack,fun_Bioinstall)  #用sapply批量安装R包

# 批量library
fun_library <- function(x){library(x, character.only = T)}
sapply(Biopack,fun_library)

}

pack_install()

#--------------kegg all pathways---------------------#
prot_path_all <- protpaths(protein_data)

#------------kegg metabolism pathways----------------#
meta_path <- prot_meta_path(prot_path_all)

#------------kegg signaling pathways----------------#
signal_path <- prot_signal_path(prot_path_all)

#--------------up-regulated pathways----------------#
up_path <- prot_up_path(prot_path_all)

#--------------down-regulated pathways--------------#
down_path <- prot_down_path(prot_path_all)

#-----------------up and down-regulation------------#
up_down_path <- up_down_path(prot_path_all,up_path,down_path)

#-------------GO analysis--------------------------#
prot_gO <-prot_GO(prot_path_all)


#--------------------------------------------------
print("The results can be seen in the folder of <analysis results>")

kegg_pathways


}




