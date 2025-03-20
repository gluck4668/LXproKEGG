install.packages("devtools")
library(devtools)

install_github("gluck4668/LXproKEGG")

library(LXproKEGG)

??LXproKEGG
#---------------------
data(protein_list)
data(protein_data)

#--------------------

rm(list=ls())

# devtools::load_all()

protein_data="protein_list.xlsx" # The data should be a list of protein UNIPROT,
                                 # or three list including Protein_Uniprot, FC and pvalue

LXproKEGG(protein_data)


