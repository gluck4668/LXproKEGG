install.packages("devtools")
library(devtools)

install_github("gluck4668/LXproKEGG")

library(LXproKEGG)

??LXproKEGG
#---------------------
data(protein_data)

#--------------------

rm(list=ls())

protein_data="protein_data.xlsx" # The data should be a list of protein UNIPROT, or three list including Protein_Uniprot, FC and pvalue

LXproKEGG(protein_data)
