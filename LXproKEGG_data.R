
library(openxlsx)

protein_list <- read.xlsx("protein_list.xlsx")
protein_data <- read.xlsx("protein_data.xlsx")

usethis::use_data(protein_list,overwrite = T)
usethis::use_data(protein_data,overwrite = T)

rm(list=ls())

data(protein_list)
data(protein_data)



