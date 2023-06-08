
library(openxlsx)

protein_data <- read.xlsx("protein_data.xlsx")

usethis::use_data(protein_data,overwrite = T)

rm(list=ls())

data(protein_data)


