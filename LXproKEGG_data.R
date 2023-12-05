
library(openxlsx)

protein_data_1 <- read.xlsx("protein_df.xlsx")
protein_data_2 <- read.xlsx("protein_ID.xlsx")

usethis::use_data(protein_data_1,overwrite = T)
usethis::use_data(protein_data_2,overwrite = T)

rm(list=ls())

data(protein_data_1)
data(protein_data_2)


