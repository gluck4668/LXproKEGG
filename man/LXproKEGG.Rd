\name{LXproKEGG}
\alias{LXproKEGG}
\title{Protein KEGG analysis}
\usage{
LXproKEGG(protein_data)
}
\description{
LXproKEGG can be used for protein (a uniprot format) KEGG enrichment analysis.
}
\examples{
install.packages("devtools")
library(devtools)

install_github("gluck4668/LXproKEGG")

library(LXproKEGG)

??LXproKEGG
#---------------------
data(protein_data_1)
data(protein_data_2)

#--------------------

rm(list=ls())

protein_data="protein_data.xlsx" # The data should be a list of protein UNIPROT, or three list including Protein_Uniprot, FC and pvalue

LXproKEGG(protein_data)
}
