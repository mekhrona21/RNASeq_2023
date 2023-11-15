# installing packages

# installing packages----
install.packages("beepr")
install.packages('BiocManager')
install.packages('tidyverse')
install.packages("tximport")
install.packages('ensembldb')
install.packages('rhdf5')
install.packages('datapasta')
install.packages("cowplot")


BiocManager::install('tximport')
BiocManager::install('ensembldb')
BiocManager::install('rhdf5')
BiocManager::install("edgeR")

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.18") 

BiocManager::install("EnsDb.Tomato")
library(tidyverse) # provides access to Hadley Wickham's collection of R packages for data science
library(tximport) # package to import Kallisto results into R
library(ensembldb) # helps with ensembl
library(biomaRt) # annotations
library(edgeR) # package for differential expression analysis (used for the DGEList object and for normalization methods)
library(matrixStats)
library(cowplot)
library(ensembldb)
library(readr)
library(EnsDb.Tomato) # as this package is not available for our R version, we are using different way to load data

# Loading the plants_mart for Ensembl plant genomes
myMart <- useMart(biomart = "plants_mart", host = "https://plants.ensembl.org")

listDatasets(mart = myMart) # to check if the Solanum lycopersicum is there


# Load the Ensembl annotations for Solanum lycopersicum with the correct dataset name
tomato.anno <- useMart(biomart = "plants_mart", dataset = "slycopersicum_eg_gene", host = "https://plants.ensembl.org")

#formatting
Tx.cof <- getBM(attributes=c('ensembl_transcript_id',
                             'ensembl_gene_id', 'description'),
                mart = tomato.anno)
Tx.cof <- as_tibble(Tx.cof) 
Tx.cof <- dplyr::rename(Tx.cof, target_id = ensembl_transcript_id, gene_name = ensembl_gene_id)
Tx.cof <- dplyr::select(Tx.cof, "target_id", "gene_name")

# import a study design file
targets <- read_tsv("~/Desktop/Kallisto/studydesign.txt")

# Replace "your/directory/path" with the actual path to your sample directories
setwd("~/Desktop/Kallisto")

# Construct paths relative to the working directory
path <- file.path(targets$sample, "abundance.tsv", fsep = .Platform$file.sep)

file.exists(path) # check to see if path worked, should be TRUE before proceeding
Txi_gene <- tximport(path, 
                     type = "kallisto", 
                     tx2gene = Tx.cof, # mapps of transcript IDs to gene IDs 
                     txOut = TRUE, # data represented at gene level instead of at transcript level
                     countsFromAbundance = "lengthScaledTPM",
                     ignoreTxVersion = TRUE) 

