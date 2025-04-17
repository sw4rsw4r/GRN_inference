#!/usr/bin/env Rscript

# the gene annotation in downloaded FANTOM files is based on Ensembl 92 (Gencode human  v28) for humans and Ensembl 93 (Gencode mouse vm18) for mouse
# wget http://humantfs.ccbr.utoronto.ca/download/v_1.01/TFs_Ensembl_v_1.01.txt -O data/Lambert_TFlist_ENSG.txt
library(biomaRt)
options(stringsAsFactors = F)

# define biomart objects, current and archived
mart.ensembl98 <- useMart(host = "sep2019.archive.ensembl.org", biomart = "ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl")
fname <- "data/Lambert_TFlist_ENSG.txt"
ENSG_IDs <- read.table(fname, sep = "\t", check.names = F, header = F)$V1

# query biomart for the ensembl_gene_id and get the corresponding gene_names
ensembl98.results <- getBM(attributes = c("ensembl_gene_id", "external_gene_name"), filters = "ensembl_gene_id", values = ENSG_IDs, mart = mart.ensembl98)
TF_genenames <- unique(toupper(ensembl98.results$external_gene_name))

write.table(TF_genenames, "data/human_TFs_v2.7_mapped_to_Ensembl98.txt", quote = F, row.names = F, col.names = F, sep = "\n")
