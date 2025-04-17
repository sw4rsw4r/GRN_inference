library(GENIE3)
library(igraph)

# exprMatr <- matrix(sample(1:10, 100, replace=TRUE), nrow=20)
# rownames(exprMatr) <- paste("Gene", 1:20, sep="")
# colnames(exprMatr) <- paste("Sample", 1:5, sep="")
# TFs_all <- c("Gene2", "Gene4", "Gene7")
convid <- "dermal.fibroblast_myoblast"
exprMatr <- as.matrix(read.delim(paste0("../fantom5/data/lognorm_", convid, ".txt")))
head(exprMatr)
df_TFs <- read.table("../network/data/human_TFs_v2.7_mapped_to_Ensembl98.txt")
TFs_all <- intersect(df_TFs[, 1], rownames(exprMatr))

set.seed(1) # For reproducibility of results/
weightMat <- GENIE3(exprMatr, regulators = TFs_all, treeMethod = "ET", nTrees = 1000, nCores = 4, verbose = T)

linkList <- getLinkList(weightMat)
linkList_selected <- subset(linkList, weight > quantile(linkList$weight, .7))
write.table(linkList_selected, paste0("data/GENIE3_", convid, ".txt"), quote = F, row.names = F, col.names = F, sep = "\t")

STRING <- read.table('../network/data/human_string_v2.7_withdirection_new_noexpfilter')
STRING$V4=NULL
colnames(STRING)= c('regulatoryGene', 'targetGene', 'weight')
linkList_selected2 <- subset(linkList, weight > quantile(linkList$weight, .8))
combined = rbind(STRING, linkList_selected2)
write.table(combined, paste0("data/GENIE3_STRING_", convid, ".txt"), quote = F, row.names = F, col.names = F, sep = "\t")


g1 <- igraph::graph_from_data_frame(linkList_selected, directed = T)
lst_eigen_centrality1 <- eigen_centrality(g1)$vector

c("SOX2", "NANOG", "POU5F1") %in% names(tail(sort(lst_eigen_centrality1[TFs_all]), 20))
