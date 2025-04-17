library(igraph)
library(bnlearn)
library(pROC)
path_input = '/home/seongwonhwang/Desktop/projects/git/Node_ablation_practice/GRN/data/'

######## 
network <- read.delim(file.path(path_input, "NodeAb13_control_Mesoderm_all_possible_edges.txt"), header=F)
network <- network[!duplicated(network), ]
colnames(network) <- c("from", "to")
r <- igraph::graph_from_edgelist(as.matrix(network), directed = TRUE)
ig_network <- igraph::simplify(r,
  remove.loops = FALSE,
  remove.multiple = TRUE
)
######## read expression data -> data
expt <- read.delim(file.path(path_input, "NodeAb13_control_Mesoderm_embeddings_for_training.txt"), header=F)
genenames <- read.delim(file.path(path_input, "NodeAb13_control_Mesoderm_embeddings_gene_for_training.txt"), header=F)[,1]
rownames(expt) <- genenames
colnames(expt) <- paste0('cell', 1:ncol(expt))
expt_scaled <- as.data.frame(scale(t(expt)))

g_a <- igraph::graph_from_data_frame(expand.grid(genenames, genenames))

df_auc = NULL
for (rng_seed in 1:20){
  set.seed(rng_seed)
  
  n_edges_true <- length(igraph::E(ig_network))
  n_true = n_edges_true/2
  edges_true <- sample(attr(igraph::E(ig_network), "vnames"), n_true)
  edges_false <- sample(setdiff(attr(igraph::E(g_a), "vnames"), attr(igraph::E(ig_network), "vnames")), n_true *3)
  edges_all <- c(edges_true, edges_false)

  blacklist <-
    igraph::as_edgelist(
      igraph::delete_edges(g_a, which(attr(igraph::E(g_a), "vnames") %in% edges_all))
    )

  edge_strength <- bnlearn::boot.strength(
    data = expt_scaled,
    m = floor(0.20 * nrow(expt_scaled)),
    R = 200,
    debug = FALSE,
    algorithm = "hc",
    algorithm.args = list(
      score = "ebic-g",
      gamma = 0.5,
      blacklist = blacklist
    )
  )
  edge_strength <- edge_strength[match(
    edges_all,
    paste(edge_strength$from, edge_strength$to, sep = "|")
  ), ]

  edge_strength$true_edge = paste(edge_strength$from, edge_strength$to, sep = "|") %in% edges_true
  roc_object <- with(edge_strength, roc(true_edge, strength))
  df_auc = rbind(df_auc, data.frame(seed = rng_seed, auc = auc(roc_object)))
}
write.table(df_auc, 'eBIC_auc.txt', quote=F, row.names=F, col.names = T , sep='\t')