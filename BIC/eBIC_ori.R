library(igraph)
library(bnlearn)

path_input <- "/home/seongwonhwang/Desktop/projects/mogrify/Statistical\ Consulting/"

######## read network data -> mog_network
network <- read.delim(file.path(path_input, "BIC/data/networks_anonymize.txt"))
network <- network[!duplicated(network), ]
colnames(network) <- c("TF", "target")

r <- igraph::graph_from_edgelist(as.matrix(network), directed = TRUE)
# remove (if-any) multiple edges, self-loops and nodes named as '1'
mog_network <- igraph::simplify(r,
  # Keep loops as we do in Mogrify
  remove.loops = FALSE,
  # Remove multiple edges (sanity-check)
  remove.multiple = TRUE
)
######## read expression data -> data
# Log-normalized counts
expt <- read.delim(file.path(path_input, "Bayesian_DE/iterative_test/iterative_test/TF_experiment_expression_matrix.gz"), row.names = 1)
expt$id <- NULL

cell_type <- c("FSK", "TDF1", "TDF2", "TDF3", "TDF4")
coldata <- read.delim(file.path(path_input, "Bayesian_DE/iterative_test/iterative_test/TF_experiment_metadata.gz"))
samples <- which(coldata[, "label.main"] %in% cell_type)
lc <- expt[, samples]

data <- scale(t(lc))
data <- as.data.frame(data[, !apply(is.na(data), 2, all)])
#####################################

# Restrict the scoring space to the genes we have in our scRNA-seq dataset.
candidates <- sort(intersect(
  igraph::V(mog_network)$name,
  colnames(data)
))

# Since the Hill-Climbing method is agnostic to node classes (e.g. TF or target),
# we calculate a second graph with all potential edges ('universe') between filtered nodes.
g_a <- igraph::graph_from_data_frame(expand.grid(
  candidates,
  candidates
))
g_cand <- igraph::delete.vertices(
  mog_network,
  which(!igraph::V(mog_network)$name %in% candidates)
)

# We eliminate the edges we want to learn about from this 'universe', so we could
# use it as 'blacklist' to the Hill-Climbing and increase accuracy while reduce computational cost.
blacklist <-
  igraph::as_edgelist(igraph::delete_edges(g_a, which(
    attr(igraph::E(g_a), "vnames") %in% attr(igraph::E(g_cand), "vnames")
  )))


# To calculate the strength of every edge, we bootstrap 'iters' number of times
# along the dataset, and learn/score our 'mog_network' graph using the
# Hill-Climbing approach with eBIC score for Gaussian data (scaled data)
# (see https://www.bnlearn.com/documentation/man/arc.strength.html).
edge_strength <- bnlearn::boot.strength(
  # Data to learn/score from
  data = data[, candidates],
  # Number of cells/samples to keep at every iteration
  m = floor(0.20 * nrow(data)),
  # Number of iterations (bootstrap). NOTE: TRY WITH iters>=200
  R = 200,
  # Verbose? (it is very detailed)
  debug = FALSE,
  # Hill-Climbing algorithm
  algorithm = "hc",
  # Parameters to Hill-Climbing approach
  algorithm.args = list(
    # By default, Extended Bayesian Information Criterion for Gaussian data.
    # Others might be: 'bic-g', 'bge', etc.
    # (see more about this in the Gaussian Bayesian Networks section in
    # https://www.bnlearn.com/documentation/man/network.scores.html).
    score = "ebic-g",
    # Default gamma is 0.5, which works better for
    # most large networks (see https://www.bnlearn.com/simulations/ebic-default-gamma/
    # for more information).
    gamma = 0.5,
    # Edges not to learn.
    blacklist = blacklist
  )
)

# Match to original edges.
edge_strength <-
  edge_strength[match(
    attr(igraph::E(g_cand), "vname"),
    paste(edge_strength$from, edge_strength$to, sep = "|")
  ), ]

write.table(edge_strength, "edge_strength.txt", quote = F, row.names = F, col.names = T, sep = "\t")
