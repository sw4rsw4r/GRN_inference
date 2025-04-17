# setwd("/home/seongwonhwang/Desktop/projects/git/GRN_inference_practice/analysis")
source("plot_utils.R")

# File path where GCN output stored
path_input <- "/home/seongwonhwang/Desktop/projects/mogrify/Statistical\ Consulting/"
path_tf_and_reqdgenes <- "/home/seongwonhwang/Desktop/projects/GRN_in_general/PyG/data/Anonymized_tables/Anonymized_tables/"

# get all genes
input_expr <- file.path(path_input, "Bayesian_DE/iterative_test/iterative_test/TF_experiment_expression_matrix")
df_expr_full <- read.delim(input_expr, row.names = 1)
genes_all <- df_expr_full$id

# get all TFs
external_tfs <- grep("tf", genes_all, value = T)
tfs <- get_tfs(path_tf_and_reqdgenes)
tfs_all <- unique(c(external_tfs, tfs))

# Load graph predictions
list_predicted <- list()
for (TYPE in c("GCN", "graphSAGE", "GAT")) {
    if (TYPE == "GAT") {
        TEST_ID <- "TEST3"
        lst_rng_seed <- c("", "_rng111_neg2.0", "_rng123_neg2.0", "_rng1234_neg2.0", "_rng111_neg3.0", "_rng123_neg3.0", "_rng1234_neg3.0", "_rng111_neg5.0", "_rng123_neg5.0", "_rng1234_neg5.0", "_rng111_neg3.0_super", "_rng123_neg3.0_super", "_rng1234_neg3.0_super", "_rng111_neg3.0_linear", "_rng123_neg3.0_linear", "_rng1234_neg3.0_linear", "_rng111_neg3.0_suplinear", "_rng123_neg3.0_suplinear", "_rng1234_neg3.0_suplinear")
        for (rng in lst_rng_seed) {
            ID <- paste0(TEST_ID, rng)
            list_predicted[[paste0(TEST_ID, "_", TYPE, rng)]] <- read_graph_prediction(TYPE, ID)
        }
        lst_rng_seed <- c("_rng111_neg3.0_linear", "_rng123_neg3.0_linear", "_rng1234_neg3.0_linear")
        for (TEST_ID in c("TEST4", "TEST5", "TEST6")) {
            for (rng in lst_rng_seed) {
                ID <- paste0(TEST_ID, rng)
                if (TEST_ID == "TEST5") {
                    list_predicted[[paste0(TEST_ID, "_", TYPE, rng)]] <- read_graph_prediction(TYPE, ID, .9)
                } else {
                    list_predicted[[paste0(TEST_ID, "_", TYPE, rng)]] <- read_graph_prediction(TYPE, ID)
                }
            }
        }
    } else {
        list_predicted[[TYPE]] <- read_graph_prediction(TYPE, "TEST3")
    }
}
networks <- load_all_networks(list_predicted, tfs_all)

df_all_possible_edges <- expand.grid(tfs_all, genes_all)
all_g <- igraph::graph_from_edgelist(as.matrix(df_all_possible_edges[, 1:2]))

# Reference
reference <- attr(E(simplify(all_g)), "vnames")
refers <- function(v) factor(ifelse(reference %in% attr(E(v), "vnames"), "Edge", "Not"), levels = c("Edge", "Not"))
ref_dorothea <- refers(networks$Dorothea$ig)
ref_GTRD <- refers(networks$GTRD$ig)

results_dorothea <- lapply(networks, function(x) {
    caret::confusionMatrix(refers(x$ig), ref_dorothea, mode = "everything")
})
results_gtrd <- lapply(networks, function(x) {
    caret::confusionMatrix(refers(x$ig), ref_GTRD, mode = "everything")
})

# Represent metrics
TEST_ID <- "TEST6"
plot_f1(results_dorothea, paste0(TEST_ID, "_rng_boxplot"), "Dorothea")
plot_f1(results_gtrd, paste0(TEST_ID, "_rng_boxplot"), "GTRD")

plot_score_distribution(dir_network = paste0("../GAT/data/"), TEST_ID)

plot_scores_comp(dir_network = paste0("../GAT/data/"), TEST_ID, lst_rng_seed, type = "eigen_centrality", tfs_all)
plot_edge_scores_comp(dir_network = paste0("../GAT/data/"), TEST_ID, lst_rng_seed)


# plot number of edges
# plot_n_edges(networks)
test_n_rng("TEST11v2")

TESTID <- "TEST15v4"
df_combined_edge_scores <- combine_rng_edge_scores(TESTID, 50)

p_df <- NULL
for (rng_seed in 1:50) {
    fname_log <- paste0("../GAT/data/", TESTID, "_rng", rng_seed, "_neg3.0_linear.log")
    log <- readLines(fname_log)
    idx_best_model <- rev(grep("model saved", log))[1]
    log_best_model <- log[idx_best_model + 1]
    auc_test <- as.numeric(rev(strsplit(log_best_model, "Test AUC: ")[[1]])[1])
    p_df <- rbind(p_df, data.frame(TESTID, auc_test))
}

lst_edge_scores <- apply(df_combined_edge_scores, 1, median)
df_net <- matrix(unlist(strsplit(names(lst_edge_scores), "_")), byrow = T, ncol = 2)
g1 <- igraph::graph_from_data_frame(df_net, directed = T)
E(g1)$weight <- lst_edge_scores
df_net_weight <- as_data_frame(g1, what = "edges")
g2 <- igraph::graph_from_data_frame(df_net, directed = T)

lst_eigen_centrality1 <- eigen_centrality(g1)$vector
lst_eigen_centrality2 <- eigen_centrality(g2)$vector
plot(lst_eigen_centrality1, lst_eigen_centrality2)


df_out <- ddply(df_net_weight, "from", summarise, n_target = length(weight), weight_sum = sum(weight))
df_out$centrality1 = lst_eigen_centrality1[df_out$from]
df_out$centrality2 = lst_eigen_centrality2[df_out$from]

df_deg <- read.table("../fantom5/data/deg_dermal.fibroblast_iPSC.txt")
df_out_deg <- merge(df_out, df_deg, by.x = "from", by.y = 0)
df_out_deg = transform(df_out_deg, stat_ = stat/sd(stat))
df_out_deg = transform(df_out_deg, score1 = weight_sum)

tfnames <- read.table("../network/data/human_TFs_v2.7_mapped_to_Ensembl98.txt")[, 1]
df_out_subset <- subset(df_out_deg, from %in% tfnames)
df_out_subset <- df_out_subset[order(df_out_subset$score1), ]
tail(df_out_subset, 13)


df_out_subset <- df_out_subset[order(df_out_subset$n_target), ]
tail(df_out_subset, 20)
