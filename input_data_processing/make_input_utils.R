library(R6)
library(edgeR)
library(Seurat)
library(scran)
library(DESeq2)
##########################

# Required functions for make_input.R
MakeInput <- R6Class("MakeInput",
    public = list(
        TEST_ID = NULL,
        path_tf_and_reqdgenes = NULL,
        path_output = NULL,
        df_expr = NULL,
        df_net = NULL,
        tfs_all = NULL,
        genes_of_interest = NULL,
        meta = NULL,
        # required_genes = NULL,
        initialize = function(TEST_ID = NULL,
                              path_expr = NULL, path_meta = NULL,
                              pseudobulking = T, n_cells_for_selecting = 150, column_name = "label.main",
                              is_processed = F,
                              use_deg = F, is_singlecell = T,
                              cells_to_be_removed = NULL,
                              path_genes_of_interest = NULL,
                              path_network = NULL, path_tf_and_reqdgenes = NULL,
                              path_output = NULL, rng_seed = 1234) {
            set.seed(rng_seed)
            self$TEST_ID <- TEST_ID
            self$path_tf_and_reqdgenes <- path_tf_and_reqdgenes
            self$path_output <- path_output

            # read data
            self$read_expression(path_expr, path_meta, pseudobulking, n_cells_for_selecting, column_name, is_processed, use_deg, is_singlecell, cells_to_be_removed, path_genes_of_interest)
            self$read_tfs()
            self$read_network(path_network)
        },
        make_pseudobulk = function(data, n_cells_for_selecting) {
            n_rep <- trunc(ncol(data) / n_cells_for_selecting)
            if (n_rep >= 2) {
                selected_ids <- sample(colnames(data), n_rep * n_cells_for_selecting)
                random_group <- sample(rep(1:n_rep, rep(n_cells_for_selecting, n_rep)))
                df_pseudo_counts <- t(apply(data[, selected_ids], 1, tapply, random_group, sum))
            } else {
                selected_ids <- rep(sample(colnames(data)), length.out = n_cells_for_selecting)
                df_pseudo_counts <- matrix(apply(data[, selected_ids], 1, sum), ncol = 1)
            }
            return(df_pseudo_counts)
        },
        get_deg = function(df_expr_full, column_name, cells_to_be_removed, is_singlecell) {
            if (is_singlecell) {
                expt <- CreateSeuratObject(counts = df_expr_full)
                sce <- as.SingleCellExperiment(expt)
                set.seed(1234)
                sce_clust <- scran::quickCluster(sce, min.size = 50)
                # Normalize data
                sce <- scran::computeSumFactors(sce, cluster = sce_clust, min.mean = 0.1)
                sce <- scuttle::logNormCounts(sce)

                filt_cells <- !self$meta[[column_name]] %in% cells_to_be_removed
                sce_filt <- sce[, filt_cells]
                out_wilcox <- scran::findMarkers(sce_filt,
                    assay.type = "logcounts", groups = self$meta[[column_name]][filt_cells], test.type = "wilcox",
                    pval.type = "all"
                )
                out_binom <- scran::findMarkers(sce_filt,
                    assay.type = "logcounts", groups = self$meta[[column_name]][filt_cells], test.type = "binom",
                    pval.type = "all"
                )
                diff_genes <- unique(
                    c(
                        unlist(sapply(names(out_wilcox), function(x) rownames(subset(out_wilcox[[x]], FDR < .05)))),
                        unlist(sapply(names(out_binom), function(x) rownames(subset(out_binom[[x]], FDR < .05))))
                    )
                )
            } else {
                data <- apply(df_expr_full, 2, round, 0)
                dds <- DESeqDataSetFromMatrix(
                    countData = data,
                    colData = self$meta,
                    design = ~label.main
                )
                dds <- DESeq(dds, test = "LRT", reduced = ~1)
                res <- results(dds)
                res <- res[!is.na(res$padj), ]
                diff_genes <- rownames(subset(res, padj < 1e-6))
            }
            return(diff_genes)
        },
        read_expression = function(path_expr, path_meta, pseudobulking, n_cells_for_selecting, column_name, is_processed, use_deg, is_singlecell, cells_to_be_removed, path_genes_of_interest) {
            df_expr_full <- read.delim(path_expr, row.names = 1)
            df_expr_full$id <- NULL
            df_expr_full <- as.matrix(df_expr_full)
            self$meta <- read.delim(path_meta)
            celltypes <- setdiff(self$meta[[column_name]], cells_to_be_removed)

            if (pseudobulking) {
                message("Construct pseudobulk data")
                df_combined_rep <- NULL
                for (celltype in celltypes) {
                    selected_idx <- self$meta[[column_name]] == celltype
                    df_pseudo <- self$make_pseudobulk(df_expr_full[, selected_idx], n_cells_for_selecting)
                    colnames(df_pseudo) <- paste0(celltype, "_rep", 1:ncol(df_pseudo))
                    df_combined_rep <- cbind(df_combined_rep, df_pseudo)
                }
                rownames(df_combined_rep) <- rownames(df_expr_full)
                df_combined <- t(apply(df_combined_rep, 1, tapply, sub("_rep[0-9]*$", "", colnames(df_combined_rep)), median))
            } else if (n_cells_for_selecting == 1) {
                df_combined <- df_expr_full
            } else {
                list_selected <- sapply(celltypes, function(celltype) {
                    message(celltype)
                    selected_idx <- self$meta[[column_name]] == celltype
                    selected_ids <- rep(sample(colnames(df_expr_full)[selected_idx]),
                        length.out = n_cells_for_selecting
                    )
                    df_expr_full[, selected_ids]
                }, simplify = F)
                df_combined <- as.data.frame(list_selected)
            }
            if (!is_processed) {
                if (use_deg) {
                    genes_selected <- self$get_deg(df_expr_full, column_name, cells_to_be_removed, is_singlecell)
                } else {
                    # Lenient filtering
                    genes_selected <- apply(df_combined > 2, 1, sum) > 0
                }
                dge <- DGEList(counts = df_combined[genes_selected, ])
                dge <- calcNormFactors(dge, method = "TMM")
                df_lognorm <- edgeR::cpm(dge, normalized.lib.sizes = T, log = T)
                # Filter genes
                # df_norm <- df_lognorm[apply(df_lognorm, 1, sd) > 0.8, ]
                df_scaled <- t(apply(df_lognorm, 1, scale))
                colnames(df_scaled) <- colnames(df_combined)
            } else {
                df_scaled <- df_combined
            }
            self$df_expr <- df_scaled
            if (!is.null(path_genes_of_interest) && file.exists(path_genes_of_interest)) {
                self$genes_of_interest <- read.table(path_genes_of_interest)[, 1]
            } else {
                self$genes_of_interest <- rownames(self$df_expr)
            }
        },
        read_tfs = function() {
            if (file.exists(self$path_tf_and_reqdgenes) && !file.info(self$path_tf_and_reqdgenes)$isdir) {
                tfs <- setdiff(read.delim(self$path_tf_and_reqdgenes)$Symbol, c("-", "", "Symbol"))
                self$tfs_all <- tfs
            } else {
                input_tfs <- file.path(self$path_tf_and_reqdgenes, "tfs_anonymized.txt")
                tfs <- read.table(input_tfs)[, 1]
                external_tfs <- grep("tf", rownames(self$df_expr), value = T)
                self$tfs_all <- unique(c(external_tfs, tfs))
            }
        },
        read_network = function(path_network) {
            df_net <- read.delim(path_network, header = F)
            df_net <- df_net[!grepl("target", df_net$V2), ]
            colnames(df_net) <- c("TF", "target")
            df_net_subset <- df_net[, c("TF", "target")]
            df_net_subset <- df_net_subset[!duplicated(df_net_subset), ]
            self$df_net <- df_net_subset
        },
        write_files = function(type, use_all_possible_edges_for_prediction = F, file_deg = NULL) {
            message("N genes of interest: ", length(self$genes_of_interest))
            genes_expressed <- intersect(rownames(self$df_expr), self$genes_of_interest)
            # df_net_filt <- subset(self$df_net, TF %in% intersect(self$tfs_all, genes_expressed) & target %in% genes_expressed)
            df_net_filt <- subset(self$df_net, TF %in% genes_expressed & target %in% genes_expressed)
            selected_genes <- unique(c(df_net_filt$TF, df_net_filt$target))

            df_embeddings <- self$df_expr[selected_genes, ]
            if (type == "predicting") {
                # df_net_string = read.delim(file.path("/home/seongwonhwang/Desktop/projects/mogrify/Statistical\ Consulting/", "BIC/data/networks_anonymize.txt"))
                # df_net_string <- df_net_string[!duplicated(df_net_string), ]
                # colnames(df_net_string) <- c("TF", "target")
                # df_net_string <- subset(df_net_string, TF %in% intersect(self$tfs_all, rownames(df_embeddings)) & target %in% rownames(df_embeddings))

                if (!is.null(file_deg)) {
                    df_deg <- read.delim(file_deg)
                    df_deg <- df_deg[order(df_deg$padj), ]
                    genes_sig <- head(rownames(subset(df_deg, padj < .05)), 1000)
                }
                if (use_all_possible_edges_for_prediction & is.null(file_deg)) {
                    ###############################################
                    # write all possible edges from TFs to target
                    selected_tfs <- unique(df_net_filt$TF)
                    index_all_tfs <- match(selected_tfs, rownames(df_embeddings)) - 1
                    index_all_genes <- seq_len(nrow(df_embeddings)) - 1
                    df_all_possible_edges_idx <- expand.grid(index_all_tfs, index_all_genes)
                    df_all_possible_edges <- expand.grid(selected_tfs, selected_genes)
                } else if (use_all_possible_edges_for_prediction & !is.null(file_deg)) {
                    selected_tfs <- intersect(df_net_filt$TF, genes_sig)
                    selected_diff_genes <- intersect(df_net_filt$target, genes_sig)
                    index_all_tfs <- match(selected_tfs, rownames(df_embeddings)) - 1
                    index_all_genes <- match(selected_diff_genes, rownames(df_embeddings)) - 1
                    df_all_possible_edges_idx <- expand.grid(index_all_tfs, index_all_genes)
                    df_all_possible_edges <- expand.grid(selected_tfs, selected_diff_genes)
                } else {
                    df_all_possible_edges_idx <- data.frame(
                        TF = match(df_net_filt$TF, rownames(df_embeddings)) - 1,
                        target = match(df_net_filt$target, rownames(df_embeddings)) - 1
                    )
                    df_all_possible_edges <- df_net_filt
                }

                write.table(df_all_possible_edges_idx,
                    paste0(
                        self$path_output, "/",
                        self$TEST_ID, "_edge_index_for_predicting.txt"
                    ),
                    quote = F, row.names = F, col.names = F, sep = "\t"
                )

                write.table(df_all_possible_edges,
                    paste0(
                        self$path_output, "/",
                        self$TEST_ID, "_all_possible_edges.txt"
                    ),
                    quote = F, row.names = F, col.names = F, sep = "\t"
                )
                ###############################################
            } else if (type == "training") {
                df_edge_index <- data.frame(
                    TF = match(df_net_filt$TF, rownames(df_embeddings)) - 1,
                    target = match(df_net_filt$target, rownames(df_embeddings)) - 1
                )
                write.table(df_edge_index,
                    paste0(self$path_output, "/", self$TEST_ID, "_edge_index_for_training.txt"),
                    quote = F, row.names = F, col.names = F, sep = "\t"
                )
            }
            write.table(df_embeddings,
                paste0(self$path_output, "/", self$TEST_ID, "_embeddings_for_", type, ".txt"),
                quote = F, row.names = F, col.names = F, sep = "\t"
            )
            write.table(data.frame(rownames(df_embeddings)),
                paste0(self$path_output, "/", self$TEST_ID, "_embeddings_gene_for_", type, ".txt"),
                quote = F, row.names = F, col.names = F, sep = "\t"
            )
        }
    )
)
