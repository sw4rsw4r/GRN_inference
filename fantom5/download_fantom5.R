# wget https://fantom.gsc.riken.jp/5/datafiles/reprocessed/hg38_v7/extra/CAGE_peaks_df_fantom5ession//hg38_fair+new_CAGE_peaks_phase1and2_tpm_ann.osc.txt.gz -O data/hg38_fair+new_CAGE_peaks_phase1and2_tpm_ann.osc.txt.gz
# gunzip data/hg38_fair+new_CAGE_peaks_phase1and2_tpm_ann.osc.txt.gz

# wget https://fantom.gsc.riken.jp/5/datafiles/reprocessed/hg38_v7/extra/CAGE_peaks_expression/hg38_fair+new_CAGE_peaks_phase1and2_counts_ann.osc.txt.gz -O data/hg38_fair+new_CAGE_peaks_phase1and2_counts_ann.osc.txt.gz
# gunzip data/hg38_fair+new_CAGE_peaks_phase1and2_counts_ann.osc.txt.gz
library(DESeq2)

load_data <- function(fname_path) {
  df_fantom5 <- read.table(fname_path, sep = "\t", check.names = F, row.names = 1, header = T)
  df_fantom5 <- df_fantom5[!grepl("STAT:", rownames(df_fantom5)), ]
  df_fantom5 <- subset(df_fantom5, !is.na(short_description) & grepl("@", short_description) & !grepl("chr[0-9MX]*:|chrY:", short_description))

  info <- strsplit(df_fantom5$short_description, ",")
  df_fantom5 <- df_fantom5[rep(1:nrow(df_fantom5), sapply(info, length)), ]

  geneSymbol <- toupper(gsub("^p[0-9]*@", "", unlist(info)))
  df_expr <- apply(df_fantom5[, -c(1:6)], 2, tapply, geneSymbol, sum)

  names_expr <- sapply(colnames(df_expr), function(x) sub(".hg38.nobarcode$", "", URLdecode(x)))
  ffid <- sapply(strsplit(names_expr, "\\."), function(x) x[length(x)])
  riken_id <- sapply(strsplit(names_expr, "\\."), function(x) tolower(x[length(x) - 1]))
  description <- sapply(strsplit(names_expr, "\\."), function(x) paste(x[2:(length(x) - 2)], collapse = "."))

  colnames(df_expr) <- riken_id
  df_meta <- data.frame(ffid, riken_id, description)
  rownames(df_meta) <- df_meta$riken_id
  return(list(expr = df_expr, meta = df_meta))
}
df_expr_tpm <- load_data("data/hg38_fair+new_CAGE_peaks_phase1and2_tpm_ann.osc.txt")
df_expr_counts <- load_data("data/hg38_fair+new_CAGE_peaks_phase1and2_counts_ann.osc.txt")

write.table(df_expr_counts$expr, "data/hg38_fantom_counts.txt", quote = F, row.names = T, col.names = T, sep = "\t")
write.table(df_expr_tpm$expr, "data/hg38_fantom_tpm.txt", quote = F, row.names = T, col.names = T, sep = "\t")
write.table(df_expr_tpm$meta, "data/hg38_fantom_meta.txt", quote = F, row.names = F, col.names = T, sep = "\t")

df_conv_sample_ids <- read.delim("data/hg38_conv_sample_ids.txt")
df_conv_list <- read.table("data/hg38_conv_list.txt", fill = T)

for (idx in 1:7) {
  group_source <- df_conv_list[idx, 1]
  group_target <- df_conv_list[idx, 2]
  conv_id <- paste0(group_source, "_", group_target)
  samples_selected <- subset(df_conv_sample_ids, group %in% c(group_source, group_target))$sampleID

  df_counts_subset <- df_expr_counts$expr[, samples_selected]
  df_tpm_subset <- df_expr_tpm$expr[, samples_selected]
  df_meta_subset <- df_expr_tpm$meta[samples_selected, ]
  if (group_source == "dermal.fibroblast") {
    df_meta_subset <- transform(df_meta_subset,
      label.main = ifelse(startsWith(description, "Fibroblast"), group_source, group_target)
    )
    # Run DESeq2
    dds <- DESeqDataSetFromMatrix(
      countData = df_counts_subset,
      colData = df_meta_subset,
      design = ~label.main
    )
    dds <- DESeq2::DESeq(dds)
    res <- results(dds)
    res <- res[!is.na(res$padj), ]
    genes_deg <- rownames(subset(res, padj < 0.05))

    genes_expressed_source <- apply(df_tpm_subset[, df_meta_subset$label.main == group_source, drop = F] > 10, 1, mean) == 1
    genes_expressed_target <- apply(df_tpm_subset[, df_meta_subset$label.main == group_target, drop = F] > 10, 1, mean) == 1
    # genes_selected <- apply(df_tpm_subset >= 10, 1, sum) >= 2
    genes_expressed = rownames(df_tpm_subset)[genes_expressed_source | genes_expressed_target]
    genes_selected = intersect(genes_deg, genes_expressed)
    df_filt <- df_tpm_subset[genes_selected, ]
    df_lognorm <- log2(df_filt + .1)
    df_scaled <- t(apply(df_lognorm, 1, scale))
    colnames(df_scaled) <- colnames(df_lognorm)
    write.table(df_lognorm, paste0("data/lognorm_", group_source, "_", group_target, ".txt"), quote = F, row.names = T, col.names = T, sep = "\t")
    write.table(df_scaled, paste0("data/expr_", group_source, "_", group_target, ".txt"), quote = F, row.names = T, col.names = T, sep = "\t")
    write.table(df_meta_subset, paste0("data/meta_", group_source, "_", group_target, ".txt"), quote = F, row.names = F, col.names = T, sep = "\t")
    write.table(res, paste0("data/deg_", group_source, "_", group_target, ".txt"), quote=F, row.names= T, col.names=T, sep='\t')
  }
}
