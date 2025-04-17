options(stringsAsFactors = F)
library(data.table)
library(plyr)

f_action <- fread("data/9606.protein.actions.v11.0.txt")
f_action_filt <- subset(f_action, is_directional == "t")
f_action_c <- ddply(f_action_filt, c("item_id_a", "item_id_b"), summarise,
  mode = paste(unique(sort(mode)), collapse = ","),
  action = paste(unique(sort(action)), collapse = ","),
  is_directional = paste(unique(sort(is_directional)), collapse = ","),
  a_is_acting = paste(unique(sort(a_is_acting)), collapse = ","),
  score = paste(unique(sort(score)), collapse = ",")
)
f_action_w <- subset(f_action_c, a_is_acting == "t")
write.table(f_action_w, "data/directional_string_edges.txt", quote = F, row.names = F, col.names = T, sep = "\t")


f_detail <- fread("data/9606.protein.links.detailed.v11.0.txt")
# f_detail_filt1 = subset(f_detail, experimental> 0 & combined_score> 500 )
f_detail_filt1 <- subset(f_detail, combined_score > 500)

f_detail_filt2 <- subset(f_detail_filt1, !paste(protein1, protein2) %in% with(f_action_w, paste(item_id_b, item_id_a)))
f_detail_filt2$is_direction <- with(f_detail_filt2, paste(protein1, protein2)) %in% with(f_action_w, paste(item_id_a, item_id_b))

f_detail_filt3 <- with(f_detail_filt2, data.frame(protein1, protein2, combined_score, is_direction))
write.table(f_detail_filt3, "data/v11.0_string_new.raw", quote = F, row.names = F, col.names = T, sep = " ")
