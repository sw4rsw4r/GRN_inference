# import functions
setwd("/home/seongwonhwang/Desktop/projects/git/GRN_inference_practice/input_data_processing")
source("make_input_utils.R")

MI <- MakeInput$new(
  TEST_ID = "TEST15v5",
  pseudobulking = F, n_cells_for_selecting = 1, column_name = "label.main",
  path_expr = "../fantom5/data/expr_dermal.fibroblast_iPSC.txt",
  path_meta = "../fantom5/data/meta_dermal.fibroblast_iPSC.txt",
  is_processed = T, use_deg = F, is_singlecell = F,
  cells_to_be_removed = NULL,
  path_network = "../network/data/human_string_v2.7_withdirection_new_noexpfilter",
  path_tf_and_reqdgenes = "../network/data/human_TFs_v2.7_mapped_to_Ensembl98.txt",
  path_output = "data/"
)
MI$write_files("training")
MI$write_files("predicting")