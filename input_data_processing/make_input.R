# import functions
setwd("/home/seongwonhwang/Desktop/projects/git/GRN_inference_practice/input_data_processing")
source("make_input_utils.R")

dir_data <- "/home/seongwonhwang/Desktop/projects/mogrify/Statistical\ Consulting/"
MI <- MakeInput$new(
  TEST_ID = "TEST14",
  # path_expr = file.path(dir_data, "Bayesian_DE/iterative_test/iterative_test/TF_experiment_expression_matrix"),
  # path_meta = file.path(dir_data, "Bayesian_DE/iterative_test/iterative_test/TF_experiment_metadata.gz"),
  pseudobulking = F, n_cells_for_selecting = 1, column_name = "label.main",
  path_expr = "data/bulk_anonymize/TF_experiment_expression_matrix_bulk",
  path_meta = "data/bulk_anonymize/TF_experiment_metadata_bulk",
  is_processed = F, use_deg = F, is_singlecell = F,
  cells_to_be_removed = NULL,
  # path_network = file.path(dir_data, "BIC/data/networks_anonymize.txt"),
  path_network = file.path(dir_data, "BIC/data/string_mara_anonymized.txt"),
  # path_network = 'data/jaspar_anonymized.txt',
  # path_network = "data/string_jaspar_overlap.txt",
  path_tf_and_reqdgenes = "/home/seongwonhwang/Desktop/projects/GRN_in_general/PyG/data/Anonymized_tables/Anonymized_tables/",
  path_output = "data/"
)
MI$write_files("training")
MI$write_files("predicting")
