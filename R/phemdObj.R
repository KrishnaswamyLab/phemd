##################
# Class constructor
#################
setClass("phemdObj",
         slots=c(data = "list",
                 markers = "character",
                 snames = "character",
                 data_aggregate = "matrix",
                 data_subsample_idx = "list",
                 subsampled_bool = "logical",
                 monocle_obj = "CellDataSet",
                 data_cluster_weights = "matrix",
                 emd_dist_mat = "matrix",
                 seurat_obj = "seurat",
                 experiment_ids = "character"))
