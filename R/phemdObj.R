##################
# Class constructor
#################

#' phemdObj class
#'
#' The main PhEMD class to store single-cell expression data.
#' @field data List of matrices, each of which represents a single-cell sample (num_cells x num_genes)
#' @field markers Column names (e.g. genes) for each element (i.e. data matrix) in "data"
#' @field snames Sample ID for each element in "data"
#' @field data_aggregate Numeric matrix representing (rows = markers, cols = cells)
#' @field data_subsample_idx List of vectors each representing the indices of elements in "data" that were subsampled and combined to form "data_aggregate"
#' @field subsampled_bool Boolean represent whether or not subsampling was performed in the data aggregation process
#' @field monocle_obj Data object of type "CellDataSet" that is the core Monocle data structure
#' @field data_cluster_weights Matrix representing cell subtype relative frequencies for each sample (num_samples x num_genes)
#' @field emd_dist_mat Matrix representing pairwise distances between each pair of single-cell samples
#' @field seurat_obj Object of type "seurat" that is the core Seurat data structure
#' @field experiment_ids Vector of length num_samples representing the experiment
#' @name phemdObj
#' @rdname phemdObj
#' @aliases phemdObj-class
#' @exportClass phemdObj
#' @import monocle Seurat

setClassUnion("CDSorNULL",members=c("CellDataSet", "NULL"))
#setClass('CellDataSet')
setClass("phemdObj",
         contains=c('CellDataSet', 'seurat'),
         slots=c(data = "list",
                 markers = "character",
                 snames = "character",
                 data_aggregate = "matrix",
                 data_subsample_idx = "list",
                 subsampled_bool = "logical",
                 monocle_obj = "CDSorNULL",
                 data_cluster_weights = "matrix",
                 emd_dist_mat = "matrix",
                 seurat_obj = "seurat",
                 experiment_ids = "character"))
