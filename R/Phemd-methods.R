##################
# Class constructor
#################

#' Phemd class
#'
#' The main PhEMD class to store single-cell expression data.
#' @field data List of matrices, each of which represents a single-cell sample (num_cells x num_genes)
#' @field markers Column names (e.g. genes) for each element (i.e. data matrix) in "data"
#' @field snames Sample ID for each element in "data"
#' @field data_aggregate Numeric matrix representing expression data for cells from all experimental conditions (rows = markers, cols = cells)
#' @field data_subsample_idx List of vectors each representing the indices of elements in "data" that were subsampled and combined to form "data_aggregate"
#' @field subsampled_bool Boolean represent whether or not subsampling was performed in the data aggregation process
#' @field monocle_obj Data object of type "CellDataSet" that is the core Monocle data structure
#' @field data_cluster_weights Matrix representing cell subtype relative frequencies for each sample (num_samples x num_genes)
#' @field emd_dist_mat Matrix representing pairwise distances between each pair of cell subtypes
#' @field seurat_obj Object of type "seurat" that is the core Seurat data structure
#' @field experiment_ids Vector of length num_samples representing the experiment (batch) in which the sample was profiled
#' @name Phemd
#' @rdname Phemd
#' @aliases Phemd-class
#' @exportClass Phemd
#' @importClassesFrom Seurat Seurat

setClassUnion("CDSorNULL",members=c('CellDataSet', "NULL"))
setClassUnion("SeuratorNULL",members=c('Seurat', "NULL"))
setClass("Phemd",
         contains=c('CellDataSet', 'Seurat'),
         slots=c(data = "list",
                 markers = "character",
                 snames = "character",
                 data_aggregate = "matrix",
                 data_subsample_idx = "list",
                 subsampled_bool = "logical",
                 monocle_obj = "CDSorNULL",
                 data_cluster_weights = "matrix",
                 emd_dist_mat = "matrix",
                 seurat_obj = "SeuratorNULL",
                 experiment_ids = "character", 
                 version='package_version'))

###########################
# Methods for Phemd class
##########################
#' @name Phemd-methods
#' @docType methods
#' @rdname Phemd-methods
#' 
setValidity("Phemd", function(object) {
  if(length(rawExpn(object)) < 1) {
    return('Phemd object must have at least 1 sample in rawExpn(object)')
  }
  if(length(sNames(object)) != length(rawExpn(object))) {
    return('sNames(object) must be the same length as rawExpn(object)')
  }
  if(sum(dim(pooledCells(object))) == 0 && ncol(rawExpn(object)[[1]]) != length(selectMarkers(object))) {
    return('Number of markers measured in rawExpn(object) must equal number of markers listed in selectMarkers(object)')
  }
  if(sum(dim(pooledCells(object))) != 0 && nrow(pooledCells(object)) != length(selectMarkers(object))) {
    return('Number of markers measured in pooledCells(object) must equal number of markers listed in selectMarkers(object)')
  }
  return(TRUE)  
})


##################
# Accessor functions
###################

#' Accessor function for stored multi-sample raw expression data
#' 
#' @param obj A Phemd object.
#' @return List of matrices, each of which represents a single-cell sample
#' @export
#' @examples
#' phemdObj <- createDataObj(all_expn_data, all_genes, as.character(snames_data))
#' raw_expn_data <- rawExpn(phemdObj)
#' 
rawExpn <- function(obj) {
  stopifnot(is(obj,"Phemd"))
  obj@data
}


#' Accessor function for stored Monocle object
#' 
#' @param obj A Phemd object.
#' @return An object of class 'CellDataSet' (from Monocle)
#' @export
#' @examples
#' phemdObj <- createDataObj(all_expn_data, all_genes, as.character(snames_data))
#' monocle_obj <- monocleInfo(phemdObj)
#' 
monocleInfo <- function(obj) {
  stopifnot(is(obj,"Phemd"))
  obj@monocle_obj
}

#' Accessor function for stored Seurat object within Phemd object
#' 
#' @param obj A Phemd object.
#' @return An object of class 'Seurat'
#' @export
#' @examples
#' phemdObj <- createDataObj(all_expn_data, all_genes, as.character(snames_data))
#' seurat_obj <- seuratInfo(phemdObj)
#' 
seuratInfo <- function(obj) {
  stopifnot(is(obj,"Phemd"))
  obj@seurat_obj
}

#' Accessor function for EMD ground distance matrix
#' 
#' @param obj A Phemd object
#' @return Sqaure matrix representing pairwise distances between cell subtypes
#' @export
#' @examples
#' phemdObj <- createDataObj(all_expn_data, all_genes, as.character(snames_data))
#' gdm <- GDM(phemdObj)
#' 
GDM <- function(obj) {
  stopifnot(is(obj,"Phemd"))
  obj@emd_dist_mat
}

#' Accessor function for gene/protein markers measured in experiment
#' 
#' @param obj Phemd object
#' @return Vector representing gene/protein markers corresponding to expression matrices
#' @export
#' @examples
#' phemdObj <- createDataObj(all_expn_data, all_genes, as.character(snames_data))
#' genes <- selectMarkers(phemdObj)
#' 
selectMarkers <- function(obj) {
  stopifnot(is(obj,"Phemd"))
  obj@markers
}

#' Accessor function for identifiers of all single-cell samples in experiment
#' 
#' @param obj Phemd object
#' @return Vector representing sample names corresponding to expression matrices
#' @export
#' @examples
#' phemdObj <- createDataObj(all_expn_data, all_genes, as.character(snames_data))
#' sampleIDs <- sNames(phemdObj)
#' 
sNames <- function(obj) {
  stopifnot(is(obj,"Phemd"))
  obj@snames
}


#' Accessor function for aggregated cells used for cell subtype definition
#' 
#' @param obj Phemd object
#' @return Numeric matrix representing expression data for cells from all experimental conditions (rows = markers, cols = cells)
#' @export
#' @examples
#' phemdObj <- createDataObj(all_expn_data, all_genes, as.character(snames_data))
#' cells_aggregated <- pooledCells(phemdObj)
#' 
pooledCells <- function(obj) {
  stopifnot(is(obj,"Phemd"))
  obj@data_aggregate
}

#' Accessor function for aggregated cells used for cell subtype definition
#' 
#' @param obj Phemd object
#' @return List of vectors each representing the indices of elements in rawExpn(obj) that were subsampled and combined to form "data_aggregate"
#' @export
#' @examples
#' phemdObj <- createDataObj(all_expn_data, all_genes, as.character(snames_data))
#' subsampled_idx_list <- subsampledIdx(phemdObj)
#' 
subsampledIdx <- function(obj) {
  stopifnot(is(obj,"Phemd"))
  obj@data_subsample_idx
}

#' Accessor function for whether or not cells were subsampled when aggregated for cell subtype analysis
#' 
#' @param obj Phemd object
#' @return Boolean represent whether or not subsampling was performed in the data aggregation process
#' @export
#' @examples
#' phemdObj <- createDataObj(all_expn_data, all_genes, as.character(snames_data))
#' subsampled <- subsampledBool(phemdObj)
#' 
subsampledBool <- function(obj) {
  stopifnot(is(obj,"Phemd"))
  obj@subsampled_bool
}

#' Accessor function for cell subtype distribution for each sample
#' 
#' @param obj Phemd object
#' @return Matrix representing cell subtype relative frequencies for each sample (num_samples x num_genes)
#' @export
#' @examples
#' phemdObj <- createDataObj(all_expn_data, all_genes, as.character(snames_data))
#' celltype_weights <- celltypeFreqs(phemdObj)
#' 
celltypeFreqs <- function(obj) {
  stopifnot(is(obj,"Phemd"))
  obj@data_cluster_weights
}

#' Accessor function for batch ID for each sample
#' 
#' @param obj Phemd object
#' @return Vector of length num_samples representing the experiment (batch) in which the sample was profiled
#' @export
#' @examples
#' phemdObj <- createDataObj(all_expn_data, all_genes, as.character(snames_data))
#' batch_metadata <- batchIDs(phemdObj)
#' 
batchIDs <- function(obj) {
  stopifnot(is(obj,"Phemd"))
  obj@experiment_ids
}

##################
# Setter functions
###################


#' Setter function for protein / gene markers
#' 
#' @rdname Phemd-methods
#' @docType methods
#' @return Updated Phemd object
#' @export
#' @examples
#' phemdObj <- createDataObj(all_expn_data, all_genes, as.character(snames_data))
#' new_genes <- all_genes
#' new_genes[1] <- 'IL2R'
#' selectMarkers(phemdObj) <- new_genes
#' 
setGeneric("selectMarkers<-", function(obj, value) standardGeneric("selectMarkers<-"))

#' @rdname Phemd-methods
#' @aliases Phemd,ANY,ANY-method
setMethod("selectMarkers<-", "Phemd", function(obj, value) {
  obj@markers <- value
  validObject(obj)
  obj
})

#' Setter function for stored expression data
#' 
#' @rdname Phemd-methods
#' @aliases Phemd,character,ANY-method
#' @export
#' @examples
#' phemdObj <- createDataObj(all_expn_data, all_genes, as.character(snames_data))
#' new_expn_data <- all_expn_data
#' new_expn_data <- lapply(new_expn_data, function(x) {log2(x+1)})
#' rawExpn(phemdObj) <- new_expn_data
#'
setGeneric("rawExpn<-", function(obj, value) standardGeneric("rawExpn<-"))

#' @rdname Phemd-methods
#' @aliases Phemd,ANY,ANY-method
setMethod("rawExpn<-", "Phemd", function(obj, value) {
  obj@data <- value
  validObject(obj)
  obj
})

#' Setter function for single-cell expression data aggregated from multiple samples
#' 
#' @rdname Phemd-methods
#' @docType methods
#' @return Updated Phemd object
#' @export
#' @examples
#' phemdObj <- createDataObj(all_expn_data, all_genes, as.character(snames_data))
#' aggregated_data <- t(do.call(rbind,all_expn_data))
#' pooledCells(phemdObj) <- aggregated_data
#' 
setGeneric("pooledCells<-", function(obj, value) standardGeneric("pooledCells<-"))

#' @rdname Phemd-methods
#' @aliases Phemd,ANY,ANY-method
setMethod("pooledCells<-", "Phemd", function(obj, value) {
  obj@data_aggregate <- value
  validObject(obj)
  obj
})

#' Setter function for indices of cells subsampled from each sample during aggregation
#' 
#' @rdname Phemd-methods
#' @docType methods
#' @return Updated Phemd object
#' @export
#' @examples
#' phemdObj <- createDataObj(all_expn_data, all_genes, as.character(snames_data))
#' subsampledIdxList<- rep(list(1:10), length(all_expn_data)) #subsampled cells 1-10 from each sample
#' subsampledIdx(phemdObj) <- subsampledIdxList
#' 
setGeneric("subsampledIdx<-", function(obj, value) standardGeneric("subsampledIdx<-"))

#' @rdname Phemd-methods
#' @aliases Phemd,ANY,ANY-method
setMethod("subsampledIdx<-", "Phemd", function(obj, value) {
  obj@data_subsample_idx <- value
  validObject(obj)
  obj
})

#' Setter function for boolean denoting whether cells were subsampled from each sample during aggregation
#' 
#' @rdname Phemd-methods
#' @docType methods
#' @return Updated Phemd object
#' @export
#' @examples
#' phemdObj <- createDataObj(all_expn_data, all_genes, as.character(snames_data))
#' subsampledBool(phemdObj) <- TRUE
#' 
setGeneric("subsampledBool<-", function(obj, value) standardGeneric("subsampledBool<-"))

#' @rdname Phemd-methods
#' @aliases Phemd,ANY,ANY-method
setMethod("subsampledBool<-", "Phemd", function(obj, value) {
  obj@subsampled_bool <- value
  validObject(obj)
  obj
})

#' Setter function for Monocle2 CellDataSet object for experiment
#' 
#' @rdname Phemd-methods
#' @param obj A Phemd object
#' @param value Assignment object
#' @return Updated Phemd object
#' @export
#' @examples
#' phemdObj <- createDataObj(all_expn_data, all_genes, as.character(snames_data))
#' mydata <- pooledCells(phemdObj)
#' myCellDataSet <- newCellDataSet(mydata,phenoData=NULL, expressionFamily=VGAM::negbinomial.size())
#' monocleInfo(phemdObj) <- myCellDataSet
#' 
setGeneric("monocleInfo<-", function(obj, value) standardGeneric("monocleInfo<-"))

#' @rdname Phemd-methods
#' @aliases Phemd,ANY,ANY-method
setMethod("monocleInfo<-", "Phemd", function(obj, value) {
  obj@monocle_obj <- value
  validObject(obj)
  obj
})

#' Setter function for Seurat object for experiment
#' 
#' @rdname Phemd-methods
#' @docType methods
#' @return Updated Phemd object containing Seurat object
#' @export
#' @examples
#' phemdObj <- createDataObj(all_expn_data, all_genes, as.character(snames_data))
#' my_seuratObj <- Seurat::CreateSeuratObject(raw.data = t(all_expn_data[[1]]), project = "A")
#' seuratInfo(phemdObj) <- my_seuratObj
#' 
setGeneric("seuratInfo<-", function(obj, value) standardGeneric("seuratInfo<-"))

#' @rdname Phemd-methods
#' @aliases Phemd,ANY,ANY-method
setMethod("seuratInfo<-", "Phemd", function(obj, value) {
  obj@seurat_obj <- value
  validObject(obj)
  obj
})

#' Setter function for cell subtype frequencies of each single-cell sample
#' 
#' @rdname Phemd-methods
#' @docType methods
#' @return Updated Phemd object
#' @export
#' @examples
#' phemdObj <- createDataObj(all_expn_data, all_genes, as.character(snames_data))
#' myCellTypeFreqs <- matrix(rexp(length(all_expn_data)*10, rate=.1), ncol=10)
#' myCellTypeFreqs <- apply(myCellTypeFreqs, 1, function(x) {x / sum(x)})
#' celltypeFreqs(phemdObj) <- myCellTypeFreqs
#' 
setGeneric("celltypeFreqs<-", function(obj, value) standardGeneric("celltypeFreqs<-"))

#' @rdname Phemd-methods
#' @aliases Phemd,ANY,ANY-method
setMethod("celltypeFreqs<-", "Phemd", function(obj, value) {
  obj@data_cluster_weights <- value
  validObject(obj)
  obj
})

#' Setter function for batch IDs of each single-cell sample
#' 
#' @rdname Phemd-methods
#' @docType methods
#' @return Updated Phemd object
#' @export
#' @examples
#' phemdObj <- createDataObj(all_expn_data, all_genes, as.character(snames_data))
#' my_seuratObj <- Seurat::CreateSeuratObject(raw.data = t(all_expn_data[[1]]), project = "A")
#' seuratInfo(phemdObj) <- my_seuratObj
#' batchIDs(phemdObj) <- rep('A', length(all_expn_data))
#' 
setGeneric("batchIDs<-", function(obj, value) standardGeneric("batchIDs<-"))

#' @rdname Phemd-methods
#' @aliases Phemd,ANY,ANY-method
setMethod("batchIDs<-", "Phemd", function(obj, value) {
  obj@experiment_ids <- value
  validObject(obj)
  obj
})

#' Setter function for EMD ground distance matrix
#' 
#' @rdname Phemd-methods
#' @docType methods
#' @return Updated Phemd object
#' @export
#' @examples
#' phemdObj <- createDataObj(all_expn_data, all_genes, as.character(snames_data))
#' cluster_locs <- 1:10
#' myGDM <- as.matrix(dist(cluster_locs))
#' GDM(phemdObj) <- myGDM
#' 
setGeneric("GDM<-", function(obj, value) standardGeneric("GDM<-"))

#' @rdname Phemd-methods
#' @aliases Phemd,ANY,ANY-method
setMethod("GDM<-", "Phemd", function(obj, value) {
  obj@emd_dist_mat <- value
  validObject(obj)
  obj
})

