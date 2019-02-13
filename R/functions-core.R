################
## FUNCTIONS ###
################

#########################################################
### Private methods below (not exported in namespace) ###
#########################################################

#' @title Retrieve single-cell sample sizes
#' @description Takes initial list of single-cell samples and returns vector containing number of cells in each sample.
#' @details Private method (not exported in namespace)
#' @param data_list List of length num_samples (each element has dimension num_cells x num_markers)
#' @return Vector of length num_samples representing number of cells in each sample
#' @examples
#' \dontrun{
#' sample_sizes <- getSampleSizes(all_expn_data)
#' }

getSampleSizes <- function(data_list) {
  return(vapply(data_list, nrow, integer(1L)))
}

#' @title Retrieve reference cell clusters
#' @description Takes initial Phemd struct and returns cell clusters as assigned by clustering algorithm (i.e. Monocle 2)
#' @details Private method (not exported in namespace)
#' @param obj Phemd struct containing Monocle2 object and underlying expression data
#' @param cell_model String representing data model for cell state space (Seurat or Monocle 2)
#' @param expn_type String representing whether to return raw expression values or coordinates in dimensionality-reduced, aligned feature space (only relevant for Seurat data models)
#' @param ndim Number of dimensions (e.g. CCA) to use (only relevant for Seurat data models)
#' @return List of data matrices; each list element is of size num_cells_in_cluster x num_markers and represents a distinct cell cluster
#' @examples
#' \dontrun{
#' cluster_expression_data <- retrieveRefClusters(my_phemdObj)
#' }
#' 

retrieveRefClusters <- function(obj, cell_model=c('monocle2','seurat'), 
                                expn_type='aligned', ndim=10) {
  cell_model <- match.arg(cell_model, c('monocle2','seurat'))
  if(cell_model == 'monocle2') {
    # Extract state labels from monocle data object
    monocle_obj <- monocleInfo(obj)
    labels <- pData(phenoData(monocle_obj))
    state_labels <- as.numeric(labels$State)

    # Split data frame based on cluster assignments
    mydata <- as.data.frame(t(pooledCells(obj)))
    ref_clusters <- split(mydata, state_labels)
    
  } else if(cell_model == 'seurat') {
    seurat_obj <- seuratInfo(obj)
    state_labels <- as.numeric(as.character(GetIdent(seurat_obj, uniq=FALSE)))
    if(min(state_labels) == 0) state_labels <- state_labels + 1 #ensure cluster labels are 1 indexed instead of zero indexed
    names(state_labels) <- names(GetIdent(seurat_obj, uniq=FALSE)) # label cluster assignments with cell name
    if(expn_type == 'aligned') {
      # aligned CCA expression data (num_cells x num_markers)
      mydata <- GetDimReduction(object = seurat_obj, reduction.type = 'cca.aligned',
                                           slot = "cell.embeddings")[,seq_len(ndim)]
    } else if(expn_type == 'raw') {
      mydata <- t(as.matrix(GetAssayData(seurat_obj, assay.type='RNA', slot='raw.data')))
    } else {
      stop('Error: expn_type must be either "raw" or "aligned"')
    }
    
    # Split data frame based on cluster assignments
    mydata <- as.data.frame(mydata)
    ref_clusters <- split(mydata, state_labels)
  } else {
      stop('Error: cell_model must be either "monocle2" or "seurat"')
  }
  return(ref_clusters)
}


#' @title Identify cluster centroids (cell names)
#' @description Takes initial list and returns list of cell names representing centroid of cluster
#' @details Private method (not exported in namespace)
#' @param ref_clusters list containing each cluster of interest (each list element is a matrix of dimension num_cells x num_markers)
#' @return List of names; element \var{i} represents the name of the cell in cluster \var{i} that is closest to the centroid (arithmetic mean) of cluster \var{i}
#' @examples
#' \dontrun{
#' centroid_names <- identifyCentroids(ref_clusters)
#' }

identifyCentroids <- function(ref_clusters) {
  centroids <- lapply(ref_clusters, function(cur_cluster) {
    arith_centroid <- colMeans(cur_cluster)
    curdist <- t(as.matrix(apply(cur_cluster, 1, function(x) norm(x-arith_centroid, type="2"))))
    closest_cell <- colnames(curdist)[which.min(curdist)] #closest cell to arithmetic centroid
    return(closest_cell)
  })
  return(centroids)
}


#' @title Get arithmetic centroids (coordinates)
#' @description Takes initial list and returns a matrix with row \var{i} representing the arithmetic centroid of cluster \var{i}
#' @details Private method (not exported in namespace)
#' @param ref_clusters list containing each cluster of interest (each list element is a matrix of dimension num_cells x num_markers)
#' @return Matrix of dimension num_cluster x num_markers; row \var{i} representing the arithmetic centroid of cluster \var{i}
#' @examples
#' \dontrun{
#' cluster_centroids <- getArithmeticCentroids(ref_clusters)
#' }
getArithmeticCentroids <- function(ref_clusters) {
  if(length(ref_clusters) < 1) stop('Error: input requires at least 1 reference cluster')
  
  #centroids <- matrix(0, nrow=length(ref_clusters), ncol = ncol(ref_clusters[[1]]))
  #for(i in seq_len(length(ref_clusters))) {
  #  cur_cluster <- ref_clusters[[i]]
  #  centroids[i,] <- colMeans(cur_cluster)
  #}
  
  centroids_list <- lapply(ref_clusters, colMeans)
  centroids <- do.call(rbind, centroids_list)
  return(centroids)
}

#' @title Assign cells to a reference cell subtype
#' @description Assigns each cell in \code{cur_cells} to a cluster based on nearest cell in Monocle 2 tree
#' @details Private method (not exported in namespace). Uses RANN package for fast knn search
#' @param cur_cells Matrix of cells to be assigned to clusters (Dim: \var{num_cells} x \var{num_markers})
#' @param ref_cells Matrix of cells used to build reference Monocle 2 tree (Dim: \var{num_monocle_cells} x \var{num_markers})
#' @param ref_cell_labels Vector of length \var{num_monocle_cells} containing Monocle 2 cell branch assignments
#' @param cell_model Either "monocle2" or "seurat" depending on method used to model cell state space
#' @return Vector of length \var{num_cells} representing cluster assignments for each cell in \var{cur_cells}
#' @examples
#' \dontrun{
#' cur_cells_cluster_labels <- assignCellClusterNearestNode(cur_cells_expn_data, 
#' clustered_cells_expn_data, clustered_cells_cluster_labels, cell_model='monocle2')
#' }
assignCellClusterNearestNode <- function(cur_cells, ref_cells, ref_cell_labels, cell_model=c('monocle2', 'seurat')) {
  if(nrow(ref_cells) != length(ref_cell_labels)) stop("Error: number of cells and cell labels do not match")

  closest <- RANN::nn2(data = ref_cells, query = cur_cells, k = 1) #fast nearest neighbor search
  nearest_cell <- closest$nn.idx

  cell_model <- match.arg(cell_model, c('monocle2','seurat'))
  if(cell_model == 'monocle2') {
    assigned <- ref_cell_labels[nearest_cell]
  } else if(cell_model == 'seurat') {
    nearest_cell_names <- rownames(ref_cells)[nearest_cell]
    assigned <- as.numeric(ref_cell_labels[nearest_cell_names])
  } else {
    stop('Error: cell_model must be either monocle2 or seurat')
  }

  return(assigned)
}

#' @title Models expression data using generalized linear model with Gaussian error
#' @description Useful for modeling pre-normalized single-cell expression data.
#' @details Private method (not to be called by user directly). Requires VGAM package. Obtained from VGAM v1.0-5 (https://www.rdocumentation.org/packages/VGAM/versions/1.0-5/topics/gaussianff)
#' @param dispersion Dispersion parameter. If 0, then estimate as described in VGAM 1.0-5 documentation.
#' @param parallel A logical or formula. If a formula, the response of the formula should be a logical and the terms of the formula indicates whether or not those terms are parallel.
#' @param zero An integer-valued vector specifying which linear/additive predictors are modelled as intercepts only. The values must be from the set {1...M} where Mis the number of columns of the matrix response.
#' @return Generalized linear model with Gaussian error
#' 
gaussianffLocal <- function(dispersion = 0, parallel = FALSE, zero = NULL) {
  if (!VGAM::is.Numeric(dispersion, length.arg = 1) || dispersion <
    0) {
    stop("bad input for argument 'dispersion'")
  }
  estimated.dispersion <- dispersion == 0
  cur_constraints <- expression({
    cur_constraints <- VGAM::cm.VGAM(matrix(1, M, 1),
                                 x = x, bool = parallel,
                                 constraints = constraints
    )
    cur_constraints <- VGAM::cm.zero.VGAM(constraints,
                                      x = x, zero,
                                      M = M, predictors.names = predictors.names, M1 = 1
    )
  })
  deviance_fn <- function(mu,y, w, residuals = FALSE, eta, extra = NULL) {
    M <- if (is.matrix(y)) {
      ncol(y)
    } else {
      1
    }
    n <- if (is.matrix(y)) {
      nrow(y)
    } else {
      length(y)
    }
    wz <- VGAM:::VGAM.weights.function(w = w, M = M, n = n)
    if (residuals) {
      if (M > 1) {
        U <- vchol(wz, M = M, n = n)
        temp <- mux22(U, y - mu,
                      M = M, upper = TRUE,
                      as.matrix = TRUE
        )
        dimnames(temp) <- dimnames(y)
        temp
      }
      else {
        (y - mu) * sqrt(wz)
      }
    }
    else {
      ResSS.vgam(y - mu, wz = wz, M = M)
    }
  }
  
  infos_fn <- function(...) {
    list(
      M1 = 1, Q1 = 1, expected = TRUE, multipleResponses = TRUE,
      quasi.type = TRUE, zero = zero
    )
  }
  
  init_expn <- expression({
    if (is.R()) assign("CQO.FastAlgorithm", TRUE, envir = VGAM::VGAMenv) else CQO.FastAlgorithm <<- TRUE
    if (any(function.name == c("cqo", "cao")) && (length(zero) ||
                                                  (is.logical(parallel) && parallel))) {
      stop("cannot handle non-default arguments for cqo() and cao()")
    }
    temp5 <- w.y.check(
      w = w, y = y, ncol.w.max = Inf, ncol.y.max = Inf,
      out.wy = TRUE, colsyperw = 1, maximize = TRUE
    )
    w <- temp5$w
    y <- temp5$y
    M <- if (is.matrix(y)) ncol(y) else 1
    dy <- dimnames(y)
    predictors.names <- if (!is.null(dy[[2]])) {
      dy[[2]]
    } else {
      param.names(
        "Y",
        M
      )
    }
    if (!length(etastart)) etastart <- 0 * y
  })
  
  last_expn <- expression({
    dy <- dimnames(y)
    if (!is.null(dy[[2]])) dimnames(fit$fitted.values) <- dy
    dpar <- dispersion
    if (!dpar) {
      wz <- VGAM:::VGAM.weights.function(w = w, M = M, n = n)
      temp5 <- ResSS.vgam(y - mu, wz = wz, M = M)
      dpar <- temp5 / (length(y) - (if (is.numeric(ncol(X.vlm.save))) ncol(X.vlm.save) else 0))
    }
    misc$dispersion <- dpar
    misc$default.dispersion <- 0
    misc$estimated.dispersion <- estimated.dispersion
    misc$link <- rep_len("identitylink", M)
    names(misc$link) <- predictors.names
    misc$earg <- vector("list", M)
    for (ilocal in seq_len(M)) misc$earg[[ilocal]] <- list()
    names(misc$link) <- predictors.names
    if (is.R()) {
      if (exists("CQO.FastAlgorithm", envir = VGAM::VGAMenv)) {
        rm("CQO.FastAlgorithm",
           envir = VGAM::VGAMenv
        )
      }
    } else {
      while (exists("CQO.FastAlgorithm")) remove("CQO.FastAlgorithm")
    }
    misc$expected <- TRUE
    misc$multipleResponses <- TRUE
  })
  
  loglikelihood_fn <- function(mu, y, w, residuals = FALSE,
                               eta, extra = NULL, summation = TRUE) {
    M <- if (is.matrix(y)) {
      ncol(y)
    } else {
      1
    }
    n <- if (is.matrix(y)) {
      nrow(y)
    } else {
      length(y)
    }
    wz <- VGAM:::VGAM.weights.function(w = w, M = M, n = n)
    if (residuals) {
      stop("loglikelihood residuals not implemented yet")
    }
    else {
      temp1 <- ResSS.vgam(y - mu, wz = wz, M = M)
      ll.elts <- if (M == 1 || ncol(wz) == M) {
        -0.5 * temp1 + 0.5 * (log(wz)) - n * (M / 2) *
          log(2 * pi)
      }
      else {
        if (all(wz[1, ] == apply(wz, 2, min)) && all(wz[1, ] == apply(wz, 2, max))) {
          onewz <- m2a(wz[1, , drop = FALSE], M = M)
          onewz <- onewz[, , 1]
          logdet <- determinant(onewz)$modulus
          logretval <- -0.5 * temp1 + 0.5 * n * logdet -
            n * (M / 2) * log(2 * pi)
          distval <- stop("variable 'distval' not computed yet")
          logretval <- -(ncol(onewz) * log(2 * pi) +
                           logdet + distval) / 2
          logretval
        }
        else {
          logretval <- -0.5 * temp1 - n * (M / 2) * log(2 * pi)
          for (ii in seq_len(n)) {
            onewz <- m2a(wz[ii, , drop = FALSE], M = M)
            onewz <- onewz[, , 1]
            logdet <- determinant(onewz)$modulus
            logretval <- logretval + 0.5 * logdet
          }
          logretval
        }
      }
      if (summation) {
        sum(ll.elts)
      }
      else {
        ll.elts
      }
    }
  }
  validparams_fn <- function(eta, y, extra = NULL) {
    okay1 <- all(is.finite(eta))
    okay1
  }
  
  new("vglmff",
    blurb = c(
      "Vector linear/additive model\n",
      "Links:    identitylink for Y1,...,YM"
    ), constraints = cur_constraints, deviance = deviance_fn, infos = infos_fn, 
    initialize = init_expn, linkinv = function(eta, extra = NULL) eta, 
    last = last_expn,
    loglikelihood = loglikelihood_fn, linkfun = function(mu, extra = NULL) mu, 
    vfamily = "gaussianff",
    validparams = validparams_fn, deriv = expression({
      wz <- VGAM:::VGAM.weights.function(w = w, M = M, n = n)
      mux22(cc = t(wz), xmat = y - mu, M = M, as.matrix = TRUE)
    }), weight = expression({
      wz
    })
  )
}

############################
### Public methods below ###
############################

#' @title Create 'Phemd' object
#' @description Wrapper function to create 'Phemd' object containing raw expression data and metadata
#' @details Note that each element in list can have different number of rows (i.e. number of cells in each sample can vary).
#' @param data List of length \var{num_samples} containing expression data; each element is of size \var{num_cells} x \var{num_markers}. Alternately a SingleCellExperiment object.
#' @param markers Vector containing marker names (i.e. column names of \code{all_data})
#' @param snames Vector containing sample names (i.e. names of samples contained in \code{all_data})
#' @param datatype Either "list" or "sce" (SingleCellExperiment with genes x cells)
#' @param valtype Type of assay data (i.e. "counts", "normcounts", "logcounts", "tpm", "cpm") if datatype is "sce"
#' @return 'Phemd' object containing raw multi-sample expression data and associated metadata
#' @examples
#' my_phemdObj <- createDataObj(all_expn_data, all_genes, as.character(snames_data))
#' 
createDataObj <- function(data, markers, snames, datatype='list', valtype='counts') {
  if(datatype == 'list') {
    all_data <- data
  } else if(datatype == 'sce') {
    if(valtype %in% names(assays(data))) {
      stop(sprintf('Error: %s not found in input SingleCellExperiment', valtype))
    }
    all_data <- t(assay(data, valtype)) # generally SCE objects are genes x cells; transpose
  } else {
    stop('Error: Input datatype must be either "list" or "sce" (SingleCellExperiment)')
  }

  nsample <- length(all_data)
  if(nsample == 0) stop('all_data is empty (length=0)')
  stopifnot(nsample == length(snames))
  if(nsample > 0) {
    stopifnot(ncol(all_data[[1]]) == length(markers))
  }

  # Error-checking to ensure that all samples in data list have same number of markers
  nmarker_vec <- rep(0, nsample)
  for(i in seq_len(nsample)) {
    nmarker_vec[i] <- ncol(all_data[[i]])
  }
  if(sum(nmarker_vec - nmarker_vec[1]) != 0) {
    stop(sprintf("Error: Sample %d has a different number of columns than Sample 1", which(nmarker_vec-nmarker_vec[1] != 0)))
  }
  data_obj <- new('Phemd', data = all_data, markers = markers, snames = snames, monocle_obj=NULL)

  return(data_obj)
}


#' @title Attach 'seurat' object to 'Phemd' object
#' @description Allows user to attach batch-normalized reference cell data from Seurat into 'Phemd' object containing raw expression data and metadata
#' @param phemd_obj Phemd object initialized using createDataObj
#' @param seurat_obj S4 'seurat' object containing batch-normalized reference cell data
#' @param batch.colname Name of column in Seurat object that denotes batch ID
#' @return 'Phemd' object containing with attached Seurat object
#' @examples
#' my_phemdObj <- createDataObj(all_expn_data, all_genes, as.character(snames_data))
#' my_seuratObj <- Seurat::CreateSeuratObject(raw.data = t(all_expn_data[[1]]), project = "A")
#' my_seuratObj <- Seurat::ScaleData(object = my_seuratObj, do.scale=FALSE, do.center=FALSE)
#' my_seuratObj <- Seurat::RunPCA(object = my_seuratObj, pc.genes = colnames(all_expn_data[[1]]), do.print = FALSE)
#' my_seuratObj <- Seurat::FindClusters(my_seuratObj, reduction.type = "pca", dims.use = 1:10, resolution = 0.6, print.output = 0, save.SNN = TRUE)
#' my_phemdObj <- bindSeuratObj(my_phemdObj, my_seuratObj)
#' 
bindSeuratObj <- function(phemd_obj, seurat_obj, batch.colname='plt') {
  stopifnot(is(seurat_obj,'seurat'))
  # ensure cluster names are 1-indexed
  if(min(as.numeric(as.character(GetIdent(seurat_obj, uniq=FALSE)))) == 0) {
    label_names <- names(GetIdent(seurat_obj, uniq=FALSE))
    labels_renumbered <- factor(as.numeric(as.character(GetIdent(seurat_obj, uniq=FALSE))) +1)
    names(labels_renumbered) <- label_names
    seurat_obj <- SetIdent(object=seurat_obj, ident.use=labels_renumbered)
  }
  if(batch.colname != 'plt') {
    seurat_obj@meta.data$plt <- seurat_obj@meta.data[[batch.colname]]
  }
  seuratInfo(phemd_obj) <- seurat_obj

  return(phemd_obj)
}

#' @title Remove samples with too few cells
#' @description Removes samples from Phemd that have fewer cells than \code{min_sz}
#' @details Note: If used, this function must be called before (and not after) the \code{aggregateSamples} function is called
#' @param obj 'Phemd' object containing raw expression data and associated metadata
#' @param min_sz Minimum number of cells in each sample to be retained
#' @return 'Phemd' object containing raw multi-sample expression data and associated metadata (same as input minus removed samples)
#' @examples
#' my_phemdObj <- createDataObj(all_expn_data, all_genes, as.character(snames_data))
#' my_phemdObj_lg <- removeTinySamples(my_phemdObj, 10) #removes samples with fewer than 10 cells
#' 
removeTinySamples <- function(obj, min_sz=20) {
  stopifnot(is(obj,'Phemd'))
  stopifnot(mode(min_sz) == 'numeric')
  all_data <- rawExpn(obj)
  all_snames <- sNames(obj)
  all_sample_sz <- getSampleSizes(all_data)
  to_remove_idx <- which(all_sample_sz < min_sz)
  if(length(to_remove_idx) == 0) return(obj)

  to_remove_idx <- to_remove_idx[order(to_remove_idx, decreasing=TRUE)] #remove from end to front
  for(i in to_remove_idx) {
    print(sprintf('%s removed because only contains %d cells', all_snames[i], all_sample_sz[i]))
    all_data[[i]] <- NULL
  }
  all_snames <- all_snames[-to_remove_idx]

  obj@data <- all_data
  obj@snames <- all_snames
  validObject(obj)
  return(obj)
}

#' @title Aggregate expression data from all samples
#' @description Takes initial Phemd object and returns object with additional data frame in slot @@data_aggregate containing cells aggregated from all samples (to be used for further analyses e.g. Monocle 2 trajectory building / pseudotime mapping / cell clustering)
#' @details Subsamples cells as necessary based on \code{max_cells}. If subsampling is performed, an equal number of cells are subsampled from each sample
#' @param obj 'Phemd' object containing raw expression data and associated metadata
#' @param max_cells Maximum number of cells across all samples to be included in final matrix on which Monocle 2 will be run
#' @return Same as input 'Phemd' object with additional slot 'data_aggregate' containing aggregated expression data (num_markers x num_cells)
#' @examples
#' my_phemdObj <- createDataObj(all_expn_data, all_genes, as.character(snames_data))
#' my_phemdObj_lg <- removeTinySamples(my_phemdObj, 10)
#' my_phemdObj_lg <- aggregateSamples(my_phemdObj_lg, max_cells=1000, cur_seed=112)
#' 
aggregateSamples <- function(obj, max_cells=12000) {
  stopifnot(is(obj, 'Phemd'))
  stopifnot(mode(max_cells) == 'numeric')
  
  all_data <- rawExpn(obj)
  nsample <- length(all_data)
  if(nsample == 0) return(obj)

  subsample_sz <- floor(max_cells/nsample)
  all_aggregate_data <- matrix(nrow=0,ncol=ncol(all_data[[1]]))
  all_subsample_idx <- list()
  subsample_bool = FALSE

  all_sample_sz <- getSampleSizes(all_data)
  if(sum(all_sample_sz) > max_cells) subsample_bool = TRUE

  for(i in seq_len(nsample)) {
    cur_data <- all_data[[i]]
    # take all cells unless total cells across all samples > max_cells
    if(subsample_bool) {
      cur_subsample_idx <- sample(seq_len(nrow(cur_data)), min(subsample_sz, nrow(cur_data)))
      all_subsample_idx[[i]] <- cur_subsample_idx
      all_aggregate_data <- rbind(all_aggregate_data, cur_data[cur_subsample_idx,])
    } else {
      all_aggregate_data <- rbind(all_aggregate_data, cur_data)
    }
  }
  all_aggregate_data <- t(all_aggregate_data) #rows = markers, cols = cells
  colnames(all_aggregate_data) <- seq_len(ncol(all_aggregate_data))
  rownames(all_aggregate_data) <- selectMarkers(obj)
  
  pooledCells(obj) <- as.matrix(all_aggregate_data)
  subsampledIdx(obj) <- all_subsample_idx
  subsampledBool(obj) <- subsample_bool
  return(obj)
}

#' @title Perform feature selection on aggregated data
#' @description Takes as input a Phemd object with aggregated data and returns updated object after performing feature selection on aggregated data
#' @details \code{aggregateSamples} needs to be called before running this function
#' @param obj 'Phemd' object containing aggregated data
#' @param selected_genes Vector containing names of genes to use for downstream analyses
#' @return Same as input 'Phemd' object after performing feature-selection based dimensionality reduction on aggregated expression data
#' @examples
#' 
#' my_phemdObj <- createDataObj(all_expn_data, all_genes, as.character(snames_data))
#' my_phemdObj_lg <- removeTinySamples(my_phemdObj, 10)
#' my_phemdObj_lg <- aggregateSamples(my_phemdObj_lg, max_cells=1000)
#' my_phemdObj_lg <- selectFeatures(my_phemdObj_lg, selected_genes=c('TP53', 
#' 'EGFR', 'KRAS', 'FOXP3', 'LAG3'))
#' 
selectFeatures <- function(obj, selected_genes) {
  if(isempty(pooledCells(obj))) stop('slot "data_aggregate" is empty; please call aggregateSamples() before running selectFeatures()')
  all_aggregate_data <- pooledCells(obj)
  all_genes <- selectMarkers(obj)
  selected_gene_map <- match(selected_genes, all_genes)
  #TODO: print genes in selected_genes that were unable to map to all_genes
  if(sum(!is.na(selected_gene_map)) == 0) {
    stop('None of the genes in "selected_genes" were found. Aborting feature selection.')
  }
  
  all_aggregate_data <- all_aggregate_data[selected_gene_map,]
  obj@data_aggregate <- all_aggregate_data #set slots manually due to validity constraints
  obj@markers <- all_genes[selected_gene_map]
  validObject(obj)
  return(obj)
}

#' @title Generate Monocle2 embedding
#' @description Takes as input a Phemd object with aggregated data and returns updated object with Monocle2 object in @@monocle_obj slot
#' @details Wrapper function for \code{reduceDimension} in Monocle 2 package. \code{aggregateSamples} needs to be called before running this function.
#' @param obj 'Phemd' object containing aggregated data
#' @param data_model One of the following: 'negbinomial_sz', 'negbinomial', 'tobit', 'uninormal', 'gaussianff'. See "Family Function" table at the following link for more details on selecting the proper one. \url{http://cole-trapnell-lab.github.io/monocle-release/docs/#getting-started-with-monocle}
#' @param ... Additional parameters to be passed to \code{reduceDimension} function
#' @return Same as input 'Phemd' object with additional Monocle2 object in @@monocle_obj slot
#' @examples
#' my_phemdObj <- createDataObj(all_expn_data, all_genes, as.character(snames_data))
#' my_phemdObj_lg <- removeTinySamples(my_phemdObj, 10)
#' my_phemdObj_lg <- aggregateSamples(my_phemdObj_lg, max_cells=1000)
#' my_phemdObj_lg <- embedCells(my_phemdObj_lg, data_model = 'gaussianff', sigma=0.02, maxIter=2)
embedCells <- function(obj, data_model = 'negbinomial_sz', ...) {
  extra_args <- list(...)

  if(isempty(pooledCells(obj))) stop('slot "data_aggregate" is empty; please call aggregateSamples() before running embedCells()')
  mydata <- pooledCells(obj)
  if(is.null(mydata)) stop("Error: call 'aggregateSamples' function first ")

  myFeatureData <- as.data.frame(selectMarkers(obj))
  fd <- new("AnnotatedDataFrame", data=myFeatureData)
  rownames(fd) <- selectMarkers(obj)

  if(is.null(data_model)) {
    print('Assuming data fit negative binomial pattern of expression...')
    data_model <- 'negbinomial_sz'
  }

  if(data_model == 'negbinomial_sz') {
    expression_fam_fn <- VGAM::negbinomial.size()
  } else if(data_model == 'negbinomial') {
    expression_fam_fn <- VGAM::negbinomial()
  } else if(data_model == 'tobit') {
    expression_fam_fn <- VGAM::tobit()
  } else if(data_model == 'uninormal') {
    expression_fam_fn <- VGAM::uninormal()
  } else if(data_model == 'gaussianff') {
    #expression_fam_fn <- VGAM::gaussianff()
    expression_fam_fn <- gaussianffLocal()
  } else {
    stop("Error: Invalid data_model specified")
  }

  # Helpful to use ncenter and sigma that's higher than default in reduceDimension for datasets w/ relatively more cells
  if(!('ncenter' %in% names(extra_args)) && ncol(mydata) > 3000) {
    extra_args['ncenter'] <- 750
    if(!('sigma' %in% names(extra_args))) extra_args['sigma'] <- 0.03
  }

  if(!('maxIter' %in% names(extra_args))) {
    extra_args['maxIter'] <- 12 #set maximum number of iterations to 12
  }

  monocle_obj <- newCellDataSet(mydata,phenoData=NULL,featureData=fd, 
                                expressionFamily=expression_fam_fn)
  varLabels(featureData(monocle_obj)) <- 'gene_short_name' #random formatting requirement for monocle

  monocle_obj <- estimateSizeFactors(monocle_obj)
  if(data_model == 'negbinomial_sz') {
    monocle_obj <- estimateDispersions(monocle_obj)
  }
  if(data_model == 'gaussianff') {
    extra_args['norm_method'] <- 'none'
    extra_args['scaling'] <- 'FALSE'
  }

  rd_args <- c(list(cds=monocle_obj, max_components=2, reduction_method='DDRTree'),
               extra_args[names(extra_args) %in% c("verbose", "ncenter", "norm_method", 'scaling', 'pseudo_expr', "initial_method", "maxIter", "sigma", "lambda", "param.gamma", "tol")])
  monocle_obj_red <- do.call(reduceDimension, rd_args)

  monocleInfo(obj) <- monocle_obj_red
  return(obj)
}

#' @title Compute Monocle2 cell state and pseudotime assignments
#' @description Takes as input a Phemd object with Monocle2 object and returns updated object with Monocle2 object containing cell state and pseudotime assignments
#' @details Wrapper function for \code{orderCells} in Monocle 2 package. \code{embedCells} needs to be called before calling this function.
#' @param obj 'Phemd' object containing initial Monocle 2 object
#' @param ... Additional parameters to be passed into \code{orderCells} function
#' @return Same as input 'Phemd' object with updated Monocle2 object in @@monocle_obj slot containing cell state and pseudotime assignments
#' @examples
#' 
#' my_phemdObj <- createDataObj(all_expn_data, all_genes, as.character(snames_data))
#' my_phemdObj_lg <- removeTinySamples(my_phemdObj, 10)
#' my_phemdObj_lg <- aggregateSamples(my_phemdObj_lg, max_cells=1000)
#' my_phemdObj_monocle <- embedCells(my_phemdObj_lg, data_model='gaussianff', sigma=0.02, maxIter=2)
#' my_phemdObj_monocle <- orderCellsMonocle(my_phemdObj_monocle)
orderCellsMonocle <- function(obj, ...) {
  monocle_obj <- monocleInfo(obj)
  if(ncol(monocle_obj) == 0) stop('slot "monocle_obj" is empty; please call embedCells() before running orderCellsMonocle()')

  extra_args <- list(...)
  oc_args <- c(list(cds=monocle_obj),
               extra_args[names(extra_args) %in% c("root_state", "reverse")])

  monocle_obj_ordered <- do.call(orderCells, oc_args)
  monocleInfo(obj) <- monocle_obj_ordered
  return(obj)
}

#' @title Computes cell subtype abundances for each sample
#' @description Takes as input a Phemd object with all single-cell expression data of all single-cell samples in @@data slot and Monocle2 object (already embedded and ordered) in @@monocle_obj slot. Returns updated object with cell subtype frequencies of each sample in @@data_cluster_weights slot
#' @details \code{embedCells} and \code{orderCellsMonocle} need to be called before calling this function.
#' @param obj 'Phemd' object containing single-cell expression data of all samples in @@data slot and Monocle2 object (already embedded and ordered) in @@monocle_obj slot
#' @param verbose Boolean that determines whether progress (sequential processing of samples) should be printed. FALSE by default
#' @param cell_model Either "monocle2" or "seurat" depending on method used to model cell state space
#' @return 'Phemd' object with cell subtype frequencies of each sample in @@data_cluster_weights slot
#' @examples
#' 
#' my_phemdObj <- createDataObj(all_expn_data, all_genes, as.character(snames_data))
#' my_phemdObj_lg <- removeTinySamples(my_phemdObj, 10)
#' my_phemdObj_lg <- aggregateSamples(my_phemdObj_lg, max_cells=1000)
#' my_phemdObj_monocle <- embedCells(my_phemdObj_lg, data_model = 'gaussianff', sigma=0.02, maxIter=2)
#' my_phemdObj_monocle <- orderCellsMonocle(my_phemdObj_monocle)
#' my_phemdObj_final <- clusterIndividualSamples(my_phemdObj_monocle)
#' 
clusterIndividualSamples <- function(obj, verbose=FALSE, cell_model=c('monocle2', 'seurat')) {
  stopifnot(is(obj,'Phemd'))
  all_data <- rawExpn(obj)
  cell_model <- match.arg(cell_model, c('monocle2','seurat'))
  if(cell_model == 'monocle2') {
    monocle_obj <- monocleInfo(obj)
    if(ncol(monocle_obj) == 0) stop('slot "monocle_obj" is empty; please call embedCells() and orderCellsMonocle() before calling this function')

    # Extract state labels from monocle data object
    labels <- pData(phenoData(monocle_obj))
    if(!('State' %in% names(labels))) stop('monocleInfo(obj) does not have cell state assignments; please call embedCells() and orderCellsMonocle() before calling this function')
    state_labels <- as.numeric(labels$State)
    # retrieve reference clusters
    ref_clusters <- retrieveRefClusters(obj, cell_model='monocle2')
    nclusters <- length(ref_clusters)

    if(!subsampledBool(obj)) {
      ## subsampling not performed; all cells assigned to clusters as-is

      cluster_weights <- matrix(0, nrow=length(all_data), ncol = nclusters)

      start_idx <- 1
      for(i in seq_len(length(all_data))) {
        cur_sample_sz <- nrow(all_data[[i]])
        end_idx <- start_idx + cur_sample_sz - 1
        sample_labels <- state_labels[start_idx:end_idx] #cell subtype assignments for current sample
        cur_hist <- rep(0, nclusters)
        for(j in seq_len(nclusters)) {
          cur_hist[j] <- sum(sample_labels == j)
        }
        cur_hist <- cur_hist / sum(cur_hist)
        cluster_weights[i,] <- cur_hist
        start_idx <- end_idx + 1
      }
      celltypeFreqs(obj) <- cluster_weights
    } else {
      ## subsampling performed; need to assign cells to cluster of nearest cell in embedding
      refcluster_sizes <- rep(0, length(ref_clusters))
      counter1 <- 0
      for(i in seq_len(nclusters)) {
        refcluster_sizes[i] <- nrow(ref_clusters[[i]])
        counter1 <- counter1 + nrow(ref_clusters[[i]])
      }
      print(refcluster_sizes / counter1)

      cluster_weights <- matrix(0, nrow=length(all_data), ncol=nclusters)
      for(i in seq_len(length(all_data))) {
        cur_data <- all_data[[i]]
        cur_ncells <- nrow(cur_data)
        cur_hist <- rep(0, nclusters)

        if(verbose && i %% 10 == 0) {
          print(sprintf('Processing sample number: %d', i))
        }
        
        # Use nearest-cell mapping instead of nearest-centroid mapping
        cur_cell_labels <- assignCellClusterNearestNode(cur_data, 
                                                        t(exprs(monocle_obj)), 
                                                        state_labels, 
                                                        cell_model=cell_model) 

        for(j in seq_len(nclusters)) {
          cur_hist[j] <- sum(cur_cell_labels == j)
        }
        cur_hist <- cur_hist / sum(cur_hist)
        cluster_weights[i,] <- cur_hist
      }
      print(colMeans(cluster_weights))

      celltypeFreqs(obj) <- cluster_weights
    }
  } else if(cell_model == 'seurat') {
    seurat_obj <- seuratInfo(obj)
    # retrieve reference clusters for starting estimate of cluster sizes
    ref_clusters <- retrieveRefClusters(obj, cell_model='seurat', 
                                        expn_type = 'raw')
    nclusters <- length(ref_clusters)
    ## subsampling performed; assign cells to same cluster as nearest cell 
    # in embedding
    refcluster_sizes <- rep(0, length(ref_clusters))
    counter1 <- 0
    for(i in seq_len(nclusters)) {
      refcluster_sizes[i] <- nrow(ref_clusters[[i]])
      counter1 <- counter1 + nrow(ref_clusters[[i]])
    }
    print(refcluster_sizes / counter1)

    cluster_weights <- matrix(0, nrow=length(all_data), ncol = nclusters)
    sample_batches <- batchIDs(obj)
    for(i in seq_len(length(all_data))) {
      cur_data <- all_data[[i]]
      cur_plt <- sample_batches[i]
      if(verbose && i %% 10 == 0) {
        print(sprintf('Processing sample number: %d', i))
      }

      cur_hist <- rep(0, nclusters)
      state_labels <- as.numeric(as.character(GetIdent(seurat_obj, uniq=FALSE)))
      names(state_labels) <- rownames(seurat_obj@meta.data)
      ref_data <- t(as.matrix(GetAssayData(seurat_obj, assay.type='RNA', slot='raw.data')))
      cell_idx_curplt <- which(seurat_obj@meta.data$plt == cur_plt)
      if(length(cell_idx_curplt) == 0) {
        stop(sprintf('Error: no cells in reference set match the experiment_id %s of sample %d', cur_plt, i))
      }
      ref_data <- ref_data[cell_idx_curplt,]
      state_labels <- state_labels[cell_idx_curplt]

      # Use nearest-cell mapping
      cur_cell_labels <- assignCellClusterNearestNode(cur_data, ref_data, 
                                                      state_labels, 
                                                      cell_model=cell_model) 

      for(j in seq_len(nclusters)) {
        cur_hist[j] <- sum(cur_cell_labels == j)
      }
      cur_hist <- cur_hist / sum(cur_hist);
      cluster_weights[i,] <- cur_hist
    }
    print(colMeans(cluster_weights))

    celltypeFreqs(obj) <- cluster_weights
  } else {
    stop('Error: cell_model must be either "monocle2" or "seurat"')
  }

  return(obj)
}


#' @title Computes ground distance matrix based on cell embedding
#' @description Takes as input a Phemd object with Monocle2 object (already embedded and ordered) in @@monocle_obj slot. Returns updated object with ground distance matrix representing pairwise distance between 2 cell subtypes based on cell state embedding.
#' @details \code{embedCells} and \code{orderCellsMonocle} need to be called before calling this function. Requires 'igraph' package
#' @param obj 'Phemd' object containing Monocle2 object (already embedded and ordered) in @@monocle_obj slot
#' @param cell_model Method by which cell state was modeled (either "monocle2" or "seurat")
#' @return Phemd object with ground distance matrix (to be used in EMD computation) in @@data_cluster_weights slot
#' @examples
#' 
#' my_phemdObj <- createDataObj(all_expn_data, all_genes, as.character(snames_data))
#' my_phemdObj_lg <- removeTinySamples(my_phemdObj, 10)
#' my_phemdObj_lg <- aggregateSamples(my_phemdObj_lg, max_cells=1000)
#' my_phemdObj_monocle <- embedCells(my_phemdObj_lg, data_model = 'gaussianff', sigma=0.02, maxIter=2)
#' my_phemdObj_monocle <- orderCellsMonocle(my_phemdObj_monocle)
#' my_phemdObj_final <- clusterIndividualSamples(my_phemdObj_monocle)
#' my_phemdObj_final <- generateGDM(my_phemdObj_final)
#' 
generateGDM <- function(obj, cell_model=c('monocle2', 'seurat')) {
  stopifnot(is(obj,'Phemd'))

  cell_model <- match.arg(cell_model, c('monocle2','seurat'))
  if(cell_model == 'monocle2') {
    monocle_obj <- monocleInfo(obj)
    # retrieve reference clusters
    ref_clusters <- retrieveRefClusters(obj, cell_model='monocle2')
    nclusters <- length(ref_clusters)

    # get graph underlying Monocle tree
    mst_graph <- monocle_obj@auxOrderingData[["DDRTree"]]$pr_graph_cell_proj_tree 
    centroids <- identifyCentroids(ref_clusters)
    pseudotimes <- pData(phenoData(monocle_obj))$Pseudotime
    emd_dists <- matrix(0, nrow=nclusters, ncol=nclusters)
    for(i in seq_len(nclusters)) {
      for(j in seq_len(nclusters)) {
        if(i == j) next
        path_btwn_cells <- shortest_paths(mst_graph, centroids[[i]], centroids[[j]])$vpath[[1]] #cell1 and cell2 should be strings representing vertices in mst
        pseudotime_curpath <- pseudotimes[path_btwn_cells]
        pseudotime_dist <- abs(pseudotimes[as.numeric(centroids[[i]])] - min(pseudotime_curpath)) + abs(pseudotimes[as.numeric(centroids[[j]])] - min(pseudotime_curpath))
        emd_dists[i,j] <- pseudotime_dist
      }
    }
    GDM(obj) <- emd_dists
  } else if(cell_model == 'seurat') {
    ref_clusters <- retrieveRefClusters(obj, cell_model='seurat', expn_type='aligned', ndim=8)
    seurat_obj <- seuratInfo(obj)
    centroids <- getArithmeticCentroids(ref_clusters)
    GDM(obj) <- as.matrix(dist(centroids))
  } else {
    stop('Error: cell_model must be either "monocle2" or "seurat"')
  }
  return(obj)
}

#' @title Computes EMD distance matrix representing pairwise dissimilarity between samples
#' @description Takes as input a Phemd object with cell subtype relative frequencies for each sample in @@data_cluster_weights slot and ground distance matrix (representing cell subtype pairwise dissimilarity) in @@emd_dist_mat slot. Returns distance matrix representing pairwise dissimilarity between samples
#' @details Requires 'transport' and 'pracma' packages
#' @param obj 'Phemd' object containing cell subtype relative frequencies for each sample in @@data_cluster_weights slot and ground distance matrix (representing cell subtype dissimilarity) in @@emd_dist_mat slot
#' @return Distance matrix of dimension num_samples x num_samples representing pairwise dissimilarity between samples
#' @examples
#'
#' my_phemdObj <- createDataObj(all_expn_data, all_genes, as.character(snames_data))
#' my_phemdObj_lg <- removeTinySamples(my_phemdObj, 10)
#' my_phemdObj_lg <- aggregateSamples(my_phemdObj_lg, max_cells=1000)
#' my_phemdObj_monocle <- embedCells(my_phemdObj_lg, data_model = 'gaussianff', sigma=0.02, maxIter=2)
#' my_phemdObj_monocle <- orderCellsMonocle(my_phemdObj_monocle)
#' my_phemdObj_final <- clusterIndividualSamples(my_phemdObj_monocle)
#' my_phemdObj_final <- generateGDM(my_phemdObj_final)
#' my_EMD_mat <- compareSamples(my_phemdObj_final)
#' 
compareSamples <- function(obj) {
  stopifnot(is(obj,'Phemd'))
  cluster_weights <- celltypeFreqs(obj)
  emd_dists <- GDM(obj)
  # generate inhibitor distance matrix
  Y <- rep(0, (nrow(cluster_weights)-1)*nrow(cluster_weights)/2)
  counter <- 1
  for(i in seq_len(nrow(cluster_weights))){
    cur_f1_weights <- cluster_weights[i,]
    for(j in (i+1):nrow(cluster_weights)) {
      if(j > nrow(cluster_weights)) break #this doesn't automatically happen in R
      cur_f2_weights <- cluster_weights[j,]
      # Compute EMD
      flow <- transport(cur_f1_weights, cur_f2_weights, emd_dists, 
                       method='primaldual')
      curdist <- 0
      for(k in seq_len(nrow(flow))) {
        cur_penalty <- emd_dists[flow[k,1], flow[k,2]]
        curdist <- curdist+cur_penalty*flow[k,3]
      }
      Y[counter] <- curdist
      counter <- counter + 1
    }
  }
  Y_sq <- squareform(Y)
  return(Y_sq)
}

#' @title Performs community detection on sample-sample distance matrix
#' @description Takes sample-sample distance matrix as input and returns group assignments for each sample
#' @details By default, uses 'kgs' (Kelley-Gardner-Sutcliffe) method for determining optimal number of groups. Alternatively, can take user-specified number of groups). Requires 'cluster' and 'maptree' packages.
#' @param distmat Distance matrix of dimension num_samples x num_samples representing pairwise dissimilarity between samples
#' @param distfun Method of partitioning network of samples (currently either 'hclust' or 'pam')
#' @param ncluster Optional parameter specifying total number of sample groups
#' @param method Optional parameter for hierarchical clustering (see "hclust" documentation)
#' @return Vector containing group assignments for each sample (same order as row-order of distmat) based on user-specified partitioning method (e.g. hierarchical clustering)
#' @examples
#' 
#' my_phemdObj <- createDataObj(all_expn_data, all_genes, as.character(snames_data))
#' my_phemdObj_lg <- removeTinySamples(my_phemdObj, 10)
#' my_phemdObj_lg <- aggregateSamples(my_phemdObj_lg, max_cells=1000)
#' my_phemdObj_monocle <- embedCells(my_phemdObj_lg, data_model = 'gaussianff', sigma=0.02, maxIter=2)
#' my_phemdObj_monocle <- orderCellsMonocle(my_phemdObj_monocle)
#' my_phemdObj_final <- clusterIndividualSamples(my_phemdObj_monocle)
#' my_phemdObj_final <- generateGDM(my_phemdObj_final)
#' my_EMD_mat <- compareSamples(my_phemdObj_final)
#' cluster_assignments <- groupSamples(my_EMD_mat, distfun = 'hclust', ncluster=4)
#' 
groupSamples <- function(distmat, distfun = 'hclust', ncluster=NULL, method='complete') {
  ## Identify similar groups of inhibitors (hierarchical clustering for now)
  if(nrow(distmat) != ncol(distmat)) {
    stop('Error: distmat must be a square distance matrix of dimension num_samples x num_samples')
  }
  if(distfun == 'hclust') {
    cluster_results <- hclust(as.dist(distmat), method=method)
    if(is.null(ncluster)) {
      # kgs method for determining optimal number of clusters
      op_k <- kgs(cluster_results, as.dist(distmat), maxclust = 15)
      ncluster <- op_k[which(op_k == min(op_k))]
    }
    cluster_assignments <- cutree(cluster_results, k=ncluster)
  } else if(distfun == 'pam') {
    if(is.null(ncluster)) ncluster <- 4
    cluster_results <- pam(distmat, ncluster, diss=TRUE)
    cluster_assignments=cluster_results$clustering
  } else {
    stop("Error: Please specify distfun as either 'hclust' or 'pam'")
  }
  return(cluster_assignments)
}





