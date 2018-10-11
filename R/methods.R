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
#' sample_sizes <- get_sample_sizes(all_expn_data)
#' }

get_sample_sizes <- function(data_list) {
  sample_sizes <- rep(0, length(data_list))
  for(i in 1:length(sample_sizes)) {
    sample_sizes[i] <- nrow(data_list[[i]])
  }
  return(sample_sizes)
}

#' @title Retrieve reference cell clusters
#' @description Takes initial phemdObj struct and returns cell clusters as assigned by clustering algorithm (i.e. Monocle 2)
#' @details Private method (not exported in namespace)
#' @param obj phemdObj struct containing Monocle2 object and underlying expression data
#' @param cell_model String representing data model for cell state space (Seurat or Monocle 2)
#' @param expn_type String representing whether to return raw expression values or coordinates in dimensionality-reduced, aligned feature space (only relevant for Seurat data models)
#' @param ndim Number of dimensions (e.g. CCA) to use (only relevant for Seurat data models)
#' @return List of data matrices; each list element is of size num_cells_in_cluster x num_markers and represents a distinct cell cluster
#' @examples
#' \dontrun{
#' cluster_expression_data <- retrieve_refclusters(my_phemdObj)
#' }
#' 

retrieve_refclusters <- function(obj, cell_model='monocle2', expn_type='aligned', ndim=10) {
  ref_clusters = list();
  if(cell_model == 'monocle2') {
    # Extract state labels from monocle data object
    monocle_obj <- obj@monocle_obj
    labels <- pData(phenoData(monocle_obj))
    state_labels <- as.numeric(labels$State)

    mydata <- t(obj@data_aggregate)
    for(i in 1:max(state_labels)) {
      ref_clusters[[i]] <- mydata[state_labels == i,]
    }
  } else if(cell_model == 'seurat') {
    seurat_obj <- obj@seurat_obj
    state_labels <- as.numeric(as.character(seurat_obj@ident))
    if(min(state_labels) == 0) state_labels <- state_labels + 1 #ensure cluster labels are 1 indexed instead of zero indexed
    names(state_labels) <- names(seurat_obj@ident) # label cluster assignments with cell name
    if(expn_type == 'aligned') {
      # aligned CCA expression data (ncells x nmarkers)
      mydata <- GetDimReduction(object = seurat_obj, reduction.type = 'cca.aligned',
                                           slot = "cell.embeddings")[,1:ndim]
    } else if(expn_type == 'raw') {
      mydata <- t(as.matrix(seurat_obj@raw.data))
    } else {
      stop('Error: expn_type must be either "raw" or "aligned"')
    }
    for(i in 1:max(state_labels)) {
      ref_clusters[[i]] <- mydata[state_labels == i,]
    }
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
#' centroid_names <- identify_centroids(ref_clusters)
#' }

identify_centroids <- function(ref_clusters) {
  centroids = list()
  for(i in 1:length(ref_clusters)) {
    cur_cluster = ref_clusters[[i]]
    arith_centroid = colMeans(cur_cluster)

    curdist = t(as.matrix(apply(cur_cluster, 1, function(x) norm(x-arith_centroid, type="2"))))
    closest_cell = colnames(curdist)[which.min(curdist)] #closest cell to arithmetic centroid
    #closest_vertex <- get_closest_vertex(closest_cell) #closest cell in mst
    #centroids[[i]] <- closest_vertex
    centroids[[i]] <- closest_cell
  }
  return(centroids)
}


#' @title Get arithmetic centroids (coordinates)
#' @description Takes initial list and returns a matrix with row \var{i} representing the arithmetic centroid of cluster \var{i}
#' @details Private method (not exported in namespace)
#' @param ref_clusters list containing each cluster of interest (each list element is a matrix of dimension num_cells x num_markers)
#' @return Matrix of dimension num_cluster x num_markers; row \var{i} representing the arithmetic centroid of cluster \var{i}
#' @examples
#' \dontrun{
#' cluster_centroids <- get_arithmetic_centroids(ref_clusters)
#' }

get_arithmetic_centroids <- function(ref_clusters) {
  if(length(ref_clusters) < 1) stop('Error: input requires at least 1 reference cluster')
  centroids = matrix(0, nrow=length(ref_clusters), ncol = ncol(ref_clusters[[1]]))
  for(i in 1:length(ref_clusters)) {
    cur_cluster = ref_clusters[[i]]
    centroids[i,] = colMeans(cur_cluster)
  }
  return(centroids)
}

#' @title Assign cells to a reference cell subtype
#' @description Assigns each cell in \code{cur_cells} to a cluster based on nearest cell in Monocle 2 tree
#' @details Private method (not exported in namespace). Uses RANN package for fast knn search
#' @param cur_cells Matrix of cells to be assigned to clusters (Dim: \var{num_cells} x \var{num_markers})
#' @param ref_cells Matrix of cells used to build reference Monocle 2 tree (Dim: \var{num_monocle_cells} x \var{num_markers})
#' @param ref_cell_labels Vector of length \var{num_monocle_cells} containing Monocle 2 cell branch assignments
#' @param data_model Either "monocle2" or "seurat" depending on method used to model cell state space
#' @return Vector of length \var{num_cells} representing cluster assignments for each cell in \var{cur_cells}
#' @examples
#' \dontrun{
#' cur_cells_cluster_labels <- assign_cell_cluster_nearest_node(cur_cells_expn_data, clustered_cells_expn_data, clustered_cells_cluster_labels, data_model='monocle2')
#' }
assign_cell_cluster_nearest_node <- function(cur_cells, ref_cells, ref_cell_labels, data_model='monocle2') {
  if(nrow(ref_cells) != length(ref_cell_labels)) stop("Error: number of cells and cell labels do not match")

  closest <- RANN::nn2(data = ref_cells, query = cur_cells, k = 1) #fast nearest neighbor search
  nearest_cell <- closest$nn.idx

  if(data_model == 'monocle2') {
    assigned <- ref_cell_labels[nearest_cell]
  } else if(data_model == 'seurat') {
    nearest_cell_names <- rownames(ref_cells)[nearest_cell]
    assigned <- as.numeric(ref_cell_labels[nearest_cell_names])
  } else {
    stop('Error: data_model must be either monocle2 or seurat')
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
gaussianff_local <- function(dispersion = 0, parallel = FALSE, zero = NULL) {
  if (!VGAM::is.Numeric(dispersion, length.arg = 1) || dispersion <
    0) {
    stop("bad input for argument 'dispersion'")
  }
  estimated.dispersion <- dispersion == 0
  new("vglmff",
    blurb = c(
      "Vector linear/additive model\n",
      "Links:    identitylink for Y1,...,YM"
    ), constraints = eval(substitute(expression({
      constraints <- VGAM::cm.VGAM(matrix(1, M, 1),
        x = x, bool = .parallel,
        constraints = constraints
      )
      constraints <- VGAM::cm.zero.VGAM(constraints,
        x = x, .zero,
        M = M, predictors.names = predictors.names, M1 = 1
      )
    }), list(.parallel = parallel, .zero = zero))), deviance = function(mu,
                                                                            y, w, residuals = FALSE, eta, extra = NULL) {
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
    }, infos = eval(substitute(function(...) {
      list(
        M1 = 1, Q1 = 1, expected = TRUE, multipleResponses = TRUE,
        quasi.type = TRUE, zero = .zero
      )
    }, list(.zero = zero))), initialize = eval(substitute(expression({
      if (is.R()) assign("CQO.FastAlgorithm", TRUE, envir = VGAM::VGAMenv) else CQO.FastAlgorithm <<- TRUE
      if (any(function.name == c("cqo", "cao")) && (length(.zero) ||
        (is.logical(.parallel) && .parallel))) {
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
    }), list(.parallel = parallel, .zero = zero))), linkinv = function(eta,
                                                                           extra = NULL) eta, last = eval(substitute(expression({
      dy <- dimnames(y)
      if (!is.null(dy[[2]])) dimnames(fit$fitted.values) <- dy
      dpar <- .dispersion
      if (!dpar) {
        wz <- VGAM:::VGAM.weights.function(w = w, M = M, n = n)
        temp5 <- ResSS.vgam(y - mu, wz = wz, M = M)
        dpar <- temp5 / (length(y) - (if (is.numeric(ncol(X.vlm.save))) ncol(X.vlm.save) else 0))
      }
      misc$dispersion <- dpar
      misc$default.dispersion <- 0
      misc$estimated.dispersion <- .estimated.dispersion
      misc$link <- rep_len("identitylink", M)
      names(misc$link) <- predictors.names
      misc$earg <- vector("list", M)
      for (ilocal in 1:M) misc$earg[[ilocal]] <- list()
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
    }), list(.dispersion = dispersion, .estimated.dispersion = estimated.dispersion))),
    loglikelihood = function(mu, y, w, residuals = FALSE,
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
            logretval <- -0.5 * temp1 - n * (M / 2) * log(2 *
              pi)
            for (ii in 1:n) {
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
    }, linkfun = function(mu, extra = NULL) mu, vfamily = "gaussianff",
    validparams = eval(substitute(function(eta, y, extra = NULL) {
      okay1 <- all(is.finite(eta))
      okay1
    }, list(.zero = zero))), deriv = expression({
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

#' @title Create 'phemdObj' object
#' @description Wrapper function to create 'phemdObj' object containing raw expression data and metadata
#' @details Note that each element in list can have different number of rows (i.e. number of cells in each sample can vary).
#' @param data List of length \var{num_samples} containing expression data; each element is of size \var{num_cells} x \var{num_markers}. Alternately a SingleCellExperiment object.
#' @param markers Vector containing marker names (i.e. column names of \code{all_data})
#' @param snames Vector containing sample names (i.e. names of samples contained in \code{all_data})
#' @param datatype Either "list" or "sce" (SingleCellExperiment with genes x cells)
#' @param valtype Type of assay data (i.e. "counts", "normcounts", "logcounts", "tpm", "cpm") if datatype is "sce"
#' @return 'phemdObj' object containing raw multi-sample expression data and associated metadata
#' @examples
#' my_phemdObj <- create_dataobj(all_expn_data, all_genes, as.character(snames))
#' 
create_dataobj <- function(data, markers, snames, datatype='list', valtype='counts') {
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
  data_obj <- new('phemdObj', data = all_data, markers = markers, snames = snames, monocle_obj=NULL)

  return(data_obj)
}


#' @title Attach 'seurat' object to 'phemdObj' object
#' @description Allows user to attach batch-normalized reference cell data from Seurat into 'phemdObj' object containing raw expression data and metadata
#' @param phemd_obj phemdObj struct initialized using create_dataobj
#' @param seurat_obj S4 'seurat' object containing batch-normalized reference cell data
#' @param batch.colname Name of column in Seurat object that denotes batch ID
#' @return 'phemdObj' object containing with attached Seurat object
#' @examples
#' \dontrun{
#' my_phemdObj <- attach_seuratobj(my_phemdObj, my_seuratObj)
#' }
attach_seuratobj <- function(phemd_obj, seurat_obj, batch.colname='plt') {
  stopifnot(is(seurat_obj,'seurat'))
  # ensure cluster names are 1-indexed
  if(min(as.numeric(as.character(seurat_obj@ident))) == 0) {
    label_names <- names(seurat_obj@ident)
    labels_renumbered <- factor(as.numeric(as.character(seurat_obj@ident)) +1)
    names(labels_renumbered) <- label_names
    seurat_obj@ident <- labels_renumbered
  }
  if(batch.colname != 'plt') {
    seurat_obj@meta.data$plt <- seurat_obj@meta.data[[batch.colname]]
  }
  phemd_obj@seurat_obj <- seurat_obj

  return(phemd_obj)
}

#' @title Remove samples with too few cells
#' @description Removes samples from phemdObj that have fewer cells than \code{min_sz}
#' @details Note: If used, this function must be called before (and not after) the \code{aggregate_samples} function is called
#' @param obj 'phemdObj' object containing raw expression data and associated metadata
#' @param min_sz Minimum number of cells in each sample to be retained
#' @return 'phemdObj' object containing raw multi-sample expression data and associated metadata (same as input minus removed samples)
#' @examples
#' my_phemdObj <- create_dataobj(all_expn_data, all_genes, as.character(snames))
#' my_phemdObj_lg <- remove_tiny_samples(my_phemdObj, 10) #removes samples with fewer than 10 cells
#' 
remove_tiny_samples <- function(obj, min_sz=20) {
  stopifnot(is(obj,'phemdObj'))
  stopifnot(mode(min_sz) == 'numeric')
  all_data <- obj@data
  all_snames <- obj@snames
  all_sample_sz <- get_sample_sizes(all_data)
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
  return(obj)
}

#' @title Aggregate expression data from all samples
#' @description Takes initial phemdObj and returns phemdObj with additional data frame in slot @@data_aggregate containing cells aggregated from all samples (to be used for further analyses e.g. Monocle 2 trajectory building / pseudotime mapping / cell clustering)
#' @details Subsamples cells as necessary based on \code{max_cells}. If subsampling is performed, an equal number of cells are subsampled from each sample
#' @param obj 'phemdObj' object containing raw expression data and associated metadata
#' @param max_cells Maximum number of cells across all samples to be included in final matrix on which Monocle 2 will be run
#' @return Same as input 'phemdObj' object with additional slot 'data_aggregate' containing aggregated expression data (num_markers x num_cells)
#' @examples
#' my_phemdObj <- create_dataobj(all_expn_data, all_genes, as.character(snames))
#' my_phemdObj_lg <- remove_tiny_samples(my_phemdObj, 10)
#' my_phemdObj_lg <- aggregate_samples(my_phemdObj_lg, max_cells=1000)
#' 
aggregate_samples <- function(obj, max_cells=12000) {
  stopifnot(is(obj, 'phemdObj'))
  stopifnot(mode(max_cells) == 'numeric')
  set.seed(112) # for reproducibility
  all_data <- obj@data
  nsample <- length(all_data)
  if(nsample == 0) return(obj)

  subsample_sz <- floor(max_cells/nsample)
  all_aggregate_data <- matrix(nrow=0,ncol=ncol(all_data[[1]]))
  all_subsample_idx <- list()
  subsample_bool = FALSE

  all_sample_sz <- get_sample_sizes(all_data)
  if(sum(all_sample_sz) > max_cells) subsample_bool = TRUE

  for(i in seq_len(nsample)) {
    cur_data <- all_data[[i]]
    # take all cells unless total cells across all samples > max_cells
    if(subsample_bool) {
      cur_subsample_idx <- sample(1:nrow(cur_data), min(subsample_sz, nrow(cur_data)))
      all_subsample_idx[[i]] <- cur_subsample_idx
      all_aggregate_data <- rbind(all_aggregate_data, cur_data[cur_subsample_idx,])
    } else {
      all_aggregate_data <- rbind(all_aggregate_data, cur_data)
    }
  }
  all_aggregate_data <- t(all_aggregate_data) #rows = markers, cols = cells
  colnames(all_aggregate_data) <- 1:ncol(all_aggregate_data)
  rownames(all_aggregate_data) <- obj@markers

  obj@data_aggregate <- as.matrix(all_aggregate_data)
  obj@data_subsample_idx = all_subsample_idx
  obj@subsampled_bool <- subsample_bool
  return(obj)
}

#' @title Perform feature selection on aggregated data
#' @description Takes as input an phemdObj with aggregated data and returns same phemdObj after performing feature selection on aggregated data
#' @details \code{aggregate_samples} needs to be called before running this function
#' @param obj 'phemdObj' object containing aggregated data
#' @param selected_genes Vector containing names of genes to use for downstream analyses
#' @return Same as input 'phemdObj' object after performing feature-selection based dimensionality reduction on aggregated expression data
#' @examples
#' 
#' my_phemdObj <- create_dataobj(all_expn_data, all_genes, as.character(snames))
#' my_phemdObj_lg <- remove_tiny_samples(my_phemdObj, 10)
#' my_phemdObj_lg <- aggregate_samples(my_phemdObj_lg, max_cells=1000)
#' my_phemdObj_lg <- select_features(my_phemdObj_lg, selected_genes=c('TP53', 'EGFR', 'KRAS'))
#' 
select_features <- function(obj, selected_genes) {
  if(isempty(obj@data_aggregate)) stop('slot "data_aggregate" is empty; please call aggregate_samples() before running select_features()')
  all_aggregate_data <- obj@data_aggregate
  all_genes <- obj@markers
  selected_gene_map <- match(selected_genes, all_genes)
  #TODO: print genes in selected_genes that were unable to map to all_genes

  all_aggregate_data <- all_aggregate_data[selected_gene_map,]
  obj@data_aggregate <- all_aggregate_data
  obj@markers <- all_genes[selected_gene_map]
  return(obj)
}

#' @title Generate Monocle2 embedding
#' @description Takes as input an phemdObj with aggregated data and returns same phemdObj with Monocle2 object in @@monocle_obj slot
#' @details Wrapper function for \code{reduceDimension} in Monocle 2 package. \code{aggregate_samples} needs to be called before running this function.
#' @param obj 'phemdObj' object containing aggregated data
#' @param data_model One of the following: 'negbinomial_sz', 'negbinomial', 'tobit', 'gaussianff'. See "Family Function" table at the following link for more details on selecting the proper one. \url{http://cole-trapnell-lab.github.io/monocle-release/docs/#getting-started-with-monocle}
#' @param ... Additional parameters to be passed to \code{reduceDimension} function
#' @return Same as input 'phemdObj' object with additional Monocle2 object in @@monocle_obj slot
#' @examples
#' my_phemdObj <- create_dataobj(all_expn_data, all_genes, as.character(snames))
#' my_phemdObj_lg <- remove_tiny_samples(my_phemdObj, 10)
#' my_phemdObj_lg <- aggregate_samples(my_phemdObj_lg, max_cells=1000)
#' my_phemdObj_lg <- embed_cells(my_phemdObj_lg, data_model = 'gaussianff', sigma=0.02, maxIter=2)
embed_cells <- function(obj, data_model = 'negbinomial_sz', ...) {
  extra_args <- list(...)

  if(isempty(obj@data_aggregate)) stop('slot "data_aggregate" is empty; please call aggregate_samples() before running embed_cells()')
  mydata <- obj@data_aggregate
  if(is.null(mydata)) stop("Error: call 'aggregate_samples' function first ")

  myFeatureData <- as.data.frame(obj@markers)
  fd <- new("AnnotatedDataFrame", data = myFeatureData)
  rownames(fd) <- obj@markers

  if(is.null(data_model)) {
    print('Assuming data fit negative binomial pattern of expression...')
    data_model = 'negbinomial_sz'
  }

  if(data_model == 'negbinomial_sz') {
    expression_fam_fn = VGAM::negbinomial.size()
  } else if(data_model == 'negbinomial') {
    expression_fam_fn = VGAM::negbinomial()
  } else if(data_model == 'tobit') {
    expression_fam_fn = VGAM::tobit()
  } else if(data_model == 'gaussianff') {
    #expression_fam_fn = VGAM::gaussianff()
    expression_fam_fn = gaussianff_local()
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

  monocle_obj <- newCellDataSet(mydata,phenoData=NULL,featureData = fd, expressionFamily = expression_fam_fn)
  varLabels(monocle_obj@featureData) <- 'gene_short_name' #random formatting requirement for monocle

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
  #rd_args <- c(list(cds=monocle_obj, reduction_method='tSNE')) #NEW FOR MONOCLE 3
  monocle_obj_red <- do.call(reduceDimension, rd_args)

  obj@monocle_obj <- monocle_obj_red
  return(obj)
}

#' @title Compute Monocle2 cell state and pseudotime assignments
#' @description Takes as input an phemdObj with Monocle2 object and returns same phemdObj with Monocle2 object containing cell state and pseudotime assignments
#' @details Wrapper function for \code{orderCells} in Monocle 2 package. \code{embed_cells} needs to be called before calling this function.
#' @param obj 'phemdObj' object containing initial Monocle 2 object
#' @param ... Additional parameters to be passed into \code{orderCells} function
#' @return Same as input 'phemdObj' object with updated Monocle2 object in @@monocle_obj slot containing cell state and pseudotime assignments
#' @examples
#' 
#' my_phemdObj <- create_dataobj(all_expn_data, all_genes, as.character(snames))
#' my_phemdObj_lg <- remove_tiny_samples(my_phemdObj, 10)
#' my_phemdObj_lg <- aggregate_samples(my_phemdObj_lg, max_cells=1000)
#' my_phemdObj_monocle <- embed_cells(my_phemdObj_lg, data_model = 'gaussianff', sigma=0.02, maxIter=2)
#' my_phemdObj_monocle <- order_cells(my_phemdObj_monocle)
order_cells <- function(obj, ...) {
  monocle_obj <- obj@monocle_obj
  if(ncol(monocle_obj) == 0) stop('slot "monocle_obj" is empty; please call embed_cells() before running order_cells()')

  extra_args <- list(...)
  oc_args <- c(list(cds=monocle_obj),
               extra_args[names(extra_args) %in% c("root_state", "reverse")])

  monocle_obj_ordered <- do.call(orderCells, oc_args)
  obj@monocle_obj <- monocle_obj_ordered
  return(obj)
}

#' @title Save Monocle2 cell embedding plots to folder
#' @description Takes as input an phemdObj with Monocle2 object (already embedded and ordered) and saves cell embedding plots to specified folder
#' @details \code{embed_cells} and \code{order_cells} need to be called before calling this function. Required additional package: 'RColorBrewer'
#' @param obj 'phemdObj' object containing Monocle 2 object
#' @param path Path to destination folder (must already exist)
#' @param cell_model Method by which cell state was modeled (either "monocle2" or "seurat")
#' @param cmap User-specified colormap to use to color cell state embedding (optional)
#' @param w Width of plot in inches
#' @param h Height of plot in inches
#' @param pt_sz Scalar factor for point size
#' @param ndims Number of dimensions to use for dimensionality reduction in case it hasn't been performed yet (only relevant when using Seurat data as input)
#' @return Colormap (vector of colors) used to color Monocle2 cell state embedding
#' @examples
#' my_phemdObj <- create_dataobj(all_expn_data, all_genes, as.character(snames))
#' my_phemdObj_lg <- remove_tiny_samples(my_phemdObj, 10)
#' my_phemdObj_lg <- aggregate_samples(my_phemdObj_lg, max_cells=1000)
#' my_phemdObj_monocle <- embed_cells(my_phemdObj_lg, data_model = 'gaussianff', sigma=0.02, maxIter=2)
#' my_phemdObj_monocle <- order_cells(my_phemdObj_monocle)
#' cmap <- plot_embeddings(my_phemdObj_monocle, path='.')
plot_embeddings <- function(obj, path, cell_model='monocle2', cmap=NULL, w=4, h=5, pt_sz = 1, ndims=NULL) {
  if(substr(path,nchar(path), nchar(path)) != '/') path <- paste(path, '/', sep='') #ensure path ends with a slash
  if(cell_model == 'monocle2') {
    monocle_obj <- obj@monocle_obj
    cell_embedding = reducedDimS(monocle_obj)
    mydata <- obj@data_aggregate

    # Extract state labels from monocle data object
    labels <- pData(phenoData(monocle_obj))
    state_labels <- as.numeric(labels$State)


    levels <- levels(factor(state_labels))
    levels_renamed <- sapply(levels, function(x) paste("C-", x, sep=""))

    if(is.null(cmap)) {
      getPalette <- colorRampPalette(brewer.pal(11, "Spectral"))
      cmap <- getPalette(max(state_labels))
      if("#FFFFBF" %in% cmap) cmap[which(cmap == "#FFFFBF")] <- "#D3D3D3" #replace light yellow with grey
      cmap <- sample(cmap)
    }
    palette(cmap)

    cell_embedding_t <- as.data.frame(t(cell_embedding))
    # visualize traj colored by state
    myplot <- ggplot(cell_embedding_t, aes(x=cell_embedding_t[,1], y=cell_embedding_t[,2], color=factor(state_labels))) +
      geom_point(size=0.4) +
      scale_color_manual(labels = levels_renamed,
                         values = cmap) +
      guides(colour = guide_legend(override.aes = list(size=2))) +
      labs(x="", y = "", color = "Cell subtype") +
      theme_classic() +
      theme(axis.title.x=element_blank(),
            axis.text.x=element_blank(),
            axis.ticks.x=element_blank(),
            axis.title.y=element_blank(),
            axis.text.y=element_blank(),
            axis.ticks.y=element_blank(),
            axis.line = element_line(colour = "black",
                                     size = 1, linetype = "solid"))

    ggsave(filename=paste(path, "traj_state.png", sep=""), plot=myplot, width=w, height=h, dpi=300)

    png(filename=paste(path, "traj_state_labeled.png", sep=""),
        units="in",
        width=w,
        height=h,
        res=300)
    plot(cell_embedding[1,], cell_embedding[2,], col=state_labels, pch=20, cex=1, xlab="", ylab="", xaxt='n', yaxt='n')
    ref_cluster_centroids <- matrix(0, nrow=max(state_labels), ncol=nrow(cell_embedding))
    ref_clusters = list()
    for(i in 1:max(state_labels)) {
      ref_clusters[[i]] <- t(monocle_obj[,state_labels == i]@assayData$exprs)
      ref_cluster_centroids[i,] = rowMeans(cell_embedding[,state_labels == i])
    }
    text(ref_cluster_centroids[,1], ref_cluster_centroids[,2], cex=1.5, col='black')
    dev.off()

    # visualize traj colored by pseudotime
    ncolor = 9
    palette(brewer.pal(ncolor, "Blues"))
    col.labels <- labels$Pseudotime / max(labels$Pseudotime) * (ncolor - 2) + 2 # linear interpolation between [2,ncolor]

    png(filename=paste(path, "traj_pseudotime_labeled.png", sep=""),
        units="in",
        width=w,
        height=h,
        res=300)
    plot(cell_embedding[1,], cell_embedding[2,], col=col.labels, pch=20, cex=1, xlab="", ylab="", xaxt='n', yaxt='n')
    ref_cluster_centroids <- matrix(0, nrow=max(state_labels), ncol=nrow(cell_embedding))
    ref_clusters = list()
    for(i in 1:max(state_labels)) {
      ref_clusters[[i]] <- t(mydata[,state_labels == i])
      ref_cluster_centroids[i,] = rowMeans(cell_embedding[,state_labels == i])
    }
    text(ref_cluster_centroids[,1], ref_cluster_centroids[,2], cex=1.5, col='black')
    dev.off()

    col.labels <- labels$Pseudotime

    # visualize traj colored by pseudotime
    myplot <- ggplot(cell_embedding_t, aes(x=cell_embedding_t[,1], y=cell_embedding_t[,2], color=col.labels)) +
      geom_point(size=0.4) +
      labs(x="", y = "", color = "Pseudotime") +
      theme_classic() +
      theme(axis.title.x=element_blank(),
            axis.text.x=element_blank(),
            axis.ticks.x=element_blank(),
            axis.title.y=element_blank(),
            axis.text.y=element_blank(),
            axis.ticks.y=element_blank(),
            axis.line = element_line(colour = "black",
                                     size = 1, linetype = "solid"))

    ggsave(filename=paste(path, "traj_pseudotime.png", sep=""), plot=myplot, width=w, height=h, dpi=300)
    return(cmap)
  } else if(cell_model == 'seurat') {
    # TODO: implement this
    seurat_obj <- obj@seurat_obj

    if(!'tsne' %in% names(seurat.combined@dr)) {
      print('Running t-SNE...')
      if(is.null(ndims)) ndims = 10
      seurat_obj <- RunTSNE(seurat_obj,reduction.use = "cca.aligned", dims.use = 1:ndims)
    }

    # define color map
    if(is.null(cmap)) {
      getPalette <- colorRampPalette(brewer.pal(11, "Spectral"))
      cmap <- getPalette(max(as.numeric(seurat.combined@ident)))
      if("#FFFFBF" %in% cmap) cmap[which(cmap == "#FFFFBF")] <- "#D3D3D3" #replace light yellow with grey
      cmap <- sample(cmap)
    }


    png(filename=paste(path, "cell_state_tsne.png", sep=""),
        units="in",
        width=w,
        height=h,
        res=300)
    TSNEPlot(seurat_obj, do.label = FALSE, pt.size = pt_sz, colors.use = cmap)
    dev.off()

  } else {
    stop('Error: cell_model must be either "monocle2" or "seurat"')
  }
  return(cmap)
}


#' @title Save Monocle2 heatmap plot to folder
#' @description Takes as input an phemdObj with Monocle2 object (already embedded and ordered) and saves heatmap describing cell subtypes to specified folder
#' @details \code{embed_cells} and \code{order_cells} need to be called before calling this function. Required additional package: 'pheatmap'
#' @param obj 'phemdObj' object containing Monocle 2 object
#' @param path Path to destination folder (must already exist)
#' @param cell_model Method by which cell state was modeled (either "monocle2" or "seurat")
#' @param selected_genes Vector containing gene names to include in heatmap (optional)
#' @param w Width of plot in inches
#' @param h Height of plot in inches
#' @return None
#' @examples
#' my_phemdObj <- create_dataobj(all_expn_data, all_genes, as.character(snames))
#' my_phemdObj_lg <- remove_tiny_samples(my_phemdObj, 10)
#' my_phemdObj_lg <- aggregate_samples(my_phemdObj_lg, max_cells=1000)
#' my_phemdObj_monocle <- embed_cells(my_phemdObj_lg, data_model = 'gaussianff', sigma=0.02, maxIter=2)
#' my_phemdObj_monocle <- order_cells(my_phemdObj_monocle)
#' cmap <- plot_embeddings(my_phemdObj_monocle, path='.')
#' plot_heatmaps(my_phemdObj_monocle, path='.')
plot_heatmaps <- function(obj, path, cell_model='monocle2', selected_genes=NULL, w=8, h=5) {
  if(substr(path,nchar(path), nchar(path)) != '/') path <- paste(path, '/', sep='') #ensure path ends with a slash
  if(cell_model == 'monocle2') {
    # retrieve reference clusters
    ref_clusters <- retrieve_refclusters(obj, cell_model='monocle2')
    selected_clusters <- 1:length(ref_clusters)
    myheatmap <- matrix(0, nrow=length(selected_clusters), ncol = ncol(ref_clusters[[1]]))
    for(i in 1:length(selected_clusters)) {
      cur_cluster_idx = selected_clusters[i]
      cur_cluster = ref_clusters[[i]]
      if(!is.null(cur_cluster)) myheatmap[i,] = colMeans(cur_cluster)
    }

    selected_clusters_renamed <- sapply(selected_clusters, function(x) paste("C-", x, sep=""))

    rownames(myheatmap) <- selected_clusters_renamed
    colnames(myheatmap) <- obj@markers

    if(!is.null(selected_genes)) {
      col_tokeep <- match(selected_genes, obj@markers)
      if(sum(is.na(col_tokeep)) > 0) {
        genes_not_found <- ''
        missing_idx <- which(is.na(col_tokeep))
        for(i in 1:length(missing_idx)) {
          if(i == 1) genes_not_found <- paste(genes_not_found, selected_genes[missing_idx[i]], sep='')
          else genes_not_found <- paste(genes_not_found, selected_genes[missing_idx[i]], sep=', ')
        }
        print(sprintf("Genes not found: %s", genes_not_found, sep=""))
      }
      col_tokeep <- col_tokeep[!is.na(col_tokeep)]
      myheatmap <- myheatmap[,col_tokeep]
    }

    myheatmap[is.nan(myheatmap)] <- 0 #this in the event of empty clusters

    assignInNamespace(
      x = "draw_colnames",
      value = "draw_colnames_45",
      ns = asNamespace("pheatmap")
    )

    myheatmap2 <- log2(myheatmap - min(myheatmap) + 1)

    pheatmap(myheatmap2,
             cluster_rows=FALSE,
             cluster_cols=TRUE,
             border_color=NA,
             show_colnames=TRUE,
             show_rownames=TRUE,
             fontsize_col=8,
             fontsize_row=12,
             cellwidth=10,
             filename=paste(path, "heatmap.png", sep=""),
             width=w,
             height=h
    )
  } else if(cell_model == 'seurat') {
    #TODO: implement this
    seurat_obj <- obj@seurat_obj
    state_labels <- as.numeric(as.character(seurat_obj@ident))
    names(state_labels) <- rownames(seurat_obj@meta.data)
    ref_data <- t(as.matrix(seurat_obj@raw.data))

    batches <- unique(obj@experiment_ids)
    myheatmaps_all <- list()
    for(batch_id in batches) {
      cell_idx_curplt <- which(seurat_obj@meta.data$plt == batch_id)
      if(length(cell_idx_curplt) == 0) {
        stop(sprintf('Error: no cells in reference set match the experiment_id %s. Please check phemdObj@experiment_ids.', batch_id))
      }
      cur_ref_data <- ref_data[cell_idx_curplt,]
      cur_state_labels <- state_labels[cell_idx_curplt]


      myheatmap <- matrix(0, nrow=max(state_labels), ncol = ncol(cur_ref_data))
      for(i in 1:max(state_labels)) {
        cur_idx <- which(cur_state_labels == i)
        cur_cluster = cur_ref_data[cur_idx,]
        if(length(cur_idx) > 1) myheatmap[i,] = colMeans(cur_cluster)
        if(length(cur_idx) == 1) myheatmap[i,] = cur_cluster
      }

      selected_clusters_renamed <- sapply(1:max(state_labels), function(x) paste("C-", x, sep=""))

      rownames(myheatmap) <- selected_clusters_renamed
      colnames(myheatmap) <- obj@markers

      if(!is.null(selected_genes)) {
        col_tokeep <- match(selected_genes, obj@markers)
        if(sum(is.na(col_tokeep)) > 0) {
          genes_not_found <- ''
          missing_idx <- which(is.na(col_tokeep))
          for(i in 1:length(missing_idx)) {
            if(i == 1) genes_not_found <- paste(genes_not_found, selected_genes[missing_idx[i]], sep='')
            else genes_not_found <- paste(genes_not_found, selected_genes[missing_idx[i]], sep=', ')
          }
          print(sprintf("Genes not found: %s", genes_not_found, sep=""))
        }
        col_tokeep <- col_tokeep[!is.na(col_tokeep)]
        myheatmap <- myheatmap[,col_tokeep]
      }

      myheatmap[is.nan(myheatmap)] <- 0 #this in the event of empty clusters

      assignInNamespace(
        x = "draw_colnames",
        value = "draw_colnames_45",
        ns = asNamespace("pheatmap")
      )

      myheatmap2 <- log2(myheatmap - min(myheatmap) + 1)
      myheatmaps_all[[batch_id]] <- myheatmap2
      pheatmap(myheatmap2,
               cluster_rows=FALSE,
               cluster_cols=FALSE,
               border_color=NA,
               show_colnames=TRUE,
               show_rownames=TRUE,
               fontsize_col=8,
               fontsize_row=12,
               cellwidth=10,
               filename=paste(path, sprintf("heatmap_%s.png", batch_id), sep=""),
               width=w,
               height=h
      )
    }

    myheatmaps_avg <- myheatmaps_all[[1]]
    for(i in 2:length(myheatmaps_all)) {
      myheatmaps_avg <- myheatmaps_avg + myheatmaps_all[[i]]
    }
    myheatmaps_avg <- myheatmaps_avg / length(myheatmaps_all)
    pheatmap(myheatmaps_avg,
             cluster_rows=TRUE,
             cluster_cols=FALSE,
             border_color=NA,
             show_colnames=TRUE,
             show_rownames=TRUE,
             fontsize_col=8,
             fontsize_row=12,
             cellwidth=10,
             filename=paste(path, "heatmap_avg.png", sep=""),
             width=w,
             height=h)

  } else {
    stop('Error: cell_model must be either "monocle2" or "seurat"')
  }
}


#' @title Rotates heatmap marker labels 45 degrees
#' @description Overwrites default draw_colnames in the pheatmap package
#' @details To be used with pheatmap plotting function; not to be called directly. Thanks to Josh O'Brien at http://stackoverflow.com/questions/15505607
#' @param coln Column names
#' @param gaps Spacing of labels
#' @param ... Additional parameters to be passed to \code{gpar}
#' @return Formatted marker labels in heatmap
#' @examples
#' #Not to be called directly
draw_colnames_45 <- function(coln, gaps, ...) {
  coord <- pheatmap:::find_coordinates(length(coln), gaps)
  x     <- coord$coord - 0.5 * coord$size
  res   <- grid::textGrob(
    coln, x = x, y = unit(1, "npc") - unit(3,"bigpts"),
    vjust = 0.75, hjust = 1, rot = 45, gp = grid::gpar(...)
  )
  return(res)
}

#' @title Computes cell subtype abundances for each sample
#' @description Takes as input an phemdObj with all single-cell expression data of all single-cell samples in @@data slot and Monocle2 object (already embedded and ordered) in @@monocle_obj slot. Returns same phemdObj with cell subtype frequencies of each sample in @@data_cluster_weights slot
#' @details \code{embed_cells} and \code{order_cells} need to be called before calling this function.
#' @param obj 'phemdObj' object containing single-cell expression data of all samples in @@data slot and Monocle2 object (already embedded and ordered) in @@monocle_obj slot
#' @param verbose Boolean that determines whether progress (sequential processing of samples) should be printed. FALSE by default
#' @param cell_model Either "monocle2" or "seurat" depending on method used to model cell state space
#' @return phemdObj with cell subtype frequencies of each sample in @@data_cluster_weights slot
#' @examples
#' 
#' my_phemdObj <- create_dataobj(all_expn_data, all_genes, as.character(snames))
#' my_phemdObj_lg <- remove_tiny_samples(my_phemdObj, 10)
#' my_phemdObj_lg <- aggregate_samples(my_phemdObj_lg, max_cells=1000)
#' my_phemdObj_monocle <- embed_cells(my_phemdObj_lg, data_model = 'gaussianff', sigma=0.02, maxIter=2)
#' my_phemdObj_monocle <- order_cells(my_phemdObj_monocle)
#' my_phemdObj_final <- cluster_individual_samples(my_phemdObj_monocle)
#' 
cluster_individual_samples <- function(obj, verbose=FALSE, cell_model='monocle2') {
  stopifnot(is(obj,'phemdObj'))
  all_data <- obj@data
  if(cell_model == 'monocle2') {
    monocle_obj <- obj@monocle_obj
    if(ncol(monocle_obj) == 0) stop('slot "monocle_obj" is empty; please call embed_cells() and order_cells() before calling this function')

    # Extract state labels from monocle data object
    labels <- pData(phenoData(monocle_obj))
    if(!('State' %in% names(labels))) stop('obj@monocle_obj does not have cell state assignments; please call embed_cells() and order_cells() before calling this function')
    state_labels <- as.numeric(labels$State)
    # retrieve reference clusters
    ref_clusters <- retrieve_refclusters(obj, cell_model='monocle2')
    nclusters = length(ref_clusters)

    if(!obj@subsampled_bool) {
      ## subsampling not performed; all cells assigned to clusters as-is

      cluster_weights <- matrix(0, nrow=length(all_data), ncol = nclusters)

      start_idx <- 1
      for(i in 1:length(all_data)) {
        cur_sample_sz <- nrow(all_data[[i]])
        end_idx <- start_idx + cur_sample_sz - 1
        sample_labels <- state_labels[start_idx:end_idx] #cell subtype assignments for current sample
        cur_hist = rep(0, nclusters)
        for(j in 1:nclusters) {
          cur_hist[j] = sum(sample_labels == j);
        }
        cur_hist = cur_hist / sum(cur_hist);
        cluster_weights[i,] = cur_hist;
        start_idx <- end_idx + 1
      }
      obj@data_cluster_weights <- cluster_weights
    } else {
      ## subsampling performed; need to assign cells to cluster of nearest cell in embedding
      refcluster_sizes <- rep(0, length(ref_clusters))
      counter1 = 0
      for(i in 1:nclusters) {
        refcluster_sizes[i] = nrow(ref_clusters[[i]])
        counter1 = counter1 + nrow(ref_clusters[[i]])
      }
      print(refcluster_sizes / counter1)

      cluster_weights <- matrix(0, nrow=length(all_data), ncol = nclusters)
      for(i in 1:length(all_data)) {
        cur_data = all_data[[i]]
        cur_ncells = nrow(cur_data)
        cur_hist = rep(0, nclusters)

        if(verbose && i %% 10 == 0) {
          print(sprintf('Processing sample number: %d', i))
        }
        # Consider using euclidean / mahalanobis distance to centroids of reference clusters
        #cur_cell_labels = assign_cell_cluster(cur_data, ref_clusters);

        cur_cell_labels = assign_cell_cluster_nearest_node(cur_data, t(exprs(monocle_obj)), state_labels, data_model=cell_model); # Use nearest-cell mapping instead of nearest-centroid mapping

        for(j in 1:nclusters) {
          cur_hist[j] = sum(cur_cell_labels == j);
        }
        cur_hist = cur_hist / sum(cur_hist);
        cluster_weights[i,] = cur_hist;
      }
      print(colMeans(cluster_weights))

      obj@data_cluster_weights <- cluster_weights
    }
  } else if(cell_model == 'seurat') {
    seurat_obj <- obj@seurat_obj
    # retrieve reference clusters for starting estimate of cluster sizes
    ref_clusters <- retrieve_refclusters(obj, cell_model='seurat', expn_type = 'raw')
    nclusters <- length(ref_clusters)
    ## subsampling performed; need to assign cells to cluster of nearest cell in embedding
    refcluster_sizes <- rep(0, length(ref_clusters))
    counter1 = 0
    for(i in 1:nclusters) {
      refcluster_sizes[i] = nrow(ref_clusters[[i]])
      counter1 = counter1 + nrow(ref_clusters[[i]])
    }
    print(refcluster_sizes / counter1)

    cluster_weights <- matrix(0, nrow=length(all_data), ncol = nclusters)
    for(i in 1:length(all_data)) {
      cur_data <- all_data[[i]]
      cur_plt <- obj@experiment_ids[i]
      if(verbose && i %% 10 == 0) {
        print(sprintf('Processing sample number: %d', i))
      }

      cur_hist <- rep(0, nclusters)
      state_labels <- as.numeric(as.character(seurat_obj@ident))
      names(state_labels) <- rownames(seurat_obj@meta.data)
      ref_data <- t(as.matrix(seurat_obj@raw.data))
      cell_idx_curplt <- which(seurat_obj@meta.data$plt == cur_plt)
      if(length(cell_idx_curplt) == 0) {
        stop(sprintf('Error: no cells in reference set match the experiment_id %s of sample %d', cur_plt, i))
      }
      ref_data <- ref_data[cell_idx_curplt,]
      state_labels <- state_labels[cell_idx_curplt]

      cur_cell_labels = assign_cell_cluster_nearest_node(cur_data, ref_data, state_labels, data_model=cell_model); # Use nearest-cell mapping instead of nearest-centroid mapping

      for(j in 1:nclusters) {
        cur_hist[j] = sum(cur_cell_labels == j);
      }
      cur_hist = cur_hist / sum(cur_hist);
      cluster_weights[i,] = cur_hist;
    }
    print(colMeans(cluster_weights))

    obj@data_cluster_weights <- cluster_weights
  } else {
    stop('Error: cell_model must be either "monocle2" or "seurat"')
  }

  return(obj)
}


#' @title Computes ground distance matrix based on cell embedding
#' @description Takes as input an phemdObj with Monocle2 object (already embedded and ordered) in @@monocle_obj slot. Returns same phemdObj with ground distance matrix representing pairwise distance between 2 cell subtypes based on cell state embedding.
#' @details \code{embed_cells} and \code{order_cells} need to be called before calling this function. Requires 'igraph' package
#' @param obj 'phemdObj' object containing Monocle2 object (already embedded and ordered) in @@monocle_obj slot
#' @param cell_model Method by which cell state was modeled (either "monocle2" or "seurat")
#' @return phemdObj with ground distance matrix (to be used in EMD computation) in @@data_cluster_weights slot
#' @examples
#' 
#' my_phemdObj <- create_dataobj(all_expn_data, all_genes, as.character(snames))
#' my_phemdObj_lg <- remove_tiny_samples(my_phemdObj, 10)
#' my_phemdObj_lg <- aggregate_samples(my_phemdObj_lg, max_cells=1000)
#' my_phemdObj_monocle <- embed_cells(my_phemdObj_lg, data_model = 'gaussianff', sigma=0.02, maxIter=2)
#' my_phemdObj_monocle <- order_cells(my_phemdObj_monocle)
#' my_phemdObj_final <- cluster_individual_samples(my_phemdObj_monocle)
#' my_phemdObj_final <- generate_gdm(my_phemdObj_final)
#' 
generate_gdm <- function(obj, cell_model='monocle2') {
  stopifnot(is(obj,'phemdObj'))

  if(cell_model == 'monocle2') {
    monocle_obj <- obj@monocle_obj
    # retrieve reference clusters
    ref_clusters <- retrieve_refclusters(obj, cell_model='monocle2')
    nclusters = length(ref_clusters)

    mst_graph <- monocle_obj@auxOrderingData[["DDRTree"]]$pr_graph_cell_proj_tree # get graph underlying Monocle tree
    centroids <- identify_centroids(ref_clusters)
    pseudotimes <- pData(phenoData(monocle_obj))$Pseudotime
    emd_dists <- matrix(0, nrow=nclusters, ncol=nclusters)
    for(i in 1:nclusters) {
      for(j in 1:nclusters) {
        if(i == j) next
        path_btwn_cells <- shortest_paths(mst_graph, centroids[[i]], centroids[[j]])$vpath[[1]] #cell1 and cell2 should be strings representing vertices in mst
        pseudotime_curpath <- pseudotimes[path_btwn_cells]
        pseudotime_dist <- abs(pseudotimes[as.numeric(centroids[[i]])] - min(pseudotime_curpath)) + abs(pseudotimes[as.numeric(centroids[[j]])] - min(pseudotime_curpath))
        emd_dists[i,j] <- pseudotime_dist
      }
    }
    obj@emd_dist_mat <- emd_dists
  } else if(cell_model == 'seurat') {
    # TODO: implement seurat option for computing gdm
    ref_clusters <- retrieve_refclusters(obj, cell_model='seurat', expn_type='aligned', ndim=8)
    seurat_obj <- obj@seurat_obj
    centroids <- get_arithmetic_centroids(ref_clusters)
    obj@emd_dist_mat <- as.matrix(dist(centroids))
  } else {
    stop('Error: cell_model must be either "monocle2" or "seurat"')
  }
  return(obj)
}

#' @title Computes EMD distance matrix representing pairwise dissimilarity between samples
#' @description Takes as input an phemdObj with cell subtype relative frequencies for each sample in @@data_cluster_weights slot and ground distance matrix (representing cell subtype pairwise dissimilarity) in @@emd_dist_mat slot. Returns distance matrix representing pairwise dissimilarity between samples
#' @details Requires 'transport' and 'pracma' packages
#' @param obj 'phemdObj' object containing cell subtype relative frequencies for each sample in @@data_cluster_weights slot and ground distance matrix (representing cell subtype dissimilarity) in @@emd_dist_mat slot
#' @return Distance matrix of dimension num_samples x num_samples representing pairwise dissimilarity between samples
#' @examples
#'
#' my_phemdObj <- create_dataobj(all_expn_data, all_genes, as.character(snames))
#' my_phemdObj_lg <- remove_tiny_samples(my_phemdObj, 10)
#' my_phemdObj_lg <- aggregate_samples(my_phemdObj_lg, max_cells=1000)
#' my_phemdObj_monocle <- embed_cells(my_phemdObj_lg, data_model = 'gaussianff', sigma=0.02, maxIter=2)
#' my_phemdObj_monocle <- order_cells(my_phemdObj_monocle)
#' my_phemdObj_final <- cluster_individual_samples(my_phemdObj_monocle)
#' my_phemdObj_final <- generate_gdm(my_phemdObj_final)
#' my_EMD_mat <- compare_samples(my_phemdObj_final)
#' 
compare_samples <- function(obj) {
  stopifnot(is(obj,'phemdObj'))
  cluster_weights <- obj@data_cluster_weights
  emd_dists <- obj@emd_dist_mat
  # generate inhibitor distance matrix
  Y = rep(0, (nrow(cluster_weights)-1)*nrow(cluster_weights)/2)
  counter = 1
  for(i in 1:nrow(cluster_weights)){
    cur_f1_weights = cluster_weights[i,]
    for(j in (i+1):nrow(cluster_weights)) {
      if(j > nrow(cluster_weights)) break #weird that this doesn't automatically happen in R
      cur_f2_weights = cluster_weights[j,]
      # Compute EMD
      flow = transport(cur_f1_weights, cur_f2_weights, emd_dists, method='primaldual')
      curdist=0
      for(k in 1:nrow(flow)) {
        cur_penalty = emd_dists[flow[k,1], flow[k,2]]
        curdist=curdist+cur_penalty*flow[k,3]
      }
      Y[counter] = curdist
      counter = counter + 1
    }
  }
  Y_sq = squareform(Y)
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
#' my_phemdObj <- create_dataobj(all_expn_data, all_genes, as.character(snames))
#' my_phemdObj_lg <- remove_tiny_samples(my_phemdObj, 10)
#' my_phemdObj_lg <- aggregate_samples(my_phemdObj_lg, max_cells=1000)
#' my_phemdObj_monocle <- embed_cells(my_phemdObj_lg, data_model = 'gaussianff', sigma=0.02, maxIter=2)
#' my_phemdObj_monocle <- order_cells(my_phemdObj_monocle)
#' my_phemdObj_final <- cluster_individual_samples(my_phemdObj_monocle)
#' my_phemdObj_final <- generate_gdm(my_phemdObj_final)
#' my_EMD_mat <- compare_samples(my_phemdObj_final)
#' cluster_assignments <- group_samples(my_EMD_mat, distfun = 'hclust', ncluster=4)
#' 
group_samples <- function(distmat, distfun = 'hclust', ncluster=NULL, method='complete') {
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
    if(is.null(ncluster)) ncluster = 4
    cluster_results <- pam(distmat, ncluster, diss=TRUE)
    cluster_assignments=cluster_results$clustering
  } else {
    stop("Error: Please specify distfun as either 'hclust' or 'pam'")
  }
  return(cluster_assignments)
}

#' @title Writes samples to file based on community detection group assignments
#' @description Takes vector of cluster assignments and phemdObj containing sample names and writes sample groups to file
#' @details Order of samples in obj@@snames is assumed to be the same as the order of group assignments in cluster_assignments
#' @param cluster_assignments Vector containing group assignments for each sample
#' @param obj phemdObj object containing sample names in @@snames slot
#' @param dest Path to existing directory where output should be saved
#' @return None
#' @examples
#' 
#' my_phemdObj <- create_dataobj(all_expn_data, all_genes, as.character(snames))
#' my_phemdObj_lg <- remove_tiny_samples(my_phemdObj, 10)
#' my_phemdObj_lg <- aggregate_samples(my_phemdObj_lg, max_cells=1000)
#' my_phemdObj_monocle <- embed_cells(my_phemdObj_lg, data_model = 'gaussianff', sigma=0.02, maxIter=2)
#' my_phemdObj_monocle <- order_cells(my_phemdObj_monocle)
#' my_phemdObj_final <- cluster_individual_samples(my_phemdObj_monocle)
#' my_phemdObj_final <- generate_gdm(my_phemdObj_final)
#' my_EMD_mat <- compare_samples(my_phemdObj_final)
#' cluster_assignments <- group_samples(my_EMD_mat, distfun = 'hclust', ncluster=4)
#' print_cluster_assignments(cluster_assignments, my_phemdObj_final, '.')
#' 
print_cluster_assignments <- function(cluster_assignments, obj, dest) {
  snames <- obj@snames
  unlink(paste(dest, 'sample_groups', sep=''), recursive=TRUE)
  dir.create(file.path(paste(dest, 'sample_groups', sep='')), showWarnings = FALSE) # create folder for output
  for(i in 1:max(cluster_assignments)) {
    cur_cluster_idx <- which(cluster_assignments == i)
    cur_file <- sprintf('sample_groups/scluster_%s.txt', intToUtf8(64+i));
    cur_file <- strcat(dest, cur_file)
    write(snames[cur_cluster_idx], file=cur_file, sep="\n")
  }
}

#' @title Plot diffusion map embedding of samples based on distance matrix
#' @description Visualizes diffusion map for network of samples based on square distance matrix (sample-sample pairwise dissimilarity)
#' @details Requires 'destiny' package
#' @param my_distmat phemdObj object containing sample names in @@snames slot
#' @param cluster_assignments Vector containing group assignments for each sample
#' @param dest Path to existing directory where output should be saved
#' @param pt_sz Size of points representing samples in plot (scaling factor)
#' @param n_dim Number of dimensions for embedding (either 2 or 3)
#' @param pt_label Vector of sample names corresponding to each point (same order as samples in \code{my_distmat} and \code{cluster_assignments})
#' @param cmap Vector containing colors by which points should be colored (corresponding to cluster_assignments)
#' @param w Width of plot in inches
#' @param h Height of plot in inches
#' @param scale.y Scaling factor for diffusion map y-axis
#' @param angle Rotation factor for diffusion map plot
#' @param ... Additional parameters to be passed to \code{DiffusionMap} function
#' @return DiffusionMap object containing biological sample embedding and associated metadata
#' @examples
#' 
#' my_phemdObj <- create_dataobj(all_expn_data, all_genes, as.character(snames))
#' my_phemdObj_lg <- remove_tiny_samples(my_phemdObj, 10)
#' my_phemdObj_lg <- aggregate_samples(my_phemdObj_lg, max_cells=1000)
#' my_phemdObj_monocle <- embed_cells(my_phemdObj_lg, data_model = 'gaussianff', sigma=0.02, maxIter=2)
#' my_phemdObj_monocle <- order_cells(my_phemdObj_monocle)
#' my_phemdObj_final <- cluster_individual_samples(my_phemdObj_monocle)
#' my_phemdObj_final <- generate_gdm(my_phemdObj_final)
#' my_EMD_mat <- compare_samples(my_phemdObj_final)
#' cluster_assignments <- group_samples(my_EMD_mat, distfun = 'hclust', ncluster=4)
#' print_cluster_assignments(cluster_assignments, my_phemdObj_final, '.')
#' dm <- plot_grouped_samples_dmap(my_EMD_mat, cluster_assignments, '.', pt_sz=2, pt_label = my_phemdObj_final@@snames)
#' 
plot_grouped_samples_dmap <- function(my_distmat, cluster_assignments, dest, pt_sz=1, n_dim=3, pt_label = NULL, cmap = NULL, w=8, h=5, scale.y=1, angle=40, ...) {
  extra_args <- list(...)
  #set.seed(15) #doesn't help because DiffusionMap function still produces non-reproducible results
  if(nrow(my_distmat) != ncol(my_distmat)) {
    stop('Error: my_distmat must be a square distance matrix')
  }
  if(nrow(my_distmat) != length(cluster_assignments)) {
    stop('Error: cluster_assignments must be the same length as the number of rows in my_distmat')
  }
  if(is.null(cmap)) {
    #cur_palette <- brewer.pal(max(cluster_assignments),"Set3")
    #if(length(cur_palette) > 1) cur_palette[2] <- "#FFD92F" #darker yellow
    #palette(cur_palette)

    getPalette <- colorRampPalette(brewer.pal(11, "Spectral"))
    cmap <- getPalette(max(cluster_assignments))
    if("#FFFFBF" %in% cmap) cmap[which(cmap == "#FFFFBF")] <- "#D3D3D3" #replace light yellow with grey
  }
  palette(cmap)


  # Plot inhibitor groups using diffusion map
  covars <- data.frame(covar1 = 1:nrow(my_distmat))
  if(nrow(my_distmat) < 30) {
    #dm <- DiffusionMap(covars, distance = as.dist(my_distmat), rotate=TRUE, n_local=3) # define preset n_local for small datasets e.g. melanoma dataset
    extra_args['n_local'] <- 3
  } else {
    #dm <- DiffusionMap(covars, distance = as.dist(my_distmat), rotate=TRUE)
  }
  #plot(dm, 1:3, pch=20, col=cluster_assignments, cex.symbols = pt_sz, interactive=interactive)

  dm_args <- c(list(data=covars, distance = as.dist(my_distmat)),
               extra_args[names(extra_args) %in% c("n_local", "density_norm", "rotate", "k", "sigma", "verbose")])
  dm <- do.call(DiffusionMap, dm_args)

  save(dm, file=paste(dest,'dm.RData',sep=""))

  # save embedding as png
  load(paste(dest,'dm.RData',sep=""))
  png(filename=paste(dest, "dm_inhib_embedding.png", sep=""),
      units="in",
      width=w,
      height=h,
      res=300)
  par(mar=c(2,2,1,3))
  #cluster_assignments_named <- sapply(cluster_assignments, function(x) paste("G-", x, sep=""))
  cluster_assignments_named <- sapply(cluster_assignments, function(x) intToUtf8(64+x))
  if(n_dim >= 3) {

    plot(dm, 1:3, pch=20, col=factor(cluster_assignments_named), pal=cmap, cex.symbols = pt_sz, box=FALSE, xlab="", ylab="", zlab="", y.margin.add = -0.5, draw_legend=TRUE, legend_opts = list(posx = c(0.85,0.88), posy = c(0.05, 0.7)), scale.y=scale.y, angle=angle)

  } else {
    plot(eigenvectors(dm)[,1], eigenvectors(dm)[,2], main = '', xlab = '', ylab = '', xaxt = 'n', yaxt = 'n', pch=20, col=factor(cluster_assignments_named), cex = pt_sz)
    #plot(dm, 1:2, pch=20, col=factor(cluster_assignments_named), pal=palette(), cex = pt_sz)
  }
  dev.off()

  if(!is.null(pt_label)) {
    png(filename=paste(dest, "dm_inhib_embedding_labeled.png", sep=""),
        units="px",
        width=8000,
        height=6000,
        res=300)
    cluster_assignments_named <- sapply(cluster_assignments, function(x) paste("G-", x, sep=""))
    if(n_dim >= 3) {
      s3d <- scatterplot3d(eigenvectors(dm)[,1], eigenvectors(dm)[,2], eigenvectors(dm)[,3], color=as.numeric(factor(cluster_assignments_named)), pch=20)
      s3d.coords <- s3d$xyz.convert(eigenvectors(dm)[,1], eigenvectors(dm)[,2], eigenvectors(dm)[,3])
      text(s3d.coords$x, s3d.coords$y,             # x and y coordinates
           labels=pt_label,               # text to plot
           cex=.3, pos=2)
    } else {
      plot(eigenvectors(dm)[,1], eigenvectors(dm)[,2], main = '', xlab = '', ylab = '', xaxt = 'n', yaxt = 'n', pch=20, col=factor(cluster_assignments_named), cex = pt_sz)
      if(!is.null(pt_label)) text(eigenvectors(dm)[,1:2],labels = pt_label, pos = 2, cex=0.4)
    }
    dev.off()
  }
  return(dm)
}

#' @title Plot cell subtype frequency histograms for each sample
#' @description Plots relative frequency ("weights") of cell subtypes ("bins" or "signatures") in each single-cell sample
#' @details \code{group_samples} must be called before calling this function. Saves plots in directory called "individual_inhibs"
#' @param myobj phemdObj object containing cell subtype relative frequency in @@data_cluster_weights slot
#' @param cluster_assignments Vector containing group assignments for each sample in myobj
#' @param dest Path to existing directory where output should be saved
#' @param cell_model Method by which cell state was modeled (either "monocle2" or "seurat")
#' @param cmap Vector containing colors by which histogram bars should be colored (optional)
#' @return None
#' @examples
#' 
#' my_phemdObj <- create_dataobj(all_expn_data, all_genes, as.character(snames))
#' my_phemdObj_lg <- remove_tiny_samples(my_phemdObj, 10)
#' my_phemdObj_lg <- aggregate_samples(my_phemdObj_lg, max_cells=1000)
#' my_phemdObj_monocle <- embed_cells(my_phemdObj_lg, data_model = 'gaussianff', sigma=0.02, maxIter=2)
#' my_phemdObj_monocle <- order_cells(my_phemdObj_monocle)
#' my_phemdObj_final <- cluster_individual_samples(my_phemdObj_monocle)
#' my_phemdObj_final <- generate_gdm(my_phemdObj_final)
#' my_EMD_mat <- compare_samples(my_phemdObj_final)
#' cluster_assignments <- group_samples(my_EMD_mat, distfun = 'hclust', ncluster=4)
#' print_cluster_assignments(cluster_assignments, my_phemdObj_final, '.')
#' dm <- plot_grouped_samples_dmap(my_EMD_mat, cluster_assignments, '.', pt_sz=2, pt_label = my_phemdObj_final@@snames)
#' plot_sample_histograms(my_phemdObj_final, cluster_assignments, '.')
#' 
plot_sample_histograms <- function(myobj, cluster_assignments, dest, cell_model='monocle2', cmap=NULL) {
  if(substr(dest,nchar(dest), nchar(dest)) != '/') dest <- paste(dest, '/', sep='') #ensure path ends with a slash
  unlink(paste(dest, 'individual_inhibs', sep=''), recursive=TRUE)
  dir.create(file.path(paste(dest, 'individual_inhibs', sep='')), showWarnings = FALSE) # create folder for output

  if(cell_model == 'monocle2') {
    monocle_obj <- myobj@monocle_obj
    labels <- pData(phenoData(monocle_obj))
    state_labels <- as.numeric(labels$State)
  } else if(cell_model == 'seurat') {
    seurat_obj <- myobj@seurat_obj
    state_labels <- as.numeric(seurat_obj@ident)
  } else {
    stop('Error: cell_model must either be "monocle2" or "seurat"')
  }


  cluster_weights <- myobj@data_cluster_weights
  #cmap <- rainbow(max(state_labels))
  if(is.null(cmap)) {
    getPalette <- colorRampPalette(brewer.pal(11, "Spectral"))
    cmap <- getPalette(max(state_labels))
  }

  snames <- myobj@snames

  for(i in 1:max(cluster_assignments)) {
    # Create folder "Group %s" if it doesn't already exist
    dir.create(file.path(paste(dest, sprintf('individual_inhibs/Group %s',intToUtf8(64+i)), sep='')), showWarnings = FALSE)
    cur_inhibs <- which(cluster_assignments == i)
    for(j in 1:length(cur_inhibs)) {
      cur_idx <- cur_inhibs[j]

      png(filename=paste(dest, sprintf("individual_inhibs/Group %s/%s.png",intToUtf8(64+i),snames[cur_idx]), sep=""),
          units="px",
          width=2400,
          height=1800,
          res=300)
      par(mar=c(6,6,4,2))
      barplot(cluster_weights[cur_idx,], main='', col=cmap, xlab='', ylab = "Frequency (%)", ylim = c(0, max(max(cluster_weights[cur_idx,])+0.1, 0.4)), cex.axis=1.5, cex.names = 2, cex.lab = 2.5, names.arg = 1:ncol(cluster_weights))
      title(xlab="Cell subtype", line=3.5, cex.lab=2.5)
      title(main=snames[cur_idx], line=0, cex.main=3)
      dev.off()
    }
  }
}

#' @title Plot cell subtype frequency histograms summarizing each group of samples
#' @description Plots relative frequency ("weights") of cell subtypes ("bins" or "signatures") summarizing each group of single-cell samples. Each summary histogram is computed by taking the bin-wise mean of all samples in the group
#' @details \code{group_samples} must be called before calling this function. Saves plots in directory called "summary_inhibs"
#' @param myobj phemdObj object containing cell subtype relative frequency in @@data_cluster_weights slot
#' @param cluster_assignments Vector containing group assignments for each sample in myobj
#' @param dest Path to existing directory where output should be saved
#' @param cell_model Method by which cell state was modeled (either "monocle2" or "seurat")
#' @param cmap Vector containing colors by which histogram bars should be colored (optional)
#' @return None
#' @examples
#' 
#' my_phemdObj <- create_dataobj(all_expn_data, all_genes, as.character(snames))
#' my_phemdObj_lg <- remove_tiny_samples(my_phemdObj, 10)
#' my_phemdObj_lg <- aggregate_samples(my_phemdObj_lg, max_cells=1000)
#' my_phemdObj_monocle <- embed_cells(my_phemdObj_lg, data_model = 'gaussianff', sigma=0.02, maxIter=2)
#' my_phemdObj_monocle <- order_cells(my_phemdObj_monocle)
#' my_phemdObj_final <- cluster_individual_samples(my_phemdObj_monocle)
#' my_phemdObj_final <- generate_gdm(my_phemdObj_final)
#' my_EMD_mat <- compare_samples(my_phemdObj_final)
#' cluster_assignments <- group_samples(my_EMD_mat, distfun = 'hclust', ncluster=4)
#' print_cluster_assignments(cluster_assignments, my_phemdObj_final, '.')
#' dm <- plot_grouped_samples_dmap(my_EMD_mat, cluster_assignments, '.', pt_sz=2, pt_label = my_phemdObj_final@@snames)
#' plot_summary_histograms(my_phemdObj_final, cluster_assignments, '.')
#' 
plot_summary_histograms <- function(myobj, cluster_assignments, dest, cell_model='monocle2', cmap=NULL) {
  if(substr(dest,nchar(dest), nchar(dest)) != '/') dest <- paste(dest, '/', sep='') #ensure path ends with a slash
  if(cell_model == 'monocle2') {
    monocle_obj <- myobj@monocle_obj
    labels <- pData(phenoData(monocle_obj))
    state_labels <- as.numeric(labels$State)

  } else if(cell_model == 'seurat') {
    seurat_obj <- myobj@seurat_obj
    state_labels <- as.numeric(seurat_obj@ident)

  } else {
    stop('Error: cell_model must be either "monocle2" or "seurat"')
  }


  cluster_weights <- myobj@data_cluster_weights

  if(is.null(cmap)) {
    getPalette <- colorRampPalette(brewer.pal(11, "Spectral"))
    cmap <- getPalette(max(state_labels))
  }

  unlink(paste(dest, 'summary_inhibs', sep=''), recursive=TRUE)
  dir.create(file.path(paste(dest, 'summary_inhibs', sep='')), showWarnings = FALSE)
  proto_inhibs <- matrix(0, max(cluster_assignments), ncol(cluster_weights))
  for(i in 1:max(cluster_assignments)) {
    if(sum(cluster_assignments == i) == 1) {
      proto_inhibs[i,] <- cluster_weights[which(cluster_assignments == i),]
    } else {
      proto_inhibs[i,] <- colMeans(cluster_weights[which(cluster_assignments == i),])
    }
  }

  for(i in 1:max(cluster_assignments)) {
    png(filename=paste(dest, sprintf("summary_inhibs/Group %s.png",intToUtf8(64+i)), sep=""),
        units="px",
        width=2400,
        height=1800,
        res=300)
    par(mar=c(6,6,4,2))
    if(max(proto_inhibs[i,]) > 0.4) ymax = max(proto_inhibs[i,])+0.1 else ymax = 0.4
    barplot(proto_inhibs[i,], col=cmap, main='', xlab='', ylab = "Frequency (%)", ylim = c(0, ymax), cex.axis=1.5, cex.names = 2, cex.lab = 2.5, names.arg = 1:ncol(proto_inhibs))

    title(xlab="Cell subtype", line=3.5, cex.lab=2.5)
    title(main=sprintf("Group %s", intToUtf8(64+i)), line=0, cex.main=3)

    dev.off()
  }
}

#' @title Plot cell yield of each sample as bar plot
#' @description Plots cell yield (number of viable cells) of each single-cell sample in decreasing order as horizontal bar plot
#' @param myobj phemdObj object containing expression data for each sample in 'data' slot
#' @param dest Path to existing directory where output should be saved
#' @param labels Vector containing group labels for samples (optional). If not provided, bars will be of uniform color (blue)
#' @param cmap Vector containing colors by which histogram bars should be colored (optional)
#' @param font_sz Scaling factor for font size of sample names in barplot
#' @param w Width of plot in inches
#' @param h Height of plot in inches
#' @return None
#' @examples
#' 
#' my_phemdObj <- create_dataobj(all_expn_data, all_genes, as.character(snames))
#' my_phemdObj_lg <- remove_tiny_samples(my_phemdObj, 10)
#' my_phemdObj_lg <- aggregate_samples(my_phemdObj_lg, max_cells=1000)
#' my_phemdObj_monocle <- embed_cells(my_phemdObj_lg, data_model = 'gaussianff', sigma=0.02, maxIter=2)
#' my_phemdObj_monocle <- order_cells(my_phemdObj_monocle)
#' my_phemdObj_final <- cluster_individual_samples(my_phemdObj_monocle)
#' my_phemdObj_final <- generate_gdm(my_phemdObj_final)
#' my_EMD_mat <- compare_samples(my_phemdObj_final)
#' cluster_assignments <- group_samples(my_EMD_mat, distfun = 'hclust', ncluster=4)
#' plot_cell_yield(my_phemdObj_final, '.', cluster_assignments, font_sz = 0.8)
#' 
plot_cell_yield <- function(myobj, dest, labels=NULL, cmap=NULL, font_sz = 0.6, w=8, h=9.5) {
  if(substr(dest,nchar(dest), nchar(dest)) != '/') dest <- paste(dest, '/', sep='') #ensure path ends with a slash
  nsample <- length(myobj@data)
  cell_yield <- rep(0, nsample)
  for(i in 1:nsample) {
    cell_yield[i] <- nrow(myobj@data[[i]])
  }

  order_idx <- order(cell_yield, decreasing=FALSE)
  cell_yield_ordered <- cell_yield[order_idx]
  snames_ordered <- myobj@snames[order_idx]

  png(filename=paste(dest, "cell_yield_ordered.png", sep=""),
      units="in",
      width=w,
      height=h,
      res=300)

  par(mar=c(6,6,2,2))

  if(!is.null(labels)) {
    if(length(labels) != nsample) {
      stop('Error: length of "labels" vector must be equal to length of myobj@data')
    }
    labels_ordered <- labels[order_idx]
    if(is.null(cmap)) {
      getPalette <- colorRampPalette(brewer.pal(11, "Spectral"))
      cmap <- getPalette(max(labels))
      if("#FFFFBF" %in% cmap) cmap[which(cmap == "#FFFFBF")] <- "#D3D3D3" #replace light yellow with grey
    }
    color_vec <- cmap[labels_ordered]
    xx <- barplot(cell_yield_ordered, main='', horiz=TRUE, names.arg=snames_ordered, las=1, cex.names=font_sz, col=color_vec)
  }  else {
    xx <- barplot(cell_yield_ordered, main='', horiz=TRUE, names.arg=snames_ordered, las=1, cex.names=font_sz, col='blue')
  }
  title(xlab="Cell yield (number of cells)", line=3, cex.lab=1.5)

  dev.off()
}



#' @title Prints cell yield of each sample as a table
#' @description Prints cell yield (number of viable cells) of each single-cell sample in decreasing order
#' @param myobj phemdObj object containing expression data for each sample in 'data' slot
#' @param dest Path to existing directory where output should be saved
#' @param cluster_assignments Vector of cluster assignments to be included as additional column in output table (optional)
#' @return None
#' @examples
#' 
#' my_phemdObj <- create_dataobj(all_expn_data, all_genes, as.character(snames))
#' my_phemdObj_lg <- remove_tiny_samples(my_phemdObj, 10)
#' my_phemdObj_lg <- aggregate_samples(my_phemdObj_lg, max_cells=1000)
#' my_phemdObj_monocle <- embed_cells(my_phemdObj_lg, data_model = 'gaussianff', sigma=0.02, maxIter=2)
#' my_phemdObj_monocle <- order_cells(my_phemdObj_monocle)
#' my_phemdObj_final <- cluster_individual_samples(my_phemdObj_monocle)
#' my_phemdObj_final <- generate_gdm(my_phemdObj_final)
#' my_EMD_mat <- compare_samples(my_phemdObj_final)
#' cluster_assignments <- group_samples(my_EMD_mat, distfun = 'hclust', ncluster=4)
#' print_cell_yield(my_phemdObj_final, '.', cluster_assignments)
#' 
print_cell_yield <- function(myobj, dest, cluster_assignments=NULL) {
  if(substr(dest,nchar(dest), nchar(dest)) != '/') dest <- paste(dest, '/', sep='') #ensure path ends with a slash
  nsample <- length(myobj@data)
  cell_yield <- rep(0, nsample)
  for(i in 1:nsample) {
    cell_yield[i] <- nrow(myobj@data[[i]])
  }

  order_idx <- order(cell_yield, decreasing=FALSE)
  cell_yield_ordered <- cell_yield[order_idx]
  snames_ordered <- myobj@snames[order_idx]
  cell_yield_tab <- cbind.data.frame(snames_ordered, cell_yield_ordered)
  colnames(cell_yield_tab) <- c('sample_ID', 'cell_yield')
  if(!is.null(cluster_assignments)) {
    cell_yield_tab$cluster_ID <- cluster_assignments[order_idx]
  }
  
  write.table(cell_yield_tab, file=paste(dest, 'cell_yield_tab.txt', sep=''), sep='\t', quote=FALSE, row.names=FALSE, col.names=TRUE)
}
