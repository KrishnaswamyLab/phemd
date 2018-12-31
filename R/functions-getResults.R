#' @title Writes samples to file based on community detection group assignments
#' @description Takes vector of cluster assignments and phemdObj containing sample names and writes sample groups to file
#' @details Order of samples in obj@@snames is assumed to be the same as the order of group assignments in cluster_assignments
#' @param cluster_assignments Vector containing group assignments for each sample
#' @param obj phemdObj object containing sample names in @@snames slot
#' @param dest Path to existing directory where output should be saved
#' @param overwrite Boolean representing whether or not to overwrite contents of "dest" with output of printClusterAssignments
#' @return None
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
#' printClusterAssignments(cluster_assignments, my_phemdObj_final, '.', overwrite=TRUE)
#' 
printClusterAssignments <- function(cluster_assignments, obj, dest, overwrite=FALSE) {
  snames <- sNames(obj)
  if(dir.exists(paste(dest, 'sample_groups', sep='')) && overwrite==FALSE) {
    stop('Directory "sample_groups" already exists in specified path. Set "overwrite" parameter to TRUE if you want to overwrite existing directory')
  }
  unlink(paste(dest, 'sample_groups', sep=''), recursive=TRUE)
  dir.create(file.path(paste(dest, 'sample_groups', sep='')), showWarnings = FALSE) # create folder for output
  for(i in seq_len(max(cluster_assignments))) {
    cur_cluster_idx <- which(cluster_assignments == i)
    cur_file <- sprintf('sample_groups/scluster_%s.txt', intToUtf8(64+i))
    cur_file <- strcat(dest, cur_file)
    write(snames[cur_cluster_idx], file=cur_file, sep="\n")
  }
}

#' @title Gets cell subtype frequency histograms for each sample by cluster ID
#' @description Gets relative frequency ("weights") of cell subtypes ("bins" or "signatures") in each single-cell sample
#' @details \code{groupSamples} must be called before calling this function. Saves plots in directory called "individual_inhibs"
#' @param myobj phemdObj object containing cell subtype relative frequency in @@data_cluster_weights slot
#' @param cluster_assignments Vector containing group assignments for each sample in myobj
#' @param cell_model Method by which cell state was modeled (either "monocle2" or "seurat")
#' @return List of lists, with outer list representing sample cluster ID and inner list representing cell subtype frequencies of given sample
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
#' weights_by_cluster <- getSampleHistsByCluster(my_phemdObj_final, cluster_assignments)
#' 
getSampleHistsByCluster <- function(myobj, cluster_assignments, cell_model=c('monocle2', 'seurat')) {
  cell_model <- match.arg(cell_model, c('monocle2','seurat'))
  if(cell_model == 'monocle2') {
    monocle_obj <- monocleInfo(myobj)
    labels <- pData(phenoData(monocle_obj))
    state_labels <- as.numeric(labels$State)
  } else if(cell_model == 'seurat') {
    seurat_obj <- seuratInfo(myobj)
    state_labels <- as.numeric(GetIdent(seurat_obj, uniq=FALSE))
  } else {
    stop('Error: cell_model must either be "monocle2" or "seurat"')
  }
  
  cluster_weights <- celltypeFreqs(myobj)
  snames <- sNames(myobj)
  
  weights_by_cluster <- list()
  for(i in seq_len(max(cluster_assignments))) {
    cur_inhibs <- which(cluster_assignments == i)
    for(j in seq_len(length(cur_inhibs))) {
      cur_idx <- cur_inhibs[j]
      cur_sname <- snames[cur_idx]
      if(is.null(weights_by_cluster[[intToUtf8(64+i)]])) {
        weights_by_cluster[[intToUtf8(64+i)]] <- list()
      }
      weights_by_cluster[[intToUtf8(64+i)]][[cur_sname]] <- cluster_weights[cur_idx,]
    }
  }
  return(weights_by_cluster)
}


#' @title Gets cell yield of each sample as a table
#' @description Gets cell yield (number of viable cells) of each single-cell sample in decreasing order
#' @param myobj phemdObj object containing expression data for each sample in 'data' slot
#' @param cluster_assignments Vector of cluster assignments to be included as additional column in output table (optional)
#' @return Data frame representing cell yield of each sample
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
#' getCellYield(my_phemdObj_final, cluster_assignments)
#' 
getCellYield <- function(myobj, cluster_assignments=NULL) {
  nsample <- length(rawExpn(myobj))
  cell_yield <- vapply(rawExpn(myobj), nrow, integer(1L))
  
  order_idx <- order(cell_yield, decreasing=FALSE)
  cell_yield_ordered <- cell_yield[order_idx]
  snames_ordered <- sNames(myobj)[order_idx]
  cell_yield_tab <- cbind.data.frame(snames_ordered, cell_yield_ordered)
  colnames(cell_yield_tab) <- c('sample_ID', 'cell_yield')
  if(!is.null(cluster_assignments)) {
    cluster_assignments_reordered <- cluster_assignments[order_idx]
    cell_yield_tab$cluster_ID <- vapply(cluster_assignments_reordered, function(x) intToUtf8(64+x), '')
  }
  return(cell_yield_tab)
}

#' @title Returns cell subtype distribution for each sample as a table
#' @description Returns cell subtype distribution for each single-cell sample along with (optional) final inhibitor cluster assignment
#' @param myobj phemdObj object containing expression data for each sample in 'data' slot
#' @param cluster_assignments Vector of cluster assignments to be included as additional column in output table (optional)
#' @return Data frame representing relative frequencies of each cell subtype along with (optional) final inhibitor cluster assignment for each single-cell sample
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
#' getSampleCelltypeFreqs(my_phemdObj_final, cluster_assignments)
#' 
getSampleCelltypeFreqs <- function(myobj, cluster_assignments=NULL) {
  celltype_freqs <- as.data.frame(celltypeFreqs(myobj))
  celltype_freqs <- round(celltype_freqs, digits=3)
  colnames(celltype_freqs) <- paste('C-', seq_len(ncol(celltype_freqs)), sep='')
  
  if(!is.null(cluster_assignments)) {
    celltype_freqs$Sample.Cluster.ID <- vapply(cluster_assignments, function(x) intToUtf8(64+x), '')
  }
  
  weightsTab <- cbind.data.frame(Sample.Name=sNames(myobj), celltype_freqs)
  return(weightsTab)
}
