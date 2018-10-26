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
#' my_phemdObj <- createDataObj(all_expn_data, all_genes, as.character(snames))
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
  snames <- sampleNames(obj)
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

#' @title Saves cell subtype frequency histograms for each sample
#' @description Saves relative frequency ("weights") of cell subtypes ("bins" or "signatures") in each single-cell sample
#' @details \code{groupSamples} must be called before calling this function. Saves plots in directory called "individual_inhibs"
#' @param myobj phemdObj object containing cell subtype relative frequency in @@data_cluster_weights slot
#' @param cluster_assignments Vector containing group assignments for each sample in myobj
#' @param dest Path to existing directory where output should be saved
#' @param cell_model Method by which cell state was modeled (either "monocle2" or "seurat")
#' @param cmap Vector containing colors by which histogram bars should be colored (optional)
#' @param overwrite Boolean representing whether or not to overwrite contents of "dest" with output of saveSampleHistograms
#' @return None
#' @examples
#' 
#' my_phemdObj <- createDataObj(all_expn_data, all_genes, as.character(snames))
#' my_phemdObj_lg <- removeTinySamples(my_phemdObj, 10)
#' my_phemdObj_lg <- aggregateSamples(my_phemdObj_lg, max_cells=1000)
#' my_phemdObj_monocle <- embedCells(my_phemdObj_lg, data_model = 'gaussianff', sigma=0.02, maxIter=2)
#' my_phemdObj_monocle <- orderCellsMonocle(my_phemdObj_monocle)
#' my_phemdObj_final <- clusterIndividualSamples(my_phemdObj_monocle)
#' my_phemdObj_final <- generateGDM(my_phemdObj_final)
#' my_EMD_mat <- compareSamples(my_phemdObj_final)
#' cluster_assignments <- groupSamples(my_EMD_mat, distfun = 'hclust', ncluster=4)
#' printClusterAssignments(cluster_assignments, my_phemdObj_final, '.', overwrite=TRUE)
#' dm <- plotGroupedSamplesDmap(my_EMD_mat, cluster_assignments, dest=NULL, pt_sz=2)
#' saveSampleHistograms(my_phemdObj_final, cluster_assignments, '.', overwrite=TRUE)
#' 
saveSampleHistograms <- function(myobj, cluster_assignments, dest, cell_model=c('monocle2', 'seurat'), cmap=NULL, overwrite=FALSE) {
  if(substr(dest,nchar(dest), nchar(dest)) != '/') dest <- paste(dest, '/', sep='') #ensure path ends with a slash
  if(dir.exists(paste(dest, 'individual_inhibs', sep='')) && overwrite==FALSE) {
    stop('Directory "individual_inhibs" already exists in specified path. Set "overwrite" parameter to TRUE if you want to overwrite existing directory')
  }
  unlink(paste(dest, 'individual_inhibs', sep=''), recursive=TRUE)
  dir.create(file.path(paste(dest, 'individual_inhibs', sep='')), showWarnings = FALSE) # create folder for output
  
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
  #cmap <- rainbow(max(state_labels))
  if(is.null(cmap)) {
    getPalette <- colorRampPalette(brewer.pal(11, "Spectral"))
    cmap <- getPalette(max(state_labels))
  }
  
  snames <- sampleNames(myobj)
  
  for(i in seq_len(max(cluster_assignments))) {
    # Create folder "Group %s" if it doesn't already exist
    dir.create(file.path(paste(dest, sprintf('individual_inhibs/Group %s',intToUtf8(64+i)), sep='')), showWarnings = FALSE)
    cur_inhibs <- which(cluster_assignments == i)
    for(j in seq_len(length(cur_inhibs))) {
      cur_idx <- cur_inhibs[j]
      
      png(filename=paste(dest, sprintf("individual_inhibs/Group %s/%s.png",intToUtf8(64+i),snames[cur_idx]), sep=""),
          units="px",
          width=2400,
          height=1800,
          res=300)
      par(mar=c(6,6,4,2))
      barplot(cluster_weights[cur_idx,], main='', col=cmap, xlab='', ylab = "Frequency (%)", ylim = c(0, max(max(cluster_weights[cur_idx,])+0.1, 0.4)), cex.axis=1.5, cex.names = 2, cex.lab = 2.5, names.arg = seq_len(ncol(cluster_weights)))
      title(xlab="Cell subtype", line=3.5, cex.lab=2.5)
      title(main=snames[cur_idx], line=0, cex.main=3)
      dev.off()
    }
  }
}


#' @title Prints cell yield of each sample as a table
#' @description Prints cell yield (number of viable cells) of each single-cell sample in decreasing order
#' @param myobj phemdObj object containing expression data for each sample in 'data' slot
#' @param dest Path to existing directory where output should be saved
#' @param cluster_assignments Vector of cluster assignments to be included as additional column in output table (optional)
#' @return None
#' @examples
#' 
#' my_phemdObj <- createDataObj(all_expn_data, all_genes, as.character(snames))
#' my_phemdObj_lg <- removeTinySamples(my_phemdObj, 10)
#' my_phemdObj_lg <- aggregateSamples(my_phemdObj_lg, max_cells=1000)
#' my_phemdObj_monocle <- embedCells(my_phemdObj_lg, data_model = 'gaussianff', sigma=0.02, maxIter=2)
#' my_phemdObj_monocle <- orderCellsMonocle(my_phemdObj_monocle)
#' my_phemdObj_final <- clusterIndividualSamples(my_phemdObj_monocle)
#' my_phemdObj_final <- generateGDM(my_phemdObj_final)
#' my_EMD_mat <- compareSamples(my_phemdObj_final)
#' cluster_assignments <- groupSamples(my_EMD_mat, distfun = 'hclust', ncluster=4)
#' printCellYield(my_phemdObj_final, '.', cluster_assignments)
#' 
printCellYield <- function(myobj, dest, cluster_assignments=NULL) {
  if(substr(dest,nchar(dest), nchar(dest)) != '/') dest <- paste(dest, '/', sep='') #ensure path ends with a slash
  nsample <- length(rawExpn(myobj))
  cell_yield <- vapply(rawExpn(myobj), nrow, integer(1L))
  
  order_idx <- order(cell_yield, decreasing=FALSE)
  cell_yield_ordered <- cell_yield[order_idx]
  snames_ordered <- sampleNames(myobj)[order_idx]
  cell_yield_tab <- cbind.data.frame(snames_ordered, cell_yield_ordered)
  colnames(cell_yield_tab) <- c('sample_ID', 'cell_yield')
  if(!is.null(cluster_assignments)) {
    cluster_assignments_reordered <- cluster_assignments[order_idx]
    cell_yield_tab$cluster_ID <- sapply(cluster_assignments_reordered, function(x) intToUtf8(64+x))
  }
  
  write.table(cell_yield_tab, file=paste(dest, 'cell_yield_tab.txt', sep=''), sep='\t', quote=FALSE, row.names=FALSE, col.names=TRUE)
}

#' @title Prints cell subtype distribution for each sample as a table
#' @description Prints cell subtype distribution for each single-cell sample along with (optional) final inhibitor cluster assignment
#' @param myobj phemdObj object containing expression data for each sample in 'data' slot
#' @param dest Path to existing directory where output should be saved
#' @param cluster_assignments Vector of cluster assignments to be included as additional column in output table (optional)
#' @return None
#' @examples
#' 
#' my_phemdObj <- createDataObj(all_expn_data, all_genes, as.character(snames))
#' my_phemdObj_lg <- removeTinySamples(my_phemdObj, 10)
#' my_phemdObj_lg <- aggregateSamples(my_phemdObj_lg, max_cells=1000)
#' my_phemdObj_monocle <- embedCells(my_phemdObj_lg, data_model = 'gaussianff', sigma=0.02, maxIter=2)
#' my_phemdObj_monocle <- orderCellsMonocle(my_phemdObj_monocle)
#' my_phemdObj_final <- clusterIndividualSamples(my_phemdObj_monocle)
#' my_phemdObj_final <- generateGDM(my_phemdObj_final)
#' my_EMD_mat <- compareSamples(my_phemdObj_final)
#' cluster_assignments <- groupSamples(my_EMD_mat, distfun = 'hclust', ncluster=4)
#' printSampleCelltypeFreqs(my_phemdObj_final, '.', cluster_assignments)
#' 
printSampleCelltypeFreqs <- function(myobj, dest, cluster_assignments=NULL) {
  if(substr(dest,nchar(dest), nchar(dest)) != '/') dest <- paste(dest, '/', sep='') #ensure path ends with a slash
  celltype_freqs <- as.data.frame(celltypeFreqs(myobj))
  celltype_freqs <- round(celltype_freqs, digits=3)
  colnames(celltype_freqs) <- paste('C-', seq_len(ncol(celltype_freqs)), sep='')
  
  if(!is.null(cluster_assignments)) {
    celltype_freqs$Sample.Cluster.ID <- sapply(cluster_assignments, function(x) intToUtf8(64+x))
  }
  
  outfile <- cbind.data.frame(Sample.Name=sampleNames(myobj), celltype_freqs)
  
  write.table(outfile, file=paste(dest, 'cell_subtype_distributions.txt', sep=''), sep='\t', quote=FALSE, row.names=FALSE, col.names=TRUE)
}
