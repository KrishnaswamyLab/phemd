#' @title Plots Monocle2 cell embedding plots
#' @description Takes as input a Phemd object containing either a Monocle2 object or Seurat object (already embedded and ordered) and plots cell embedding plots side by side. Optionally saves to specified folder.
#' @details \code{embedCells} and \code{orderCellsMonocle} need to be called before calling this function. Required additional packages: 'RColorBrewer', 'cowplot'
#' @param obj 'Phemd' object containing Monocle 2 object
#' @param cell_model Method by which cell state was modeled (either "monocle2" or "seurat")
#' @param cmap User-specified colormap to use to color cell state embedding (optional)
#' @param w Width of plot in inches
#' @param h Height of plot in inches
#' @param pt_sz Scalar factor for point size
#' @param ndims Number of dimensions to use for dimensionality reduction in case it hasn't been performed yet (only relevant when using Seurat data as input)
#' @return Colormap (vector of colors) used to color Monocle2 cell state embedding
#' @examples
#' my_phemdObj <- createDataObj(all_expn_data, all_genes, as.character(snames_data))
#' my_phemdObj_lg <- removeTinySamples(my_phemdObj, 10)
#' my_phemdObj_lg <- aggregateSamples(my_phemdObj_lg, max_cells=1000)
#' my_phemdObj_monocle <- embedCells(my_phemdObj_lg, data_model='gaussianff', sigma=0.02, maxIter=2)
#' my_phemdObj_monocle <- orderCellsMonocle(my_phemdObj_monocle)
#' cmap <- plotEmbeddings(my_phemdObj_monocle)
plotEmbeddings <- function(obj, cell_model=c('monocle2', 'seurat'), cmap=NULL, w=4, h=5, pt_sz=1, ndims=NULL) {
  cell_model <- match.arg(cell_model, c('monocle2','seurat'))
  if(cell_model == 'monocle2') {
    monocle_obj <- monocleInfo(obj)
    cell_embedding <- reducedDimS(monocle_obj)
    mydata <- pooledCells(obj)
    
    # Extract state labels from monocle data object
    labels <- pData(phenoData(monocle_obj))
    state_labels <- as.numeric(labels$State)
    
    levels <- levels(factor(state_labels))
    levels_renamed <- vapply(levels, function(x) paste("C-", x, sep=""), "")
    
    if(is.null(cmap)) {
      getPalette <- colorRampPalette(brewer.pal(11, "Spectral"))
      cmap <- getPalette(max(state_labels))
      if("#FFFFBF" %in% cmap) cmap[which(cmap == "#FFFFBF")] <- "#D3D3D3" #replace light yellow with grey
      cmap <- sample(cmap)
    }
    palette(cmap)
    
    cell_embedding_t <- as.data.frame(t(cell_embedding))
    # visualize traj colored by state
    myplot_state <- ggplot(cell_embedding_t, aes(x=cell_embedding_t[,1], y=cell_embedding_t[,2], color=factor(state_labels))) +
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
    
    # visualize traj colored by pseudotime
    ncolor <- 9
    palette(brewer.pal(ncolor, "Blues"))
    
    col.labels <- labels$Pseudotime
    
    # visualize traj colored by pseudotime
    myplot_pt <- ggplot(cell_embedding_t, aes(x=cell_embedding_t[,1], y=cell_embedding_t[,2], color=col.labels)) +
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
    
    print(plot_grid(myplot_state, myplot_pt, ncol=2))
    return(cmap)
  } else if(cell_model == 'seurat') {
    seurat_obj <- seuratInfo(obj)
    
    if(!'tsne' %in% names(seurat_obj@dr)) {
      print('Running t-SNE...')
      if(is.null(ndims)) ndims <- 10
      seurat_obj <- RunTSNE(seurat_obj, reduction.use="cca.aligned", dims.use=seq_len(ndims))
    }
    
    # define color map
    if(is.null(cmap)) {
      getPalette <- colorRampPalette(brewer.pal(11, "Spectral"))
      cmap <- getPalette(max(as.numeric(seurat_obj@ident)))
      if("#FFFFBF" %in% cmap) cmap[which(cmap == "#FFFFBF")] <- "#D3D3D3" #replace light yellow with grey
      cmap <- sample(cmap)
    }
    
    TSNEPlot(seurat_obj, do.label=FALSE, pt.size=pt_sz, colors.use=cmap)
    
  } else {
    stop('Error: cell_model must be either "monocle2" or "seurat"')
  }
  return(cmap)
}

#' @title Plot heatmap of cell subtypes
#' @description Takes as input a Phemd object containing either a Monocle2 or Seurat object (already embedded and ordered) and saves heatmap describing cell subtypes to specified folder
#' @details \code{embedCells} and \code{orderCellsMonocle} need to be called before calling this function. Required additional package: 'pheatmap'
#' @param obj 'Phemd' object containing Monocle 2 object
#' @param cell_model Method by which cell state was modeled (either "monocle2" or "seurat")
#' @param selected_genes Vector containing gene names to include in heatmap (optional)
#' @param w Width of plot in inches
#' @param h Height of plot in inches
#' @param ... Additional parameters to be passed on to pheatmap function
#' @return Heatmap containing expression values for each cell subtype. If cell_model is 'seurat', then returns a list of heatmaps (1 for each batch) that may be subsequently plotted individually
#' @examples
#' 
#' my_phemdObj <- createDataObj(all_expn_data, all_genes, as.character(snames_data))
#' my_phemdObj_lg <- removeTinySamples(my_phemdObj, 10)
#' my_phemdObj_lg <- aggregateSamples(my_phemdObj_lg, max_cells=1000)
#' my_phemdObj_lg <- selectFeatures(my_phemdObj_lg, selected_genes)
#' my_phemdObj_monocle <- embedCells(my_phemdObj_lg, data_model = 'gaussianff', 
#' pseudo_expr=0, sigma=0.02, maxIter=2)
#' my_phemdObj_monocle <- orderCellsMonocle(my_phemdObj_monocle)
#' myheatmap <- plotHeatmaps(my_phemdObj_monocle, cell_model='monocle2')
#' 
plotHeatmaps <- function(obj, cell_model=c('monocle2','seurat'), selected_genes=NULL, w=8, h=5, ...) {
  cell_model <- match.arg(cell_model, c('monocle2','seurat'))
  if(cell_model == 'monocle2') {
    # retrieve reference clusters
    ref_clusters <- retrieveRefClusters(obj, cell_model='monocle2')
    selected_clusters <- seq_len(length(ref_clusters))
    myheatmap <- matrix(0, nrow=length(selected_clusters), ncol=ncol(ref_clusters[[1]]))
    for(i in selected_clusters) {
      cur_cluster <- ref_clusters[[i]]
      if(!is.null(cur_cluster)) { #at least 1 cell
        if(nrow(cur_cluster) > 1) {
          myheatmap[i,] <- colMeans(cur_cluster)
        } else {
          myheatmap[i,] <- cur_cluster #only 1 cell
        }
      } 
    }
    
    selected_clusters_renamed <- vapply(names(ref_clusters), function(x) paste("C-", x, sep=""), "")
    
    rownames(myheatmap) <- selected_clusters_renamed
    colnames(myheatmap) <- selectMarkers(obj)
    
    if(!is.null(selected_genes)) {
      col_tokeep <- match(selected_genes, selectMarkers(obj))
      if(sum(is.na(col_tokeep)) > 0) {
        genes_not_found <- ''
        missing_idx <- which(is.na(col_tokeep))
        for(i in seq_len(length(missing_idx))) {
          if(i == 1) genes_not_found <- paste(genes_not_found, selected_genes[missing_idx[i]], sep='')
          else genes_not_found <- paste(genes_not_found, selected_genes[missing_idx[i]], sep=', ')
        }
        print(sprintf("Genes not found: %s", genes_not_found, sep=""))
      }
      col_tokeep <- col_tokeep[!is.na(col_tokeep)]
      myheatmap <- myheatmap[,col_tokeep]
    }
    
    myheatmap[is.nan(myheatmap)] <- 0 #this in the event of empty clusters
    
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
             width=w,
             height=h,
             ...
    )
    return(myheatmap2)
  } else if(cell_model == 'seurat') {
    seurat_obj <- seuratInfo(obj)
    state_labels <- as.numeric(as.character(GetIdent(seurat_obj, uniq=FALSE)))
    names(state_labels) <- rownames(seurat_obj@meta.data)
    ref_data <- t(as.matrix(GetAssayData(seurat_obj, assay.type='RNA', slot='raw.data')))
    
    batches <- unique(batchIDs(obj))
    myheatmaps_all <- list()
    for(batch_id in batches) {
      cell_idx_curplt <- which(seurat_obj@meta.data$plt == batch_id)
      if(length(cell_idx_curplt) == 0) {
        stop(sprintf('Error: no cells in reference set match the experiment_id %s. Please check batchIDs(phemdObj).', batch_id))
      }
      cur_ref_data <- ref_data[cell_idx_curplt,]
      cur_state_labels <- state_labels[cell_idx_curplt]
      
      myheatmap <- matrix(0, nrow=max(state_labels), ncol=ncol(cur_ref_data))
      for(i in seq_len(max(state_labels))) {
        cur_idx <- which(cur_state_labels == i)
        cur_cluster <- cur_ref_data[cur_idx,]
        if(length(cur_idx) > 1) myheatmap[i,] <- colMeans(cur_cluster)
        if(length(cur_idx) == 1) myheatmap[i,] <- cur_cluster
      }
      
      selected_clusters_renamed <- vapply(seq_len(max(state_labels)), function(x) paste("C-", x, sep=""), "")
      
      rownames(myheatmap) <- selected_clusters_renamed
      colnames(myheatmap) <- selectMarkers(obj)
      
      if(!is.null(selected_genes)) {
        col_tokeep <- match(selected_genes, selectMarkers(obj))
        if(sum(is.na(col_tokeep)) > 0) {
          genes_not_found <- ''
          missing_idx <- which(is.na(col_tokeep))
          for(i in seq_len(length(missing_idx))) {
            if(i == 1) genes_not_found <- paste(genes_not_found, selected_genes[missing_idx[i]], sep='')
            else genes_not_found <- paste(genes_not_found, selected_genes[missing_idx[i]], sep=', ')
          }
          print(sprintf("Genes not found: %s", genes_not_found, sep=""))
        }
        col_tokeep <- col_tokeep[!is.na(col_tokeep)]
        myheatmap <- myheatmap[,col_tokeep]
      }
      
      myheatmap[is.nan(myheatmap)] <- 0 #this in the event of empty clusters
      
      myheatmap2 <- log2(myheatmap - min(myheatmap) + 1)
      myheatmaps_all[[batch_id]] <- myheatmap2
    }
    
    for(i in seq_len(length(myheatmaps_all))) {
      if(!exists('myheatmaps_avg')) myheatmaps_avg <- myheatmaps_all[[i]]
      else myheatmaps_avg <- myheatmaps_avg + myheatmaps_all[[i]]
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
             width=w,
             height=h,
             ...)
    return(myheatmaps_all)
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
drawColnames45 <- function(coln, gaps, ...) {
  coord <- pheatmap:::find_coordinates(length(coln), gaps)
  x     <- coord$coord - 0.5 * coord$size
  res   <- grid::textGrob(
    coln, x = x, y = unit(1, "npc") - unit(3,"bigpts"),
    vjust = 0.75, hjust = 1, rot = 45, gp = grid::gpar(...)
  )
  return(res)
}


#' @title Plot diffusion map embedding of samples based on distance matrix
#' @description Visualizes diffusion map for network of samples based on square distance matrix (sample-sample pairwise dissimilarity)
#' @details Requires 'destiny' package
#' @param my_distmat phemdObj object containing sample names in @@snames slot
#' @param cluster_assignments Vector containing group assignments for each sample
#' @param pt_sz Size of points representing samples in plot (scaling factor)
#' @param n_dim Number of dimensions for embedding (either 2 or 3)
#' @param pt_label Vector of sample names corresponding to each point (same order as samples in \code{my_distmat} and \code{cluster_assignments})
#' @param cmap Vector containing colors by which points should be colored (corresponding to cluster_assignments)
#' @param w Width of plot in inches
#' @param h Height of plot in inches
#' @param scale.y Scaling factor for diffusion map y-axis
#' @param angle Rotation factor for diffusion map plot
#' @param autosave Boolean denoting whether or not to save output diffusion map
#' @param ... Additional parameters to be passed to \code{DiffusionMap} function
#' @return DiffusionMap object containing biological sample embedding and associated metadata
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
#' dm <- plotGroupedSamplesDmap(my_EMD_mat, cluster_assignments, pt_sz=2)
#' 
plotGroupedSamplesDmap <- function(my_distmat, cluster_assignments, pt_sz=1, n_dim=3, pt_label = NULL, cmap = NULL, w=8, h=5, scale.y=1, angle=40, autosave=FALSE, ...) {
  extra_args <- list(...)
  if(nrow(my_distmat) != ncol(my_distmat)) {
    stop('Error: my_distmat must be a square distance matrix')
  }
  if(nrow(my_distmat) != length(cluster_assignments)) {
    stop('Error: cluster_assignments must be the same length as the number of rows in my_distmat')
  }
  if(is.null(cmap)) {
    getPalette <- colorRampPalette(brewer.pal(11, "Spectral"))
    cmap <- getPalette(max(c(cluster_assignments),3)) # min palette size = 3
    if("#FFFFBF" %in% cmap) cmap[which(cmap == "#FFFFBF")] <- "#D3D3D3" #replace light yellow with grey
  }
  if(length(cmap) > 1) palette(cmap)
  
  # Plot inhibitor groups using diffusion map
  covars <- data.frame(covar1 = seq_len(nrow(my_distmat)))
  if(nrow(my_distmat) < 30) {
    extra_args['n_local'] <- 3
  }
  
  dm_args <- c(list(data=covars, distance = as.dist(my_distmat)),
               extra_args[names(extra_args) %in% c("n_local", "density_norm", "rotate", "k", "sigma", "verbose")])
  dm <- do.call(DiffusionMap, dm_args)
  
  
  cluster_assignments_named <- vapply(cluster_assignments, function(x) intToUtf8(64+x), "")
  if(n_dim >= 3) {
    
    plot(dm, c(1,2,3), pch=20, col=factor(cluster_assignments_named), pal=cmap, cex.symbols = pt_sz, box=FALSE, xlab="", ylab="", zlab="", y.margin.add = -0.5, draw_legend=TRUE, legend_opts = list(posx = c(0.85,0.88), posy = c(0.05, 0.7)), scale.y=scale.y, angle=angle)
    
  } else {
    plot(eigenvectors(dm)[,1], eigenvectors(dm)[,2], main = '', xlab = '', ylab = '', xaxt = 'n', yaxt = 'n', pch=20, col=factor(cluster_assignments_named), cex = pt_sz)
  }
  
  if(!is.null(pt_label)) {
    cluster_assignments_named <- vapply(cluster_assignments, function(x) paste("G-", x, sep=""), "")
    if(n_dim >= 3) {
      s3d <- scatterplot3d(eigenvectors(dm)[,1], eigenvectors(dm)[,2], eigenvectors(dm)[,3], color=as.numeric(factor(cluster_assignments_named)), pch=20)
      s3d.coords <- s3d$xyz.convert(eigenvectors(dm)[,1], eigenvectors(dm)[,2], eigenvectors(dm)[,3])
      text(s3d.coords$x, s3d.coords$y,             # x and y coordinates
           labels=pt_label,               # text to plot
           cex=.3, pos=2)
    } else {
      plot(eigenvectors(dm)[,1], eigenvectors(dm)[,2], main = '', xlab = '', ylab = '', xaxt = 'n', yaxt = 'n', pch=20, col=factor(cluster_assignments_named), cex = pt_sz)
      text(eigenvectors(dm)[,c(1,2)],labels = pt_label, pos = 2, cex=0.4)
    }
  }
  return(dm)
}

#' @title Plots cell subtype frequency histograms summarizing each group of samples
#' @description Visualizes plots of relative frequency ("weights") of cell subtypes ("bins" or "signatures") summarizing each group of single-cell samples. Each summary histogram is computed by taking the bin-wise mean of all samples in the group
#' @details \code{groupSamples} must be called before calling this function. Saves plots in directory called "summary_inhibs"
#' @param myobj Phemd object containing cell subtype relative frequency in @@data_cluster_weights slot
#' @param cluster_assignments Vector containing group assignments for each sample in myobj
#' @param cell_model Method by which cell state was modeled (either "monocle2" or "seurat")
#' @param cmap Vector containing colors by which histogram bars should be colored (optional)
#' @param ncol.plot Number of columns to use to plot multi-panel histogram plot
#' @param ax.lab.sz Scaling factor for axis labels (default 2.5)
#' @param title.sz Scaling factor for plot title (default 3)
#' @return None
#' @examples
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
#' dm <- plotGroupedSamplesDmap(my_EMD_mat, cluster_assignments, '.', pt_sz=2, pt_label = sNames(my_phemdObj_final))
#' plotSummaryHistograms(my_phemdObj_final, cluster_assignments, cell_model='monocle2)
#' 
plotSummaryHistograms <- function(myobj, cluster_assignments, cell_model=c('monocle2','seurat'), cmap=NULL, ncol.plot=4, ax.lab.sz=2.5, title.sz=3) {
  cell_model <- match.arg(cell_model, c('monocle2','seurat'))
  if(cell_model == 'monocle2') {
    monocle_obj <- monocleInfo(myobj)
    labels <- pData(phenoData(monocle_obj))
    state_labels <- as.numeric(labels$State)
    
  } else if(cell_model == 'seurat') {
    seurat_obj <- seuratInfo(myobj)
    state_labels <- as.numeric(GetIdent(seurat_obj, uniq=FALSE))
    
  } else {
    stop('Error: cell_model must be either "monocle2" or "seurat"')
  }
  
  
  cluster_weights <- celltypeFreqs(myobj)
  
  if(is.null(cmap)) {
    getPalette <- colorRampPalette(brewer.pal(11, "Spectral"))
    cmap <- getPalette(max(state_labels))
  }
  
  proto_inhibs <- matrix(0, max(cluster_assignments), ncol(cluster_weights))
  for(i in seq_len(max(cluster_assignments))) {
    if(sum(cluster_assignments == i) == 1) {
      proto_inhibs[i,] <- cluster_weights[which(cluster_assignments == i),]
    } else {
      proto_inhibs[i,] <- colMeans(cluster_weights[which(cluster_assignments == i),])
    }
  }
  
  nrow.plot <- ceiling(max(cluster_assignments) / ncol.plot)
  par(mfrow=c(nrow.plot,ncol.plot))
  for(i in seq_len(max(cluster_assignments))) {
    if(max(proto_inhibs[i,]) > 0.4) ymax <- max(proto_inhibs[i,])+0.1 
    else ymax <- 0.4
    barplot(proto_inhibs[i,], col=cmap, main='', xlab='', ylab = "Frequency (%)", ylim = c(0, ymax), cex.axis=1.5, cex.names = 2, cex.lab = ax.lab.sz, names.arg = seq_len(ncol(proto_inhibs)))
    
    title(xlab="Cell subtype", line=3.5, cex.lab=ax.lab.sz)
    title(main=sprintf("Group %s", intToUtf8(64+i)), line=0, cex.main=title.sz)
  }
}

#' @title Plot cell yield of each sample as bar plot
#' @description Plots cell yield (number of viable cells) of each single-cell sample in decreasing order as horizontal bar plot
#' @param myobj Phmed object containing expression data for each sample in 'data' slot
#' @param labels Vector containing group labels for samples (optional). If not provided, bars will be of uniform color (blue)
#' @param cmap Vector containing colors by which histogram bars should be colored (optional)
#' @param font_sz Scaling factor for font size of sample names in barplot
#' @param w Width of plot in inches
#' @param h Height of plot in inches
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
#' plotCellYield(my_phemdObj_final, labels=cluster_assignments, font_sz = 0.8)
#' 
plotCellYield <- function(myobj, labels=NULL, cmap=NULL, font_sz = 0.6, w=8, h=9.5) {
  nsample <- length(rawExpn(myobj))
  cell_yield <- vapply(rawExpn(myobj), nrow, integer(1L))
  
  order_idx <- order(cell_yield, decreasing=FALSE)
  cell_yield_ordered <- cell_yield[order_idx]
  snames_ordered <- sNames(myobj)[order_idx]
  
  
  par(mar=c(6,6,2,2))
  
  if(!is.null(labels)) {
    if(length(labels) != nsample) {
      stop('Error: length of "labels" vector must be equal to length of rawExpn(myobj)')
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

}


