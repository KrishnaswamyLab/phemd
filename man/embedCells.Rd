% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/functions-core.R
\name{embedCells}
\alias{embedCells}
\title{Generate cell-state embedding}
\usage{
embedCells(
  obj,
  cell_model = c("monocle2", "seurat", "phate"),
  data_model = "negbinomial_sz",
  phate_ncluster = 8,
  phate_cluster_seed = NULL,
  ...
)
}
\arguments{
\item{obj}{'Phemd' object containing aggregated data}

\item{cell_model}{Method to use to generate cell-state embedding. Currently supports "phate" and "monocle2". If using the Seurat to model the cell-state space, please identify cell subtypes as outlined in the Seurat software package and then use the \code{bindSeuratObj} function.}

\item{data_model}{Only relevant if cell_model = "monocle2". One of the following: 'negbinomial_sz', 'negbinomial', 'tobit', 'uninormal', 'gaussianff'. See "Family Function" table at the following link for more details on selecting the proper one. \url{http://cole-trapnell-lab.github.io/monocle-release/docs/#getting-started-with-monocle}}

\item{phate_ncluster}{Only relevant if cell_model = "phate". Number of cell state clusters to return when using PHATE}

\item{phate_cluster_seed}{Only relevant if cell_model = "phate". Seed to use when performing cell state clustering (optional)}

\item{...}{Additional parameters to be passed to \code{reduceDimension} function for Monocle or \code{phate} function for PHATE}
}
\value{
Same as input 'Phemd' object containing additional cell-state embedding object
}
\description{
Takes as input a Phemd object with aggregated data and returns updated object containing cell-state embedding
}
\details{
\code{aggregateSamples} needs to be called before running this function.
}
\examples{
my_phemdObj <- createDataObj(all_expn_data, all_genes, as.character(snames_data))
my_phemdObj_lg <- removeTinySamples(my_phemdObj, 10)
my_phemdObj_lg <- aggregateSamples(my_phemdObj_lg, max_cells=1000)
my_phemdObj_lg <- embedCells(my_phemdObj_lg, cell_model='monocle2', data_model = 'gaussianff', sigma=0.02, maxIter=2)
}
