\name{all_expn_data}
\alias{all_expn_data}
\docType{data}
\title{Single-cell RNA-seq expression data for melanoma samples}
\description{This dataset contains normalized single-cell RNA-seq expression data for 19 melanoma samples (immune cells).}
\usage{
data(all_expn_data)
}
\format{
  A list of length 19 with each element representing a distinct sample. Each list element (sample) is a matrix with dimension num_genes x num_cells.
}
\source{https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE72056}
\references{Tirosh, I. et al. Dissecting the multicellular ecosystem of metastatic melanoma by single-cell RNA-seq. Science 352, 189–196 (2016)}
\keyword{datasets}
