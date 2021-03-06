% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/calc_specificity.R
\name{calc_specificity}
\alias{calc_specificity}
\title{Calculate cell type specificity of genes}
\usage{
calc_specificity(datExpr, cell_type)
}
\arguments{
\item{datExpr}{expression data (numeric matrix of read counts),
a matrix or data frame with row and column names, in which rows are genes and columns are samples.}

\item{cell_type}{a vector of cell type annotation, whose dimension equals to ncol(datExpr).}
}
\value{
A dataframe including the specificity of each gene for each cell type, rows denoting genes.
}
\description{
Calculation of gene specificity for each cell type.
The inputs include a gene expression dataset and a cell type annotation information for cells.
}
\examples{
{}
}
