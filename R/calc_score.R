#'@title  Calculate cell type score of genes
#'
#'@description Calculation of score for each gene using cell type specificity.
#'
#'@param specificity a matrix gotten from calc_specificity function (without NA) with row names,
#'rows denoting genes.
#'@param cell_type a vector of cell type annotation, whose dimension equals to ncol(datExpr).
#'@return A dataframe including cell type score for each gene, row denoting genes.
#'@export
#'@examples {}


calc_score <- function(specificity, cell_type) {
    spe = specificity
    index <- apply(spe, 1, function(x) {
        return(length(which(as.character(x) == "NaN")) == 0)
    })
    specificity = spe[index, ]
    if (nrow(specificity) < nrow(spe)) 
        warning("Remove genes whose specificity values are NaN in any cell type.")
    Ej_all = apply(specificity, 1, median)
    iqr = apply(specificity, 1, IQR)
    unique_cell = factor(unique(cell_type), levels = unique(cell_type))
    genes = rownames(specificity)
    score = matrix(nrow = length(genes), ncol = length(unique_cell))
    for (i in 1:length(unique_cell)) {
        cell <- unique_cell[i]
        Ej <- specificity[, cell]
        score[, i] = (Ej - Ej_all)/iqr
    }
    score = as.data.frame(score)
    dimnames(score) = list(genes, unique_cell)
    return(score)
}


