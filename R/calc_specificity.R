#'@title Calculate cell type specificity of genes
#'
#'@description Calculation of gene specificity for each cell type.
#'The inputs include a gene expression dataset and a cell type annotation information for cells.
#'
#'@param datExpr expression data (numeric matrix of read counts),
#'a matrix or data frame with row and column names, in which rows are genes and columns are samples.
#'@param cell_type a vector of cell type annotation, whose dimension equals to ncol(datExpr).
#'@return A dataframe including the specificity of each gene for each cell type, rows denoting genes.
#'@importFrom edgeR DGEList calcNormFactors cpm
#'@export
#'@examples {}

calc_specificity <- function(datExpr, cell_type) {
    if (length(cell_type) != ncol(datExpr)) 
        stop("Length of cell_type must be equal to ncol(datExpr).")
    d <- DGEList(datExpr)
    d <- calcNormFactors(d, method = "TMM")
    filter_count <- cpm(d, prior.count = 10, normalized.lib.sizes = TRUE)
    unique_cell = factor(unique(cell_type), levels = unique(cell_type))
    minfc <- matrix(nrow = nrow(datExpr), ncol = length(unique_cell))
    genes <- rownames(datExpr)
    for (i in 1:length(unique_cell)) {
        cell <- unique_cell[i]
        list_other_cells <- unique_cell[which(unique_cell != cell)]
        # compute fold-changes for CTOI vs others
        fc_df <- data.frame(rep(1:nrow(filter_count)))
        for (j in 1:length(list_other_cells)) {
            fc <- rowMeans(filter_count[, which(cell_type %in% cell)])/rowMeans(filter_count[, which(cell_type %in% 
                list_other_cells[[j]])])
            if (!exists("fc_df")) {
                fc_df <- data.frame(fc)
            } else {
                fc_df[, j] <- fc
            }
        }
        fc_min <- apply(fc_df, 1, min)
        minfc[, i] <- fc_min
    }
    minfc <- as.data.frame(minfc)
    dimnames(minfc) = list(genes, unique_cell)
    return(minfc)
}

