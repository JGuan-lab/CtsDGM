#'@title  Construct cell type-specific gene network
#'
#'@description   Construction of cell type-specific gene network from a referenced gene-gene interaction network,
#'               using cell type scores of genes.
#'
#'@param net a referenced gene-gene interaction network, a two-column matrix, character or numeric, each row defining one edge.
#'@param genescore  a dataframe of cell type scores of genes, gotten from calc_score function.
#'@param sthre a cutoff of cell type scores used for getting cell type-specific genes (by default 0).
#'@return A list containing specific_net, cell type-specific gene network each row defining one edge;
#' and gene_score, cell type score for each gene in the cell type-specific gene network.
#'@export
#'@examples {}

constr_specific_net <- function(net, genescore, sthre = 0) {
    if (ncol(net) == 2) 
        colnames(net) = c("V1", "V2")
    if (!is.null(colnames(genescore))) 
        cells <- colnames(genescore)
    if (!is.null(rownames(genescore))) 
        genes <- rownames(genescore)
    out = list()
    for (cell in cells) {
        score = genescore[, cell]
        spe_gene <- genes[score > sthre]
        ix <- (net$V1 %in% spe_gene) & (net$V2 %in% spe_gene)
        net_spe = net[ix, ]
        uni_g <- unique(c(net_spe$V1, net_spe$V2))
        single_node = spe_gene[!spe_gene %in% uni_g]
        single_node = data.frame(V1 = single_node, V2 = rep("", times = length(single_node)))
        spe_net = rbind(net_spe, single_node)
        if (nrow(spe_net) == 0) 
            warning(paste0("The given sthre cause the specific net of ", cell, " to contain no genes.Try a smaller threshold.\n"))
        score = data.frame(gene = spe_gene, score = score[score > sthre])
        out[[cell]][["specific_net"]] = spe_net
        out[[cell]][["gene_score"]] = score
    }
    return(out)
}

