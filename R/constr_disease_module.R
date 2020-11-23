#'@title  Construct cell type-specific disease gene module
#'
#'@description   Identification of cell type-specific disease gene module from cell type-specific gene network,
#'               using cell type scores of genes and disease-associated genes.
#'
#'@param Cts_list a list, gotten from constr_specific_net function, containing cell type-specific gene network, a two-column matrix, character or numeric, with each row defining one edge,
#'and genescore, a dataframe of cell type scores of genes in the cell type-specific gene network.
#'@param disease_gene a vector of disease-associated genes. The keytype of disease_gene must be consistent with the one of network.
#'@param itera a number (a positive integer) of resampling iterations (by default 1000).
#'@param returnAllmodules a boolean value determining if all candidate cell type-specific disease gene modules should be return or not (by default FALSE, only return the significant ones).
#'@param padjCutoff a cutoff of fdr-adjusted P value.
#'@return A list containing disease_module, cell type-specific disease gene module, each row defining one edge;
#' gene_score, cell type scores of genes in the cell type-specific disease gene module; Pvalue, significance of the disease module
#' and fdr-adjusted P value.
#'@importFrom igraph graph_from_edgelist components subcomponent
#'@export
#'@examples {}
constr_disease_module <- function(Cts_list, disease_gene, itera = 1000, returnAllmodules = FALSE, padjCutoff = 0.1) {
    y <- disease_gene
    cells <- names(Cts_list)
    n = length(cells)
    out = list()
    # identify all candidate disease gene modules
    p_value = list()
    for (i in 1:n) {
        cell = cells[i]
        out_cell = list()
        net <- Cts_list[[cell]]$specific_net
        score <- Cts_list[[cell]]$gene_score
        spe_gene <- score$gene
        inters_g <- intersect(spe_gene, y)
        ntotal <- length(inters_g)
        ix <- (net$V1 %in% inters_g) & (net$V2 %in% inters_g)
        net_disease = net[ix, ]
        if (nrow(net_disease) > 0) {
            el <- as.matrix(net_disease, nc = 2, byrow = TRUE)
            g = igraph::graph_from_edgelist(el)
            clu <- igraph::components(g)
            modules = unique(clu$membership)
            for (index in modules) {
                s_ob = length(which(clu$membership == index))
                subg <- igraph::subcomponent(g, igraph::V(g)[which(clu$membership == index)])
                disemoduleGene <- igraph::as_ids(subg)
                sub_net = net_disease[(net_disease$V1 %in% disemoduleGene) & (net_disease$V2 %in% disemoduleGene),
                  ]
                gene_score = score[score$gene %in% disemoduleGene, ]
                ## calculate disease module significance
                srand = c()
                for (r in 1:itera) {
                  gene_random = sample(spe_gene, ntotal)
                  inModule_rand = (net$V1 %in% gene_random) & (net$V2 %in% gene_random)
                  TOM_rand = net[inModule_rand, ]
                  if (nrow(TOM_rand) > 0) {
                    el_rand <- as.matrix(TOM_rand, nc = 2, byrow = TRUE)
                    g_rand = igraph::graph_from_edgelist(el_rand)
                    clu_rand <- igraph::components(g_rand)
                    s_rand = max(clu_rand$csize)
                  } else {
                    s_rand = 0
                  }
                  srand[r] = s_rand
                }  #end itera
                p_value[[paste0(cell, index)]] <- length(which(srand >= s_ob))/itera
                out_cell[[as.character(index)]][["disease_module"]] = sub_net
                out_cell[[as.character(index)]][["gene_score"]] = gene_score
                out_cell[[as.character(index)]][["Pvalue"]] = length(which(srand >= s_ob))/itera
            }  #end modules
        } else {
            warning(paste0("The given sthre cause the disease module of ", cell, " to contain no genes.Try a smaller threshold.\n"))
        }
        out[[cell]] = out_cell
    }
    p.adj = data.frame(cell = names(p_value), p_adj = p.adjust(as.vector(unlist(p_value)), "fdr"), stringsAsFactors = FALSE)
    res = list()
    for (cell in names(out)) {
        temp = out[[cell]]
        for (name in names(temp)) {
            module = paste0(cell, name)
            padj = p.adj$p_adj[which(p.adj$cell == module)]
            out[[cell]][[name]][["p.adjust"]] = padj
            if (padj < padjCutoff)
                {
                  res[[cell]][[name]] = out[[cell]][[name]]
                }  #end if
        }  #end for name
    }  #end for cell

    # only significant ones can be returned.
    if (returnAllmodules == FALSE) {
        if (length(res) == 0)
            warning("No significant cell type-specific disease gene modules returned. One can set returnAllmodules=TRUE to check all candidate cell type-specific disease gene modules.")
        return(res)
    }
    return(out)
}


