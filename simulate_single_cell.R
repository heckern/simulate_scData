## functions for simulating single cell data

get_neg_binomial <- function(mu, n, size=5) {
    vals <- round(rnbinom(n, size=size, mu=mu), 0)
    return(vals)    
}


init_count_matrix <- function(n_samples, cell_tag, list_genes, list_mu_expression, size=5) {

    if(length(list_genes) != length(list_mu_expression)) {
        stop('Size of `list_genes` does not match size of `list_mu_expression`')
    }

    list_mu_expression[ list_mu_expression < 0] <- 0
    
    cell_labels <- paste0(cell_tag, as.character(seq(1:n_samples)))
    mat_counts <- t(sapply(list_mu_expression, get_neg_binomial, n=n_samples, size=5))

    rownames(mat_counts) <- list_genes
    colnames(mat_counts) <- cell_labels

    return(mat_counts)
}
