#' Compute gamma distribution parameters in plasma samples
#'
#' Given the estimated gamma distribution parameters for the intensity of contamination in negative control samples, \code{alphaBetaContInPlasma}
#' computes the
#' gamma distribution parameters in the plasma sample using the scaling
#' property of the gamma distribution.
#'
#' @param psPlByBlock A list of phyloseq objects with plasma samples.
#' @param psallzeroInNC A list of phyloseq objects with all taxa with no
#'        prevalence in negative control samples.
#' @param blk Numeric, which block samples to be used.
#' @param alphaBetaNegControl A list of estimated distribution parameters for the intensity of contamination in negative control samples returned by \code{alphaBetaNegControl}.
#' @return A list of estimated distribution parameters for the intensity of contamination in the plasma samples in a selected block.
#' @import phyloseq
#' @export
#' @section Note:
#' If all the samples are from one block, create a variable with one block
#' level.
#' @seealso \code{\link{alphaBetaNegControl}}

alphaBetaContInPlasma <- function(psPlByBlock, psallzeroInNC, blk, alphaBetaNegControl) {
    ps_blk <- psPlByBlock[[blk]]
    if (dim(otu_table(ps_blk))[1] != ntaxa(ps_blk)) {
        otu_table(ps_blk) = t(otu_table(ps_blk))
    }

    compAlphaBeta <- alphaBetaNegControl
    compAlphaBeta_blk <- compAlphaBeta[[blk]]
    samples_in_block <- as.list(1:nsamples(ps_blk))
    gammaPrior_c <- lapply(samples_in_block, function(x) {
        kij <- otu_table(ps_blk)[, x]
        S_j <- sum(kij)
        S_j0 <- compAlphaBeta_blk$S_j0

        mu_ij_0 <- compAlphaBeta_blk$mu_ij_0
        ind_not_na_of_mu_ij_0 <- which(!is.infinite(mu_ij_0))

        alpha_ij_c <- rep(1e-04, length(mu_ij_0))
        beta_ij_c <- rep(1, length(mu_ij_0))

        alpha_ij_c[ind_not_na_of_mu_ij_0] <- S_j/S_j0 * compAlphaBeta_blk$alpha_ij_0[ind_not_na_of_mu_ij_0]
        beta_ij_c[ind_not_na_of_mu_ij_0] <- compAlphaBeta_blk$beta_ij_0[ind_not_na_of_mu_ij_0]

        species_name <- compAlphaBeta_blk$species_name

        all_zero_neg_controls <- ifelse(species_name %in% taxa_names(psallzeroInNC[[blk]]),
            "Yes", "No")

        rt <- list(alpha_ij_c, beta_ij_c, kij, species_name, all_zero_neg_controls,
            S_j, S_j0)
        names(rt) <- c("alpha_ij_c", "beta_ij_c", "kij", "species_name", "all_zero_neg_controls",
            "S_j", "S_j0")

        return(rt)
    })

    return(gammaPrior_c)
}



