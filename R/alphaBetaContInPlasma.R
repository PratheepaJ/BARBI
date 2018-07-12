#' Compute gamma distribution parameters
#'
#' Given the estimated gamma distribution parameters of contaminant
#' intensity in negative control samples, \code{alphaBetaContInPlasma}
#' computes the
#' gamma distribution parameters in the plasma sample using the scaling
#' property of the gamma distribution.
#'
#' @param psPlByBlock A list of phyloseq objects with plasma samples.
#' @param psallzeroInNC A list of phyloseq objects with all taxa with no
#'        prevalence in negative control samples.
#' @param blk Numeric, which block samples to be used.
#' @param alphaBetaNegControl A list of estimated contamination
#'        distribution parameters using all negative control samples.

#' @return A list of estimated distribution parameters of contamination
#'         intensity in the plasma samples in a particular block.
#' @import phyloseq
#' @export
#' @section Note:
#' If all the samples are from one block, create a variable with one block
#' level.
#' @seealso \code{\link{alphaBetaNegControl}}

alphaBetaContInPlasma = function(psPlByBlock,psallzeroInNC,blk,alphaBetaNegControl){
        if (dim(otu_table(psPlByBlock[[blk]]))[1] != ntaxa(psPlByBlock[[blk]])){otu_table(psPlByBlock[[blk]]) = t(otu_table(psPlByBlock[[blk]]))}

        compAlphaBeta = alphaBetaNegControl
        samples_in_block = as.list(1:nsamples(psPlByBlock[[blk]]))
        gammaPrior_c = lapply(samples_in_block, function(x){
                kj = phyloseq::otu_table(psPlByBlock[[blk]])[,x]
                S_j = sum(kj)
                S_j0 = compAlphaBeta[[blk]]$S_j0

                mu_ij_0 = compAlphaBeta[[blk]]$mu_ij_0
                gamma_ij_0 = compAlphaBeta[[blk]]$gamma_ij_0

                alpha_ij_c = rep(.0001, length(mu_ij_0))
                beta_ij_c = rep(1, length(mu_ij_0))

                ind_not_na_of_mu_ij_0 = which(!is.na(mu_ij_0))

                alpha_ij_c[ind_not_na_of_mu_ij_0     ] = S_j/S_j0*(1/gamma_ij_0[ind_not_na_of_mu_ij_0     ])
                beta_ij_c[ind_not_na_of_mu_ij_0     ] = 1/(gamma_ij_0[ind_not_na_of_mu_ij_0     ]*mu_ij_0[ind_not_na_of_mu_ij_0     ])

                #       muhat.pl = alpha_ij_c/beta_ij_c
                species_name = compAlphaBeta[[blk]]$species_name

                all_zero_neg_controls = ifelse(species_name%in%phyloseq::taxa_names(psallzeroInNC[[blk]]),"Yes","No")

                rt = list(alpha_ij_c,beta_ij_c,kj,species_name,all_zero_neg_controls,S_j,S_j0)
                names(rt) = c("alpha_ij_c","beta_ij_c","kj","species_name","all_zero_neg_controls","S_j","S_j0")

                return(rt)
        })

        return(gammaPrior_c)
}
