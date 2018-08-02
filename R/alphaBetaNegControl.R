#' Compute gamma distirbution parameters
#'
#' Computes the gamma distribution parameters for the intensity of contamination using negative control samples in each block.
#'
#' @param psNCbyBlock A list of phyloseq objects with all negative control
#'        samples in each block.
#' @return A list of estimated distribution parameters for the intensity of contamination in negative control samples.
#' @examples
#' \dontrun{
#' alphaBetaNegControl(psNCByBlock)
#' }
#' @import DESeq2
#' @import SummarizedExperiment
#' @import Biobase
#' @importFrom stats median var runif dgamma rgamma
#' @export
alphaBetaNegControl <- function(psNCbyBlock, add1 = TRUE){

        compAlphaBeta <- lapply(psNCbyBlock, function(x){

                if(add1==TRUE){
                        otu_table(x) <- otu_table(x)+1

                        psTodq <- phyloseq_to_deseq2(x, design = ~1)

                        dq <- DESeq(psTodq, fitType = "local")

                        library_size_norm <- sizeFactors(dq)

                        ot.tab <- t(t(otu_table(x))/library_size_norm)

                        S_j0 <- round(median(colSums(ot.tab)), digits = 0)

                        mu_ij_0_all <- assays(dq)[["mu"]]

                        mu_ij_0 <- rowMedians(mu_ij_0_all)

                        gamma_ij_0 <- dispersions(dq)

                        species_name <- taxa_names(x)

                        sample_mean <- apply(ot.tab, 1, mean)

                        sample_var <- apply(ot.tab, 1, var)

                        alpha_ij_0 <- rep(.0001, length(mu_ij_0))
                        beta_ij_0 <- rep(1, length(mu_ij_0))

                        ind_not_na_of_mu_ij_0 <- which(!is.na(mu_ij_0))

                        alpha_ij_0[ind_not_na_of_mu_ij_0] <- 1/gamma_ij_0[ind_not_na_of_mu_ij_0]
                        beta_ij_0[ind_not_na_of_mu_ij_0] <- 1/(gamma_ij_0[ind_not_na_of_mu_ij_0]*mu_ij_0[ind_not_na_of_mu_ij_0])


                }else{
                        psTodq <- phyloseq_to_deseq2(x, design = ~1)

                        gm_mean <-  function(x, na.rm=TRUE){
                                exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
                        }

                        geoMeans <-  apply(counts(psTodq), 1, gm_mean)

                        psTodq <- estimateSizeFactors(psTodq, geoMeans = geoMeans)

                        # dq = DESeq2::DESeq(psTodq,fitType = "local",sfType = "poscounts")
                        dq <- DESeq(psTodq, fitType = "local")

                        library_size_norm <- sizeFactors(dq)

                        ot.tab <- t(t(otu_table(x))/library_size_norm)

                        S_j0 <- round(median(colSums(ot.tab)), digits = 0)

                        mu_ij_0_all <- assays(dq)[["mu"]]

                        mu_ij_0 <- rowMedians(mu_ij_0_all)

                        gamma_ij_0 <- dispersions(dq)

                        species_name <- taxa_names(x)

                        sample_mean <- apply(ot.tab, 1, mean)

                        sample_var <- apply(ot.tab, 1, var)

                        alpha_ij_0 <- rep(.0001, length(mu_ij_0))
                        beta_ij_0 <- rep(1, length(mu_ij_0))

                        ind_not_na_of_mu_ij_0 <- which(!is.na(mu_ij_0))

                        alpha_ij_0[ind_not_na_of_mu_ij_0] <- 1/gamma_ij_0[ind_not_na_of_mu_ij_0]
                        beta_ij_0[ind_not_na_of_mu_ij_0] <- 1/(gamma_ij_0[ind_not_na_of_mu_ij_0]*mu_ij_0[ind_not_na_of_mu_ij_0])

                }


                out <- list(mu_ij_0,
                            gamma_ij_0,S_j0,
                            species_name,sample_mean,
                            sample_var,
                            alpha_ij_0,
                            beta_ij_0)
                names(out) = c("mu_ij_0",
                               "gamma_ij_0",
                               "S_j0",
                               "species_name",
                               "sample_mean",
                               "sample_var",
                               "alpha_ij_0",
                               "beta_ij_0")
                return(out)
        })
        return(compAlphaBeta)
}
