#' Compute gamma distirbution parameters in negative control samples.
#'
#' Computes the gamma distribution parameters for the intensity of contamination using negative control samples in each block.
#'
#' @param psNCbyBlock A list of phyloseq objects with all negative control
#'        samples in each block.
#' @param add1 logical, whether adding one to the otu table in the phyloseq
#' @param stringent Logical, TRUE/FALSE: (i) TRUE use the sample mean of non-zero abundances in the negative control samples for the taxa with mean less than 1, and (ii) FALSE use the less than one mean.
#'
#' @return A list of estimated distribution parameters for the intensity of contamination in negative control samples.
#' @examples
#' \dontrun{
#' alphaBetaNegControl(psNCByBlock)
#' }
#' @import DESeq2
#' @importFrom SummarizedExperiment assays
#' @importFrom Biobase rowMedians
#' @importFrom stats median var runif dgamma rgamma
#' @export
alphaBetaNegControl <- function(psNCbyBlock, stringent = FALSE) {

    compAlphaBeta <- lapply(psNCbyBlock, function(x) {

            ps.to.dq <- phyloseq_to_deseq2(x, design = ~1)

            geo.mean <- function(y) {
                    if(all(y == 0)){
                            val <- 0
                    }else{
                            val <- exp(sum(log(y[y > 0]))/length(y))
                    }
                    return(val)
            }

            geom.mean.row <- apply(counts(ps.to.dq), 1, FUN = geo.mean)

            ps.to.dq <- estimateSizeFactors(ps.to.dq, geoMeans = geom.mean.row)

            dq <- DESeq(ps.to.dq, fitType = "local", minReplicatesForReplace= Inf)

            library.size.norm <- sizeFactors(dq)

            ot.tab <- t(t(otu_table(x))/library.size.norm)

            S_j0 <- round(min(colSums(ot.tab)), digits = 0)

            mu_ij_0_all <- assays(dq)[["mu"]]

            mu_ij_0 <- apply(mu_ij_0_all, 1, function(x){
                    max(x, na.rm = T)
            })

            gamma_ij_0 <- dispersions(dq)

            species_name <- taxa_names(x)

            sample_mean <- apply(ot.tab, 1, function(y){
                    if(all(y == 0)){
                            0
                    }else{
                            mean(y[y > 0])
                    }
            })

            sample_var <- apply(ot.tab, 1, function(y){
                    if(all(y == 0)){
                            0
                    }else{
                            var(y)
                    }
            })

            disp <- numeric(0)

            for(i in 1:length(sample_var)){
                    disp[i] <- (sample_var[i] - sample_mean[i])/(sample_mean[i])^2
            }

            alpha_ij_0 <- rep(1e-04, length(mu_ij_0))
            beta_ij_0 <- rep(1, length(mu_ij_0))

            ind_not_na_of_mu_ij_0 <- which(!is.infinite(mu_ij_0))#max produce NA for all zeros

            ind_less_one_mu_ij_0 <- which(abs(mu_ij_0) < 1)

            if(stringent){
                    alpha_ij_0[ind_not_na_of_mu_ij_0] <- 1/gamma_ij_0[ind_not_na_of_mu_ij_0]
                    beta_ij_0[ind_not_na_of_mu_ij_0] <- 1/(gamma_ij_0[ind_not_na_of_mu_ij_0] *
                                    mu_ij_0[ind_not_na_of_mu_ij_0])

                    beta_ij_0[ind_less_one_mu_ij_0] <- 1/(gamma_ij_0[ind_less_one_mu_ij_0]*sample_mean[ind_less_one_mu_ij_0])

            }else{
                    alpha_ij_0[ind_not_na_of_mu_ij_0] <- 1/gamma_ij_0[ind_not_na_of_mu_ij_0]
                    beta_ij_0[ind_not_na_of_mu_ij_0] <- 1/(gamma_ij_0[ind_not_na_of_mu_ij_0] *
                                    mu_ij_0[ind_not_na_of_mu_ij_0])
            }

        out <- list(mu_ij_0, gamma_ij_0, S_j0, species_name, sample_mean,
            sample_var, alpha_ij_0, beta_ij_0)
        names(out) = c("mu_ij_0", "gamma_ij_0", "S_j0", "species_name", "sample_mean",
            "sample_var", "alpha_ij_0", "beta_ij_0")
        return(out)
    })
    return(compAlphaBeta)
}

