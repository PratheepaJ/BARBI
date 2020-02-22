#' Bayesian latent Dirichlet allocation for 16S and metagenomics data.
#'
#' @param stan_data A list of stan data.
#' @param ... Arguments passed to `rstan::sampling` (e.g. iter, chains).
#'
#' @return An object of class `stanfit` returned by `rstan::sampling`
#' @export
#' @import Rcpp
#' @import methods
#' @importFrom rstan sampling
#' @useDynLib BARBI, .registration=TRUE
LDAtopicmodel = function(stan_data, ...){
        stan.data = stan_data
        stan.fit = rstan::sampling(stanmodels$lda,
                            data = stan.data, ...)
        return(stan.fit)
}


