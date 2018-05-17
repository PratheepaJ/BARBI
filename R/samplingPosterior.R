#' samplingPosterior
#'
#' @param psPlByBlock list, phyloseq objects with plasma samples.
#' @param blk numeric, which block samples to be used
#' @param gammPrior list, contamination distribution in the plasma samples
#' @param iter numeric, number of MCMC
#'
#' @return list, estiamted reference posterior for real reads in all samples
#' @export
samplingPosterior <- function(psPlByBlock,blk,gammPrior,iter=10000){
        doParallel::registerDoParallel(parallel::detectCores())
        BiocParallel::register(BiocParallel::DoparParam())
        taxa.post <- list()
        sam <- seq(1,nsamples(psPlByBlock[[blk]]))
        sam <- as.list(sam)
        taxa_post_all_sam <- bplapply(sam,function(x){
                sam <- x
                for(taxa in 1:ntaxa(psPlByBlock[[blk]])){
                        chain <- MH_MCMC(iterations = iter,k=as.numeric(gammPrior[[sam]][[3]][taxa]),al=gammPrior[[sam]][[1]][taxa],be=gammPrior[[sam]][[2]][taxa])
                        taxa.post[[taxa]] <- chain
                }
                return(taxa.post)
        })
        return(taxa_post_all_sam)
}
