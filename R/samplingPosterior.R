#' samplingPosterior
#'
#' @param psPlByBlock list, phyloseq objects with plasma samples.
#' @param blk numeric, which block samples to be used
#' @param gammaPrior_Cont list, contamination distribution in the plasma samples
#' @param itera numeric, number of MCMC
#'
#' @return list, estiamted reference posterior for real reads in all samples
#' @export
samplingPosterior <- function(psPlByBlock,blk,gammaPrior_Cont,itera=10000){

        doParallel::registerDoParallel(parallel::detectCores())
        BiocParallel::register(BiocParallel::DoparParam())


        sampleLst <- seq(1,phyloseq::nsamples(psPlByBlock[[blk]]))
        sampleLst <- as.list(sampleLst)

        taxa_post_all_sam <- BiocParallel::bplapply(sampleLst,function(x){
                #taxa_post <- list()
                sam <- x
                taxa_list <- as.list(1:ntaxa(psPlByBlock[[blk]]))

                taxa_post <- lapply(taxa_list, function(taxa){
                        chain <- MH_MCMC(iterations = itera,k = as.numeric(gammaPrior_Cont[[sam]]$kj[taxa]),al_c = gammaPrior_Cont[[sam]]$alpha_ij_c[taxa],be_c = gammaPrior_Cont[[sam]]$beta_ij_c[taxa],startvalue_lamda_r=0)
                        return(chain)
                })

                # for(taxa in 1:phyloseq::ntaxa(psPlByBlock[[blk]])){
                #         chain <- MH_MCMC(iterations = itera,k = as.numeric(gammaPrior_Cont[[sam]][[3]][taxa]),al = gammaPrior_Cont[[sam]][[1]][taxa],be = gammaPrior_Cont[[sam]][[2]][taxa])
                #         taxa_post[[taxa]] <- chain
                # }
                return(taxa_post)
        })

        return(taxa_post_all_sam)
}
