#' Sampling posterior
#'
#' Sampling from the marginal reference posterior for the intensity of true signal
#'
#' @param gammaPrior_Cont A list of estimated distribution parameters of
#'        contamination intensity in the plasma samples in a particular
#'        block.
#' @param itera Numeric, number of MCMC samples.
#' @inheritParams alphaBetaContInPlasma
#' @inheritParams MH_MCMC
#' @importFrom doParallel registerDoParallel
#' @importFrom parallel detectCores
#' @importFrom BiocParallel register DoparParam bplapply
#' @return A list of estimated marginal reference posterior for the intensity of true signal.
#' @export
samplingPosterior <- function(psPlByBlock, blk, gammaPrior_Cont, itera = 10000, 
    ncores = ncores) {
    
    # registerDoParallel(detectCores()) register(DoparParam())
    
    sampleLst <- seq(1, nsamples(psPlByBlock[[blk]]))
    sampleLst <- as.list(sampleLst)
    
    sub_sampling_pos <- function(x, gammaPrior_Cont, blk, itera) {
        sam <- x
        taxa_list <- as.list(1:ntaxa(psPlByBlock[[blk]]))
        
        taxa_post <- lapply(taxa_list, function(taxa) {
            chain = MH_MCMC(itera = itera, k = as.numeric(gammaPrior_Cont[[sam]]$kij[taxa]), 
                al_c = gammaPrior_Cont[[sam]]$alpha_ij_c[taxa], be_c = gammaPrior_Cont[[sam]]$beta_ij_c[taxa], 
                startvalue_lamda_r = 0)
            return(chain)
        })
        
        return(taxa_post)
    }
    
    
    # taxa_post_all_sam <- bplapply(sampleLst, FUN=sub_sampling_pos,
    # gammaPrior_Cont=gammaPrior_Cont, blk=blk, itera=itera, BPPARAM =
    # MulticoreParam(workers = ncores, tasks = ncores))
    
    taxa_post_all_sam <- mclapply(sampleLst, FUN = sub_sampling_pos, gammaPrior_Cont = gammaPrior_Cont, 
        blk = blk, itera = itera, mc.cores = ncores)
    
    return(taxa_post_all_sam)
}
