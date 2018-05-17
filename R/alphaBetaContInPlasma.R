#' alphaBetaContInPlasma
#'
#' @param psPlByBlock list, phyloseq objects with plasma samples.
#' @param psallzeroInNC list, phyloseq objects with all taxa with no prevalence in negative controls.
#' @param blk numeric, which block samples to be used.
#' @param alphaBetaNegControl list, estimated contamination distribution parameters using all negative control samples.
#'
#' @return list, estimated contamination distribution parameters for the plasma samples.
#'
#' @export

alphaBetaContInPlasma <- function(psPlByBlock,psallzeroInNC,blk,alphaBetaNegControl){
        gammPrior <- list()
        compAlphaBeta <- alphaBetaNegControl
        for(j in 1:phyloseq::nsamples(psPlByBlock[[blk]])){
                Np <- colSums(otu_table(psPlByBlock[[blk]]))[j]
                Nc <- compAlphaBeta[[blk]][[3]]
                mhat <- compAlphaBeta[[blk]][[1]]
                dishat <- compAlphaBeta[[blk]][[2]]
                ahat <- Np/Nc*(1/dishat) ## account for dependency of plasma sample's library size and negative control samples
                bhat <- 1/(dishat*mhat)

                muhat.pl <- ahat/bhat
                # ind <- which(muhat.pl < 1)# if muhat.pl less than 1, set it to negative control muhat. i.e not considering library size of plasma sample
                # ahat[ind] <- 1/dishat[ind]
                # bhat[ind] <- 1/(dishat[ind]*mhat[ind])
                xj <- phyloseq::otu_table(psPlByBlock[[blk]])[,j]
                txname <- phyloseq::taxa_names(psPlByBlock[[blk]])
                allzero <- ifelse(txname%in%phyloseq::taxa_names(psallzeroInNC[[blk]]),"Yes","No")
                gammPrior[[j]] <- list(ahat,bhat,xj,txname,allzero,muhat.pl)
        }

        return(gammPrior)
}
