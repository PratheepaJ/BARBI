alphaBetaContInPlasma <- function(psPlByBlock,psallzeroInNC,blk,alphaBetaNegControl){
        gammPrior <- list()
        compAlphaBeta <- alphaBetaNegControl
        for(j in 1:nsamples(psPlByBlock[[blk]])){
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

                xj <- otu_table(psPlByBlock[[blk]])[,j]
                txname <- taxa_names(psPlByBlock[[blk]])
                allzero <- ifelse(txname%in%taxa_names(psallzeroInNC[[blk]]),"Yes","No")
                gammPrior[[j]] <- list(ahat,bhat,xj,txname,allzero,muhat.pl)
        }

        return(gammPrior)
}
