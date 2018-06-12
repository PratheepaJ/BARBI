#' alphaBetaNegControl
#'
#' @param psNCbyBlock list, phyloseq objects with all negative control samples in each block.
#'
#' @return list, estimated contamination density parameters using all negative control samples.
#' @export
alphaBetaNegControl <- function(psNCbyBlock){
        compAlphaBeta <- lapply(psNCbyBlock,function(x){
                # if(add1==TRUE){
                #         otu_table(x) <- otu_table(x)+1
                # }
                psTodq <- phyloseq::phyloseq_to_deseq2(x,design = ~1)

                # dq <- DESeq2::DESeq(psTodq,fitType = "local",sfType = "poscounts")
                dq <- DESeq2::DESeq(psTodq,fitType = "local")

                library_size_norm <- DESeq2::sizeFactors(dq)

                ot.tab <- t(t(otu_table(x))/library_size_norm)

                S_j0 <- round(median(colSums(ot.tab)),digits = 0)

                mu_ij_0_all <- SummarizedExperiment::assays(dq)[["mu"]]

                mu_ij_0 <- Biobase::rowMedians(mu_ij_0_all)

                gamma_ij_0 <- DESeq2::dispersions(dq)

                species_name <- phyloseq::taxa_names(x)

                sample_mean <- apply(ot.tab,1,mean)

                sample_var <- apply(ot.tab,1,var)

                alpha_ij_0 <- rep(.0001, length(mu_ij_0))
                beta_ij_0 <- rep(1, length(mu_ij_0))

                ind_not_na_of_mu_ij_0 <- which(!is.na(mu_ij_0))

                alpha_ij_0[ind_not_na_of_mu_ij_0] <- 1/gamma_ij_0[ind_not_na_of_mu_ij_0]
                beta_ij_0[ind_not_na_of_mu_ij_0] <- 1/(gamma_ij_0[ind_not_na_of_mu_ij_0]*mu_ij_0[ind_not_na_of_mu_ij_0])

                out <- list(mu_ij_0,gamma_ij_0,S_j0,species_name,sample_mean,sample_var,alpha_ij_0,beta_ij_0)
                names(out) <- c("mu_ij_0","gamma_ij_0","S_j0","species_name","sample_mean","sample_var","alpha_ij_0","beta_ij_0")
                return(out)
        })
        return(compAlphaBeta)

}
