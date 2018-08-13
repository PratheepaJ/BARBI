#' Make summaires from the BARBI results
#'
#' @param burnIn numeric, number of states from MCMC to discard.
#' @param num_blks numeric, number of blocks.
#' @param cov.pro numeric, coverage probability to construct the highest posterior density interval.
#' @param psByBlock list, phyloseq objects by blocks.
#' @param gammaPrior_all_blks list, estimated distribution parameters for the intensity of contamination in the plasma samples in all blocks.
#' @param post_all_blocks list, posterior results
#' @param mak_tab logical, whether make a table for each sample with true signal.
#' @inheritParams MH_MCMC
#'
#' @return tables with true signal and a list of true signal
#' @export

makeTable <- function(itera, burnIn, num_blks, cov.pro, psByBlock, gammaPrior_all_blks, post_all_blocks, mak_tab=TRUE){
        all_real_taxa_lt <- list()
        for(blk in 1:num_blks){

                taxa_post_all_sam <- post_all_blocks[[blk]]
                gammPrior <- gammaPrior_all_blks[[blk]]

                total_summary_table <- NULL

                all_real_taxa <- character()

                for(sam in 1:nsamples(psPlByBlock[[blk]])){

                        taxa_post <- taxa_post_all_sam[[sam]]
                        acceptance <- list()
                        exp_post_s <- list()
                        lower.s <- list()
                        upper.s <- list()
                        lower.b <- list()
                        upper.b <- list()
                        all.zero.nc <- list()

                        for(taxa in 1:length(taxa_post)){
                                burnIn  <- burnIn
                                acceptance[[taxa]]  <-  1-mean(duplicated(taxa_post[[taxa]][-(1:burnIn),]))

                                exp_post_s[[taxa]] <- mean(taxa_post[[taxa]][-(1:burnIn),])

                                hdi.v <- hdi(taxa_post[[taxa]][-(1:burnIn),], credMass = cov.pro)
                                lower.s[[taxa]] <- round(hdi.v[1], digits = 0)
                                upper.s[[taxa]] <- round(hdi.v[2], digits = 0)

                                hdi.b <- hdi(rgamma((itera-burnIn+1), shape=gammPrior[[sam]][[1]][taxa], rate = gammPrior[[sam]][[2]][taxa]), credMass = cov.pro)
                                lower.b[[taxa]] <- round(hdi.b[1], digits = 0)
                                upper.b[[taxa]] <- round(hdi.b[2], digits = 0)


                                all.zero.nc[[taxa]] <-  gammPrior[[sam]][[5]][taxa]
                        }


                        df <- data.frame(Species=taxa_names(psPlByBlock[[blk]]),
                                         xj=as.numeric(gammPrior[[sam]][[3]]),
                                         l.s=unlist(lower.s),
                                         u.s=unlist(upper.s),
                                         l.b=unlist(lower.b),
                                         u.b=unlist(upper.b),
                                         all.zero.nc=unlist(all.zero.nc)
                        )

                        df <- arrange(filter(df,(l.s>u.b)&(l.s>0)), desc(xj))


                        if(dim(df)[1]==0){
                                df <- data.frame(Species="Negative",
                                                 xj="Negative",
                                                 l.s="Negative",
                                                 u.s="Negative",
                                                 l.b ="Negative",
                                                 u.b="Negative",
                                                 all.zero.nc = "Negative")
                        }


                        all_real_taxa <- c(all_real_taxa,
                                           as.character(df$Species))
                        if(mak_tab){
                                filname <- paste("./",
                                                 sample_names(psPlByBlock[[blk]])[sam],
                                                 ".png",
                                                 sep="")

                                png(filname, height = 600, width = 750)

                                df.p <- tableGrob(df)
                                title <- textGrob(sample_names(psPlByBlock[[blk]])[sam], gp = gpar(fontsize = 12))

                                padding <- unit(0.5,"line")

                                df.p <- gtable_add_rows(
                                        df.p, heights = grobHeight(title) + padding, pos = 0
                                )

                                df.p <- gtable_add_grob(
                                        df.p, list(title),
                                        t = 1, l = 1, r = ncol(df.p)
                                )

                                grid.newpage()
                                grid.draw(df.p)
                                dev.off()
                        }


                        all_real_taxa <- unique(all_real_taxa)
                }

                all_real_taxa_lt[[blk]] <- all_real_taxa
        }

        return(all_real_taxa_lt)
}
