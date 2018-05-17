#' psBlockResults
#'
#' @param ps phyloseq object
#' @param sampleTypeVar factor, with two levels "Plasma" and "Control" samples
#' @param blockVar factor, block of the samples
#'
#' @return list, a list of phyloseq by block, phyloseq of negative control samples by block, phyloseq of taxa with prevalence zero in negative control samples by block, phyloseq of plasma samples by block
#' @export
psBlockResults <- function(ps,sampleTypeVar="Sample_Type",blockVar="block"){

        names(phyloseq::sample_data(ps))[names(phyloseq::sample_data(ps))==sampleTypeVar] <- "Sample_Type"
        names(phyloseq::sample_data(ps))[names(phyloseq::sample_data(ps))==blockVar] <- "block"

        samdf <- phyloseq::sample_data(ps)
        g <- samdf$block
        samdfL <- split(samdf,g)
        blockSamples <- lapply(samdfL,function(x){rownames(x)})


        psByBlock <- list()

        for(i in 1:length(blockSamples)){
                psByBlock[[i]] <- phyloseq::prune_samples(blockSamples[[i]],ps)
        }

        psByBlock <- lapply(psByBlock,function(x){
                allzero <- apply(phyloseq::otu_table(phyloseq::subset_samples(x,Sample_Type=="Plasma")),1,function(x){all(x<1)})
                x <- phyloseq::prune_taxa(!allzero,x)
                return(x)
        })

        psNCbyBlock <- lapply(psByBlock,function(x){
                x <- phyloseq::subset_samples(x,Sample_Type=="Control")
                return(x)
        })


        psallzeroInNC <- lapply(psNCbyBlock,function(x){
                allzero <- apply(otu_table(x),1,function(y){all(y<1)})
                z <- phyloseq::prune_taxa(allzero,x)
                return(z)
        })

        psPlByBlock <- lapply(psByBlock,function(x){
                x <- phyloseq::subset_samples(x,Sample_Type=="Plasma")
                return(x)
        })


        rt <- list(psByBlock,psNCbyBlock,psallzeroInNC,psPlByBlock)
        return(rt)
}
