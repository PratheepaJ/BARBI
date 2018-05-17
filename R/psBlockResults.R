#' psBlockResults
#'
#' @param ps phyloseq object
#' @param sampleTypeVar factor, with two levels "Plasma" and "Control" samples
#' @param blockVar factor, block of the samples
#'
#' @return list, a list of phyloseq by block, phyloseq of negative control samples by block, phyloseq of taxa with prevalence zero in negative control samples by block, phyloseq of plasma samples by block
#' @export
psBlockResults <- function(ps,sampleTypeVar="Sample_Type",blockVar="block"){

        names(sample_data(ps))[names(sample_data(ps))==sampleTypeVar] <- "Sample_Type"
        names(sample_data(ps))[names(sample_data(ps))==blockVar] <- "block"

        samdf <- sample_data(ps)
        g <- samdf$block
        samdfL <- split(samdf,g)
        blockSamples <- lapply(samdfL,function(x){rownames(x)})


        psByBlock <- list()

        for(i in 1:length(blockSamples)){
                psByBlock[[i]] <- prune_samples(blockSamples[[i]],ps)
        }

        psByBlock <- lapply(psByBlock,function(x){
                allzero <- apply(otu_table(subset_samples(x,Sample_Type=="Plasma")),1,function(x){all(x<1)})
                x <- prune_taxa(!allzero,x)
                return(x)
        })

        psNCbyBlock <- lapply(psByBlock,function(x){
                x <- subset_samples(x,Sample_Type=="Control")
                return(x)
        })


        psallzeroInNC <- lapply(psNCbyBlock,function(x){
                allzero <- apply(otu_table(x),1,function(y){all(y<1)})
                z <- prune_taxa(allzero,x)
                return(z)
        })

        psPlByBlock <- lapply(psByBlock,function(x){
                x <- subset_samples(x,Sample_Type=="Plasma")
                return(x)
        })


        rt <- list(psByBlock,psNCbyBlock,psallzeroInNC,psPlByBlock)
        return(rt)
}
