#' psBlockResults
#'
#' Prepare the phyloseq object for the Bayesian inference
#'
#' @param ps phyloseq object.
#' @param sampleTypeVar character, variable defines the 'Plasma' and 'Control' samples.
#' @param caselevels character vector, levels of the plasma samples, e.g. 'Healthy_Plasma' and 'Sick_Plasma'.
#' @param controllevel character, level of the control samples, e.g. 'Control'
#' @param sampleName character, variable defines the sample names.
#' @param blockVar character, variable defines the block of the samples. If there is only one block, create a variable with only one level.
#'
#' @return A list of phyloseq objects, phyloseq of negative control samples by block, phyloseq of taxa with prevalence zero in negative control samples by block, phyloseq of case samples by block.
#' @import phyloseq
#' @export
psBlockResults <- function(ps, sampleTypeVar = "Sample_Type", caselevels = "Plasma", 
    controllevel = "Control", sampleName = "SampleCode", blockVar = "block") {
    if (dim(otu_table(ps))[1] != ntaxa(ps)) {
        otu_table(ps) <- t(otu_table(ps))
    }
    
    Sample_Type <- NULL
    names(sample_data(ps))[names(sample_data(ps)) == sampleTypeVar] <- "Sample_Type"
    
    names(sample_data(ps))[names(sample_data(ps)) == blockVar] <- "block"
    
    names(sample_data(ps))[names(sample_data(ps)) == sampleName] <- "SampleCode"
    
    samdf <- sample_data(ps)
    gr <- samdf$block
    samdfL <- split(samdf, gr)
    blockSamples <- lapply(samdfL, function(x) {
        as.character(x$SampleCode)
    })
    
    psByBlock <- list()
    
    for (i in 1:length(blockSamples)) {
        psByBlock[[i]] <- prune_samples(blockSamples[[i]], ps)
    }
    
    psByBlock <- lapply(psByBlock, function(y) {
        sb_samples <- as.character(sample_data(y)$SampleCode)[sample_data(y)$Sample_Type %in% 
            caselevels]
        y <- prune_taxa(taxa_sums(prune_samples(sb_samples, y)) > 0, y)
        return(y)
    })
    
    psNCbyBlock <- lapply(psByBlock, function(z) {
        ct_samples <- as.character(sample_data(z)$SampleCode[sample_data(z)$Sample_Type %in% 
            controllevel])
        z <- prune_samples(ct_samples, z)
        return(z)
    })
    
    psallzeroInNC <- lapply(psNCbyBlock, function(m) {
        pm <- prune_taxa(taxa_sums(m) == 0, m)
        return(pm)
    })
    
    psPlByBlock <- lapply(psByBlock, function(x) {
        sb_samples <- as.character(sample_data(x)$SampleCode[sample_data(x)$Sample_Type %in% 
            caselevels])
        x <- prune_samples(sb_samples, x)
        return(x)
    })
    
    rt <- list(psByBlock, psNCbyBlock, psallzeroInNC, psPlByBlock)
    return(rt)
}
