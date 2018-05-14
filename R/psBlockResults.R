psBlockResults <- function(ps,sampleTypeVar="Sample_Type",blockVar="block"){
        
        names(sample_data(ps))[names(sample_data(ps))==sampleTypeVar] <- "Sample_Type"
        names(sample_data(ps))[names(sample_data(ps))==blockVar] <- "block"
        
        nblocks <- as.character(unique(sample_data(ps)$block))
        
        psByBlock <- list()
        
        for(i in 1:length(nblocks)){
                psByBlock[[i]] <- subset_samples(ps,block==nblocks[i])
        }
        
        psByBlock <- lapply(psByBlock,function(x){
                allzero <- apply(otu_table(subset_samples(x,Sample_Type=="Plasma")),1,function(x){all(x<1)})
                x <- prune_taxa(!allzero,x)
                return(x)
        })
        
        psNCbyBlock <- lapply(psByBlock,function(x){
                subset_samples(x,Sample_Type=="Control")
        })
        
        
        psallzeroInNC <- lapply(psNCbyBlock,function(x){
                allzero <- apply(otu_table(x),1,function(y){all(y<1)})
                z <- prune_taxa(allzero,x)
                return(z)
        })
        
        psPlByBlock <- lapply(psByBlock,function(x){
                subset_samples(x,Sample_Type=="Plasma")
        })
        
        
        rt <- as.list(psByBlock,psNCbyBlock,psallzeroInNC,psPlByBlock)
        return(rt)
}
