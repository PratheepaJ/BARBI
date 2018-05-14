alphaBetaNegControl <- function(psNCbyBlock){
        compAlphaBeta <- lapply(psNCbyBlock,function(x){
                otu_table(x) <- otu_table(x)+1
                psTodq <- phyloseq_to_deseq2(x,design = ~1)
                dq <- DESeq(psTodq,fitType = "local")##Account for the library sizes of negative control samples
                sizefac <- sizeFactors(dq)
                ot.tab <- t(t(otu_table(x))/sizefac)
                Nc <- round(min(colSums(ot.tab)),digits = 0)
                muhat.all <- assays(dq)[["mu"]]
                muhat <- rowMax(muhat.all)
                dis <- dispersions(dq)
                txname <- taxa_names(x)
                sam.mean <- apply(ot.tab,1,mean)
                sam.sd <- apply(ot.tab,1,sd)
                out <- list(muhat,dis,Nc,txname,sam.mean,sam.sd)
                names(out) <- c("muhat","dis","Nc","txname","sam.mean","sam.sd")
                return(out)
        })
        return(compAlphaBeta)
        
}
