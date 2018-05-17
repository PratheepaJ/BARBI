MH_MCMC <- function(iterations,k,al,be){
        startvalue <- 0
        chain <- array(dim = c(iterations+1,1))
        chain[1,] <- startvalue
        probab <- numeric()

        for (mc in 1:iterations){
                proposal <- proposalFunction(param=chain[mc,],al.p=(k+.5),be.p=1)

                pro <- exp(posterior(proposal,k=k,al=al,be=be) - posterior(chain[mc,],k=k,al=al,be=be))

                probab[mc] <- min(c(1,pro))

                if (runif(1) < probab[mc]){
                        chain[mc+1,] = proposal
                }else{
                        chain[mc+1,] = chain[mc,]
                }
        }
        return(chain)
}
