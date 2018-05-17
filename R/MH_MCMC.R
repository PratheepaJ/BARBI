#' MH_MCMC
#'
#' @param iterations numeric, number of MCMC
#' @param k numeric vector, observed reads for all taxa in a sample
#' @param al numeric vector, estimated shape parameter for the contamination distribution
#' @param be numeric vector, estimated scale parameter for the contamination distribution
#'
#' @return array
#' @export
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
