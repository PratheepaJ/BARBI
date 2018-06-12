#' MH_MCMC
#'
#' @param iterations numeric, number of MCMC
#' @param k numeric vector, observed reads for all taxa in a sample
#' @param al_c numeric vector, estimated shape parameter for the contamination distribution
#' @param be_c numeric vector, estimated scale parameter for the contamination distribution
#'
#' @return array
#' @export
MH_MCMC <- function(iterations,k,al_c,be_c,startvalue_lamda_r=0){
        startvalue_lamda_r <- startvalue_lamda_r
        chain_lamda_r <- array(dim = c(iterations+1,1))
        chain_lamda_r[1,] <- startvalue_lamda_r

        probab <- numeric()

        #diff_shape_para_prop <- k
        diff_shape_para_prop <- (k-al_c)
        if(diff_shape_para_prop <= 0){ diff_shape_para_prop <- (k+.5)}

        for (mc in 1:iterations){
                proposal <- proposalFunction(al_p= diff_shape_para_prop,be_p = 1)

                pro <- exp(posterior(lambda_r=proposal,k=k,al_c=al_c,be_c=be_c) - posterior(lambda_r=chain_lamda_r[mc,],k=k,al_c=al_c,be_c=be_c))

                probab[mc] <- min(c(1,pro))

                if (runif(1) < probab[mc]){
                        chain_lamda_r[mc+1,] = proposal
                }else{
                        chain_lamda_r[mc+1,] = chain_lamda_r[mc,]
                }
        }
        return(chain_lamda_r)
}
