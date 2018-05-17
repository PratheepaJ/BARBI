#' posterior
#'
#' @param s numeric, intensity of real read.
#' @param k numeric, observed read for a taxa in a sample.
#' @param al numeric, estimated shape parameter for the contamination distribution for each taxa.
#' @param be numeric, estimated scale parameter for the contamination distribution for each taxa.
#'
#' @return numeric, posterior density
#' @export
posterior <- function(s,k,al,be){
        ga.mu <- al/be

        if(k==0){
                pos.v <- dgamma(s,shape=.5,rate = 1)

        }else if(k==1){
                pos.v <- exp(-log(factorial(.5))+.5*log(s+ga.mu)-(s+ga.mu))
        }else{
                c.m <- numeric()
                c.m[1] <- 0

                for(m in 2:k){
                        c.m[m] <- log(m)+c.m[m-1]
                }

                pos.v <- exp(-c.m[(k-1)]+(k-.5)*log(s+ga.mu)-(s+ga.mu))
        }

        return(pos.v)
}
