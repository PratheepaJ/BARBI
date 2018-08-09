#' proposalFunction
#'
#' Proposed gamma density function
#'
#' @param al_p numeric, shape parameter for the proposed density function
#' @param be_p numeric, scale parameter for the proposed density function
#'
#' @return numeric, proposed density function
#' @export
proposalFunction = function(al_p,
                            be_p){
        return(rgamma(1, shape=al_p, rate=be_p))
}
