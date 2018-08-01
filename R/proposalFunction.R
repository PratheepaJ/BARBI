#' proposalFunction
#'
#' @param al_p numeric, shape parameter for the proposal function
#' @param be_p numeric, scale parameter for the proposal function
#'
#' @return numeric, proposal density
#' @export
proposalFunction = function(al_p,
                            be_p){
        return(rgamma(1, shape=al_p, rate=be_p))
}
