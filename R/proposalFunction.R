#' proposalFunction
#'
#' @param param numeric, a proposal value
#' @param al.p numeric, shape parameter of the proposal function
#' @param be.p numeric, scale parameter of the proposal function
#'
#' @return numeric, proposal density
#' @export
proposalFunction <- function(al_p,be_p){
        return(rgamma(1,shape=al_p,rate=be_p))
}
