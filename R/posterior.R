#' posterior
#'
#' Posterior density function
#'
#' @param lambda_r numeric, true intensity.
#' @param k numeric, observed count for a taxa in a sample.
#' @param al_c numeric, estimated shape parameter for the contaminant intensity.
#' @param be_c numeric, estimated scale parameter for the contaminant intensity.
#'
#' @return numeric, posterior density function for the true intensity (up to a constant)
#' @export
posterior = function(lambda_r, k, al_c, be_c) {

    ga_mu = al_c/be_c

    if (k == 0) {
        pos_v = dgamma(lambda_r, shape = 0.5, rate = 1)

    } else if (k == 1) {
        pos_v = exp(-log(factorial(0.5)) + 0.5 * log(lambda_r + ga_mu) -
            (lambda_r + ga_mu))
    } else {
        c_m = numeric()
        c_m[1] = 0

        for (m in 2:k) {
            c_m[m] = log(m) + c_m[m - 1]
        }

        pos_v = exp(-c_m[(k - 1)] + (k - 0.5) * log(lambda_r + ga_mu) - (lambda_r +
            ga_mu))
    }

    return(pos_v)
}
