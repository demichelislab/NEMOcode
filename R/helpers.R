

#' Probit transformation
#'
#' @param x beta values
#' @param multiplier the number that multiplies -x at the exponent (default is 1)
#'
#' @return probit transformed value
#' @export
.probit_trans = function(x, multiplier = 1) {
    res = 1 / (1 + exp(-x * multiplier))
    return(res)
}



#' Clamp values
#'
#' @param vec a vector to clamp
#' @param lower upper bound
#' @param upper lower bound
#'
#' @return a clamped numeric vector
.clamp_val = function(vec, lower, upper){
    cla = vec
    cla[which(cla < lower)] = lower
    cla[which(cla > upper)] = upper
    return(cla)
}



#' Convert to Beta
#'
#' @param mval mval vector
#'
#' @return converted values in Beta scale
.to_beta <- function(mval) {
    bb <- 2^mval / (2^mval + 1)
    return(bb)
}


#' Convert to MVal
#'
#' @param beta beta vector
#'
#' @return converted values in Mval scale
.to_mval <- function(beta) {
    ## Extreme values should be converted to avoid errors
    beta[which(beta < 0.01)] <- 0.01
    beta[which(beta > 0.99)] <- .99
    mm <- log2(beta / (1 - beta))
    return(mm)
}



#' Normalize to 1
#'
#' @param x numeric vector
#'
#' @return a normalized numeric vector
norm_one <- function(x){
    x/sum(x)
}
