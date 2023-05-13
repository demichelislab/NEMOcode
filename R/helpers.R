#' Clamp values
#'
#' @param vec a vector to clamp
#' @param lower upper bound
#' @param upper lower bound
#'
#' @return a clamped numeric vector
#' @export
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
#' @export
.to_beta <- function(mval) {
    bb <- 2^mval / (2^mval + 1)
    return(bb)
}


#' Convert to MVal
#'
#' @param beta beta vector
#'
#' @return converted values in Mval scale
#' @export
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
#' @export
norm_one <- function(x){
    x/sum(x)
}


#' Beta sine transform
#'
#' @param beta numeric vector of values to transform
#'
#' @return transformed vector
#' @export
siTr = function(beta){
    return(asin((2*beta) - 1))
}


#' Inverse Beta-sine transform
#'
#' @param x numeric vetor of values to transform
#'
#' @return transformed vector
#' @export
invSi = function(x){
    return((sin(x)+1)/2)
}
