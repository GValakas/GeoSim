#' Transform a Variable into Normal Scores
#'
#' This function transforms a variable into normal scores.
#'
#' @param values A vector of values to be transformed.
#' @details The z-values weights are taken as equal.
#' Ensure that the provided data values follow a normal distribution for accurate transformation into normal scores.
#'
#' @export

normal.scores <- function(values){
#library(pracma)
y <- values
z <- unique(y)
n <- length(z)
w <- c(0,rep(1,n) / n) # weights are taken as equal

suppressWarnings(cdfzfile <- cumsum(0.5 * w[1:n] + 0.5 * w[2:(n + 1)]))
z <- sort(z)
cdfyinterp <- approx(z, cdfzfile, xout = y, method = "linear")
cdfyinterp <- cdfyinterp$y
#cdfyinterp <- interp1(z,cdfzfile,y,"linear");

mean_par <- 0
var_par <- 1
# zgaussian <- sqrt(2) * erfinv(2 * cdfyinterp-1) # or 
zgaussian <- sqrt(2) * qnorm((1 + 2 * cdfyinterp - 1) / 2) / sqrt(2)
index_inf <- which(is.infinite(zgaussian))
if (length(index_inf) > 0){
for (i in 1:length(index_inf)){
if (zgaussian[index_inf[i]] > 0){
zgaussian[index_inf[i]] <- 8.3
}
else{
zgaussian[index_inf[i]] <- -8.3
} 
} 
}

zgaussian <- zgaussian * sqrt(var_par) + mean_par
return(zgaussian)
}

