crps <- function(y, y.hat, y.var){
    sd <- sqrt(y.var)
    y.std <- (y-y.hat)/sd
    mean(-sd*(1/sqrt(pi) - 2*dnorm(y.std) - y.std*(2*pnorm(y.std) - 1)))
}

rmspe <- function(y, y.hat){
    sqrt(mean((y-y.hat)^2))
}

