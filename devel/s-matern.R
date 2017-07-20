rm(list=ls())
setwd("/home/andy/Rdevel/R-NNGP/devel")
source("spNNGP.R")
source("mkMatUtil.R")
dyn.load("sNNGP.so")
source("util.R")
library(coda)

rmvn <- function(n, mu=0, V = matrix(1)){

  p <- length(mu)

  if(any(is.na(match(dim(V),p))))

    stop("Dimension problem!")

  D <- chol(V)

  t(matrix(rnorm(n*p), ncol=p)%*%D + rep(mu,rep(n,p)))

}

##Make some data

set.seed(1)

n <- 5000

coords <- cbind(runif(n,0,1), runif(n,0,1))

x <- cbind(1, rnorm(n))

 

B <- as.matrix(c(1,5))

sigma.sq <- 5

tau.sq <- 1

phi <- 3/0.5

D <- as.matrix(dist(coords))

R <- exp(-phi*D)

w <- rmvn(1, rep(0,n), sigma.sq*R)

y <- rnorm(n, x%*%B + w, sqrt(tau.sq))

##Fit a Response and Sequential NNGP model

n.samples <- 5000

starting <- list("phi"=phi, "sigma.sq"=5, "tau.sq"=1)

tuning <- list("phi"=0, "sigma.sq"=0.1, "tau.sq"=0.1)

priors <- list("phi.Unif"=c(3/1, 3/0.01), "sigma.sq.IG"=c(2, 5), "tau.sq.IG"=c(2, 1))

starting$nu=0.5

priors$"nu.Unif"=c(0.1,4)

tuning$nu=0

dyn.load("sNNGP.so")
set.seed(1)
m.s.2 <- spNNGP(y~x-1, coords=coords, starting=starting, method="sequential", n.neighbors=10,

  tuning=tuning, priors=priors, cov.model="matern", search.type = "brute",

  n.samples=1, n.omp.threads=1, return.neighbors=TRUE)


## source("spCor.R")

matern <- function(d, phi, nu){
    (d*phi)^nu/(2^(nu-1)*gamma(nu))*besselK(x=d*phi, nu=nu)
}

d <- c(0.000308567, 0.00239029, 0.00415754, 0.00816737, 0.0100943, 0.0103185, 0.011832, 0.0119166, 0.0208336, 0.0260726)

phi <- 6
nu <- 1.9

c <- matern(d, phi, nu)

coords.ord <- coords[order(coords[,1]),]

n.indx <- m.s.2$n.indx

s <- 2949

library(spBayes)

d <- iDist(coords.ord[s,,drop=FALSE], coords.ord[n.indx[[s]],])
D <- iDist(coords.ord[n.indx[[s]],], coords.ord[n.indx[[s]],])

c <- matern(d, phi, nu)
C <- matern(D, phi, nu); diag(C) <- 1

c%*%chol2inv(chol(C))%*%t(c)


##formatC(c, digits = 8, format = "f")
