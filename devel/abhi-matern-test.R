rm(list=ls())
library(spNNGP)
library(spBayes)

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
n.samples <- 1000
starting <- list("phi"=phi, "sigma.sq"=5, "tau.sq"=1)
tuning <- list("phi"=0.1, "sigma.sq"=0.1, "tau.sq"=0.1)
priors <- list("phi.Unif"=c(3/1, 3/0.01), "sigma.sq.IG"=c(2, 5), "tau.sq.IG"=c(2, 1))

###Set nu to 1 and set tuning to 0 to get it to run so we can get the nearest neighbor index list
starting$nu=1
priors$"nu.Unif"=c(1,4)
tuning$nu=0

m.s.2 <- spNNGP(y~x-1, coords=coords, starting=starting, method="sequential", n.neighbors=10, tuning=tuning, priors=priors, cov.model="matern",
  n.samples=n.samples, n.omp.threads=2, return.neighbors=TRUE)


##Okay, things start going bad for nu=1.9 (and higher) for location 2949 in the update w loop.
s <- 2949
nu <- 1.9

matern <- function(d, phi, nu){
    (d*phi)^nu/(2^(nu-1)*gamma(nu))*besselK(x=d*phi, nu=nu)
}

coords.ord <- m.s.2$coords
n.indx <- m.s.2$n.indx

d <- iDist(coords.ord[s,,drop=FALSE], coords.ord[n.indx[[s]],])
D <- iDist(coords.ord[n.indx[[s]],], coords.ord[n.indx[[s]],])

c <- matern(d, phi, nu)
C <- matern(D, phi, nu); diag(C) <- 1

c%*%chol2inv(chol(C))%*%t(c)
