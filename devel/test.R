rm(list=ls())
dyn.load("rNNGP.so")
dyn.load("sNNGP.so")
source("spNNGP.R")
source("mkMatUtil.R")
library(Formula)
library(coda)
library(RANN)

rmvn <- function(n, mu=0, V = matrix(1)){
  p <- length(mu)
  if(any(is.na(match(dim(V),p))))
    stop("Dimension problem!")
  D <- chol(V)
  t(matrix(rnorm(n*p), ncol=p)%*%D + rep(mu,rep(n,p)))
}

set.seed(1)

n <- 2000
coords <- cbind(runif(n,0,1), runif(n,0,1))

x <- as.matrix(rep(1, n))

B <- as.matrix(c(-15))

sigma.sq <- 5
tau.sq <- 1
phi <- 3/0.7

D <- as.matrix(dist(coords))
R <- exp(-phi*D)
w <- rmvn(1, rep(0,n), sigma.sq*R)
y <- rnorm(n, x%*%B + w, sqrt(tau.sq))

set.seed(1)
ho <- sample(1:n, 1000)

y.ho <- y[ho]
x.ho <- x[ho,,drop=FALSE]
w.ho <- w[ho]
coords.ho <- coords[ho,]

y <- y[-ho]
x <- x[-ho,,drop=FALSE]
w <- w[-ho,,drop=FALSE]
coords <- coords[-ho,]

##########
n.samples <- 500

starting <- list("phi"=phi, "sigma.sq"=5, "tau.sq"=1)

tuning <- list("phi"=0.1, "sigma.sq"=0.1, "tau.sq"=0.1)

priors <- list("phi.Unif"=c(3/1, 3/0.01), "sigma.sq.IG"=c(2, 5), "tau.sq.IG"=c(2, 1))

cov.model <- "exponential"

n.report <- 500
verbose <- TRUE

m.1 <- spNNGP(y~x-1, coords=coords, starting=starting, method="sequential", n.neighbors=10,
              tuning=tuning, priors=priors, cov.model=cov.model,
              n.samples=n.samples, n.omp.threads=2, verbose=verbose, n.report=n.report)

round(summary(m.1$p.beta.samples)$quantiles[c(3,1,5)],2)
round(summary(m.1$p.theta.samples)$quantiles[,c(3,1,5)],2)

m.2 <- spNNGP(y~x-1, coords=coords, starting=starting, method="response", n.neighbors=10,
              tuning=tuning, priors=priors, cov.model=cov.model,
              n.samples=n.samples, n.omp.threads=2, verbose=verbose, n.report=n.report)

round(summary(m.2$p.beta.samples)$quantiles[c(3,1,5)],2)
round(summary(m.2$p.theta.samples)$quantiles[,c(3,1,5)],2)

source("spPredict.R")

dyn.load("sNNGPPredict.so")
set.seed(1)
c <- spPredict(m.1, X.0 = x.ho, coords.0 = coords.ho, n.omp.threads=4)

plot(apply(a$p.y0, 1, mean), y.ho)

dyn.load("sNNGPPredict.so")
set.seed(1)
b <- spPredict(m.1, X.0 = x.ho, coords.0 = coords.ho, n.omp.threads=1)

plot(apply(b$p.y0, 1, mean), y.ho)


dyn.load("rNNGPPredict.so")
b <- spPredict(m.2, X.0 = x.ho, coords.0 = coords.ho, n.omp.threads=1)

plot(y.ho, apply(b$p.y0, 1, mean))
