rm(list=ls())
setwd("/home/andy/Rdevel/R-NNGP/devel")
dyn.load("cNNGP.so")
source("spConjNNGP.R")
source("mkMatUtil.R")
source("util.R")
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

n <- 10000
coords <- cbind(runif(n,0,1), runif(n,0,1))

x <- as.matrix(cbind(1, rnorm(n)))

B <- as.matrix(c(1,5))

sigma.sq <- 5
tau.sq <- 0.0001
phi <- 3/0.5

D <- as.matrix(dist(coords))
R <- exp(-phi*D)
w <- rmvn(1, rep(0,n), sigma.sq*R)
y <- rnorm(n, x%*%B + w, sqrt(tau.sq))

##########
sigma.sq.IG <- c(2, sigma.sq)

cov.model <- "matern"#"exponential"

n.report <- 25000
verbose <- TRUE

#theta.alpha <- c(3/0.5, tau.sq/sigma.sq, 2)
#names(theta.alpha) <- c("phi", "alpha", "nu")
g <- 10
theta.alpha <- cbind(seq(phi,30,length.out=g), seq(tau.sq/sigma.sq,5,length.out=g), seq(0.5,2,length.out=g))
colnames(theta.alpha) <- c("phi", "alpha", "nu")

source("spConjNNGP.R")
m.1 <- spConjNNGP(y~x-1, coords=coords, n.neighbors = 10,
                  X.0 = x, coords.0 = coords,
                  k.fold = 5, score.rule = "rmspe",
                  n.omp.threads = 4,
                  theta.alpha = theta.alpha, sigma.sq.IG = sigma.sq.IG, cov.model = cov.model)

m.1$beta.hat
m.1$theta.alpha.sigmaSq
m.1$k.fold.scores
