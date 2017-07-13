rm(list=ls())

library(spNNGP)

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
tau.sq <- 0.1
phi <- 3/0.5

D <- as.matrix(dist(coords))
R <- exp(-phi*D)
w <- rmvn(1, rep(0,n), sigma.sq*R)
y <- rnorm(n, x%*%B + w, sqrt(tau.sq))

sigma.sq.IG <- c(2, sigma.sq)

cov.model <- "matern"

#theta.alpha <- c(3/0.5, tau.sq/sigma.sq, 2)
#names(theta.alpha) <- c("phi", "alpha", "nu")

grd.res <- 5
theta.alpha <- as.matrix(expand.grid(seq(0.01/3,1/3,length.out=grd.res),
                                     seq(0.001,1,length.out=grd.res),
                                     seq(0.1,2,length.out=grd.res)))

colnames(theta.alpha) <- c("phi", "alpha", "nu")

m.1 <- spConjNNGP(y~x-1, coords=coords, n.neighbors = 10,
                  X.0 = x, coords.0 = coords,
                  k.fold = 5, score.rule = "crps",
                  n.omp.threads = 2,
                  theta.alpha = theta.alpha, sigma.sq.IG = sigma.sq.IG, cov.model = cov.model)

m.1$beta.hat
m.1$theta.alpha.sigmaSq
m.1$k.fold.scores

library(MBA)
library(fields)

z.lim <- range(c(s11[["z"]], s12[["z"]], s22[["z"]]))

s11 <- mba.surf(m.1$k.fold.scores[,c("alpha","phi","rmspe")], no.X=100, no.Y=100, extend=TRUE)$xyz.est

s12 <- mba.surf(m.1$k.fold.scores[,c("alpha","nu","rmspe")], no.X=100, no.Y=100, extend=TRUE)$xyz.est

s22 <- mba.surf(m.1$k.fold.scores[,c("phi","nu","rmspe")], no.X=100, no.Y=100, extend=TRUE)$xyz.est

layout(matrix(c(1,2,3,4), 2, 2, byrow = FALSE))

image(s11, xlab="alpha", ylab="phi", zlim=z.lim)
image(s12, xlab="alpha", ylab="nu", zlim=z.lim)
image.plot(legend.only=TRUE, zlim=z.lim) 
image(s22, xlab="phi", ylab="nu", zlim=z.lim)

