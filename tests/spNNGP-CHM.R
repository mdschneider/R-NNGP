rm(list=ls())
library(spNNGP)

rmvn <- function(n, mu=0, V = matrix(1)){
  p <- length(mu)
  if(any(is.na(match(dim(V),p))))
    stop("Dimension problem!")
  D <- chol(V)
  t(matrix(rnorm(n*p), ncol=p)%*%D + rep(mu,rep(n,p)))
}

##Make some data
set.seed(1)
n <- 20000
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

ho <- sample(1:n, 5000)

y.ho <- y[ho]
x.ho <- x[ho,,drop=FALSE]
w.ho <- w[ho]
coords.ho <- coords[ho,]

y <- y[-ho]
x <- x[-ho,,drop=FALSE]
w <- w[-ho,,drop=FALSE]
coords <- coords[-ho,]

##Fit a Response and Sequential NNGP model
n.samples <- 10000

starting <- list("phi"=phi, "sigma.sq"=5, "tau.sq"=1)

tuning <- list("phi"=0.1, "sigma.sq"=0.1, "tau.sq"=0.1)

priors <- list("phi.Unif"=c(3/1, 3/0.01), "sigma.sq.IG"=c(2, 5), "tau.sq.IG"=c(2, 1))

cov.model <- "exponential"

n.report <- 500
verbose <- TRUE

##sequential
set.seed(1)
m.s.b <- spNNGP(y~x-1, coords=coords, starting=starting, method="sequential", n.neighbors=15,
                tuning=tuning, priors=priors, cov.model=cov.model,
                n.samples=n.samples, search.type="brute", return.neighbors = TRUE, n.omp.threads=10, verbose=verbose, n.report=n.report)

set.seed(1)
m.s.t <- spNNGP(y~x-1, coords=coords, starting=starting, method="sequential", n.neighbors=15,
                tuning=tuning, priors=priors, cov.model=cov.model,
                n.samples=n.samples, search.type="tree", return.neighbors = TRUE, n.omp.threads=2, verbose=verbose, n.report=n.report)

identical(m.s.b$n.indx,m.s.t$n.indx)

##response
set.seed(1)
m.r.b <- spNNGP(y~x-1, coords=coords, starting=starting, method="response", n.neighbors=4,
                tuning=tuning, priors=priors, cov.model=cov.model,
                n.samples=n.samples, search.type="brute", n.omp.threads=2, verbose=verbose, n.report=n.report)

set.seed(1)
m.r.t <- spNNGP(y~x-1, coords=coords, starting=starting, method="response", n.neighbors=4,
                tuning=tuning, priors=priors, cov.model=cov.model,
                n.samples=n.samples, search.type="tree", n.omp.threads=2, verbose=verbose, n.report=n.report)

identical(m.s.b$n.indx,m.s.t$n.indx)

##conj
sigma.sq.IG <- c(2, sigma.sq)

g <- 10
theta.alpha <- cbind(seq(phi,30,length.out=g), seq(tau.sq/sigma.sq,5,length.out=g))

colnames(theta.alpha) <- c("phi", "alpha")

set.seed(1)
m.c.b <- spConjNNGP(y~x-1, coords=coords, n.neighbors = 15,
                    X.0 = x.ho, coords.0 = coords.ho,
                    k.fold = 5, score.rule = "crps",
                    n.omp.threads = 2, return.neighbors = TRUE, 
                    theta.alpha = theta.alpha, sigma.sq.IG = sigma.sq.IG, cov.model = cov.model, search.type="brute")

set.seed(1)
m.c.t <- spConjNNGP(y~x-1, coords=coords, n.neighbors = 15,
                    X.0 = x.ho, coords.0 = coords.ho,
                    k.fold = 5, score.rule = "crps",
                    n.omp.threads = 2, return.neighbors = TRUE, 
                    theta.alpha = theta.alpha, sigma.sq.IG = sigma.sq.IG, cov.model = cov.model, search.type="tree")

identical(m.c.b$n.indx, m.c.t$n.indx)

m.c.b$beta.hat
m.c.b$theta.alpha.sigmaSq
m.c.b$k.fold.scores

m.c.t$beta.hat
m.c.t$theta.alpha.sigmaSq
m.c.t$k.fold.scores

##n.indx
ord.coords <- m.c.t$coords

for(i in 1:n){
    plot(ord.coords)
    points(ord.coords[i,,drop=FALSE], pch=19)
    points(ord.coords[m.c.t$n.indx[[i]],,drop=FALSE],pch=19, col="blue")
    readline(prompt = "Pause. Press <Enter> to continue...")
}
