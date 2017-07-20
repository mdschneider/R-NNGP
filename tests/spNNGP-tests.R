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
n <- 200
coords <- cbind(runif(n,0,1), runif(n,0,1))

x <- cbind(1, rnorm(n), rbinom(n, 1, 0.5))

B <- as.matrix(c(1,5,-1))

sigma.sq <- 5
tau.sq <- 1
phi <- 3/0.5

D <- as.matrix(dist(coords))
R <- exp(-phi*D)
w <- rmvn(1, rep(0,n), sigma.sq*R)
y <- rnorm(n, x%*%B + w, sqrt(tau.sq))

ho <- sample(1:n, 100)

y.ho <- y[ho]
x.ho <- x[ho,,drop=FALSE]
w.ho <- w[ho]
coords.ho <- coords[ho,]

y <- y[-ho]
x <- x[-ho,,drop=FALSE]
w <- w[-ho,,drop=FALSE]
coords <- coords[-ho,]

##conj 
sigma.sq.IG <- c(2, sigma.sq)
cov.model <- "exponential"
g <- 20
theta.alpha <- cbind(seq(phi,30,length.out=g), seq(tau.sq/sigma.sq,5,length.out=g))
colnames(theta.alpha) <- c("phi", "alpha")

set.seed(1)
m.c.b <- spConjNNGP(y~x-1, coords=coords, n.neighbors = 15,
                    X.0 = x.ho, coords.0 = coords.ho,
                    k.fold = 20, score.rule = "rmspe",
                    n.omp.threads = 2, return.neighbors = TRUE, 
                    theta.alpha = theta.alpha, sigma.sq.IG = sigma.sq.IG, cov.model = cov.model, search.type="brute")

set.seed(1)
m.c.t <- spConjNNGP(y~x-1, coords=coords, n.neighbors = 15,
                    X.0 = x.ho, coords.0 = coords.ho,
                    k.fold = 20, score.rule = "rmspe",
                    n.omp.threads = 2, return.neighbors = TRUE, 
                    theta.alpha = theta.alpha, sigma.sq.IG = sigma.sq.IG, cov.model = cov.model, search.type="tree")

for(i in 1:length(m.c.b)){
    if(!identical(m.c.b[[i]], m.c.t[[i]])){
        print(paste(names(m.c.b)[i], "not identical"))    
    }
}


##conj
cov.model <- "matern"
theta.alpha <- cbind(seq(phi,30,length.out=g), seq(0.1, 2, length.out=g), seq(tau.sq/sigma.sq,5,length.out=g))
colnames(theta.alpha) <- c("phi", "nu", "alpha")

set.seed(1)
m.c.b <- spConjNNGP(y~x-1, coords=coords, n.neighbors = 15,
                    X.0 = x.ho, coords.0 = coords.ho,
                    k.fold = 20, score.rule = "rmspe",
                    n.omp.threads = 2, return.neighbors = TRUE, 
                    theta.alpha = theta.alpha, sigma.sq.IG = sigma.sq.IG, cov.model = cov.model, search.type="brute")

set.seed(1)
m.c.t <- spConjNNGP(y~x-1, coords=coords, n.neighbors = 15,
                    X.0 = x.ho, coords.0 = coords.ho,
                    k.fold = 20, score.rule = "rmspe",
                    n.omp.threads = 2, return.neighbors = TRUE, 
                    theta.alpha = theta.alpha, sigma.sq.IG = sigma.sq.IG, cov.model = cov.model, search.type="tree")

for(i in 1:length(m.c.b)){
    if(!identical(m.c.b[[i]], m.c.t[[i]])){
        print(paste(names(m.c.b)[i], "not identical"))    
    }
}

##conj
cov.model <- "matern"
theta.alpha <- c(phi, 0.5, tau.sq/sigma.sq)
names(theta.alpha) <- c("phi", "nu", "alpha")

set.seed(1)
m.c.b <- spConjNNGP(y~x-1, coords=coords, n.neighbors = 15,
                    X.0 = x.ho, coords.0 = coords.ho,
                    k.fold = 20, score.rule = "rmspe",
                    n.omp.threads = 2, return.neighbors = TRUE, 
                    theta.alpha = theta.alpha, sigma.sq.IG = sigma.sq.IG, cov.model = cov.model, search.type="brute")

set.seed(1)
m.c.t <- spConjNNGP(y~x-1, coords=coords, n.neighbors = 15,
                    X.0 = x.ho, coords.0 = coords.ho,
                    k.fold = 20, score.rule = "rmspe",
                    n.omp.threads = 2, return.neighbors = TRUE, 
                    theta.alpha = theta.alpha, sigma.sq.IG = sigma.sq.IG, cov.model = cov.model, search.type="tree")

for(i in 1:length(m.c.b)){
    if(!identical(m.c.b[[i]], m.c.t[[i]])){
        print(paste(names(m.c.b)[i], "not identical"))    
    }
}

plot(m.c.b$y.0.hat, y.ho)

##spnngp

rmvn <- function(n, mu=0, V = matrix(1)){
  p <- length(mu)
  if(any(is.na(match(dim(V),p))))
    stop("Dimension problem!")
  D <- chol(V)
  t(matrix(rnorm(n*p), ncol=p)%*%D + rep(mu,rep(n,p)))
}

##Make some data
set.seed(1)
n <- 2000
coords <- cbind(runif(n,0,1), runif(n,0,1))

x <- cbind(1, rnorm(n), rbinom(n, 1, 0.5))

B <- as.matrix(c(1,5,-1))

sigma.sq <- 5
tau.sq <- 1
phi <- 3/0.5

D <- as.matrix(dist(coords))
R <- exp(-phi*D)
w <- rmvn(1, rep(0,n), sigma.sq*R)
y <- rnorm(n, x%*%B + w, sqrt(tau.sq))

ho <- sample(1:n, 1000)

y.ho <- y[ho]
x.ho <- x[ho,,drop=FALSE]
w.ho <- w[ho]
coords.ho <- coords[ho,]

y <- y[-ho]
x <- x[-ho,,drop=FALSE]
w <- w[-ho,,drop=FALSE]
coords <- coords[-ho,]

##Fit a Response and Sequential NNGP model
n.samples <- 1000

starting <- list("phi"=phi, "sigma.sq"=5, "tau.sq"=1)

tuning <- list("phi"=0.1, "sigma.sq"=0.1, "tau.sq"=0.1)

priors <- list("phi.Unif"=c(3/1, 3/0.01), "sigma.sq.IG"=c(2, 5), "tau.sq.IG"=c(2, 1))

cov.model <- "exponential"

n.report <- 500

##Predict for holdout set using both models
m.s <- spNNGP(y~x-1, coords=coords, starting=starting, method="sequential", n.neighbors=10,
              tuning=tuning, priors=priors, cov.model=cov.model,
              n.samples=n.samples, n.omp.threads=2, n.report=n.report, return.neighbors=TRUE)

m.r <- spNNGP(y~x-1, coords=coords, starting=starting, method="response", n.neighbors=10,
              tuning=tuning, priors=priors, cov.model=cov.model,
              n.samples=n.samples, n.omp.threads=2, n.report=n.report)

round(summary(m.s$p.beta.samples)$quantiles[,c(3,1,5)],2)
round(summary(m.s$p.theta.samples)$quantiles[,c(3,1,5)],2)
plot(apply(m.s$p.w.samples, 1, median), w)

round(summary(m.r$p.beta.samples)$quantiles[,c(3,1,5)],2)
round(summary(m.r$p.theta.samples)$quantiles[,c(3,1,5)],2)


##Prediction for holdout data
p.s <- spPredict(m.s, X.0 = x.ho, coords.0 = coords.ho, n.omp.threads=2)

plot(apply(p.s$p.w.0, 1, mean), w.ho)
plot(apply(p.s$p.y.0, 1, mean), y.ho)

p.r <- spPredict(m.r, X.0 = x.ho, coords.0 = coords.ho, n.omp.threads=2)

plot(apply(p.r$p.y.0, 1, mean), y.ho)

plot(apply(p.r$p.y.0, 1, mean), apply(p.s$p.y.0, 1, mean))

##threads
##response
m.r.1 <- spNNGP(y~x-1, coords=coords, starting=starting, method="response", n.neighbors=10,
                tuning=tuning, priors=priors, cov.model=cov.model,
                n.samples=n.samples, n.omp.threads=2, verbose=FALSE, n.report=n.report)

m.r.2 <- spNNGP(y~x-1, coords=coords, starting=starting, method="response", n.neighbors=10,
                tuning=tuning, priors=priors, cov.model=cov.model,
                n.samples=n.samples, n.omp.threads=2, verbose=FALSE, n.report=n.report)

m.r.1$run.time
m.r.2$run.time


m.s.1 <- spNNGP(y~x-1, coords=coords, starting=starting, method="sequential", n.neighbors=10,
                tuning=tuning, priors=priors, cov.model=cov.model,
                n.samples=n.samples, n.omp.threads=2, verbose=FALSE, n.report=n.report)

m.s.2 <- spNNGP(y~x-1, coords=coords, starting=starting, method="sequential", n.neighbors=10,
                tuning=tuning, priors=priors, cov.model=cov.model,
                n.samples=n.samples, n.omp.threads=2, verbose=FALSE, n.report=n.report)

m.s.1$run.time
m.s.2$run.time
