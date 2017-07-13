pkgname <- "spNNGP"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
library('spNNGP')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
cleanEx()
nameEx("spConjNNGP")
### * spConjNNGP

flush(stderr()); flush(stdout())

### Name: spConjNNGP
### Title: Function for fitting univariate Bayesian conjugate spatial
###   regression models
### Aliases: spConjNNGP
### Keywords: model

### ** Examples

## Not run: 
##D rmvn <- function(n, mu=0, V = matrix(1)){
##D   p <- length(mu)
##D   if(any(is.na(match(dim(V),p))))
##D     stop("Dimension problem!")
##D   D <- chol(V)
##D   t(matrix(rnorm(n*p), ncol=p)%*%D + rep(mu,rep(n,p)))
##D }
##D 
##D ##Make some data
##D set.seed(1)
##D n <- 2000
##D coords <- cbind(runif(n,0,1), runif(n,0,1))
##D 
##D x <- cbind(1, rnorm(n))
##D 
##D B <- as.matrix(c(1,5))
##D 
##D sigma.sq <- 5
##D tau.sq <- 1
##D phi <- 3/0.5
##D 
##D D <- as.matrix(dist(coords))
##D R <- exp(-phi*D)
##D w <- rmvn(1, rep(0,n), sigma.sq*R)
##D y <- rnorm(n, x%*%B + w, sqrt(tau.sq))
##D 
##D ho <- sample(1:n, 1000)
##D 
##D y.ho <- y[ho]
##D x.ho <- x[ho,,drop=FALSE]
##D w.ho <- w[ho]
##D coords.ho <- coords[ho,]
##D 
##D y <- y[-ho]
##D x <- x[-ho,,drop=FALSE]
##D w <- w[-ho,,drop=FALSE]
##D coords <- coords[-ho,]
##D 
##D ##Fit a Conjugate NNGP model and predict for the holdout
##D sigma.sq.IG <- c(2, sigma.sq)
##D 
##D cov.model <- "exponential"
##D 
##D g <- 10
##D theta.alpha <- cbind(seq(phi,30,length.out=g), seq(tau.sq/sigma.sq,5,length.out=g))
##D 
##D colnames(theta.alpha) <- c("phi", "alpha")
##D 
##D m.c <- spConjNNGP(y~x-1, coords=coords, n.neighbors = 10,
##D                   X.0 = x.ho, coords.0 = coords.ho,
##D                   k.fold = 5, score.rule = "crps",
##D                   n.omp.threads = 2,
##D                   theta.alpha = theta.alpha, sigma.sq.IG = sigma.sq.IG, cov.model = cov.model)
##D 
##D m.c$beta.hat
##D m.c$theta.alpha.sigmaSq
##D m.c$k.fold.scores
##D 
##D plot(m.c$y.0.hat, y.ho)
## End(Not run)



cleanEx()
nameEx("spNNGP")
### * spNNGP

flush(stderr()); flush(stdout())

### Name: spNNGP
### Title: Function for fitting univariate Bayesian spatial regression
###   models
### Aliases: spNNGP
### Keywords: model

### ** Examples

## Not run: 
##D 
##D rmvn <- function(n, mu=0, V = matrix(1)){
##D   p <- length(mu)
##D   if(any(is.na(match(dim(V),p))))
##D     stop("Dimension problem!")
##D   D <- chol(V)
##D   t(matrix(rnorm(n*p), ncol=p)%*%D + rep(mu,rep(n,p)))
##D }
##D 
##D ##Make some data
##D set.seed(1)
##D n <- 2000
##D coords <- cbind(runif(n,0,1), runif(n,0,1))
##D 
##D x <- cbind(1, rnorm(n))
##D 
##D B <- as.matrix(c(1,5))
##D 
##D sigma.sq <- 5
##D tau.sq <- 1
##D phi <- 3/0.5
##D 
##D D <- as.matrix(dist(coords))
##D R <- exp(-phi*D)
##D w <- rmvn(1, rep(0,n), sigma.sq*R)
##D y <- rnorm(n, x%*%B + w, sqrt(tau.sq))
##D 
##D ho <- sample(1:n, 1000)
##D 
##D y.ho <- y[ho]
##D x.ho <- x[ho,,drop=FALSE]
##D w.ho <- w[ho]
##D coords.ho <- coords[ho,]
##D 
##D y <- y[-ho]
##D x <- x[-ho,,drop=FALSE]
##D w <- w[-ho,,drop=FALSE]
##D coords <- coords[-ho,]
##D 
##D ##Fit a Response and Sequential NNGP model
##D n.samples <- 1000
##D 
##D starting <- list("phi"=phi, "sigma.sq"=5, "tau.sq"=1)
##D 
##D tuning <- list("phi"=0.1, "sigma.sq"=0.1, "tau.sq"=0.1)
##D 
##D priors <- list("phi.Unif"=c(3/1, 3/0.01), "sigma.sq.IG"=c(2, 5), "tau.sq.IG"=c(2, 1))
##D 
##D cov.model <- "exponential"
##D 
##D n.report <- 500
##D verbose <- TRUE
##D 
##D m.s <- spNNGP(y~x-1, coords=coords, starting=starting, method="sequential", n.neighbors=10,
##D               tuning=tuning, priors=priors, cov.model=cov.model,
##D               n.samples=n.samples, n.omp.threads=2, verbose=verbose, n.report=n.report)
##D 
##D round(summary(m.s$p.beta.samples)$quantiles[,c(3,1,5)],2)
##D round(summary(m.s$p.theta.samples)$quantiles[,c(3,1,5)],2)
##D 
##D m.r <- spNNGP(y~x-1, coords=coords, starting=starting, method="response", n.neighbors=10,
##D               tuning=tuning, priors=priors, cov.model=cov.model,
##D               n.samples=n.samples, n.omp.threads=2, verbose=verbose, n.report=n.report)
##D 
##D round(summary(m.r$p.beta.samples)$quantiles[,c(3,1,5)],2)
##D round(summary(m.r$p.theta.samples)$quantiles[,c(3,1,5)],2)
##D 
## End(Not run)



cleanEx()
nameEx("spPredict")
### * spPredict

flush(stderr()); flush(stdout())

### Name: spPredict
### Title: Function for prediction at new locations using 'spNNGP' models.
### Aliases: spPredict
### Keywords: model

### ** Examples

## Not run: 
##D 
##D rmvn <- function(n, mu=0, V = matrix(1)){
##D   p <- length(mu)
##D   if(any(is.na(match(dim(V),p))))
##D     stop("Dimension problem!")
##D   D <- chol(V)
##D   t(matrix(rnorm(n*p), ncol=p)%*%D + rep(mu,rep(n,p)))
##D }
##D 
##D ##Make some data
##D set.seed(1)
##D n <- 2000
##D coords <- cbind(runif(n,0,1), runif(n,0,1))
##D 
##D x <- cbind(1, rnorm(n))
##D 
##D B <- as.matrix(c(1,5))
##D 
##D sigma.sq <- 5
##D tau.sq <- 1
##D phi <- 3/0.5
##D 
##D D <- as.matrix(dist(coords))
##D R <- exp(-phi*D)
##D w <- rmvn(1, rep(0,n), sigma.sq*R)
##D y <- rnorm(n, x%*%B + w, sqrt(tau.sq))
##D 
##D ho <- sample(1:n, 1000)
##D 
##D y.ho <- y[ho]
##D x.ho <- x[ho,,drop=FALSE]
##D w.ho <- w[ho]
##D coords.ho <- coords[ho,]
##D 
##D y <- y[-ho]
##D x <- x[-ho,,drop=FALSE]
##D w <- w[-ho,,drop=FALSE]
##D coords <- coords[-ho,]
##D 
##D ##Fit a Response and Sequential NNGP model
##D n.samples <- 1000
##D 
##D starting <- list("phi"=phi, "sigma.sq"=5, "tau.sq"=1)
##D 
##D tuning <- list("phi"=0.1, "sigma.sq"=0.1, "tau.sq"=0.1)
##D 
##D priors <- list("phi.Unif"=c(3/1, 3/0.01), "sigma.sq.IG"=c(2, 5), "tau.sq.IG"=c(2, 1))
##D 
##D cov.model <- "exponential"
##D 
##D n.report <- 500
##D verbose <- TRUE
##D 
##D m.s <- spNNGP(y~x-1, coords=coords, starting=starting, method="sequential", n.neighbors=10,
##D               tuning=tuning, priors=priors, cov.model=cov.model,
##D               n.samples=n.samples, n.omp.threads=2, verbose=verbose, n.report=n.report)
##D 
##D round(summary(m.s$p.beta.samples)$quantiles[,c(3,1,5)],2)
##D round(summary(m.s$p.theta.samples)$quantiles[,c(3,1,5)],2)
##D 
##D m.r <- spNNGP(y~x-1, coords=coords, starting=starting, method="response", n.neighbors=10,
##D               tuning=tuning, priors=priors, cov.model=cov.model,
##D               n.samples=n.samples, n.omp.threads=2, verbose=verbose, n.report=n.report)
##D 
##D round(summary(m.r$p.beta.samples)$quantiles[,c(3,1,5)],2)
##D round(summary(m.r$p.theta.samples)$quantiles[,c(3,1,5)],2)
##D 
##D ##Prediction for holdout data
##D p.s <- spPredict(m.s, X.0 = x.ho, coords.0 = coords.ho, n.omp.threads=2)
##D 
##D plot(apply(p.s$p.y.0, 1, mean), y.ho)
##D plot(apply(p.s$p.w.0, 1, mean), w.ho)
##D 
##D p.r <- spPredict(m.r, X.0 = x.ho, coords.0 = coords.ho, n.omp.threads=2)
##D 
##D plot(apply(p.r$p.y.0, 1, mean), y.ho)
## End(Not run)



### * <FOOTER>
###
options(digits = 7L)
base::cat("Time elapsed: ", proc.time() - base::get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
