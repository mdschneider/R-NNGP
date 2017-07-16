pkgname <- "spNNGP"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
base::assign(".ExTimings", "spNNGP-Ex.timings", pos = 'CheckExEnv')
base::cat("name\tuser\tsystem\telapsed\n", file=base::get(".ExTimings", pos = 'CheckExEnv'))
base::assign(".format_ptime",
function(x) {
  if(!is.na(x[4L])) x[1L] <- x[1L] + x[4L]
  if(!is.na(x[5L])) x[2L] <- x[2L] + x[5L]
  options(OutDec = '.')
  format(x[1L:3L], digits = 7L)
},
pos = 'CheckExEnv')

### * </HEADER>
library('spNNGP')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
cleanEx()
nameEx("spConjNNGP")
### * spConjNNGP

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
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
##D ##one thread
##D m.c <- spConjNNGP(y~x-1, coords=coords, n.neighbors = 10,
##D                   X.0 = x.ho, coords.0 = coords.ho,
##D                   k.fold = 5, score.rule = "crps",
##D                   n.omp.threads = 1,
##D                   theta.alpha = theta.alpha, sigma.sq.IG = sigma.sq.IG, cov.model = cov.model)
##D 
##D m.c$beta.hat
##D m.c$theta.alpha.sigmaSq
##D m.c$k.fold.scores
##D 
##D ##two threads
##D m.c <- spConjNNGP(y~x-1, coords=coords, n.neighbors = 10,
##D                   X.0 = x.ho, coords.0 = coords.ho,
##D                   k.fold = 5, score.rule = "crps",
##D                   n.omp.threads = 2,
##D                   theta.alpha = theta.alpha, sigma.sq.IG = sigma.sq.IG, cov.model = cov.model)
##D 
##D m.c$beta.hat
##D m.c$theta.alpha.sigmaSq
##D m.c$k.fold.scores
## End(Not run)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("spConjNNGP", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("spNNGP")
### * spNNGP

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: spNNGP
### Title: Function for fitting univariate Bayesian spatial regression
###   models
### Aliases: spNNGP
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
##D n <- 1000
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
##D ##Fit a Response and Sequential NNGP model
##D n.samples <- 500
##D 
##D starting <- list("phi"=phi, "sigma.sq"=5, "tau.sq"=1)
##D 
##D tuning <- list("phi"=0.1, "sigma.sq"=0.1, "tau.sq"=0.1)
##D 
##D priors <- list("phi.Unif"=c(3/1, 3/0.01), "sigma.sq.IG"=c(2, 5), "tau.sq.IG"=c(2, 1))
##D 
##D cov.model <- "exponential"
##D 
##D m.s <- spNNGP(y~x-1, coords=coords, starting=starting, method="sequential", n.neighbors=10,
##D               tuning=tuning, priors=priors, cov.model=cov.model,
##D               n.samples=n.samples, n.omp.threads=2)
##D 
##D round(summary(m.s$p.beta.samples)$quantiles[,c(3,1,5)],2)
##D round(summary(m.s$p.theta.samples)$quantiles[,c(3,1,5)],2)
##D 
##D m.r <- spNNGP(y~x-1, coords=coords, starting=starting, method="response", n.neighbors=10,
##D               tuning=tuning, priors=priors, cov.model=cov.model,
##D               n.samples=n.samples, n.omp.threads=2)
##D 
##D round(summary(m.r$p.beta.samples)$quantiles[,c(3,1,5)],2)
##D round(summary(m.r$p.theta.samples)$quantiles[,c(3,1,5)],2)
## End(Not run)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("spNNGP", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("spPredict")
### * spPredict

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: spPredict
### Title: Function for prediction at new locations using 'spNNGP' models.
### Aliases: spPredict
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
##D n <- 1000
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
##D ho <- sample(1:n, 500)
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
##D n.samples <- 500
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
##D 
##D ##Predict for holdout set using both models
##D m.s <- spNNGP(y~x-1, coords=coords, starting=starting, method="sequential", n.neighbors=10,
##D               tuning=tuning, priors=priors, cov.model=cov.model,
##D               n.samples=n.samples, n.omp.threads=2, n.report=n.report)
##D 
##D m.r <- spNNGP(y~x-1, coords=coords, starting=starting, method="response", n.neighbors=10,
##D               tuning=tuning, priors=priors, cov.model=cov.model,
##D               n.samples=n.samples, n.omp.threads=2, n.report=n.report)
##D 
##D ##Prediction for holdout data
##D p.s <- spPredict(m.s, X.0 = x.ho, coords.0 = coords.ho, n.omp.threads=2)
##D 
##D plot(apply(p.s$p.w.0, 1, mean), w.ho)
##D plot(apply(p.s$p.y.0, 1, mean), y.ho)
##D 
##D p.r <- spPredict(m.r, X.0 = x.ho, coords.0 = coords.ho, n.omp.threads=2)
##D 
##D plot(apply(p.r$p.y.0, 1, mean), y.ho)
## End(Not run)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("spPredict", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
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
