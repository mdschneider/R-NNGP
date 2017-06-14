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
nameEx("rNNGP")
### * rNNGP

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: rNNGP
### Title: Function for fitting univariate Bayesian spatial regression
###   models
### Aliases: rNNGP
### Keywords: model

### ** Examples

## Not run: 
##D 
##D rmvn <- function(n, mu=0, V = matrix(1)){
##D   p <- length(mu)
##D   if(any(is.na(match(dim(V),p))))
##D     stop("Dimension problem!")
##D   D <- chol(V)
##D   t(matrix(rnorm(n*p), ncol=p)##D 
##D }
##D 
##D set.seed(1)
##D 
##D n <- 5000
##D coords <- cbind(runif(n,0,1), runif(n,0,1))
##D X <- as.matrix(cbind(1, rnorm(n)))
##D 
##D B <- as.matrix(c(1,5))
##D p <- length(B)
##D 
##D sigma.sq <- 2
##D tau.sq <- 1
##D phi <- 3/0.5
##D 
##D D <- as.matrix(dist(coords))
##D R <- exp(-phi*D)
##D w <- rmvn(1, rep(0,n), sigma.sq*R)
##D y <- rnorm(n, X##D 
##D 
##D n.samples <- 2000
##D 
##D starting <- list("phi"=3/0.5, "sigma.sq"=5, "tau.sq"=1)
##D 
##D tuning <- list("phi"=0.1, "sigma.sq"=0.1, "tau.sq"=0.1)
##D 
##D priors <- list("phi.Unif"=c(3/1, 3/0.01), "sigma.sq.IG"=c(2, 1), "tau.sq.IG"=c(2, 1))
##D 
##D cov.model <- "exponential"
##D 
##D n.report <- 500
##D verbose <- TRUE
##D 
##D m.1 <- rNNGP(y~X-1, coords=coords, starting=starting, n.neighbors=10,
##D               tuning=tuning, priors=priors, cov.model=cov.model,
##D               n.samples=n.samples, n.omp.threads=2, verbose=verbose, n.report=n.report)
##D 
##D plot(mcmc(m.1$p.beta.samples), density=FALSE)
##D 
##D plot(mcmc(m.1$p.theta.samples), density=FALSE)
##D 
## End(Not run)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("rNNGP", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
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
