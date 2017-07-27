library(spNNGP)

rmvn <- function(n, mu=0, V = matrix(1)){
  p <- length(mu)
  if(any(is.na(match(dim(V),p))))
    stop("Dimension problem!")
  D <- chol(V)
  t(matrix(rnorm(n*p), ncol=p)%*%D + rep(mu,rep(n,p)))
}

##Make some data
set.seed(2)
n <- 50
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
n.samples <- 500

starting <- list("phi"=phi, "sigma.sq"=5, "tau.sq"=1)

tuning <- list("phi"=0.1, "sigma.sq"=0.1, "tau.sq"=0.1)

priors <- list("phi.Unif"=c(3/1, 3/0.01), "sigma.sq.IG"=c(2, 5), "tau.sq.IG"=c(2, 1))

cov.model <- "exponential"

m.s <- spNNGP(y~x-1, coords=coords, starting=starting, method="sequential", n.neighbors=5,
              tuning=tuning, priors=priors, cov.model=cov.model,
              n.samples=n.samples, n.omp.threads=2, return.neighbors=TRUE)
pdf(file="neighbs-for-location-18.pdf")
i <- 18
plot(m.s$coords.ord, cex=2, xlab="Easting", ylab="Northing")
abline(v=m.s$coords.ord[i,1,drop=FALSE], lty=3, lwd=2)
points(m.s$coords.ord[i,,drop=FALSE], col="blue", pch=19, cex=2)
points(m.s$coords.ord[m.s$n.indx[[i]],,drop=FALSE], col="red", pch=19, cex=2)
dev.off()

pdf(file="neighbs-for-location-46.pdf")
i <- 46
plot(m.s$coords.ord, cex=2, xlab="Easting", ylab="Northing")
abline(v=m.s$coords.ord[i,1,drop=FALSE], lty=3, lwd=2)
points(m.s$coords.ord[i,,drop=FALSE], col="blue", pch=19, cex=2)
points(m.s$coords.ord[m.s$n.indx[[i]],,drop=FALSE], col="red", pch=19, cex=2)
dev.off()
