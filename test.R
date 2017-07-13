library(spNNGP)

coords.mod <- as.matrix(read.table("coords.mod", header=F))
y.mod <- as.matrix(read.table("y.mod", header=F)[,1])
x.mod <- as.matrix(read.table("x.coords.mod", header=F))

##Fit a Response and Sequential NNGP model
n.samples <- 1000

starting <- list("phi"=10, "sigma.sq"=5, "tau.sq"=1)

tuning <- list("phi"=0.1, "sigma.sq"=0.1, "tau.sq"=0.1)

priors <- list("phi.Unif"=c(0.6, 30), "sigma.sq.IG"=c(2, 5), "tau.sq.IG"=c(2, 0.001))

cov.model <- "exponential"

n.report <- 500
verbose <- TRUE

m.s <- spNNGP(y.mod~x.mod-1, coords=coords.mod, starting=starting, method="sequential", n.neighbors=10,
              tuning=tuning, priors=priors, cov.model=cov.model,
              n.samples=n.samples, n.omp.threads=2, verbose=verbose, n.report=n.report)

round(summary(m.s$p.beta.samples)$quantiles[,c(3,1,5)],2)
round(summary(m.s$p.theta.samples)$quantiles[,c(3,1,5)],2)

m.r <- spNNGP(y~x-1, coords=coords, starting=starting, method="response", n.neighbors=10,
              tuning=tuning, priors=priors, cov.model=cov.model,
              n.samples=n.samples, n.omp.threads=2, verbose=verbose, n.report=n.report)

round(summary(m.r$p.beta.samples)$quantiles[,c(3,1,5)],2)
round(summary(m.r$p.theta.samples)$quantiles[,c(3,1,5)],2)

##Prediction for holdout data
p.s <- spPredict(m.s, X.0 = x.ho, coords.0 = coords.ho, n.omp.threads=2)

plot(apply(p.s$p.y.0, 1, mean), y.ho)
plot(apply(p.s$p.w.0, 1, mean), w.ho)

p.r <- spPredict(m.r, X.0 = x.ho, coords.0 = coords.ho, n.omp.threads=2)

plot(apply(p.r$p.y.0, 1, mean), y.ho)
