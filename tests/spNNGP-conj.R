rm(list=ls())
library(spNNGP)
library(geoR)

data(CHM)
colnames(CHM)
dim(CHM)

##remove zeros
CHM <- CHM[CHM[,3]>0,]

set.seed(1)
ho <- sample(1:nrow(CHM), 718498)

y <- CHM[,"CHM"]
coords <- CHM[,c("x","y")]
x <- coords

y.ho <- y[ho]
x.ho <- x[ho,,drop=FALSE]
coords.ho <- coords[ho,]

y.mod <- y[-ho]
x.mod <- x[-ho,,drop=FALSE]
coords.mod <- coords[-ho,]

##conj
#sub.samp <- sample(1:nrow(coords.mod), 10000)

#d.max <- max(as.matrix(dist(coords.mod[sub.samp,])))

#v <- variog(coords=coords.mod[sub.samp,],
#            data=resid(lm(y.mod[sub.samp] ~ x.mod[sub.samp,])),
#            uvec=(seq(0, 0.25*d.max, length=25)))
#plot(v)

##set up x-val
sigma.sq.IG <- c(2, 10)

#g <- 5
#theta.alpha <- as.matrix(expand.grid(seq(3/1000, 3/100, length.out=g), seq(1/10, 1, length.out=g)))
#colnames(theta.alpha) <- c("phi", "alpha")

theta.alpha <- c(0.03, 0.1)
names(theta.alpha) <- c("phi", "alpha")

set.seed(1)
m.c.b <- spConjNNGP(y.mod~x.mod, coords=coords.mod, n.neighbors = 10,
                    X.0 = cbind(1,x.ho), coords.0 = coords.ho,
                    k.fold = 10, score.rule = "crps",
                    n.omp.threads = 2, return.neighbors = TRUE, 
                    theta.alpha = theta.alpha, sigma.sq.IG = sigma.sq.IG, cov.model = "exponential", search.type="tree")

names(m.c.b)

m.c.b$theta.alpha
m.c.b$sigma.sq.hat
m.c.b$sigma.sq.var
m.c.b$beta.hat
m.c.b$beta.var
m.c.b$k.fold.scores
plot(y.ho, m.c.b$y.0.hat)

##
n.samples <- 100

starting <- list("phi"=0.3, "sigma.sq"=15, "tau.sq"=5)

tuning <- list("phi"=0.1, "sigma.sq"=0.1, "tau.sq"=0.1)

priors <- list("phi.Unif"=c(3/3000, 3/1), "sigma.sq.IG"=c(2, 15), "tau.sq.IG"=c(2, 5))

cov.model <- "exponential"

n.report <- 10

##Predict for holdout set using both models
m.s <- spNNGP(y.mod~x.mod, coords=coords.mod, starting=starting, method="sequential", n.neighbors=5,
              tuning=tuning, priors=priors, cov.model=cov.model,
              n.samples=n.samples, n.omp.threads=2, n.report=n.report)

m.r <- spNNGP(y.mod~x.mod, coords=coords.mod, starting=starting, method="response", n.neighbors=10,
              tuning=tuning, priors=priors, cov.model=cov.model,
              n.samples=n.samples, n.omp.threads=2, n.report=n.report)

round(summary(m.s$p.beta.samples)$quantiles[,c(3,1,5)],2)
round(summary(m.s$p.theta.samples)$quantiles[,c(3,1,5)],2)

summary(m.r$p.beta.samples)$quantiles[,c(3,1,5)]
round(summary(m.r$p.theta.samples)$quantiles[,c(3,1,5)],2)
