rm(list=ls())
library(yaImpute)
library(spBayes)
library(plyr)

y <- as.matrix(read.table(paste("large-syn-2/y.mod",sep="")))[,1]
X <- as.matrix(read.table(paste("large-syn-2/x.mod",sep="")))
coords <- as.matrix(read.table(paste("large-syn-2/coords.mod",sep="")))
w <- as.matrix(read.table(paste("large-syn-2/w.mod",sep="")))[,1]

pred.X <- as.matrix(read.table(paste("large-syn-2/x.ho",sep="")))
pred.coords <- as.matrix(read.table(paste("large-syn-2/coords.ho",sep="")))
pred.y <- as.matrix(read.table(paste("large-syn-2/y.ho",sep="")))[,1]
pred.w <- as.matrix(read.table(paste("large-syn-2/w.ho",sep="")))[,1]
n.pred <- nrow(pred.coords)
n <- nrow(pred.coords)

m <- 15

n.samples <- 25000
set.seed(1)
sub <- sample(floor(0.9*n.samples):n.samples, 100)

samps <- cbind(t(matrix(scan(paste("syn-large-chains/m",m,"-chain-2-beta",sep="")), nrow=2, byrow=TRUE)),
               t(matrix(scan(paste("syn-large-chains/m",m,"-chain-2-theta",sep="")), nrow=3, byrow=TRUE)))[sub,]

colnames(samps) <- c("b0","b1","sigma.sq", "tau.sq", "phi")

B <- samps[,c("b0","b1")]
sigma.sq <- samps[,"sigma.sq"]
tau.sq <- samps[,"tau.sq"]
phi <- samps[,"phi"]

n.samples <- nrow(samps)

y.samps <- matrix(0, n.pred, n.samples)

nn <- ann(coords, pred.coords, k=m)$knnIndexDist
nn.indx <- nn[,1:m]
#nn.dist <- sqrt(nn[,(m+1):ncol(nn)])

for(s in 1:n.samples){
   
    y.samps[,s] <- apply(cbind(1:n,nn), 1, function(i){

        loc <- i[1]
        nn.indx <- i[2:(m+1)]
        nn.dist <- sqrt(i[(m+2):length(i)])

        c <- sigma.sq[s]*exp(-phi[s]*nn.dist)
        
        C <- sigma.sq[s]*exp(-phi[s]*iDist(coords[nn.indx,])) + diag(tau.sq[s],m)

        a <- t(c)%*%chol2inv(chol(C))
        v <- a%*%(y[nn.indx] - X[nn.indx,]%*%B[s,])
        F <- sigma.sq[s] + tau.sq[s] - a%*%c

        rnorm(1, pred.X[loc,]%*%B[s,]+v, sqrt(F))
        
    })

    
    print(s)
}

plot(pred.y, apply(y.samps, 1, mean))

rmse <- function(a,b){
    sqrt(mean((a-b)^2))
}

rmse(pred.y, apply(y.samps, 1, mean))

write.table(y.samps, paste("samples-large-sym-m",m,"-pred",sep=""))
