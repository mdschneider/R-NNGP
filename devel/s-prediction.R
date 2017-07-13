rm(list=ls())
library(yaImpute)
library(spBayes)
library(plyr)

y <- as.matrix(read.table(paste("../model1/large-syn-2/y.mod",sep="")))[,1]
X <- as.matrix(read.table(paste("../model1/large-syn-2/x.mod",sep="")))
coords <- as.matrix(read.table(paste("../model1/large-syn-2/coords.mod",sep="")))
w <- as.matrix(read.table(paste("../model1/large-syn-2/w.mod",sep="")))[,1]
n <- nrow(coords)

pred.X <- as.matrix(read.table(paste("../model1/large-syn-2/x.ho",sep="")))
pred.coords <- as.matrix(read.table(paste("../model1/large-syn-2/coords.ho",sep="")))
pred.y <- as.matrix(read.table(paste("../model1/large-syn-2/y.ho",sep="")))[,1]
pred.w <- as.matrix(read.table(paste("../model1/large-syn-2/w.ho",sep="")))[,1]
n.pred <- nrow(pred.coords)

m <- 15

samps <- cbind(t(matrix(scan(paste("nngp-phi-sig-blk/syn-large-chains/m",m,"-chain-1-beta",sep="")), nrow=2, byrow=TRUE)),
               t(matrix(scan(paste("nngp-phi-sig-blk/syn-large-chains/m",m,"-chain-1-theta",sep="")), nrow=3, byrow=TRUE)))

w.samps <- matrix(scan(paste("nngp-phi-sig-blk/syn-large-chains/m",m,"-chain-1-w",sep="")), nrow=n, byrow=TRUE)

colnames(samps) <- c("b0","b1","sigma.sq", "tau.sq", "phi")

samps <- samps[20001:nrow(samps),]
w.samps <- w.samps[,20001:nrow(samps)]

set.seed(1)
sub <- sample(1:nrow(samps), 100)
n.samples <- length(sub)

samps <- samps[sub,]
w.samps <- w.samps[,sub]


B <- samps[,c("b0","b1")]
sigma.sq <- samps[,"sigma.sq"]
tau.sq <- samps[,"tau.sq"]
phi <- samps[,"phi"]

y.hat <- matrix(0, n.pred, n.samples)
w.hat <- matrix(0, n.pred, n.samples)

nn <- ann(coords, pred.coords, k=m)$knnIndexDist
nn.indx <- nn[,1:m]
nn.dist <- sqrt(nn[,(m+1):ncol(nn)])

for(s in 1:n.samples){
   
    w.hat[,s] <- apply(nn, 1, function(i){

        nn.indx <- i[1:m]
        nn.dist <- sqrt(i[(m+1):ncol(nn)])

        c <- sigma.sq[s]*exp(-phi[s]*nn.dist)
        
        C <- sigma.sq[s]*exp(-phi[s]*iDist(coords[nn.indx,]))

        a <- t(c)%*%chol2inv(chol(C))
        B <- a%*%w.samps[nn.indx,s]
        F <- sigma.sq[s] - a%*%c
        rnorm(1, B, sqrt(F))
        
    })

    y.hat[,s] <- rnorm(nrow(pred.X), pred.X%*%B[s,] + w.hat[,s], sqrt(tau.sq[s]))
    
    print(s)
}


plot(pred.w, apply(w.hat, 1, mean))

plot(pred.y, apply(y.hat, 1, mean))

rmse <- function(a,b){
    sqrt(mean((a-b)^2))
}

rmse(pred.y, apply(y.hat, 1, mean))

write.table(y.hat, paste("samples-large-sym-m",m,"-pred",sep=""))
write.table(w.hat, paste("samples-large-sym-m",m,"-w-pred",sep=""))
