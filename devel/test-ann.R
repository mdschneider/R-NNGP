library(RANN)
library(yaImpute)

set.seed(1)

n <- 10
coords.1 <- cbind(runif(n,0,1), runif(n,0,1))
coords.2 <- cbind(runif(n,0,1), runif(n,0,1))

k <- 3
a <- ann(coords.1, coords.2, k=k)$knnIndexDist[,(k+1):(2*k)]
a

b <- nn2(coords.1, coords.2, k=k)$nn.dist
b

max(sqrt(a)-b)
