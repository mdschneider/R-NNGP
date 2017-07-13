\name{spConjNNGP}
\alias{spConjNNGP}
\title{Function for fitting univariate Bayesian conjugate spatial regression
  models}

\description{
  The function \code{spConjNNGP} fits Gaussian univariate Bayesian conjugate spatial
  regression models using Nearest Neighbor Gaussian Processes (NNGP).
}

\usage{
    spConjNNGP(formula, data = parent.frame(), coords, n.neighbors = 15,
               theta.alpha, sigma.sq.IG, cov.model = "exponential",
               k.fold, score.rule,
               X.0, coords.0, 
               n.omp.threads = 1, verbose=TRUE, ...)
}

\arguments{
  \item{formula}{a symbolic description of the regression model to be
    fit. See example below. }
  \item{data}{an optional data frame containing the variables in the
    model. If not found in data, the variables are taken from
    \code{environment(formula)}, typically the environment from which \code{spConjNNGP} is called.  }
  \item{coords}{an \eqn{n \times 2}{n x 2} matrix of the observation coordinates
    in \eqn{R^2}{R^2} (e.g., easting and northing). }
  \item{n.neighbors}{number of neighbors used in the NNGP.}
  \item{theta.alpha}{a vector or matrix of parameter values for \code{phi}, \code{nu}, and \code{alpha}, where \eqn{\alpha=\tau^2/sigma^2}{alpha=tau^2/sigma^2}\code{nu} is only required if \code{cov.model} equals \code{"matern"}. A vector is passed if you want to run the model using one set of parameters. The vector elements must be named and hold values for \code{phi}, \code{nu}, and \code{alpha}. If a matrix is passed, columns must be named and hold values for \code{phi}, \code{alpha}, and \code{nu}. Each row in the matrix defines each set of parameters for which the the model will be run. }
  \item{sigma.sq.IG}{a vector of length two that holds the hyperparameters, \emph{shape} and \emph{scale} respectively, for the inverse-Gamma prior on \eqn{\sigma^2}{sigma^2}. }
  \item{cov.model}{a quoted keyword that specifies the covariance
    function used to model the spatial dependence structure among the
    observations.  Supported covariance model key words are:
    \code{"exponential"}, \code{"matern"}, \code{"spherical"}, and
    \code{"gaussian"}.  See below for details.}
  \item{k.fold}{an optional argument used to specify the number of \emph{k} folds for cross-validation. In k-fold cross-validation, the data specified in \code{model} is randomly partitioned into \emph{k} equal sized subsamples. Of the \emph{k} subsamples, \emph{k}-1 subsamples are used to fit the model and the remaining \emph{k} samples are used for prediction. The cross-validation process is repeated \emph{k} times (the folds). Root mean squared prediction error (RMSPE) and continuous ranked probability score (CRPS; Gneiting and Raftery, 2007) rules are averaged over the \emph{k} fold prediction results and reported for the parameter set(s) defined by \code{theta.alpha}. The parameter set that yields the \emph{best} performance based on the scoring rule defined by \code{score.rule} is used fit the final model that uses all the data and make predictions if \code{X.0} and \code{coords.0} are specified.  Results from the k-fold cross-validation are returned in the \code{k.fold.scores} matrix. }
      \item{score.rule}{a quoted keyword "rmspe" or "crps" that specifies the scoring rule used to select the \emph{best} parameter set, see argument definition for \code{k.fold} for more details.}
      \item{X.0}{the design matrix for prediction locations. An
  intercept should be provided in the first column if one is specified in \code{model}.}
      \item{coords.0}{the spatial coordinates corresponding to \code{X.0}. }
  \item{n.omp.threads}{a positive integer indicating
the number of threads to use for SMP parallel processing. The package must
be compiled for OpenMP support. For most Intel-based machines, we recommend setting
n.omp.threads to up to the number of hyperthreaded cores.}
  \item{verbose}{if \code{TRUE}, model specification and progress of the
    sampler is printed to the screen. Otherwise, nothing is printed to
    the screen.}
  \item{...}{currently no additional arguments.  }
}

\value{
  An object of class \code{cNNGP}, which is a list comprising:
  \item{beta.hat}{a matrix of regression coefficient estimates corresponding to parameter set(s) defined in \code{theta.alpha}. }
  \item{theta.alpha.sigmaSq}{the \code{theta.alpha} matrix with
    \eqn{\sigma^2}{sigma^2} estimates appended. }
\item{run.time}{execution time for building the nearest neighbor index
  and parameter estimation reported using \code{proc.time()}.}
  The return object will include additional data used for subsequent
  prediction and/or model fit evaluation.
}

\references{
  Datta, A., S. Banerjee, A.O. Finley, and A.E. Gelfand. (2016) Hierarchical Nearest-Neighbor Gaussian process models for large geostatistical datasets. Journal of the American Statistical Association, 111:800-812.

  Gneiting, T and A.E. Raftery. (2007) Strictly proper scoring rules, prediction, and estimation. Journal of the American Statistical Association, 102:359-378.
}

\author{
  Andrew O. Finley \email{finleya@msu.edu}, \cr
  Abhirup Datta \email{abhidatta@jhu.edu}, \cr
  Sudipto Banerjee \email{sudipto@ucla.edu}
}

\examples{
\dontrun{
rmvn <- function(n, mu=0, V = matrix(1)){
  p <- length(mu)
  if(any(is.na(match(dim(V),p))))
    stop("Dimension problem!")
  D <- chol(V)
  t(matrix(rnorm(n*p), ncol=p)%*%D + rep(mu,rep(n,p)))
}

set.seed(1)

n <- 1000
coords <- cbind(runif(n,0,1), runif(n,0,1))

x <- as.matrix(cbind(1, rnorm(n)))

B <- as.matrix(c(1,5))

sigma.sq <- 5
tau.sq <- 0.1
phi <- 3/0.5

D <- as.matrix(dist(coords))
R <- exp(-phi*D)
w <- rmvn(1, rep(0,n), sigma.sq*R)
y <- rnorm(n, x\%*\%B + w, sqrt(tau.sq))

sigma.sq.IG <- c(2, sigma.sq)

cov.model <- "matern"

#theta.alpha <- c(3/0.5, tau.sq/sigma.sq, 2)
#names(theta.alpha) <- c("phi", "alpha", "nu")

n.test.sets <- 10
theta.alpha <- cbind(seq(phi,30,length.out=n.test.sets),
                     seq(tau.sq/sigma.sq,5,length.out=n.test.sets),
                     seq(0.5,2,length.out=n.test.sets))

colnames(theta.alpha) <- c("phi", "alpha", "nu")

m.1 <- spConjNNGP(y~x-1, coords=coords, n.neighbors = 10,
                  X.0 = x, coords.0 = coords,
                  k.fold = 5, score.rule = "crps",
                  n.omp.threads = 2,
                  theta.alpha = theta.alpha, sigma.sq.IG = sigma.sq.IG, cov.model = cov.model)

m.1$beta.hat
m.1$theta.alpha.sigmaSq
m.1$k.fold.scores

}
}

\keyword{model}
