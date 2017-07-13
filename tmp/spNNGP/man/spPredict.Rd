\name{spPredict}
\alias{spPredict}
\title{Function for prediction at new locations using \code{spNNGP} models.}

\description{
  The function \code{spPredict} collects posterior predictive samples
  for a set of new locations given a \code{\link{spNNGP}} object.
}

\usage{
spPredict(sp.obj, X.0, coords.0, start=1, end, thin=1, n.omp.threads =
1, verbose=TRUE, n.report=100, ...)
}

\arguments{
  \item{sp.obj}{an object returned by \code{\link{spNNGP}}. }
  \item{X.0}{the design matrix for prediction locations. An
    intercept should be provided in the first column if one is specified
    in \code{sp.obj} model.}
  \item{coords.0}{the spatial coordinates corresponding to
    \code{X.0}. }
  \item{start}{specifies the first sample included in the composition sampling.}
  \item{end}{specifies the last sample included in the composition.
    The default is to use all posterior samples in \code{sp.obj}. }
  \item{thin}{a sample thinning factor.  The default of 1 considers all
    samples between \code{start} and \code{end}.  For example, if \code{thin = 10} then 1 in 10 samples are considered between \code{start} and
    \code{end}. }
  \item{n.omp.threads}{a positive integer indicating
    the number of threads to use for SMP parallel processing. The package must
    be compiled for OpenMP support. For most Intel-based machines, we recommend setting
    n.omp.threads to up to the number of hyperthreaded cores.}
  \item{verbose}{if \code{TRUE}, model specification and progress of the
    sampler is printed to the screen. Otherwise, nothing is printed to
    the screen.  }
  \item{n.report}{the interval to report sampling progress. }  
  \item{...}{currently no additional arguments.  }
}

\value{
  A list comprising:
  \item{p.y0}{a matrix that holds the response variable posterior
    predictive samples where rows are location and columns are samples. }
  \item{p.w0}{a matrix that holds the random effect posterior
    predictive samples where rows are location and columns are
  samples. This is only returned if \code{\link{spNNGP}} \code{method}
  is "sequential".}
\item{run.time}{execution time reported using \code{proc.time()}.}
  The return object will include additional data used for subsequent
  prediction and/or model fit evaluation.
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

}
}

\keyword{model}