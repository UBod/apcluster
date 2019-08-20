\name{apcluster}
\docType{methods}
\alias{apcluster}
\alias{apcluster-methods}
\alias{apcluster,matrix,missing-method}
\alias{apcluster,dgTMatrix,missing-method}
\alias{apcluster,sparseMatrix,missing-method}
\alias{apcluster,Matrix,missing-method}
\alias{apcluster,character,ANY-method}
\alias{apcluster,function,ANY-method}
\title{Affinity Propagation}
\description{
   Runs affinity propagation clustering
}
\usage{
\S4method{apcluster}{matrix,missing}(s, x, p=NA, q=NA, maxits=1000,
    convits=100, lam=0.9, includeSim=FALSE, details=FALSE,
    nonoise=FALSE, seed=NA)
\S4method{apcluster}{dgTMatrix,missing}(s, x, p=NA, q=NA, maxits=1000,
    convits=100, lam=0.9, includeSim=FALSE, details=FALSE,
    nonoise=FALSE, seed=NA)
\S4method{apcluster}{sparseMatrix,missing}(s, x, ...)
\S4method{apcluster}{Matrix,missing}(s, x, ...)
\S4method{apcluster}{character,ANY}(s, x, p=NA, q=NA, maxits=1000,
    convits=100, lam=0.9, includeSim=TRUE, details=FALSE,
    nonoise=FALSE, seed=NA, ...)
\S4method{apcluster}{function,ANY}(s, x, p=NA, q=NA, maxits=1000,
    convits=100, lam=0.9, includeSim=TRUE, details=FALSE,
    nonoise=FALSE, seed=NA, ...)
}
\arguments{
  \item{s}{an \eqn{l\times l}{l x l} similarity matrix or a similarity
    function either specified as the name of a package-provided
    similarity function as character string or a user provided
    function object. \code{s} may also be a sparse matrix according to
    the \pkg{Matrix} package. Internally, \code{apcluster} uses the
    \code{\linkS4class{dgTMatrix}} class; all other sparse matrices are
    cast to this class (if possible, otherwise the function quits with
    an error). If \code{s} is any other object of class
    \code{\linkS4class{Matrix}}, \code{s} is cast to a regular matrix
    internally (if possible, otherwise the function quits with
    an error).} 
  \item{x}{input data to be clustered; if \code{x} is a matrix or data
    frame, rows are interpreted as samples and columns are 
    interpreted as features; apart from matrices or data frames, 
    \code{x} may be any other structured data type that
    contains multiple data items - provided that an appropriate
    \code{\link[base:length]{length}} function is available that
    returns the number of items}
  \item{p}{input preference; can be a vector that specifies
    individual preferences for each data point. If scalar,
    the same value is used for all data points. If \code{NA},
    exemplar preferences are initialized according to the
    distribution of non-Inf values in \code{s}. How this
    is done is controlled by the parameter \code{q}.}
  \item{q}{if \code{p=NA}, exemplar preferences are initialized
    according to the distribution of non-Inf values in \code{s}.
    If \code{q=NA}, exemplar preferences are set to the median
    of non-Inf values in \code{s}. If \code{q} is a value
    between 0 and 1, the sample quantile with threshold
    \code{q} is used, whereas \code{q=0.5} again results in
    the median.}
  \item{maxits}{maximal number of iterations that should be executed}
  \item{convits}{the algorithm terminates if the examplars have not
    changed for \code{convits} iterations}
  \item{lam}{damping factor; should be a value in the range [0.5, 1);
    higher values correspond to heavy damping which may be
    needed if oscillations occur}
  \item{includeSim}{if \code{TRUE}, the similarity matrix (either computed
    internally or passed via the \code{s} argument) is stored to the
    slot \code{sim} of the returned
    \code{\linkS4class{APResult}} object. The default is \code{FALSE}
    if \code{apcluster} has been called for a similarity matrix,
    otherwise the default is \code{TRUE}.}
  \item{details}{if \code{TRUE}, more detailed information about the
    algorithm's progress is stored in the output object
    (see \code{\linkS4class{APResult}})}
  \item{nonoise}{\code{apcluster} adds a small amount of noise to
    \code{s} to prevent degenerate cases; if \code{TRUE},
    this is disabled}
  \item{seed}{for reproducibility, the seed of the random number
    generator can be set to a fixed value before
    adding noise (see above), if \code{NA}, the seed remains
    unchanged}
  \item{...}{for the methods with signatures \code{character,ANY} and
    \code{function,ANY}, all other arguments are passed to the selected 
    similarity function as they are; for the methods with signatures
    \code{Matrix,missing} and \code{sparseMatrix,missing}, further
    arguments are passed on to the \code{apcluster} methods with
    signatures \code{Matrix,missing} and \code{dgTMatrix,missing},
    respectively.}
}
\details{Affinity Propagation clusters data using a set of
  real-valued pairwise data point similarities as input. Each cluster
  is represented by a cluster center data point (the so-called exemplar). 
  The method is iterative and searches for clusters maximizing
  an objective function called net similarity.

  When called with a similarity matrix as input (which may also be a
  sparse matrix according to the \pkg{Matrix} package), the function performs
  AP clustering. When called with the name of a package-provided
  similarity function or a user-provided similarity function object and 
  input data, the function first computes the similarity matrix before
  performing AP clustering. The similarity
  matrix is returned for later use as part of the
  \code{\linkS4class{APResult}}
  object depending on whether \code{includeSim} was set to \code{TRUE} (see
  argument description above).
 
  Apart from minor adaptations and optimizations, the AP
  clustering functionality of the function \code{apcluster} is 
  largely analogous to Frey's and Dueck's Matlab code
  (see \url{http://genes.toronto.edu/affinitypropagation/}).

  The new argument \code{q} allows for better controlling the number of
  clusters without knowing the distribution of similarity
  values. A meaningful range for the parameter \code{p} can be determined
  using the function \code{\link{preferenceRange}}. Alternatively, a
  certain fixed number of clusters may be desirable. For this purpose,
  the function \code{\link{apclusterK}} is available.
}
\value{
  Upon successful completion, the function returns an
  \code{\linkS4class{APResult}} object.
}
\author{Ulrich Bodenhofer, Andreas Kothmeier & Johannes Palme
\email{apcluster@bioinf.jku.at}}
\references{\url{http://www.bioinf.jku.at/software/apcluster}

Frey, B. J. and Dueck, D. (2007) Clustering by passing messages
between data points. \emph{Science} \bold{315}, 972-976.
DOI: \href{http://dx.doi.org/10.1126/science.1136800}{10.1126/science.1136800}.

Bodenhofer, U., Kothmeier, A., and Hochreiter, S. (2011)
APCluster: an R package for affinity propagation clustering.
\emph{Bioinformatics} \bold{27}, 2463-2464.
DOI: \href{http://dx.doi.org/10.1093/bioinformatics/btr406}{10.1093/bioinformatics/btr406}.
}
\seealso{\code{\link{APResult}}, \code{\link{show-methods}},
  \code{\link{plot-methods}}, \code{\link{labels-methods}},
  \code{\link{preferenceRange}}, \code{\link{apclusterL-methods}}, 
  \code{\link{apclusterK}}}
\examples{
## create two Gaussian clouds
cl1 <- cbind(rnorm(100, 0.2, 0.05), rnorm(100, 0.8, 0.06))
cl2 <- cbind(rnorm(50, 0.7, 0.08), rnorm(50, 0.3, 0.05))
x <- rbind(cl1, cl2)

## compute similarity matrix and run affinity propagation 
## (p defaults to median of similarity)
apres <- apcluster(negDistMat(r=2), x, details=TRUE)

## show details of clustering results
show(apres)

## plot clustering result
plot(apres, x)

## plot heatmap
heatmap(apres)

## run affinity propagation with default preference of 10% quantile
## of similarities; this should lead to a smaller number of clusters
## reuse similarity matrix from previous run
apres <- apcluster(s=apres@sim, q=0.1)
show(apres)
plot(apres, x)

## now try the same with RBF kernel
sim <- expSimMat(x, r=2)
apres <- apcluster(s=sim, q=0.2)
show(apres)
plot(apres, x)

## create sparse similarity matrix
cl1 <- cbind(rnorm(20, 0.2, 0.05), rnorm(20, 0.8, 0.06))
cl2 <- cbind(rnorm(20, 0.7, 0.08), rnorm(20, 0.3, 0.05))
x <- rbind(cl1, cl2)

sim <- negDistMat(x, r=2)
ssim <- as.SparseSimilarityMatrix(sim, lower=-0.2)

## run apcluster() on the sparse similarity matrix
apres <- apcluster(ssim, q=0)
apres
}
\keyword{cluster}
\keyword{methods}