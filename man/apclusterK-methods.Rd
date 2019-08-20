\name{apclusterK}
\docType{methods}
\alias{apclusterK}
\alias{apclusterK-methods}
\alias{apclusterK,matrix,missing-method}
\alias{apclusterK,Matrix,missing-method}
\alias{apclusterK,dgTMatrix,missing-method}
\alias{apclusterK,sparseMatrix,missing-method}
\alias{apclusterK,character,ANY-method}
\alias{apclusterK,function,ANY-method}
\title{Affinity Propagation for Pre-defined Number of Clusters}
\description{
  Runs affinity propagation clustering for a given similarity matrix
  adjusting input preferences iteratively in order to achieve a desired
  number of clusters
}
\usage{
\S4method{apclusterK}{matrix,missing}(s, x, K, prc=10, bimaxit=20, exact=FALSE,
     maxits=1000, convits=100, lam=0.9, includeSim=FALSE, details=FALSE,
     nonoise=FALSE, seed=NA, verbose=TRUE)
\S4method{apclusterK}{Matrix,missing}(s, x, K, ...)
\S4method{apclusterK}{dgTMatrix,missing}(s, x, K, prc=10, bimaxit=20,
     exact=FALSE, maxits=1000, convits=100, lam=0.9, includeSim=FALSE,
     details=FALSE, nonoise=FALSE, seed=NA, verbose=TRUE)
\S4method{apclusterK}{sparseMatrix,missing}(s, x, K, ...)
\S4method{apclusterK}{function,ANY}(s, x, K, prc=10, bimaxit=20, exact=FALSE,
     maxits=1000, convits=100, lam=0.9, includeSim=TRUE, details=FALSE,
     nonoise=FALSE, seed=NA, verbose=TRUE, ...)
\S4method{apclusterK}{character,ANY}(s, x, K, prc=10, bimaxit=20, exact=FALSE,
     maxits=1000, convits=100, lam=0.9, includeSim=TRUE, details=FALSE,
     nonoise=FALSE, seed=NA, verbose=TRUE, ...)
}
\arguments{
  \item{s}{an \eqn{l\times l}{lxl} similarity matrix in sparse or dense
    format or a similarity
    function either specified as the name of a package-provided
    similarity function as character string or a user provided
    function object.} 
  \item{x}{input data to be clustered; if \code{x} is a matrix or data
    frame, rows are interpreted as samples and columns are 
    interpreted as features; apart from matrices or data frames, 
    \code{x} may be any other structured data type that
    contains multiple data items - provided that an appropriate
    \code{\link[base:length]{length}} function is available that
    returns the number of items}
  \item{K}{desired number of clusters}
  \item{prc}{the algorithm stops if the number of clusters does not
    deviate more than prc percent from desired value K; set to 0 if you
    want to have exactly K clusters}
  \item{bimaxit}{maximum number of bisection steps to perform; note that
    no warning is issued if the number of clusters is still not in
    the desired range}
  \item{exact}{flag indicating whether or not to compute the initial
    preference range exactly (see \code{\link{preferenceRange}})}
  \item{maxits}{maximal number of iterations that \code{\link{apcluster}}
    should execute}
  \item{convits}{\code{\link{apcluster}} terminates if the examplars have not
                 changed for \code{convits} iterations}
  \item{lam}{damping factor for  \code{\link{apcluster}}; should be a value in
             the range [0.5, 1);
             higher values correspond to heavy damping which may be
             needed if oscillations occur}
  \item{includeSim}{if \code{TRUE}, the similarity matrix (either computed
    internally or passed via the \code{s} argument) is stored to the
    slot \code{sim} of the returned
    \code{\linkS4class{APResult}} object. The default is \code{FALSE}
    if \code{apclusterK} has been called for a similarity matrix,
    otherwise the default is \code{TRUE}.}
  \item{details}{if \code{TRUE}, more detailed information about the
                 algorithm's progress is stored in the output object
                 (see \code{\linkS4class{APResult}})}
  \item{nonoise}{\code{\link{apcluster}} adds a small amount of noise to
                 \code{s} to prevent degenerate cases; if \code{TRUE},
                 this is disabled}
  \item{seed}{for reproducibility, the seed of the random number
    generator can be set to a fixed value, if \code{NA},
    the seed remains unchanged}
  \item{verbose}{flag indicating whether status information should be
	         displayed during bisection}
  \item{...}{for the methods with signatures \code{character,ANY} and
    \code{function,ANY}, all other arguments are passed to the selected 
    similarity function as they are; for the methods with signatures
    \code{Matrix,missing} and \code{sparseMatrix,missing}, further
    arguments are passed on to the \code{apclusterK} methods with
    signatures \code{Matrix,missing} and \code{dgTMatrix,missing},
    respectively.}
}
\details{
  \code{apclusterK} first runs \code{\link{preferenceRange}} to determine
  the range of meaningful choices of the input preference \code{p}. Then
  it decreases \code{p} exponentially for a few iterations to obtain a
  good initial guess for \code{p}. If the number of clusters is still
  too far from the desired goal, bisection is applied.

  When called with a similarity matrix as input, the function performs
  the procedure described above. When called with the name of a package-provided
  similarity function or a user-provided similarity function object and 
  input data, the function first computes the similarity matrix before
  running \code{apclusterK} on this similarity matrix. The similarity
  matrix is returned for later use as part of the APResult object
  depending on whether \code{includeSim} was set to \code{TRUE} (see
  argument description above).
 
  Apart from minor adaptations and optimizations, the implementation is
  largely analogous to Frey's and Dueck's Matlab code
  (see \url{http://genes.toronto.edu/affinitypropagation/}).
}
\value{
  Upon successful completion, the function returns a
  \code{\linkS4class{APResult}} object.
}
\author{Ulrich Bodenhofer & Andreas Kothmeier
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
\seealso{\code{\link{apcluster}}, \code{\link{preferenceRange}},
         \code{\linkS4class{APResult}}}
\examples{
## create three Gaussian clouds
cl1 <- cbind(rnorm(70, 0.2, 0.05), rnorm(70, 0.8, 0.06))
cl2 <- cbind(rnorm(50, 0.7, 0.08), rnorm(50, 0.3, 0.05))
cl3 <- cbind(rnorm(60, 0.8, 0.04), rnorm(60, 0.8, 0.05))
x <- rbind(cl1, cl2, cl3)

## run affinity propagation such that 3 clusters are obtained
apres <- apclusterK(negDistMat(r=2), x, K=3)

## show details of clustering results
show(apres)

## plot clustering result
plot(apres, x)

## create sparse similarity matrix
cl1 <- cbind(rnorm(20, 0.2, 0.05), rnorm(20, 0.8, 0.06))
cl2 <- cbind(rnorm(20, 0.7, 0.08), rnorm(20, 0.3, 0.05))
x <- rbind(cl1, cl2)

sim <- negDistMat(x, r=2)
ssim <- as.SparseSimilarityMatrix(sim, lower=-0.2)

## run apcluster() on the sparse similarity matrix
apres <- apclusterK(ssim, K=2)
apres
}
\keyword{cluster}

