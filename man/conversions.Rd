\name{conversions}
\alias{conversions}
\alias{as.SparseSimilarityMatrix}
\alias{as.DenseSimilarityMatrix}
\alias{as.SparseSimilarityMatrix-methods}
\alias{as.DenseSimilarityMatrix-method}
\alias{as.SparseSimilarityMatrix,matrix-method}
\alias{as.DenseSimilarityMatrix,matrix-method}
\alias{as.SparseSimilarityMatrix,sparseMatrix-method}
\alias{as.DenseSimilarityMatrix,sparseMatrix-method}
\alias{as.SparseSimilarityMatrix,Matrix-method}
\alias{as.DenseSimilarityMatrix,Matrix-method}
\title{Conversions Between Dense and Sparse Similarity Matrices}
\description{
   Converts a dense similarity matrix into a sparse one or vice versa
}
\usage{
\S4method{as.SparseSimilarityMatrix}{matrix}(s, lower=-Inf)
\S4method{as.SparseSimilarityMatrix}{Matrix}(s, lower=-Inf)
\S4method{as.SparseSimilarityMatrix}{sparseMatrix}(s, lower=-Inf)
\S4method{as.DenseSimilarityMatrix}{matrix}(s, fill=-Inf)
\S4method{as.DenseSimilarityMatrix}{Matrix}(s, fill=-Inf)
\S4method{as.DenseSimilarityMatrix}{sparseMatrix}(s, fill=-Inf)
}

\arguments{
  \item{s}{a similarity matrix in sparse or dense format (see details
    below)}
  \item{lower}{cut-off threshold to apply when converting similarity
    matrices into sparse format. All similarities lower than or equal to
    \code{lower} will be omitted from the result. The default is
    \code{-Inf}), i.e. only \code{-Inf} values are removed.}
  \item{fill}{value to fill in for entries that are missing from sparse
	   similarity matrix 's' (defaults to \code{-Inf}).}
}
\details{
  The function \code{as.SparseSimilarityMatrix} takes a matrix argument,
  removes all diagonal elements and all values that are lower than or
  equal to the cut-off threshold \code{lower} and returns a sparse
  matrix of class \code{\linkS4class{dgTMatrix}}.

  If the function \code{as.DenseSimilarityMatrix} is called for a
  sparse matrix (class \code{\linkS4class{sparseMatrix}} or any
  class derived from this class), a dense matrix is returned, where all
  values that were missing in the sparse matrix are replaced with
  \code{fill}.

  \code{as.DenseSimilarityMatrix} can also be called for dense
  \code{\link{matrix}} and \code{\linkS4class{Matrix}} objects.
  In this case, \code{as.DenseSimilarityMatrix} assumes that the
  matrices have three columns that encode for a sparse matrix
  in the same way as the Matlab implementation of Frey's and Dueck's
  sparse affinity propagation accepts it:
  the first column contains 1-based row indices, the second column
  contains 1-based column indices, and the third column contains the
  similarity values. The same format is also accepted by
  \code{as.SparseSimilarityMatrix} to convert a sparse similarity matrix
  of this format into a \code{\linkS4class{dgTMatrix}} object.
  Note that, for matrices of this format,
  \code{as.DenseSimilarityMatrix} replaces the deprectated function
  \code{sparseToFull} (see \code{\link{apcluster-deprecated}}).

  Note that \code{as.SparseSimilarityMatrix} and
  \code{as.DenseSimilarityMatrix} are no S4 coercion methods.
  There are no classes named \code{SparseSimilarityMatrix}
  or \code{DenseSimilarityMatrix}.
}
\value{returns a square similarity matrix in sparse format (class
  \code{\linkS4class{dgTMatrix}} or in dense format (standard class
  \code{\link{matrix}}).}
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
\examples{
## create similarity matrix in sparse format according to Frey and Dueck
sp <- matrix(c(1, 2, 0.5, 3, 1, 0.2, 5, 4, -0.2, 3, 4, 1.2), 4, 3, byrow=TRUE)
sp

## perform conversions
as.DenseSimilarityMatrix(sp, fill=0)
as.SparseSimilarityMatrix(sp)

## create dense similarity matrix
cl1 <- cbind(rnorm(20, 0.2, 0.05), rnorm(20, 0.8, 0.06))
cl2 <- cbind(rnorm(20, 0.7, 0.08), rnorm(20, 0.3, 0.05))
x <- rbind(cl1, cl2)

sim <- negDistMat(x, r=2)
ssim <- as.SparseSimilarityMatrix(sim, lower=-0.2)

## run apcluster() on the sparse similarity matrix
apres <- apcluster(ssim, q=0)
apres
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{cluster}

