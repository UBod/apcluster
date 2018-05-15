\name{cutree-methods}
\docType{methods}
\alias{cutree}
\alias{cutree-methods}
\alias{cutree,AggExResult-method}
\alias{cutree,APResult-method}
\title{Cut Out Clustering Level from Cluster Hierarchy}
\description{
  Cut out a clustering level from a cluster hierarchy
}
\usage{
\S4method{cutree}{AggExResult}(tree, k, h)
\S4method{cutree}{APResult}(tree, k, h)
}
\arguments{
  \item{tree}{an object of class \code{\linkS4class{AggExResult}}
    containing a cluster hierarchy; may also be an object of
    class \code{\linkS4class{APResult}}}
  \item{k}{the level (i.e. the number of clusters) to be selected}
  \item{h}{alternatively, the level can be selected by specifying a
    cut-off for the merging objective}
}
\details{
  The function \code{cutree} extracts a clustering level from a
  cluster hierarchy stored in an \code{\linkS4class{AggExResult}}
  object. Which level is selected can be determined by one of the
  two arguments \code{k} and \code{h} (see above). If both \code{k} and
  \code{h} are specified, \code{k} overrides \code{h}. This is
  done largely analogous to the standard function 
  \code{\link[stats:cutree]{cutree}}. The differences are (1) that
  only one level can be extracted at a time and (2) that an
  \code{\linkS4class{ExClust}} is returned instead of an index list.

  The function \code{cutree} may further be used to convert an
  \code{\linkS4class{APResult}} object into an
  \code{\linkS4class{ExClust}} object. In this case, the arguments
  \code{k} and \code{h} are ignored.
}
\value{
  returns an object of class \code{\linkS4class{ExClust}}
}
\author{Ulrich Bodenhofer & Andreas Kothmeier
\email{apcluster@bioinf.jku.at}}
\references{\url{http://www.bioinf.jku.at/software/apcluster}

Bodenhofer, U., Kothmeier, A., and Hochreiter, S. (2011)
APCluster: an R package for affinity propagation clustering.
\emph{Bioinformatics} \bold{27}, 2463-2464.
DOI: \href{http://dx.doi.org/10.1093/bioinformatics/btr406}{10.1093/bioinformatics/btr406}.
}
\seealso{\code{\linkS4class{AggExResult}},
  \code{\linkS4class{ExClust}}}
\examples{
## create two simple clusters
x <- c(1, 2, 3, 7, 8, 9)
names(x) <- c("a", "b", "c", "d", "e", "f")

## compute similarity matrix (negative squared distance)
sim <- negDistMat(x, r=2)

## run affinity propagation
aggres <- aggExCluster(sim)

## show details of clustering results
show(aggres)

## retrieve clustering with 2 clusters
cutree(aggres, 2)

## retrieve clustering with cut-off h=-1
cutree(aggres, h=-1)
}
\keyword{cluster}
\keyword{methods}
