\name{labels-methods}
\docType{methods}
\alias{labels}
\alias{labels-methods}
\alias{labels,APResult-method}
\alias{labels,ExClust-method}
\title{Generate label vector from clustering result}
\description{
  Generate a label vector from an clustering result
}
\usage{
\S4method{labels}{ExClust}(object, type="names")
}
\arguments{
  \item{object}{object of class \code{\linkS4class{APResult}}
    or \code{\linkS4class{ExClust}}}
  \item{type}{specifies which kind of label vector should be
    created, see details below}
}
\details{
  The function \code{labels} creates a label vector from a clustering
  result. Which kind of labels are produced is controlled by the
  argument \code{type}:
  \describe{
    \item{\dQuote{names}}{(default) returns the name of the exemplar
    to which each
    data sample belongs to; if no names are available, the function
    stops with an error;}
    \item{\dQuote{enum}}{returns the index of the cluster to which
    each data sample belongs to, where clusters are enumerated
    consecutively from 1 to the number of clusters (analogous to
    other clustering methods like \code{\link{kmeans}});}
    \item{\dQuote{exemplars}}{returns the index of the exemplar to
    which each data sample belongs to, where indices of exemplars are
    within the original data, which is nothing else but the slot
    \code{object@idx} with attributes removed.}}
}
\value{
  returns a label vector as long as the number of samples in the
  original data set
}
\author{Ulrich Bodenhofer & Andreas Kothmeier
\email{apcluster@bioinf.jku.at}}
\references{\url{http://www.bioinf.jku.at/software/apcluster}

Bodenhofer, U., Kothmeier, A., and Hochreiter, S. (2011)
APCluster: an R package for affinity propagation clustering.
\emph{Bioinformatics} \bold{27}, 2463-2464.
DOI: \href{http://dx.doi.org/10.1093/bioinformatics/btr406}{10.1093/bioinformatics/btr406}.
}
\seealso{\code{\linkS4class{APResult}},
  \code{\linkS4class{ExClust}}, \code{\link{cutree}}}
\examples{
## create two simple clusters
x <- c(1, 2, 3, 7, 8, 9)
names(x) <- c("a", "b", "c", "d", "e", "f")

## compute similarity matrix (negative squared distance)
sim <- negDistMat(x, r=2)

## run affinity propagation
apres <- apcluster(sim)

## show details of clustering results
show(apres)

## label vector (names of exemplars)
labels(apres)

## label vector (consecutive index of exemplars)
labels(apres, type="enum")

## label vector (index of exemplars within original data set)
labels(apres, type="exemplars")

## now with agglomerative clustering
aggres <- aggExCluster(sim)

## label (names of exemplars)
labels(cutree(aggres, 2))
}
\keyword{cluster}
\keyword{methods}
