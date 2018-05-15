\name{show-methods}
\docType{methods}
\alias{show}
\alias{show-methods}
\alias{show,APResult-method}
\alias{show,ExClust-method}
\alias{show,AggExResult-method}
\title{Display Clustering Result Objects}
\description{
  Display methods for S4 classes \code{\linkS4class{APResult}},
  \code{\linkS4class{ExClust}}, and \code{\linkS4class{AggExResult}}
}
\usage{
\S4method{show}{APResult}(object)
\S4method{show}{ExClust}(object)
\S4method{show}{AggExResult}(object)
}
\arguments{
  \item{object}{an object of class
    \code{\linkS4class{APResult}}, \code{\linkS4class{ExClust}}, or
    \code{\linkS4class{AggExResult}}
  }
}
\details{
  \code{show} displays the most important information stored in
  \code{object}.

  For \code{\linkS4class{APResult}} objects,
  the number of data samples, the number of clusters, the number of
  iterations, the input preference, the final objective
  function values, the vector of exemplars, the list of clusters and
  for leveraged clustering the selected sample subset are printed.  

  For \code{\linkS4class{ExClust}} objects,
  the number of data samples, the number of clusters,
  the vector of exemplars, and list of clusters are printed.

  For \code{\linkS4class{AggExResult}} objects,
  only the number of data samples and the maximum
  number of clusters are printed. For retrieving a particular
  clustering level, use the function \code{\link{cutree}}.

  For accessing more detailed information, it is necessary to
  access the slots of \code{object} directly. Use
  \code{\link{str}} to get a compact overview of all slots of an object.
}
\value{
  \code{show} returns an invisible \code{NULL}
}
\author{Ulrich Bodenhofer, Andreas Kothmeier & Johannes Palme
\email{apcluster@bioinf.jku.at}}
\references{\url{http://www.bioinf.jku.at/software/apcluster}

Bodenhofer, U., Kothmeier, A., and Hochreiter, S. (2011)
APCluster: an R package for affinity propagation clustering.
\emph{Bioinformatics} \bold{27}, 2463-2464.
DOI: \href{http://dx.doi.org/10.1093/bioinformatics/btr406}{10.1093/bioinformatics/btr406}.
}
\seealso{\code{\linkS4class{APResult}},
  \code{\linkS4class{ExClust}}, \code{\linkS4class{AggExResult}},
  \code{\link{cutree-methods}}
}
\examples{
## create two Gaussian clouds
cl1 <- cbind(rnorm(100, 0.2, 0.05), rnorm(100, 0.8, 0.06))
cl2 <- cbind(rnorm(50, 0.7, 0.08), rnorm(50, 0.3, 0.05))
x <- rbind(cl1, cl2)

## compute similarity matrix (negative squared Euclidean)
sim <- negDistMat(x, r=2)

## run affinity propagation
apres <- apcluster(sim)

## show details of clustering results
show(apres)

## apply agglomerative clustering to apres
aggres <- aggExCluster(sim, apres)

## display overview of result
show(aggres)

## show clustering level with two clusters
show(cutree(aggres, 2))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{cluster}
\keyword{methods}
