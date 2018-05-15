\name{coerce-methods}
\docType{methods}
\alias{as.hclust}
\alias{as.dendrogram}
\alias{as.hclust,AggExResult-method}
\alias{as.hclust,ExClust-method}
\alias{as.dendrogram,AggExResult-method}
\alias{as.dendrogram,ExClust-method}
\title{Coercion of cluster hierarchies}
\description{
  Functions for coercing clustering object to hclust and dendrogram objects
}
\usage{
\S4method{as.hclust}{AggExResult}(x, base=0.05)
\S4method{as.hclust}{ExClust}(x, base=0.05, ...)
\S4method{as.dendrogram}{AggExResult}(object, base=0.05, useNames=TRUE)
\S4method{as.dendrogram}{ExClust}(object, base=0.05, useNames=TRUE, ...)
}
\arguments{
  \item{x}{a clustering result object of class
     \code{\linkS4class{APResult}}, \code{\linkS4class{ExClust}},
     or \code{\linkS4class{AggExResult}}}
  \item{object}{a clustering result object of class
     \code{\linkS4class{APResult}}, \code{\linkS4class{ExClust}},
     or \code{\linkS4class{AggExResult}}}
  \item{base}{fraction of height used for the very first join; defaults
    to 0.05, i.e. the first join appears at 5\% of the total height of the
    dendrogram (see details below).}
  \item{useNames}{if \code{TRUE} (default), the labels of the dendrogram
    are the sample/cluster names (if available); otherwise, the labels
    are indices.}
  \item{...}{all other arguments are passed on to
    \code{\link{aggExCluster}} (see details below).}
}
\details{
  If called for an \code{\linkS4class{AggExResult}} object,
  \code{as.hclust} creates an \code{\link[stats:hclust]{hclust}} object.
  The heights are transformed to the interval from \code{base} (height
  of lowest join) to 1 (height of highest join).
  If called for an \code{\linkS4class{ExClust}} or
  \code{\linkS4class{APResult}} object, \code{\link{aggExCluster}} is
  called internally to create a cluster hierarchy first. This is only
  possible if the pairwise similarities are included in the \code{sim}
  slot of \code{x} (see \code{\link{aggExCluster}} on how to ensure
  this).

  If \code{x} is an \code{\linkS4class{AggExResult}} object obtained by
  clustering an entire data set, \code{as.hclust} produces a complete
  hierarchy. If, however, \code{x} is an \code{\linkS4class{ExClust}} (or
  \code{\linkS4class{APResult}}) object or an
  \code{\linkS4class{AggExResult}} obtained by running
  \code{\link{aggExCluster}} on an \code{\linkS4class{ExClust}} or
  \code{\linkS4class{APResult}} object, then \code{as.hclust} produces
  a hierarchy of clusters, not of samples.

  If called for an \code{\linkS4class{AggExResult}} object,
  \code{as.dendrogram} creates an
  \code{\link[stats:dendrogram]{dendrogram}} object.
  Analogously to \code{as.hclust}, the heights are transformed to the
  interval ranging from \code{base} (height
  of lowest join) to 1 (height of highest join). So, any information
  about heights of merges is lost. If the original join heights are
  relevant, call \code{\link[apcluster:plot]{plot}} on the original
  \code{\linkS4class{AggExResult}} object directly without coercing it
  to a \code{\link[stats:dendrogram]{dendrogram}} object first.
  If called for an \code{\linkS4class{ExClust}} or
  \code{\linkS4class{APResult}} object, \code{\link{aggExCluster}} is
  called first to create a cluster hierarchy. Again this is only
  possible if the pairwise similarities are included in the \code{sim}
  slot of \code{object}.

  If \code{object} is an \code{\linkS4class{AggExResult}} object obtained by
  clustering an entire data set, \code{as.dendrogram} produces a complete
  dendrogram. If \code{object} is an \code{\linkS4class{ExClust}} (or
  \code{\linkS4class{APResult}}) object or an
  \code{\linkS4class{AggExResult}} obtained by previously running
  \code{\link{aggExCluster}} on an \code{\linkS4class{ExClust}} or
  \code{\linkS4class{APResult}} object, then \code{as.dendrogram} produces
  a complete dendrogram of all samples, too, but with the difference
  that entire clusters of the previous \code{\linkS4class{ExClust}} or
  \code{\linkS4class{APResult}} object are not further split up
  hierarchically.
  Consequently, if \code{x} is not a complete cluster hierarchy, but a
  hierarchy of clusters, \code{as.dendrogram(as.hclust(x))} produces a
  dendrogram of clusters, whereas \code{as.dendrogram(x)} in any case
  produces a dendrogram of samples (with the special property mentioned
  above).
}
\value{
  see details above
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
  \code{\linkS4class{AggExResult}}, \code{\linkS4class{ExClust}},
  \code{\link{heatmap-methods}}, \code{\link{apcluster}}, 
  \code{\link{apclusterL}}, \code{\link{aggExCluster}}, 
  \code{\link{cutree-methods}}}
\examples{
## create two Gaussian clouds
cl1 <- cbind(rnorm(20, 0.2, 0.05), rnorm(20, 0.8, 0.06))
cl2 <- cbind(rnorm(20, 0.7, 0.08), rnorm(20, 0.3, 0.05))
x <- rbind(cl1, cl2)

## run affinity propagation
apres <- apcluster(negDistMat(r=2), x, q=0.7, details=TRUE)

## perform agglomerative clustering of affinity propagation clusters
aggres1 <- aggExCluster(x=apres)

## compute and plot dendrogram
dend1 <- as.dendrogram(aggres1)
dend1
plot(dend1)

## compute and show dendrogram computed from hclust object
dend2 <- as.dendrogram(as.hclust(aggres1))
dend2
plot(dend2)

## perform agglomerative clustering of whole data set
aggres2 <- aggExCluster(negDistMat(r=2), x)

## compute and plot dendrogram
dend3 <- as.dendrogram(aggres2)
dend3
plot(dend3)
}
\keyword{cluster}
\keyword{methods}
