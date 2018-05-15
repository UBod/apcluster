\name{apcluster-package}
\alias{apcluster-package}
\docType{package}
\title{APCluster Package}
\description{
The apcluster package implements affinity propagation according to
Frey and Dueck and a method for exemplar-based agglomerative clustering.
It further offers various functions for plotting clustering results.
}
\details{The central function is \code{\link{apcluster}}. It runs affinity
propagation on a given similarity matrix or it creates a similarity matrix
for a given data set and similarity measure and runs affinity propagation
on this matrix. The function returns an \code{\linkS4class{APResult}} 
object from which the clustering itself and information about the affinity 
propagation run can be obtained. Leveraged affinity propagation clustering 
\code{\link{apclusterL}} allows efficient clustering of large datasets by 
using only a subset of the similarities. The package further implements 
an exemplar-based agglomerative clustering method \code{\link{aggExCluster}} 
that can be used for computing a complete cluster hierarchy, but also for 
joining fine-grained clusters previously obtained by affinity propagation 
clustering. Further functions are implemented to visualize the
results and to create distance matrices.
}
\author{Ulrich Bodenhofer, Andreas Kothmeier & Johannes Palme
\email{apcluster@bioinf.jku.at}}
\references{\url{http://www.bioinf.jku.at/software/apcluster}

Bodenhofer, U., Kothmeier, A., and Hochreiter, S. (2011)
APCluster: an R package for affinity propagation clustering.
\emph{Bioinformatics} \bold{27}, 2463-2464.
DOI: \href{http://dx.doi.org/10.1093/bioinformatics/btr406}{10.1093/bioinformatics/btr406}.

Frey, B. J. and Dueck, D. (2007) Clustering by passing messages
between data points. \emph{Science} \bold{315}, 972-976.
DOI: \href{http://dx.doi.org/10.1126/science.1136800}{10.1126/science.1136800}.
}
\keyword{package}
\examples{
## create two Gaussian clouds
cl1 <- cbind(rnorm(100, 0.2, 0.05), rnorm(100, 0.8, 0.06))
cl2 <- cbind(rnorm(100, 0.7, 0.08), rnorm(100, 0.3, 0.05))
x <- rbind(cl1, cl2)

## compute similarity matrix (negative squared Euclidean)
sim <- negDistMat(x, r=2)

## run affinity propagation
apres <- apcluster(sim, details=TRUE)

## show details of clustering results
show(apres)

## plot information about clustering run
plot(apres)

## plot clustering result
plot(apres, x)

## employ agglomerative clustering to join clusters
aggres <- aggExCluster(sim, apres)

## show information
show(aggres)
show(cutree(aggres, 2))

## plot dendrogram
plot(aggres)

## plot clustering result for k=2 clusters
plot(aggres, x, k=2)

## plot heatmap
heatmap(apres, sim)

## leveraged apcluster
apresL <- apclusterL(s=negDistMat(r=2), x=x,  frac=0.2, sweeps=3)

## show details of clustering results
show(apresL)

## plot clustering result
plot(apresL, x)
}
