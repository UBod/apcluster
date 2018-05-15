\name{aggExCluster}
\docType{methods}
\alias{aggExCluster}
\alias{aggexcluster}
\alias{aggExCluster-methods}
\alias{aggExCluster,matrix,missing-method}
\alias{aggExCluster,matrix,ExClust-method}
\alias{aggExCluster,Matrix,missing-method}
\alias{aggExCluster,Matrix,ExClust-method}
\alias{aggExCluster,missing,ExClust-method}
\alias{aggExCluster,function,ANY-method}
\alias{aggExCluster,character,ANY-method}
\title{Exemplar-based Agglomerative Clustering}
\description{
  Runs exemplar-based agglomerative clustering
}
\usage{
\S4method{aggExCluster}{matrix,missing}(s, x, includeSim=FALSE)
\S4method{aggExCluster}{matrix,ExClust}(s, x, includeSim=FALSE)
\S4method{aggExCluster}{Matrix,missing}(s, x, includeSim=FALSE)
\S4method{aggExCluster}{Matrix,ExClust}(s, x, includeSim=FALSE)
\S4method{aggExCluster}{missing,ExClust}(s, x, includeSim=TRUE)
\S4method{aggExCluster}{function,ANY}(s, x, includeSim=TRUE, ...)
\S4method{aggExCluster}{character,ANY}(s, x, includeSim=TRUE, ...)
}
\arguments{
  \item{s}{an \eqn{l\times l}{lxl} similarity matrix or a similarity
  	      function either specified as the name of a package-provided
  	      similarity function as character string or a user provided
  	      function object} 
  \item{x}{either a prior clustering of class \code{\linkS4class{ExClust}} (or
    \code{\linkS4class{APResult}}) or, if called with \code{s} being a
    function or function name, input data to be clustered (see
    \code{\link{apcluster}} for a detailed specification)}
  \item{includeSim}{if \code{TRUE}, the similarity matrix (either computed
    internally or passed via the \code{s} argument) is stored to the
    slot \code{sim} of the returned
    \code{\linkS4class{AggExResult}} object. The default is \code{FALSE}
    if \code{aggExCluster} has been called for a similarity matrix,
    otherwise the default is \code{TRUE}.}
  \item{...}{all other arguments are passed to the selected 
             similarity function as they are.}
}
\details{\code{aggExCluster} performs agglomerative clustering.
  Unlike other methods, e.g., the ones implemented in \code{\link{hclust}},
  \code{aggExCluster} is computing exemplars for each cluster and
  its merging objective is geared towards the identification of
  meaningful exemplars, too.

  For each pair of clusters, the merging objective is computed as
  follows:
  \enumerate{
    \item{An intermediate cluster is created as the union
      of the two clusters.}
    \item{The potential exemplar is selected from the intermediate
      cluster as the sample that has the largest average similarity
      to all other samples in the intermediate cluster.}
    \item{Then the average similarity of the exemplar with all
      samples in the first cluster and the average similarity with
      all samples in the second cluster is computed. These two values
      measure how well the joint exemplar describes the samples in the
      two clusters.}
    \item{The merging objective is finally computed as the
      average of the two measures above. Hence, we can consider the
      merging objective as some kind of \dQuote{balanced average
      similarity to the joint exemplar}.}
  }
  In each step, all pairs of clusters are considered and
  the pair with the largest merging objective is actually merged.
  The joint exemplar is then chosen as the exemplar of the merged
  cluster.

  \code{aggExCluster} can be used in two ways, either by performing
  agglomerative clustering of an entire data set or by performing
  agglomerative clustering of data previously clustered by
  affinity propagation or another clustering algorithm.
  \enumerate{
  \item{Agglomerative clustering of an entire data set can be
    accomplished either by calling \code{aggExCluster} on a
    quadratic similarity matrix without further argument or by
    calling \code{aggExCluster} for a function or function name
    along with data to be clustered (as argument \code{x}).
    A full agglomeration run is performed that starts from \code{l}
    clusters (all samples in separate one-element clusters) and ends
    with one cluster (all samples in one single cluster).}  
  \item{Agglomerative clustering starting from a given clustering
    result can be accomplished by calling \code{aggExCluster} for
    an \code{\linkS4class{APResult}} or \code{\linkS4class{ExClust}}
    object passed as parameter \code{x}. The similarity matrix
    can either be passed as argument \code{s} or, if missing,
    \code{aggExCluster} looks if the similarity matrix is
    included in the clustering object \code{x}. A cluster hierarchy
    with numbers of clusters ranging from the
    number of clusters in \code{x} down to 1 is created.}  
  }

  The result is stored in an \code{\linkS4class{AggExResult}} object.
  The slot \code{height} is filled with the merging
  objective of each of the \code{maxNoClusters-1} merges. The slot
  \code{order} contains a permutation of the samples/clusters for
  dendrogram plotting. The algorithm for computing this permutation
  is the same as the one used in \code{\link{hclust}}. If \code{aggExCluster}
  was called for an entire data set, the slot \code{label}
  contains the names of the objects to be clustered (if available,
  otherwise the indices are used). If \code{aggExCluster} was called
  for a prior clustering, then labels are set to \sQuote{Cluster 1},
  \sQuote{Cluster 2}, etc.
}
\note{Similarity matrices can be supplied in dense or sparse
  format. Note, however, that sparse matrices are converted to full
  dense matrices before clustering which may lead to memory and/or
  performance bottlenecks for larger data sets.}
\value{
  Upon successful completion, the function returns an
  \code{\linkS4class{AggExResult}} object.
}
\author{Ulrich Bodenhofer & Johannes Palme
\email{apcluster@bioinf.jku.at}}
\references{\url{http://www.bioinf.jku.at/software/apcluster}

Bodenhofer, U., Kothmeier, A., and Hochreiter, S. (2011)
APCluster: an R package for affinity propagation clustering.
\emph{Bioinformatics} \bold{27}, 2463-2464.
DOI: \href{http://dx.doi.org/10.1093/bioinformatics/btr406}{10.1093/bioinformatics/btr406}.
}
\seealso{\code{\linkS4class{AggExResult}}, \code{\link{apcluster-methods}},  
  \code{\link{plot-methods}}, \code{\link{heatmap-methods}},
  \code{\link{cutree-methods}}}
\examples{
## create two Gaussian clouds
cl1 <- cbind(rnorm(50, 0.2, 0.05), rnorm(50, 0.8, 0.06))
cl2 <- cbind(rnorm(50, 0.7, 0.08), rnorm(50, 0.3, 0.05))
x <- rbind(cl1, cl2)

## compute agglomerative clustering from scratch
aggres1 <- aggExCluster(negDistMat(r=2), x)

## show results
show(aggres1)

## plot dendrogram
plot(aggres1)

## plot heatmap along with dendrogram
heatmap(aggres1)

## plot level with two clusters
plot(aggres1, x, k=2)

## run affinity propagation
apres <- apcluster(negDistMat(r=2), x, q=0.7)

## create hierarchy of clusters determined by affinity propagation
aggres2 <- aggExCluster(x=apres)

## show results
show(aggres2)

## plot dendrogram
plot(aggres2)
plot(aggres2, showSamples=TRUE)

## plot heatmap
heatmap(aggres2)

## plot level with two clusters
plot(aggres2, x, k=2)
}
\keyword{cluster}
\keyword{methods}
