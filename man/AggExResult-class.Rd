\name{AggExResult-class}
\docType{class}
\alias{AggExResult-class}
\alias{AggExResult}
\alias{aggexresult}
\alias{[,AggExResult,index,missing,missing-method}
\alias{[[,AggExResult,index,missing-method}
\alias{length,AggExResult-method}
\alias{similarity,AggExResult-method}

\title{Class "AggExResult"}
\description{S4 class for storing results of exemplar-based
  agglomerative clustering}
\section{Objects}{
Objects of this class can be created by calling \code{\link{aggExCluster}}
for a given similarity matrix.
}
\section{Slots}{
The following slots are defined for \link{AggExResult} objects:
    \describe{
    \item{\code{l}:}{number of samples in the data set}
    \item{\code{sel}:}{subset of samples used for leveraged clustering 
    	               (empty for normal clustering)}
    \item{\code{maxNoClusters}:}{maximum number of clusters in the 
      cluster hierarchy, i.e. it
      contains clusterings with 1 - \code{maxNoClusters} clusters.}
    \item{\code{exemplars}:}{list of length \code{maxNoClusters};
      the \code{i}-th component of the list is a vector of \code{i}
      exemplars (corresponding to the level with \code{i} clusters).}
    \item{\code{clusters}:}{list of length \code{maxNoClusters};
      the \code{i}-th component of \code{clusters} is a list of \code{i}
      clusters, each of which is a vector of sample indices.}
    \item{\code{merge}:}{a \code{maxNoClusters-1} by 2 matrix that
      contains the merging hierarchy; fully analogous to the
      slot \code{merge} in the class \code{\link{hclust}}.}
    \item{\code{height}:}{a vector of length \code{maxNoClusters-1} that
      contains the merging objective of each merge; largely analogous to
      the slot \code{height} in the class \code{\link{hclust}} except
      that the slot \code{height} in \code{AggExResult} objects is
      supposed to be non-increasing, since \code{\link{aggExCluster}}
      is based on similarities, whereas \code{\link{hclust}} uses
      dissimilarities.}
    \item{\code{order}:}{a vector containing a permutation of indices
      that can be used for plotting proper dendrograms without crossing
      branches; fully analogous to the
      slot \code{order} in the class \code{\link{hclust}}.}
    \item{\code{labels}:}{a character vector containing labels of
      clustered objects used for plotting dendrograms.} 
    \item{\code{sim}:}{similarity matrix; only available if
		       \code{\link{aggExCluster}} was called with similarity
		       function and \code{includeSim=TRUE}.}
    \item{\code{call}:}{method call used to produce this
    	                clustering result}
  }
}
\section{Methods}{
  \describe{
    \item{plot}{\code{signature(x="AggExResult")}: see
       \code{\link{plot-methods}}}
    \item{plot}{\code{signature(x="AggExResult", y="matrix")}: see
       \code{\link{plot-methods}}}
    \item{heatmap}{\code{signature(x="AggExResult")}: see
       \code{\link{heatmap-methods}}}
    \item{heatmap}{\code{signature(x="AggExResult", y="matrix")}: see
       \code{\link{heatmap-methods}}}
    \item{show}{\code{signature(object="AggExResult")}: see
       \code{\link{show-methods}}}
    \item{cutree}{\code{signature(object="AggExResult", k="ANY",
	h="ANY")}: see \code{\link{cutree-methods}}}
    \item{length}{\code{signature(x="AggExResult")}: gives the number of
      clustering levels in the clustering result.}
    \item{as.hclust}{\code{signature(x="AggExResult")}: see
       \code{\link{coerce-methods}}}
    \item{as.dendrogram}{\code{signature(object="AggExResult")}: see
       \code{\link{coerce-methods}}}
    }
}
\section{Accessors}{
  In the following code snippets, \code{x} is an \code{AggExResult} object.
  \describe{
    \item{}{\code{x[[i]]}: Returns an object of class
      \code{\linkS4class{ExClust}} corresponding to the clustering level
      with \code{i} clusters; synonymous to \code{\link{cutree}(x, i)}.
    }
    \item{}{\code{x[i]}: Returns a list of \code{\linkS4class{ExClust}}
      objects with all clustering levels specified in vector \code{i}.
      So, the list has as many components as the argument \code{i} has
      elements. A list is returned even if \code{i} is a single level.
    }
    \item{similarity}{\code{signature(x="AggExResult")}: gives the similarity
      matrix.
    }
  }
}
\author{Ulrich Bodenhofer, Andreas Kothmeier & Johannes Palme
\email{apcluster@bioinf.jku.at}}
\references{\url{http://www.bioinf.jku.at/software/apcluster}

Bodenhofer, U., Kothmeier, A., and Hochreiter, S. (2011)
APCluster: an R package for affinity propagation clustering.
\emph{Bioinformatics} \bold{27}, 2463-2464.
DOI: \href{http://dx.doi.org/10.1093/bioinformatics/btr406}{10.1093/bioinformatics/btr406}.
}
\seealso{\code{\link{aggExCluster}}, \code{\link{show-methods}},
  \code{\link{plot-methods}}, \code{\link{cutree-methods}}}
\examples{
## create two Gaussian clouds
cl1 <- cbind(rnorm(50, 0.2, 0.05), rnorm(50, 0.8, 0.06))
cl2 <- cbind(rnorm(50, 0.7, 0.08), rnorm(50, 0.3, 0.05))
x <- rbind(cl1, cl2)

## compute similarity matrix (negative squared Euclidean)
sim <- negDistMat(x, r=2)

## compute agglomerative clustering from scratch
aggres1 <- aggExCluster(sim)

## show results
show(aggres1)

## plot dendrogram
plot(aggres1)

## plot heatmap along with dendrogram
heatmap(aggres1, sim)

## plot level with two clusters
plot(aggres1, x, k=2)

## run affinity propagation
apres <- apcluster(sim, q=0.7)

## create hierarchy of clusters determined by affinity propagation
aggres2 <- aggExCluster(sim, apres)

## show results
show(aggres2)

## plot dendrogram
plot(aggres2)

## plot heatmap
heatmap(aggres2, sim)

## plot level with two clusters
plot(aggres2, x, k=2)
}
\keyword{classes}
