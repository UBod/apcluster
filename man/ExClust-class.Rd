\name{ExClust-class}
\docType{class}
\alias{ExClust-class}
\alias{ExClust}
\alias{exclust}
\alias{[,ExClust,index,missing,missing-method}
\alias{[[,ExClust,index,missing-method}
\alias{length,ExClust-method}
\alias{similarity,ExClust-method}

\title{Class "ExClust"}
\description{S4 class for storing exemplar-based clusterings}
\section{Objects}{
  Objects of this class can be created by calling \code{\link{cutree}}
  to cut out a clustering level from a cluster hierarchy
  of class \code{\linkS4class{AggExResult}}. Moreover,
  \code{\link{cutree}} can also be used to convert an object of
  class \code{\linkS4class{APResult}} to class \code{ExClust}.
}
\section{Slots}{
The following slots are defined for \link{ExClust} objects:
    \describe{
    \item{\code{l}:}{number of samples in the data set}
    \item{\code{sel}:}{subset of samples used for leveraged clustering}
    \item{\code{exemplars}:}{vector containing indices of exemplars}
    \item{\code{clusters}:}{list containing the clusters; the i-th
                            component is a vector of indices of
                            data points belonging to the i-th
                            exemplar (including the exemplar itself)}
    \item{\code{idx}:}{vector of length \code{l} realizing a
                       sample-to-exemplar mapping; the i-th entry
                       contains the index of the exemplar the i-th
                       sample belongs to}
    \item{\code{sim}:}{similarity matrix; only available if
      the preceding clustering method was called with
      \code{includeSim=TRUE}.}
    \item{\code{call}:}{method call of the preceding clustering method}
  }
}
\section{Methods}{
  \describe{
    \item{plot}{\code{signature(x="ExClust")}: see
       \code{\link{plot-methods}}}
    \item{plot}{\code{signature(x="ExClust", y="matrix")}: see
       \code{\link{plot-methods}}}
    \item{heatmap}{\code{signature(x="ExClust")}: see
       \code{\link{heatmap-methods}}}
    \item{heatmap}{\code{signature(x="ExClust", y="matrix")}: see
       \code{\link{heatmap-methods}}}
    \item{show}{\code{signature(object="ExClust")}: see
       \code{\link{show-methods}}}
    \item{labels}{\code{signature(object="ExClust")}: see
       \code{\link{labels-methods}}}
    \item{cutree}{\code{signature(object="ExClust", k="ANY", h="ANY")}: see
       \code{\link{cutree-methods}}}
    \item{length}{\code{signature(x="ExClust")}: gives the number of
      clusters.}
    \item{sort}{\code{signature(x="ExClust")}: see
       \code{\link{sort-methods}}}
    \item{as.hclust}{\code{signature(x="ExClust")}: see
       \code{\link{coerce-methods}}}
    \item{as.dendrogram}{\code{signature(object="ExClust")}: see
       \code{\link{coerce-methods}}}
    }
}
\section{Accessors}{
  In the following code snippets, \code{x} is an \code{ExClust} object.
  \describe{
    \item{}{\code{x[[i]]}: Returns the i-th cluster as a list of indices
      of samples belonging to the i-th cluster.
    }
    \item{}{\code{x[i]}: Returns a list of integer vectors with the
      indices  of samples belonging to this cluster. The list has as
      many components as the argument \code{i} has elements. A list is
      returned even if \code{i} is a single integer.
    }
    \item{similarity}{\code{signature(x="ExClust")}: gives the similarity
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
  \code{\link{plot-methods}}, \code{\link{labels-methods}},
  \code{\link{cutree-methods}}, \code{\linkS4class{AggExResult}},
  \code{\linkS4class{APResult}}}
\examples{
## create two Gaussian clouds
cl1 <- cbind(rnorm(20, 0.2, 0.05), rnorm(20, 0.8, 0.06))
cl2 <- cbind(rnorm(25, 0.7, 0.08), rnorm(25, 0.3, 0.05))
x <- rbind(cl1, cl2)

## compute similarity matrix (negative squared Euclidean)
sim <- negDistMat(x, r=2)

## run affinity propagation
aggres <- aggExCluster(sim)

## extract level with two clusters
excl <- cutree(aggres, k=2)

## show details of clustering results
show(excl)

## plot information about clustering run
plot(excl, x)
}
\keyword{classes}
