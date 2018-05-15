\name{APResult-class}
\docType{class}
\alias{APResult-class}
\alias{APResult}
\alias{apresult}
\alias{similarity}
\alias{[,APResult,index,missing,missing-method}
\alias{[[,APResult,index,missing-method}
\alias{length,APResult-method}
\alias{similarity,APResult-method}

\title{Class "APResult"}
\description{S4 class for storing results of affinity propagation
clustering. It extends the class \code{\linkS4class{ExClust}}.}
\section{Objects}{
Objects of this class can be created by calling \code{\link{apcluster}}
or \code{\link{apclusterL}} for a given similarity matrix or calling
one of these procedures with a data set and a similarity measure.
}
\section{Slots}{
The following slots are defined for \link{APResult} objects. Most names
are taken from Frey's and Dueck's original Matlab package:
    \describe{
    \item{\code{sweeps}:}{number of times leveraged clustering ran with 
    	                  different subsets of samples}
    \item{\code{it}:}{number of iterations the algorithm ran}
    \item{\code{p}:}{input preference (either set by user or
                     computed by \code{\link{apcluster}} or
                     \code{\link{apclusterL}})}
    \item{\code{netsim}:}{final total net similarity, defined as the
                          sum of \code{expref} and \code{dpsim}
                          (see below)}
    \item{\code{dpsim}:}{final sum of similarities of data points to
                         exemplars}
    \item{\code{expref}:}{final sum of preferences of the identified
                          exemplars}
    \item{\code{netsimLev}:}{total net similarity of the individual 
    	                     sweeps for leveraged clustering; only
    	                     available for leveraged clustering}
    \item{\code{netsimAll}:}{vector containing the total net similarity
                             for each iteration; only available if
                             \code{\link{apcluster}} was called with
                             \code{details=TRUE}}
    \item{\code{exprefAll}:}{vector containing the sum of preferences
                             of the identified exemplars
                             for each iteration; only available if
                             \code{\link{apcluster}} was called with
                             \code{details=TRUE}}
    \item{\code{dpsimAll}:}{vector containing the sum of similarities
                            of data points to exemplars
                            for each iteration; only available if
                            \code{\link{apcluster}} was called with
                            \code{details=TRUE}}
    \item{\code{idxAll}:}{matrix with sample-to-exemplar indices
                          for each iteration; only available if
                          \code{\link{apcluster}} was called with
                          \code{details=TRUE}}
  }
}
\section{Extends}{
Class \code{"ExClust"}, directly.
}
\section{Methods}{
  \describe{
    \item{plot}{\code{signature(x="APResult")}: see
       \code{\link{plot-methods}}}
    \item{plot}{\code{signature(x="ExClust", y="matrix")}: see
       \code{\link{plot-methods}}}
    \item{heatmap}{\code{signature(x="ExClust")}: see
       \code{\link{heatmap-methods}}}
    \item{heatmap}{\code{signature(x="ExClust", y="matrix")}: see
       \code{\link{heatmap-methods}}}
    \item{show}{\code{signature(object="APResult")}: see
       \code{\link{show-methods}}}
    \item{labels}{\code{signature(object="APResult")}: see
       \code{\link{labels-methods}}}
    \item{cutree}{\code{signature(object="APResult")}: see
       \code{\link{cutree-methods}}}
    \item{length}{\code{signature(x="APResult")}: gives the number of
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
  In the following code snippets, \code{x} is an \code{APResult} object.
  \describe{
    \item{}{\code{x[[i]]}: Returns the i-th cluster as a list of indices
      of samples belonging to the i-th cluster.
    }
    \item{}{\code{x[i]}: Returns a list of integer vectors with the
      indices  of samples belonging to this cluster. The list has as
      many components as the argument \code{i} has elements. A list is
      returned even if \code{i} is a single integer.
    }
    \item{similarity}{\code{signature(x="APResult")}: gives the similarity
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

Frey, B. J. and Dueck, D. (2007) Clustering by passing messages
between data points. \emph{Science} \bold{315}, 972-976.
}
\seealso{\code{\link{apcluster}}, \code{\link{apclusterL}}, 
  \code{\link{show-methods}}, \code{\link{plot-methods}}, 
  \code{\link{labels-methods}}, \code{\link{cutree-methods}}}
\examples{
## create two Gaussian clouds
cl1 <- cbind(rnorm(100, 0.2, 0.05), rnorm(100, 0.8, 0.06))
cl2 <- cbind(rnorm(50, 0.7, 0.08), rnorm(50, 0.3, 0.05))
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

## plot heatmap
heatmap(apres, sim)
}
\keyword{classes}
