\name{sort-methods}
\docType{methods}
\alias{sort}
\alias{sort-methods}
\alias{sort,APResult-method}
\alias{sort,ExClust-method}
\title{Sort clusters}
\description{
  Rearrange clusters according to sort criterion
}
\usage{
\S4method{sort}{ExClust}(x, decreasing=FALSE,
                   sortBy=c("aggExCluster", "size",
                            "nameExemplar", "noExemplar"), ...)
}
\arguments{
  \item{x}{object of class \code{\linkS4class{APResult}}
    or \code{\linkS4class{ExClust}}}
  \item{decreasing}{logical indicating if sorting should be done in
    decreasing order, see details below}
  \item{sortBy}{sort criterion, see details below}
  \item{...}{further arguments are ignored; only defined for S3 method
    consistency}
}
\details{
  The function \code{sort} takes an \code{\linkS4class{APResult}}
  or \code{\linkS4class{ExClust}} clustering object \code{x} and creates
  a new clustering object of the same class, but with clusters arranged
  according to the sort criterion passed as argument \code{sortBy}:
  \describe{
    \item{\dQuote{aggExCluster}}{(default) order clusters as they
      would appear in the dendrogram produced by
      \code{\link{aggExCluster}}. This is also the same ordering in
      which the clusters are arranged by \code{\link{heatmap}}.
      Note that this only works if the similarity matrix is included
      in the input object \code{x}, otherwise an error message is
      produced.}
    \item{\dQuote{size}}{sorts clusters according to their size (from
      small to large).}
    \item{\dQuote{nameExemplar}}{sorts clusters according to the names
      of the examplars (if available, otherwise an error is produced).}
    \item{\dQuote{noExemplar}}{sorts clusters according to the indices
      of the examplars.}
  }
  If \code{decreasing} is \code{TRUE}, the order is reversed and, for
  example, \code{sortBy="size"} sorts clusters with such that the larger
  clusters come first.

  Note that the cluster numbers of \code{x} are not preserved by
  \code{sort}, i.e. the cluster no. 1 of the object returned by
  \code{sort} is the one that has been ranked first by \code{sort},
  which may not necessarily coincide with cluster no. 1 of the original
  clustering object \code{x}.

  Note that this is an S3 method (whereas all other methods in this
  package are S4 methods). This inconsistency has been introduced in
  order to avoid interoperability problems with the \pkg{BiocGenerics}
  package which may overwrite the definition of the \code{sort} generic
  if it is loaded after the \pkg{apcluster} package.
}
\value{
  returns a copy of \code{x}, but with slots \code{exemplars} and
  \code{clusters} (see \code{\linkS4class{APResult}}
  or \code{\linkS4class{ExClust}}) reordered.
}
\author{Ulrich Bodenhofer & Andreas Kothmeier
\email{apcluster@bioinf.jku.at}}
\references{\url{http://www.bioinf.jku.at/software/apcluster}

Bodenhofer, U., Kothmeier, A., and Hochreiter, S. (2011)
APCluster: an R package for affinity propagation clustering.
\emph{Bioinformatics} \bold{27}, 2463-2464.
DOI: \href{http://dx.doi.org/10.1093/bioinformatics/btr406}{10.1093/bioinformatics/btr406}.
}
\seealso{\code{\linkS4class{APResult}}, \code{\linkS4class{ExClust}}}
\examples{
## create two Gaussian clouds
cl1 <- cbind(rnorm(50,0.2,0.05),rnorm(50,0.8,0.06))
cl2 <- cbind(rnorm(50,0.7,0.08),rnorm(50,0.3,0.05))
x <- rbind(cl1,cl2)

## run affinity propagation
apres <- apcluster(negDistMat(r=2), x, q=0.7)

show(apres)
## show dendrogram
plot(aggExCluster(x=apres))

## default sort order: like in heatmap or dendrogram
show(sort(apres))

## show dendrogram (note the different cluster numbers!)
plot(aggExCluster(x=sort(apres)))

## sort by size
show(sort(apres, decreasing=TRUE, sortBy="size"))
}
\keyword{cluster}
\keyword{methods}
