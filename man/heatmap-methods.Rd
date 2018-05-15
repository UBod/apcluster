\name{heatmap}
\docType{methods}
\alias{heatmap}
\alias{heatmap-methods}
\alias{heatmap,ExClust,missing-method}
\alias{heatmap,ExClust,matrix-method}
\alias{heatmap,ExClust,Matrix-method}
\alias{heatmap,ExClust,sparseMatrix-method}
\alias{heatmap,AggExResult,missing-method}
\alias{heatmap,AggExResult,matrix-method}
\alias{heatmap,missing,matrix-method}
\alias{heatmap,matrix,missing-method}
\title{Plot Heatmap}
\description{
  Functions for Plotting of Heatmap
}
\usage{
\S4method{heatmap}{ExClust,missing}(x, y, ...)
\S4method{heatmap}{ExClust,matrix}(x, y, ...)
\S4method{heatmap}{ExClust,Matrix}(x, y, ...)
\S4method{heatmap}{ExClust,sparseMatrix}(x, y, ...)
\S4method{heatmap}{AggExResult,missing}(x, y, ...)
\S4method{heatmap}{AggExResult,matrix}(x, y, Rowv=TRUE, Colv=TRUE,
    sideColors=NULL, col=heat.colors(12),
    base=0.05,  add.expr, margins=c(5, 5, 2),
    cexRow=max(min(35 / nrow(y), 1), 0.1),
    cexCol=max(min(35 / ncol(y), 1), 0.1), main=NULL, dendScale=1,
    barScale=1, legend=c("none", "col"), ...)
\S4method{heatmap}{matrix,missing}(x, y, ...)
\S4method{heatmap}{missing,matrix}(x, y, ...)
}
\arguments{
  \item{x}{a clustering result object of class
     \code{\linkS4class{APResult}}, \code{\linkS4class{ExClust}},
     or \code{\linkS4class{AggExResult}}; for compatibility, \code{x}
     may also be a similarity matrix (see details below).}
   \item{y}{a similarity matrix}
   \item{sideColors}{character vector of colors to be used for plotting color
     bars that visualize clusters of the finest clustering level in
     \code{x}. This is done in a fashion similar to using
     \code{RowSideColors} or \code{ColSideColors} in the standard
     \code{\link[stats:heatmap]{heatmap}} function. However, color bars
     are plotted either on both sides or not at all. The
     \code{sideColors} argument determines the coloring of both
     horizontal and vertical bars.
     If \code{sideColors} is shorter than the number of clusters in the finest
     clustering level, \code{sideColors} is recycled. In any case, a
     minimum number of two colors (two elements of sideColors) is
     required. If \code{NA}, no color bars are plotted.
     If \code{NULL} (default), color bars are only plotted if
     the finest cluster level does not only consist of single samples.
     In this case, the \code{\link{rainbow}} function is used to
     compute the vector of colors which is shuffled such that dissimilar
     colors are placed next to each other in the color bar.}
   \item{col}{color ramp used for the heatmap image; see
     \code{\link[graphics:image]{image}}}
   \item{Rowv}{determines whether or not a row dendrogram should be
     plotted. If \code{FALSE} or \code{NA}, no row dendrogram is
     plotted. In any other
     case, a row dendrogram is plotted unless the number of clusters in
     the finest clustering level is less than 2. Note that, in the
     latter case, the actual values in \code{Rowv} are ignored, so this
     argument cannot be used to supply a previously computed dendrogram
     or re-ordering of elements as in the standard
     \code{\link[stats:heatmap]{heatmap}} function.}
   \item{Colv}{determines whether or not a column dendrogram should be
     plotted. Fully analogous to \code{Rowv}, except that column
     dendrograms are never plotted if the similarity  matrix \code{y}
     is non-quadratic.}
  \item{base}{fraction of height used for the very first join in
    dendrograms; see \code{\link[apcluster:coerce-methods]{coerce-methods}}.}
  \item{add.expr,margins,cexRow,cexCol,main}{largely analogous to the standard
    \code{\link[stats:heatmap]{heatmap}} function; to omit row/column
    labeling, set \code{cexRow}/\code{cexCol} to \code{0}.
    The default for \code{margins} is a vector of length 3, where the
    third element is the right-hand side margin for the color legend
    (see \code{legend} argument). It remains unused (and can also be
    omitted) if no color legend is plotted.}
  \item{dendScale}{factor scaling the width of vertical and height of
    horizontal dendrograms; values have to be larger than 0 and no
    larger than 2. The default is 1 which corresponds to the same size
    as the dendrograms plot by the standard
    \code{\link[stats:heatmap]{heatmap}} function}
  \item{barScale}{factor scaling the width of color bars; values have to
    be larger than 0 and no larger than 4. The default is 1 which
    corresponds to half the width of the color bars plot by the standard
    \code{\link[stats:heatmap]{heatmap}} function}
  \item{legend}{if \code{"col"}, then a color legend similar to 
    \code{\link[graphics:filled.contour]{filled.contour}} is added on
    the right-hand side of the heatmap plot; if \code{"none"} (default),
    no such legend is added.}
  \item{...}{see details below}
}
\details{
  The \code{heatmap} functions provide plotting of heatmaps from several
  different types of input object. The implementation is similar to the standard
  graphics function \code{\link[stats:heatmap]{heatmap}}. 
  Plotting heatmaps via the \code{\link{plot}} command as available in previous 
  versions of this package is still available for backward 
  compatibility.

  If \code{heatmap} is called for objects of classes
  \code{\linkS4class{APResult}} or \code{\linkS4class{ExClust}},
  a heatmap of the similarity matrix in slot \code{sim} of the parameter
  \code{x} is created with clusters grouped together and highlighted in
  different colors. The order of clusters is determined by running 
  \code{\link{aggExCluster}} on the clustering result \code{x}. This
  variant of \code{heatmap} returns an  invisible
  \code{\linkS4class{AggExResult}} object.

  If \code{heatmap} is called for an \code{\linkS4class{AggExResult}}
  object that contains all levels of clustering, the heatmap is
  displayed with the corresponding clustering dendrogram. If the
  \code{\linkS4class{AggExResult}} object is the result of running
  \code{\link{aggExCluster}} on a prior clustering result, the same heatmap
  plot is produced as if \code{heatmap} had been called on this
  prior clustering result, however, returning the cluster hierarchy's
  \code{\link{dendrogram}}. In the latter case, color bars are plotted
  to visualize the prior clustering result (see description of
  argument \code{sideColors} above).

  All variants described above only work if the input object \code{x}
  contains a slot \code{sim} with the similarity matrix (which is only
  the case if the preceding clustering method has been called with
  \code{includeSim=TRUE}). In case the slot \code{sim} of \code{x} does not
  contain the similarity matrix, the similarity matrix must be supplied
  as second argument \code{y}.

  All variants described above internally use \code{heatmap} with signature
  \code{AggExResult,matrix}, so all arguments list above can be used for
  all variants, as they are passed through using the \code{...}
  argument. All other arguments, analogously to the standard 
  \code{\link[stats:heatmap]{heatmap}} function, are passed on to the
  standard function \code{\link[graphics:image]{image}}. This is
  particularly useful for using alternative color schemes via the
  \code{col} argument.

  The two variants with one of the two arguments being a matrix and one
  being missing are just wrappers around the standard
  \code{\link[stats:heatmap]{heatmap}} function with the aim to provide
  compatibility with this standard case.
}
\note{Similarity matrices can be supplied in dense or sparse
  format. Note, however, that sparse matrices are converted to full
  dense matrices before plotting heatmaps which may lead to memory and/or
  performance bottlenecks for larger data sets.}
\value{
  see details above
}
\author{Ulrich Bodenhofer, Andreas Kothmeier & Johannes Palme
\email{apcluster@bioinf.jku.at}}
\references{\url{http://www.bioinf.jku.at/software/apcluster}

Bodenhofer, U., Kothmeier, A., and Hochreiter, S. (2011)
APCluster: an R package for affinity propagation clustering.
\emph{Bioinformatics} \bold{27}, 2463-2464.
DOI: \href{http://dx.doi.org/10.1093/bioinformatics/btr406}{10.1093/bioinformatics/btr406}.}
\seealso{\code{\linkS4class{APResult}},
  \code{\linkS4class{AggExResult}}, \code{\linkS4class{ExClust}},
  \code{\link{apcluster}}, \code{\link{apclusterL}}, 
  \code{\link{aggExCluster}}, \code{\link{cutree-methods}},
  \code{\link{plot-methods}}}
\examples{
## create two Gaussian clouds
cl1 <- cbind(rnorm(50, 0.2, 0.05), rnorm(50, 0.8, 0.06))
cl2 <- cbind(rnorm(50, 0.7, 0.08), rnorm(50, 0.3, 0.05))
x <- rbind(cl1, cl2)

## run affinity propagation using negative squared Euclidean
apres <- apcluster(negDistMat(r=2), x, p=-0.1)

## plot heatmap clustering run
heatmap(apres)

## rerun affinity propagation
## reuse similarity matrix from previous run
apres2 <- apcluster(s=apres@sim, q=0.6)

## plot heatmap of second run
heatmap(apres2, apres@sim)

## with alternate heatmap coloring, alternating color bars, and no dendrograms
heatmap(apres2, apres@sim, Rowv=NA, Colv=NA,
        sideColors=c("darkgreen", "yellowgreen"), col=terrain.colors(12))

## perform agglomerative clustering of affinity propagation clusters
aggres1 <- aggExCluster(apres@sim, apres2)

## plot heatmap
heatmap(cutree(aggres1, 2), apres@sim)

## perform agglomerative clustering of whole data set
aggres2 <- aggExCluster(negDistMat(r=2), x)

## show heatmap along with dendrogram
heatmap(aggres2)
}
\keyword{cluster}
\keyword{methods}
