\name{plot}
\docType{methods}
\alias{plot}
\alias{plot-methods}
\alias{plot,APResult,missing-method}
\alias{plot,ExClust,missing-method}
\alias{plot,ExClust,matrix-method}
\alias{plot,ExClust,data.frame-method}
\alias{plot,AggExResult,missing-method}
\alias{plot,AggExResult,matrix-method}
\alias{plot,AggExResult,data.frame-method}
\title{Plot Clustering Results}
\description{
  Functions for Visualizing Clustering Results
}
\usage{
\S4method{plot}{APResult,missing}(x, y, type=c("netsim", "dpsim", "expref"),
    xlab="# Iterations", ylab="Similarity", ...)
\S4method{plot}{ExClust,matrix}(x, y, connect=TRUE, xlab="", ylab="",
labels=NA, limitNo=15, ...)
\S4method{plot}{ExClust,data.frame}(x, y, connect=TRUE, xlab="",
ylab="", labels=NA, limitNo=15, ...)
\S4method{plot}{AggExResult,missing}(x, y, main="Cluster dendrogram",
    xlab="", ylab="", ticks=4, digits=2, base=0.05, showSamples=FALSE,
    horiz=FALSE, ...)
\S4method{plot}{AggExResult,matrix}(x, y, k=NA, h=NA, ...)
\S4method{plot}{AggExResult,data.frame}(x, y, k=NA, h=NA, ...)
}
\arguments{
  \item{x}{a clustering result object of class
     \code{\linkS4class{APResult}}, \code{\linkS4class{ExClust}},
     or \code{\linkS4class{AggExResult}}}
  \item{y}{a matrix or data frame (see details below)}
  \item{type}{a string or array of strings indicating which
              performance measures should be plotted; valid values are
              \code{"netsim"}, \code{"dpsim"}, and \code{"expref"}
              which can be used in any combination or order;
              all other strings are ignored (for the meaning see
              \code{\linkS4class{APResult}})}
  \item{xlab, ylab}{labels for axes of 2D plots; ignored if \code{y} has
  more than two columns}
  \item{labels}{names used for variables in scatter plot matrix
    (displayed if \code{y} has more than two columns). If \code{NA}
    (default), column names are used. If no column names are available,
    labels such as \code{x[, 2]} are displayed.}
  \item{limitNo}{if the number of columns/features in \code{y} is too large,
    problems may occur when attempting to plot a scatter plot matrix.
    To avoid problems, the \code{plot} method throws an error if the
    number of columns exceeds \code{limitNo}. For special applications,
    users can increase the value (15 by default). If \code{limitNo}
    is set to \code{NA} or any other non-numeric value, the limit is
    ignored entirely. Please note that attempting to plot scatter plot
    matrices with too many features may corrupt the graphics device.
    So users are making changes at their own risk. If plotting of many
    features is necessary, make sure that the graphics device is large
    enough to accommodate the plot (e.g. by using a sufficiently large
    graphics file device).}
  \item{connect}{used only if clustering is plotted on original data,
    ignored otherwise. If \code{connect} is \code{TRUE}, lines are
    drawn connecting exemplars with their cluster members.}
  \item{main}{title of plot}
  \item{ticks}{number of ticks used for the axis on the left side of the
    plot (applies to dendrogram plots only, see below)}
  \item{digits}{number of digits used for the axis tickmarks on the left
    side of the plot (applies to dendrogram plots only, see below)}
  \item{base}{fraction of height used for the very first join; defaults
    to 0.05, i.e. the first join appears at 5\% of the total height of the
    dendrogram.}
  \item{showSamples}{if \code{TRUE}, a complete cluster hierarchy is
    shown, otherwise, in case that \code{x} is a hierarchy of clusters,
    the dendrogram of clusters is shown. For backward compatibility, the
    default is \code{FALSE}.}
  \item{horiz}{if \code{TRUE}, the dendrogram is plotted horizontally
    (analogous to \code{\link[stats:plot.dendrogram]{plot.dendrogram}}).
    The default is \code{FALSE}.}
  \item{k}{level to be selected when plotting a single clustering
    level of cluster hierarchy (i.e. the number of clusters; see
    \code{\link{cutree-methods}})}
  \item{h}{cut-off to be used  when plotting a single clustering
    level of cluster hierarchy (see \code{\link{cutree-methods}})}
  \item{...}{all other arguments are passed to the plotting command that
    are used internally, \code{\link[graphics:plot]{plot}} or
    \code{\link[stats:heatmap]{heatmap}}.}
}
\details{
  If \code{plot} is called for an \code{\linkS4class{APResult}} object
  without specifying the second argument \code{y},
  a plot is created that displays graphs of performance
  measures over execution time of the affinity propagation run.
  This only works if \code{\link{apcluster}} was called with
  \code{details=TRUE}.

  If \code{plot} is called for an \code{\linkS4class{APResult}} object
  along with a matrix or data frame as argument \code{y}, then the dimensions of
  the matrix determine the behavior of \code{plot}:
  \enumerate{
    \item{If the matrix \code{y} has two columns, \code{y} is
      interpreted as the original data set. Then a plot of
      the clustering result superimposed on the original
      data set is created. Each cluster is displayed in a
      different color. The exemplar of each cluster is highlighted
      by a black square. If \code{connect} is \code{TRUE}, lines
      connecting the cluster members to their exemplars are drawn.
      This variant of \code{plot} does not return any value.}
    \item{If \code{y} has more than two columns, clustering results are
      superimposed in a sort of scatter plot matrix. The variant
      that \code{y} is interpreted as similarity matrix if it is
      quadratic has been removed in version 1.3.2. Use
      \code{\link{heatmap}} instead.}
    \item{If \code{y} has only one column, an error is displayed.}
  }
  
  If \code{plot} is called for an \code{\linkS4class{ExClust}} object
  along with a matrix or data frame as argument \code{y}, then
  \code{plot} behaves exactly the same as described in the previous
  paragraph.
  
  If \code{plot} is called for an \code{\linkS4class{AggExResult}} object
  without specifying the second argument \code{y}, then a dendrogram
  plot is drawn. This variant returns an invisible
  \code{\link{dendrogram}} object. The \code{showSamples} argument
  determines whether a complete dendrogram or a dendrogram of clusters
  is plotted (see above). If the option \code{horiz=TRUE} is used, the
  dendrogram is rotated. Note that, in this case, the margin to the
  right of the plot may not be wide enough to accommodate long
  cluster/sample labels. In such a case, the figure margins have to
  be widened before \code{plot} is called.

  If \code{plot} is called for an \code{\linkS4class{AggExResult}} object
  along with a matrix or data frame \code{y}, \code{y} is
  again interpreted
  as original data set. If one of the two arguments \code{k} or
  \code{h} is present, a clustering is cut out from the cluster hierarchy 
  using \code{\link{cutree}} and this clustering is displayed with the
  original data set as described above. This variant of
  \code{plot} returns an invisible \code{\linkS4class{ExClust}} object
  containing the extracted clustering.
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
DOI: \href{http://dx.doi.org/10.1093/bioinformatics/btr406}{10.1093/bioinformatics/btr406}.}
\seealso{\code{\linkS4class{APResult}},
  \code{\linkS4class{AggExResult}}, \code{\linkS4class{ExClust}},
  \code{\link{heatmap-methods}}, \code{\link{apcluster}}, 
  \code{\link{apclusterL}}, \code{\link{aggExCluster}}, 
  \code{\link{cutree-methods}}}
\examples{
## create two Gaussian clouds
cl1 <- cbind(rnorm(50, 0.2, 0.05), rnorm(50, 0.8, 0.06))
cl2 <- cbind(rnorm(50, 0.7, 0.08), rnorm(50, 0.3, 0.05))
x <- rbind(cl1, cl2)

## run affinity propagation
apres <- apcluster(negDistMat(r=2), x, q=0.7, details=TRUE)

## plot information about clustering run
plot(apres)

## plot clustering result
plot(apres, x)

## perform agglomerative clustering of affinity propagation clusters
aggres1 <- aggExCluster(x=apres)

## show dendrograms
plot(aggres1)
plot(aggres1, showSamples=TRUE)

## show clustering result for 4 clusters
plot(aggres1, x, k=4)

## perform agglomerative clustering of whole data set
aggres2 <- aggExCluster(negDistMat(r=2), x)

## show dendrogram
plot(aggres2)

## show heatmap along with dendrogram
heatmap(aggres2)

## show clustering result for 2 clusters
plot(aggres2, x, k=2)

## cluster iris data set
data(iris)
apIris <- apcluster(negDistMat(r=2), iris, q=0)
plot(apIris, iris)
}
\keyword{cluster}
\keyword{methods}
