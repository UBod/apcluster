\name{similarities}
\alias{similarities}
\alias{negDistMat}
\alias{expSimMat}
\alias{linSimMat}
\alias{corSimMat}
\alias{linKernel}
\title{Methods for Computing Similarity Matrices}
\description{
   Compute similarity matrices from data set
}
\usage{
negDistMat(x, sel=NA, r=1, method="euclidean", p=2)
expSimMat(x, sel=NA, r=2, w=1, method="euclidean", p=2)
linSimMat(x, sel=NA, w=1, method="euclidean", p=2)
corSimMat(x, sel=NA, r=1, signed=TRUE, method="pearson")
linKernel(x, sel=NA, normalize=FALSE)
}
\arguments{
  \item{x}{input data to be clustered; if \code{x} is a vector, it is
    interpreted as a list of scalar values that are to be clustered;
    if \code{x} is a matrix or data frame, rows are interpreted as
    samples and columns are interpreted as features; in the case
    that \code{x} is a data frame, only numerical columns/features are
    taken into account, whereas categorical features are neglected. If
    \code{x} is missing, all functions return a function that
    can be used as similarity measure, in particular, as
    \code{s} argument for \code{\link{apclusterL}}.}
  \item{sel}{selected samples subset; vector of row indices
  	      for x in increasing order (see details below)}
  \item{r}{exponent (see details below)}
  \item{w}{radius (see details below)}
  \item{signed}{take sign of correlation into account (see details below)}
  \item{normalize}{see details below}
  \item{method}{type of distance measure to be used; for \code{negDistMat},
    \code{expSimMat} and \code{linSimMat}, this argument is analogous to the
    \code{method} argument of \code{\link[stats:dist]{dist}}.
    For \code{corSimMat}, this argument is analogous to the \code{method}
    argument of \code{\link[stats:cor]{cor}}.}
  \item{p}{exponent for Minkowski distance; only used for
    \code{method="minkowski"}, otherwise ignored. See
    \code{\link[stats:dist]{dist}}.}
}
\details{\code{negDistMat} creates a square matrix of mutual
pairwise similarities of data vectors as negative distances. The
argument \code{r} (default is 1) is used to transform the resulting
distances by computing the r-th power (use \code{r=2} to obtain
negative squared distances as in Frey's and Dueck's demos), i.e.,
given a distance d, the resulting similarity is computed as
\eqn{s=-d^r}. With the parameter \code{sel} a subset of samples
can be specified for distance calculation. In this case not the 
full distance matrix is computed but a rectangular similarity matrix
of all samples (rows) against the subset (cols) as needed for 
leveraged clustering. Internally, the computation of distances is 
done using an internal method derived from 
\code{\link[stats:dist]{dist}}. All options of this function except 
\code{diag} and \code{upper} can be used, especially \code{method} 
which allows for selecting different distance measures.
Note that, since version 1.4.4. of the package, there is an additional
method \code{"discrepancy"} that implements Weyl's discrepancy measure.

\code{expSimMat} computes similarities in a way similar to
\code{negDistMat}, but the transformation of distances to similarities
is done in the following way:
\deqn{s=\exp\left(-\left(\frac{d}{w}\right)^r\right)}{s=exp(-(d/w)^r)}
The parameter \code{sel} allows the creation of a rectangular
similarity matrix. As above, r is an exponent. The parameter w controls 
the speed of descent. \code{r=2} in conjunction with Euclidean 
distances corresponds to the well-known Gaussian/RBF kernel, 
whereas \code{r=1} corresponds to the Laplace kernel. Note that these
similarity measures can also be understood as fuzzy equality relations.

\code{linSimMat} provides another way of transforming distances
into similarities by applying the following transformation to a
distance d:
\deqn{s=\max\left(0,1-\frac{d}{w}\right)}{s=max(0,1-d/w)}
Thw parameter \code{sel} is used again for creation of a rectangular
similarity matrix. Here \code{w} corresponds to a maximal radius of 
interest. Note that this is a fuzzy equality relation with respect to 
the Lukasiewicz t-norm.

Unlike the above three functions, \code{linKernel} computes pairwise
similarities as scalar products of data vectors, i.e. it corresponds,
as the name suggests, to the \dQuote{linear kernel}. Use parameter 
\code{sel} to compute only a submatrix of the full kernel matrix as 
described above. If \code{normalize=TRUE}, the values are scaled to 
the unit sphere in the following way (for two samples \code{x} and 
\code{y}:
\deqn{s=\frac{\vec{x}^T\vec{y}}{\|\vec{x}\| \|\vec{y}\|}}{s=(x^T y)/(|x| |y|)}

The function \code{corSimMat} computes pairwise similarities as
correlations. It uses \code{link[stats:cor]{cor}} internally.
The \code{method} argument is passed on to \code{link[stats:cor]{cor}}.
The argument \code{r} serves as an exponent with which the correlations
can be transformed. If \code{signed=TRUE} (default), negative correlations are
taken into account, i.e. two samples are maximally dissimilar if they
are negatively correlated. If \code{signed=FALSE}, similarities are
computed as absolute values of correlations, i.e. two samples are
maximally similar if they are positively or negatively correlated and
the two samples are maximally dissimilar if they are uncorrelated.

Note that the naming of the argument \code{p} has been chosen for
consistency with \code{\link[stats:dist]{dist}} and previous versions 
of the package. When using leveraged AP in
conjunction with the Minkowski distance, this leads to conflicts with
the input preference parameter \code{p} of
\code{\link{apclusterL}}. In order to avoid that, use the above
functions without \code{x} argument to create a custom similarity
measure with fixed parameter \code{p} (see example below).
}
\value{
   All functions listed above return square or rectangular matrices 
   of similarities.
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

Micchelli, C. A. (1986) Interpolation of scattered data: distance
matrices and conditionally positive definite functions.
\emph{Constr. Approx.} \bold{2}, 11-20.

De Baets, B. and Mesiar, R. (1997) Pseudo-metrics and T-equivalences.
\emph{J. Fuzzy Math.} \bold{5}, 471-481.

Bauer, P., Bodenhofer, U., and Klement, E. P. (1996)
A fuzzy algorithm for pixel classification based on the discrepancy
norm.
In \emph{Proc. 5th IEEE Int. Conf. on Fuzzy Systems}, volume III, pages
2007--2012, New Orleans, LA.
DOI: \href{http://dx.doi.org/10.1109/FUZZY.1996.552744}{10.1109/FUZZY.1996.552744}.
}
\seealso{\code{\link[stats:dist]{dist}},
\code{\link{apcluster}}, \code{\link{apclusterL}}}
\examples{
## create two Gaussian clouds
cl1 <- cbind(rnorm(100, 0.2, 0.05), rnorm(100, 0.8, 0.06))
cl2 <- cbind(rnorm(100, 0.7, 0.08), rnorm(100, 0.3, 0.05))
x <- rbind(cl1, cl2)

## create negative distance matrix (default Euclidean)
sim1 <- negDistMat(x)

## compute similarities as squared negative distances
## (in accordance with Frey's and Dueck's demos)
sim2 <- negDistMat(x, r=2)

## compute RBF kernel
sim3 <- expSimMat(x, r=2)

## compute similarities as squared negative distances
## all samples versus a randomly chosen subset 
## of 50 samples (for leveraged AP clustering)
sel <- sort(sample(1:nrow(x), nrow(x)*0.25)) 
sim4 <- negDistMat(x, sel, r=2)


## example of leveraged AP using Minkowski distance with non-default
## parameter p
cl1 <- cbind(rnorm(150, 0.2, 0.05), rnorm(150, 0.8, 0.06))
cl2 <- cbind(rnorm(100, 0.7, 0.08), rnorm(100, 0.3, 0.05))
x <- rbind(cl1, cl2)

apres <- apclusterL(s=negDistMat(method="minkowski", p=2.5, r=2),
                       x, frac=0.2, sweeps=3, p=-0.2)
show(apres)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{cluster}

