\name{apclusterL}
\docType{methods}
\alias{apclusterL}
\alias{apclusterL-methods}
\alias{apclusterL,matrix,missing-method}
\alias{apclusterL,character,ANY-method}
\alias{apclusterL,function,ANY-method}
\title{Leveraged Affinity Propagation}
\description{
   Runs leveraged affinity propagation clustering}
\usage{
\S4method{apclusterL}{matrix,missing}(s, x,
          sel, p=NA, q=NA, maxits=1000, convits=100, lam=0.9,
          includeSim=FALSE, nonoise=FALSE, seed=NA)
\S4method{apclusterL}{character,ANY}(s, x, 
          frac, sweeps, p=NA, q=NA, maxits=1000, convits=100, lam=0.9,
          includeSim=TRUE, nonoise=FALSE, seed=NA, ...)
\S4method{apclusterL}{function,ANY}(s, x,
          frac, sweeps, p=NA, q=NA, maxits=1000, convits=100, lam=0.9,
          includeSim=TRUE, nonoise=FALSE, seed=NA, ...)
}
\arguments{
  \item{s}{an \eqn{l \times length(sel)}{l x length(sel)} similarity 
  	matrix or a similarity function either specified as the name of 
  	a package provided similarity function as character string or a 
  	user provided function object for similarity calculation.
	If \code{s} is supplied as a similarity matrix, the columns
	must correspond to the same sub-selection of samples as
	specified in the \code{sel} argument and must be in the same
	increasing order.
  	For a package- or user-defined similarity function, additional 
  	parameters can be specified as appropriate for the chosen method 
  	and are passed on to the similarity function via the \code{...} 
  	argument (see below). See the package vignette for a non-trivial 
  	example or supplying a user-defined similarity measure.}
  \item{x}{input data to be clustered; if \code{x} is a matrix or data
    frame, rows are interpreted as samples and columns are interpreted 
    as features; apart from matrices or data frames, \code{x} may be 
    any other structured data type that contains multiple data items - 
    provided that an appropriate \code{\link[base:length]{length}} 
    function is available that returns the number of items}
  \item{frac}{fraction of samples that should be used for leveraged 
  	clustering. The similarity matrix will be generated for
  	all samples against a random fraction of the samples as 
  	specified by this parameter.}
  \item{sweeps}{number of sweeps of leveraged clustering performed 
  	with changing randomly selected subset of samples.}
  \item{sel}{selected sample indices; a vector containing the 
  	sample indices of the sample subset used for leveraged 
  	AP clustering in increasing order.}
  \item{p}{input preference; can be a vector that specifies
    individual preferences for each data point. If scalar,
    the same value is used for all data points. If \code{NA},
    exemplar preferences are initialized according to the
    distribution of non-Inf values in \code{s}. How this
    is done is controlled by the parameter \code{q}. See also
	\code{\link{apcluster}}.}
  \item{q}{if \code{p=NA}, exemplar preferences are initialized
    according to the distribution of non-Inf values in \code{s}.
    If \code{q=NA}, exemplar preferences are set to the median
    of non-Inf values in \code{s}. If \code{q} is a value
    between 0 and 1, the sample quantile with threshold
    \code{q} is used, whereas \code{q=0.5} again results in
    the median. See also \code{\link{apcluster}}.}
  \item{maxits}{maximal number of iterations that should be executed}
  \item{convits}{the algorithm terminates if the examplars have not
    changed for \code{convits} iterations}
  \item{lam}{damping factor; should be a value in the range [0.5, 1);
    higher values correspond to heavy damping which may be needed 
    if oscillations occur}
  \item{includeSim}{if \code{TRUE}, the similarity matrix (either computed
    internally or passed via the \code{s} argument) is stored to the
    slot \code{sim} of the returned
    \code{\linkS4class{APResult}} object. The default is \code{FALSE}
    if \code{apclusterL} has been called for a similarity matrix,
    otherwise the default is \code{TRUE}.}
  \item{nonoise}{\code{apcluster} adds a small amount of noise to
    \code{s} to prevent degenerate cases; if \code{TRUE},
    this is disabled}
  \item{seed}{for reproducibility, the seed of the random number
    generator can be set to a fixed value before
    adding noise (see above), if \code{NA}, the seed remains
    unchanged}
  \item{...}{all other arguments are passed to the selected 
    similarity function as they are; note that possible name conflicts between
    arguments of \code{apcluster} and arguments of the similarity
    function may occur; therefore, we recommend to write user-defined
    similarity functions without additional parameters or to use
    closures to fix parameters (such as, in the example below);}
}
\details{Affinity Propagation clusters data using a set of
real-valued pairwise similarities as input. Each cluster
is represented by a representative cluster center (the so-called exemplar). 
The method is iterative and searches for clusters maximizing
an objective function called net similarity.

Leveraged Affinity Propagation reduces dynamic and static load for 
large datasets. Only a subset of the samples are considered
in the clustering process assuming that they provide already enough
information about the cluster structure.

When called with input data and the name of a package provided or a user 
provided similarity function the function selects a random sample subset
according to the \code{frac} parameter, calculates a rectangular 
similarity matrix of all samples against this subset and repeats 
affinity propagation \code{sweep} times. A new sample subset is used 
for each repetition. The clustering result of the sweep with the highest 
net similarity is returned. Any parameters specific to the chosen 
method of similarity calculation can be passed to \code{apcluster} 
in addition to the parameters described above. The similarity matrix 
for the best trial is also returned in the result object when requested 
by the user (argument \code{includeSim}).

When called with a rectangular similarity matrix (which represents a 
column subset of the full similarity matrix) the function performs 
AP clustering on this similarity matrix. The information 
about the selected samples is passed to clustering with the 
parameter \code{sel}. This function is only needed when the user needs full
control of distance calculation or sample subset selection.
   
Apart from minor adaptations and optimizations, the implementation 
of the function \code{apclusterL}
is largely analogous to Frey's and Dueck's Matlab code
(see \url{http://genes.toronto.edu/affinitypropagation/}).}
\value{
  Upon successful completion, both functions returns an
  \code{\linkS4class{APResult}} object.
}
\author{Ulrich Bodenhofer, Andreas Kothmeier & Johannes Palme
\email{apcluster@bioinf.jku.at}}
\references{\url{http://www.bioinf.jku.at/software/apcluster}

Frey, B. J. and Dueck, D. (2007) Clustering by passing messages
between data points. \emph{Science} \bold{315}, 972-976.
DOI: \href{http://dx.doi.org/10.1126/science.1136800}{10.1126/science.1136800}.

Bodenhofer, U., Kothmeier, A., and Hochreiter, S. (2011)
APCluster: an R package for affinity propagation clustering.
\emph{Bioinformatics} \bold{27}, 2463-2464.
DOI: \href{http://dx.doi.org/10.1093/bioinformatics/btr406}{10.1093/bioinformatics/btr406}.
}
\seealso{\code{\linkS4class{APResult}}, \code{\link{show-methods}},
  \code{\link{plot-methods}}, \code{\link{labels-methods}},
  \code{\link{preferenceRange}}, \code{\link{apcluster-methods}},
  \code{\link{apclusterK}}}
\examples{
## create two Gaussian clouds
cl1 <- cbind(rnorm(150, 0.2, 0.05), rnorm(150, 0.8, 0.06))
cl2 <- cbind(rnorm(100, 0.7, 0.08), rnorm(100, 0.3, 0.05))
x <- rbind(cl1, cl2)

## leveraged apcluster
apres <- apclusterL(negDistMat(r=2), x, frac=0.2, sweeps=3, p=-0.2)

## show details of leveraged clustering results
show(apres)

## plot leveraged clustering result
plot(apres, x)

## plot heatmap of clustering result
heatmap(apres)

## show net similarities of single sweeps
apres@netsimLev

## show samples on which best sweep was based
apres@sel
}
\keyword{cluster}
