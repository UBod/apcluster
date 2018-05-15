\name{apclusterDemo}
\alias{apclusterDemo}
\title{Affinity Propagation Demo}
\description{
  Runs affinity propagation demo for randomly generated data set
  according to Frey and Dueck
}
\usage{
apclusterDemo(l=100, d=2, seed=NA, ...)
}
\arguments{
  \item{l}{number of data points to be generated}
  \item{d}{dimension of data to be created}
  \item{seed}{for reproducibility, the seed of the random number
              generator can be set to a fixed value; if \code{NA},
              the seed remains unchanged}
  \item{...}{all other arguments are passed on to 
             \code{\link{apcluster}}}
}
\details{\code{apclusterDemo} creates \code{l} \code{d}-dimensional
data points that are uniformly distributed in \eqn{[0,1]^d}. Affinity
propagation is executed for this data set with default parameters.
Alternative settings can be passed to \code{\link{apcluster}} with
additional arguments. After completion of affinity propagation,
the results are shown and the performance measures are plotted.

This function corresponds to the demo function in the original
Matlab code of Frey and Dueck. We warn the user, however, that
uniformly distributed data are not necessarily ideal for demonstrating
clustering, as there can never be real clusters in uniformly
distributed data - all clusters found must be random artefacts.
}
\value{
  Upon successful completion, the function returns an invisible list
  with three components. The first is the data set that has been
  created, the second is the similarity matrix, and the third is an
  \code{\linkS4class{APResult}} object with the clustering results (see
  examples below).
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
\seealso{\code{\linkS4class{APResult}}, \code{\link{plot-methods}},
  \code{\link{apcluster}}, \code{\link{apclusterL}}}
\examples{
## create random data set and run affinity propagation
apd <- apclusterDemo()

## plot clustering result along with data set
plot(apd[[3]], apd[[1]])
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{cluster}

