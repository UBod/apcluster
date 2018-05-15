\name{ch22Promoters}
\alias{ch22Promoters}
\title{ch22Promoters Data Set}
\description{Toy data set consisting of promoter regions of 150
  random genes from the human chromosome 22 (obsolete)}

\usage{ch22Promoters}
\format{A data set with 150 DNA sequences. Each string is
  a nucleotide sequence that corresponds to
  the promoter region of a gene from the human chromosome no. 22
  (according to the human genome assembly hg18). The sequences
  start 999 bases upstream of the transcription start site (TSS)
  and end with the TSS itself.
  The names attribute contains the RefSeq IDs of the genes.

  In previous version of the \pkg{apcluster} package, this was an
  R object that can be loaded via \code{data(ch22Promoters)}. For
  better compatibility with the \pkg{kebabs} package, the data set
  has been moved to a plain text file (in FASTA format)
  that can be loaded from \code{inst/examples} (see examples below).
}
\examples{
## if Biostrings package is available
if (requireNamespace("Biostrings", quietly=TRUE))
{
    ## load data set
    filepath <- system.file("examples", "ch22Promoters.fasta",
                            package="apcluster")
    ch22Promoters <- Biostrings::readDNAStringSet(filepath)

    ## display sequences
    ch22Promoters
}
}
\keyword{datasets}

