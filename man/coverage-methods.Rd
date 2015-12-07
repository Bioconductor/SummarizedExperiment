\name{coverage-methods}

\alias{coverage-methods}

\alias{coverage}
\alias{coverage,RangedSummarizedExperiment-method}


\title{Coverage of a RangedSummarizedExperiment object}

\description{
  This man page documents the \code{coverage} method for
  \link{RangedSummarizedExperiment} objects.
}

\usage{
\S4method{coverage}{RangedSummarizedExperiment}(x, shift=0L, width=NULL, weight=1L,
            method=c("auto", "sort", "hash"))
}

\arguments{
  \item{x}{
    A \link{RangedSummarizedExperiment} object.
  }
  \item{shift, width, weight, method}{
    See \code{?\link[GenomicRanges]{coverage}} in the \pkg{GenomicRanges}
    package.
  }
}

\details{
  This method operates on the \code{rowRanges} component of the
  \link{RangedSummarizedExperiment} object, which can be a
  \link[GenomicRanges]{GenomicRanges} or \link[GenomicRanges]{GRangesList}
  object.

  More precisely, on \link{RangedSummarizedExperiment} object \code{x},
  \code{coverage(x, ...)} is equivalent to \code{coverage(rowRanges(x), ...)}.

  See \code{?\link[GenomicRanges]{coverage}} in the \pkg{GenomicRanges}
  package for the details of how \code{coverage} operates on a
  \link[GenomicRanges]{GenomicRanges} or \link[GenomicRanges]{GRangesList}
  object.
}

\value{
  See \code{?\link[GenomicRanges]{coverage}} in the \pkg{GenomicRanges}
  package.
}

\seealso{
  \itemize{
    \item \link{RangedSummarizedExperiment} objects.

    \item The \link[GenomicRanges]{coverage} man page in the
          \pkg{GenomicRanges} package where the \code{coverage} methods
          for \link[GenomicRanges]{GenomicRanges} and
          \link[GenomicRanges]{GRangesList} objects are documented.
  }
}

\examples{
nrows <- 20; ncols <- 6
counts <- matrix(runif(nrows * ncols, 1, 1e4), nrows)
rowRanges <- GRanges(rep(c("chr1", "chr2"), c(5, 15)),
                     IRanges(sample(1000L, 20), width=100),
                     strand=Rle(c("+", "-"), c(12, 8)),
                     seqlengths=c(chr1=1800, chr2=1300))
colData <- DataFrame(Treatment=rep(c("ChIP", "Input"), 3),
                     row.names=LETTERS[1:6])
rse <- SummarizedExperiment(assays=SimpleList(counts=counts),
                            rowRanges=rowRanges, colData=colData)

cvg <- coverage(rse)
cvg
stopifnot(identical(cvg, coverage(rowRanges(rse))))
}

\keyword{methods}
\keyword{utilities}
