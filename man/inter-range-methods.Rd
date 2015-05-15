\name{inter-range-methods}

\alias{inter-range-methods}

\alias{isDisjoint}
\alias{isDisjoint,RangedSummarizedExperiment-method}

\alias{disjointBins}
\alias{disjointBins,RangedSummarizedExperiment-method}


\title{Inter range transformations of a RangedSummarizedExperiment object}

\description{
  This man page documents \emph{inter range transformations} of a
  \link{RangedSummarizedExperiment} object. These transformations
  operate on the \code{rowRanges} component of the object, which can be a
  \link[GenomicRanges]{GenomicRanges} or \link[GenomicRanges]{GRangesList}
  object.

  See \code{?\link[GenomicRanges]{isDisjoint}} in the \pkg{GenomicRanges}
  package for the details of how these transformations operate on a
  \link[GenomicRanges]{GenomicRanges} or \link[GenomicRanges]{GRangesList}
  object.
}

\usage{
\S4method{isDisjoint}{RangedSummarizedExperiment}(x, ignore.strand=FALSE)

\S4method{disjointBins}{RangedSummarizedExperiment}(x, ignore.strand=FALSE)
}

\arguments{
  \item{x}{
    A \link{RangedSummarizedExperiment} object.
  }
  \item{ignore.strand}{
    See \code{?\link[GenomicRanges]{isDisjoint}} in the
    \pkg{GenomicRanges} package.
  }
}

\details{
  Any of these functions performs the following transformation on
  \link{RangedSummarizedExperiment} object \code{x}:
\preformatted{f(rowRanges(x), ...)
}
  where \code{f} is the name of the function and \code{...} any additional
  arguments passed to it.
}

\value{
  See \code{?\link[GenomicRanges]{isDisjoint}} in the
  \pkg{GenomicRanges} package.
}

\seealso{
  \itemize{
    \item \link{RangedSummarizedExperiment} objects.

    \item The \link[GenomicRanges]{isDisjoint} man page in the
          \pkg{GenomicRanges} package where \emph{inter range transformations}
          of a \link[GenomicRanges]{GenomicRanges} or
          \link[GenomicRanges]{GRangesList} object are documented.
  }
}

\examples{
nrows <- 20; ncols <- 6
counts <- matrix(runif(nrows * ncols, 1, 1e4), nrows)
rowRanges <- GRanges(rep(c("chr1", "chr2"), c(5, 15)),
                     IRanges(sample(1000L, 20), width=100),
                     strand=Rle(c("+", "-"), c(12, 8)))
colData <- DataFrame(Treatment=rep(c("ChIP", "Input"), 3),
                     row.names=LETTERS[1:6])
se0 <- SummarizedExperiment(assays=SimpleList(counts=counts),
                            rowRanges=rowRanges, colData=colData)
se1 <- shift(se0, 99*start(se0))

isDisjoint(se0)  # FALSE
isDisjoint(se1)  # TRUE

bins0 <- disjointBins(se0)
bins0
stopifnot(identical(bins0, disjointBins(rowRanges(se0))))

bins1 <- disjointBins(se1)
bins1
stopifnot(all(bins1 == bins1[1]))
}

\keyword{methods} 
\keyword{utilities}
