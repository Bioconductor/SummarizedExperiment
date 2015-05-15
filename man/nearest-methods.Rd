\name{nearest-methods}

\alias{nearest-methods}

\alias{precede}
\alias{precede,RangedSummarizedExperiment,ANY-method}
\alias{precede,ANY,RangedSummarizedExperiment-method}
\alias{precede,RangedSummarizedExperiment,RangedSummarizedExperiment-method}

\alias{follow}
\alias{follow,RangedSummarizedExperiment,ANY-method}
\alias{follow,ANY,RangedSummarizedExperiment-method}
\alias{follow,RangedSummarizedExperiment,RangedSummarizedExperiment-method}

\alias{nearest}
\alias{nearest,RangedSummarizedExperiment,ANY-method}
\alias{nearest,ANY,RangedSummarizedExperiment-method}
\alias{nearest,RangedSummarizedExperiment,RangedSummarizedExperiment-method}

\alias{distance}
\alias{distance,RangedSummarizedExperiment,ANY-method}
\alias{distance,ANY,RangedSummarizedExperiment-method}
\alias{distance,RangedSummarizedExperiment,RangedSummarizedExperiment-method}

\alias{distanceToNearest}
\alias{distanceToNearest,RangedSummarizedExperiment,ANY-method}
\alias{distanceToNearest,ANY,RangedSummarizedExperiment-method}
\alias{distanceToNearest,RangedSummarizedExperiment,RangedSummarizedExperiment-method}


\title{Finding the nearest range neighbor in RangedSummarizedExperiment objects}

\description{
  This man page documents the \code{nearest} methods and family (i.e.
  \code{precede}, \code{follow}, \code{distance}, and \code{distanceToNearest}
  methods) for \link{RangedSummarizedExperiment} objects.
}

\usage{
\S4method{precede}{RangedSummarizedExperiment,ANY}(x, subject, select=c("arbitrary", "all"),
        ignore.strand=FALSE)
\S4method{precede}{ANY,RangedSummarizedExperiment}(x, subject, select=c("arbitrary", "all"),
        ignore.strand=FALSE)

\S4method{follow}{RangedSummarizedExperiment,ANY}(x, subject, select=c("arbitrary", "all"),
        ignore.strand=FALSE)
\S4method{follow}{ANY,RangedSummarizedExperiment}(x, subject, select=c("arbitrary", "all"),
        ignore.strand=FALSE)

\S4method{nearest}{RangedSummarizedExperiment,ANY}(x, subject, select=c("arbitrary", "all"),
        algorithm=c("nclist", "intervaltree"), ignore.strand=FALSE)
\S4method{nearest}{ANY,RangedSummarizedExperiment}(x, subject, select=c("arbitrary", "all"),
        algorithm=c("nclist", "intervaltree"), ignore.strand=FALSE)

\S4method{distance}{RangedSummarizedExperiment,ANY}(x, y, ignore.strand=FALSE, ...)
\S4method{distance}{ANY,RangedSummarizedExperiment}(x, y, ignore.strand=FALSE, ...)

\S4method{distanceToNearest}{RangedSummarizedExperiment,ANY}(x, subject,
        algorithm=c("nclist", "intervaltree"), ignore.strand=FALSE, ...)
\S4method{distanceToNearest}{ANY,RangedSummarizedExperiment}(x, subject,
        algorithm=c("nclist", "intervaltree"), ignore.strand=FALSE, ...)
}

\arguments{
  \item{x, subject}{
    One of these two arguments must be a \link{RangedSummarizedExperiment}
    object.
  }
  \item{select, ignore.strand, algorithm}{
    See \code{?\link[GenomicRanges]{nearest}} in the \pkg{GenomicRanges}
    package.
  }
  \item{y}{
    For the \code{distance} methods, one of \code{x} or \code{y} must be a
    \link{RangedSummarizedExperiment} object.
  }
  \item{...}{Additional arguments for methods.}
}

\details{
  These methods operate on the \code{rowRanges} component of the
  \link{RangedSummarizedExperiment} object, which can be a
  \link[GenomicRanges]{GenomicRanges} or \link[GenomicRanges]{GRangesList}
  object.

  More precisely, if any of the above functions is passed a
  \link{RangedSummarizedExperiment} object thru the \code{x}, \code{subject},
  and/or \code{y} argument, then it behaves as if \code{rowRanges(x)},
  \code{rowRanges(subject)}, and/or \code{rowRanges(y)} had been passed
  instead.

  See \code{?\link[GenomicRanges]{nearest}} in the \pkg{GenomicRanges}
  package for the details of how \code{nearest} and family operate on
  \link[GenomicRanges]{GenomicRanges} and \link[GenomicRanges]{GRangesList}
  objects.
}

\value{
  See \code{?\link[GenomicRanges]{nearest}} in the \pkg{GenomicRanges}
  package.
}

\seealso{
  \itemize{
    \item \link{RangedSummarizedExperiment} objects.

    \item The \link[GenomicRanges]{nearest} man page in the
          \pkg{GenomicRanges} package where the \code{nearest} family
          of methods for \link[GenomicRanges]{GenomicRanges} and
          \link[GenomicRanges]{GRangesList} objects is documented.
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
se1 <- shift(se0, 100)

res <- nearest(se0, se1)
res
stopifnot(identical(res, nearest(rowRanges(se0), rowRanges(se1))))
stopifnot(identical(res, nearest(se0, rowRanges(se1))))
stopifnot(identical(res, nearest(rowRanges(se0), se1)))

res <- nearest(se0)  # missing subject
res
stopifnot(identical(res, nearest(rowRanges(se0))))

hits <- nearest(se0, se1, select="all")
hits
stopifnot(identical(
  hits,
  nearest(rowRanges(se0), rowRanges(se1), select="all")
))
stopifnot(identical(
  hits,
  nearest(se0, rowRanges(se1), select="all")
))
stopifnot(identical(
  hits,
  nearest(rowRanges(se0), se1, select="all")
))
}

\keyword{methods} 
\keyword{utilities}
