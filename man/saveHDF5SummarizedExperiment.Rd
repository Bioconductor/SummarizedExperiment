\name{saveHDF5SummarizedExperiment}

\alias{saveHDF5SummarizedExperiment}
\alias{loadHDF5SummarizedExperiment}

\title{Save/load a HDF5-based SummarizedExperiment object}

\description{
  \code{saveHDF5SummarizedExperiment} and \code{loadHDF5SummarizedExperiment}
  can be used to save/load a HDF5-based \link{SummarizedExperiment} object
  to/from disk.
}

\usage{
saveHDF5SummarizedExperiment(x, dir="myse")
loadHDF5SummarizedExperiment(dir="myse")
}

\arguments{
  \item{x}{
    A \link{SummarizedExperiment} object.
  }
  \item{dir}{
    The path (as a single string) to the directory where to save/load the
    HDF5-based \link{SummarizedExperiment} object.
    The directory should not exist yet and will be created when saving.
  }
}

\details{
  These functions use functionalities from the \pkg{rhdf5} and \pkg{HDF5Array}
  packages internally so they will fail if these packages are not installed.
}

\author{Hervé Pagès}

\seealso{
  \itemize{
    \item \link{SummarizedExperiment} and \link{RangedSummarizedExperiment}
          objects.

    \item The \code{\link[HDF5Array]{writeHDF5Array}} function in the
          \pkg{HDF5Array} package, which \code{saveHDF5SummarizedExperiment}
          is based on.
  }
}

\examples{
nrows <- 200; ncols <- 6
counts <- matrix(runif(nrows * ncols, 1, 1e4), nrows)
colData <- DataFrame(Treatment=rep(c("ChIP", "Input"), 3),
                     row.names=LETTERS[1:6])
se0 <- SummarizedExperiment(assays=SimpleList(counts=counts),
                            colData=colData)

dir <- sub("file", "myse", tempfile())
saveHDF5SummarizedExperiment(se0, dir)
se1 <- loadHDF5SummarizedExperiment(dir)
stopifnot(is(assay(se1, withDimnames=FALSE), "HDF5Matrix"))
}
