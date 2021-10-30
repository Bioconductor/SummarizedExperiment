\name{readKallisto}
\alias{readKallisto}
\alias{readKallistoBootstrap}
\alias{KALLISTO_ASSAYS}

\title{
  Input kallisto or kallisto bootstrap results.
}

\description{
  WARNING: \code{readKallisto()} is defunct in BioC 3.15. Please use
  \code{tximeta()} from the \pkg{tximeta} package instead.

  \code{readKallisto} inputs several kallisto output files into a single
  \code{SummarizedExperiment} instance, with rows corresponding to
  estimated transcript abundance and columns to
  samples. \code{readKallistoBootstrap} inputs kallisto bootstrap
  replicates of a single sample into a matrix of transcript x bootstrap
  abundance estimates.
}

\usage{
readKallisto(files,
    json = file.path(dirname(files), "run_info.json"), 
    h5 = any(grepl("\\\\.h5$", files)), what = KALLISTO_ASSAYS,
    as = c("SummarizedExperiment", "list", "matrix"))

readKallistoBootstrap(file, i, j)
}

\arguments{

  \item{files}{character() paths to kallisto \sQuote{abundance.tsv}
    output files. The assumption is that files are organized in the way
    implied by kallisto, with each sample in a distinct directory, and
    the directory containing files abundance.tsv, run_info.json, and
    perhaps abundance.h5.}

  \item{json}{character() vector of the same length as \code{files}
    specifying the location of JSON files produced by kallisto and
    containing information on the run. The default assumes that json
    files are in the same directory as the corresponding abundance
    file.}

  \item{h5}{character() vector of the same length as \code{files}
    specifying the location of HDF5 files produced by kallisto and
    containing bootstrap estimates. The default assumes that HDF5 files
    are in the same directory as the corresponding abundance file.}

  \item{what}{character() vector of kallisto per-sample outputs to be
    input. See KALLISTO_ASSAYS for available values.}

  \item{as}{character(1) specifying the output format. See \code{Value}
    for additional detail.}

  \item{file}{character(1) path to a single HDF5 output file.}

  \item{i, j}{integer() vector of row (\code{i}) and column (\code{j})
    indexes to input.}

}

\value{

  A \code{SummarizedExperiment}, \code{list}, or \code{matrix},
  depending on the value of argument \code{as}; by default a
  \code{SummarizedExperiment}. The \code{as="SummarizedExperiment"}
  \code{rowData(se)} the length of each transcript;
  \code{colData(se)} includes summary information on each sample,
  including the number of targets and bootstraps, the kallisto and index
  version, the start time and operating system call used to create the
  file. \code{assays()} contains one or more transcript x sample
  matrices of parameters estimated by kallisto (see
  \code{KALLISTO_ASSAYS}).

  \code{as="list"} return value contains information simillar to
  \code{SummarizedExperiment} with row, column and assay data as
  elements of the list without coordination of row and column
  annotations into an integrated data container. \code{as="matrix"}
  returns the specified assay as a simple \emph{R} matrix.

}

\references{
  \url{http://pachterlab.github.io/kallisto} software for quantifying
  transcript abundance.
}

\author{
  Martin Morgan \url{martin.morgan@roswellpark.org}
}

\examples{
\dontrun{
outputs <- system.file(package="SummarizedExperiment", "extdata",
    "kallisto")
files <- dir(outputs, pattern="abundance.tsv", full=TRUE, recursive=TRUE)
stopifnot(all(file.exists(files)))

## default: input 'est_counts'
(se <- readKallisto(files, as="SummarizedExperiment"))
str(readKallisto(files, as="list"))
str(readKallisto(files, as="matrix"))

## available assays
KALLISTO_ASSAYS
## one or more assay
readKallisto(files, what=c("tpm", "eff_length"))

## alternatively: read hdf5 files
files <- sub(".tsv", ".h5", files, fixed=TRUE)
readKallisto(files)

## input all bootstraps
xx <- readKallistoBootstrap(files[1])
ridx <- head(which(rowSums(xx) != 0), 3)
cidx <- c(1:5, 96:100)
xx[ridx, cidx]

## selective input of rows (transcripts) and/or bootstraps
readKallistoBootstrap(files[1], i=c(ridx, rev(ridx)), j=cidx)
}
}

\keyword{file}
