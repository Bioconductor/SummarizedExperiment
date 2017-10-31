\name{makeSummarizedExperimentFromLoom}

\alias{makeSummarizedExperimentFromLoom}

\title{Make a SummarizedExperiment from a '.loom' hdf5 file}

\description{
    \code{makeSummarizedExperimentFromLoom} represents a '.loom' file as
    a \code{SummarizedExperiment}. The \code{'/matrix'} and
    \code{'/layers'} are represented as \code{HDF5Array} objects; row
    and column attributes are parsed to \code{DataFrame}. Optionally,
    row or column attributes can be specified as row and and column
    names.
}

\usage{
makeSummarizedExperimentFromLoom(file,
                                 rownames_attr = NULL,
                                 colnames_attr = NULL)
}
\arguments{
    \item{file}{
        The path (as a single character string) to the HDF5 file where
        the dataset is located.
    }
    \item{rownames_attr}{
        The name of the row attribute to be used as row names.
    }
    \item{colnames_attr}{
        The name of the column attribute to be used as column names.
    }
}
\value{
    A \link{SummarizedExperiment} object with row and column data and
    one or more assays.
}
\author{
    Martin Morgan
}
\seealso{
    \url{http://loompy.org/loompy-docs/format/index.html} for a
    specification of the .loom format.
}
\examples{
## ---------------------------------------------------------------------
## BASIC EXAMPLE
## ---------------------------------------------------------------------

file <- system.file(
    package="SummarizedExperiment", "extdata", "example.loom"
)
se <- makeSummarizedExperimentFromLoom(file)
se
assay(se)
metadata(se)
}
