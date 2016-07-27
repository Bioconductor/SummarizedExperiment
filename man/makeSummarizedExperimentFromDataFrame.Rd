\name{makeSummarizedExperimentFromDataFrame}

\alias{makeSummarizedExperimentFromDataFrame}

\title{Make a RangedSummarizedExperiment from a data.frame or DataFrame}

\description{
    \code{makeSummarizedExperimentFromDataFrame} uses \code{data.frame}
    or \code{DataFrame} column names to create a \link{GRanges} object for the
    \code{rowRanges} of the resulting \link{SummarizedExperiment} object.
    It requires that non-range data columns be coercible into a \code{numeric}
    \code{matrix} for the \link{SummarizedExperiment} constructor. All columns
    that are not part of the row ranges attribute are assumed to be experiment
    data; thus, keeping metadata columns will not be supported. Note that this
    function only returns \link{SummarizedExperiment} objects with a single
    assay.

    If metadata columns are to be kept, one can first construct the row ranges
    attribute by using the \link{makeGRangesFromDataFrame} function and
    subsequently creating the \link{SummarizedExperiment}.
}

\usage{
makeSummarizedExperimentFromDataFrame(df,
                                    ...,
                                    seqinfo = NULL,
                                    starts.in.df.are.0based = FALSE)
}
\arguments{
    \item{df}{
        A data.frame or \link[S4Vectors]{DataFrame} object. If not, then
        the function first tries to turn \code{df} into a data frame with
        \code{as.data.frame(df)}.
    }
    \item{...}{
        Additional arguments passed on to \link{makeGRangesFromDataFrame}
    }
    \item{seqinfo}{
        Either \code{NULL}, or a \link{Seqinfo} object, or a character vector
        of seqlevels, or a named numeric vector of sequence lengths.
        When not \code{NULL}, it must be compatible with the genomic ranges
        in \code{df} i.e. it must include at least the sequence levels
        represented in \code{df}.
    }
    \item{starts.in.df.are.0based}{
        \code{TRUE} or \code{FALSE} (the default).
        If \code{TRUE}, then the start positions of the genomic ranges in
        \code{df} are considered to be \emph{0-based} and are converted to
        \emph{1-based} in the returned \link{GRanges} object.
        This feature is intended to make it more convenient to handle input
        that contains data obtained from resources using the "0-based
        start" convention. A notorious example of such resource is the UCSC
        Table Browser (\url{http://genome.ucsc.edu/cgi-bin/hgTables}).
    }
}
\value{
A \link{RangedSummarizedExperiment} object with rowRanges and a single assay
}
\author{
    M. Ramos
}
\seealso{
  \itemize{
      \item \link{makeGRangesFromDataFrame} 
  }
}

\examples{
## ---------------------------------------------------------------------
## BASIC EXAMPLES
## ---------------------------------------------------------------------

# Note that rownames of the data.frame are also rownames of the result
df <- data.frame(chr="chr2", start = 11:15, end = 12:16,
                 strand = c("+", "-", "+", "*", "."), expr0 = 3:7,
                 expr1 = 8:12, expr2 = 12:16,
                 row.names = paste0("GENE", letters[5:1]))
df

exRSE <- makeSummarizedExperimentFromDataFrame(df)

exRSE

assay(exRSE)

rowRanges(exRSE)
}
