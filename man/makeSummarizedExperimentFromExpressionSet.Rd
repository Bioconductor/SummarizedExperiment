\name{makeSummarizedExperimentFromExpressionSet}

\alias{makeSummarizedExperimentFromExpressionSet}
\alias{naiveRangeMapper}
\alias{probeRangeMapper}
\alias{geneRangeMapper}
\alias{coerce,ExpressionSet,RangedSummarizedExperiment-method}
\alias{coerce,RangedSummarizedExperiment,ExpressionSet-method}


\title{Make a RangedSummarizedExperiment object from an ExpressionSet}

\description{
  \code{makeSummarizedExperimentFromExpressionSet} takes an
  \link[Biobase]{ExpressionSet} object as input and a \emph{range mapping
  function} that maps the features to ranges.  It then returns a
  \link{RangedSummarizedExperiment} object.
}

\usage{
makeSummarizedExperimentFromExpressionSet(from,
                                          mapFun=naiveRangeMapper)

## range mapping functions
naiveRangeMapper(from)
probeRangeMapper(from)
geneRangeMapper(txDbPackage, key = "ENTREZID")
}

\arguments{
  \item{from}{
    An \link[Biobase]{ExpressionSet} object.
  }
  \item{mapFun}{
    A function which takes an \link[Biobase]{ExpressionSet} object and
    returns a \link{GRanges}, or \link{GRangesList} object which
    corresponds to the genomic ranges used in the ExpressionSet.  The
    \link[base]{rownames} of the returned \link[GenomicRanges]{GRanges}
    are used to match the \link[Biobase]{featureNames} of the
    \link[Biobase]{ExpressionSet}.

    The \code{naiveRangeMapper} function is used by default.
  }
  \item{txDbPackage}{
    A character string with the Transcript Database to use for the mapping.
  }
  \item{key}{
    A character string with the Gene key to use for the mapping.
  }
}

\value{
  \code{makeSummarizedExperimentFromExpressionSet} returns a
  \link{RangedSummarizedExperiment} object that corresponds to the input.

  The range mapping functions return a \link{GRanges} object, with the
  \code{rownames} corresponding to the \link[Biobase]{featureNames} of
  the \link[Biobase]{ExpressionSet} object.
}

\author{Jim Hester, \url{james.f.hester@gmail.com}}

\seealso{
  \itemize{
    \item \link{RangedSummarizedExperiment} objects.

    \item \link[Biobase]{ExpressionSet} objects in the \pkg{Biobase} package.

    \item \link[GenomicFeatures]{TxDb} objects in the \pkg{GenomicFeatures}
          package.
  }
}

\examples{
data(sample.ExpressionSet, package="Biobase")

# using the naive coercion
makeSummarizedExperimentFromExpressionSet(sample.ExpressionSet)

# using probe range mapper
makeSummarizedExperimentFromExpressionSet(sample.ExpressionSet, probeRangeMapper)

# using the gene range mapper
makeSummarizedExperimentFromExpressionSet(sample.ExpressionSet,
                                          geneRangeMapper("TxDb.Hsapiens.UCSC.hg19.knownGene"))
}

\keyword{manip}
