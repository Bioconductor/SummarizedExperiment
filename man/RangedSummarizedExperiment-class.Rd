\name{RangedSummarizedExperiment-class}
\docType{class}

% Class
\alias{class:RangedSummarizedExperiment}
\alias{RangedSummarizedExperiment-class}
\alias{RangedSummarizedExperiment}

\alias{parallel_slot_names,RangedSummarizedExperiment-method}

% Coercion
\alias{coerce,RangedSummarizedExperiment,SummarizedExperiment-method}
\alias{coerce,SummarizedExperiment,RangedSummarizedExperiment-method}

% Accessors
\alias{rowRanges}
\alias{rowRanges,SummarizedExperiment-method}
\alias{rowRanges,RangedSummarizedExperiment-method}
\alias{rowRanges<-}
\alias{rowRanges<-,SummarizedExperiment,NULL-method}
\alias{rowRanges<-,RangedSummarizedExperiment,NULL-method}
\alias{rowRanges<-,SummarizedExperiment,GenomicRanges-method}
\alias{rowRanges<-,SummarizedExperiment,GRangesList-method}
\alias{names,RangedSummarizedExperiment-method}
\alias{names<-,RangedSummarizedExperiment-method}
\alias{dimnames<-,RangedSummarizedExperiment,list-method}

% GenomicRanges compatibility methods
\alias{Compare,ANY,RangedSummarizedExperiment-method}
\alias{Compare,RangedSummarizedExperiment,ANY-method}
\alias{Compare,RangedSummarizedExperiment,RangedSummarizedExperiment-method}
\alias{pcompare,ANY,RangedSummarizedExperiment-method}
\alias{pcompare,RangedSummarizedExperiment,ANY-method}
\alias{pcompare,RangedSummarizedExperiment,RangedSummarizedExperiment-method}
\alias{duplicated,RangedSummarizedExperiment-method}
\alias{elementMetadata,RangedSummarizedExperiment-method}
\alias{elementMetadata<-,RangedSummarizedExperiment-method}
\alias{end,RangedSummarizedExperiment-method}
\alias{end<-,RangedSummarizedExperiment-method}
\alias{granges,RangedSummarizedExperiment-method}
\alias{is.unsorted,RangedSummarizedExperiment-method}
\alias{mcols,RangedSummarizedExperiment-method}
\alias{mcols<-,RangedSummarizedExperiment-method}
\alias{order,RangedSummarizedExperiment-method}
\alias{ranges,RangedSummarizedExperiment-method}
\alias{ranges<-,RangedSummarizedExperiment-method}
\alias{rank,RangedSummarizedExperiment-method}
\alias{seqinfo,RangedSummarizedExperiment-method}
\alias{seqinfo<-,RangedSummarizedExperiment-method}
\alias{seqnames,RangedSummarizedExperiment-method}
\alias{sort,RangedSummarizedExperiment-method}
\alias{split,RangedSummarizedExperiment-method}
\alias{split,RangedSummarizedExperiment,ANY-method}
\alias{start,RangedSummarizedExperiment-method}
\alias{start<-,RangedSummarizedExperiment-method}
\alias{strand,RangedSummarizedExperiment-method}
\alias{strand<-,RangedSummarizedExperiment,ANY-method}
\alias{subset,RangedSummarizedExperiment-method}
\alias{width,RangedSummarizedExperiment-method}
\alias{width<-,RangedSummarizedExperiment-method}

% updateObject
\alias{updateObject,RangedSummarizedExperiment-method}

\title{RangedSummarizedExperiment objects}

\description{

  The RangedSummarizedExperiment class is a matrix-like container where rows
  represent ranges of interest (as a \link{GRanges} or \link{GRangesList}
  object) and columns represent samples (with sample data summarized as a
  \link{DataFrame}). A RangedSummarizedExperiment object contains one or more
  assays, each represented by a matrix-like object of numeric or other mode.

  RangedSummarizedExperiment is a subclass of \link{SummarizedExperiment} and,
  as such, all the methods documented in \code{?SummarizedExperiment}
  also work on a RangedSummarizedExperiment object or any
  \link{SummarizedExperiment} derivative. The methods documented below
  are additional methods that are specific to RangedSummarizedExperiment
  objects.
}

\usage{

## Constructor

# See ?SummarizedExperiment for the constructor function.

## Accessors

rowRanges(x, ...)
rowRanges(x, ...) <- value

## Subsetting

\S4method{subset}{RangedSummarizedExperiment}(x, subset, select, ...)

## rowRanges access
## see 'GRanges compatibility', below
}

\arguments{

  \item{x}{A RangedSummarizedExperiment object or derivative. The
    \code{rowRanges} setter will also accept a \link{SummarizedExperiment}
    \emph{instance} and will first coerce it to RangedSummarizedExperiment
    before it sets \code{value} on it.}

  \item{...}{Further arguments to be passed to or from other methods.}

  \item{value}{A \link[GenomicRanges]{GRanges} or
    \link[GenomicRanges]{GRangesList} object.}

  \item{subset}{An expression which, when evaluated in the
    context of \code{rowRanges(x)}, is a logical vector indicating
    elements or rows to keep: missing values are taken as false.}

  \item{select}{An expression which, when evaluated in the
    context of \code{colData(x)}, is a logical vector indicating
    elements or rows to keep: missing values are taken as false.}

}

\details{

  The rows of a RangedSummarizedExperiment object represent ranges
  (in genomic coordinates) of interest. The ranges of interest are
  described by a \link{GRanges} or a \link{GRangesList} object, accessible
  using the \code{rowRanges} function, described below. The \link{GRanges}
  and \link{GRangesList} classes contains sequence (e.g., chromosome) name,
  genomic coordinates, and strand information. Each range can be
  annotated with additional data; this data might be used to describe
  the range or to summarize results (e.g., statistics of differential
  abundance) relevant to the range. Rows may or may not have row names;
  they often will not.

}

\section{Constructor}{

  RangedSummarizedExperiment instances are constructed using the
  \code{SummarizedExperiment()} function documented in
  \code{?\link{SummarizedExperiment}}.

}

\section{Accessors}{

  In the code snippets below, \code{x} is a RangedSummarizedExperiment object
  or derivative (e.g. a \link[SingleCellExperiment]{SingleCellExperiment}
  object).

  \describe{

    \item{\code{rowRanges(x)}, \code{rowRanges(x) <- value}:}{Get or set the
      row data. \code{value} is a \code{GenomicRanges} object. Row
      names of \code{value} must be NULL or consistent with the existing
      row names of \code{x}.}
   
  }
}

\section{GRanges compatibility (rowRanges access)}{

  Many \link{GRanges} and \link{GRangesList} operations are supported on
  RangedSummarizedExperiment objects, using \code{rowRanges}.

  Supported operations include: \code{\link{pcompare}},
  \code{\link{duplicated}}, \code{\link{end}}, \code{\link{end<-}}, 
  \code{\link{granges}}, \code{\link{is.unsorted}}, \code{\link{match}},
  \code{\link{mcols}}, \code{\link{mcols<-}}, \code{\link{order}},
  \code{\link{ranges}}, \code{\link{ranges<-}}, \code{\link{rank}},
  \code{\link{seqinfo}}, \code{\link{seqinfo<-}}, \code{\link{seqnames}},
  \code{\link{sort}}, \code{\link{start}}, \code{\link{start<-}},
  \code{\link{strand}}, \code{\link{strand<-}},
  \code{\link{width}}, \code{\link{width<-}}.

  See also \code{?\link[SummarizedExperiment]{shift}},
  \code{?\link[SummarizedExperiment]{isDisjoint}},
  \code{?\link[SummarizedExperiment]{coverage}},
  \code{?\link[SummarizedExperiment]{findOverlaps}}, and
  \code{?\link[SummarizedExperiment]{nearest}} for more
  \emph{GRanges compatibility methods}.

  Not all \link{GRanges} operations are supported, because
  they do not make sense for RangedSummarizedExperiment objects
  (e.g., length, name, as.data.frame, c, splitAsList), involve
  non-trivial combination or splitting of rows (e.g., disjoin, gaps,
  reduce, unique), or have not yet been implemented (Ops, map, window,
  window<-).

}

\section{Subsetting}{

  In the code snippets below, \code{x} is a RangedSummarizedExperiment object
  or derivative (e.g. a \link[SingleCellExperiment]{SingleCellExperiment}
  object).

  \describe{

    \item{\code{subset(x, subset, select)}:}{Create a subset of \code{x}
      using an expression \code{subset} referring to columns of
      \code{rowRanges(x)} (including \sQuote{seqnames}, \sQuote{start},
      \sQuote{end}, \sQuote{width}, \sQuote{strand}, and
      \code{names(rowData(x))}) and / or \code{select} referring to
      column names of \code{colData(x)}.}

  }

}

\section{Extension}{

  RangedSummarizedExperiment is implemented as an S4 class, and can be
  extended in the usual way, using \code{contains="RangedSummarizedExperiment"}
  in the new class definition.

  See the \link[SingleCellExperiment]{SingleCellExperiment} class defined
  in the \pkg{SingleCellExperiment} package for an example of such extension.

}

\author{Martin Morgan, \url{mtmorgan@fhcrc.org}}

\seealso{
  \itemize{
    \item \link{SummarizedExperiment} for the parent class of
          RangedSummarizedExperiment and the
          RangedSummarizedExperiment/SummarizedExperiment constructor function.

    \item \link[SummarizedExperiment]{shift},
          \link[SummarizedExperiment]{isDisjoint},
          \link[SummarizedExperiment]{coverage},
          \link[SummarizedExperiment]{findOverlaps}, and
          \link[SummarizedExperiment]{nearest} for more
          \emph{GRanges compatibility methods}.

    \item \link[GenomicRanges]{GRanges} objects in the \pkg{GenomicRanges}
          package.

    \item The \link[SingleCellExperiment]{SingleCellExperiment} defined
          in the \pkg{SingleCellExperiment} package, a subclass of
          RangedSummarizedExperiment specifically designed to represent
          single-cell sequencing data.

  }
}

\examples{
nrows <- 200; ncols <- 6
counts <- matrix(runif(nrows * ncols, 1, 1e4), nrows)
rowRanges <- GRanges(rep(c("chr1", "chr2"), c(50, 150)),
                     IRanges(floor(runif(200, 1e5, 1e6)), width=100),
                     strand=sample(c("+", "-"), 200, TRUE),
                     feature_id=sprintf("ID\%03d", 1:200))
colData <- DataFrame(Treatment=rep(c("ChIP", "Input"), 3),
                     row.names=LETTERS[1:6])
rse <- SummarizedExperiment(assays=SimpleList(counts=counts),
                            rowRanges=rowRanges, colData=colData)
rse
dim(rse)
dimnames(rse)
assayNames(rse)
head(assay(rse))
assays(rse) <- endoapply(assays(rse), asinh)
head(assay(rse))

rowRanges(rse)
rowData(rse)  # same as 'mcols(rowRanges(rse))'
colData(rse)

rse[ , rse$Treatment == "ChIP"]

## cbind() combines objects with the same ranges but different samples:
rse1 <- rse
rse2 <- rse1[ , 1:3]
colnames(rse2) <- letters[1:ncol(rse2)] 
cmb1 <- cbind(rse1, rse2)
dim(cmb1)
dimnames(cmb1)

## rbind() combines objects with the same samples but different ranges:
rse1 <- rse
rse2 <- rse1[1:50, ]
rownames(rse2) <- letters[1:nrow(rse2)] 
cmb2 <- rbind(rse1, rse2)
dim(cmb2)
dimnames(cmb2)

## Coercion to/from SummarizedExperiment:
se0 <- as(rse, "SummarizedExperiment")
se0

as(se0, "RangedSummarizedExperiment")

## Setting rowRanges on a SummarizedExperiment object turns it into a
## RangedSummarizedExperiment object:
se <- se0
rowRanges(se) <- rowRanges
se  # RangedSummarizedExperiment

## Sanity checks:
stopifnot(identical(assays(se0), assays(rse)))
stopifnot(identical(dim(se0), dim(rse)))
stopifnot(identical(dimnames(se0), dimnames(rse)))
stopifnot(identical(rowData(se0), rowData(rse)))
stopifnot(identical(colData(se0), colData(rse)))
}
