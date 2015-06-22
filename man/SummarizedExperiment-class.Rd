\name{RangedSummarizedExperiment-class}
\docType{class}

% Classes
\alias{SummarizedExperiment0-class}
\alias{SummarizedExperiment0}
\alias{RangedSummarizedExperiment-class}
\alias{RangedSummarizedExperiment}

% Constructor
\alias{SummarizedExperiment}
\alias{SummarizedExperiment,missing-method}
\alias{SummarizedExperiment,matrix-method}
\alias{SummarizedExperiment,list-method}
\alias{SummarizedExperiment,SimpleList-method}

% Coercion
\alias{coerce,RangedSummarizedExperiment,SummarizedExperiment0-method}

% Accessors
\alias{length,SummarizedExperiment0-method}
\alias{names,SummarizedExperiment0-method}
\alias{names<-,SummarizedExperiment0-method}
\alias{metadata<-,SummarizedExperiment0-method}
\alias{colData}
\alias{colData,SummarizedExperiment0-method}
\alias{colData<-}
\alias{colData<-,SummarizedExperiment0,DataFrame-method}
\alias{assayNames}
\alias{assayNames,SummarizedExperiment0-method}
\alias{assayNames<-}
\alias{assayNames<-,SummarizedExperiment0,character-method}
\alias{assays}
\alias{assays,SummarizedExperiment0-method}
\alias{assays<-}
\alias{assays<-,SummarizedExperiment0,SimpleList-method}
\alias{assays<-,SummarizedExperiment0,list-method}
\alias{assay}
\alias{assay,SummarizedExperiment0,missing-method}
\alias{assay,SummarizedExperiment0,numeric-method}
\alias{assay,SummarizedExperiment0,character-method}
\alias{assay<-}
\alias{assay<-,SummarizedExperiment0,missing,matrix-method}
\alias{assay<-,SummarizedExperiment0,numeric,matrix-method}
\alias{assay<-,SummarizedExperiment0,character,matrix-method}
\alias{dim,SummarizedExperiment0-method}
\alias{dimnames,SummarizedExperiment0-method}

\alias{rowRanges}
\alias{rowRanges,RangedSummarizedExperiment-method}
\alias{rowRanges<-}
\alias{rowRanges<-,RangedSummarizedExperiment,GenomicRanges-method}
\alias{rowRanges<-,RangedSummarizedExperiment,GRangesList-method}
\alias{names,RangedSummarizedExperiment-method}
\alias{names<-,RangedSummarizedExperiment-method}
\alias{dimnames,RangedSummarizedExperiment-method}
\alias{dimnames<-,RangedSummarizedExperiment,list-method}
\alias{dimnames<-,RangedSummarizedExperiment,NULL-method}

% Subsetting
\alias{[,RangedSummarizedExperiment-method}
\alias{[,RangedSummarizedExperiment,ANY-method}
\alias{[<-,RangedSummarizedExperiment,ANY,ANY,RangedSummarizedExperiment-method}

% Quick colData access
\alias{[[,SummarizedExperiment0,ANY,missing-method}
\alias{[[<-,SummarizedExperiment0,ANY,missing-method}
\alias{$,SummarizedExperiment0-method}
\alias{$<-,SummarizedExperiment0-method}

% Display
\alias{show,SummarizedExperiment0-method}

% Combine 
\alias{rbind,RangedSummarizedExperiment-method}
\alias{cbind,RangedSummarizedExperiment-method}

% Backward compatibility with "classic" SummarizedExperiment objects
\alias{exptData,SummarizedExperiment0-method}
\alias{exptData<-,SummarizedExperiment0-method}
\alias{updateObject,RangedSummarizedExperiment-method}

% GenomicRanges compatibility methods
\alias{Compare,ANY,RangedSummarizedExperiment-method}
\alias{Compare,RangedSummarizedExperiment,ANY-method}
\alias{Compare,RangedSummarizedExperiment,RangedSummarizedExperiment-method}
\alias{compare,ANY,RangedSummarizedExperiment-method}
\alias{compare,RangedSummarizedExperiment,ANY-method}
\alias{compare,RangedSummarizedExperiment,RangedSummarizedExperiment-method}
\alias{duplicated,RangedSummarizedExperiment-method}
\alias{elementMetadata,RangedSummarizedExperiment-method}
\alias{elementMetadata<-,RangedSummarizedExperiment-method}
\alias{end,RangedSummarizedExperiment-method}
\alias{end<-,RangedSummarizedExperiment-method}
\alias{granges,RangedSummarizedExperiment-method}
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
\alias{values,RangedSummarizedExperiment-method}
\alias{values<-,RangedSummarizedExperiment-method}
\alias{width,RangedSummarizedExperiment-method}
\alias{width<-,RangedSummarizedExperiment-method}

\title{RangedSummarizedExperiment objects}

\description{

  The RangedSummarizedExperiment class is a matrix-like container where rows
  represent ranges of interest (as a \code{\linkS4class{GRanges} or
  \linkS4class{GRangesList}-class}) and columns represent samples (with
  sample data summarized as a \code{\linkS4class{DataFrame}-class}). A
  \code{RangedSummarizedExperiment} contains one or more assays, each
  represented by a matrix-like object of numeric or other mode.

}

\usage{

## Constructors

SummarizedExperiment(assays, ...)
\S4method{SummarizedExperiment}{SimpleList}(assays, rowRanges=GRangesList(),
    colData=DataFrame(), metadata=list(), exptData=SimpleList())
\S4method{SummarizedExperiment}{missing}(assays, ...)
\S4method{SummarizedExperiment}{list}(assays, ...)
\S4method{SummarizedExperiment}{matrix}(assays, ...)

## Accessors

assayNames(x, ...)
assayNames(x, ...) <- value
assays(x, ..., withDimnames=TRUE)
assays(x, ..., withDimnames=TRUE) <- value
assay(x, i, ...)
assay(x, i, ...) <- value
rowRanges(x, ...)
rowRanges(x, ...) <- value
colData(x, ...)
colData(x, ...) <- value
\S4method{dim}{SummarizedExperiment0}(x)
\S4method{dimnames}{RangedSummarizedExperiment}(x)
\S4method{dimnames}{RangedSummarizedExperiment,NULL}(x) <- value
\S4method{dimnames}{RangedSummarizedExperiment,list}(x) <- value

## Quick colData access

\S4method{$}{SummarizedExperiment0}(x, name)
\S4method{$}{SummarizedExperiment0}(x, name) <- value
\S4method{[[}{SummarizedExperiment0,ANY,missing}(x, i, j, ...)
\S4method{[[}{SummarizedExperiment0,ANY,missing}(x, i, j, ...) <- value

## rowRanges access
## see 'GRanges compatibility', below

## Subsetting

\S4method{[}{RangedSummarizedExperiment}(x, i, j, ..., drop=TRUE)
\S4method{[}{RangedSummarizedExperiment,ANY,ANY,RangedSummarizedExperiment}(x, i, j) <- value
\S4method{subset}{RangedSummarizedExperiment}(x, subset, select, ...)

## Combining 

\S4method{cbind}{RangedSummarizedExperiment}(..., deparse.level=1)
\S4method{rbind}{RangedSummarizedExperiment}(..., deparse.level=1)

}

\arguments{

  \item{assays}{A \code{list} or \code{SimpleList} of matrix elements,
    or a \code{matrix}. All elements of the list must have the same
    dimensions, and dimension names (if present) must be consistent
    across elements and with the row names of \code{rowRanges} and
    \code{colData}.}

  \item{rowRanges}{A \code{GRanges} or \code{GRangesList} instance describing 
    the ranges of interest. Row names, if present, become the row names of the
    \code{RangedSummarizedExperiment}. The length of the \code{GRanges} or the
    \code{GRangesList} must equal the number of rows of the 
    matrices in \code{assays}.}

  \item{colData}{An optional \code{DataFrame} describing the
    samples. Row names, if present, become the column names of the
    \code{RangedSummarizedExperiment}.}

  \item{metadata}{An optional \code{list} of arbitrary content
    describing the overall experiment.}

  \item{exptData}{Deprecated. Please use the \code{metadata} argument instead
    (see above).}

  \item{...}{For \code{SummarizedExperiment}, S4 methods \code{list}
    and \code{matrix}, arguments identical to those of the
    \code{SimpleList} method.

    For \code{assay}, \code{...} may contain \code{withDimnames}, which is
    forwarded to \code{assays}.

    For \code{cbind, rbind}, \code{...} contains
    \code{RangedSummarizedExperiment} objects to be combined.

    For other accessors, ignored.

  }

  \item{x}{An instance of \code{RangedSummarizedExperiment}-class.}

  \item{i, j}{For \code{assay}, \code{assay<-}, \code{i} is a integer or
    numeric scalar; see \sQuote{Details} for additional constraints.

    For \code{[,RangedSummarizedExperiment},
    \code{[,RangedSummarizedExperiment<-}, \code{i}, \code{j} are instances
    that can act to subset the underlying \code{rowRanges}, \code{colData},
    and \code{matrix} elements of \code{assays}.

    For \code{[[,RangedSummarizedExperiment},
    \code{[[<-,RangedSummarizedExperiment}, \code{i} is a scalar index (e.g.,
    \code{character(1)} or \code{integer(1)}) into a column of
    \code{colData}.

  }

  \item{subset}{An expression which, when evaluated in the
    context of \code{rowRanges(x)}, is a logical vector indicating
    elements or rows to keep: missing values are taken as false.}

  \item{select}{An expression which, when evaluated in the
    context of \code{colData(x)}, is a logical vector indicating
    elements or rows to keep: missing values are taken as false.}

  \item{name}{A symbol representing the name of a column of
    \code{colData}.}

  \item{withDimnames}{A \code{logical(1)}, indicating whether dimnames
    should be applied to extracted assay elements. Setting
    \code{withDimnames=FALSE} increases the speed and memory efficiency
    with which assays are extracted. \code{withDimnames=TRUE} in the
    getter \code{assays<-} allows efficient complex assignments (e.g.,
    updating names of assays, \code{names(assays(x, withDimnames=FALSE))
    = ...} is more efficient than \code{names(assays(x)) = ...}); it
    does not influence actual assignment of dimnames to assays.}

  \item{drop}{A \code{logical(1)}, ignored by these methods.}

  \item{value}{An instance of a class specified in the S4 method
    signature or as outlined in \sQuote{Details}.}

  \item{deparse.level}{See \code{?base::\link[base]{cbind}} for a description of
    this argument.}

}

\details{

  The \code{RangedSummarizedExperiment} class is meant for numeric and other
  data types derived from a sequencing experiment. The structure is
  rectangular like a \code{matrix}, but with additional annotations on
  the rows and columns, and with the possibility to manage several
  assays simultaneously.

  The rows of a \code{RangedSummarizedExperiment} instance represent ranges
  (in genomic coordinates) of interest. The ranges of interest are
  described by a \code{\linkS4class{GRanges}-class} or a
  \code{\linkS4class{GRangesList}-class} instance, accessible using the
  \code{rowRanges} function, described below. The \code{GRanges} and
  \code{GRangesList} classes contains sequence (e.g., chromosome) name,
  genomic coordinates, and strand information. Each range can be
  annotated with additional data; this data might be used to describe
  the range or to summarize results (e.g., statistics of differential
  abundance) relevant to the range. Rows may or may not have row names;
  they often will not.

  Each column of a \code{RangedSummarizedExperiment} instance represents a
  sample.  Information about the samples are stored in a
  \code{\linkS4class{DataFrame}-class}, accessible using the function
  \code{colData}, described below. The \code{DataFrame} must have as
  many rows as there are columns in the \code{RangedSummarizedExperiment},
  with each row of the \code{DataFrame} providing information on the
  sample in the corresponding column of the
  \code{RangedSummarizedExperiment}. Columns of the \code{DataFrame} represent
  different sample attributes, e.g., tissue of origin, etc. Columns of
  the \code{DataFrame} can themselves be annotated (via the
  \code{\link{mcols}} function). Column names typically provide a short
  identifier unique to each sample.

  A \code{RangedSummarizedExperiment} can also contain information about the
  overall experiment, for instance the lab in which it was conducted,
  the publications with which it is associated, etc. This information is
  stored as a \code{list} object, accessible using the \code{metadata}
  function. The form of the data associated with the experiment is left to
  the discretion of the user.

  The \code{RangedSummarizedExperiment} is appropriate for matrix-like
  data. The data are accessed using the \code{assays} function,
  described below. This returns a \code{SimpleList}-class instance. Each
  element of the list must itself be a matrix (of any mode) and must
  have dimensions that are the same as the dimensions of the
  \code{RangedSummarizedExperiment} in which they are stored. Row and column
  names of each matrix must either be \code{NULL} or match those of the
  \code{RangedSummarizedExperiment} during construction. It is convenient for
  the elements of \code{SimpleList} of assays to be named.

  The \code{RangedSummarizedExperiment} class has the following slots; this
  detail of class structure is not relevant to the user.
  \describe{
    \item{\code{metadata}}{A \code{list} containing information about the
      overall experiment.}
    \item{\code{rowRanges}}{A \link{GRanges}-class instance defining the
      ranges of interest and associated metadata.}
    \item{\code{colData}}{A \link{DataFrame}-class instance describing
      the samples and associated metadata.}
    \item{\code{assays}}{A \link{SimpleList}-class instance, each
      element of which is a matrix summarizing data associated with the
      corresponding range and sample.}
  }
}

\section{Constructor}{

  Instances are constructed using the \code{SummarizedExperiment}
  function with arguments outlined above.

}

\section{Accessors}{

  In the following code snippets, \code{x} is a
  \code{RangedSummarizedExperiment} instance.

  \describe{

    \item{\code{assays(x)}, \code{assays(x) <- value}:}{Get or set the
      assays. \code{value} is a \code{list} or \code{SimpleList}, each
      element of which is a matrix with the same dimensions as
      \code{x}.}

    \item{\code{assay(x, i)}, \code{assay(x, i) <- value}:}{A convenient
      alternative (to \code{assays(x)[[i]]}, \code{assays(x)[[i]] <-
      value}) to get or set the \code{i}th (default first) assay
      element. \code{value} must be a matrix of the same dimension as
      \code{x}, and with dimension names \code{NULL} or consistent with
      those of \code{x}.}

    \item{\code{assayNames(x)}, \code{assayNames(x) <- value}:}{Get or
     set the names of \code{assay()} elements.}

    \item{\code{rowRanges(x)}, \code{rowRanges(x) <- value}:}{Get or set the
      row data. \code{value} is a \code{GenomicRanges} instance. Row
      names of \code{value} must be NULL or consistent with the existing
      row names of \code{x}.}
   
    \item{\code{colData(x)}, \code{colData(x) <- value}:}{Get or set the
      column data. \code{value} is a \code{DataFrame} instance. Row
      names of \code{value} must be NULL or consistent with the existing
      column names of \code{x}.}

    \item{\code{metadata(x)}, \code{metadata(x) <- value}:}{Get or set
      the experiment data. \code{value} is a \code{list} with arbitrary
      content.}

    \item{\code{dim(x)}:}{Get the dimensions (ranges x samples) of the
      \code{RangedSummarizedExperiment}.}

    \item{\code{dimnames(x)}, \code{dimnames(x) <- value}:}{Get or set
      the dimension names. \code{value} is usually a list of length 2,
      containing elements that are either \code{NULL} or vectors of
      appropriate length for the corresponding dimension. \code{value}
      can be \code{NULL}, which removes dimension names. This method
      implies that \code{rownames}, \code{rownames<-}, \code{colnames},
      and \code{colnames<-} are all available.}

  }
}

\section{GRanges compatibility (rowRanges access)}{

  Many \code{\linkS4class{GRanges}-class} and
  \code{\linkS4class{GRangesList}-class} operations are supported on
  \sQuote{RangedSummarizedExperiment} and derived instances, using
  \code{rowRanges}.

  Supported operations include: \code{\link{compare}},
  \code{\link{duplicated}}, \code{\link{end}}, \code{\link{end<-}}, 
  \code{\link{granges}}, \code{\link{match}}, \code{\link{mcols}},
  \code{\link{mcols<-}}, \code{\link{order}}, \code{\link{ranges}},
  \code{\link{ranges<-}}, \code{\link{rank}}, \code{\link{seqinfo}},
  \code{\link{seqinfo<-}}, \code{\link{seqnames}}, \code{\link{sort}}, 
  \code{\link{start}}, \code{\link{start<-}},
  \code{\link{strand}}, \code{\link{strand<-}},
  \code{\link{width}}, \code{\link{width<-}}.

  See also \code{?\link[SummarizedExperiment]{shift}},
  \code{?\link[SummarizedExperiment]{isDisjoint}},
  \code{?\link[SummarizedExperiment]{coverage}},
  \code{?\link[SummarizedExperiment]{findOverlaps}}, and
  \code{?\link[SummarizedExperiment]{nearest}} for more
  \emph{GRanges compatibility methods}.

  Not all \code{\link{GRanges}-class} operations are supported, because
  they do not make sense for \sQuote{RangedSummarizedExperiment} objects
  (e.g., length, name, as.data.frame, c, splitAsList), involve
  non-trivial combination or splitting of rows (e.g., disjoin, gaps,
  reduce, unique), or have not yet been implemented (Ops, map, window,
  window<-).

}

\section{Subsetting}{

  In the code snippets below, \code{x} is a \code{RangedSummarizedExperiment}
  instance.

  \describe{

    \item{\code{x[i,j]}, \code{x[i,j] <- value}:}{Create or replace a
      subset of \code{x}. \code{i}, \code{j} can be \code{numeric},
      \code{logical}, \code{character}, or \code{missing}. \code{value}
      must be a \code{RangedSummarizedExperiment} instance with dimensions,
      dimension names, and assay elements consistent with the subset
      \code{x[i,j]} being replaced.}

    \item{\code{subset(x, subset, select)}:}{Create a subset of \code{x}
      using an expression \code{subset} referring to columns of
      \code{rowRanges(x)} (including \sQuote{seqnames}, \sQuote{start},
      \sQuote{end}, \sQuote{width}, \sQuote{strand}, and
      \code{names(mcols(x))}) and / or \code{select} referring to
      column names of \code{colData(x)}.}

  }

  Additional subsetting accessors provide convenient access to
  \code{colData} columns

  \describe{

    \item{\code{x$name}, \code{x$name <- value}}{Access or replace
    column \code{name} in \code{x}.}

    \item{\code{x[[i, ...]]}, \code{x[[i, ...]] <- value}}{Access or
    replace column \code{i} in \code{x}.}

  }

}

\section{Combining}{

  In the code snippets below, \code{...} are \code{RangedSummarizedExperiment}
  instances to be combined.

  \describe{

    \item{\code{cbind(...)}, \code{rbind(...)}:}{
      \code{cbind} combines objects with identical ranges (\code{rowRanges}) 
      but different samples (columns in \code{assays}). The colnames in
      \code{colData} must match or an error is thrown. Duplicate columns
      of \code{mcols(rowRanges(RangedSummarizedExperiment))} must contain the
      same data. Data in \code{assays} are combined by name matching; if all
      names are NULL matching is by position. A mixture of names and NULL
      throws an error.

      \code{rbind} combines objects with different ranges (\code{rowRanges}) 
      and the same subjects (columns in \code{assays}). Duplicate columns
      of \code{colData} must contain the same data. 

      \code{metadata} from all objects are combined into a
      \code{list} with no name checking.

    }

  }

}

\section{Implementation and Extension}{

  This section contains advanced material meant for package developers.

  \code{RangedSummarizedExperiment} is implemented as an S4 class, and can be
  extended in the usual way, using
  \code{contains="RangedSummarizedExperiment"} in the new class definition.

  In addition, the representation of the \code{assays} slot of
  \code{RangedSummarizedExperiment} is as a virtual class \code{Assays}. This
  allows derived classes (\code{contains="Assays"}) to easily implement
  alternative requirements for the assays, e.g., backed by file-based
  storage like NetCDF or the \code{ff} package, while re-using the
  existing \code{RangedSummarizedExperiment} class without modification. The
  requirements on \code{Assays} are list-like semantics (e.g.,
  \code{sapply}, \code{[[} subsetting, \code{names}) with elements
  having matrix- or array-like semantics (e.g., \code{dim},
  \code{dimnames}). These requirements can be made more precise if
  developers express interest.

  The current \code{assays} slot is implemented as a reference class
  that has copy-on-change semantics. This means that modifying non-assay
  slots does not copy the (large) assay data, and at the same time the
  user is not surprised by reference-based semantics. Updates to
  non-assay slots are very fast; updating the assays slot itself can be
  5x or more faster than with an S4 instance in the slot. One useful
  technique when working with \code{assay} or \code{assays} function is
  use of the \code{withDimnames=FALSE} argument, which benefits speed
  and memory use by not copying dimnames from the row- and colData
  elements to each assay.

  In a little more detail, a small reference class hierarchy (not
  exported from the GenomicRanges name space) defines a reference class
  \code{ShallowData} with a single field \code{data} of type \code{ANY},
  and a derived class \code{ShallowSimpleListAssays} that specializes
  the type of \code{data} as \code{SimpleList}, and
  \code{contains=c("ShallowData", "Assays")}. The assays slot contains
  an instance of \code{ShallowSimpleListAssays}. Invoking
  \code{assays()} on a \code{RangedSummarizedExperiment} re-dispatches from
  the \code{assays} slot to retrieve the \code{SimpleList} from the
  field of the reference class. This was achieved by implementing a
  generic (not exported) \code{value(x, name, ...)}, with a method
  implemented on \code{RangedSummarizedExperiment} that retrieves a slot when
  \code{name} is a slot containing an S4 object in \code{x}, and a field
  when \code{name} is a slot containing a \code{ShallowData} instance in
  \code{x}. Copy-on-change semantics is maintained by implementing the
  \code{clone} method (\code{clone} methods are supposed to do a deep
  copy, \code{update} methods a shallow copy; the \code{clone} generic
  is introduced, and not exported, in the GenomicRanges package). The
  \sQuote{getter} and \sQuote{setter} code for methods implemented on
  \code{RangedSummarizedExperiment} use \code{value} for slot access, and
  \code{clone} for replacement. This makes it easy to implement
  \code{ShallowData} instances for other slots if the need arises.

}

\author{Martin Morgan, \url{mtmorgan@fhcrc.org}}

\seealso{
  \itemize{
    \item \link[SummarizedExperiment]{shift},
          \link[SummarizedExperiment]{isDisjoint},
          \link[SummarizedExperiment]{coverage},
          \link[SummarizedExperiment]{findOverlaps}, and
          \link[SummarizedExperiment]{nearest} for more
          \emph{GRanges compatibility methods}.

    \item \link[GenomicRanges]{GRanges} objects in the \pkg{GenomicRanges}
          package.

    \item \link[S4Vectors]{DataFrame}, \link[S4Vectors]{SimpleList}, and
          \link[S4Vectors]{Annotated} objects in the \pkg{S4Vectors} package.
  }
}

\examples{
  nrows <- 200; ncols <- 6
  counts <- matrix(runif(nrows * ncols, 1, 1e4), nrows)
  rowRanges <- GRanges(rep(c("chr1", "chr2"), c(50, 150)),
                     IRanges(floor(runif(200, 1e5, 1e6)), width=100),
                     strand=sample(c("+", "-"), 200, TRUE))
  colData <- DataFrame(Treatment=rep(c("ChIP", "Input"), 3),
                       row.names=LETTERS[1:6])
  sset <- SummarizedExperiment(assays=SimpleList(counts=counts),
                 rowRanges=rowRanges, colData=colData)
  sset
  assayNames(sset)
  assays(sset) <- endoapply(assays(sset), asinh)
  head(assay(sset))

  sset[, sset$Treatment == "ChIP"]

  ## cbind combines objects with the same ranges and different samples.
  se1 <- sset
  se2 <- se1[,1:3]
  colnames(se2) <- letters[1:ncol(se2)] 
  cmb1 <- cbind(se1, se2)

  ## rbind combines objects with the same samples and different ranges.
  se1 <- sset
  se2 <- se1[1:50,]
  rownames(se2) <- letters[1:nrow(se2)] 
  cmb2 <- rbind(se1, se2)

}
