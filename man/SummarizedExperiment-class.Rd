\name{SummarizedExperiment-class}
\docType{class}

% Class
\alias{class:SummarizedExperiment}
\alias{SummarizedExperiment-class}

\alias{class:SummarizedExperiment0}
\alias{SummarizedExperiment0-class}
\alias{SummarizedExperiment0}

% Accessors
\alias{length,SummarizedExperiment-method}
\alias{names,SummarizedExperiment-method}
\alias{names<-,SummarizedExperiment-method}
\alias{exptData}
\alias{exptData,SummarizedExperiment-method}
\alias{exptData<-}
\alias{exptData<-,SummarizedExperiment-method}
\alias{rowData}
\alias{rowData,SummarizedExperiment-method}
\alias{rowData<-}
\alias{rowData<-,SummarizedExperiment-method}
\alias{colData}
\alias{colData,SummarizedExperiment-method}
\alias{colData<-}
\alias{colData<-,SummarizedExperiment,DataFrame-method}
\alias{assays}
\alias{assays,SummarizedExperiment-method}
\alias{assays<-}
\alias{assays<-,SummarizedExperiment,SimpleList-method}
\alias{assays<-,SummarizedExperiment,list-method}
\alias{assay}
\alias{assay,SummarizedExperiment,missing-method}
\alias{assay,SummarizedExperiment,numeric-method}
\alias{assay,SummarizedExperiment,character-method}
\alias{assay<-}
\alias{assay<-,SummarizedExperiment,missing,matrix-method}
\alias{assay<-,SummarizedExperiment,numeric,matrix-method}
\alias{assay<-,SummarizedExperiment,character,matrix-method}
\alias{assayNames}
\alias{assayNames,SummarizedExperiment-method}
\alias{assayNames<-}
\alias{assayNames<-,SummarizedExperiment,character-method}
\alias{dim,SummarizedExperiment-method}
\alias{dimnames,SummarizedExperiment-method}
\alias{dimnames<-,SummarizedExperiment,list-method}
\alias{dimnames<-,SummarizedExperiment,NULL-method}

% Subsetting
\alias{[,SummarizedExperiment-method}
\alias{[,SummarizedExperiment,ANY-method}
\alias{[<-,SummarizedExperiment,ANY,ANY,SummarizedExperiment-method}
\alias{extractROWS,SummarizedExperiment,ANY-method}
\alias{replaceROWS,SummarizedExperiment-method}

% Quick colData access
\alias{[[,SummarizedExperiment,ANY,missing-method}
\alias{[[<-,SummarizedExperiment,ANY,missing-method}
\alias{$,SummarizedExperiment-method}
\alias{$<-,SummarizedExperiment-method}

% Display
\alias{show,SummarizedExperiment-method}

% Combine 
\alias{rbind,SummarizedExperiment-method}
\alias{cbind,SummarizedExperiment-method}

\title{SummarizedExperiment objects}

\description{

  The SummarizedExperiment class is a matrix-like container where rows
  represent features of interest (e.g. genes, transcripts, exons, etc...)
  and columns represent samples (with sample data summarized as a
  \link{DataFrame}). A SummarizedExperiment object contains one or more
  assays, each represented by a matrix-like object of numeric or other mode.

  Note that SummarizedExperiment is the parent of the
  \link{RangedSummarizedExperiment} class which means that all the methods
  documented below also work on a \link{RangedSummarizedExperiment} object.
}

\usage{

## Constructor

# See ?RangedSummarizedExperiment for the constructor function.

## Accessors

assayNames(x, ...)
assayNames(x, ...) <- value
assays(x, ..., withDimnames=TRUE)
assays(x, ..., withDimnames=TRUE) <- value
assay(x, i, ...)
assay(x, i, ...) <- value
rowData(x, ...)
rowData(x, ...) <- value
colData(x, ...)
colData(x, ...) <- value
#dim(x)
#dimnames(x)
#dimnames(x) <- value

## Quick colData access

\S4method{$}{SummarizedExperiment}(x, name)
\S4method{$}{SummarizedExperiment}(x, name) <- value
\S4method{[[}{SummarizedExperiment,ANY,missing}(x, i, j, ...)
\S4method{[[}{SummarizedExperiment,ANY,missing}(x, i, j, ...) <- value

## Subsetting

\S4method{[}{SummarizedExperiment}(x, i, j, ..., drop=TRUE)
\S4method{[}{SummarizedExperiment,ANY,ANY,SummarizedExperiment}(x, i, j) <- value

## Combining 

\S4method{cbind}{SummarizedExperiment}(..., deparse.level=1)
\S4method{rbind}{SummarizedExperiment}(..., deparse.level=1)

}

\arguments{

  \item{x}{A SummarizedExperiment object.}

  \item{...}{
    For \code{assay}, \code{...} may contain \code{withDimnames}, which is
    forwarded to \code{assays}.

    For \code{rowData}, arguments passed thru \code{...} are forwarded to
    \code{\link[S4Vectors]{mcols}}.

    For \code{cbind}, \code{rbind}, \code{...} contains SummarizedExperiment
    objects to be combined.

    For other accessors, ignored.
  }

  \item{i, j}{
    For \code{assay}, \code{assay<-}, \code{i} is an integer or
    numeric scalar; see \sQuote{Details} for additional constraints.

    For \code{[,SummarizedExperiment},
    \code{[,SummarizedExperiment<-}, \code{i}, \code{j} are subscripts
    that can act to subset the rows and columns of \code{x}, that is the
    \code{matrix} elements of \code{assays}.

    For \code{[[,SummarizedExperiment},
    \code{[[<-,SummarizedExperiment}, \code{i} is a scalar index (e.g.,
    \code{character(1)} or \code{integer(1)}) into a column of
    \code{colData}.
  }

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

  \item{value}{An object of a class specified in the S4 method
    signature or as outlined in \sQuote{Details}.}

  \item{deparse.level}{See \code{?base::\link[base]{cbind}} for a description of
    this argument.}

}

\details{

  The SummarizedExperiment class is meant for numeric and other
  data types derived from a sequencing experiment. The structure is
  rectangular like a \code{matrix}, but with additional annotations on
  the rows and columns, and with the possibility to manage several
  assays simultaneously.

  The rows of a SummarizedExperiment object represent features
  of interest. Information about these features is stored in a
  \link{DataFrame} object, accessible using the function
  \code{rowData}. The \link{DataFrame} must have as many rows
  as there are rows in the SummarizedExperiment object, with each row
  of the \link{DataFrame} providing information on the feature in the
  corresponding row of the SummarizedExperiment object. Columns of the
  \link{DataFrame} represent different attributes of the features
  of interest, e.g., gene or transcript IDs, etc.

  Each column of a SummarizedExperiment object represents a sample.
  Information about the samples are stored in a \link{DataFrame},
  accessible using the function \code{colData}, described below.
  The \link{DataFrame} must have as many rows as there are
  columns in the SummarizedExperiment object, with each row of the
  \link{DataFrame} providing information on the sample in the
  corresponding column of the SummarizedExperiment object.
  Columns of the \link{DataFrame} represent different sample
  attributes, e.g., tissue of origin, etc. Columns of the
  \link{DataFrame} can themselves be annotated (via the
  \code{\link[S4Vectors]{mcols}} function). Column names typically
  provide a short identifier unique to each sample.

  A SummarizedExperiment object can also contain information about
  the overall experiment, for instance the lab in which it was conducted,
  the publications with which it is associated, etc. This information is
  stored as a \code{list} object, accessible using the \code{metadata}
  function. The form of the data associated with the experiment is left to
  the discretion of the user.

  The SummarizedExperiment container is appropriate for matrix-like
  data. The data are accessed using the \code{assays} function,
  described below. This returns a \link{SimpleList} object. Each
  element of the list must itself be a matrix (of any mode) and must
  have dimensions that are the same as the dimensions of the
  SummarizedExperiment in which they are stored. Row and column
  names of each matrix must either be \code{NULL} or match those of the
  SummarizedExperiment during construction. It is convenient for
  the elements of \link{SimpleList} of assays to be named.

}

\section{Constructor}{

  SummarizedExperiment instances are constructed using the
  \code{SummarizedExperiment} function documented in
  \code{?\link{RangedSummarizedExperiment}}.

}

\section{Accessors}{

  In the following code snippets, \code{x} is a SummarizedExperiment
  object.

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

    \item{\code{rowData(x)}, \code{rowData(x) <- value}:}{Get or set the
      row data. \code{value} is a \link{DataFrame} object. Row
      names of \code{value} must be NULL or consistent with the existing
      row names of \code{x}.}

    \item{\code{colData(x)}, \code{colData(x) <- value}:}{Get or set the
      column data. \code{value} is a \link{DataFrame} object. Row
      names of \code{value} must be NULL or consistent with the existing
      column names of \code{x}.}

    \item{\code{metadata(x)}, \code{metadata(x) <- value}:}{Get or set
      the experiment data. \code{value} is a \code{list} with arbitrary
      content.}

    \item{\code{dim(x)}:}{Get the dimensions (features of interest x samples)
      of the SummarizedExperiment.}

    \item{\code{dimnames(x)}, \code{dimnames(x) <- value}:}{Get or set
      the dimension names. \code{value} is usually a list of length 2,
      containing elements that are either \code{NULL} or vectors of
      appropriate length for the corresponding dimension. \code{value}
      can be \code{NULL}, which removes dimension names. This method
      implies that \code{rownames}, \code{rownames<-}, \code{colnames},
      and \code{colnames<-} are all available.}

  }
}

\section{Subsetting}{

  In the code snippets below, \code{x} is a SummarizedExperiment object.

  \describe{

    \item{\code{x[i,j]}, \code{x[i,j] <- value}:}{Create or replace a
      subset of \code{x}. \code{i}, \code{j} can be \code{numeric},
      \code{logical}, \code{character}, or \code{missing}. \code{value}
      must be a SummarizedExperiment object with dimensions,
      dimension names, and assay elements consistent with the subset
      \code{x[i,j]} being replaced.}

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

  In the code snippets below, \code{...} are SummarizedExperiment objects
  to be combined.

  \describe{

    \item{\code{cbind(...)}:}{
      \code{cbind} combines objects with the same features of interest
      but different samples (columns in \code{assays}).
      The colnames in \code{colData(SummarizedExperiment)} must match or
      an error is thrown.
      Duplicate columns of \code{rowData(SummarizedExperiment)} must
      contain the same data.

      Data in \code{assays} are combined by name matching; if all assay
      names are NULL matching is by position. A mixture of names and NULL
      throws an error.

      \code{metadata} from all objects are combined into a \code{list}
      with no name checking.
    }

    \item{\code{rbind(...)}:}{
      \code{rbind} combines objects with the same samples
      but different features of interest (rows in \code{assays}).
      The colnames in \code{rowData(SummarizedExperiment)} must match or
      an error is thrown.
      Duplicate columns of \code{colData(SummarizedExperiment)} must
      contain the same data. 

      Data in \code{assays} are combined by name matching; if all assay
      names are NULL matching is by position. A mixture of names and NULL
      throws an error.

      \code{metadata} from all objects are combined into a \code{list}
      with no name checking.
    }

  }

}

\section{Implementation and Extension}{

  This section contains advanced material meant for package developers.

  SummarizedExperiment is implemented as an S4 class, and can be extended in
  the usual way, using \code{contains="SummarizedExperiment"} in the new
  class definition.

  In addition, the representation of the \code{assays} slot of
  SummarizedExperiment is as a virtual class Assays. This
  allows derived classes (\code{contains="Assays"}) to easily implement
  alternative requirements for the assays, e.g., backed by file-based
  storage like NetCDF or the \code{ff} package, while re-using the existing
  SummarizedExperiment class without modification.
  See \link{Assays} for more information.

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
}

\author{Martin Morgan, \url{mtmorgan@fhcrc.org}}

\seealso{
  \itemize{
    \item \link{RangedSummarizedExperiment} objects.

    \item \link[S4Vectors]{DataFrame}, \link[S4Vectors]{SimpleList}, and
          \link[S4Vectors]{Annotated} objects in the \pkg{S4Vectors} package.

    \item The \code{\link[S4Vectors]{metadata}} and
          \code{\link[S4Vectors]{mcols}} accessors in the \pkg{S4Vectors}
          package.
  }
}

\examples{
nrows <- 200; ncols <- 6
counts <- matrix(runif(nrows * ncols, 1, 1e4), nrows)
colData <- DataFrame(Treatment=rep(c("ChIP", "Input"), 3),
                     row.names=LETTERS[1:6])
se0 <- SummarizedExperiment(assays=SimpleList(counts=counts),
                            colData=colData)
se0
dim(se0)
dimnames(se0)
assayNames(se0)
head(assay(se0))
assays(se0) <- endoapply(assays(se0), asinh)
head(assay(se0))

rowData(se0)
colData(se0)

se0[, se0$Treatment == "ChIP"]

## cbind() combines objects with the same features of interest
## but different samples:
se1 <- se0
se2 <- se1[,1:3]
colnames(se2) <- letters[1:ncol(se2)] 
cmb1 <- cbind(se1, se2)
dim(cmb1)
dimnames(cmb1)

## rbind() combines objects with the same samples but different
## features of interest:
se1 <- se0
se2 <- se1[1:50,]
rownames(se2) <- letters[1:nrow(se2)] 
cmb2 <- rbind(se1, se2)
dim(cmb2)
dimnames(cmb2)
}
