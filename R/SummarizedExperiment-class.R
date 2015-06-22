### =========================================================================
### SummarizedExperiment0 and RangedSummarizedExperiment objects
### -------------------------------------------------------------------------
###

### TODO: Rename this class SummarizedExperiment once the "old"
### SummarizedExperiment class in GenomicRanges is gone.
setClass("SummarizedExperiment0",
    contains="Vector",
    representation(
        colData="DataFrame",              # columns and their annotations
        assays="Assays",                  # Data -- e.g., list of matricies
        NAMES="characterORNULL",
        elementMetadata="DataFrame"
    ),
    prototype(
        assays=GenomicRanges:::.ShallowSimpleListAssays(data=SimpleList())
    )
)

### Combine the new parallel slots with those of the parent class. Make sure
### to put the new parallel slots *first*.
setMethod("parallelSlotNames", "SummarizedExperiment0",
    function(x) c("NAMES", callNextMethod())
)

### The 'elementMetadata' slot must contain a zero-column DataFrame at all time
### (this is checked by the validity method). The top-level mcols are stored on
### the rowRanges component.
setClass("RangedSummarizedExperiment",
    contains="SummarizedExperiment0", 
    representation(
        rowRanges="GenomicRangesORGRangesList"
    ),
    prototype(
        rowRanges=GRanges()
    )
)

### Combine the new parallel slots with those of the parent class. Make sure
### to put the new parallel slots *first*.
setMethod("parallelSlotNames", "RangedSummarizedExperiment",
    function(x) c("rowRanges", callNextMethod())
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Validity.
###

.valid.SummarizedExperiment0.assays_current <- function(x)
{
    if (!is(slot(x, "assays"), "Assays"))
        return("'assays' is out-of-date; use updateObject()")
    NULL
}

.valid.SummarizedExperiment0.assays_class <- function(x)
{
    ok <- sapply(assays(x, withDimnames=FALSE), function(cl) {
        (!is.null(dim(cl))) && (length(dim(cl)) >= 2L)
    })
    if (!all(ok))
        return("'assays' must be matrix-like with 2 (or more?) dimensions")
    NULL
}

.valid.SummarizedExperiment0.assays_nrow <- function(x)
{
    if (!all(sapply(assays(x, withDimnames=FALSE), nrow) ==
             length(x)))
        return("'mcols' nrow differs from 'assays' nrow")
    NULL
}

.valid.SummarizedExperiment0.assays_ncol <- function(x)
{
    if (!all(sapply(assays(x, withDimnames=FALSE), ncol) ==
             nrow(colData(x))))
        return("'colData' nrow differs from 'assays' ncol")
    NULL
}

.valid.SummarizedExperiment0.assays_dim <- function(x)
{
    c(.valid.SummarizedExperiment0.assays_nrow(x),
      .valid.SummarizedExperiment0.assays_ncol(x))
}

.valid.SummarizedExperiment0 <- function(x)
{
    c(.valid.SummarizedExperiment0.assays_current(x),
      msg <- .valid.SummarizedExperiment0.assays_class(x),
      if (is.null(msg)) {
          .valid.SummarizedExperiment0.assays_dim(x)
      } else NULL)
}

setValidity2("SummarizedExperiment0", .valid.SummarizedExperiment0)

.valid.RangedSummarizedExperiment <- function(x)
{
    if (ncol(x@elementMetadata) != 0L)
        return("'elementMetadata' slot must contain a zero-column DataFrame")
    NULL
}

setValidity2("RangedSummarizedExperiment", .valid.RangedSummarizedExperiment)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Constructor.
###

.new_RangedSummarizedExperiment <- function(rowRanges, colData, assays,
                                            metadata)
{
    elementMetadata <- new("DataFrame", nrows=length(rowRanges))
    new("RangedSummarizedExperiment", rowRanges=rowRanges,
                                      colData=colData,
                                      assays=assays,
                                      elementMetadata=elementMetadata,
                                      metadata=as.list(metadata))
}

.GRangesList_assays <-
    function(assays)
{
    m <- assays[[1]]
    n <- nrow(m)
    names <- rownames(m)
    relist(GRanges(), PartitioningByEnd(integer(n), names=names))
}

setMethod(SummarizedExperiment, "SimpleList",
   function(assays, rowRanges=GRangesList(), colData=DataFrame(),
            metadata=list(), exptData=SimpleList())
{
    if (missing(rowRanges) && 0L != length(assays)) {
        rowRanges <- .GRangesList_assays(assays)
    }

    if (missing(colData) && 0L != length(assays)) {
        nms <- colnames(assays[[1]])
        if (is.null(nms) && 0L != ncol(assays[[1]]))
            stop("'SummarizedExperiment' assay colnames must not be NULL")
        colData <- DataFrame(row.names=nms)
    }

    FUN <- function(x) {
        exp <- list(names(rowRanges), rownames(colData))
        ## dimnames as NULL or list(NULL, NULL)
        all(sapply(dimnames(x), is.null)) ||
            ## or consistent with row / colData
            identical(dimnames(x)[1:2], exp)
    }
    if (!all(sapply(assays, FUN)))
        assays <- endoapply(assays, unname)
    if (!is(assays, "Assays"))
        assays <- GenomicRanges:::.ShallowSimpleListAssays(data=assays)

    ## For backward compatibility with "classic" SummarizedExperiment objects.
    if (!missing(exptData)) {
        if (!missing(metadata))
            stop("only one of 'metadata' and 'exptData' can be ",
                 "specified, but not both")
        msg <- c("the 'exptData' argument is deprecated, ",
                 "please use 'metadata' instead")
        .Deprecated(msg=msg)
        metadata <- exptData
    }
    .new_RangedSummarizedExperiment(rowRanges, colData, assays, metadata)
})

setMethod(SummarizedExperiment, "missing",
    function(assays, ...)
{
    SummarizedExperiment(SimpleList(), ...)
})

setMethod(SummarizedExperiment, "list",
    function(assays, ...)
{
    SummarizedExperiment(do.call(SimpleList, assays), ...)
})

setMethod(SummarizedExperiment, "matrix",
    function(assays, ...)
{
    if (is.list(assays))
        ## special case -- matrix of lists
        assays <- list(assays)
    SummarizedExperiment(SimpleList(assays), ...)
})


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Coercion.
###
### See makeSummarizedExperimentFromExpressionSet.R for coercion back and
### forth between SummarizedExperiment and ExpressionSet.
###

.from_RangedSummarizedExperiment_to_SummarizedExperiment0 <- function(from)
{
    ans_NAMES <- names(from@rowRanges)
    ans_elementMetadata <- mcols(from@rowRanges)
    new("SummarizedExperiment0", colData=from@colData,
                                 assays=from@assays,
                                 NAMES=ans_NAMES,
                                 elementMetadata=ans_elementMetadata,
                                 metadata=from@metadata)
}

setAs("RangedSummarizedExperiment", "SummarizedExperiment0",
    .from_RangedSummarizedExperiment_to_SummarizedExperiment0
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Clone.
###

setMethod(GenomicRanges:::clone, "SummarizedExperiment0",  # not exported
    function(x, ...)
{
    ## S4Vectors:::extraArgsAsList would prevent using clone on
    ## subclasses
    args <- list(...)
    firstTime <- TRUE
    for (nm in names(args)) {
        s <- slot(x, nm)
        v <- args[[nm]]
        if (is(s, "ShallowData"))
            v <- GenomicRanges:::clone(s, data=v)
        if (firstTime) {
            slot(x, nm, FALSE) <- v
            firstTime <- FALSE
        } else {
            `slot<-`(x, nm, FALSE, v)
        }
    }
    x
})


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### SummarizedExperiment0 getters and setters.
###

setMethod("length", "SummarizedExperiment0",
    function(x) nrow(x@elementMetadata)
)

setMethod("names", "SummarizedExperiment0", function(x) x@NAMES)

setReplaceMethod("names", "SummarizedExperiment0",
    function(x, value)
    {
        value <- S4Vectors:::normalize_names_replacement_value(value, x)
        GenomicRanges:::clone(x, NAMES=value)
    }
)

### We override the default "metadata<-" method (defined in S4Vectors for
### Annotated objects) with more memory efficient version that uses clone()
### in order to minimize data copy (ShallowData trick).
setReplaceMethod("metadata", "SummarizedExperiment0",
    function(x, value)
{
    value <- as.list(value)
    GenomicRanges:::clone(x, metadata=value)
})

setGeneric("value",  # not exported
    function(x, name, ...) standardGeneric("value"),
    signature = "x")

setMethod("value", "SummarizedExperiment0",  # not exported
    function(x, name, ...)
{
    s <- slot(x, name)
    if (is(s, "ShallowData"))
        s <- s$data
    s
})

setMethod(colData, "SummarizedExperiment0",
    function(x, ...) value(x, "colData"))

setReplaceMethod("colData", c("SummarizedExperiment0", "DataFrame"),
    function(x, ..., value)
{
    if (nrow(value) != ncol(x))
        stop("nrow of supplied 'colData' must equal ncol of object")
    GenomicRanges:::clone(x, colData=value)
})

setMethod("assayNames", "SummarizedExperiment0",
    function(x, ...)
{
    names(assays(x, withDimnames=FALSE))
})

setMethod("assayNames<-", c("SummarizedExperiment0", "character"),
    function(x, ..., value)
{
    names(assays(x, withDimnames=FALSE)) <- value
    x
})

setMethod(assays, "SummarizedExperiment0",
    function(x, ..., withDimnames=TRUE)
{
    if (withDimnames)
        endoapply(value(x, "assays"), "dimnames<-", dimnames(x))
    else
        value(x, "assays")
})

.SummarizedExperiment.assays.replace <-
    function(x, ..., withDimnames=TRUE, value)
{
    ## withDimnames arg allows names(assays(se, withDimnames=FALSE)) <- value
    ok <- vapply(value, function(elt, xdimnames) {
        e <- dimnames(elt)
        (is.null(e[[1]]) || identical(e[[1]], xdimnames[[1]])) &&
             (is.null(e[[2]]) || identical(e[[2]], xdimnames[[2]]))
    }, logical(1), xdimnames=dimnames(x))
    if (!all(ok))
        stop("current and replacement dimnames() differ")
    x <- GenomicRanges:::clone(x, assays=value)
    ## validObject(x) should be called below because it would then fully
    ## re-validate objects that derive from SummarizedExperiment0 (e.g.
    ## DESeqDataSet objects) after the user sets the assays slot with
    ## assays(x) <- value. For example the assays slot of a DESeqDataSet
    ## object must contain a matrix named 'counts' and calling validObject(x)
    ## would check that but .valid.SummarizedExperiment0(x) doesn't.
    ## The FourC() constructor function defined in the FourCSeq package
    ## actually takes advantage of the incomplete validation below to
    ## purposedly return invalid FourC objects!
    msg <- .valid.SummarizedExperiment0(x)
    if (!is.null(msg)) 
        stop(msg)
    x
}

setReplaceMethod("assays", c("SummarizedExperiment0", "SimpleList"),
    .SummarizedExperiment.assays.replace)

setReplaceMethod("assays", c("SummarizedExperiment0", "list"),
    .SummarizedExperiment.assays.replace)

## convenience for common use case
setMethod(assay, c("SummarizedExperiment0", "missing"),
    function(x, i, ...)
{
    assays <- assays(x, ...)
    if (0L == length(assays))
        stop("'assay(<", class(x), ">, i=\"missing\", ...) ",
             "length(assays(<", class(x), ">)) is 0'")
    assays[[1]]
})

setMethod(assay, c("SummarizedExperiment0", "numeric"),
    function(x, i, ...)
{
    tryCatch({
        assays(x, ...)[[i]]
    }, error=function(err) {
        stop("'assay(<", class(x), ">, i=\"numeric\", ...)' ",
             "invalid subscript 'i'\n", conditionMessage(err))
    })
})

setMethod(assay, c("SummarizedExperiment0", "character"),
    function(x, i = names(x)[1], ...)
{
    msg <- paste0("'assay(<", class(x), ">, i=\"character\", ...)' ",
                  "invalid subscript 'i'")
    res <- tryCatch({
        assays(x, ...)[[i]]
    }, error=function(err) {
        stop(msg, "\n", conditionMessage(err))
    })
    if (is.null(res))
        stop(msg, "\n'i' not in names(assays(<", class(x), ">))")
    res
})

setReplaceMethod("assay", c("SummarizedExperiment0", "missing", "matrix"),
    function(x, i, ..., value)
{
    if (0L == length(assays(x)))
        stop("'assay(<", class(x), ">) <- value' ", "length(assays(<",
             class(x), ">)) is 0")
    assays(x)[[1]] <- value
    x
})

setReplaceMethod("assay",
    c("SummarizedExperiment0", "numeric", "matrix"),
    function(x, i = 1, ..., value)
{
    assays(x, ...)[[i]] <- value
    x
})

setReplaceMethod("assay",
    c("SummarizedExperiment0", "character", "matrix"),
    function(x, i = names(x)[1], ..., value)
{
    assays(x, ...)[[i]] <- value
    x
})

## cannonical location for dim, dimnames; dimnames should be checked
## for consistency (if non-null) and stripped from assays on
## construction, or added from assays if row/col names are NULL in
## <SummarizedExperiment0> but not assays. dimnames need to be added on
## to assays when assays() invoked
setMethod(dim, "SummarizedExperiment0",
    function(x)
{
    c(length(x), nrow(colData(x)))
})

setMethod(dimnames, "SummarizedExperiment0",
    function(x)
{
    list(names(x), rownames(colData(x)))
})


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### RangedSummarizedExperiment getters and setters.
###

setMethod(rowRanges, "RangedSummarizedExperiment",
    function(x, ...) value(x, "rowRanges"))

.RangedSummarizedExperiment.rowRanges.replace <-
    function(x, ..., value)
{
    x <- GenomicRanges:::clone(
            x, ...,
            rowRanges=value,
            elementMetadata=new("DataFrame", nrows=length(value)))
    msg <- .valid.SummarizedExperiment0.assays_nrow(x)
    if (!is.null(msg))
        stop(msg)
    x
}

setReplaceMethod("rowRanges", c("RangedSummarizedExperiment", "GenomicRanges"),
    .RangedSummarizedExperiment.rowRanges.replace)

setReplaceMethod("rowRanges", c("RangedSummarizedExperiment", "GRangesList"),
    .RangedSummarizedExperiment.rowRanges.replace)

setMethod("names", "RangedSummarizedExperiment",
    function(x) names(rowRanges(x))
)

setReplaceMethod("names", "RangedSummarizedExperiment",
    function(x, value)
{
    rowRanges <- rowRanges(x)
    names(rowRanges) <- value
    GenomicRanges:::clone(x, rowRanges=rowRanges)
})

setMethod(dimnames, "RangedSummarizedExperiment",
    function(x)
{
    list(names(x), rownames(colData(x)))
})

setReplaceMethod("dimnames", c("RangedSummarizedExperiment", "list"),
    function(x, value)
{
    rowRanges <- rowRanges(x)
    names(rowRanges) <- value[[1]]
    colData <- colData(x)
    rownames(colData) <- value[[2]]
    GenomicRanges:::clone(x, rowRanges=rowRanges, colData=colData)
})

setReplaceMethod("dimnames", c("RangedSummarizedExperiment", "NULL"),
    function(x, value)
{
    dimnames(x) <- list(NULL, NULL)
    x
})


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Subsetting.
###

.SummarizedExperiment.charbound <-
    function(idx, txt, fmt)
{
    orig <- idx
    idx <- match(idx, txt)
    if (any(bad <- is.na(idx))) {
        msg <- paste(BiocGenerics:::selectSome(orig[bad]), collapse=" ")
        stop(sprintf(fmt, msg))
    }
    idx
}

.SummarizedExperiment.assays.subset <- function(x, i, j)
{
    ## need to expand Rle's for subsetting standard matrix
    if (!missing(i) && !missing(j)) {
        fun <- function(x) {
            switch(length(dim(x)),
                   stop("'[' on assays() with 1 dimension not supported"),
                   x[i, j, drop=FALSE],
                   x[i, j, , drop=FALSE],
                   x[i, j, , , drop=FALSE],
                   stop("'[' on assays() with >4 dimensions not supported"))
        }
    } else if (!missing(i)) {
        fun <- function(x) {
            switch(length(dim(x)),
                   stop("'[' on assays() with 1 dimension not supported"),
                   x[i, , drop=FALSE],
                   x[i, , , drop=FALSE],
                   x[i, , , , drop=FALSE],
                   stop("'[' on assays() with >4 dimensions not supported"))
        }
    } else if (!missing(j)) {
        fun <- function(x) {
            switch(length(dim(x)),
                   stop("'[' on assays() with 1 dimension not supported"),
                   x[, j, drop=FALSE],
                   x[, j, , drop=FALSE],
                   x[, j, , , drop=FALSE],
                   stop("'[' on assays() with >4 dimensions not supported"))
        }
    }
    endoapply(assays(x, withDimnames=FALSE), fun)
}

setMethod("[", c("RangedSummarizedExperiment", "ANY", "ANY"),
    function(x, i, j, ..., drop=TRUE)
{
    if (1L != length(drop) || (!missing(drop) && drop))
        warning("'drop' ignored '[,", class(x), ",ANY,ANY-method'")

    if (missing(i) && missing(j))
        return(x)

    if (!missing(i)) {
        if (is.character(i)) {
            fmt <- paste0("<", class(x), ">[i,] index out of bounds: %s")
            i <- .SummarizedExperiment.charbound(i, rownames(x), fmt)
        }
        ans_elementMetadata <- x@elementMetadata[i, , drop=FALSE]
        ans_rowRanges <- rowRanges(x)[i]
        ii <- as.vector(i)
    }
    if (!missing(j)) {
        if (is.character(j)) {
            fmt <- paste0("<", class(x), ">[,j] index out of bounds: %s")
            j <- .SummarizedExperiment.charbound(j, colnames(x), fmt)
        }
        ans_colData <- x@colData[j, , drop=FALSE]
        jj <- as.vector(j)
    }

    if (missing(i)) {
        ans_assays <- .SummarizedExperiment.assays.subset(x, j=jj)
        ans <- GenomicRanges:::clone(x, ...,
                                        colData=ans_colData,
                                        assays=ans_assays)
    } else if (missing(j)) {
        ans_assays <- .SummarizedExperiment.assays.subset(x, i=ii)
        ans <- GenomicRanges:::clone(x, ...,
                                        elementMetadata=ans_elementMetadata,
                                        rowRanges=ans_rowRanges,
                                        assays=ans_assays)
    } else {
        ans_assays <- .SummarizedExperiment.assays.subset(x, ii, jj)
        ans <- GenomicRanges:::clone(x, ...,
                                        elementMetadata=ans_elementMetadata,
                                        rowRanges=ans_rowRanges,
                                        colData=ans_colData,
                                        assays=ans_assays)
    }
    ans
})

.SummarizedExperiment.assays.subsetgets <- function(x, i, j, value)
{
    ## need to expand Rle's for subsetting standard matrix
    if (!missing(i) && !missing(j)) {
        fun <- function(x, value) {
            switch(length(dim(x)),
                   stop("'[<-' on assays() with 1 dimension not supported"),
                   x[i, j] <- value,
                   x[i, j, ] <- value,
                   x[i, j, , ] <- value,
                   stop("'[<-' on assays() with >4 dimensions not supported"))
            x
        }
    } else if (!missing(i)) {
        fun <- function(x, value) {
            switch(length(dim(x)),
                   stop("'[<-' on assays() with 1 dimension not supported"),
                   x[i, ] <- value,
                   x[i, , ] <- value,
                   x[i, , , ] <- value,
                   stop("'[<-' on assays() with >4 dimensions not supported"))
            x
        }
    } else if (!missing(j)) {
        fun <- function(x, value) {
            switch(length(dim(x)),
                   stop("'[<-' on assays() with 1 dimension not supported"),
                   x[, j] <- value,
                   x[, j, ] <- value,
                   x[, j, , ] <- value,
                   stop("'[<-' on assays() with >4 dimensions not supported"))
            x
        }
    }
    a <- assays(x, withDimnames=FALSE)
    v <- assays(value, withDimnames=FALSE)
    mendoapply(fun, x=a, value=v)
}

setReplaceMethod("[",
    c("RangedSummarizedExperiment", "ANY", "ANY", "RangedSummarizedExperiment"),
    function(x, i, j, ..., value)
{
    if (missing(i) && missing(j))
        return(value)

    ans_metadata <- c(metadata(x), metadata(value))

    if (!missing(i)) {
        if (is.character(i)) {
            fmt <- paste0("<", class(x), ">[i,] index out of bounds: %s")
            i <- .SummarizedExperiment.charbound(i, rownames(x), fmt)
        }
        ii <- as.vector(i)
        ans_elementMetadata <- local({
            emd <- x@elementMetadata
            emd[i,] <- value@elementMetadata
            emd
        })
        ans_rowRanges <- local({
            r <- rowRanges(x)
            r[i] <- rowRanges(value)
            names(r)[ii] <- names(rowRanges(value))
            r
        })
    }

    if (!missing(j)) {
        if (is.character(j)) {
            fmt <- paste0("<", class(x), ">[,j] index out of bounds: %s")
            j <- .SummarizedExperiment.charbound(j, colnames(x), fmt)
        }
        jj <- as.vector(j)
        ans_colData <- local({
            c <- x@colData
            c[j,] <- value@colData
            rownames(c)[jj] <- rownames(value@colData)
            c
        })
    }

    if (missing(i)) {
        ans_assays <- .SummarizedExperiment.assays.subsetgets(x, j=jj,
                                                              value=value)
        ans <- GenomicRanges:::clone(x, ...,
                                        metadata=ans_metadata,
                                        colData=ans_colData,
                                        assays=ans_assays)
        msg <- .valid.SummarizedExperiment0.assays_ncol(ans)
    } else if (missing(j)) {
        ans_assays <- .SummarizedExperiment.assays.subsetgets(x, i=ii,
                                                              value=value)
        ans <- GenomicRanges:::clone(x, ...,
                                        metadata=ans_metadata,
                                        elementMetadata=ans_elementMetadata,
                                        rowRanges=ans_rowRanges,
                                        assays=ans_assays)
        msg <- .valid.SummarizedExperiment0.assays_nrow(ans)
    } else {
        ans_assays <- .SummarizedExperiment.assays.subsetgets(x, ii, jj,
                                                              value=value)
        ans <- GenomicRanges:::clone(x, ...,
                                        metadata=ans_metadata,
                                        elementMetadata=ans_elementMetadata,
                                        rowRanges=ans_rowRanges,
                                        colData=ans_colData,
                                        assays=ans_assays)
        msg <- .valid.SummarizedExperiment0.assays_dim(ans)
    }
    if (!is.null(msg))
        stop(msg)
    ans
})


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Quick colData access.
###

setMethod("[[", c("SummarizedExperiment0", "ANY", "missing"),
    function(x, i, j, ...)
{
    colData(x)[[i, ...]]
})

setReplaceMethod("[[", c("SummarizedExperiment0", "ANY", "missing"),
    function(x, i, j, ..., value)
{
    colData(x)[[i, ...]] <- value
    x
})

.DollarNames.SummarizedExperiment0 <- function(x, pattern)
    grep(pattern, names(colData(x)), value=TRUE)

setMethod("$", "SummarizedExperiment0",
    function(x, name)
{
    colData(x)[[name]]
})

setReplaceMethod("$", "SummarizedExperiment0",
    function(x, name, value)
{
    colData(x)[[name]] <- value
    x
})


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Display.
###

setMethod(show, "SummarizedExperiment0",
    function(object)
{
    selectSome <- BiocGenerics:::selectSome
    scat <- function(fmt, vals=character(), exdent=2, ...)
    {
        vals <- ifelse(nzchar(vals), vals, "''")
        lbls <- paste(BiocGenerics:::selectSome(vals), collapse=" ")
        txt <- sprintf(fmt, length(vals), lbls)
        cat(strwrap(txt, exdent=exdent, ...), sep="\n")
    }

    cat("class:", class(object), "\n")
    cat("dim:", dim(object), "\n")

    ## metadata()
    expt <- names(metadata(object))
    if (is.null(expt))
        expt <- character(length(metadata(object)))
    scat("metadata(%d): %s\n", expt)

    ## assays()
    nms <- names(assays(object, withDimnames=FALSE))
    if (is.null(nms))
        nms <- character(length(assays(object, withDimnames=FALSE)))
    scat("assays(%d): %s\n", nms)

    ## rownames()
    dimnames <- dimnames(object)
    dlen <- sapply(dimnames, length)
    if (dlen[[1]]) scat("rownames(%d): %s\n", dimnames[[1]])
    else scat("rownames: NULL\n")

    ## mcols()
    mcolnames <- names(mcols(object))
    fmt <- "metadata column names(%d): %s\n"
    if (is(object, "RangedSummarizedExperiment"))
        fmt <- paste("rowRanges", fmt)
    scat("rowRanges metadata column names(%d): %s\n", mcolnames)

    ## colnames()
    if (dlen[[2]]) scat("colnames(%d): %s\n", dimnames[[2]])
    else cat("colnames: NULL\n")

    ## colData()
    scat("colData names(%d): %s\n", names(colData(object)))
})


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Combine.
###

### Appropriate for objects with different ranges and same samples.
setMethod("rbind", "RangedSummarizedExperiment",
    function(..., deparse.level=1)
{
    args <- unname(list(...))
    .rbind.RangedSummarizedExperiment(args)
})

.rbind.RangedSummarizedExperiment <- function(args)
{
    if (!.compare(lapply(args, colnames)))
            stop("'...' objects must have the same colnames")
    if (!.compare(lapply(args, ncol)))
            stop("'...' objects must have the same number of samples")

    rowRanges <- do.call(c, lapply(args, rowRanges))
    colData <- .cbind.DataFrame(args, colData, "colData")
    assays <- .bind.arrays(args, rbind, "assays")
    elementMetadata <- do.call(rbind, lapply(args, slot, "elementMetadata"))
    metadata <- do.call(c, lapply(args, metadata))

    initialize(args[[1L]],
               rowRanges=rowRanges, colData=colData, assays=assays,
               elementMetadata=elementMetadata, metadata=metadata)
}

### Appropriate for objects with same ranges and different samples.
setMethod("cbind", "RangedSummarizedExperiment",
    function(..., deparse.level=1)
{
    args <- unname(list(...))
    .cbind.RangedSummarizedExperiment(args)
})

.cbind.RangedSummarizedExperiment <- function(args)
{
    if (!.compare(lapply(args, rowRanges), TRUE))
        stop("'...' object ranges (rows) are not compatible")

    rowRanges <- rowRanges(args[[1L]])
    mcols(rowRanges) <- .cbind.DataFrame(args, mcols, "mcols")
    colData <- do.call(rbind, lapply(args, colData))
    assays <- .bind.arrays(args, cbind, "assays")
    metadata <- do.call(c, lapply(args, metadata))

    initialize(args[[1L]],
               rowRanges=rowRanges, colData=colData, assays=assays,
               metadata=metadata)
}

.compare <- function(x, GenomicRanges=FALSE)
{
    x1 <- x[[1]]
    if (GenomicRanges) {
        if (is(x1, "GRangesList")) {
            x <- lapply(x, unlist)
            x1 <- x[[1]]
        }
        for (i in seq_along(x)[-1]) {
            if (length(x1) != length(x[[i]]))
                return(FALSE)
            ok <- x1 == x[[i]]
            if (!all(ok))
                return(FALSE)
        }
        return(TRUE)
    } else {
        all(sapply(x[-1],
            function(xelt) all(identical(xelt, x[[1]]))))
    }
}

.cbind.DataFrame <- function(args, accessor, accessorName)
{
    lst <- lapply(args, accessor)
    if (!.compare(lst)) {
        nms <- lapply(lst, names)
        nmsv <- unlist(nms, use.names=FALSE)
        names(nmsv) <- rep(seq_along(nms), elementLengths(nms))
        dups <- duplicated(nmsv)
        ## no duplicates
        if (!any(dups))
            return(do.call(cbind, lst))
        ## confirm duplicates are the same
        lapply(nmsv[duplicated(nmsv)], function(d) {
            if (!.compare(lapply(lst, "[", d)))
                stop("column(s) '", unname(d),
                     "' in ", sQuote(accessorName),
                     " are duplicated and the data do not match")})
        ## remove duplicates
        do.call(cbind, lst)[,!dups]
    } else {
        lst[[1]]
    }
}

.bind.array.elements <- function(index, lst, bind) {
    e1 <- lapply(lst, "[[", index)
    dim <- .get.assay.dimension(e1, bind)
    if (is.na(dim[3])) {
        do.call(bind, e1)
    } else {
        e2 <- lapply(1:dim[3], function(i) {
            do.call(bind, lapply(e1, "[", ,,i))
        })
        array(do.call(c, e2), dim=dim)
    }
}

.bind.arrays <- function(args, bind, accessor)
{
    lst <- lapply(args, accessor)
    if (!length(elts <- unique(elementLengths(lst))) == 1L)
        stop("elements in ", sQuote(accessor),
             " must have the same length")
    if (elts == 0L)
        return(GenomicRanges:::.ShallowSimpleListAssays(data=SimpleList()))
    var <- lapply(lst,  names)
    if (is.null(uvar <- unique(unlist(var)))) {
        ## no names, match by position
        res <- lapply(seq_along(elts), .bind.array.elements, lst=lst, bind=bind)
    } else {
        ## match by name
        if (!.compare(var))
            stop("elements in ", sQuote(accessor),
                 " must have the same names")
        res <- lapply(uvar, .bind.array.elements, lst=lst, bind=bind)
        names(res) <- uvar
    }
    GenomicRanges:::.ShallowSimpleListAssays(data=SimpleList(res))
}

.get.assay.dimension <- function(lst, bind)
{
    dim <- lapply(lst, dim)
    if (!.compare(lapply(dim, "[", 3)))
        stop("elements in assays must have the same dimension")
    if (identical(bind, cbind))
        c(dim[[1]][1], do.call(sum, lapply(dim, "[", 2)), dim[[1]][3])
    else
        c(do.call(sum, lapply(dim, "[", 1)), dim[[1]][2], dim[[1]][3])
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### To facilitate transition from "classic" SummarizedExperiment objects to
### new SummarizedExperiment0/RangedSummarizedExperiment objects.
###

### For backward compatibility with "classic" SummarizedExperiment objects.
setMethod("exptData", "SummarizedExperiment0",
    function(x, ...)
{
    .Deprecated("metadata")
    SimpleList(metadata(x, ...))
})

setReplaceMethod("exptData", "SummarizedExperiment0",
    function(x, ..., value)
{
    .Deprecated("metadata<-")
    `metadata<-`(x, ..., value=value)
})

### Update "classic" SummarizedExperiment objects.

.has_SummarizedExperiment_internal_structure <- function(object)
    all(sapply(slotNames("SummarizedExperiment"), .hasSlot, object=object))

### Used in GenomicRanges!
.from_SummarizedExperiment_to_RangedSummarizedExperiment <- function(from)
    .new_RangedSummarizedExperiment(from@rowData,
                                    from@colData,
                                    from@assays,
                                    from@exptData)

### Should work on any object that (1) belongs to a class that now derives
### from RangedSummarizedExperiment and (2) was created at a time when it was
### deriving from the old SummarizedExperiment class defined in GenomicRanges.
### For example: DESeqDataSet objects created before the definition of class
### DESeqDataSet was modified to extend RangedSummarizedExperiment instead of
### SummarizedExperiment.  
setMethod(updateObject, "RangedSummarizedExperiment",
    function(object, ..., verbose=FALSE)
{
    if (.has_SummarizedExperiment_internal_structure(object)) {
        rse <- .from_SummarizedExperiment_to_RangedSummarizedExperiment(object)
    } else if (!(.hasSlot(object, "NAMES") &&
                 .hasSlot(object, "elementMetadata"))) {
        rse <- .new_RangedSummarizedExperiment(object@rowRanges,
                                               object@colData,
                                               object@assays,
                                               object@metadata)
    } else {
        return(object)
    }

    xslotnames <- setdiff(slotNames(class(object)), slotNames(class(rse)))
    xslots <- attributes(object)[xslotnames]
    #do.call("new", c(list(Class=class(object), rse), xslots))

    ## The line above doesn't work because of a bug in R (see
    ## https://stat.ethz.ch/pipermail/r-devel/2015-May/071130.html),
    ## so we use the workaround below.
    rse_slots <- attributes(rse)[slotNames(class(rse))]

    ## Because of another bug in R, rse_slots$NAMES is broken when the NAMES
    ## slot is NULL so we repair it (see
    ## https://bugs.r-project.org/bugzilla/show_bug.cgi?id=16428).
    if (is.name(rse_slots$NAMES))
        rse_slots$NAMES <- NULL

    do.call("new", c(list(Class=class(object)), rse_slots, xslots))
})

