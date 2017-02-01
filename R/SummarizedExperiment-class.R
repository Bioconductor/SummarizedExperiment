### =========================================================================
### SummarizedExperiment objects
### -------------------------------------------------------------------------
###


setClass("SummarizedExperiment",
    contains="Vector",
    representation(
        colData="DataFrame",              # columns and their annotations
        assays="Assays",                  # Data -- e.g., list of matricies
        NAMES="character_OR_NULL",
        elementMetadata="DataFrame"
    ),
    prototype(
        assays=Assays()
    )
)

### Combine the new parallel slots with those of the parent class. Make sure
### to put the new parallel slots *first*.
setMethod("parallelSlotNames", "SummarizedExperiment",
    function(x) c("NAMES", callNextMethod())
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Validity.
###

.valid.SummarizedExperiment.assays_nrow <- function(x)
{
    if (length(x@assays) == 0L)
        return(NULL)
    assays_nrow <- nrow(x@assays)
    rowData_nrow <- length(x)
    if (assays_nrow != rowData_nrow) {
        txt <- sprintf(
            "\n  nb of rows in 'assay' (%d) must equal nb of rows in 'rowData' (%d)",
            assays_nrow, rowData_nrow)
        return(txt)
    }
    NULL
}

.valid.SummarizedExperiment.assays_ncol <- function(x)
{
    if (length(x@assays) == 0L)
        return(NULL)
    assays_ncol <- ncol(x@assays)
    colData_nrow <- nrow(colData(x))
    if (assays_ncol != colData_nrow) {
        txt <- sprintf(
            "\n  nb of cols in 'assay' (%d) must equal nb of rows in 'colData' (%d)",
            assays_ncol, colData_nrow)
        return(txt)
    }
    NULL
}

.valid.SummarizedExperiment.assays_dim <- function(x)
{
    c(.valid.SummarizedExperiment.assays_nrow(x),
      .valid.SummarizedExperiment.assays_ncol(x))
}

.valid.SummarizedExperiment <- function(x)
{
    .valid.SummarizedExperiment.assays_dim(x)
}

setValidity2("SummarizedExperiment", .valid.SummarizedExperiment)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Low-level constructor (not exported).
###

new_SummarizedExperiment <- function(assays, names, rowData, colData,
                                     metadata)
{
    if (!is(assays, "Assays"))
        assays <- Assays(assays)
    if (is.null(rowData)) {
        if (is.null(names))
            nrow <- nrow(assays)
        else
            nrow <- length(names)
        rowData <- S4Vectors:::make_zero_col_DataFrame(nrow)
    } else {
        rownames(rowData) <- NULL
    }
    new("SummarizedExperiment", NAMES=names,
                                elementMetadata=rowData,
                                colData=colData,
                                assays=assays,
                                metadata=as.list(metadata))
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Getters and setters.
###

setMethod("length", "SummarizedExperiment",
    function(x) nrow(x@elementMetadata)
)

setMethod("names", "SummarizedExperiment", function(x) x@NAMES)

setReplaceMethod("names", "SummarizedExperiment",
    function(x, value)
    {
        NAMES <- S4Vectors:::normalize_names_replacement_value(value, x)
        BiocGenerics:::replaceSlots(x, NAMES=NAMES, check=FALSE)
    }
)

## rowData, colData seem too vague, but from eSet derived classes wanted to
## call the rows / cols something different from 'features' or 'samples', so
## might as well avoid the issue

setGeneric("rowData", function(x, ...) standardGeneric("rowData"))

setMethod("rowData", "SummarizedExperiment",
    function(x, ...) mcols(x, ...)
)

setGeneric("rowData<-",
    function(x, ..., value) standardGeneric("rowData<-"))

setReplaceMethod("rowData", "SummarizedExperiment",
    function(x, ..., value) `mcols<-`(x, ..., value=value)
)

setGeneric("colData", function(x, ...) standardGeneric("colData"))

setMethod("colData", "SummarizedExperiment", function(x, ...) x@colData)

setGeneric("colData<-",
    function(x, ..., value) standardGeneric("colData<-"))

setReplaceMethod("colData", c("SummarizedExperiment", "DataFrame"),
    function(x, ..., value)
{
    if (nrow(value) != ncol(x))
        stop("nrow of supplied 'colData' must equal ncol of object")
    BiocGenerics:::replaceSlots(x, colData=value, check=FALSE)
})

setGeneric("assays",
    function(x, ..., withDimnames=TRUE) standardGeneric("assays"),
    signature="x")

setMethod("assays", "SummarizedExperiment",
    function(x, ..., withDimnames=TRUE)
{
    assays <- as(x@assays, "SimpleList")
    if (withDimnames) {
        assays <- endoapply(assays,
            function(assay) {
                dimnames(assay)[1:2] <- dimnames(x)
                assay
            }
        )
    }
    assays
})

setGeneric("assays<-",
    function(x, ..., withDimnames=TRUE, value) standardGeneric("assays<-"),
    signature=c("x", "value"))

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
    x <- BiocGenerics:::replaceSlots(x, assays=Assays(value), check=FALSE)
    ## validObject(x) should be called below because it would then fully
    ## re-validate objects that derive from SummarizedExperiment (e.g.
    ## DESeqDataSet objects) after the user sets the assays slot with
    ## assays(x) <- value. For example the assays slot of a DESeqDataSet
    ## object must contain a matrix named 'counts' and calling validObject(x)
    ## would check that but .valid.SummarizedExperiment(x) doesn't.
    ## The FourC() constructor function defined in the FourCSeq package
    ## actually takes advantage of the incomplete validation below to
    ## purposedly return invalid FourC objects!
    msg <- .valid.SummarizedExperiment(x)
    if (!is.null(msg)) 
        stop(msg)
    x
}

setReplaceMethod("assays", c("SummarizedExperiment", "SimpleList"),
    .SummarizedExperiment.assays.replace)

setReplaceMethod("assays", c("SummarizedExperiment", "list"),
    .SummarizedExperiment.assays.replace)

setGeneric("assay", function(x, i, ...) standardGeneric("assay"))

## convenience for common use case
setMethod("assay", c("SummarizedExperiment", "missing"),
    function(x, i, ...)
{
    assays <- assays(x, ...)
    if (0L == length(assays))
        stop("'assay(<", class(x), ">, i=\"missing\", ...) ",
             "length(assays(<", class(x), ">)) is 0'")
    assays[[1]]
})

setMethod("assay", c("SummarizedExperiment", "numeric"),
    function(x, i, ...)
{
    tryCatch({
        assays(x, ...)[[i]]
    }, error=function(err) {
        stop("'assay(<", class(x), ">, i=\"numeric\", ...)' ",
             "invalid subscript 'i'\n", conditionMessage(err))
    })
})

setMethod("assay", c("SummarizedExperiment", "character"),
    function(x, i, ...)
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

setGeneric("assay<-", signature=c("x", "i"),
    function(x, i, ..., value) standardGeneric("assay<-"))

setReplaceMethod("assay", c("SummarizedExperiment", "missing"),
    function(x, i, ..., value)
{
    if (0L == length(assays(x)))
        stop("'assay(<", class(x), ">) <- value' ", "length(assays(<",
             class(x), ">)) is 0")
    assays(x)[[1]] <- value
    x
})

setReplaceMethod("assay", c("SummarizedExperiment", "numeric"),
    function(x, i = 1, ..., value)
{
    assays(x, ...)[[i]] <- value
    x
})

setReplaceMethod("assay", c("SummarizedExperiment", "character"),
    function(x, i, ..., value)
{
    assays(x, ...)[[i]] <- value
    x
})

setGeneric("assayNames", function(x, ...) standardGeneric("assayNames"))

setMethod("assayNames", "SummarizedExperiment",
    function(x, ...)
{
    names(assays(x, withDimnames=FALSE))
})

setGeneric("assayNames<-",
    function(x, ..., value) standardGeneric("assayNames<-"))

setReplaceMethod("assayNames", c("SummarizedExperiment", "character"),
    function(x, ..., value)
{
    names(assays(x, withDimnames=FALSE)) <- value
    x
})

## cannonical location for dim, dimnames; dimnames should be checked
## for consistency (if non-null) and stripped from assays on
## construction, or added from assays if row/col names are NULL in
## <SummarizedExperiment> but not assays. dimnames need to be added on
## to assays when assays() invoked
setMethod("dim", "SummarizedExperiment",
    function(x)
{
    c(length(x), nrow(colData(x)))
})

setMethod("dimnames", "SummarizedExperiment",
    function(x)
{
    list(names(x), rownames(colData(x)))
})

setReplaceMethod("dimnames", c("SummarizedExperiment", "list"),
    function(x, value)
{
    NAMES <- S4Vectors:::normalize_names_replacement_value(value[[1]], x)
    colData <- colData(x)
    rownames(colData) <- value[[2]]
    BiocGenerics:::replaceSlots(x, NAMES=NAMES, colData=colData, check=FALSE)
})

setReplaceMethod("dimnames", c("SummarizedExperiment", "NULL"),
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
        msg <- paste(S4Vectors:::selectSome(orig[bad]), collapse=" ")
        stop(sprintf(fmt, msg))
    }
    idx
}

setMethod("[", c("SummarizedExperiment", "ANY", "ANY"),
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
        ii <- as.vector(i)
        ans_elementMetadata <- x@elementMetadata[i, , drop=FALSE]
        if (is(x, "RangedSummarizedExperiment")) {
            ans_rowRanges <- x@rowRanges[i]
        } else {
            ans_NAMES <- x@NAMES[ii]
        }
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
        ans_assays <- x@assays[ , jj]
        ans <- BiocGenerics:::replaceSlots(x, ...,
                       colData=ans_colData,
                       assays=ans_assays,
                       check=FALSE)
    } else if (missing(j)) {
        ans_assays <- x@assays[ii, ]
        if (is(x, "RangedSummarizedExperiment")) {
            ans <- BiocGenerics:::replaceSlots(x, ...,
                       elementMetadata=ans_elementMetadata,
                       rowRanges=ans_rowRanges,
                       assays=ans_assays,
                       check=FALSE)
        } else {
            ans <- BiocGenerics:::replaceSlots(x, ...,
                       elementMetadata=ans_elementMetadata,
                       NAMES=ans_NAMES,
                       assays=ans_assays,
                       check=FALSE)
        }
    } else {
        ans_assays <- x@assays[ii, jj]
        if (is(x, "RangedSummarizedExperiment")) {
            ans <- BiocGenerics:::replaceSlots(x, ...,
                       elementMetadata=ans_elementMetadata,
                       rowRanges=ans_rowRanges,
                       colData=ans_colData,
                       assays=ans_assays,
                       check=FALSE)
        } else {
            ans <- BiocGenerics:::replaceSlots(x, ...,
                       elementMetadata=ans_elementMetadata,
                       NAMES=ans_NAMES,
                       colData=ans_colData,
                       assays=ans_assays,
                       check=FALSE)
        }
    }
    ans
})

setReplaceMethod("[",
    c("SummarizedExperiment", "ANY", "ANY", "SummarizedExperiment"),
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
        if (is(x, "RangedSummarizedExperiment")) {
            ans_rowRanges <- local({
                r <- x@rowRanges
                r[i] <- value@rowRanges
                names(r)[ii] <- names(value@rowRanges)
                r
            })
        } else {
            ans_NAMES <- local({
                nms <- x@NAMES
                nms[ii] <- value@NAMES
                nms
            })
        }
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
        ans_assays <- local({
            a <- x@assays
            a[ , jj] <- value@assays
            a
        })
        ans <- BiocGenerics:::replaceSlots(x, ...,
                   metadata=ans_metadata,
                   colData=ans_colData,
                   assays=ans_assays,
                   check=FALSE)
        msg <- .valid.SummarizedExperiment.assays_ncol(ans)
    } else if (missing(j)) {
        ans_assays <- local({
            a <- x@assays
            a[ii, ] <- value@assays
            a
        })
        if (is(x, "RangedSummarizedExperiment")) {
            ans <- BiocGenerics:::replaceSlots(x, ...,
                       metadata=ans_metadata,
                       elementMetadata=ans_elementMetadata,
                       rowRanges=ans_rowRanges,
                       assays=ans_assays,
                       check=FALSE)
        } else {
            ans <- BiocGenerics:::replaceSlots(x, ...,
                       metadata=ans_metadata,
                       elementMetadata=ans_elementMetadata,
                       NAMES=ans_NAMES,
                       assays=ans_assays,
                       check=FALSE)
        }
        msg <- .valid.SummarizedExperiment.assays_nrow(ans)
    } else {
        ans_assays <- local({
            a <- x@assays
            a[ii, jj] <- value@assays
            a
        })
        if (is(x, "RangedSummarizedExperiment")) {
            ans <- BiocGenerics:::replaceSlots(x, ...,
                       metadata=ans_metadata,
                       elementMetadata=ans_elementMetadata,
                       rowRanges=ans_rowRanges,
                       colData=ans_colData,
                       assays=ans_assays,
                       check=FALSE)
        } else {
            ans <- BiocGenerics:::replaceSlots(x, ...,
                       metadata=ans_metadata,
                       elementMetadata=ans_elementMetadata,
                       NAMES=ans_NAMES,
                       colData=ans_colData,
                       assays=ans_assays,
                       check=FALSE)
        }
        msg <- .valid.SummarizedExperiment.assays_dim(ans)
    }
    if (!is.null(msg))
        stop(msg)
    ans
})

setMethod("extractROWS", "SummarizedExperiment",
    function(x, i)
    {
        i <- normalizeSingleBracketSubscript(i, x)
        x[i, ]
    }
)

setMethod("replaceROWS", "SummarizedExperiment",
    function(x, i, value)
    {
        i <- normalizeSingleBracketSubscript(i, x)
        x[i, ] <- value
        x
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Quick colData access.
###

setMethod("[[", c("SummarizedExperiment", "ANY", "missing"),
    function(x, i, j, ...)
{
    colData(x)[[i, ...]]
})

setReplaceMethod("[[", c("SummarizedExperiment", "ANY", "missing"),
    function(x, i, j, ..., value)
{
    colData(x)[[i, ...]] <- value
    x
})

.DollarNames.SummarizedExperiment <- function(x, pattern)
    grep(pattern, names(colData(x)), value=TRUE)

setMethod("$", "SummarizedExperiment",
    function(x, name)
{
    colData(x)[[name]]
})

setReplaceMethod("$", "SummarizedExperiment",
    function(x, name, value)
{
    colData(x)[[name]] <- value
    x
})


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Display.
###

setMethod("show", "SummarizedExperiment",
    function(object)
{
    selectSome <- S4Vectors:::selectSome
    scat <- function(fmt, vals=character(), exdent=2, ...)
    {
        vals <- ifelse(nzchar(vals), vals, "''")
        lbls <- paste(S4Vectors:::selectSome(vals), collapse=" ")
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
    nms <- assayNames(object)
    if (is.null(nms))
        nms <- character(length(assays(object, withDimnames=FALSE)))
    scat("assays(%d): %s\n", nms)

    ## rownames()
    dimnames <- dimnames(object)
    dlen <- sapply(dimnames, length)
    if (dlen[[1]]) scat("rownames(%d): %s\n", dimnames[[1]])
    else scat("rownames: NULL\n")

    ## rowData()
    scat("rowData names(%d): %s\n", names(rowData(object)))

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
setMethod("rbind", "SummarizedExperiment",
    function(..., deparse.level=1)
{
    args <- unname(list(...))
    .rbind.SummarizedExperiment(args)
})

.rbind.SummarizedExperiment <- function(args)
{
    if (!.compare(lapply(args, colnames)))
            stop("'...' objects must have the same colnames")
    if (!.compare(lapply(args, ncol)))
            stop("'...' objects must have the same number of samples")

    if (is(args[[1L]], "RangedSummarizedExperiment")) {
        rowRanges <- do.call(c, lapply(args, rowRanges))
    } else {
        ## Code below taken from combine_GAlignments_objects() from the
        ## GenomicAlignments package.

        ## Combine "NAMES" slots.
        NAMES_slots <- lapply(args, function(x) x@NAMES)
        ## TODO: Use elementIsNull() here when it becomes available.
        has_no_names <- sapply(NAMES_slots, is.null, USE.NAMES=FALSE)
        if (all(has_no_names)) {
            NAMES <- NULL
        } else {
            noname_idx <- which(has_no_names)
            if (length(noname_idx) != 0L)
                NAMES_slots[noname_idx] <-
                    lapply(elementNROWS(args[noname_idx]), character)
            NAMES <- unlist(NAMES_slots, use.names=FALSE)
        }
    }
    colData <- .cbind.DataFrame(args, colData, "colData")
    assays <- do.call(rbind, lapply(args, slot, "assays"))
    elementMetadata <- do.call(rbind, lapply(args, slot, "elementMetadata"))
    metadata <- do.call(c, lapply(args, metadata))

    if (is(args[[1L]], "RangedSummarizedExperiment")) {
        BiocGenerics:::replaceSlots(args[[1L]],
            rowRanges=rowRanges, colData=colData, assays=assays,
            elementMetadata=elementMetadata, metadata=metadata)
    } else {
        BiocGenerics:::replaceSlots(args[[1L]],
            NAMES=NAMES, colData=colData, assays=assays,
            elementMetadata=elementMetadata, metadata=metadata)
    }
}

### Appropriate for objects with same ranges and different samples.
setMethod("cbind", "SummarizedExperiment",
    function(..., deparse.level=1)
{
    args <- unname(list(...))
    .cbind.SummarizedExperiment(args)
})

.cbind.SummarizedExperiment <- function(args)
{
    if (is(args[[1L]], "RangedSummarizedExperiment")) {
        if (!.compare(lapply(args, rowRanges), TRUE))
            stop("'...' object ranges (rows) are not compatible")
        rowRanges <- rowRanges(args[[1L]])
        mcols(rowRanges) <- .cbind.DataFrame(args, mcols, "mcols")
    } else {
        elementMetadata <- .cbind.DataFrame(args, mcols, "mcols")
    }
    colData <- do.call(rbind, lapply(args, colData))
    assays <- do.call(cbind, lapply(args, slot, "assays"))
    metadata <- do.call(c, lapply(args, metadata))

    if (is(args[[1L]], "RangedSummarizedExperiment")) {
        BiocGenerics:::replaceSlots(args[[1L]],
            rowRanges=rowRanges,
            colData=colData, assays=assays, metadata=metadata)
    } else {
        BiocGenerics:::replaceSlots(args[[1L]],
            elementMetadata=elementMetadata,
            colData=colData, assays=assays, metadata=metadata)
    }
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
            if (!identicalVals(x1, x[[i]]))
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
        names(nmsv) <- rep(seq_along(nms), elementNROWS(nms))
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


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### identicalVals()
###
### Internal generic and methods (i.e. not exported).
### Provides a fast implementation of 'length(x) == length(y) && all(x == y)'
### for various kinds of vector-like objects.
### TODO: Move this to S4Vectors (for the generic and methods for factor and
### Rle objects) and IRanges (for the method for Ranges objects).
###

setGeneric("identicalVals", function(x, y) standardGeneric("identicalVals"))

### Semantically equivalent to identical(as.character(x), as.character(y))
### but avoids turning the 2 factor objects into character vectors so is more
### efficient.
setMethod("identicalVals", c("factor", "factor"),
    function(x, y)
    {
        m <- match(levels(y), levels(x), nomatch=0L)
        identical(as.integer(x), m[y])
    }
)

### Only support factor-Rle objects at the moment!
### Semantically equivalent to identical(as.character(x), as.character(y))
### but avoids turning the 2 factor-Rle objects into character vectors so is
### more efficient.
setMethod("identicalVals", c("Rle", "Rle"),
    function(x, y) identical(runLength(x), runLength(y)) &&
                   identicalVals(runValue(x), runValue(y))
)

setMethod("identicalVals", c("Ranges", "Ranges"),
    function(x, y) identical(start(x), start(y)) &&
                   identical(width(x), width(y))
)

### Like 'x == y' this method ignores circularity of the underlying sequences
### e.g. ranges [1, 10] and [101, 110] represent the same position on a
### circular sequence of length 100 so should be considered equal. However
### for 'x == y' and the method below, they are not.
### TODO: Take circularity of the underlying sequences into account.
setMethod("identicalVals", c("GenomicRanges", "GenomicRanges"),
    function(x, y)
    {
        ## Trying to merge 'seqinfo(x)' and 'seqinfo(y)' will raise an error
        ## if 'x' and 'y' are not based on the same reference genome. This is
        ## the standard way to check that 'x' and 'y' are based on the same
        ## reference genome.
        merge(seqinfo(x), seqinfo(y))  # we ignore the returned value

        identicalVals(seqnames(x), seqnames(y)) &&
            identicalVals(ranges(x), ranges(y)) &&
            identicalVals(strand(x), strand(y))
    }
)

