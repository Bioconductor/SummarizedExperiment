### =========================================================================
### SummarizedExperiment objects
### -------------------------------------------------------------------------
###


setClassUnion("Assays_OR_NULL", c("Assays", "NULL"))

setClass("SummarizedExperiment",
    contains=c("RectangularData", "Vector"),
    representation(
        colData="DataFrame",            # columns and their annotations
        assays="Assays_OR_NULL",        # Data -- e.g., list of matrices
        NAMES="character_OR_NULL",
        elementMetadata="DataFrame"
    ),
    prototype(
        colData=new("DFrame"),
        elementMetadata=new("DFrame")
    )
)

### Combine the new "parallel slots" with those of the parent class. Make
### sure to put the new parallel slots **first**. See R/Vector-class.R file
### in the S4Vectors package for what slots should or should not be considered
### "parallel".
setMethod("parallel_slot_names", "SummarizedExperiment",
    function(x) c("assays", "NAMES", callNextMethod())
)

setMethod("vertical_slot_names", "SummarizedExperiment",
    function(x) parallel_slot_names(x)
)

### Like parallel_slot_names() methods, horizontal_slot_names() methods for
### SummarizedExperiment derivatives should be defined in an incremental
### fashion, that is, they should only explicitly list the new "horizontal
### slots" (i.e. the horizontal slots that they add to their parent class).
### See R/RectangularData-class.R file in the S4Vectors package for what
### slots should or should not be considered "horizontal".
setMethod("horizontal_slot_names", "SummarizedExperiment",
    function(x) "colData"
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Validity
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
### Low-level constructor (not exported)
###

new_SummarizedExperiment <- function(assays, names, rowData, colData,
                                     metadata)
{
    assays <- Assays(assays, as.null.if.no.assay=TRUE)
    if (is.null(rowData)) {
        if (!is.null(names)) {
            nrow <- length(names)
        } else if (!is.null(assays)) {
            nrow <- nrow(assays)
        } else {
            nrow <- 0L
        }
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
### Accessors
###

setMethod("length", "SummarizedExperiment",
    function(x) nrow(x@elementMetadata)
)

setMethod("names", "SummarizedExperiment", function(x) x@NAMES)

setReplaceMethod("names", "SummarizedExperiment",
    function(x, value)
    {
        NAMES <- S4Vectors:::normarg_names(value, class(x), length(x))
        BiocGenerics:::replaceSlots(x, NAMES=NAMES, check=FALSE)
    }
)

## rowData, colData seem too vague, but from eSet derived classes wanted to
## call the rows / cols something different from 'features' or 'samples', so
## might as well avoid the issue

setGeneric("rowData", signature="x",
    function(x, use.names=TRUE, ...) standardGeneric("rowData")
)

### Fix old DataFrame instances on-the-fly (mcols() does it).
setMethod("rowData", "SummarizedExperiment",
    function(x, use.names=TRUE, ...) mcols(x, use.names=use.names, ...)
)

setGeneric("rowData<-",
    function(x, ..., value) standardGeneric("rowData<-"))

setReplaceMethod("rowData", "SummarizedExperiment",
    function(x, ..., value) `mcols<-`(x, ..., value=value)
)

setGeneric("colData", function(x, ...) standardGeneric("colData"))

### Fix old DataFrame instances on-the-fly.
setMethod("colData", "SummarizedExperiment",
    function(x, ...) updateObject(x@colData, check=FALSE)
)

setGeneric("colData<-",
    function(x, ..., value) standardGeneric("colData<-"))

setReplaceMethod("colData", c("SummarizedExperiment", "DataFrame"),
    function(x, ..., value)
{
    if (nrow(value) != ncol(x))
        stop("nrow of supplied 'colData' must equal ncol of object")
    x <- updateObject(x, check=FALSE)
    BiocGenerics:::replaceSlots(x, colData=value, check=FALSE)
})

setReplaceMethod("colData", c("SummarizedExperiment", "NULL"),
    function(x, ..., value)
{
    value <- new2("DFrame", nrows=ncol(x), rownames=colnames(x), check=FALSE)
    x <- updateObject(x, check=FALSE)
    BiocGenerics:::replaceSlots(x, colData=value, check=FALSE)
})

setGeneric("assays", signature="x",
    function(x, withDimnames=TRUE, ...) standardGeneric("assays")
)

setMethod("assays", "SummarizedExperiment",
    function(x, withDimnames=TRUE, ...)
{
    if (!isTRUEorFALSE(withDimnames))
        stop(wmsg("'withDimnames' must be TRUE or FALSE"))
    assays <- as(x@assays, "SimpleList")
    if (withDimnames) {
        x_dimnames <- dimnames(x)
        if (is.null(x_dimnames))
            x_dimnames <- list(NULL, NULL)
        assays <- endoapply(assays,
            function(a) {
                a_dimnames <- dimnames(a)
                a_dimnames[1:2] <- x_dimnames
                a_dimnames <- DelayedArray:::simplify_NULL_dimnames(a_dimnames)
                DelayedArray:::set_dimnames(a, a_dimnames)
            }
        )
    }
    assays
})

setGeneric("assays<-", signature=c("x", "value"),
    function(x, withDimnames=TRUE, ..., value) standardGeneric("assays<-"),
)

### For assays with more than 2 dimensions, **only** the first 2
### dimnames components (i.e. rownames and colnames) are compared
### with 'expected_dimnames'.
### If 'strict' is FALSE, then a NULL on either side of the comparison
### is enough for the comparison to be considered succesful.
assays_have_expected_dimnames <- function(assays, expected_dimnames,
                                          strict=TRUE)
{
    if (length(assays) == 0L)
        return(TRUE)
    ## Returns the rownames and colnames in a list of length 2. Also list
    ## elements in the returned list that are equal to character(0) are
    ## replaced with NULLs. This is to accommodate the fact that, unlike
    ## an ordinary vector or data.frame, an ordinary matrix or array object
    ## will never hold a character(0) along a dimension of zero extend.
    get_rownames_and_colnames <- function(dn) {
        if (is.null(dn))
            return(vector("list", length=2L))
        ans <- dn[1:2]
        if (length(ans[[1L]]) == 0L)
            ans[1L] <- list(NULL)
        if (length(ans[[2L]]) == 0L)
            ans[2L] <- list(NULL)
        ans
    }
    expected_dimnames <- get_rownames_and_colnames(expected_dimnames)
    if (!strict) {
        expected_rownames_is_NULL <- is.null(expected_dimnames[[1L]])
        expected_colnames_is_NULL <- is.null(expected_dimnames[[2L]])
        if (expected_rownames_is_NULL && expected_colnames_is_NULL)
            return(TRUE)
    }
    ok <- vapply(seq_along(assays),
        function(i) {
            a <- getListElement(assays, i)
            a_dimnames <- get_rownames_and_colnames(dimnames(a))
            a_rownames <- a_dimnames[[1L]]
            a_colnames <- a_dimnames[[2L]]
            if (strict) {
                ok1 <- identical(a_rownames, expected_dimnames[[1L]])
                ok2 <- identical(a_colnames, expected_dimnames[[2L]])
            } else {
                ## Comparing the rownames/colnames with 'identical(x, y)'
                ## would cause something like
                ##   m <- matrix(1:12, ncol=3,
                ##               dimnames=list(NULL, c(A="a", B="b", C="c")))
                ##   SummarizedExperiment(m)
                ## to fail (because of the names on 'colnames(m)'). So we
                ## use less stringent 'length(x) == length(y) && all(x == y)'
                ## instead.
                ok1 <- is.null(a_rownames) ||
                       expected_rownames_is_NULL ||
                       (length(a_rownames) == length(expected_dimnames[[1L]])
                        && all(a_rownames == expected_dimnames[[1L]]))
                ok2 <- is.null(a_colnames) ||
                       expected_colnames_is_NULL ||
                       (length(a_colnames) == length(expected_dimnames[[2L]])
                        && all(a_colnames == expected_dimnames[[2L]]))
            }
            ok1 && ok2
        },
        logical(1)
    )
    all(ok)
}

.SummarizedExperiment.assays.replace <-
    function(x, withDimnames=TRUE, ..., value)
{
    if (!isTRUEorFALSE(withDimnames))
        stop(wmsg("'withDimnames' must be TRUE or FALSE"))
    value <- normarg_assays(value, as.null.if.no.assay=TRUE)
    ## By default the dimnames on the supplied assays must be identical to
    ## the dimnames on 'x'. The user must use 'withDimnames=FALSE' if it's
    ## not the case. This is for symetry with the behavior of the getter.
    ## See https://github.com/Bioconductor/SummarizedExperiment/issues/35
    if (withDimnames && !assays_have_expected_dimnames(value, dimnames(x)))
        stop(wmsg("please use 'assay(x, withDimnames=FALSE)) <- value' ",
                  "or 'assays(x, withDimnames=FALSE)) <- value' when ",
                  "the rownames or colnames of the supplied assay(s) ",
                  "are not identical to those of the receiving ",
                  class(x), " object 'x'"))
    new_assays <- Assays(value, as.null.if.no.assay=TRUE)
    x <- BiocGenerics:::replaceSlots(x, assays=new_assays, check=FALSE)
    ## validObject(x) should NOT be called below because it would then
    ## fully re-validate objects that derive from SummarizedExperiment
    ## (e.g. DESeqDataSet objects) after the user sets the assays slot with
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

setGeneric("assay", signature=c("x", "i"),
    function(x, i, withDimnames=TRUE, ...) standardGeneric("assay")
)

## convenience for common use case
setMethod("assay", c("SummarizedExperiment", "missing"),
    function(x, i, withDimnames=TRUE, ...)
{
    assays <- assays(x, withDimnames=withDimnames, ...)
    if (0L == length(assays))
        stop("'assay(<", class(x), ">, i=\"missing\", ...) ",
             "length(assays(<", class(x), ">)) is 0'")
    assays[[1]]
})

setMethod("assay", c("SummarizedExperiment", "numeric"),
    function(x, i, withDimnames=TRUE, ...)
{
    tryCatch({
        assays(x, withDimnames=withDimnames, ...)[[i]]
    }, error=function(err) {
        stop("'assay(<", class(x), ">, i=\"numeric\", ...)' ",
             "invalid subscript 'i'\n", conditionMessage(err))
    })
})

setMethod("assay", c("SummarizedExperiment", "character"),
    function(x, i, withDimnames=TRUE, ...)
{
    msg <- paste0("'assay(<", class(x), ">, i=\"character\", ...)' ",
                  "invalid subscript 'i'")
    res <- tryCatch({
        assays(x, withDimnames=withDimnames, ...)[[i]]
    }, error=function(err) {
        stop(msg, "\n", conditionMessage(err))
    })
    if (is.null(res))
        stop(msg, "\n'", i, "' not in names(assays(<", class(x), ">))")
    res
})

setGeneric("assay<-", signature=c("x", "i"),
    function(x, i, withDimnames=TRUE, ..., value) standardGeneric("assay<-"))

setReplaceMethod("assay", c("SummarizedExperiment", "missing"),
    function(x, i, withDimnames=TRUE, ..., value)
{
    if (0L == length(assays(x, withDimnames=FALSE)))
        stop("'assay(<", class(x), ">) <- value' ", "length(assays(<",
             class(x), ">)) is 0")
    assays(x, withDimnames=withDimnames, ...)[[1]] <- value
    x
})

setReplaceMethod("assay", c("SummarizedExperiment", "numeric"),
    function(x, i, withDimnames=TRUE, ..., value)
{
    assays(x, withDimnames=withDimnames, ...)[[i]] <- value
    x
})

setReplaceMethod("assay", c("SummarizedExperiment", "character"),
    function(x, i, withDimnames=TRUE, ..., value)
{
    assays(x, withDimnames=withDimnames, ...)[[i]] <- value
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

setMethod("nrow", "SummarizedExperiment", function(x) length(x))
setMethod("ncol", "SummarizedExperiment", function(x) nrow(colData(x)))

setMethod("rownames", "SummarizedExperiment", function(x) names(x))
setMethod("colnames", "SummarizedExperiment", function(x) rownames(colData(x)))

setReplaceMethod("dimnames", c("SummarizedExperiment", "list"),
    function(x, value)
{
    NAMES <- S4Vectors:::normarg_names(value[[1L]], class(x), length(x))
    colData <- colData(x)
    rownames(colData) <- value[[2L]]
    BiocGenerics:::replaceSlots(x, NAMES=NAMES, colData=colData, check=FALSE)
})

setReplaceMethod("dimnames", c("SummarizedExperiment", "NULL"),
    function(x, value)
{
    dimnames(x) <- list(NULL, NULL)
    x
})


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Subsetting
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

### TODO: Refactor this to use extractROWS().
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

### TODO: Refactor this to use replaceROWS().
setReplaceMethod("[",
    c("SummarizedExperiment", "ANY", "ANY", "SummarizedExperiment"),
    function(x, i, j, ..., value)
{
    if (missing(i) && missing(j))
        return(value)

    ans_metadata <- metadata(x)

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

setMethod("subset", "SummarizedExperiment",
    function(x, subset, select, ...)
{
    i <- S4Vectors:::evalqForSubset(subset, rowData(x, use.names=FALSE), ...)
    j <- S4Vectors:::evalqForSubset(select, colData(x), ...)
    x[i, j]
})


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Quick colData access
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

.DollarNames.SummarizedExperiment <- function(x, pattern = "")
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
### Show
###

setMethod("show", "SummarizedExperiment",
    function(object)
{
    cat("class:", class(object), "\n")
    cat("dim:", dim(object), "\n")

    ## metadata()
    expt <- names(metadata(object))
    if (is.null(expt))
        expt <- character(length(metadata(object)))
    coolcat("metadata(%d): %s\n", expt)

    ## assays()
    nms <- assayNames(object)
    if (is.null(nms))
        nms <- character(length(assays(object, withDimnames=FALSE)))
    coolcat("assays(%d): %s\n", nms)

    ## rownames()
    rownames <- rownames(object)
    if (!is.null(rownames)) coolcat("rownames(%d): %s\n", rownames)
    else cat("rownames: NULL\n")

    ## rowData`()
    coolcat("rowData names(%d): %s\n", names(rowData(object, use.names=FALSE)))

    ## colnames()
    colnames <- colnames(object)
    if (!is.null(colnames)) coolcat("colnames(%d): %s\n", colnames)
    else cat("colnames: NULL\n")

    ## colData()
    coolcat("colData names(%d): %s\n", names(colData(object)))
})

setMethod("showAsCell", "SummarizedExperiment",
    function(object) rep.int("####", NROW(object))
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Combine
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
    lst <- lapply(args, accessor, use.names=FALSE)
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
### Rle objects) and IRanges (for the method for IntegerRanges objects).
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

setMethod("identicalVals", c("IntegerRanges", "IntegerRanges"),
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


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### On-disk realization.
###

setMethod("realize", "SummarizedExperiment",
    function(x, BACKEND=getAutoRealizationBackend())
    {
        for (i in seq_along(assays(x))) {
            ## We drop the dimnames of the individual assays for 2 reasons:
            ##   1) These dimnames are kind of irrelevant. The dimnames that
            ##      really matter are 'dimnames(x)' and they are stored
            ##      somewhere else in 'x'. So we don't loose them by not
            ##      realizing the assay dimnames on disk. As a little extra
            ##      bonus, this actually saves a little bit of time and disk
            ##      space.
            ##   2) Using the HDF5Array backend to realize an array-like object
            ##      on disk doesn't store the dimnames in the HDF5 file at the
            ##      moment.
            a <- assay(x, i, withDimnames=FALSE)
            a <- DelayedArray:::set_dimnames(a, NULL)
            assay(x, i, withDimnames=FALSE) <- realize(a, BACKEND=BACKEND)
        }
        x
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### updateObject()
###

.updateObject_SummarizedExperiment <- function(object, ..., verbose=FALSE)
{
    object@assays <- updateObject(object@assays, ..., verbose=verbose)
    object@colData <- updateObject(object@colData, ..., verbose=verbose)
    callNextMethod()  # call method for Vector objects
}

setMethod("updateObject", "SummarizedExperiment",
    .updateObject_SummarizedExperiment
)

