### =========================================================================
### RangedSummarizedExperiment objects
### -------------------------------------------------------------------------
###


### The 'elementMetadata' slot must contain a zero-column DataFrame at all time
### (this is checked by the validity method). The top-level mcols are stored on
### the rowRanges component.
setClass("RangedSummarizedExperiment",
    contains="SummarizedExperiment",
    representation(
        rowRanges="GenomicRanges_OR_GRangesList"
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

### The names and mcols of a RangedSummarizedExperiment must be set on its
### rowRanges slot, not in its NAMES and elementMetadata slots!
.valid.RangedSummarizedExperiment <- function(x)
{
    if (!is.null(x@NAMES))
        return("'NAMES' slot must be set to NULL at all time")
    if (ncol(x@elementMetadata) != 0L)
        return(wmsg("'elementMetadata' slot must contain a zero-column ",
                    "DataFrame at all time"))
    rowRanges_len <- length(x@rowRanges)
    x_nrow <- length(x)
    if (rowRanges_len != x_nrow) {
        txt <- sprintf(
            "\n  length of 'rowRanges' (%d) must equal nb of rows in 'x' (%d)",
            rowRanges_len, x_nrow)
        return(txt)
    }
    NULL
}

setValidity2("RangedSummarizedExperiment", .valid.RangedSummarizedExperiment)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Constructor.
###

.new_RangedSummarizedExperiment <- function(assays, rowRanges, colData,
                                            metadata)
{
    elementMetadata <- S4Vectors:::make_zero_col_DataFrame(length(rowRanges))
    if (!is(assays, "Assays"))
        assays <- Assays(assays)
    new("RangedSummarizedExperiment", rowRanges=rowRanges,
                                      colData=colData,
                                      assays=assays,
                                      elementMetadata=elementMetadata,
                                      metadata=as.list(metadata))
}

.get_colnames_from_assays <- function(assays)
{
    if (length(assays) == 0L)
        return(NULL)
    colnames(assays[[1L]])
}

.get_rownames_from_assays <- function(assays)
{
    if (length(assays) == 0L)
        return(NULL)
    rownames(assays[[1L]])
}

setGeneric("SummarizedExperiment",
    function(assays, ...) standardGeneric("SummarizedExperiment"))

setMethod("SummarizedExperiment", "SimpleList",
   function(assays, rowData=NULL, rowRanges=GRangesList(), colData=DataFrame(),
            metadata=list())
{
    if (missing(colData) && 0L != length(assays)) {
        assay <- assays[[1]]
        nms <- colnames(assay)
        colData <- DataFrame(x=seq_len(ncol(assay)), row.names=nms)[, FALSE]
    } else if (!missing(colData)) {
        if (!is(colData, "DataFrame"))
            colData <- as(colData, "DataFrame")
        if (is.null(rownames(colData)))
            rownames(colData) <- .get_colnames_from_assays(assays)
    }
    ans_colnames <- rownames(colData)

    if (is.null(rowData)) {
        if (missing(rowRanges)) {
            ans_rownames <- .get_rownames_from_assays(assays)
        } else {
            if (is.null(names(rowRanges)))
                names(rowRanges) <- .get_rownames_from_assays(assays)
            ans_rownames <- names(rowRanges)
        }
    } else {
        if (!missing(rowRanges))
            stop("only one of 'rowData' and 'rowRanges' can be specified")
        if (is(rowData, "GenomicRanges_OR_GRangesList")) {
            rowRanges <- rowData
            if (is.null(names(rowRanges)))
                names(rowRanges) <- .get_rownames_from_assays(assays)
            ans_rownames <- names(rowRanges)
        } else {
            if (!is(rowData, "DataFrame"))
                rowData <- as(rowData, "DataFrame")
            ans_rownames <- rownames(rowData)
            if (is.null(ans_rownames))
                ans_rownames <- .get_rownames_from_assays(assays)
        }
    }

    ## validate the assay rownames and colnames
    .validate_names <- function(nms, expected_nms, what1, what2)
    {
        if (is.null(nms))
            return()
        if (!is.character(nms))
            stop(wmsg(what1, " must be NULL or a character vector"))
        if (!is.null(attributes(nms)))
            stop(wmsg(what1, " must be NULL or a character vector",
                             " with no attributes"))
        if (!identical(nms, expected_nms))
            stop(wmsg(what1, " must be NULL or identical to ", what2))
    }

    for(i in seq_along(assays)) {
        a <- assays[[i]]
        rownames <- rownames(a)
        .validate_names(rownames, ans_rownames,
                        "assay rownames()",
                        "rowData rownames() / rowRanges names()")
        colnames <- colnames(a)
        .validate_names(colnames, ans_colnames,
                        "assay colnames()",
                        "colData rownames()")
    }

    assays <- Assays(assays)

    if (missing(rowRanges) && !is(rowData, "GenomicRanges_OR_GRangesList")) {
        new_SummarizedExperiment(assays, ans_rownames, rowData, colData,
                                 metadata)
    } else {
        .new_RangedSummarizedExperiment(assays, rowRanges, colData, metadata)
    }
})

setMethod("SummarizedExperiment", "ANY",
    function(assays, ...)
{
    if (is.matrix(assays) && is.list(assays))
        ## special case -- matrix of lists
        assays <- list(assays)
    SummarizedExperiment(SimpleList(assays), ...)
})

setMethod("SummarizedExperiment", "list",
    function(assays, ...)
{
    SummarizedExperiment(do.call(SimpleList, assays), ...)
})

setMethod("SummarizedExperiment", "missing",
    function(assays, ...)
{
    SummarizedExperiment(SimpleList(), ...)
})


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Coercion.
###
### See makeSummarizedExperimentFromExpressionSet.R for coercion back and
### forth between SummarizedExperiment and ExpressionSet.
###

.from_RangedSummarizedExperiment_to_SummarizedExperiment <- function(from)
{
    new_SummarizedExperiment(from@assays,
                             names(from@rowRanges),
                             mcols(from@rowRanges, use.names=FALSE),
                             from@colData,
                             from@metadata)
}

setAs("RangedSummarizedExperiment", "SummarizedExperiment",
    .from_RangedSummarizedExperiment_to_SummarizedExperiment
)

.from_SummarizedExperiment_to_RangedSummarizedExperiment <- function(from)
{
    partitioning <- PartitioningByEnd(integer(length(from)), names=names(from))
    rowRanges <- relist(GRanges(), partitioning)
    mcols(rowRanges) <- mcols(from, use.names=FALSE)
    .new_RangedSummarizedExperiment(from@assays,
                                    rowRanges,
                                    from@colData,
                                    from@metadata)
}

setAs("SummarizedExperiment", "RangedSummarizedExperiment",
    .from_SummarizedExperiment_to_RangedSummarizedExperiment
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Getters and setters.
###

### The rowRanges() generic is defined in the DelayedArray package.
setMethod("rowRanges", "SummarizedExperiment",
    function(x, ...) NULL
)

### Fix old GRanges instances on-the-fly.
setMethod("rowRanges", "RangedSummarizedExperiment",
    function(x, ...) updateObject(x@rowRanges, check=FALSE)
)

setGeneric("rowRanges<-",
    function(x, ..., value) standardGeneric("rowRanges<-"))

### No-op.
setReplaceMethod("rowRanges", c("SummarizedExperiment", "NULL"),
    function(x, ..., value) x
)

### Degrade 'x' to SummarizedExperiment instance.
setReplaceMethod("rowRanges", c("RangedSummarizedExperiment", "NULL"),
    function(x, ..., value) as(x, "SummarizedExperiment", strict=TRUE)
)

.SummarizedExperiment.rowRanges.replace <-
    function(x, ..., value)
{
    if (!is(x, "RangedSummarizedExperiment"))
        x <- as(x, "RangedSummarizedExperiment")
    x <- BiocGenerics:::replaceSlots(x, ...,
             rowRanges=value,
             elementMetadata=S4Vectors:::make_zero_col_DataFrame(length(value)),
             check=FALSE)
    msg <- .valid.SummarizedExperiment.assays_nrow(x)
    if (!is.null(msg))
        stop(msg)
    x
}

setReplaceMethod("rowRanges", c("SummarizedExperiment", "GenomicRanges"),
    .SummarizedExperiment.rowRanges.replace)

setReplaceMethod("rowRanges", c("SummarizedExperiment", "GRangesList"),
    .SummarizedExperiment.rowRanges.replace)

setMethod("names", "RangedSummarizedExperiment",
    function(x) names(rowRanges(x))
)

setReplaceMethod("names", "RangedSummarizedExperiment",
    function(x, value)
{
    rowRanges <- rowRanges(x)
    names(rowRanges) <- value
    BiocGenerics:::replaceSlots(x, rowRanges=rowRanges, check=FALSE)
})

setMethod("dimnames", "RangedSummarizedExperiment",
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
    BiocGenerics:::replaceSlots(x,
        rowRanges=rowRanges,
        colData=colData,
        check=FALSE)
})


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Subsetting.
###

.DollarNames.RangedSummarizedExperiment <-
    .DollarNames.SummarizedExperiment

setMethod("subset", "RangedSummarizedExperiment",
    function(x, subset, select, ...)
{
    i <- S4Vectors:::evalqForSubset(subset, rowRanges(x), ...)
    j <- S4Vectors:::evalqForSubset(select, colData(x), ...)
    x[i, j]
})


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
## colData-as-GRanges compatibility: allow direct access to GRanges /
## GRangesList colData for select functions

## Not supported:
## 
## Not consistent SummarizedExperiment structure: length, names,
##   as.data.frame, c.
## Length-changing endomorphisms: disjoin, gaps, reduce, unique.
## 'legacy' data types / functions: as "RangedData", as "IntegerRangesList",
##   renameSeqlevels, keepSeqlevels.
## Possile to implement, but not yet: Ops, map, window, window<-

## mcols
setMethod("mcols", "RangedSummarizedExperiment",
    function(x, use.names=TRUE, ...)
{
    mcols(rowRanges(x), use.names=use.names, ...)
})

setReplaceMethod("mcols", "RangedSummarizedExperiment",
    function(x, ..., value)
{
    BiocGenerics:::replaceSlots(x,
        rowRanges=local({
            r <- rowRanges(x)
            mcols(r) <- value
            r
        }),
        check=FALSE)
})

### mcols() is the recommended way for accessing the metadata columns.
### Use of values() or elementMetadata() is discouraged.

setMethod("elementMetadata", "RangedSummarizedExperiment",
    function(x, use.names=FALSE, ...)
{
    elementMetadata(rowRanges(x), use.names=use.names, ...)
})

setReplaceMethod("elementMetadata", "RangedSummarizedExperiment",
    function(x, ..., value)
{
    elementMetadata(rowRanges(x), ...) <- value
    x
})

## Single dispatch, generic signature fun(x, ...)
local({
    .funs <-
        c("duplicated", "end", "end<-", "ranges", "seqinfo", "seqnames",
          "start", "start<-", "strand", "width", "width<-")

    endomorphisms <- .funs[grepl("<-$", .funs)]

    tmpl <- function() {}
    environment(tmpl) <- parent.frame(2)
    for (.fun in .funs) {
        generic <- getGeneric(.fun)
        formals(tmpl) <- formals(generic)
        fmls <- as.list(formals(tmpl))
        fmls[] <- sapply(names(fmls), as.symbol)
        fmls[[generic@signature]] <- quote(rowRanges(x))
        if (.fun %in% endomorphisms)
            body(tmpl) <- substitute({
                rowRanges(x) <- do.call(FUN, ARGS)
                x
            }, list(FUN=.fun, ARGS=fmls))
        else
            body(tmpl) <-
                substitute(do.call(FUN, ARGS),
                           list(FUN=as.symbol(.fun), ARGS=fmls))
        setMethod(.fun, "RangedSummarizedExperiment", tmpl)
    }
})

setMethod("granges", "RangedSummarizedExperiment",
    function(x, use.mcols=FALSE, ...)
{
    if (!identical(use.mcols, FALSE))
        stop("\"granges\" method for RangedSummarizedExperiment objects ",
             "does not support the 'use.mcols' argument")
    rowRanges(x)
})

## 2-argument dispatch:
## pcompare / Compare 
## 
.RangedSummarizedExperiment.pcompare <-
    function(x, y)
{
    if (is(x, "RangedSummarizedExperiment"))
        x <- rowRanges(x)
    if (is(y, "RangedSummarizedExperiment"))
        y <- rowRanges(y)
    pcompare(x, y)
}

.RangedSummarizedExperiment.Compare <-
    function(e1, e2)
{
    if (is(e1, "RangedSummarizedExperiment"))
        e1 <- rowRanges(e1)
    if (is(e2, "RangedSummarizedExperiment"))
        e2 <- rowRanges(e2)
    callGeneric(e1=e1, e2=e2)
}

local({
    .signatures <- list(
        c("RangedSummarizedExperiment", "ANY"),
        c("ANY", "RangedSummarizedExperiment"),
        c("RangedSummarizedExperiment", "RangedSummarizedExperiment"))

    for (.sig in .signatures) {
        setMethod("pcompare", .sig, .RangedSummarizedExperiment.pcompare)
        setMethod("Compare", .sig, .RangedSummarizedExperiment.Compare)
    }
})

## additional getters / setters

setReplaceMethod("strand", "RangedSummarizedExperiment",
    function(x, ..., value)
{
    strand(rowRanges(x)) <- value
    x
})

setReplaceMethod("ranges", "RangedSummarizedExperiment",
    function(x, ..., value)
{
    ranges(rowRanges(x)) <- value
    x
})

## order, rank, sort

setMethod("is.unsorted", "RangedSummarizedExperiment",
    function(x, na.rm = FALSE, strictly = FALSE, ignore.strand = FALSE)
{
    x <- rowRanges(x)
    if (!is(x, "GenomicRanges"))
        stop("is.unsorted() is not yet supported when 'rowRanges(x)' is a ",
             class(x), " object")
    callGeneric()
})

setMethod("order", "RangedSummarizedExperiment",
    function(..., na.last=TRUE, decreasing=FALSE,
             method=c("auto", "shell", "radix"))
{
    args <- lapply(list(...), rowRanges)
    do.call("order", c(args, list(na.last=na.last,
                                  decreasing=decreasing,
                                  method=method)))
})

setMethod("rank", "RangedSummarizedExperiment",
    function(x, na.last = TRUE,
        ties.method = c("average", "first", "last", "random", "max", "min"))
{
    ties.method <- match.arg(ties.method)
    rank(rowRanges(x), na.last=na.last, ties.method=ties.method)
})

setMethod("sort", "RangedSummarizedExperiment",
    function(x, decreasing = FALSE, ignore.strand = FALSE)
{
    x_rowRanges <- rowRanges(x)
    if (!is(x_rowRanges, "GenomicRanges"))
        stop("sort() is not yet supported when 'rowRanges(x)' is a ",
             class(x_rowRanges), " object")
    oo <- GenomicRanges:::order_GenomicRanges(x_rowRanges,
                                              decreasing = decreasing,
                                              ignore.strand = ignore.strand)
    x[oo]
})

## seqinfo (also seqlevels, genome, seqlevels<-, genome<-), seqinfo<-

setMethod("seqinfo", "RangedSummarizedExperiment",
    function(x)
{
    seqinfo(x@rowRanges)
})

setReplaceMethod("seqinfo", "RangedSummarizedExperiment",
    function (x, new2old= NULL,
              pruning.mode=c("error", "coarse", "fine", "tidy"),
              value)
{
    if (!is(value, "Seqinfo")) 
        stop("the supplied 'seqinfo' must be a Seqinfo object")
    dangling_seqlevels <-
        GenomeInfoDb:::getDanglingSeqlevels(x@rowRanges, new2old=new2old,
                                            pruning.mode=pruning.mode,
                                            seqlevels(value))
    if (length(dangling_seqlevels) != 0L) 
        x <- x[!(seqnames(x) %in% dangling_seqlevels)]
    x@rowRanges <-
        update(x@rowRanges,
               seqnames = GenomeInfoDb:::makeNewSeqnames(x, new2old,
                                                         seqlevels(value)),
               seqinfo = value)
    if (is.character(msg <- .valid.RangedSummarizedExperiment(x)))
        stop(msg)
    x
})

setMethod("split", "RangedSummarizedExperiment",
    function(x, f, drop=FALSE, ...) 
{
    splitAsList(x, f, drop=drop)
})


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### updateObject()
###

.updateObject_RangedSummarizedExperiment <- function(object, ..., verbose=FALSE)
{
    object <- callNextMethod()  # call method for SummarizedExperiment objects
    object@rowRanges <- updateObject(object@rowRanges, ..., verbose=verbose)
    object
}

setMethod("updateObject", "RangedSummarizedExperiment",
    .updateObject_RangedSummarizedExperiment
)

