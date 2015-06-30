### =========================================================================
### RangedSummarizedExperiment objects
### -------------------------------------------------------------------------
###


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

### The names and mcols of a RangedSummarizedExperiment must be set on its
### rowRanges slot, not in its NAMES and elementMetadata slots!
.valid.RangedSummarizedExperiment <- function(x)
{
    if (!is.null(x@NAMES))
        return("'NAMES' slot must be set to NULL at all time")
    if (ncol(x@elementMetadata) != 0L)
        return(wmsg("'elementMetadata' slot must contain a zero-column ",
                    "DataFrame at all time"))
    NULL
}

setValidity2("RangedSummarizedExperiment", .valid.RangedSummarizedExperiment)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Constructor.
###

.new_RangedSummarizedExperiment <- function(assays, rowRanges, colData,
                                            metadata)
{
    elementMetadata <- new("DataFrame", nrows=length(rowRanges))
    if (!is(assays, "Assays"))
        assays <- Assays(assays)
    new("RangedSummarizedExperiment", rowRanges=rowRanges,
                                      colData=colData,
                                      assays=assays,
                                      elementMetadata=elementMetadata,
                                      metadata=as.list(metadata))
}

get_rownames_from_assays <- function(assays)
{
    if (length(assays) == 0L)
        return(NULL)
    rownames(assays[[1L]])
}

setMethod(SummarizedExperiment, "SimpleList",
   function(assays, rowRanges=GRangesList(), colData=DataFrame(),
            metadata=list(), exptData=SimpleList())
{
    if (missing(colData) && 0L != length(assays)) {
        nms <- colnames(assays[[1]])
        if (is.null(nms) && 0L != ncol(assays[[1]]))
            stop("'SummarizedExperiment' assay colnames must not be NULL")
        colData <- DataFrame(row.names=nms)
    }

    if (missing(rowRanges)) {
        ans_rownames <- get_rownames_from_assays(assays)
    } else {
        ans_rownames <- names(rowRanges)
    }
    ans_colnames <- rownames(colData)
    ans_dimnames <- list(ans_rownames, ans_colnames)
    FUN <- function(x) {
        ## dimnames(x) as NULL or list(NULL, NULL)
        all(sapply(dimnames(x), is.null)) ||
            ## or consistent with 'ans_dimnames'
            identical(dimnames(x)[1:2], ans_dimnames)
    }
    if (!all(sapply(assays, FUN)))
        assays <- endoapply(assays, unname)
    assays <- Assays(assays)

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

    if (missing(rowRanges)) {
        new_SummarizedExperiment0(assays, ans_rownames, NULL, colData, metadata)
    } else {
        .new_RangedSummarizedExperiment(assays, rowRanges, colData, metadata)
    }
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
    new_SummarizedExperiment0(from@assays,
                              names(from@rowRanges),
                              mcols(from@rowRanges),
                              from@colData,
                              from@metadata)
}

setAs("RangedSummarizedExperiment", "SummarizedExperiment0",
    .from_RangedSummarizedExperiment_to_SummarizedExperiment0
)

.from_SummarizedExperiment0_to_RangedSummarizedExperiment <- function(from)
{
    partitioning <- PartitioningByEnd(integer(length(from)), names=names(from))
    rowRanges <- relist(GRanges(), partitioning)
    .new_RangedSummarizedExperiment(from@assays,
                                    rowRanges,
                                    from@colData,
                                    from@metadata)
}

setAs("SummarizedExperiment0", "RangedSummarizedExperiment",
    .from_SummarizedExperiment0_to_RangedSummarizedExperiment
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Getters and setters.
###

setMethod(rowRanges, "RangedSummarizedExperiment",
    function(x, ...) x@rowRanges
)

.SummarizedExperiment0.rowRanges.replace <-
    function(x, ..., value)
{
    if (!is(x, "RangedSummarizedExperiment"))
        x <- as(x, "RangedSummarizedExperiment")
    x <- BiocGenerics:::replaceSlots(x, ...,
             rowRanges=value,
             elementMetadata=new("DataFrame", nrows=length(value)),
             check=FALSE)
    msg <- .valid.SummarizedExperiment0.assays_nrow(x)
    if (!is.null(msg))
        stop(msg)
    x
}

setReplaceMethod("rowRanges", c("SummarizedExperiment0", "GenomicRanges"),
    .SummarizedExperiment0.rowRanges.replace)

setReplaceMethod("rowRanges", c("SummarizedExperiment0", "GRangesList"),
    .SummarizedExperiment0.rowRanges.replace)

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
    BiocGenerics:::replaceSlots(x,
        rowRanges=rowRanges,
        colData=colData,
        check=FALSE)
})


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Subsetting.
###

setMethod("subset", "RangedSummarizedExperiment",
    function(x, subset, select, ...)
{
    i <- S4Vectors:::evalqForSubset(subset, rowRanges(x), ...)
    j <- S4Vectors:::evalqForSubset(select, colData(x), ...)
    x[i, j]
})


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### To facilitate transition from "classic" SummarizedExperiment objects to
### new RangedSummarizedExperiment objects.
###

.has_SummarizedExperiment_internal_structure <- function(object)
    all(sapply(slotNames("SummarizedExperiment"), .hasSlot, object=object))

### Used in GenomicRanges!
.from_SummarizedExperiment_to_RangedSummarizedExperiment <- function(from)
    .new_RangedSummarizedExperiment(from@assays,
                                    from@rowData,
                                    from@colData,
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
        rse <- .new_RangedSummarizedExperiment(object@assays,
                                               object@rowRanges,
                                               object@colData,
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


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
## colData-as-GRanges compatibility: allow direct access to GRanges /
## GRangesList colData for select functions

## Not supported:
## 
## Not consistent SummarizedExperiment structure: length, names,
##   as.data.frame, c.
## Length-changing endomorphisms: disjoin, gaps, reduce, unique.
## 'legacy' data types / functions: as "RangedData", as "RangesList",
##   renameSeqlevels, keepSeqlevels.
## Possile to implement, but not yet: Ops, map, window, window<-

## mcols
setMethod(mcols, "RangedSummarizedExperiment",
    function(x, use.names=FALSE, ...)
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

setMethod(elementMetadata, "RangedSummarizedExperiment",
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
## compare / Compare 
## 
.RangedSummarizedExperiment.compare <-
    function(x, y)
{
    if (is(x, "RangedSummarizedExperiment"))
        x <- rowRanges(x)
    if (is(y, "RangedSummarizedExperiment"))
        y <- rowRanges(y)
    compare(x, y)
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
        setMethod("compare", .sig, .RangedSummarizedExperiment.compare)
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

setMethod("order", "RangedSummarizedExperiment",
    function(..., na.last = TRUE, decreasing = FALSE)
{
    args <- lapply(list(...), rowRanges)
    do.call("order",
            c(args, list(na.last=na.last, decreasing=decreasing)))
})

setMethod("rank", "RangedSummarizedExperiment",
    function (x, na.last = TRUE,
              ties.method = c("average", "first", "random", "max", "min"))
{
    ties.method <- match.arg(ties.method)
    rank(rowRanges(x), na.last=na.last, ties.method=ties.method)
})

setMethod("sort", "RangedSummarizedExperiment",
    function(x, decreasing = FALSE, ...)
{
    x[order(rowRanges(x), decreasing=decreasing),]
})

## seqinfo (also seqlevels, genome, seqlevels<-, genome<-), seqinfo<-

setMethod(seqinfo, "RangedSummarizedExperiment",
    function(x)
{
    seqinfo(x@rowRanges)
})

setReplaceMethod("seqinfo", "RangedSummarizedExperiment",
    function (x, new2old = NULL, force = FALSE, value)
{
    if (!is(value, "Seqinfo")) 
        stop("the supplied 'seqinfo' must be a Seqinfo object")
    dangling_seqlevels <-
        GenomeInfoDb:::getDanglingSeqlevels(x@rowRanges, new2old = new2old,
                                            force = force, seqlevels(value))
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

setMethod(split, "RangedSummarizedExperiment",
    function(x, f, drop=FALSE, ...) 
{
    splitAsList(x, f, drop=drop)
})
