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
    GenomicRanges:::clone(x, rowRanges=local({
        r <- rowRanges(x)
        mcols(r) <- value
        r
    }))
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

setMethod(values, "RangedSummarizedExperiment",
    function(x, ...)
{
    values(rowRanges(x), ...)
})

setReplaceMethod("values", "RangedSummarizedExperiment",
    function(x, ..., value)
{
    values(rowRanges(x), ...) <- value
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

## subset

setMethod("subset", "RangedSummarizedExperiment",
    function(x, subset, select, ...)
{
    i <- S4Vectors:::evalqForSubset(subset, rowRanges(x), ...)
    j <- S4Vectors:::evalqForSubset(select, colData(x), ...)
    x[i, j]
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

## extractROWS, replaceROWS

setMethod("extractROWS", "RangedSummarizedExperiment",
    function(x, i)
{
    ridx <- extractROWS(seq_len(nrow(x)), i)
    x[ridx,]
})

setMethod("replaceROWS", "RangedSummarizedExperiment",
    function (x, i, value)
{
    ridx <- extractROWS(seq_len(nrow(x)), i)
    x[ridx,] <- value
    x
})

setMethod(split, "RangedSummarizedExperiment",
    function(x, f, drop=FALSE, ...) 
{
    splitAsList(x, f, drop=drop)
})
