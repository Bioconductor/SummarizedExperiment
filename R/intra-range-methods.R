### =========================================================================
### Intra-range methods
### -------------------------------------------------------------------------
###


setMethod("shift", "RangedSummarizedExperiment",
    function(x, shift=0L, use.names=TRUE)
    {
        x0 <- x
        x <- rowRanges(x)
        rowRanges(x0) <- callGeneric()
        x0
    }
)

setMethod("narrow", "RangedSummarizedExperiment",
    function(x, start=NA, end=NA, width=NA, use.names=TRUE)
    {
        x0 <- x
        x <- rowRanges(x)
        rowRanges(x0) <- callGeneric()
        x0
    }
)

setMethod("resize", "RangedSummarizedExperiment",
    function(x, width, fix="start", use.names=TRUE, ignore.strand=FALSE)
    {
        x0 <- x
        x <- rowRanges(x)
        rowRanges(x0) <- callGeneric()
        x0
    }
)

setMethod("flank", "RangedSummarizedExperiment", 
    function(x, width, start=TRUE, both=FALSE, use.names=TRUE,
             ignore.strand=FALSE)
    {
        x0 <- x
        x <- rowRanges(x)
        rowRanges(x0) <- callGeneric()
        x0
    }
)

setMethod("promoters", "RangedSummarizedExperiment",
    function(x, upstream=2000, downstream=200)
    {
        x0 <- x
        x <- rowRanges(x)
        rowRanges(x0) <- callGeneric()
        x0
    }
)

### Because 'keep.all.ranges' is FALSE by default, it will break if some
### ranges are dropped.
setMethod("restrict", "RangedSummarizedExperiment",
    function(x, start=NA, end=NA, keep.all.ranges=FALSE, use.names=TRUE)
    {
        x0 <- x
        x <- rowRanges(x)
        rowRanges(x0) <- callGeneric()
        x0
    }
)

setMethod("trim", "RangedSummarizedExperiment",
    function(x, use.names=TRUE)
    {
        x0 <- x
        x <- rowRanges(x)
        rowRanges(x0) <- callGeneric()
        x0
    }
)

