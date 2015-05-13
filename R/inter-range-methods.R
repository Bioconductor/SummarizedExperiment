### =========================================================================
### Inter-range methods
### -------------------------------------------------------------------------
###


setMethod("isDisjoint", "RangedSummarizedExperiment",
    function(x, ignore.strand=FALSE)
    {
        x <- rowRanges(x)
        callGeneric()
    }
)

setMethod("disjointBins", "RangedSummarizedExperiment",
    function(x, ignore.strand = FALSE)
    {
        x <- rowRanges(x)
        callGeneric()
    }
)

