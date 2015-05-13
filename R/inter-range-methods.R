### =========================================================================
### Inter-range methods
### -------------------------------------------------------------------------
###


setMethod("isDisjoint", "GenomicRanges",
    function(x, ignore.strand=FALSE)
    {
        x <- rowRanges(x)
        callGeneric()
    }
)

setMethod("disjointBins", "GenomicRanges",
    function(x, ignore.strand = FALSE)
    {
        x <- rowRanges(x)
        callGeneric()
    }
)

