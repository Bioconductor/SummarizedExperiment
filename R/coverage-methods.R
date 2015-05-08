### =========================================================================
### "coverage" method
### -------------------------------------------------------------------------
###


setMethod("coverage", "RangedSummarizedExperiment",
    function(x, shift=0L, width=NULL, weight=1L,
                method=c("auto", "sort", "hash"))
    {
        x <- rowRanges(x)
        callGeneric()
    }
)

