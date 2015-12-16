### =========================================================================
### nearest (and related) methods
### -------------------------------------------------------------------------
###


### precede & follow

for (f in c("precede", "follow")) {
    setMethod(f, c("RangedSummarizedExperiment", "ANY"),
        function(x, subject, select=c("arbitrary", "all"), ignore.strand=FALSE)
        {
            x <- rowRanges(x)
            callGeneric()
        }
    )
    setMethod(f, c("ANY", "RangedSummarizedExperiment"),
        function(x, subject, select=c("arbitrary", "all"), ignore.strand=FALSE)
        {
            subject <- rowRanges(subject)
            callGeneric()
        }
    )
    setMethod(f, c("RangedSummarizedExperiment", "RangedSummarizedExperiment"),
        function(x, subject, select=c("arbitrary", "all"), ignore.strand=FALSE)
        {
            x <- rowRanges(x)
            subject <- rowRanges(subject)
            callGeneric()
        }
    )
}

### nearest

setMethod("nearest", c("RangedSummarizedExperiment", "ANY"),
    function(x, subject, select=c("arbitrary", "all"), ignore.strand=FALSE)
    {
        x <- rowRanges(x)
        callGeneric()
    }
)

setMethod("nearest", c("ANY", "RangedSummarizedExperiment"),
    function(x, subject, select=c("arbitrary", "all"), ignore.strand=FALSE)
    {
        subject <- rowRanges(subject)
        callGeneric()
    }
)

setMethod("nearest", c("RangedSummarizedExperiment",
                       "RangedSummarizedExperiment"),
    function(x, subject, select=c("arbitrary", "all"), ignore.strand=FALSE)
    {
        x <- rowRanges(x)
        subject <- rowRanges(subject)
        callGeneric()
    }
)

### distance

setMethod("distance", c("RangedSummarizedExperiment", "ANY"),
    function(x, y, ignore.strand=FALSE, ...)
    {
        x <- rowRanges(x)
        callGeneric()
    }
)

setMethod("distance", c("ANY", "RangedSummarizedExperiment"),
    function(x, y, ignore.strand=FALSE, ...)
    {
        y <- rowRanges(y)
        callGeneric()
    }
)

setMethod("distance", c("RangedSummarizedExperiment",
                        "RangedSummarizedExperiment"),
    function(x, y, ignore.strand=FALSE, ...)
    {
        x <- rowRanges(x)
        y <- rowRanges(y)
        callGeneric()
    }
)

### distanceToNearest

setMethod("distanceToNearest", c("RangedSummarizedExperiment", "ANY"),
    function(x, subject, ignore.strand=FALSE, ...)
    {
        x <- rowRanges(x)
        callGeneric()
    }
)

setMethod("distanceToNearest", c("ANY", "RangedSummarizedExperiment"),
    function(x, subject, ignore.strand=FALSE, ...)
    {
        subject <- rowRanges(subject)
        callGeneric()
    }
)

setMethod("distanceToNearest", c("RangedSummarizedExperiment",
                                 "RangedSummarizedExperiment"),
    function(x, subject, ignore.strand=FALSE, ...)
    {
        x <- rowRanges(x)
        subject <- rowRanges(subject)
        callGeneric()
    }
)

