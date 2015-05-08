### =========================================================================
### nearest (and related) methods
### -------------------------------------------------------------------------
###


### precede & follow

for (f in c("precede", "follow")) {
    setMethod(f, c("RangedSummarizedExperiment", "ANY"),
        function(x, subject, select=c("arbitrary", "all"), ignore.strand=FALSE)
        {
            query <- rowRanges(query)
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
            query <- rowRanges(query)
            subject <- rowRanges(subject)
            callGeneric()
        }
    )
}

### nearest

setMethod("nearest", c("RangedSummarizedExperiment", "ANY"),
    function(x, subject, select=c("arbitrary", "all"),
             algorithm=c("nclist", "intervaltree"), ignore.strand=FALSE)
    {
        query <- rowRanges(query)
        callGeneric()
    }
)

setMethod("nearest", c("ANY", "RangedSummarizedExperiment"),
    function(x, subject, select=c("arbitrary", "all"),
             algorithm=c("nclist", "intervaltree"), ignore.strand=FALSE)
    {
        subject <- rowRanges(subject)
        callGeneric()
    }
)

setMethod("nearest", c("RangedSummarizedExperiment",
                       "RangedSummarizedExperiment"),
    function(x, subject, select=c("arbitrary", "all"),
             algorithm=c("nclist", "intervaltree"), ignore.strand=FALSE)
    {
        query <- rowRanges(query)
        subject <- rowRanges(subject)
        callGeneric()
    }
)

### distance

setMethod("distance", c("RangedSummarizedExperiment", "ANY"),
    function(x, y, ignore.strand=FALSE, ...)
    {
        query <- rowRanges(query)
        callGeneric()
    }
)

setMethod("distance", c("ANY", "RangedSummarizedExperiment"),
    function(x, y, ignore.strand=FALSE, ...)
    {
        subject <- rowRanges(subject)
        callGeneric()
    }
)

setMethod("distance", c("RangedSummarizedExperiment",
                        "RangedSummarizedExperiment"),
    function(x, y, ignore.strand=FALSE, ...)
    {
        query <- rowRanges(query)
        subject <- rowRanges(subject)
        callGeneric()
    }
)

### distanceToNearest

setMethod("distanceToNearest", c("RangedSummarizedExperiment", "ANY"),
    function(x, subject, algorithm=c("nclist", "intervaltree"),
             ignore.strand=FALSE, ...)
    {
        query <- rowRanges(query)
        callGeneric()
    }
)

setMethod("distanceToNearest", c("ANY", "RangedSummarizedExperiment"),
    function(x, subject, algorithm=c("nclist", "intervaltree"),
             ignore.strand=FALSE, ...)
    {
        subject <- rowRanges(subject)
        callGeneric()
    }
)

setMethod("distanceToNearest", c("RangedSummarizedExperiment",
                                 "RangedSummarizedExperiment"),
    function(x, subject, algorithm=c("nclist", "intervaltree"),
             ignore.strand=FALSE, ...)
    {
        query <- rowRanges(query)
        subject <- rowRanges(subject)
        callGeneric()
    }
)

