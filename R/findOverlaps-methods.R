### =========================================================================
### findOverlaps methods
### -------------------------------------------------------------------------


### findOverlaps

setMethod("findOverlaps", c("RangedSummarizedExperiment", "Vector"),
    function(query, subject, maxgap=0L, minoverlap=1L,
             type=c("any", "start", "end", "within", "equal"),
             select=c("all", "first", "last", "arbitrary"),
             algorithm=c("nclist", "intervaltree"),
             ignore.strand=FALSE)
    {
        query <- rowRanges(query)
        callGeneric()
    }
)

setMethod("findOverlaps", c("Vector", "RangedSummarizedExperiment"),
    function(query, subject, maxgap=0L, minoverlap=1L,
             type=c("any", "start", "end", "within", "equal"),
             select=c("all", "first", "last", "arbitrary"),
             algorithm=c("nclist", "intervaltree"),
             ignore.strand=FALSE)
    {
        subject <- rowRanges(subject)
        callGeneric()
    }
)

setMethod("findOverlaps", c("RangedSummarizedExperiment",
                            "RangedSummarizedExperiment"),
    function(query, subject, maxgap=0L, minoverlap=1L,
             type=c("any", "start", "end", "within", "equal"),
             select=c("all", "first", "last", "arbitrary"),
             algorithm=c("nclist", "intervaltree"),
             ignore.strand=FALSE)
    {
        query <- rowRanges(query)
        subject <- rowRanges(subject)
        callGeneric()
    }
)

### countOverlaps

setMethod("countOverlaps", c("RangedSummarizedExperiment", "Vector"),
    function(query, subject, maxgap=0L, minoverlap=1L,
             type=c("any", "start", "end", "within", "equal"),
             algorithm=c("nclist", "intervaltree"),
             ignore.strand=FALSE)
    {
        query <- rowRanges(query)
        callGeneric()
    }
)

setMethod("countOverlaps", c("Vector", "RangedSummarizedExperiment"),
    function(query, subject, maxgap=0L, minoverlap=1L,
             type=c("any", "start", "end", "within", "equal"),
             algorithm=c("nclist", "intervaltree"),
             ignore.strand=FALSE)
    {
        subject <- rowRanges(subject)
        callGeneric()
    }
)

setMethod("countOverlaps", c("RangedSummarizedExperiment",
                             "RangedSummarizedExperiment"),
    function(query, subject, maxgap=0L, minoverlap=1L,
             type=c("any", "start", "end", "within", "equal"),
             algorithm=c("nclist", "intervaltree"),
             ignore.strand=FALSE)
    {
        query <- rowRanges(query)
        subject <- rowRanges(subject)
        callGeneric()
    }
)

### overlapsAny & subsetByOverlaps

.signatures2 <- list(
    c("RangedSummarizedExperiment", "Vector"),
    c("Vector", "RangedSummarizedExperiment"),
    c("RangedSummarizedExperiment", "RangedSummarizedExperiment")
)

setMethods("overlapsAny", .signatures2,
    GenomicRanges:::overlapsAny.definition
)

setMethods("subsetByOverlaps", .signatures2,
    GenomicRanges:::subsetByOverlaps.definition1
)

