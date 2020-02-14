###

M1 <- matrix(1, 5, 3)
M2 <- matrix(1, 3, 3)
assaysList <- list(gr=SimpleList(m=M1), grl=SimpleList(m=M2))
rowRangesList <-
    list(gr=GRanges("chr1", IRanges(1:5, 10)),
         grl=split(GRanges("chr1", IRanges(1:5, 10)), c(1,1,2,2,3)))
names(rowRangesList[["grl"]]) <- NULL
colData <- DataFrame(x=letters[1:3])

## a list of one SE with GRanges and one with GRangesList
rseList <-
    list(SummarizedExperiment(
           assays=assaysList[["gr"]],
           rowRanges=rowRangesList[["gr"]],
           colData=colData),
         SummarizedExperiment(
           assays=assaysList[["grl"]],
           rowRanges=rowRangesList[["grl"]],
           colData=colData))


test_interfaces <- function()
{
    fun <- "findOverlaps"
    signatures <- list(
        c("RangedSummarizedExperiment", "Vector"),
        c("Vector", "RangedSummarizedExperiment"),
        c("RangedSummarizedExperiment", "RangedSummarizedExperiment")
    )
    generic <- getGeneric(fun)
    for (sig in signatures) {
        method <- getMethod(fun, sig)
        checkIdentical(c("query", "subject"), generic@signature)
        checkIdentical(formals(generic@.Data), formals(method@.Data))
    }
}

test_findOverlaps_methods <- function()
{
    identical_SummarizedExperiment <- function(x, y) {
        x@assays <- as(assays(x), "SimpleAssays")
        y@assays <- as(assays(y), "SimpleAssays")
        identical(x, y)
    }
    for (i in 1:2) {
        x <- rseList[[i]]
        for (j in 1:2) {
            y <- rseList[[j]]

            ## findOverlaps
            target <- findOverlaps(rowRanges(x), rowRanges(y))
            current <- findOverlaps(x, rowRanges(y))
            checkIdentical(target, current)
            current <- findOverlaps(rowRanges(x), y)
            checkIdentical(target, current)
            current <- findOverlaps(x, y)
            checkIdentical(target, current)

            ## countOverlaps
            target <- countOverlaps(rowRanges(x), rowRanges(y))
            current <- countOverlaps(x, rowRanges(y))
            checkIdentical(target, current)
            current <- countOverlaps(rowRanges(x), y)
            checkIdentical(target, current)
            current <- countOverlaps(x, y)
            checkIdentical(target, current)

            ## overlapsAny
            target <- overlapsAny(rowRanges(x), rowRanges(y))
            current <- overlapsAny(x, rowRanges(y))
            checkIdentical(target, current)
            current <- overlapsAny(rowRanges(x), y)
            checkIdentical(target, current)
            current <- overlapsAny(x, y)
            checkIdentical(target, current)

            ## subsetByOverlaps
            target <- subsetByOverlaps(x, rowRanges(y))
            current <- subsetByOverlaps(x, rowRanges(y))
            checkTrue(identical_SummarizedExperiment(target, current))
            current <- subsetByOverlaps(x, y)
            checkTrue(identical_SummarizedExperiment(target, current))

            target <- subsetByOverlaps(rowRanges(x), rowRanges(y))
            current <- subsetByOverlaps(rowRanges(x), y)
            checkIdentical(target, current)
        }
    }
}

