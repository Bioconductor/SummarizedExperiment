###

M1 <- matrix(1, 5, 3)
M2 <- matrix(1, 3, 3)
assaysList <- list(gr=SimpleList(m=M1), grl=SimpleList(m=M2))
rowRangesList <-
    list(gr=GRanges("chr1", IRanges(1:5, 10), Rle(c("+", "-"), 3:2)),
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
    generic_functions <- c("shift", "narrow", "resize",
                           "flank", "promoters",
                           "restrict", "trim")
    for (fun in generic_functions) {
        generic <- getGeneric(fun)
        method <- getMethod(fun, "RangedSummarizedExperiment")
        checkIdentical("x", generic@signature)
        checkIdentical(formals(generic@.Data), formals(method@.Data))
    }
}

test_intra_range_methods <- function()
{
    identical_SummarizedExperiment <- function(x, y) {
        x@assays <- as(assays(x), "SimpleAssays")
        y@assays <- as(assays(y), "SimpleAssays")
        identical(x, y)
    }
    #for (i in 1:2) {
    for (i in 1L) {
        ## shift
        target <- rseList[[i]]
        rowRanges(target) <- shift(rowRanges(target), 50)
        current <- shift(rseList[[i]], 50)
        checkTrue(identical_SummarizedExperiment(target, current))

        ## narrow
        target <- rseList[[i]]
        rowRanges(target) <- narrow(rowRanges(target), 2, -2)
        current <- narrow(rseList[[i]], 2, -2)
        checkTrue(identical_SummarizedExperiment(target, current))

        ## resize
        target <- rseList[[i]]
        rowRanges(target) <- resize(rowRanges(target), 8)
        current <- resize(rseList[[i]], 8)
        checkTrue(identical_SummarizedExperiment(target, current))

        ## flank
        target <- rseList[[i]]
        rowRanges(target) <- flank(rowRanges(target), 5, both=TRUE)
        current <- flank(rseList[[i]], 5, both=TRUE)
        checkTrue(identical_SummarizedExperiment(target, current))

        ## promoters
        target <- rseList[[i]]
        rowRanges(target) <- promoters(rowRanges(target),
                                       upstream=20, downstream=5)
        current <- promoters(rseList[[i]], upstream=20, downstream=5)
        checkTrue(identical_SummarizedExperiment(target, current))

        ## restrict
        target <- rseList[[i]]
        rowRanges(target) <- restrict(rowRanges(target), start=2, end=3,
                                      keep.all.ranges=TRUE)
        current <- restrict(rseList[[i]], start=2, end=3,
                            keep.all.ranges=TRUE)
        checkTrue(identical_SummarizedExperiment(target, current))

        ## trim
        suppressWarnings(seqlengths(rseList[[i]]) <- 8)
        target <- rseList[[i]]
        rowRanges(target) <- trim(rowRanges(target))
        current <- trim(rseList[[i]])
        checkTrue(identical_SummarizedExperiment(target, current))
        seqlengths(rseList[[i]]) <- NA
    }
}

