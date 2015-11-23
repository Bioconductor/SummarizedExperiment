###

M1 <- matrix(1, 5, 3, dimnames=list(NULL, NULL))
M2 <- matrix(1, 3, 3, dimnames=list(NULL, NULL))
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


test_intra_range_methods <- function()
{
    #for (i in 1:2) {
    for (i in 1L) {
        ## shift
        target <- rseList[[i]]
        rowRanges(target) <- shift(rowRanges(target), 50)
        current <- shift(rseList[[i]], 50)
        checkIdentical(target, current)

        ## narrow
        target <- rseList[[i]]
        rowRanges(target) <- narrow(rowRanges(target), 2, -2)
        current <- narrow(rseList[[i]], 2, -2)
        checkIdentical(target, current)

        ## resize
        target <- rseList[[i]]
        rowRanges(target) <- resize(rowRanges(target), 8)
        current <- resize(rseList[[i]], 8)
        checkIdentical(target, current)

        ## flank
        target <- rseList[[i]]
        rowRanges(target) <- flank(rowRanges(target), 5, both=TRUE)
        current <- flank(rseList[[i]], 5, both=TRUE)
        checkIdentical(target, current)

        ## promoters
        target <- rseList[[i]]
        rowRanges(target) <- promoters(rowRanges(target),
                                       upstream=20, downstream=5)
        current <- promoters(rseList[[i]], upstream=20, downstream=5)
        checkIdentical(target, current)

        ## restrict
        target <- rseList[[i]]
        rowRanges(target) <- restrict(rowRanges(target), start=2, end=3,
                                      keep.all.ranges=TRUE)
        current <- restrict(rseList[[i]], start=2, end=3,
                            keep.all.ranges=TRUE)
        checkIdentical(target, current)

        ## trim
        suppressWarnings(seqlengths(rseList[[i]]) <- 8)
        target <- rseList[[i]]
        rowRanges(target) <- trim(rowRanges(target))
        current <- trim(rseList[[i]])
        checkIdentical(target, current)
        seqlengths(rseList[[i]]) <- NA
    }
}

