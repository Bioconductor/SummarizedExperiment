###

m <- matrix(1, 5, 3, dimnames=list(NULL, NULL))
mlst <- matrix(1, 3, 3, dimnames=list(NULL, NULL))
assaysList <- list(gr=SimpleList(m=m), grl=SimpleList(m=mlst))
rowRangesList <-
    list(gr=GRanges("chr1", IRanges(1:5, 10), Rle(c("+", "-"), 3:2)),
         grl=split(GRanges("chr1", IRanges(1:5, 10)), c(1,1,2,2,3)))
names(rowRangesList[["grl"]]) <- NULL
colData <- DataFrame(x=letters[1:3])

## a list of one SE with GRanges and one with GRangesList
ssetList <-
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
        target <- ssetList[[i]]
        rowRanges(target) <- shift(rowRanges(target), 50)
        current <- shift(ssetList[[i]], 50)
        checkIdentical(target, current)

        ## narrow
        target <- ssetList[[i]]
        rowRanges(target) <- narrow(rowRanges(target), 2, -2)
        current <- narrow(ssetList[[i]], 2, -2)
        checkIdentical(target, current)

        ## resize
        target <- ssetList[[i]]
        rowRanges(target) <- resize(rowRanges(target), 8)
        current <- resize(ssetList[[i]], 8)
        checkIdentical(target, current)

        ## flank
        target <- ssetList[[i]]
        rowRanges(target) <- flank(rowRanges(target), 5, both=TRUE)
        current <- flank(ssetList[[i]], 5, both=TRUE)
        checkIdentical(target, current)

        ## promoters
        target <- ssetList[[i]]
        rowRanges(target) <- promoters(rowRanges(target),
                                       upstream=20, downstream=5)
        current <- promoters(ssetList[[i]], upstream=20, downstream=5)
        checkIdentical(target, current)

        ## restrict
        target <- ssetList[[i]]
        rowRanges(target) <- restrict(rowRanges(target), start=2, end=3,
                                      keep.all.ranges=TRUE)
        current <- restrict(ssetList[[i]], start=2, end=3,
                            keep.all.ranges=TRUE)
        checkIdentical(target, current)

        ## trim
        suppressWarnings(seqlengths(ssetList[[i]]) <- 8)
        target <- ssetList[[i]]
        rowRanges(target) <- trim(rowRanges(target))
        current <- trim(ssetList[[i]])
        checkIdentical(target, current)
        seqlengths(ssetList[[i]]) <- NA
    }
}

