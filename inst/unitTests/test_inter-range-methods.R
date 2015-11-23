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


test_inter_range_methods <- function()
{
    #for (i in 1:2) {
    for (i in 1L) {
        ## isDisjoint
        x <- rseList[[i]]
        target <- isDisjoint(rowRanges(x))
        current <- isDisjoint(x)
        checkIdentical(target, current)

        ## disjointBins
        x <- rseList[[i]]
        target <- disjointBins(rowRanges(x))
        current <- disjointBins(x)
        checkIdentical(target, current)
    }
}

