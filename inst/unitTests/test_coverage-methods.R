###

M1 <- matrix(1, 5, 3, dimnames=list(NULL, NULL))
M2 <- matrix(1, 3, 3, dimnames=list(NULL, NULL))
mList <- list(M1, M2)
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

test_coverage_RangedSummarizedExperiment <- function()
{
    generic <- getGeneric("coverage")
    method <- getMethod("coverage", "RangedSummarizedExperiment")
    checkIdentical("x", generic@signature)
    checkIdentical(formals(generic@.Data), formals(method@.Data))

    for (i in 1:2) {
        target <- coverage(rowRanges(rseList[[i]]))
        current <- coverage(rseList[[i]])
        checkIdentical(target, current)

        weight <- runif(length(rseList[[i]]))
        target <- coverage(rowRanges(rseList[[i]]), weight=weight)
        current <- coverage(rseList[[i]], weight=weight)
        checkIdentical(target, current)
    }
}

