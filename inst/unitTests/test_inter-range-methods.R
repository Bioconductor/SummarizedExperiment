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
    generic_functions <- c("isDisjoint", "disjointBins")
    for (fun in generic_functions) {
        generic <- getGeneric(fun)
        method <- getMethod(fun, "RangedSummarizedExperiment")
        checkIdentical("x", generic@signature)
        checkIdentical(formals(generic@.Data), formals(method@.Data))
    }
}

test_inter_range_methods <- function()
{
    #for (i in 1:2) {
    for (i in 1L) {
        x <- rseList[[i]]

        ## isDisjoint
        target <- isDisjoint(rowRanges(x))
        current <- isDisjoint(x)
        checkIdentical(target, current)

        ## disjointBins
        target <- disjointBins(rowRanges(x))
        current <- disjointBins(x)
        checkIdentical(target, current)
    }
}

