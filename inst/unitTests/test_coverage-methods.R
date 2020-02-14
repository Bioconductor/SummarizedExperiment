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
    generic_functions <- "coverage"
    for (fun in generic_functions) {
        generic <- getGeneric(fun)
        method <- getMethod(fun, "RangedSummarizedExperiment")
        checkIdentical("x", generic@signature)
        checkIdentical(formals(generic@.Data), formals(method@.Data))
    }
}

test_coverage_RangedSummarizedExperiment <- function()
{
    for (i in 1:2) {
        x <- rseList[[i]]

        target <- coverage(rowRanges(x))
        current <- coverage(x)
        checkIdentical(target, current)

        weight <- runif(length(x))
        ## Issues a warning (in BioC 3.3) when rowRanges(x) is a GRangesList
        ## object, which reveals a problem with how the "coverage" method for
        ## GRangesList objects handles the 'weight' argument. The warning is
        ## expected and healthy, don't try to suppress it here. It will go
        ## away when we fix the "coverage" method for GRangesList objects
        ## (defined in the GenomicRanges package). 
        target <- coverage(rowRanges(x), weight=weight)
        current <- coverage(x, weight=weight)
        checkIdentical(target, current)
    }
}

