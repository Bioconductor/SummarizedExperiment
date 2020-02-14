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


.GENERIC_SIGNATURES <- list(
    precede=c("x", "subject"),
    follow=c("x", "subject"),
    nearest=c("x", "subject"),
    distance=c("x", "y"),
    distanceToNearest=c("x", "subject")
)

test_interfaces <- function()
{
    method_signatures <- list(
        c("RangedSummarizedExperiment", "ANY"),
        c("ANY", "RangedSummarizedExperiment"),
        c("RangedSummarizedExperiment", "RangedSummarizedExperiment")
    )
    for (fun in names(.GENERIC_SIGNATURES)) {
        generic <- getGeneric(fun)
        checkIdentical(.GENERIC_SIGNATURES[[fun]], generic@signature)
        for (sig in method_signatures) {
            method <- getMethod(fun, sig)
            checkIdentical(formals(generic@.Data), formals(method@.Data))
        }
    }
}

test_nearest_methods <- function()
{
    #for (i in 1:2) {
    for (i in 1L) {
        x <- rseList[[i]]
        #for (j in 1:2) {
        for (j in 1L) {
            y <- rseList[[j]]
            for (fun in names(.GENERIC_SIGNATURES)) {
                fun <- get(fun)
                target <- fun(rowRanges(x), rowRanges(y))
                current <- fun(x, rowRanges(y))
                checkIdentical(target, current)
                current <- fun(rowRanges(x), y)
                checkIdentical(target, current)
                current <- fun(x, y)
                checkIdentical(target, current)
            }
        }
    }
}

