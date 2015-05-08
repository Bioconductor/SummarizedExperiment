library(digest)

.singleDispatch <-
    c("duplicated", "end", "end<-", "granges", "ranges", 
      "seqinfo", "seqinfo<-", "seqnames", "start", "start<-",
      "strand", "width", "width<-")

.twoDispatch <- c("compare", "Compare")

.otherFuns <- c("order", "rank", "sort")

m <- matrix(1, 5, 3, dimnames=list(NULL, NULL))
mlst <- matrix(1, 3, 3, dimnames=list(NULL, NULL))
assaysList <- list(gr=SimpleList(m=m), grl=SimpleList(m=mlst))
rowRangesList <- 
    list(gr=GRanges("chr1", IRanges(1:5, 10)), 
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

test_RangedSummarizedExperiment_GRanges_API <- function()
{
    ## are we targetting the correct API? signature for
    ## RangedSummarizedExperiment method should match signature for
    ## GenomicRanges or similar, as in each test below

    for (.fun in .singleDispatch) {
        generic <- getGeneric(.fun)
        method <- getMethod(.fun, "RangedSummarizedExperiment")
        checkIdentical("x", generic@signature)
        checkIdentical(formals(generic@.Data), formals(method@.Data))
    }

    ## FIXME: compare, Compare

    .sig <- "RangedSummarizedExperiment"
    for (.fun in .otherFuns) {
        generic <- getGeneric(.fun)
        method <- getMethod(.fun, "RangedSummarizedExperiment")
        checkIdentical(formals(generic@.Data), formals(method@.Data))
    }        
}

test_RangedSummarizedExperiment_GRanges_values <- function()
{
    x <- ssetList[[1]]
    isAssign <- grep("<-$", .singleDispatch, value=TRUE)
    .funs <- setdiff(.singleDispatch, isAssign)
    ## 'exp' created after manual inspection of results
    exp <- setNames(c("02dde", "80339", "49a3f", "86757", "77198",
                      "ec53a", "35e2c", "625d9", "3c90a"), .funs)
    obs <- sapply(.funs, function(.fun) {
        substr(digest(getGeneric(.fun)(x)), 1, 5)
    })
    checkIdentical(exp, obs)

    .funs <- isAssign
    .gets <- sub("<-$", "", isAssign)
    for (i in seq_along(isAssign)) {
        ## self-assignment isomorphism
        value <- getGeneric(.gets[[i]])(x)
        x1 <- do.call(isAssign[[i]], list(x, value=value))
        checkIdentical(x, x1)
    }
}

test_RangedSummarizedExperiment_split <- function()
{
    gr <- GRanges(Rle(c("A", "B"), c(2, 3)), IRanges(1:5, 10))
    se <- SummarizedExperiment(m, rowRanges=gr, colData=colData)
    ## FIXME: unname should not be necessary
    obs <- split(se, seqnames(se))
    exp <- SimpleList(A=se[1:2], B=se[3:5])
    checkEquals(obs, exp)
}

