library(digest)

.singleDispatch <-
    c("duplicated", "end", "end<-", "granges", "ranges", 
      "seqinfo", "seqinfo<-", "seqnames", "start", "start<-",
      "strand", "width", "width<-")

.twoDispatch <- c("pcompare", "Compare")

.otherFuns <- c("is.unsorted", "order", "rank", "sort")

M1 <- matrix(1, 5, 3)
M2 <- matrix(1, 3, 3)
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


test_RangedSummarizedExperiment_construction <- function()
{
    ## empty-ish
    m1 <- matrix(0, 0, 0)
    checkTrue(validObject(new("RangedSummarizedExperiment")))

    ## substance
    for (i in seq_along(rseList)) {
        rse <- rseList[[i]]
        checkTrue(validObject(rse))
        checkIdentical(SimpleList(m=mList[[i]]), assays(rse))
        checkIdentical(rowRangesList[[i]], rowRanges(rse))
        checkIdentical(DataFrame(x=letters[1:3]), colData(rse))
    }

    ## array in assays slot
    ss <- rseList[[1]]
    assays(ss) <- SimpleList(array(1:5, c(5,3,2)))
    checkTrue(validObject(ss))
    checkTrue(all(dim(assays(ss[1:3,1:2])[[1]]) == c(3, 2, 2)))
}

test_RangedSummarizedExperiment_getters <- function()
{
    for (i in seq_along(rseList)) {
        rse <- rseList[[i]]
        rowRanges <- rowRangesList[[i]]

        ## dim, dimnames
        checkIdentical(c(length(rowRanges), nrow(colData)), dim(rse))
        checkIdentical(NULL, dimnames(rse))

        ## row / col / metadata
        checkIdentical(rowRanges, rowRanges(rse))
        checkIdentical(colData, colData(rse))
        checkIdentical(list(), metadata(rse))
    }
}

test_RangedSummarizedExperiment_setters <- function()
{
    for (i in seq_along(rseList)) {
        rse <- rseList[[i]]
        rowRanges <- rowRangesList[[i]]

        ## row / col / metadata<-
        ss1 <- rse
        revData <- rev(rowRanges)
        rowRanges(ss1) <- revData
        checkIdentical(revData, rowRanges(ss1))
        checkException(rowRanges(ss1) <- rowRanges(rse)[1:2,,drop=FALSE],
                       "incorrect row dimensions", TRUE)
        revData <- colData[rev(seq_len(nrow(colData))),,drop=FALSE]
        colData(ss1) <- revData
        checkIdentical(revData, colData(ss1))
        checkException(colData(ss1) <- colData(rse)[1:2,,drop=FALSE],
                       "incorrect col dimensions", TRUE)
        lst <- list("foo", "bar")
        metadata(ss1) <- lst
        checkIdentical(lst, metadata(ss1))

        ## assay / assays
        ss1 <- rse
        assay(ss1) <- assay(ss1)+1
        checkIdentical(assay(rse)+1, assay(ss1))
        ss1 <- rse
        assay(ss1, 1) <- assay(ss1, 1) + 1
        checkIdentical(assay(rse, "m") + 1, assay(ss1, "m"))
        ss1 <- rse
        assay(ss1, "m") <- assay(ss1, "m") + 1
        checkIdentical(assay(rse, "m")+1, assay(ss1, "m"))

        ## dimnames<-
        ss1 <- rse
        dimnames <- list(letters[seq_len(nrow(ss1))],
                         LETTERS[seq_len(ncol(ss1))])
        rownames(ss1) <- dimnames[[1]]
        colnames(ss1) <- dimnames[[2]]
        checkIdentical(dimnames, dimnames(ss1))
        rowRanges1 <- rowRanges
        names(rowRanges1) <- dimnames[[1]]
        checkIdentical(rowRanges1, rowRanges(ss1))
        colData1 <- colData
        row.names(colData1) <- dimnames[[2]]
        checkIdentical(colData1, colData(ss1))
        ss1 <- rse
        dimnames(ss1) <- dimnames
        checkIdentical(dimnames, dimnames(ss1))
        dimnames(ss1) <- NULL
        checkIdentical(NULL, dimnames(ss1))
    }
}

test_RangedSummarizedExperiment_subset <- function()
{
    for (i in seq_along(rseList)) {
        rse <- rseList[[i]]
        rowRanges <- rowRangesList[[i]]

        ## numeric
        ss1 <- rse[2:3,]
        checkIdentical(c(2L, ncol(rse)), dim(ss1))
        checkIdentical(rowRanges(ss1), rowRanges(rse)[2:3,])
        checkIdentical(colData(ss1), colData(rse))
        ss1 <- rse[,2:3]
        checkIdentical(c(nrow(rse), 2L), dim(ss1))
        checkIdentical(rowRanges(ss1), rowRanges(rse))
        checkIdentical(colData(ss1), colData(rse)[2:3,,drop=FALSE])
        ss1 <- rse[2:3, 2:3]
        checkIdentical(c(2L, 2L), dim(ss1))
        checkIdentical(rowRanges(ss1), rowRanges(rse)[2:3,,drop=FALSE])
        checkIdentical(colData(ss1), colData(rse)[2:3,,drop=FALSE])

        ## character
        ss1 <- rse
        dimnames(ss1) <- list(LETTERS[seq_len(nrow(ss1))],
                              letters[seq_len(ncol(ss1))])
        ridx <- c("B", "C")
        checkIdentical(rowRanges(ss1[ridx,]), rowRanges(ss1)[ridx,])
        checkIdentical(rowRanges(ss1["C",]), rowRanges(ss1)["C",,drop=FALSE])
        checkException(ss1[LETTERS,], "i-index out of bounds", TRUE)
        cidx <- c("b", "c")
        checkIdentical(colData(ss1[,cidx]), colData(ss1)[cidx,,drop=FALSE])
        checkIdentical(colData(ss1[,"a"]), colData(ss1)["a",,drop=FALSE])
        checkException(ss1[,letters], "j-index out of bounds", TRUE)

        ## logical
        ss1 <- rse
        dimnames(ss1) <- list(LETTERS[seq_len(nrow(ss1))],
                              letters[seq_len(ncol(ss1))])
        checkEquals(ss1, ss1[TRUE,])
        checkIdentical(c(0L, ncol(ss1)), dim(ss1[FALSE,]))
        checkEquals(ss1, ss1[,TRUE])
        checkIdentical(c(nrow(ss1), 0L), dim(ss1[,FALSE]))
        idx <- c(TRUE, FALSE)               # recycling
        ss2 <- ss1[idx,]
        checkIdentical(rowRanges(ss1)[idx,,drop=FALSE], rowRanges(ss2))
        ss2 <- ss1[,idx]
        checkIdentical(colData(ss1)[idx,,drop=FALSE], colData(ss2))

        ## Rle
        ss1 <- rse
        rle <- rep(c(TRUE, FALSE), each=3, length.out=nrow(ss1))
        checkIdentical(rowRanges(ss1[rle]), rowRanges(ss1[Rle(rle)]))
        checkIdentical(assays(ss1[rle]), assays(ss1[Rle(rle)]))
    }

    ## 0 columns
    se <- SummarizedExperiment(rowRanges=GRanges("chr1", IRanges(1:10, width=1)))
    checkIdentical(dim(se[1:5, ]), c(5L, 0L))
    ## 0 rows 
    se <- SummarizedExperiment(colData=DataFrame(samples=1:10))
    checkIdentical(dim(se[ ,1:5]), c(0L, 5L))
}

test_RangedSummarizedExperiment_subsetassign <- function()
{
    for (i in seq_along(rseList)) {
        rse <- rseList[[i]]
        dimnames(rse) <- list(LETTERS[seq_len(nrow(rse))],
                              letters[seq_len(ncol(rse))])
        ## rows
        ss1 <- rse
        ss1[1:2,] <- ss1[2:1,]
        checkIdentical(rowRanges(rse)[2:1,], rowRanges(ss1)[1:2,])
        checkIdentical(rowRanges(rse[-(1:2),]), rowRanges(ss1)[-(1:2),])
        checkIdentical(colData(rse), colData(ss1))
        checkIdentical(c(metadata(rse), metadata(rse)), metadata(ss1))
        ## Rle
        ss1rle <- ss1Rle <- rse
        rle <- rep(c(TRUE, FALSE), each=3, length.out=nrow(ss1))
        ss1rle[rle,] <- ss1rle[rle,]
        ss1Rle[Rle(rle),] <- ss1Rle[Rle(rle),]
        checkIdentical(rowRanges(ss1rle), rowRanges(ss1Rle))
        checkIdentical(assays(ss1rle), assays(ss1Rle))
        ## cols
        ss1 <- rse
        ss1[,1:2] <- ss1[,2:1,drop=FALSE]
        checkIdentical(colData(rse)[2:1,,drop=FALSE],
                       colData(ss1)[1:2,,drop=FALSE])
        checkIdentical(colData(rse)[-(1:2),,drop=FALSE],
                       colData(ss1)[-(1:2),,drop=FALSE])
        checkIdentical(rowRanges(rse), rowRanges(ss1))
        checkIdentical(c(metadata(rse), metadata(rse)), metadata(ss1))
    }

    ## full replacement
    ss1 <- ss2 <- rseList[[1]]
    rowRanges(ss2) <- rev(rowRanges(ss2))
    ss1[,] <- ss2
    checkIdentical(ss1, ss2)
}

quiet <- suppressWarnings
test_RangedSummarizedExperiment_cbind <- function()
## requires matching ranges
{
    ## empty
    se <- SummarizedExperiment()
    empty <- cbind(se, se)
    checkTrue(all.equal(se, empty))

    ## different ranges 
    se1 <- rseList[[1]]
    se2 <- se1[2:4]
    rownames(se2) <- month.name[seq_len(nrow(se2))]
    checkException(quiet(cbind(se1, se2)), silent=TRUE)

    ## same ranges 
    se1 <- rseList[[1]]
    se2 <- se1[,1:2]
    colnames(se2) <- month.name[seq_len(ncol(se2))]
    res <- cbind(se1, se2)
    checkTrue(nrow(res) == 5)
    checkTrue(ncol(res) == 5)
    ## rowRanges
    rowData(se1) <- DataFrame("one"=1:5)
    rowData(se2) <- DataFrame("two"=6:10)
    res <- quiet(cbind(se1, se2))
    checkIdentical(names(mcols(rowRanges(res))), c("one", "two"))
    rowData(se2) <- DataFrame("one"=6:10, "two"=6:10)
    checkException(cbind(se1, se2), silent=TRUE)
    ## colData
    checkTrue(nrow(colData(res)) == 5)
    ## assays 
    se1 <- rseList[[1]]
    se2 <- se1[,1:2]
    assays(se1) <- SimpleList("m"=matrix(rep("m", 15), nrow=5),
                              "a"=array(rep("a", 30), c(5,3,2)))
    assays(se2) <- SimpleList("m"=matrix(LETTERS[1:10], nrow=5),
                              "a"=array(LETTERS[1:20], c(5,2,2)))
    res <- cbind(se1, se2) ## same variables
    checkTrue(nrow(res) == 5)
    checkTrue(ncol(res) == 5)
    checkTrue(all.equal(dim(assays(res)$m), c(5L, 5L)))
    checkTrue(all.equal(dim(assays(res)$a), c(5L, 5L, 2L)))
    names(assays(se1)) <- c("mm", "aa")
    checkException(cbind(se1, se2), silent=TRUE) ## different variables
}

test_RangedSummarizedExperiment_rbind <- function()
## requires matching samples 
{
    ## empty
    se <- SummarizedExperiment()
    empty <- rbind(se, se)
    checkTrue(all.equal(se, empty))

    ## different samples 
    se1 <- rseList[[1]]
    se2 <- se1[,1]
    checkException(quiet(rbind(se1, se2)), silent=TRUE)

    ## same samples 
    se1 <- rseList[[1]]
    se2 <- se1
    rownames(se2) <- LETTERS[seq_len(nrow(se2))]
    res <- rbind(se1, se2)
    checkTrue(nrow(res) == 10)
    checkTrue(ncol(res) == 3)
    ## rowRanges
    rowData(se1) <- DataFrame("one"=1:5)
    rowData(se2) <- DataFrame("two"=6:10)
    checkIdentical(
        rbind(
            cbind(rowData(se1), two = NA_integer_),
            cbind(one = NA_integer_, rowData(se2))
        ),
        rowData(rbind(se1, se2), use.names = FALSE)
    )
    ## colDat
    se1 <- rseList[[1]]
    se2 <- se1
    colData(se2) <- DataFrame("one"=1:3, "two"=4:6)
    res <- quiet(rbind(se1, se2))
    checkTrue(ncol(colData(res)) == 3)
    ## assays 
    se1 <- rseList[[1]]
    se2 <- se1
    assays(se1) <- SimpleList("m"=matrix(rep("m", 15), nrow=5),
                              "a"=array(rep("a", 30), c(5,3,2)))
    assays(se2) <- SimpleList("m"=matrix(LETTERS[1:15], nrow=5),
                              "a"=array(LETTERS[1:30], c(5,3,2)))
    res <- rbind(se1, se2) ## same variables
    checkTrue(nrow(res) == 10)
    checkTrue(ncol(res) == 3)
    checkTrue(all.equal(dim(assays(res)$m), c(10L, 3L)))
    checkTrue(all.equal(dim(assays(res)$a), c(10L, 3L, 2L)))
    names(assays(se1)) <- c("mm", "aa")
    checkException(rbind(se1, se2), silent=TRUE) ## different variables
}

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

    ## FIXME: pcompare, Compare

    .sig <- "RangedSummarizedExperiment"
    for (.fun in .otherFuns) {
        generic <- getGeneric(.fun)
        method <- getMethod(.fun, "RangedSummarizedExperiment")
        checkIdentical(formals(generic@.Data), formals(method@.Data))
    }        
}

test_RangedSummarizedExperiment_GRanges_values <- function()
{
    x <- rseList[[1]]
    isAssign <- grep("<-$", .singleDispatch, value=TRUE)
    .funs <- setdiff(.singleDispatch, isAssign)
    ## 'exp' created after manual inspection of results
    exp <- setNames(c("02dde", "80339", "df2bc", "410ea", "77198",
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
    se <- SummarizedExperiment(M1, rowRanges=gr, colData=colData)
    ## FIXME: unname should not be necessary
    obs <- split(se, seqnames(se))
    exp <- SimpleList(A=se[1:2], B=se[3:5])
    checkEquals(obs, exp)
}

test_RangedSummarizedExperiment_NULL_rowRanges <- function()
{
    se <- SummarizedExperiment(M1, colData=colData)
    rse <- rseList[[1L]]
    rowRanges(rse) <- NULL
    checkTrue(identical(rowRanges(rse), NULL))
    checkTrue(is(rse, "SummarizedExperiment") &&
        !is(rse, "RangedSummarizedExperiment"))
    checkTrue(identical(rowRanges(se), NULL))
}
