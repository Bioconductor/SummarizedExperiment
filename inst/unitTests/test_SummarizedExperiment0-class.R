M1 <- matrix(1, 5, 3, dimnames=list(NULL, NULL))
M2 <- matrix(1, 3, 3, dimnames=list(NULL, NULL))
mList <- list(M1, M2)
assaysList <- list(M1=SimpleList(m=M1), M2=SimpleList(m=M2))
colData <- DataFrame(x=letters[1:3])

se0List <-
    list(SummarizedExperiment(
           assays=assaysList[["M1"]], 
           colData=colData),
         SummarizedExperiment(
           assays=assaysList[["M2"]], 
           colData=colData))


test_SummarizedExperiment0_construction <- function()
{
    ## empty-ish
    m1 <- matrix(0, 0, 0)
    checkTrue(validObject(new("SummarizedExperiment0")))
    checkTrue(validObject(SummarizedExperiment()))
    checkTrue(validObject(
        SummarizedExperiment(assays=SimpleList(m1))))
    checkException(
        SummarizedExperiment(assays=SimpleList(matrix())),
        "assays dim mismatch", TRUE)
    checkException(
        SummarizedExperiment(assays=SimpleList(m1, matrix())),
        "assays dim mismatch", TRUE)
    checkException(
        SummarizedExperiment(assays=SimpleList(character())),
        "assays class", TRUE)

    ## substance
    for (i in seq_along(se0List)) {
        se0 <- se0List[[i]] 
        checkTrue(validObject(se0))
        checkIdentical(SimpleList(m=mList[[i]]), assays(se0))
        checkIdentical(DataFrame(x=letters[1:3]), colData(se0))
    }

    ## array in assays slot
    ss <- se0List[[1]]
    assays(ss) <- SimpleList(array(1:5, c(5,3,2)))
    checkTrue(validObject(ss))
    checkTrue(all(dim(assays(ss[1:3,1:2])[[1]]) == c(3, 2, 2)))

    ## matrix-of-list in assay slot
    m <- matrix(list(), 2, 3, dimnames=list(LETTERS[1:2], letters[1:3]))
    checkTrue(validObject(se <- SummarizedExperiment(m)))
    checkIdentical(m, assay(se))
    checkIdentical(m[,1:2], assay(se[,1:2]))

    ## DataFrame in assay slot
    df <- DataFrame(a=1:3, b=1:3, row.names=LETTERS[1:3])
    checkTrue(validObject(SummarizedExperiment(list(df))))
}

test_SummarizedExperiment0_getters <- function()
{
    for (i in seq_along(se0List)) {
        se0 <- se0List[[i]] 

        ## dim, dimnames
        checkIdentical(c(nrow(mList[[i]]), nrow(colData)), dim(se0))
        checkIdentical(list(NULL, NULL), dimnames(se0))

        ## col / metadata
        checkIdentical(colData, colData(se0))
        checkIdentical(list(), metadata(se0))
    }

    ## assays
    m0 <- matrix(0L, 0, 0, dimnames=list(NULL, NULL))
    m1 <- matrix(0, 0, 0, dimnames=list(NULL, NULL))
    a <- SimpleList(a=m0, b=m1)
    checkIdentical(a, assays(SummarizedExperiment(assays=a)))
    ## assay
    checkException(
        assay(SummarizedExperiment()), "0-length assay", TRUE)
    checkIdentical(m0,
        assay(SummarizedExperiment(assays=a)), "default assay")
    checkIdentical(m1,
        assay(SummarizedExperiment(assays=a), 2),
        "assay, numeric index")
    checkException(
        assay(SummarizedExperiment(assays=a), 3),
        "invalid assay index", TRUE)
    checkIdentical(m1,
        assay(SummarizedExperiment(assays=a), "b"),
        "assay, character index")
    checkException(
        assay(SummarizedExperiment(assays=a), "c"),
        "invalid assay name", TRUE)
}

test_SummarizedExperiment0_setters <- function()
{
    for (i in seq_along(se0List)) {
        se0 <- se0List[[i]] 

        ## row / col / metadata<-
        ss1 <- se0
        revData <- colData[rev(seq_len(nrow(colData))),,drop=FALSE]
        colData(ss1) <- revData
        checkIdentical(revData, colData(ss1))
        checkException(colData(ss1) <- colData(se0)[1:2,,drop=FALSE],
                       "incorrect col dimensions", TRUE)
        lst <- list("foo", "bar")
        metadata(ss1) <- lst
        checkIdentical(lst, metadata(ss1))

        ## assay / assays
        ss1 <- se0
        assay(ss1) <- assay(ss1)+1
        checkIdentical(assay(se0)+1, assay(ss1))
        ss1 <- se0
        assay(ss1, 1) <- assay(ss1, 1) + 1
        checkIdentical(assay(se0, "m") + 1, assay(ss1, "m"))
        ss1 <- se0
        assay(ss1, "m") <- assay(ss1, "m") + 1
        checkIdentical(assay(se0, "m")+1, assay(ss1, "m"))

        ## dimnames<-
        ss1 <- se0
        dimnames <- list(letters[seq_len(nrow(ss1))],
                         LETTERS[seq_len(ncol(ss1))])
        rownames(ss1) <- dimnames[[1]]
        colnames(ss1) <- dimnames[[2]]
        checkIdentical(dimnames, dimnames(ss1))
        colData1 <- colData
        row.names(colData1) <- dimnames[[2]]
        checkIdentical(colData1, colData(ss1))
        ss1 <- se0
        dimnames(ss1) <- dimnames
        checkIdentical(dimnames, dimnames(ss1))
        dimnames(ss1) <- NULL
        checkIdentical(list(NULL, NULL), dimnames(ss1))
    }
}

test_SummarizedExperiment0_subset <- function()
{
    for (i in seq_along(se0List)) {
        se0 <- se0List[[i]] 

        ## numeric
        ss1 <- se0[2:3,]
        checkIdentical(c(2L, ncol(se0)), dim(ss1))
        checkIdentical(colData(ss1), colData(se0))
        ss1 <- se0[,2:3]
        checkIdentical(c(nrow(se0), 2L), dim(ss1))
        checkIdentical(colData(ss1), colData(se0)[2:3,,drop=FALSE])
        ss1 <- se0[2:3, 2:3]
        checkIdentical(c(2L, 2L), dim(ss1))
        checkIdentical(colData(ss1), colData(se0)[2:3,,drop=FALSE])

        ## character
        ss1 <- se0
        dimnames(ss1) <- list(LETTERS[seq_len(nrow(ss1))],
                               letters[seq_len(ncol(ss1))])
        ridx <- c("B", "C")
        checkException(ss1[LETTERS,], "i-index out of bounds", TRUE)
        cidx <- c("b", "c")
        checkIdentical(colData(ss1[,cidx]), colData(ss1)[cidx,,drop=FALSE])
        checkIdentical(colData(ss1[,"a"]), colData(ss1)["a",,drop=FALSE])
        checkException(ss1[,letters], "j-index out of bounds", TRUE)

        ## logical
        ss1 <- se0
        dimnames(ss1) <- list(LETTERS[seq_len(nrow(ss1))],
                               letters[seq_len(ncol(ss1))])
        checkEquals(ss1, ss1[TRUE,])
        checkIdentical(c(0L, ncol(ss1)), dim(ss1[FALSE,]))
        checkEquals(ss1, ss1[,TRUE])
        checkIdentical(c(nrow(ss1), 0L), dim(ss1[,FALSE]))
        idx <- c(TRUE, FALSE)               # recycling
        ss2 <- ss1[idx,]
        ss2 <- ss1[,idx]
        checkIdentical(colData(ss1)[idx,,drop=FALSE], colData(ss2))

        ## Rle
        ss1 <- se0
        rle <- rep(c(TRUE, FALSE), each=3, length.out=nrow(ss1))
        checkIdentical(assays(ss1[rle]), assays(ss1[Rle(rle)]))
    }

    ## 0 columns
    se <- SummarizedExperiment(matrix(integer(0), nrow=5))
    checkIdentical(dim(se[1:5, ]), c(5L, 0L)) 
    ## 0 rows 
    se <- SummarizedExperiment(colData=DataFrame(samples=1:10))
    checkIdentical(dim(se[ ,1:5]), c(0L, 5L)) 
}

test_SummarizedExperiment0_subsetassign <- function()
{
    for (i in seq_along(se0List)) {
        se0 <- se0List[[i]] 
        dimnames(se0) <- list(LETTERS[seq_len(nrow(se0))],
                               letters[seq_len(ncol(se0))])
        ## rows
        ss1 <- se0
        ss1[1:2,] <- ss1[2:1,]
        checkIdentical(colData(se0), colData(ss1))
        checkIdentical(c(metadata(se0), metadata(se0)), metadata(ss1))
        ## Rle
        ss1rle <- ss1Rle <- se0
        rle <- rep(c(TRUE, FALSE), each=3, length.out=nrow(ss1))
        ss1rle[rle,] <- ss1rle[rle,]
        ss1Rle[Rle(rle),] <- ss1Rle[Rle(rle),]
        checkIdentical(assays(ss1rle), assays(ss1Rle))
        ## cols
        ss1 <- se0
        ss1[,1:2] <- ss1[,2:1,drop=FALSE]
        checkIdentical(colData(se0)[2:1,,drop=FALSE],
                       colData(ss1)[1:2,,drop=FALSE])
        checkIdentical(colData(se0)[-(1:2),,drop=FALSE],
                       colData(ss1)[-(1:2),,drop=FALSE])
        checkIdentical(c(metadata(se0), metadata(se0)), metadata(ss1))
    }

    ## full replacement
    ss1 <- ss2 <- se0List[[1]]
    ss1[,] <- ss2
    checkIdentical(ss1, ss2)
}

test_SummarizedExperiment0_assays_4d <- function()
{
    ## [
    a <- array(0, c(3, 3, 3, 3), list(LETTERS[1:3], letters[1:3], NULL, NULL))
    assays <- SimpleList(a=a)
    se <- SummarizedExperiment(assays)
    checkIdentical(assays(se[1,])[[1]], a[1,,,,drop=FALSE])

    ## [<-
    a1 <- a; a1[1,,,] <- a[1,,,,drop=FALSE] + 1
    assays(se[1,])[[1]] <- 1 + assays(se[1,])[[1]]
    checkIdentical(assays(se)[[1]], a1)

    ## [, [<- don't support more than 4 dimensions
    a <- array(0, c(3, 3, 3, 3, 3),
               list(LETTERS[1:3], letters[1:3], NULL, NULL, NULL))
    assays <- SimpleList(a=a)
    se <- SummarizedExperiment(assays)
    checkException(se[1,], silent=TRUE)
}

