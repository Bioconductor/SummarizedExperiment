M1 <- matrix(1, 5, 3)
M2 <- matrix(1, 3, 3)
mList <- list(M1, M2)
assaysList <- list(M1=SimpleList(m=M1), M2=SimpleList(m=M2))
rowData1 <- DataFrame(id1=LETTERS[1:5])
rowData2 <- S4Vectors:::make_zero_col_DataFrame(3L)
rowDataList <- list(rowData1, rowData2)
colData0 <- DataFrame(x=letters[1:3])

se0List <-
    list(SummarizedExperiment(
           assays=assaysList[["M1"]],
           rowData=rowData1,
           colData=colData0),
         SummarizedExperiment(
           assays=assaysList[["M2"]],
           colData=colData0))


test_SummarizedExperiment_construction <- function()
{
    ## empty-ish
    m1 <- matrix(0, 0, 0)
    checkTrue(validObject(new("SummarizedExperiment")))
    checkTrue(validObject(SummarizedExperiment()), "empty constructor")
    checkTrue(validObject(SummarizedExperiment(SimpleList())))
    checkTrue(validObject(SummarizedExperiment(assays=SimpleList(m1))),
              "0x0 constructor")
    checkException(SummarizedExperiment(assays=SimpleList(m1, matrix())),
                   "assays dim mismatch", TRUE)

    ## substance
    for (i in seq_along(se0List)) {
        se0 <- se0List[[i]]
        checkTrue(validObject(se0))
        checkIdentical(SimpleList(m=mList[[i]]), assays(se0))
        checkIdentical(rowDataList[[i]], rowData(se0))
        checkIdentical(colData0, colData(se0))
    }

    ## array in assays slot
    ss <- se0List[[1]]
    assays(ss, withDimnames=FALSE) <- SimpleList(array(1:5, c(5,3,2)))
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

test_SummarizedExperiment_construction_dimnames <- function()
{
    do_tests <- function(m, rowData, colData) {
        ## Some quick tests to make sure that the rowData(), colData(), and
        ## assay() getters handle the rownames and colnames as expected.
        ## Getters are tested more thoroughly in dedicated unit
        ## test_SummarizedExperiment_getters().
        test_rowData_colData_assay <- function(se) {
            checkIdentical(rownames(rowData(se)), rownames(se))
            checkIdentical(rownames(colData(se)), colnames(se))
            a1 <- assay(se, withDimnames=FALSE)
            checkIdentical(dimnames(m), dimnames(a1))
            a1 <- assay(se)
            checkIdentical(rownames(se), rownames(a1))
            checkIdentical(colnames(se), colnames(a1))
        }
        target_rownames <- function() {
            rownames <- rownames(rowData)
            if (is.null(rownames)) rownames(m) else rownames
        }
        target_colnames <- function() {
            colnames <- rownames(colData)
            if (is.null(colnames)) colnames(m) else colnames
        }

        se <- SummarizedExperiment(m)
        checkTrue(validObject(se))
        checkIdentical(rownames(m), rownames(se))
        checkIdentical(colnames(m), colnames(se))
        test_rowData_colData_assay(se)

        se <- SummarizedExperiment(m, rowData=rowData)
        checkTrue(validObject(se))
        checkIdentical(target_rownames(), rownames(se))
        checkIdentical(colnames(m), colnames(se))
        test_rowData_colData_assay(se)

        se <- SummarizedExperiment(m, colData=colData)
        checkTrue(validObject(se))
        checkIdentical(rownames(m), rownames(se))
        checkIdentical(target_colnames(), colnames(se))
        test_rowData_colData_assay(se)

        se <- SummarizedExperiment(m, rowData=rowData, colData=colData)
        checkTrue(validObject(se))
        checkIdentical(target_rownames(), rownames(se))
        checkIdentical(target_colnames(), colnames(se))
        test_rowData_colData_assay(se)
    }

    m <- matrix(0, nrow=4, ncol=3)
    rowData <- DataFrame(stuff=11:14)  # no rownames
    colData <- DataFrame(row.names=letters[1:3])
    do_tests(m, rowData, colData)

    rownames(m) <- paste0("ROW", 1:4)
    do_tests(m, rowData, colData)

    colnames(m) <- letters[1:3]
    do_tests(m, rowData, colData)

    dimnames(m) <- NULL
    rownames(rowData) <- LETTERS[1:4]
    do_tests(m, rowData, colData)

    rownames(m) <- LETTERS[1:4]
    do_tests(m, rowData, colData)

    ## Empty strings in the dimnames are supported.
    dimnames(m) <- list(c("A", "", "", "D"), c("a", "b", ""))
    rowData <- DataFrame(row.names=rownames(m))
    colData <- DataFrame(row.names=colnames(m))
    do_tests(m, rowData, colData)

    ## NAs in the dimnames.
    ## WARNINGS:
    ## - NAs in the **rownames** are tolerated but will cause problems
    ##   downstream e.g. they break the rowData() getter unless
    ##   'use.names=FALSE' is used.
    ## - NAs in the **colnames** are not and cannot be supported at the moment!
    ##   Right now they break the SummarizedExperiment() constructor in an ugly
    ##   way (error message not super helpful):
    ##       > SummarizedExperiment(m)
    ##       Error in DataFrame(x = seq_len(ncol(a1)), row.names = nms) :
    ##         missing values in 'row.names'
    ##   This should be improved.
    ## - At the root of these problems is the fact that at the moment
    ##   DataFrame objects do not support NAs in their rownames (this is
    ##   a of BioC 3.16).
    ## Bottom line: NAs in the dimnames of a SummarizedExperiment object
    ## should be avoided at all cost. One way to deal with them is to
    ## replace them with empty strings (""). It's not clear whether the
    ## SummarizedExperiment() constructor should automatically take care
    ## of that (with a warning), or if it should just fail with an informative
    ## error message asking the user to "fix" the dimnames on the assay
    ## **before** calling the constructor. One possible downside of the latter
    ## is that there is no guarantee that the dimnames of an assay can be
    ## modified in general (think on-disk assay), and, even if they are, doing
    ## so could be costly (e.g. trigger a copy of a huge object). But if the
    ## SummarizedExperiment() constructor were to take care of the fix, it
    ## wouldn't need to modify the dimnames of the supplied assay.
    rownames(m)[c(1L, 3L)] <- NA
    se <- SummarizedExperiment(m)
    checkTrue(validObject(se))
}

test_SummarizedExperiment_getters <- function()
{
    for (i in seq_along(se0List)) {
        se0 <- se0List[[i]]

        ## dim, dimnames
        checkIdentical(c(nrow(mList[[i]]), nrow(colData0)), dim(se0))
        checkIdentical(NULL, dimnames(se0))

        ## col / metadata
        checkIdentical(rowDataList[[i]], rowData(se0))
        checkIdentical(colData0, colData(se0))
        checkIdentical(list(), metadata(se0))
    }

    ## assays
    m0 <- matrix(0L, 0, 0)
    m1 <- matrix(0, 0, 0)
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

test_SummarizedExperiment_setters <- function()
{
    for (i in seq_along(se0List)) {
        se0 <- se0List[[i]]

        ## row / col / metadata<-
        se1 <- se0
        rowData <- rowDataList[[i]]

        rowData <- rowData[rev(seq_len(nrow(rowData))),,drop=FALSE]
        rowData(se1) <- rowData
        checkIdentical(rowData, rowData(se1))

        colData <- colData0[rev(seq_len(nrow(colData0))),,drop=FALSE]
        colData(se1) <- colData
        checkIdentical(colData, colData(se1))

        ## The rowData (alias for mcols) setter recycles the supplied
        ## DataFrame. This is consistent with what the mcols/elementMetadata
        ## setter does on Vector objects in general.
        rowData(se1) <- rowData(se0)[1:2,,drop=FALSE]
        idx <- rep(1:2, length.out=length(se1))
        target_se1_rowData <- rowData(se0)[idx,,drop=FALSE]
        checkIdentical(target_se1_rowData, rowData(se1))

        ## The colData setter does NOT recycle the supplied DataFrame.
        checkException(colData(se1) <- colData(se0)[1:2,,drop=FALSE],
                       "incorrect col dimensions", TRUE)

        lst <- list("foo", "bar")
        metadata(se1) <- lst
        checkIdentical(lst, metadata(se1))

        ## assay / assays
        se1 <- se0
        assay(se1) <- assay(se1)+1
        checkIdentical(assay(se0)+1, assay(se1))
        se1 <- se0
        assay(se1, 1) <- assay(se1, 1) + 1
        checkIdentical(assay(se0, "m") + 1, assay(se1, "m"))
        se1 <- se0
        assay(se1, "m") <- assay(se1, "m") + 1
        checkIdentical(assay(se0, "m")+1, assay(se1, "m"))

        ## dimnames<-
        se1 <- se0
        dimnames <- list(letters[seq_len(nrow(se1))],
                         LETTERS[seq_len(ncol(se1))])
        rownames(se1) <- dimnames[[1]]
        colnames(se1) <- dimnames[[2]]
        checkIdentical(dimnames, dimnames(se1))
        colData1 <- colData0
        row.names(colData1) <- dimnames[[2]]
        checkIdentical(colData1, colData(se1))
        se1 <- se0
        dimnames(se1) <- dimnames
        checkIdentical(dimnames, dimnames(se1))
        dimnames(se1) <- NULL
        checkIdentical(NULL, dimnames(se1))
    }

    ## With empty strings in the dimnames.
    m <- matrix(1:12, nrow=4, dimnames=list(c(a="A", b="B", c="", d="D"),
                                            c(X="x", Y="y", Z="")))
    se <- se0 <- SummarizedExperiment(m)

    assay(se) <- m  # should be a no-op
    checkIdentical(se0, se)

    dimnames(se) <- dimnames(m)  # should be a no-op
    checkIdentical(se0, se)
}

test_SummarizedExperiment_subset <- function()
{
    for (i in seq_along(se0List)) {
        se0 <- se0List[[i]]

        ## numeric
        se1 <- se0[2:3,]
        checkIdentical(c(2L, ncol(se0)), dim(se1))
        checkIdentical(rowData(se1), rowData(se0)[2:3,,drop=FALSE])
        checkIdentical(colData(se1), colData(se0))
        se1 <- se0[,2:3]
        checkIdentical(c(nrow(se0), 2L), dim(se1))
        checkIdentical(rowData(se1), rowData(se0))
        checkIdentical(colData(se1), colData(se0)[2:3,,drop=FALSE])
        se1 <- se0[2:3, 2:3]
        checkIdentical(c(2L, 2L), dim(se1))
        checkIdentical(colData(se1), colData(se0)[2:3,,drop=FALSE])

        ## character
        se1 <- se0
        dimnames(se1) <- list(LETTERS[seq_len(nrow(se1))],
                               letters[seq_len(ncol(se1))])
        ridx <- c("B", "C")
        checkException(se1[LETTERS,], "i-index out of bounds", TRUE)
        cidx <- c("b", "c")
        checkIdentical(colData(se1[,cidx]), colData(se1)[cidx,,drop=FALSE])
        checkIdentical(colData(se1[,"a"]), colData(se1)["a",,drop=FALSE])
        checkException(se1[,letters], "j-index out of bounds", TRUE)

        ## logical
        se1 <- se0
        dimnames(se1) <- list(LETTERS[seq_len(nrow(se1))],
                               letters[seq_len(ncol(se1))])
        checkEquals(se1, se1[TRUE,])
        checkIdentical(c(0L, ncol(se1)), dim(se1[FALSE,]))
        checkEquals(se1, se1[,TRUE])
        checkIdentical(c(nrow(se1), 0L), dim(se1[,FALSE]))
        idx <- c(TRUE, FALSE)               # recycling
        se2 <- se1[idx,]
        se2 <- se1[,idx]
        checkIdentical(colData(se1)[idx,,drop=FALSE], colData(se2))

        ## Rle
        se1 <- se0
        rle <- rep(c(TRUE, FALSE), each=3, length.out=nrow(se1))
        checkIdentical(assays(se1[rle]), assays(se1[Rle(rle)]))
    }

    ## 0 columns
    se <- SummarizedExperiment(matrix(integer(0), nrow=5))
    checkIdentical(dim(se[1:5, ]), c(5L, 0L))
    ## 0 rows
    se <- SummarizedExperiment(colData=DataFrame(samples=1:10))
    checkIdentical(dim(se[ ,1:5]), c(0L, 5L))
}

test_SummarizedExperiment_subsetassign <- function()
{
    for (i in seq_along(se0List)) {
        se0 <- se0List[[i]]
        dimnames(se0) <- list(LETTERS[seq_len(nrow(se0))],
                               letters[seq_len(ncol(se0))])
        ## rows
        se1 <- se0
        se1[1:2,] <- se1[2:1,]
        checkIdentical(colData(se0), colData(se1))
        checkIdentical(c(metadata(se0), metadata(se0)), metadata(se1))
        ## Rle
        se1rle <- se1Rle <- se0
        rle <- rep(c(TRUE, FALSE), each=3, length.out=nrow(se1))
        se1rle[rle,] <- se1rle[rle,]
        se1Rle[Rle(rle),] <- se1Rle[Rle(rle),]
        checkIdentical(assays(se1rle), assays(se1Rle))
        ## cols
        se1 <- se0
        se1[,1:2] <- se1[,2:1,drop=FALSE]
        checkIdentical(colData(se0)[2:1,,drop=FALSE],
                       colData(se1)[1:2,,drop=FALSE])
        checkIdentical(colData(se0)[-(1:2),,drop=FALSE],
                       colData(se1)[-(1:2),,drop=FALSE])
        checkIdentical(c(metadata(se0), metadata(se0)), metadata(se1))
    }

    ## full replacement
    se1 <- se2 <- se0List[[1]]
    se1[,] <- se2
    checkIdentical(se1, se2)
}

test_SummarizedExperiment_assays_4d <- function()
{
    ## construction/validation
    A <- array(0, c(3, 2, 5, 4), list(c("x1", "x2", "x3"),
                                      c("y1", "y2"),
                                      NULL,
                                      c("t1", "t2", "t3", "t4")))
    B <- array(0, c(3, 2, 6),    list(c("x1", "x2", "x3"),
                                      c("y1", "y2"),
                                      NULL))
    assays0 <- SimpleList(A=A, B=B)
    checkTrue(validObject(SummarizedExperiment(assays0)))

    dimnames(B)[[2]] <- c("y1", "oops")
    assays0 <- SimpleList(A=A, B=B)
    checkException(SummarizedExperiment(assays0))

    dimnames(B)[1:2] <- dimnames(A)[1:2]
    C <- array(0, c(3, 2, 4),    list(NULL,
                                      c("y1", "y2"),
                                      c("z1", "z2", "z3", "z4")))
    D <- array(0, c(3, 2, 7, 2), list(NULL,
                                      NULL,
                                      NULL,
                                      c("t1", "t2")))
    E <- array(0, c(3, 2, 0))

    assays0 <- SimpleList(A=A, B=B, C=C, D=D, E=E)
    se0 <- se1 <- SummarizedExperiment(assays0)
    dimnames(se0) <- NULL
    checkTrue(validObject(se0, complete=TRUE))
    checkTrue(validObject(se1, complete=TRUE))

    ## dimnames
    checkIdentical(NULL, dimnames(se0))
    checkIdentical(dimnames(A)[1:2], dimnames(se1))

    ## assays
    for (i in seq_along(assays0)) {
        target <- assays0[[i]]
        checkIdentical(target, assay(se0, i, withDimnames=FALSE))
        checkIdentical(target, assay(se1, i, withDimnames=FALSE))
    }

    target0 <- target1 <- dimnames(assays0[[1]])
    target0[1:2] <- list(NULL)
    checkIdentical(target0, dimnames(assay(se0, 1)))
    checkIdentical(target1, dimnames(assay(se1, 1)))

    checkIdentical(NULL, dimnames(assay(se0, 2)))
    checkIdentical(dimnames(assays0[[2]]), dimnames(assay(se1, 2)))

    target0 <- target1 <- dimnames(assays0[[3]])
    target0[1:2] <- list(NULL)
    checkIdentical(target0, dimnames(assay(se0, 3)))
    target1[1:2] <- dimnames(se1)
    checkIdentical(target1, dimnames(assay(se1, 3)))

    target0 <- target1 <- dimnames(assays0[[4]])
    checkIdentical(target0, dimnames(assay(se0, 4)))
    target1[1:2] <- dimnames(se1)
    checkIdentical(target1, dimnames(assay(se1, 4)))

    checkIdentical(NULL, dimnames(assay(se0, 5)))
    target1 <- c(dimnames(se1), list(NULL))
    checkIdentical(target1, dimnames(assay(se1, 5)))

    ## [
    se2 <- se1[3:2, ]
    checkIdentical(A[3:2, , , , drop=FALSE],
                   assay(se2, 1, withDimnames=FALSE))
    checkIdentical(B[3:2, , , drop=FALSE],
                   assay(se2, 2, withDimnames=FALSE))
    checkIdentical(C[3:2, , , drop=FALSE],
                   assay(se2, 3, withDimnames=FALSE))

    ## [<-
    A1 <- A; A1[1, , , ] <- A[1, , , , drop=FALSE] + 1
    assays(se1[1, ], withDimnames=FALSE)[[1]] <-
        1 + assays(se1[1, ], withDimnames=FALSE)[[1]]
    checkIdentical(A1, assays(se1, withDimnames=FALSE)[[1]])
}

