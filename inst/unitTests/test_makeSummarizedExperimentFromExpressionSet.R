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


test_SummarizedExperiment_GenomicRanges_coercion <- function()
{
        eset1 <- ExpressionSet()

        checkTrue(validObject(eset1))

        se1 <- as(eset1, "RangedSummarizedExperiment")

        checkTrue(validObject(se1))

        data("sample.ExpressionSet", package = "Biobase")

        eset2 <- sample.ExpressionSet
        checkTrue(validObject(eset2))

        se2 <- as(eset2, "RangedSummarizedExperiment")

        checkTrue(validObject(se2))

        checkIdentical(experimentData(eset2),
                       metadata(se2)$experimentData)

        checkIdentical(annotation(eset2),
                       metadata(se2)$annotation)

        checkIdentical(protocolData(eset2),
                       metadata(se2)$protocolData)

        eset2Assays <- SimpleList(as.list(assayData(eset2)))
        se2Assays <- assays(se2)
        checkIdentical(eset2Assays$exprs, se2Assays$exprs)
        checkIdentical(eset2Assays$se.exprs, se2Assays$se.exprs)

        checkIdentical(featureNames(eset2), rownames(se2))

        checkIdentical(sampleNames(eset2), colnames(se2))
}

test_GenomicRanges_SummarizedExperiment_coercion <- function()
{
        ## empty SE
        simpleSE <- SummarizedExperiment()

        eset1 <- as(simpleSE, "ExpressionSet")

        checkTrue(validObject(eset1))

        ## Back and forth empty ES
        simpleES <- ExpressionSet()

        simpleES2 <- as(as(simpleES, "RangedSummarizedExperiment"),
                        "ExpressionSet")

        checkTrue(validObject(simpleES2))

        checkEquals(as.list(assayData(simpleES)),
                    as.list(assayData(simpleES2)))

        ## Simple SE
        simpleSE <- rseList[[1]]
        assayNames(simpleSE) <- "exprs" # No warning 'No assay named exprs..."
        eset2 <- as(simpleSE, "ExpressionSet")
        checkTrue(validObject(eset2))

        ## The ExpressionSet features should have the data from the
        ## SummarizedExperiment rows if they are from GRanges.
        checkIdentical(pData(featureData(eset2)),
                       as.data.frame(rowRanges(rseList[[1]])))

        # the rowRanges are retained if the object has them to begin with.
        se2_2 <- as(eset2, "RangedSummarizedExperiment")
        rr_se2_2 <- unname(rowRanges(se2_2))
        rr_eset2 <- rowRanges(rseList[[1]])
        checkEquals(rr_se2_2, rr_eset2)

        simpleSE <- rseList[[2]]
        assayNames(simpleSE) <- "exprs" # No warning 'No assay named exprs..."
        eset3 <- as(simpleSE, "ExpressionSet")
        checkTrue(validObject(eset3))

        ## The ExpressionSet features should not have the data from the
        ## SummarizedExperiment rows if they are from GRangesList, but they
        ## should be empty and the same length as the number of ranges.
        checkEquals(unname(NROW(featureData(eset3))),
                    unname(length(rowRanges(rseList[[2]]))))

        data("sample.ExpressionSet", package = "Biobase")
        eset4 <- sample.ExpressionSet

        eset5 <- as(as(eset4, "RangedSummarizedExperiment"), "ExpressionSet")

        checkTrue(validObject(eset5))

        ## this is necessary because the order in environments is undefined.
        compareLists <- function(x, y) {
            nmsX <- names(x)
            nmsY <- names(y)

            reorderY <- match(nmsY, nmsX)

            checkIdentical(x, y[reorderY])
        }

        compareLists(as.list(assayData(eset4)),
                       as.list(assayData(eset5)))

        checkIdentical(experimentData(eset4),
                       experimentData(eset5))

        checkIdentical(annotation(eset4),
                       annotation(eset5))

        checkIdentical(protocolData(eset4),
                       protocolData(eset5))

        checkIdentical(featureNames(eset4),
                       featureNames(eset5))

        checkIdentical(sampleNames(eset4),
                       sampleNames(eset5))
}

test_GenomicRanges_SummarizedExperiment_coercion_lockedEnvironment <- function()
{
    ## https://github.com/Bioconductor/SummarizedExperiment/issues/43
    se = SummarizedExperiment(list(exprs = matrix(1:10, 5)))
    es1 = es2 = as(se, "ExpressionSet")
    original <- exprs(es2)
    checkIdentical(original, exprs(es2))
    exprs(es1)[1, 1] = 2
    checkTrue(!identical(original, exprs(es1)))
    checkIdentical(original, exprs(es2))
}
    
test_GenomicRanges_SummarizedExperiment_coercion_mappingFunctions <- function()
{
    ## naiveRangeMapper
    ## valid object from empty object
    checkTrue(validObject(makeSummarizedExperimentFromExpressionSet(ExpressionSet())))

    ## valid object from sample ExpressionSet
    data("sample.ExpressionSet", package = "Biobase")
    eset1 <- sample.ExpressionSet
    checkTrue(validObject(makeSummarizedExperimentFromExpressionSet(eset1)))

    ## makeSummarizedExperimentFromExpressionSet should be the same as `as`
    ## with default args
    checkEquals(makeSummarizedExperimentFromExpressionSet(eset1),
                as(eset1, "RangedSummarizedExperiment"))

    ## probeRangeMapper
    ## valid object from empty object
    checkTrue(validObject(
            makeSummarizedExperimentFromExpressionSet(ExpressionSet(),
                probeRangeMapper)))

    ## valid object from sample ExpressionSet
    se1 <- makeSummarizedExperimentFromExpressionSet(eset1, probeRangeMapper)
    checkTrue(validObject(se1))

    ## Granges returned have rownames that were from the featureNames
    checkTrue(all(rownames(rowRanges(se1)) %in% featureNames(eset1)))

    ## geneRangeMapper
    ## valid object from empty object
    checkTrue(validObject(
            makeSummarizedExperimentFromExpressionSet(ExpressionSet(),
                geneRangeMapper(NULL))))

    ## valid object from sample ExpressionSet
    se2 <- makeSummarizedExperimentFromExpressionSet(eset1,
        geneRangeMapper("TxDb.Hsapiens.UCSC.hg19.knownGene"))
    checkTrue(validObject(se2))

    ## Granges returned have rownames that were from the featureNames
    checkTrue(all(rownames(rowRanges(se2)) %in% featureNames(eset1)))
}

