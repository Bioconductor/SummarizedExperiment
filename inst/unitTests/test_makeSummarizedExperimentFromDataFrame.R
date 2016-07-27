##

rowNames <- paste0("GENE", letters[5:1])

range_info <- list(chr="chr2", start = 11:15, end = 12:16,
                      strand = c("+", "-", "+", "*", "."))
expr_info <- list(expr0 = 3:7, expr1 = 8:12, expr2 = 12:16)

df <- as.data.frame(c(range_info, expr_info), row.names = rowNames)
DF <- DataFrame(c(range_info, expr_info), row.names = rowNames)

test_makeSummarizedExperimentFromDataFrame <- function()
{
    validObject(makeSummarizedExperimentFromDataFrame(df))
    validObject(makeSummarizedExperimentFromDataFrame(DF))

    rangesA <- GRanges(as.data.frame(range_info, row.names = rowNames))
    rangesB <- rowRanges(makeSummarizedExperimentFromDataFrame(df))
    # Check rowRanges to be identical
    checkIdentical(rangesA, rangesB)
    # Check assay matrix and expr_info matrix are identical
    checkIdentical(assay(makeSummarizedExperimentFromDataFrame(df)),
                   as.matrix(as.data.frame(expr_info, row.names = rowNames)))
    checkIdentical(assay(makeSummarizedExperimentFromDataFrame(DF)),
                   as.matrix(as.data.frame(expr_info, row.names = rowNames)))

    checkEquals(makeSummarizedExperimentFromDataFrame(df),
                makeSummarizedExperimentFromDataFrame(DF))

    checkException(
        makeSummarizedExperimentFromDataFrame(
            cbind(df, expr3 = letters[seq_len(nrow(df))])))

    checkException(
        makeSummarizedExperimentFromDataFrame(
            cbind(DF, DataFrame(expr3 = letters[seq_len(nrow(df))]))))

    checkIdentical(nrow(df),
                   length(rowRanges(
                       makeSummarizedExperimentFromDataFrame(df))))

    checkIdentical(nrow(DF),
                   length(rowRanges(
                       makeSummarizedExperimentFromDataFrame(DF))))

    checkIdentical(colnames(makeSummarizedExperimentFromDataFrame(df)),
                   names(expr_info))
    checkIdentical(rownames(makeSummarizedExperimentFromDataFrame(df)),
                   rowNames)

    checkIdentical(colnames(makeSummarizedExperimentFromDataFrame(DF)),
                   names(expr_info))
    checkIdentical(rownames(makeSummarizedExperimentFromDataFrame(DF)),
                   rowNames)
}

