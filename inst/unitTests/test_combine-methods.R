test_combineRows_unnamed <- function() {
    se <- SummarizedExperiment(list(counts=matrix(rpois(1000, 10), ncol=10)))
    colData(se)$A <- 1
    rowData(se)$A <- 1

    se2 <- SummarizedExperiment(list(counts=matrix(rpois(1000, 10), ncol=10), 
        normalized=matrix(rnorm(1000), ncol=10)))
    colData(se2)$B <- 2
    rowData(se2)$B <- "B"

    stuff <- combineRows(se, se2, use.names=FALSE)

    # Column data is correctly combined.
    checkIdentical(stuff$A, rep(1, ncol(stuff)))
    checkIdentical(stuff$B, rep(2, ncol(stuff)))

    # Row data is correctly combined.
    checkIdentical(rowData(stuff)$A, rep(c(1, NA), c(nrow(se), nrow(se2))))
    checkIdentical(rowData(stuff)$B, rep(c(NA, "B"), c(nrow(se), nrow(se2))))

    # Assay data is correctly combined.
    checkIdentical(as.matrix(assay(stuff)), rbind(assay(se), assay(se2)))
    checkIdentical(as.matrix(assay(stuff, 2)), rbind(matrix(NA, nrow(se), ncol(se)), assay(se2, 2)))

    # Unary methods work as expected.
    checkIdentical(se, combineRows(se, delayed=FALSE, use.names=FALSE))
}

test_combineRows_named <- function() {
    se <- SummarizedExperiment(list(counts=matrix(rpois(1000, 10), ncol=10)))
    colData(se)$A <- 1
    rowData(se)$A <- 1

    se2 <- SummarizedExperiment(list(counts=matrix(rpois(1000, 10), ncol=20), 
        normalized=matrix(rnorm(1000), ncol=20)))
    colData(se2)$B <- 2
    rowData(se2)$B <- "B"

    # This fails, because we expect matching numbers of columns when use.names=TRUE.
    checkException(combineRows(se, se2), silent=TRUE)

    colnames(se) <- letters[1:10]
    colnames(se2) <- letters[3:22]
    stuff <- combineRows(se, se2)

    # Column data is correctly combined
    checkIdentical(colnames(stuff), letters[1:22])
    checkIdentical(stuff$A, rep(c(1, NA), c(ncol(se), 12)))
    checkIdentical(stuff$B, rep(c(NA, 2), c(2, ncol(se2))))

    # Row data is correctly combined.
    checkIdentical(rowData(stuff)$A, rep(c(1, NA), c(nrow(se), nrow(se2))))
    checkIdentical(rowData(stuff)$B, rep(c(NA, "B"), c(nrow(se), nrow(se2))))

    # Assay data is correctly combined.
    mat <- as.matrix(assay(stuff))
    ref <- rbind(
        cbind(assay(se), matrix(NA, nrow(se), ncol=12)),
        cbind(NA, NA, assay(se2))
    )
    colnames(ref) <- letters[1:22]
    checkIdentical(mat, ref)

    mat <- as.matrix(assay(stuff, 2))
    ref <- rbind(
        matrix(NA, nrow(se), ncol(stuff)),
        cbind(NA, NA, assay(se2, 2))
    )
    colnames(ref) <- letters[1:22]
    checkIdentical(mat, ref)

    # Unary methods work as expected.
    checkIdentical(se, combineRows(se, delayed=FALSE))
}

test_combineRows_assays <- function() {
    # Deep dive into correct assay name behavior.
    se <- SummarizedExperiment(list(matrix(rpois(1000, 10), ncol=10)))
    se2 <- SummarizedExperiment(list(matrix(rpois(1000, 10), ncol=10), 
        matrix(rnorm(1000), ncol=10)))
    colnames(se) <- letters[1:10]
    colnames(se2) <- letters[15:24]
    rownames(se) <- paste0("GENE_", 1:100)
    rownames(se2) <- paste0("SPIKE_", 1:100)

    # This should fail due to differences in the number of _unnamed_ assays.
    checkException(combineRows(se, se2), silent=TRUE)

    # Either all assays are named, or all are unnamed.
    assays(se) <- assays(se)[c(1, 1)]
    assayNames(se2) <- c("WHEE", "BLAH")
    checkException(combineRows(se, se2), silent=TRUE)

    assays(se2) <- unname(assays(se2))
    out <- combineRows(se, se2)
    checkIdentical(colnames(out), letters[c(1:10, 15:24)])
    checkIdentical(rownames(out), c(rownames(se), rownames(se2)))
}

test_combineRows_ranges <- function() {
    se <- SummarizedExperiment(list(counts=matrix(rpois(1000, 10), ncol=10)))
    se2 <- SummarizedExperiment(list(counts=matrix(rpois(1000, 10), ncol=10)))
    rownames(se) <- paste0("GENE_", 1:100)
    rownames(se2) <- paste0("SPIKE_", 1:100)

    out <- combineRows(se, se2, use.names=FALSE)
    checkIdentical(as.character(class(out)), "SummarizedExperiment")
    checkIdentical(rownames(out), c(rownames(se), rownames(se2)))

    rowRanges(se) <- GRanges("chrA", IRanges(1, 1:100))
    rowRanges(se2) <- GRanges("chrB", IRanges(1, 1:100))
    suppressWarnings(out <- combineRows(se, se2, use.names=FALSE))
    checkIdentical(rowRanges(out), suppressWarnings(c(rowRanges(se), rowRanges(se2))))

    rowRanges(se2) <- NULL
    suppressWarnings(out <- combineRows(se, se2, use.names=FALSE))
    checkTrue(is(rowRanges(out), "GRangesList"))
    checkIdentical(rownames(out), c(rownames(se), rownames(se2)))

    # Order doesn't matter.
    suppressWarnings(out <- combineRows(se2, se, use.names=FALSE))
    checkTrue(is(rowRanges(out), "GRangesList"))
}
