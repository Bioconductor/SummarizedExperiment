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

test_combineRows_ranges_named <- function() {
    se <- SummarizedExperiment(list(counts=matrix(rpois(1000, 10), ncol=10)))
    se2 <- SummarizedExperiment(list(counts=matrix(rpois(1000, 10), ncol=10)))
    rownames(se) <- paste0("GENE_", 1:100)
    rownames(se2) <- paste0("SPIKE_", 1:100)

    # Returns a vanilla SE.
    out <- combineRows(se, se2, use.names=FALSE)
    checkIdentical(as.character(class(out)), "SummarizedExperiment")
    checkIdentical(rownames(out), c(rownames(se), rownames(se2)))

    # Returns a GRanges.
    replace <- GRanges("chrA", IRanges(1, 1:100))
    names(replace) <- rownames(se)
    rowRanges(se) <- replace

    replace2 <- GRanges("chrB", IRanges(1, 1:100))
    names(replace2) <- rownames(se2)
    rowRanges(se2) <- replace2

    suppressWarnings(out <- combineRows(se, se2, use.names=FALSE))
    checkIdentical(rowRanges(out), suppressWarnings(c(replace, replace2)))

    # Testing different objects.
    se3 <- se2
    rowRanges(se3) <- NULL
    rownames(se3) <- rownames(se2)
    suppressWarnings(out <- combineRows(se, se3, use.names=FALSE))
    checkTrue(is(rowRanges(out), "GRangesList"))
    checkIdentical(unname(lengths(rowRanges(out))), rep(c(1L, 0L), c(nrow(se), nrow(se3))))

    se4 <- se2
    rowRanges(se4) <- as(rowRanges(se4), "GRangesList")
    suppressWarnings(out <- combineRows(se, se4, use.names=FALSE))
    expected <- suppressWarnings(as(c(replace, replace2), "GRangesList"))
    checkIdentical(rowRanges(out), expected)

    # Order doesn't affect conversion to GRL.
    suppressWarnings(out <- combineRows(se4, se, use.names=FALSE))
    expected <- suppressWarnings(as(c(replace2, replace), "GRangesList"))
    checkIdentical(rowRanges(out), expected)

    suppressWarnings(combined <- rowRanges(combineRows(se, se3, se4, use.names=FALSE)))
    checkIdentical(unname(lengths(combined)), rep(c(1L, 0L, 1L), c(nrow(se), nrow(se3), nrow(se4))))
    suppressWarnings(combined <- rowRanges(combineRows(se3, se, se4, use.names=FALSE)))
    checkIdentical(unname(lengths(combined)), rep(c(0L, 1L, 1L), c(nrow(se3), nrow(se), nrow(se4))))
}

test_combineRows_ranges_unnamed <- function() {
    # Repeating the same suite of tests for SEs without rownames.
    # This checks the correctness of some edge-case behaviors.

    se <- SummarizedExperiment(list(counts=matrix(rpois(1000, 10), ncol=10)))
    se2 <- SummarizedExperiment(list(counts=matrix(rpois(1000, 10), ncol=10)))

    # Returns a vanilla SE.
    out <- combineRows(se, se2, use.names=FALSE)
    checkIdentical(as.character(class(out)), "SummarizedExperiment")
    checkIdentical(nrow(out), nrow(se) + nrow(se2))

    # Returns a GRanges.
    rowRanges(se) <- GRanges("chrA", IRanges(1, 1:100))
    rowRanges(se2) <- GRanges("chrB", IRanges(1, 1:100))
    suppressWarnings(out <- combineRows(se, se2, use.names=FALSE))
    checkIdentical(rowRanges(out), suppressWarnings(c(rowRanges(se), rowRanges(se2))))

    # Testing different objects.
    se3 <- se2
    rowRanges(se3) <- NULL
    suppressWarnings(out <- combineRows(se, se3, use.names=FALSE))
    checkTrue(is(rowRanges(out), "GRangesList"))
    checkIdentical(unname(lengths(rowRanges(out))), rep(c(1L, 0L), c(nrow(se), nrow(se3))))

    se4 <- se2
    rowRanges(se4) <- as(rowRanges(se4), "GRangesList")
    suppressWarnings(out <- combineRows(se, se4, use.names=FALSE))
    expected <- suppressWarnings(as(c(rowRanges(se), rowRanges(se2)), "GRangesList"))
    checkIdentical(rowRanges(out), expected)

    suppressWarnings(combined <- rowRanges(combineRows(se, se3, se4, use.names=FALSE)))
    checkIdentical(unname(lengths(combined)), rep(c(1L, 0L, 1L), c(nrow(se), nrow(se3), nrow(se4))))

    # Handles partial row names.
    rownames(se) <- paste0("GENE_", 1:100)
    suppressWarnings(out <- combineRows(se, se2, use.names=FALSE))
    checkIdentical(rownames(out), c(rownames(se), character(nrow(se2))))
}

test_combineCols_unnamed <- function() {
    se <- SummarizedExperiment(list(counts=matrix(rpois(1000, 10), ncol=10)))
    colData(se)$A <- 1L
    rowData(se)$A <- 1

    se2 <- SummarizedExperiment(list(counts=matrix(rpois(1000, 10), ncol=10), 
        normalized=matrix(rnorm(1000), ncol=10)))
    colData(se2)$A <- 2L
    colData(se2)$B <- 3
    rowData(se2)$B <- "B"

    stuff <- combineCols(se, se2, use.names=FALSE)

    # Column data is correctly combined.
    checkIdentical(stuff$A, rep(1:2, each=10))
    checkIdentical(stuff$B, rep(c(NA, 3), each=10))

    # Row data is correctly combined.
    checkIdentical(rowData(stuff)$A, rep(1, nrow(se)))
    checkIdentical(rowData(stuff)$B, rep("B", nrow(se)))

    # Assay data is correctly combined.
    checkIdentical(as.matrix(assay(stuff)), cbind(assay(se), assay(se2)))
    checkIdentical(as.matrix(assay(stuff, 2)), cbind(matrix(NA, nrow(se), ncol(se)), assay(se2, 2)))

    # Unary methods work as expected.
    checkIdentical(se, combineCols(se, delayed=FALSE, use.names=FALSE))
}

test_combineCols_named <- function() {
    se <- SummarizedExperiment(list(counts=matrix(rpois(1000, 10), ncol=100)))
    colData(se)$A <- 1L
    rowData(se)$A <- 1

    se2 <- SummarizedExperiment(list(counts=matrix(rpois(1000, 10), ncol=50), 
        normalized=matrix(rnorm(1000), ncol=50)))
    colData(se2)$A <- 2L
    colData(se2)$B <- 3
    rowData(se2)$B <- "B"

    # This fails, because we expect matching numbers of columns when use.names=TRUE.
    checkException(combineCols(se, se2), silent=TRUE)

    rownames(se) <- letters[1:10]
    rownames(se2) <- letters[3:22]
    stuff <- combineCols(se, se2)

    # Column data is correctly combined
    checkIdentical(rownames(stuff), letters[1:22])
    checkIdentical(stuff$A, rep(1:2, c(ncol(se), ncol(se2))))
    checkIdentical(stuff$B, rep(c(NA, 3), c(ncol(se), ncol(se2))))

    # Row data is correctly combined.
    checkIdentical(rowData(stuff)$A, rep(c(1, NA), c(nrow(se), 12)))
    checkIdentical(rowData(stuff)$B, rep(c(NA, "B"), c(2, nrow(se2))))

    # Assay data is correctly combined.
    mat <- as.matrix(assay(stuff))
    ref <- cbind(
        rbind(assay(se), matrix(NA, 12, ncol(se))),
        rbind(NA, NA, assay(se2))
    )
    rownames(ref) <- letters[1:22]
    checkIdentical(mat, ref)

    mat <- as.matrix(assay(stuff, 2))
    ref <- cbind(
        matrix(NA, nrow(stuff), ncol(se)),
        rbind(NA, NA, assay(se2, 2))
    )
    rownames(ref) <- letters[1:22]
    checkIdentical(mat, ref)

    # Unary methods work as expected.
    checkIdentical(se, combineCols(se, delayed=FALSE))
}

test_combineCols_assays <- function() {
    # Deep dive into correct assay name behavior.
    se <- SummarizedExperiment(list(matrix(rpois(1000, 10), ncol=10)))
    se2 <- SummarizedExperiment(list(matrix(rpois(1000, 10), ncol=10), 
        matrix(rnorm(1000), ncol=10)))
    colnames(se) <- letters[1:10]
    colnames(se2) <- letters[15:24]
    rownames(se) <- paste0("GENE_", 1:100)
    rownames(se2) <- paste0("SPIKE_", 1:100)

    # This should fail due to differences in the number of _unnamed_ assays.
    checkException(combineCols(se, se2), silent=TRUE)

    # Either all assays are named, or all are unnamed.
    assays(se) <- assays(se)[c(1, 1)]
    assayNames(se2) <- c("WHEE", "BLAH")
    checkException(combineCols(se, se2), silent=TRUE)

    assays(se2) <- unname(assays(se2))
    out <- combineCols(se, se2)
    checkIdentical(colnames(out), letters[c(1:10, 15:24)])
    checkIdentical(rownames(out), c(rownames(se), rownames(se2)))
}

test_combineCols_ranges_named <- function() {
    se <- SummarizedExperiment(list(counts=matrix(rpois(1000, 10), ncol=10)))
    se2 <- SummarizedExperiment(list(counts=matrix(rpois(1000, 10), ncol=10)))
    rownames(se) <- paste0("GENE_", 1:100)
    rownames(se2) <- paste0("GENE_", 21:120)

    # Checking that an SE is returned.
    out <- combineCols(se, se2, use.names=FALSE)
    checkIdentical(as.character(class(out)), "SummarizedExperiment")
    checkIdentical(rownames(out), rownames(se)) # ignoring other row names when use.names=FALSE.

    out <- combineCols(se, se2)
    checkIdentical(as.character(class(out)), "SummarizedExperiment")
    checkIdentical(rownames(out), union(rownames(se), rownames(se2))) 

    # Checking that an RSE is returned.
    ref <- GRanges("chrA", IRanges(1, 1:120), seqinfo=Seqinfo(seqlengths=c(chrA=1000)))
    names(ref) <- paste0("GENE_", 1:120)
    rowRanges(se) <- ref[1:100]
    rowRanges(se2) <- ref[21:120]

    suppressWarnings(out <- combineCols(se, se2, use.names=FALSE)) # should have a warning here due to differences in values.
    checkIdentical(as.character(class(out)), "RangedSummarizedExperiment")
    checkIdentical(rowRanges(out), rowRanges(se)) 

    out <- combineCols(se, se2)
    checkIdentical(rowRanges(out), ref)

    # Checking that it works with mixtures of object classes in rowRanges.
    se3 <- se2
    rowRanges(se3) <- NULL
    rownames(se3) <- rownames(se2)

    out <- combineCols(se, se3)
    checkTrue(is(rowRanges(out), "GRangesList"))
    checkIdentical(rownames(out), paste0("GENE_", 1:120)) 
    checkIdentical(unname(lengths(rowRanges(out))), rep(1:0, c(100, 20))) 

    out2 <- combineCols(se3, se) # flipping the order.
    checkIdentical(rownames(out2), paste0("GENE_", c(21:120, 1:20)))
    checkIdentical(unname(lengths(rowRanges(out2))), rep(c(1L,0L,1L), c(80, 20, 20)))

    out3 <- combineCols(se, se2, se3) 
    checkIdentical(rowRanges(out3), ref) # avoid unnecessary conversion to a GRL.

    se4 <- se2
    rowRanges(se4) <- as(rowRanges(se4), "GRangesList")
    out <- combineCols(se, se4)
    checkIdentical(rowRanges(out), as(ref, "GRangesList"))

    # Checking that we get the same object class, regardless of ordering of inputs.
    checkIdentical(rowRanges(out), rowRanges(combineCols(se, se3, se4)))
    checkIdentical(rowRanges(out), rowRanges(combineCols(se3, se, se4))[rownames(out)]) 
    checkIdentical(rowRanges(out), rowRanges(combineCols(se3, se4, se))[rownames(out)])

    # Handles conflicting features correctly.
    se5 <- se2
    strand(rowRanges(se5)[1]) <- "+"
    suppressWarnings(out <- combineCols(se, se5)) # this should emit a warning.
    checkTrue(is(rowRanges(out), "GRangesList"))
    checkIdentical(unname(lengths(rowRanges(out))), rep(1:0, c(100, 20))) 
}

test_combineCols_ranges_unnamed <- function() {
    # Repeating the same suite of tests for SEs without rownames.
    # This checks the correctness of some edge-case behaviors.

    se <- SummarizedExperiment(list(counts=matrix(rpois(1000, 10), ncol=10)))
    se2 <- SummarizedExperiment(list(counts=matrix(rpois(1000, 10), ncol=10)))

    # Checking that an SE is returned.
    out <- combineCols(se, se2, use.names=FALSE)
    checkIdentical(as.character(class(out)), "SummarizedExperiment")
    checkIdentical(nrow(out), nrow(se)) 
    checkException(combineCols(se, se2), silent=TRUE)

    # Checking that an RSE is returned.
    ref <- GRanges("chrA", IRanges(1, 1:120), seqinfo=Seqinfo(seqlengths=c(chrA=1000)))
    rowRanges(se) <- ref[1:100]
    rowRanges(se2) <- ref[21:120]

    suppressWarnings(out <- combineCols(se, se2, use.names=FALSE)) # should have a warning here due to differences in values.
    checkIdentical(as.character(class(out)), "RangedSummarizedExperiment")
    checkIdentical(rowRanges(out), rowRanges(se)) 

    # Checking that mixtures of objects work.
    se3 <- se2
    rowRanges(se3) <- NULL

    out <- combineCols(se, se3, use.names=FALSE) # no warning.
    checkIdentical(rowRanges(out), rowRanges(se)) 
    out2 <- combineCols(se3, se, use.names=FALSE)
    checkIdentical(rowRanges(out2), rowRanges(se)) 

    se4 <- se2
    rowRanges(se4) <- as(rowRanges(se4), "GRangesList")
    suppressWarnings(out <- combineCols(se, se4, use.names=FALSE)) # has warning.
    checkIdentical(rowRanges(out), as(rowRanges(se), "GRangesList"))

    se5 <- se
    rowRanges(se5) <- as(rowRanges(se5), "GRangesList")
    out <- combineCols(se, se5, use.names=FALSE) # no warning.
    checkIdentical(rowRanges(out), as(rowRanges(se), "GRangesList"))

    multi.com <- suppressWarnings(combineCols(se, se3, se4, use.names=FALSE))
    checkIdentical(rowRanges(out), rowRanges(multi.com))
}
