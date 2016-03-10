.TEST_matrices <- list(
    matrix(1:15, nrow=3, ncol=5,
           dimnames=list(NULL, paste0("M1y", 1:5))),
    matrix(101:135, nrow=7, ncol=5,
           dimnames=list(paste0("M2x", 1:7), paste0("M2y", 1:5))),
    matrix(1001:1025, nrow=5, ncol=5,
           dimnames=list(paste0("M3x", 1:5), NULL))
)

.TEST_arrays <- list(
    array(1:60, c(3, 5, 4),
           dimnames=list(NULL, paste0("M1y", 1:5), NULL)),
    array(101:240, c(7, 5, 4),
           dimnames=list(paste0("M2x", 1:7), paste0("M2y", 1:5), NULL)),
    array(10001:10100, c(5, 5, 4),
           dimnames=list(paste0("M3x", 1:5), NULL, paste0("M3z", 1:4)))
)


test_arbind.default <- function()
{
    arbind.default <- SummarizedExperiment:::arbind.default

    ## on matrices
    target <- do.call(rbind, .TEST_matrices)
    current <- do.call(arbind.default, .TEST_matrices)
    checkIdentical(target, current)

    ## on arrays
    current <- do.call(arbind.default, .TEST_arrays)
    for (k in 1:4) {
        target <- do.call(rbind, lapply(.TEST_arrays, `[`, , , k))
        checkIdentical(target, current[ , , k])
    }
}

test_acbind.default <- function()
{
    acbind.default <- SummarizedExperiment:::acbind.default

    ## on matrices
    matrices <- lapply(.TEST_matrices, t)
    target <- do.call(cbind, matrices)
    current <- do.call(acbind.default, matrices)
    checkIdentical(target, current)

    ## on arrays

    ## transpose the 1st 2 dimensions
    arrays <- lapply(.TEST_arrays,
        function(a) {
            a_dimnames <- dimnames(a)
            dim(a)[1:2] <- dim(a)[2:1]
            a_dimnames[1:2] <- a_dimnames[2:1]
            dimnames(a) <- a_dimnames
            a
    })
    current <- do.call(acbind.default, arrays)
    for (k in 1:4) {
        target <- do.call(cbind, lapply(arrays, `[`, , , k))
        checkIdentical(target, current[ , , k])
    }
}

