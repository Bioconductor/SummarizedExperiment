### =========================================================================
### Save/load a HDF5-based SummarizedExperiment object
### -------------------------------------------------------------------------


### Save all the assays in HDF5 format, including in-memory assays.
### Delayed assays with delayed operations on them are realized while they
### are written to disk..
saveHDF5SummarizedExperiment <- function(x, dir="myse", verbose=FALSE)
{
    if (!is(x, "SummarizedExperiment"))
        stop("'x' must be a SummarizedExperiment object")
    if (!isTRUEorFALSE(verbose))
        stop("'verbose' must be TRUE or FALSE")
    library(HDF5Array)  # for writeHDF5Array()
    if (!isSingleString(dir))
        stop(wmsg("'dir' must be a single string specifying the path ",
                  "to the directory where to save the ", class(x),
                  " object (the directory will be created)"))
    if (!suppressWarnings(dir.create(dir)))
        stop("cannot create dir \"", dir, "\"")
    h5_file <- file.path(dir, "se.h5")
    nassay <- length(assays(x))
    for (i in seq_len(nassay)) {
        a <- assay(x, i, withDimnames=FALSE)
        h5_name <- sprintf("assay%03d", i)
        if (verbose)
            message("Start writing assay ", i, "/", nassay, " to '",
                    h5_file, "':")
        a <- HDF5Array::writeHDF5Array(a, h5_file, h5_name, verbose=verbose)
        if (verbose)
            message("Finished writing assay ", i, "/", nassay, " to '",
                    h5_file, "'.")
        a@seed@file <- basename(a@seed@file)
        assay(x, i) <- a
    }
    rds_file <- file.path(dir, "se.rds")
    saveRDS(x, file=rds_file)
    invisible(x)
}

.THE_EXPECTED_STUFF <- c(
    "a HDF5-based SummarizedExperiment object previously ",
    "saved with saveHDF5SummarizedExperiment()"
)

.stop_if_bad_dir <- function(dir)
    stop(wmsg("'", dir, "' dir does not seem to contain ",
              .THE_EXPECTED_STUFF))

### Does a lot of checking and tries to fail graciously if the content
### of 'dir' doesn't look as expected.
loadHDF5SummarizedExperiment <- function(dir="myse")
{
    library(rhdf5)  # for h5ls()
    library(HDF5Array)  # for the HDF5Array class
    if (!isSingleString(dir))
        stop(wmsg("'dir' must be a single string specifying the path ",
                  "to the directory containing ", .THE_EXPECTED_STUFF))
    h5_file <- file.path(dir, "se.h5")
    rds_file <- file.path(dir, "se.rds")
    if (!file.exists(h5_file) || !file.exists(rds_file))
        .stop_if_bad_dir(dir)
    h5_content <- try(rhdf5::h5ls(h5_file), silent=TRUE)
    if (inherits(h5_content, "try-error"))
        .stop_if_bad_dir(dir)
    h5_datasets <- h5_content[ , "name"]
    ans <- readRDS(rds_file)
    if (!is(ans, "SummarizedExperiment"))
        .stop_if_bad_dir(dir)
    for (i in seq_along(assays(ans))) {
        a <- assay(ans, i, withDimnames=FALSE)
        if (!is(a, "HDF5Array") || !identical(a@seed@file, "se.h5") ||
            !(a@seed@name %in% h5_datasets))
            .stop_if_bad_dir(dir)
        a@seed@file <- file.path(dir, a@seed@file)
        assay(ans, i) <- a
    }
    ans
}

