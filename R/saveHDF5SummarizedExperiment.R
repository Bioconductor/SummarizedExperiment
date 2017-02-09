### =========================================================================
### Save/load a HDF5-based SummarizedExperiment object
### -------------------------------------------------------------------------


.write_h5_assays <- function(assays, h5_file, verbose)
{
    nassay <- length(assays)
    for (i in seq_len(nassay)) {
        a <- assays[[i]]
        h5_name <- sprintf("assay%03d", i)
        if (verbose)
            message("Start writing assay ", i, "/", nassay, " to '",
                    h5_file, "':")
        a <- HDF5Array::writeHDF5Array(a, h5_file, h5_name, verbose=verbose)
        if (verbose)
            message("Finished writing assay ", i, "/", nassay, " to '",
                    h5_file, "'.")
        assays[[i]] <- a
    }
    assays
}

.shorten_assay_h5_paths <- function(assays)
{
    nassay <- length(assays)
    for (i in seq_len(nassay)) {
        a <- assays[[i]]
        a@seed@file <- basename(a@seed@file)
        assays[[i]] <- a
    }
    assays
}

### Save all the assays in HDF5 format, including in-memory assays.
### Delayed assays with delayed operations on them are realized while they
### are written to disk..
saveHDF5SummarizedExperiment <- function(x, dir="my_h5_se", verbose=FALSE)
{
    if (!is(x, "SummarizedExperiment"))
        stop("'x' must be a SummarizedExperiment object")
    if (!isTRUEorFALSE(verbose))
        stop("'verbose' must be TRUE or FALSE")
    if (!isSingleString(dir))
        stop(wmsg("'dir' must be a single string specifying the path ",
                  "to the directory where to save the ", class(x),
                  " object (the directory will be created)"))

    ## Try library(HDF5Array) before creating directory 'dir'. That way if
    ## HDF5Array is not installed then saveHDF5SummarizedExperiment() will
    ## stop without having created an empty directory.
    library(HDF5Array)  # for writeHDF5Array()
    if (!suppressWarnings(dir.create(dir)))
        stop("cannot create dir \"", dir, "\"")

    h5_file <- file.path(dir, "assays.h5")
    x@assays <- .write_h5_assays(x@assays, h5_file, verbose)

    rds_file <- file.path(dir, "se.rds")
    ans <- x
    x@assays <- .shorten_assay_h5_paths(x@assays)
    saveRDS(x, file=rds_file)

    invisible(ans)
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
loadHDF5SummarizedExperiment <- function(dir="my_h5_se")
{
    library(rhdf5)  # for h5ls()
    library(HDF5Array)  # for the HDF5Array class
    if (!isSingleString(dir))
        stop(wmsg("'dir' must be a single string specifying the path ",
                  "to the directory containing ", .THE_EXPECTED_STUFF))
    h5_file <- file.path(dir, "assays.h5")
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
        if (!is(a, "HDF5Array") || !identical(a@seed@file, "assays.h5") ||
            !(a@seed@name %in% h5_datasets))
            .stop_if_bad_dir(dir)
        a@seed@file <- file.path(dir, a@seed@file)
        assay(ans, i) <- a
    }
    ans
}

