.import_loom_matrix <-
    function(file, name)
{
    t(HDF5Array::HDF5Array(file, name))
}

.import_loom_DataFrame <-
    function(file, name, rowname)
{
    df <- DataFrame(rhdf5::h5read(file, name))
    df[] <- lapply(df, as.vector)
    if (!is.null(rowname)) {
        rownames(df) <- df[[rowname]]
        df <- df[, -match(rowname, colnames(df)), drop = FALSE]
    }
    if (nrow(df) == 0L)
        df <- NULL
    df
}

makeSummarizedExperimentFromLoom <-
    function(file, rownames_attr = "rownames", colnames_attr = "colnames")
{
    stopifnot(file.exists(file))

    ls <- rhdf5::h5ls(file)
    rowColnames <- ls[ls$group == "/row_attrs", "name", drop=TRUE]
    colColnames <- ls[ls$group == "/col_attrs", "name", drop=TRUE]
    if (missing(rownames_attr) && "rownames" %in% rowColnames)
        rownames_attr <- "rownames"
    if (missing(colnames_attr) && "colnames" %in% colColnames)
        colnames_attr <- "colnames"
    stopifnot(
        is.null(rownames_attr) || rownames_attr %in% rowColnames,
        is.null(colnames_attr) || colnames_attr %in% colColnames
    )

    assay <- .import_loom_matrix(file, "/matrix")
    layerNames <- ls[ls$group == "/layers", "name", drop = TRUE]
    layers <- lapply(setNames(layerNames, layerNames), function(layer) {
        layer <- paste0("/layers/", layer)
        .import_loom_matrix(file, layer)
    })
    assays <- c(list(matrix = assay), layers)

    rowData <- .import_loom_DataFrame(file, "row_attrs", rownames_attr)
    colData <- .import_loom_DataFrame(file, "col_attrs", colnames_attr)

    se <- SummarizedExperiment(assays, rowData = rowData, colData = colData)
    metadata(se) <- rhdf5::h5readAttributes(file, "/")
    se
}

setGeneric(
    ".export_loom",
    function(object, file = tempfile(), ...) standardGeneric(".export_loom"),
    signature = "object"
)

setMethod(".export_loom", "matrix",
    function(object, file, name)
{
    object <- t(object)
    tryCatch({
        rhdf5::h5write(object, file, name)
    }, error = function(err) {
        warning(conditionMessage(err))
        1L
    })
})

setMethod(".export_loom", "DelayedArray",
    function(object, file, name)
{
    HDF5Array::writeHDF5Array(t(object), file, name)
    0L
})

setMethod(".export_loom", "data.frame",
    function(object, file, name, rowname_attr)
{
    if (!is.null(rowname_attr))
        object[[rowname_attr]] <- rownames(object)

    is.factor <- vapply(object, is, logical(1), "factor")
    if (any(is.factor))
        warning(
            "'.export_loom()' coerced 'factor' column(s) to character:\n  ",
            paste(sQuote(names(object)[is.factor]), collapse=", "),
            call. = FALSE
        )
    object[is.factor] <- lapply(object[is.factor], as.character)

    names <- sprintf("/%s/%s", name, names(object))
    tryCatch({
        Map(rhdf5::h5write, object, names, MoreArgs = list(file = file))
    }, error = function(err) {
        warning(conditionMessage(err))
        1L
    })
})

setMethod(".export_loom", "DataFrame",
    function(object, file, name, rowname_attr)
{
    object <- as.data.frame(object)
    .export_loom(object, file, name, rowname_attr)
})

setMethod(".export_loom", "GenomicRanges",
    function(object, file, name, rowname_attr)
{
    object <- as.data.frame(rowRanges(object))
    .export_loom(object, file, name, rowname_attr)
})    

setMethod(".export_loom", "GenomicRangesList",
    function(object, file, name, rowname_attr)
{
    warning(
        "'.export_loom()' does not support '", class(object),
        "'; using mcols()",
        call. = FALSE
    )
    object <- mcols(object, use.names=TRUE)
    .export_loom(object, file, name, rowname_attr)
})    

setMethod(".export_loom", "SummarizedExperiment",
    function(object, file,
             matrix = assayNames(object)[1],
             rownames_attr = "rownames", colnames_attr = "colnames")
{
    stopifnot(
        !file.exists(file),
        is.character(matrix), length(matrix) == 1L, !is.na(matrix),
        matrix %in% assayNames(object),
        is.character(rownames_attr), length(rownames_attr) == 1L,
        !is.na(rownames_attr),
        is.character(colnames_attr), length(colnames_attr) == 1L,
        !is.na(colnames_attr)
    )
    if (!is.null(rownames(object)) && rownames_attr %in% names(rowData(object)))
        stop("'rownames_attr' must not be in names(rowData())")
    if (!is.null(colnames(object)) && colnames_attr %in% names(colData(object)))
        stop("'colnames_attr()' must not be in names(colData())")

    rhdf5::h5createFile(file)

    assays <- assays(object, withDimnames = FALSE)
    layers <- setNames(paste0("/layers/", names(assays)), names(assays))
    layers[matrix] <- "/matrix"

    if (length(layers) > 1L)
        rhdf5::h5createGroup(file, "/layers")
    success <- unlist(Map(
        .export_loom, assays, layers, MoreArgs = list(file = file)
    ))
    if (!all(success == 0L))
        stop(
            "'.export_loom()' failed to write assay(s)\n  ",
            paste0(sQuote(names(layers)[success != 0]), collapse = ", ")
        )

    rhdf5::h5createGroup(file, "/col_attrs")
    .export_loom(colData(object), file, "col_attrs", colnames_attr)
    rhdf5::h5createGroup(file, "/row_attrs")
    if (is(object, "RangedSummarizedExperiment"))
        rowData <- rowRanges(object)
    else
        rowData <- rowData(object)
    .export_loom(rowData, file, "row_attrs", rownames_attr)

    rhdf5::h5createGroup(file, "col_graphs")
    rhdf5::h5createGroup(file, "row_graphs")

    invisible(file)
})
