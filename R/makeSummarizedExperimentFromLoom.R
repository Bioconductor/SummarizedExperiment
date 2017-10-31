.loom_make_rownames <- function(df, colname) {
    rownames(df) <- df[[colname]]
    df[, -match(colname, colnames(df)), drop = FALSE]
}

makeSummarizedExperimentFromLoom <-
    function(file, rownames_attr = NULL, colnames_attr = NULL)
{
    stopifnot(file.exists(file))

    ls <- rhdf5::h5ls(file)
    rowColnames <- ls[ls$group == "/row_attrs", "name", drop=TRUE]
    colColnames <- ls[ls$group == "/col_attrs", "name", drop=TRUE]
    stopifnot(
        is.null(rownames_attr) || rownames_attr %in% rowColnames,
        is.null(colnames_attr) || colnames_attr %in% colColnames
    )

    assay <- t(HDF5Array::HDF5Array(file, "/matrix"))
    layerNames <- ls[ls$group == "/layers", "name", drop = TRUE]
    layers <- lapply(setNames(layerNames, layerNames), function(layer) {
        layer <- paste0("/layers/", layer)
        t(HDF5Array::HDF5Array(file, layer))
    })
    assays <- c(list(matrix = assay), layers)

    rowData <- DataFrame(rhdf5::h5read(file, "row_attrs"))
    if (!is.null(rownames_attr))
        rowData <- .loom_make_rownames(rowData, rownames_attr)

    colData <- DataFrame(rhdf5::h5read(file, "col_attrs"))
    if (!is.null(colnames_attr))
        colData <- .loom_make_rownames(colData, colnames_attr)

    se <- SummarizedExperiment(assays, rowData = rowData, colData = colData)
    metadata(se) <- rhdf5::h5readAttributes(file, "/")
    se
}
