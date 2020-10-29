setGeneric("mergeRows", function(...) standardGeneric("mergeRows"))

setMethod("mergeRows", 'SummarizedExperiment', function(...) {
    args <- list(...)

    # column names should be all equal or all NULL.
    all.cn <- lapply(args, colnames)
    if (length(unique(all.cn))!=1L) {
        stop("all entries of '...' must have the same colnames")
    }
    if (is.null(all.cn[[1]])) {
        all.nc <- lapply(args, ncol)
        if (length(unique(all.nc))!=1L) {
            stop("all entries of '...' must have the same number of columns")
        }
    }

    all.assays <- lapply(args, slot, name="assays")
    final.assays <- .relaxed_combine_assays(all.assays, byrow=TRUE)
    
    # Man, what a pain.
    is.rse <- vapply(args, is, class="RangedSummarizedExperiment", FUN.VALUE=TRUE)
    final.rr <- NULL
    if (all(is.rse)) {
        all.rr <- lapply(args, rowRanges)
        for (i in seq_along(all.rr)) {
            mcols(all.rr[[i]]) <- NULL
        }
        try(final.rr <- do.call(c, all.rr))
    }

    all.rd <- lapply(args, rowData)
    final.rd <- .relaxed_rbind_data_frames(all.rd)

    all.cd <- lapply(args, slot, name="colData")
    final.cd <- .relaxed_cbind_data_frames(all.cd)

    all.meta <- lapply(args, slot, name="metadata")
    final.meta <- do.call(c, all.meta)

    if (!all(is.rse)) {
        output <- BiocGenerics:::replaceSlots(args[[1]], assays=final.assays,
            elementMetadata=final.rd, colData=final.cd, metadata=final.meta, check=FALSE)
        if (!isTRUE(validObject(output, test=TRUE))) {
            output <- new("SummarizedExperiment", output)
        }
    } else {
        if (is.null(final.rr)) {
            final.rr <- GRangesList(rep(list(GRanges())))
            final.rr <- final.rr[rep(1, nrow(final.rd))]
            names(final.rr) <- rownames(final.rd)
        }
        mcols(final.rr) <- final.rd

        output <- BiocGenerics:::replaceSlots(args[[1]], assays=final.assays,
            rowRanges=final.rr, elementMetadata=final.rd[,0],
            colData=final.cd, metadata=final.meta, check=FALSE)

        if (!isTRUE(validObject(output, test=TRUE))) {
            output <- new("RangedSummarizedExperiment", output)
        }
    }

    output
})

.relaxed_combine_assays <- function(assay.list, byrow=FALSE) {
    # Cleaning out all unshared assays. Unnamed assays only supported
    # if all entries have the same number of them.

    all.an <- lapply(assay.list, names)
    common.an <- Reduce(intersect, all.an)
    common.an <- setdiff(common.an, "")

    if (length(common.an)) {
        assay.list <- lapply(assay.list, function(x) {
            x@data <- x@data[common.an]
            x
        })
    } else {
        if (length(unique(lengths(all.an)))!=1L || any(unlist(all.an)!="")) {
            stop("mix of named and unnamed assays not supported")
        }
    }

    if (byrow) {
        do.call(rbind, assay.list)
    } else {
        do.call(cbind, assay.list)
    }
}

# NOTE: I'd like this to be exported for SCE's use.
.relaxed_rbind_data_frames <- function(df.list) {
    all.names <- lapply(df.list, colnames)
    common.rd <- Reduce(intersect, all.names)
    collated.rd <- vector("list", length(common.rd))
    names(collated.rd) <- common.rd

    for (x in common.rd) {
        copy <- df.list
        for (y in seq_along(copy)) {
            copy[[y]] <- df.list[[y]][,x,drop=FALSE]
        }
        try({
            collated.rd[[x]] <- do.call(rbind, copy)
        })
    }

    failed <- vapply(collated.rd, is.null, TRUE)
    output <- collated.rd[!failed]

    if (length(output)==0) {
        output <- lapply(df.list, function(x) x[,0,drop=FALSE])
    } 
    do.call(rbind, output)
}

# NOTE: I'd also like this to be exported for SCE's use.
.relaxed_cbind_data_frames <- function(df.list) {
    all.names <- lapply(df.list, colnames)
    by.names <- split(rep(seq_along(all.names), lengths(all.names)), unlist(all.names))

    okay <- character(0)
    for (field in names(by.names)) {
        possibles <- by.names[[field]]
        ref <- df.list[[possibles[1]]][,field]

        failed <- FALSE
        for (i in seq_along(possibles)[-1]) {
            if (!identical(ref, df.list[[possibles[i]]][,field])) {
                failed <- TRUE
            }
        }

        if (!failed) {
            okay <- c(okay, field)
        }
    }

    to.keep <- lapply(by.names[okay], head, n=1)
    for (i in seq_along(to.keep)) {
        to.keep[[i]] <- df.list[[i]][,names(to.keep)[i],drop=FALSE]
    }

    if (length(to.keep)) {
        do.call(cbind, to.keep)
    } else {
        df.list[[1]][,0]
    }
}
