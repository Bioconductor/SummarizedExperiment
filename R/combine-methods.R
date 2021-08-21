# Contains methods for combineRows and combineCols. These serve as more
# fault-tolerant relaxed counterparts to rbind and cbind, respectively.

setMethod("combineRows", "SummarizedExperiment", function(x, ..., delayed=TRUE, fill=NA, use.names=TRUE) {
    all.se <- list(x, ...)

    # Combining the rowData.
    all.rd <- lapply(all.se, rowData)
    tryCatch({
        com.rd <- do.call(combineRows, all.rd)
    }, error=function(e) {
        stop(paste0("failed to combine rowData of SummarizedExperiment objects:\n  ", conditionMessage(e)))
    })

    # Combining the colData. This constructs mappings of the columns for each
    # SE to the columns of the final object.
    all.cd <- lapply(all.se, colData)
    tryCatch({
        com.cd <- do.call(combineUniqueCols, c(all.cd, list(use.names=use.names)))
    }, error=function(e) {
        stop(paste0("failed to combine colData of SummarizedExperiment objects:\n  ", conditionMessage(e)))
    })

    if (use.names) {
        # If combineUniqueCols succeeded, all SE's should have valid column names.
        all.names <- rownames(com.cd)
        mappings <- vector("list", length(all.se))
        for (i in seq_along(mappings)) {
            mappings[[i]] <- match(all.names, rownames(all.cd[[i]]))
        }
    } else {
        mappings <- NULL
    }

    args <- list(
        assays=combine_assays_by(all.se, mappings, delayed=delayed, fill=fill, by.row=TRUE),
        colData=com.cd,
        metadata=unlist(lapply(all.se, metadata), recursive=FALSE, use.names=FALSE)
    )

    # Finally, filling in the rowRanges. Rows for SummarizedExperiment
    # inputs are filled in with empty GRangesList objects.
    extracted <- extract_granges_from_se(all.se)

    if (!is.null(extracted)) {
        filled.ranges <- FALSE
        for (i in seq_along(extracted)) {
            if (is.null(extracted[[i]])) {
                filled.ranges <- TRUE
                cur.se <- all.se[[i]]
                levels <- rownames(cur.se)
                if (is.null(levels)) {
                    levels <- seq_len(nrow(cur.se))
                }
                rr <- splitAsList(GRanges(), factor(character(0), levels))
                if (is.null(rownames(cur.se))) {
                    names(rr) <- NULL
                }
                extracted[[i]] <- rr
            }
        }

        if (filled.ranges) {
            for (s in seq_along(extracted)) {
                extracted[[s]] <- as(extracted[[s]], "GRangesList")
            }
        }

        com.rr <- do.call(c, extracted)
        mcols(com.rr) <- com.rd
        args$rowRanges <- com.rr
    } else {
        args$rowData <- com.rd
    }

    # Assembling the SE.
    do.call(SummarizedExperiment, args)
})

combine_assays_by <- function(all.se, mappings, delayed, fill, by.row) {
    all.assays <- lapply(all.se, assays, withDimnames=FALSE)
    each.assay.names <- lapply(all.assays, names)
    no.assay.names <- vapply(each.assay.names, is.null, TRUE)

    if (by.row) {
        INFLATE <- inflate_matrix_by_column
    } else {
        INFLATE <- inflate_matrix_by_row
    }

    if (any(no.assay.names)) {
        if (!all(no.assay.names)) {
            stop("named and unnamed assays cannot be mixed")
        }
        n.assays <- unique(lengths(all.assays))
        if (length(n.assays)!=1L) {
            stop("all SummarizedExperiments should have the same number of unnamed assays")
        }

        for (s in seq_along(all.se)) {
            all.assays[[s]] <- lapply(all.assays[[s]], FUN=INFLATE,
                                      idx=mappings[[s]], delayed=delayed, fill=fill)
        }
    } else {
        all.assay.names <- Reduce(union, each.assay.names)
        for (s in seq_along(all.se)) {
            cur.se <- all.se[[s]]
            cur.assays <- all.assays[[s]]
            idx <- mappings[[s]]

            # Filling in all missing assay names and columns.
            for (a in all.assay.names) {
                if (a %in% names(cur.assays)) {
                    mat <- INFLATE(cur.assays[[a]], idx, delayed=delayed, fill=fill)
                } else {
                    nr <- nrow(cur.se)
                    nc <- ncol(cur.se)
                    if (!is.null(idx)) {
                        if (by.row) {
                            nc <- length(idx)
                        } else {
                            nr <- length(idx)
                        }
                    }
                    mat <- create_dummy_matrix(nr, nc, delayed=delayed, fill=fill)
                }
                cur.assays[[a]] <- mat
            }
            all.assays[[s]] <- cur.assays
        }
    }

    # Re-use assay r/cbind'ing machinery.
    all.assays <- lapply(all.assays, Assays)
    if (by.row) {
        combined <- do.call(rbind, all.assays)
    } else {
        combined <- do.call(cbind, all.assays)
    }

    combined <- as(combined, "SimpleList")

    # Prevent the SummarizedExperiment() constructor function from choking
    # on the rownames or colnames of the combined assays.
    if (by.row) {
        endoapply(combined, `colnames<-`, NULL)
    } else {
        endoapply(combined, `rownames<-`, NULL)
    }
}

create_dummy_matrix <- function(nr, nc, delayed, fill) {
    if (!delayed) {
        array(c(nr, nc), data=fill)
    } else {
        ConstantArray(c(nr, nc), value=fill)
    }
}

inflate_matrix_by_column <- function(mat, idx, delayed, fill) {
    if (delayed) {
        mat <- DelayedArray(mat)
    }
    if (!is.null(idx)) {
        absent <- is.na(idx)
        if (any(absent)) {
            idx[absent] <- ncol(mat)+1L
            mat <- cbind(mat, create_dummy_matrix(nrow(mat), 1L, delayed, fill))
        }
        mat <- mat[,idx,drop=FALSE]
    }
    mat
}

inflate_matrix_by_row <- function(mat, idx, delayed, fill) {
    if (delayed) {
        mat <- DelayedArray(mat)
    }
    if (!is.null(idx)) {
        absent <- is.na(idx)
        if (any(absent)) {
            idx[absent] <- nrow(mat)+1L
            mat <- rbind(mat, create_dummy_matrix(1L, ncol(mat), delayed, fill))
        }
        mat <- mat[idx,,drop=FALSE]
    }
    mat
}

extract_granges_from_se <- function(all.se) {
    has.ranges <- vapply(all.se, is, class2="RangedSummarizedExperiment", FUN.VALUE=TRUE)
    if (!any(has.ranges)) {
        return(NULL)
    }

    final.rr <- vector("list", length(all.se))
    for (s in which(has.ranges)) {
        cur.se <- all.se[[s]]
        rr <- rowRanges(cur.se)
        mcols(rr) <- NULL
        names(rr) <- rownames(cur.se)
        final.rr[[s]] <- rr
    }

    # Coercing everyone to a GRL if anyone is a GRL. Note that we don't fill in
    # NULLs with GRLs yet, to give a chance for the caller to decide how to
    # handle them (e.g., fill in combineRows or merge in combineCols).
    is.grl <- vapply(final.rr[has.ranges], function(x) is(x, "GRangesList"), TRUE)
    if (any(is.grl)) {
        for (s in which(has.ranges)) {
            final.rr[[s]] <- as(final.rr[[s]], "GRangesList")
        }
    }

    final.rr
}

setMethod("combineCols", "SummarizedExperiment", function(x, ..., delayed=TRUE, fill=NA, use.names=TRUE) {
    all.se <- list(x, ...)

    # Combining the rowData. This constructs mappings of the rows for each
    # SE to the columns of the final object.
    all.rd <- lapply(all.se, rowData)
    tryCatch({
        com.rd <- do.call(combineUniqueCols, c(all.rd, list(use.names=use.names)))
    }, error=function(e) {
        stop(paste0("failed to combine rowData of SummarizedExperiment objects:\n  ", conditionMessage(e)))
    })

    # Combining the colData.
    all.cd <- lapply(all.se, colData)
    tryCatch({
        com.cd <- do.call(combineRows, all.cd)
    }, error=function(e) {
        stop(paste0("failed to combine colData of SummarizedExperiment objects:\n  ", conditionMessage(e)))
    })

    if (use.names) {
        # If combineUniqueCols succeeded for the rowData, all SE's should have valid row names.
        all.names <- rownames(com.rd)
        mappings <- vector("list", length(all.se))
        for (i in seq_along(mappings)) {
            mappings[[i]] <- match(all.names, rownames(all.rd[[i]]))
        }
    } else {
        mappings <- NULL
    }

    args <- list(
        assays=combine_assays_by(all.se, mappings, delayed=delayed, fill=fill, by.row=FALSE),
        colData=com.cd,
        metadata=unlist(lapply(all.se, metadata), recursive=FALSE, use.names=FALSE)
    )

    com.rr <- merge_granges_from_se(all.se, mappings)
    if (!is.null(com.rr)) {
        mcols(com.rr) <- com.rd
        names(com.rr) <- rownames(com.rd)
        args$rowRanges <- com.rr
    } else {
        args$rowData <- com.rd
    }

    # Assembling the SE.
    do.call(SummarizedExperiment, args)
})

merge_granges_from_se <- function(all.se, mappings) {
    extracted <- extract_granges_from_se(all.se)
    if (is.null(extracted)) {
        return(NULL)
    }

    has.ranges <- which(!vapply(extracted, is.null, FALSE))
    extracted.ranges <- extracted[has.ranges]

    if (!is.null(mappings)) {
        # We concatenate everything to automatically merge the seqinfo. We
        # then create a container for the filling process.
        temp.rr <- do.call(c, extracted.ranges)
        names(temp.rr) <- NULL
        nentries <- lengths(extracted.ranges)
        starts <- cumsum(c(0L, nentries))

        com.rr <- temp.rr
        if (length(temp.rr) > 0) {
            com.rr <- temp.rr[rep(1L, length(mappings[[1]]))]
        }
        filled <- logical(length(com.rr))

        for (i in seq_along(extracted.ranges)) {
            idx <- mappings[[has.ranges[i]]]
            candidates <- temp.rr[starts[i] + seq_len(nentries[i])]

            available <- which(!is.na(idx))
            new.rr <- candidates[idx[available]]
            old.rr <- com.rr[available]

            existing <- filled[available]
            if (!identical(new.rr[existing], old.rr[existing])) {
                warning(wmsg("different 'rowRanges' for shared rows across SummarizedExperiment objects, ",
                             "ignoring 'rowRanges' for ", class(all.se[[i]])[1], " ", i))
                next
            }

            com.rr[available[!existing]] <- new.rr[!existing]
            filled[available[!existing]] <- TRUE
        }

        if (!all(filled)) {
            com.rr <- as(com.rr, "GRangesList")
            com.rr[!filled] <- GRangesList(GRanges())
        }
    } else {
        if (length(unique(extracted.ranges)) > 1) {
            warning(wmsg("'rowRanges' are not identical across input objects, ",
                          "using 'rowRanges' from the first object only"))
        }
        com.rr <- extracted.ranges[[1]]
    }

    com.rr
}
