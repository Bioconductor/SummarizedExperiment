# Contains methods for combineRows and combineCols. These serve as more
# fault-tolerant relaxed counterparts to rbind and cbind, respectively.

setMethod("combineRows", "SummarizedExperiment", function(x, y, ..., use.names=TRUE, delayed=TRUE, fill=NA) {
    all.se <- list(x, y, ...)

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
        assays=combine_assays_by_row(all.se, mappings, delayed=delayed, fill=fill),
        colData=com.cd,
        metadata=unlist(lapply(all.se, metadata), recursive=FALSE, use.names=FALSE)
    ) 
    
    com.rr <- combine_granges_from_se(all.se)
    if (!is.null(com.rr)) {
        mcols(com.rr) <- com.rd
        args$rowRanges <- com.rr
    } else {
        args$rowData <- com.rd
    }

    # Assembling the SE.
    do.call(SummarizedExperiment, args)
})

combine_assays_by_row <- function(all.se, mappings, delayed, fill) {
    all.assays <- lapply(all.se, assays, withDimnames=FALSE)
    each.assay.names <- lapply(all.assays, names)
    no.assay.names <- vapply(each.assay.names, is.null, TRUE)
    
    if (any(no.assay.names)) {
        if (!all(no.assay.names)) {
            stop("named and unnamed assays cannot be mixed")
        }
        n.assays <- unique(lengths(all.assays))
        if (length(n.assays)!=1L) {
            stop("all SummarizedExperiments should have the same number of unnamed assays")
        }
        for (s in seq_along(all.se)) {
            all.assays[[s]] <- lapply(all.assays[[s]], FUN=inflate_matrix_by_column, 
                                      j=mappings[[s]], delayed=delayed, fill=fill)
        }
    } else {
        all.assay.names <- Reduce(union, each.assay.names)
        for (s in seq_along(all.se)) {
            cur.se <- all.se[[s]]
            cur.assays <- all.assays[[s]]
            j <- mappings[[s]]

            # Filling in all missing assay names and columns.
            for (a in all.assay.names) {
                if (a %in% names(cur.assays)) {
                    mat <- inflate_matrix_by_column(cur.assays[[a]], j, delayed=delayed, fill=fill)
                } else {
                    if (is.null(j)) {
                        nc <- ncol(cur.se)
                    } else {
                        nc <- length(j)
                    }
                    mat <- create_dummy_matrix(nrow(cur.se), nc, delayed=delayed, fill=fill)
                }
                cur.assays[[a]] <- mat
            }
            all.assays[[s]] <- cur.assays
        }
    }

    # Re-use assay rbind'ing machinery.
    all.assays <- lapply(all.assays, Assays)
    combined <- do.call(rbind, all.assays)
    as(combined, "SimpleList") 
}

create_dummy_matrix <- function(nr, nc, delayed, fill) {
    if (!delayed) {
        array(c(nr, nc), data=fill)
    } else {
        ConstantArray(c(nr, nc), value=fill)
    }
}

inflate_matrix_by_column <- function(mat, j, delayed, fill) {
    if (delayed) {
        mat <- DelayedArray(mat)
    }
    if (!is.null(j)) {
        absent <- is.na(j)
        if (any(absent)) {
            j[absent] <- ncol(mat)+1L
            mat <- cbind(mat, create_dummy_matrix(nrow(mat), 1L, delayed, fill))
        }
        mat <- mat[,j,drop=FALSE]
    }
    mat
}

combine_granges_from_se <- function(all.se) {
    has.ranges <- vapply(all.se, is, class2="RangedSummarizedExperiment", FUN.VALUE=TRUE)
    if (!any(has.ranges)) {
        return(NULL)
    }

    as.grl <- !all(has.ranges)
    final.rr <- vector("list", length(all.se))
    for (s in seq_along(all.se)) {
        cur.se <- all.se[[s]]
        if (is(cur.se, "RangedSummarizedExperiment")) {
            rr <- rowRanges(cur.se)
            mcols(rr) <- NULL
            if (as.grl) {
                rr <- as(rr, "GRangesList")
            }
        } else {
            rr <- GRangesList(rep(list(GRanges()), nrow(cur.se)))
            rownames(rr) <- rownames(cur.se)
        }
        final.rr[[s]] <- rr
    }

    do.call(c, final.rr)
}
